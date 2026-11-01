"""Discover and check the original four-strand X_3 cofactor over Q(q,t).

Run locally with Singular 4.4 and freegb.lib installed::

    python examples/x3_singular_cofactor.py --singular Singular

All 31 original Bigelow presentation relations are used. An interpreted
completion and tracing script uses Singular's free-polynomial arithmetic,
stops when X_3 reduces to zero, and expands that proof into an explicit sum
of multiples of the original relations.
The saved cofactor is reloaded and checked in the untruncated free algebra.
No existing basis, cofactor, or support is an input. A resource-limit failure
is inconclusive. This research computation is deliberately excluded from CI.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from zipper_algebra import Presentation
from zipper_algebra.algebra import AlgebraElement, AlgebraPresentation, CofactorIdentity, FreeAlgebra
from zipper_algebra.coefficients import SymbolicRing
from zipper_algebra.cofactor.general import expand_derivation, verify_cofactor
from zipper_algebra.execution.artifacts import prepare_output_root, validate_source_snapshot
from zipper_algebra.provenance import capture_provenance, file_identity, write_bundle_provenance
from zipper_algebra.verification.derivations import DerivationGraph, DerivationTerm


HELPER = Path(__file__).with_name('four_strand_square_zero.py')
ENGINE = HELPER.with_suffix('.sing')


@lru_cache(maxsize=1)
def load_helper():
    spec = importlib.util.spec_from_file_location('x3_singular_trace_helper', HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def presentation():
    """Return the literal default four-strand presentation and its X_3 target."""
    helper = load_helper()
    names, relations, sliders, tags = Presentation(4).build(SymbolicRing())
    field = helper.FIELD
    algebra = FreeAlgebra(field, tuple(names))
    q, t = (field.parameter(name) for name in ('q', 't'))

    def convert(polynomial):
        return algebra.element({
            word: sum((c*q**i*t**j for (i, j), c in coefficient.v.items()), field.zero)
            for word, coefficient in polynomial.terms.items()})

    full = AlgebraPresentation(algebra, tuple(map(convert, relations)),
                               tuple(tag.label(4) for tag in tags))
    return full, convert(sliders[(3, 1)])


def target_engine():
    """Adapt the interpreted completion/tracer to stop at a proved target node."""
    source = ENGINE.read_text(encoding='utf-8')
    check = '''
int targetDone=0;
proc checkTarget()
{
 if(reduce(Target,B)!=0){return(0);}
 list terms;list res=traceNF(Target,terms);
 if(res[1]!=0){return(0);}
 emitDerived(Target,res[2],nodeCount+1);
 print("TARGET_DONE|"+string(nodeCount+1));
 return(1);
}
'''
    replacements = {
        'proc completeToHeads()': check+'\nproc completeToHeads()',
        'while(size(Missing)>0 && nodeCount<maxNodes)':
            'targetDone=checkTarget();\nwhile(targetDone==0 && nodeCount<maxNodes)',
        'if(result[1]!=0){insertNode(result[1],result[2]);}':
            'if(result[1]!=0){insertNode(result[1],result[2]);targetDone=checkTarget();}',
    }
    for old, new in replacements.items():
        if source.count(old) != 1:
            raise ValueError('Singular trace engine differs from the supported source shape')
        source = source.replace(old, new)
    return source


def source_for(full, target, library, *, degree, max_nodes):
    """Supply every original relation, keeping the target outside the input ideal."""
    helper = load_helper()
    library_name = str(library)
    if any(character in library_name for character in ('\\', '"', '\n', '\r')):
        raise ValueError('unsupported Singular library path')
    return (f'LIB "{library_name}";\n'
            'ring r=(0,q,t),(y3,y2,y1,x3,x2,x1),(c,Dp);\n'
            f'ring R=freeAlgebra(r,{degree});\nshort=0;\nideal I='+
            ',\n'.join(helper.polynomial_source(dict(relation.terms))
                        for relation in full.relations)+';\npoly Target='+
            helper.polynomial_source(dict(target.terms))+';\n'+
            f'int maxNodes={max_nodes};int maxDegree={degree};\n'
            'option(noredSB);option(noredTail);\n'+target_engine()+
            '\ncompleteToHeads();quit;\n')


def parse_target(full, target, text, *, max_nodes):
    """Accept one framed target trace; its expansion must still be checked."""
    markers = re.findall(r'^TARGET_DONE\|([1-9][0-9]*)$', text, re.MULTILINE)
    if len(markers) != 1 or sum(line.startswith('TARGET_DONE') for line in text.splitlines()) != 1:
        raise RuntimeError('Singular did not produce one target completion marker within its limits')
    polynomials, derivations = load_helper().parse_output(text, inputs=len(full.relations))
    if len(polynomials) > max_nodes+1 or set(derivations) != set(polynomials):
        raise ValueError('Singular proof exceeds its node bound or has missing derivations')
    nodes = tuple(full.algebra.element(polynomials[index]) for index in range(1, len(polynomials)+1))
    rows = tuple(tuple(DerivationTerm(c, left, 'input' if kind == '0' else 'node', reference-1, right)
                       for c, left, kind, reference, right in derivations[index])
                 for index in range(1, len(polynomials)+1))
    node = int(markers[0])-1
    if node != len(nodes)-1 or nodes[node] != target:
        raise ValueError('Singular target marker does not identify the final target polynomial')
    return DerivationGraph(full.content_hash(), nodes, rows), node


def produce(binary, output, *, timeout=120, memory_mib=4096, max_nodes=1000,
            degree=24, max_operations=10000000, max_contexts=100000):
    """Discover, expand, publish, reload and independently check one cofactor.

    Timeout and virtual-memory caps apply to each native stage. Python expansion
    and verification have explicit scalar-operation and retained-context bounds.
    The degree bound applies to Singular's search, never the final exact check.
    """
    integers = (memory_mib, max_nodes, degree, max_operations, max_contexts)
    if (type(timeout) not in (int, float) or not 0 < timeout < float('inf')
            or any(type(value) is not int for value in integers)
            or memory_mib < 256 or max_nodes < 31 or not 12 <= degree <= 24
            or min(max_operations, max_contexts) <= 0):
        raise ValueError('invalid resource limits')
    binary, output = Path(binary).resolve(), Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError('output directory must be empty; use a fresh run directory')
    started = time.monotonic()
    prepare_output_root(output)
    producer = capture_provenance()
    paths = (Path(__file__), HELPER, ENGINE, binary)
    snapshots = {path: file_identity(path) for path in paths}
    helper = load_helper()
    full, target = presentation()
    (output/'presentation.json').write_text(json.dumps(full.to_json(), indent=2)+'\n', encoding='utf-8')
    (output/'target.json').write_text(json.dumps(target.to_json(), indent=2)+'\n', encoding='utf-8')
    library, runtime, environment_seconds = helper.solver_environment(
        binary, output, timeout=timeout, memory_mib=memory_mib)
    snapshots[library] = file_identity(library)
    source = source_for(full, target, library, degree=degree, max_nodes=max_nodes)
    text, native_seconds = helper.run_singular(binary, source, output, 'producer',
                                               timeout=timeout, memory_mib=memory_mib)
    graph, node = parse_target(full, target, text, max_nodes=max_nodes)
    (output/'graph.json').write_text(json.dumps(graph.to_json(), separators=(',', ':'))+'\n',
                                     encoding='utf-8')
    before = time.monotonic()
    identity = expand_derivation(full, graph, node, max_operations=max_operations,
                                  max_contexts=max_contexts, max_nodes=max_nodes+1)
    expansion_seconds = time.monotonic()-before
    candidate = output/'.cofactor-candidate.json'
    candidate.write_text(json.dumps(identity.to_json(), separators=(',', ':'))+'\n', encoding='utf-8')
    before = time.monotonic()
    reloaded = CofactorIdentity.from_json(json.loads(candidate.read_text(encoding='utf-8')))
    saved_presentation = AlgebraPresentation.from_json(json.loads((output/'presentation.json').read_text()))
    saved_target = AlgebraElement.from_json(json.loads((output/'target.json').read_text()))
    if (saved_presentation != full or saved_target != target or reloaded.target != saved_target
            or not verify_cofactor(saved_presentation, reloaded, max_operations=max_operations)):
        raise ValueError('published original-input cofactor failed exact verification')
    verification_seconds = time.monotonic()-before
    validate_source_snapshot(producer)
    if any(file_identity(path) != snapshot for path, snapshot in snapshots.items()):
        raise RuntimeError('producer source, solver, or library changed during the run')
    candidate.replace(output/'cofactor.json')
    report = {
        'schema': 'zipper-x3-singular-cofactor-run-v1', 'presentation_hash': full.content_hash(),
        'target': 'X_3', 'field': full.algebra.field.to_json(), 'runtime': runtime,
        'input_relations': len(full.relations), 'input_relation_labels': list(full.labels),
        'method': 'interpreted completion in Singular stopping at target reducibility, with original-input proof expansion',
        'proof_nodes': len(graph.nodes), 'target_node': node,
        'local_proof_terms': sum(map(len, graph.derivations)), 'cofactor_terms': len(reloaded.terms),
        'verification': 'saved cofactor reloaded and checked exactly against all original inputs',
        'limits': {'timeout_seconds_per_native_stage': timeout, 'native_virtual_memory_mib': memory_mib,
                   'maximum_search_nodes': max_nodes, 'maximum_search_degree': degree,
                   'maximum_scalar_operations_per_exact_stage': max_operations,
                   'maximum_retained_contexts': max_contexts},
        'seconds': {'environment': environment_seconds, 'native': native_seconds,
                    'expansion_and_check': expansion_seconds, 'reload_and_check': verification_seconds,
                    'total_before_provenance': time.monotonic()-started}}
    (output/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    inputs = [{'name': path.name, **snapshot} for path, snapshot in snapshots.items()]
    write_bundle_provenance(output, sorted(path.name for path in output.iterdir() if path.is_file()),
                            producer=producer, inputs=inputs)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--singular', default='Singular', help='path or command for Singular 4.4')
    parser.add_argument('--output', type=Path, default=Path('runs/x3-singular-cofactor'))
    parser.add_argument('--timeout', type=float, default=120, help='seconds per native solver stage')
    parser.add_argument('--memory-mib', type=int, default=4096, help='native virtual-memory cap')
    parser.add_argument('--max-nodes', type=int, default=1000)
    parser.add_argument('--degree', type=int, default=24)
    parser.add_argument('--max-operations', type=int, default=10000000)
    parser.add_argument('--max-contexts', type=int, default=100000)
    args = parser.parse_args(argv)
    executable = shutil.which(args.singular)
    if executable is None:
        parser.error(f'Singular executable not found: {args.singular}')
    try:
        report = produce(Path(executable), args.output, timeout=args.timeout, memory_mib=args.memory_mib,
                         max_nodes=args.max_nodes, degree=args.degree,
                         max_operations=args.max_operations, max_contexts=args.max_contexts)
    except (OSError, ValueError, TypeError, RuntimeError, ArithmeticError, subprocess.SubprocessError) as error:
        parser.exit(1, f'{error}\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
