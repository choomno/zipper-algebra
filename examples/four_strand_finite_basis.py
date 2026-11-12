"""Produce finite-word-basis and generator-action certificates over Q(q,t).

First run four_strand_square_zero.py, then run this local computation::

    python examples/four_strand_finite_basis.py \
        --source-run runs/four-strand-square-zero \
        --output runs/four-strand-finite-basis

The source run supplies literal presentation and proved-polynomial data. This
producer appends compact ideal identities for every generator action. Its JSON
output is a candidate WordBasisCertificate; verify it before using its dimension
or action matrices as established results. Full research runs stay outside CI.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil

from zipper_algebra.algebra import AlgebraPresentation
from zipper_algebra.execution.artifacts import prepare_output_root, validate_source_snapshot
from zipper_algebra.provenance import (OUTPUT_SCHEMA, capture_provenance, file_identity,
                                       validate_provenance, write_bundle_provenance)
from zipper_algebra.verification.derivations import DerivationGraph, DerivationTerm
from zipper_algebra.verification.groebner import GroebnerCertificate


ENGINE = Path(__file__).with_suffix('.sing')
HELPER = Path(__file__).with_name('four_strand_square_zero.py')
SOURCE_NAMES = ('presentation.json', 'certificate.json', 'proof.txt')
INPUT_BUDGET = 32*1024*1024


def load_helper():
    """Load only the adjacent, named reproduction recipe, without path searches."""
    spec = importlib.util.spec_from_file_location('four_strand_square_zero_recipe', HELPER)
    if spec is None or spec.loader is None:
        raise RuntimeError('the adjacent four-strand producer helper could not be loaded')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate JSON key')
        result[key] = value
    return result


def strict_loads(text):
    def nonfinite(value):
        raise ValueError('nonfinite JSON values are not permitted')
    return json.loads(text, object_pairs_hook=_unique_object, parse_constant=nonfinite)


def read_input_bytes(path):
    """Enforce the input limit on the bytes read, including a growing file."""
    with path.open('rb') as stream:
        content = stream.read(INPUT_BUDGET+1)
    if len(content) > INPUT_BUDGET:
        raise ValueError('source input exceeds the reproduction input budget')
    return content


def read_json(path):
    return strict_loads(read_input_bytes(path).decode('utf-8'))


def input_identity(path, *, bounded=False):
    if not bounded:
        return file_identity(path)
    content = read_input_bytes(path)
    return {'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()}


def validate_inputs(snapshots, *, absent_provenance=None):
    for path, identity, bounded in snapshots:
        try:
            current = input_identity(path, bounded=bounded)
        except OSError as error:
            raise ValueError(f'input became unreadable during the computation: {path.name}') from error
        if current != identity:
            raise ValueError(f'input changed during the computation: {path.name}')
    if absent_provenance is not None and absent_provenance.exists():
        raise ValueError('source provenance appeared during the computation')


def source_provenance_status(source_run, snapshots, *, names=SOURCE_NAMES):
    """Verify byte links to consumed inputs, without implying mathematical validity."""
    path = source_run/'provenance.json'
    if not path.is_file():
        return 'absent'
    record = read_json(path)
    if (not isinstance(record, dict) or set(record) != {'schema', 'producer', 'outputs', 'inputs'}
            or record['schema'] != OUTPUT_SCHEMA or not isinstance(record['outputs'], dict)
            or not isinstance(record['inputs'], list)):
        raise ValueError('invalid source output provenance')
    validate_provenance(record['producer'])
    identities = {item.name: identity for item, identity, _ in snapshots if item.parent == source_run}
    for name in names:
        if record['outputs'].get(name) != identities[name]:
            raise ValueError(f'source provenance does not match consumed file: {name}')
    return 'present_verified'


def _index(text, bound):
    if not re.fullmatch(r'[1-9][0-9]*', text) or int(text) > bound:
        raise ValueError('framed index is outside its declared range')
    return int(text)-1


def load_source(directory, helper):
    """Bind all three source encodings to the same full generic presentation."""
    full = AlgebraPresentation.from_json(read_json(directory/'presentation.json'))
    if full != helper.presentation():
        raise ValueError('source presentation is not the full generic four-strand square-zero presentation')
    certificate = GroebnerCertificate.from_json(read_json(directory/'certificate.json'))
    if certificate.presentation_hash != full.content_hash():
        raise ValueError('source certificate presentation hash mismatch')
    if any(polynomial.algebra != full.algebra for polynomial in certificate.basis):
        raise ValueError('source certificate uses a different free algebra')
    path = directory/'proof.txt'
    lines = read_input_bytes(path).decode('utf-8').splitlines()
    if not lines or not re.fullmatch(r'META\|[1-9][0-9]*\|[1-9][0-9]*\|[1-9][0-9]*', lines[0]):
        raise ValueError('invalid source proof framing')
    counts = tuple(int(value) for value in lines[0].split('|')[1:])
    inputs, candidates, nodes = counts
    if (inputs != len(full.relations) or nodes != len(certificate.basis)
            or candidates > nodes or lines[-1] != 'END|'+ '|'.join(map(str, counts))):
        raise ValueError('source proof counts or terminal frame disagree')
    polynomials = {'I': defaultdict(dict), 'G': defaultdict(dict)}
    derivations = defaultdict(list)
    for line in lines[1:-1]:
        parts = line.split('|')
        if len(parts) != (7 if parts[0] == 'D' else 6) or parts[0] not in ('I', 'G', 'D'):
            raise ValueError('unexpected source proof record')
        kind = parts[0]
        index = _index(parts[1], inputs if kind == 'I' else nodes)
        coefficient = helper.parse_coefficient(parts[2])
        if coefficient.is_zero():
            raise ValueError('source proof contains a zero sparse term')
        left = helper.parse_word(parts[3])
        if kind in ('I', 'G'):
            target = polynomials[kind][index]
            if parts[4:] != ['0', ''] or left in target:
                raise ValueError('invalid or duplicate source polynomial term')
            target[left] = coefficient
        else:
            if parts[4] not in ('0', '1'):
                raise ValueError('invalid source derivation reference kind')
            reference = _index(parts[5], inputs if parts[4] == '0' else nodes)
            right = helper.parse_word(parts[6])
            derivations[index].append((coefficient, left, 'input' if parts[4] == '0' else 'basis', reference, right))
    for kind, expected in (('I', full.relations), ('G', certificate.basis)):
        if set(polynomials[kind]) != set(range(len(expected))):
            raise ValueError('source proof has missing polynomials')
        if any(full.algebra.element(polynomials[kind][index]) != polynomial
               for index, polynomial in enumerate(expected)):
            raise ValueError('source proof polynomials disagree with source JSON')
    for index, row in enumerate(certificate.derivations):
        expected = [(term.coefficient, term.left, term.kind, term.index, term.right) for term in row]
        if derivations[index] != expected:
            raise ValueError('source proof derivations disagree with source JSON')
    return full, certificate, candidates


def standard_words(candidate, *, maximum=1000):
    """Enumerate the candidate's irreducible words in breadth-first degree-lex order."""
    if not candidate:
        raise ValueError('a nonempty candidate is required')
    size = len(candidate[0].algebra.generators)
    heads = [max((word for word, _ in polynomial.terms), key=lambda word: (len(word), word))
             for polynomial in candidate]
    if () in heads:
        raise ValueError('the four-strand example requires a nonzero quotient')
    words = [()]
    queue = deque([()])
    while queue:
        word = queue.popleft()
        for letter in range(1, size+1):
            child = word+(letter,)
            if any(any(child[i:i+len(head)] == head for i in range(len(child)-len(head)+1)) for head in heads):
                continue
            words.append(child)
            if len(words) > maximum:
                raise ValueError('candidate standard words exceed the enumeration budget')
            queue.append(child)
    return tuple(words)


def action_source(full, source, candidate_count, words, helper, library, runtime, *, degree):
    """Use fresh candidate data and cached leading words in the Singular engine."""
    prefix, _ = helper.header(full, degree=degree, library=library,
                              library_version=runtime['freegb_version'])
    prefix += 'ideal B='+',\n'.join(helper.polynomial_source(dict(polynomial.terms))
                                     for polynomial in source.basis[:candidate_count])+';\n'
    encoded = ','.join('"'+''.join(str(7-letter) for letter in word)+'"' for word in words)
    prefix += f'list Standards=list({encoded});int baseNodes={len(source.basis)};\n'
    return prefix+ENGINE.read_text(encoding='utf-8')


def parse_actions(text, full, source, candidate_count, words, helper):
    """Validate framing and bind every exported action target to its sparse column."""
    from zipper_algebra.verification.finite import GeneratorAction, WordBasisCertificate

    dimension, generators, base = len(words), len(full.algebra.generators), len(source.basis)
    expected = generators*dimension
    total = base+expected
    word_indices = {word: index for index, word in enumerate(words)}
    markers = set()
    targets, rows, columns = defaultdict(dict), defaultdict(list), defaultdict(dict)
    zero = set()
    for line in text.splitlines():
        if line.startswith('TRACE_MISMATCH|'):
            raise ValueError('Singular traced normal form disagrees with builtin reduction: '+line)
        if not line.startswith(('ACTION|', 'G|', 'D|', 'Z|', 'A|')):
            continue
        parts = line.split('|')
        kind = parts[0]
        if kind in ('ACTION', 'A'):
            if len(parts) != (4 if kind == 'ACTION' else 5):
                raise ValueError('malformed generator-action record')
            generator, column = _index(parts[1], generators), _index(parts[2], dimension)
            node = base+generator*dimension+column
            if kind == 'ACTION':
                if _index(parts[3], total) != node or (generator, column) in markers:
                    raise ValueError('duplicate or mismatched action node marker')
                markers.add((generator, column))
            else:
                coefficient = helper.parse_coefficient(parts[3])
                word = helper.parse_word(parts[4])
                if coefficient.is_zero() or word not in word_indices:
                    raise ValueError('action coefficient is zero or its output word is outside the basis')
                row = word_indices[word]
                if row in columns[generator, column]:
                    raise ValueError('duplicate sparse action entry')
                columns[generator, column][row] = coefficient
            continue
        if len(parts) != {'G': 6, 'D': 7, 'Z': 2}[kind]:
            raise ValueError('malformed action proof record')
        node = _index(parts[1], total)
        if node < base:
            raise ValueError('action stage attempted to modify a source proof node')
        if kind == 'Z':
            if node in zero:
                raise ValueError('duplicate zero target marker')
            zero.add(node)
            continue
        coefficient = helper.parse_coefficient(parts[2])
        left = helper.parse_word(parts[3])
        if coefficient.is_zero():
            raise ValueError('action proof sparse coefficients must be nonzero')
        if kind == 'G':
            if parts[4:] != ['0', ''] or left in targets[node]:
                raise ValueError('invalid or duplicate action target term')
            targets[node][left] = coefficient
        else:
            if parts[4] != '1':
                raise ValueError('action traces must reference proved candidate rows')
            reference = _index(parts[5], candidate_count)
            right = helper.parse_word(parts[6])
            rows[node].append(DerivationTerm(coefficient, left, 'node', reference, right))
    completion = [line for line in text.splitlines() if line.startswith('ACTIONS_DONE|')]
    if len(markers) != expected or completion != [f'ACTIONS_DONE|{expected}']:
        raise ValueError('action stage is missing a column or its completion marker')
    if zero & set(targets) or zero | set(targets) != set(range(base, total)):
        raise ValueError('each action target needs exactly one polynomial or zero marker')
    nodes = list(source.basis)
    derivations = [tuple(DerivationTerm(term.coefficient, term.left,
                                      'node' if term.kind == 'basis' else 'input', term.index, term.right)
                         for term in row) for row in source.derivations]
    actions = []
    for generator in range(generators):
        action_columns = []
        for column, word in enumerate(words):
            node = base+generator*dimension+column
            entries = columns[generator, column]
            normal = full.algebra.element({words[row]: coefficient for row, coefficient in entries.items()})
            target = full.algebra.element({(generator+1,)+word: 1})-normal
            if full.algebra.element(targets[node]) != target:
                raise ValueError('action target disagrees with its literal generator column')
            nodes.append(target)
            derivations.append(tuple(rows[node]))
            action_columns.append(tuple(sorted(entries.items())))
        actions.append(GeneratorAction(full.algebra.field, tuple(action_columns)))
    graph = DerivationGraph(full.content_hash(), tuple(nodes), tuple(derivations))
    action_nodes = tuple(tuple(base+generator*dimension+column for column in range(dimension))
                         for generator in range(generators))
    return WordBasisCertificate(full.content_hash(), words, tuple(actions), graph, action_nodes, None)


def produce(binary, source_run, output, *, timeout=180, memory_mib=6144, degree=24):
    """Extend existing generic proof data; verification is a separate bounded step."""
    if timeout <= 0 or memory_mib < 256 or not 12 <= degree <= 48:
        raise ValueError('invalid resource limits')
    if output.exists() and any(output.iterdir()):
        raise ValueError('output directory must be empty; use a fresh run directory')
    producer = capture_provenance()
    paths = [Path(__file__), ENGINE, HELPER, binary, *(source_run/name for name in SOURCE_NAMES)]
    provenance_path = source_run/'provenance.json'
    supplied_provenance = provenance_path.is_file()
    if supplied_provenance:
        paths.append(provenance_path)
    snapshots = [(path, input_identity(path, bounded=path.parent == source_run), path.parent == source_run)
                 for path in paths]
    inputs = [{'name': path.name, **identity} for path, identity, _ in snapshots]
    source_status = source_provenance_status(source_run, snapshots)
    absent_provenance = None if supplied_provenance else provenance_path
    helper = load_helper()
    full, source, candidate_count = load_source(source_run, helper)
    words = standard_words(source.basis[:candidate_count])
    if max(map(len, words))+1 > degree:
        raise ValueError('action words exceed the configured letterplace degree')
    prepare_output_root(output)
    library, runtime, setup_seconds = helper.solver_environment(
        binary, output, timeout=timeout, memory_mib=memory_mib)
    inputs.append(runtime['freegb_source'])
    snapshots.append((library, {key: runtime['freegb_source'][key] for key in ('bytes', 'sha256')}, False))
    source_text = action_source(full, source, candidate_count, words, helper, library, runtime, degree=degree)
    text, seconds = helper.run_singular(binary, source_text, output, 'actions',
                                       timeout=timeout, memory_mib=memory_mib)
    certificate = parse_actions(text, full, source, candidate_count, words, helper)
    # Exercise the strict wire encoding before committing the generated file.
    wire = json.dumps(certificate.to_json(), separators=(',', ':'), allow_nan=False)+'\n'
    if type(certificate).from_json(strict_loads(wire)) != certificate:
        raise ArithmeticError('finite-basis certificate JSON round-trip failed')
    validate_source_snapshot(producer)
    validate_inputs(snapshots, absent_provenance=absent_provenance)
    (output/'presentation.json').write_text(json.dumps(full.to_json(), indent=2)+'\n', encoding='utf-8')
    (output/'word-basis-certificate.json').write_text(wire, encoding='utf-8')
    frames = [f'W|{index}|'+''.join(f'{letter},' for letter in word) for index, word in enumerate(words, 1)]
    frames.extend(line for line in text.splitlines() if line.startswith(('ACTION|', 'G|', 'D|', 'Z|', 'A|', 'ACTIONS_DONE|')))
    (output/'actions.txt').write_text('\n'.join(frames)+'\n', encoding='utf-8')
    report = {'schema': 'zipper-four-strand-finite-basis-run-v1',
              'presentation_hash': full.content_hash(), 'runtime': runtime,
              'source_provenance_status': source_status,
              'candidate_basis_words': len(words), 'source_proof_nodes': len(source.basis),
              'action_columns': len(full.algebra.generators)*len(words),
              'matrix_nonzero_entries': sum(len(column) for action in certificate.actions for column in action.columns),
              'proof_nodes': len(certificate.proof.nodes),
              'zero_action_targets': sum(node.is_zero() for node in certificate.proof.nodes[len(source.basis):]),
              'local_action_proof_terms': sum(len(row) for row in certificate.proof.derivations[len(source.basis):]),
              'verification': 'candidate produced; finite-basis verification required',
              'limits': {'timeout_seconds_per_stage': timeout, 'memory_mib': memory_mib, 'maximum_word_degree': degree},
              'seconds': {'environment': setup_seconds, 'actions': seconds}}
    (output/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    validate_source_snapshot(producer)
    validate_inputs(snapshots, absent_provenance=absent_provenance)
    write_bundle_provenance(output, sorted(path.name for path in output.iterdir() if path.is_file()),
                            producer=producer, inputs=inputs)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run', type=Path, default=Path('runs/four-strand-square-zero'))
    parser.add_argument('--singular', default='Singular')
    parser.add_argument('--output', type=Path, default=Path('runs/four-strand-finite-basis'))
    parser.add_argument('--timeout', type=int, default=180)
    parser.add_argument('--memory-mib', type=int, default=6144)
    parser.add_argument('--degree', type=int, default=24)
    args = parser.parse_args(argv)
    binary = shutil.which(args.singular)
    if binary is None:
        parser.error(f'Singular executable not found: {args.singular}')
    try:
        report = produce(Path(binary).resolve(), args.source_run, args.output,
                         timeout=args.timeout, memory_mib=args.memory_mib, degree=args.degree)
    except (ValueError, RuntimeError, ArithmeticError, OSError, UnicodeError) as error:
        parser.exit(1, f'{error}\n')
    print(f"Produced {report['action_columns']} action columns and {report['proof_nodes']} proof nodes "
          f"in {args.output}; run finite-basis verification.")
    return report


if __name__ == '__main__':
    main()
