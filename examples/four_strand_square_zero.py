"""Produce an exact generic four-strand square-zero certificate over Q(q,t).

Run locally, with Singular 4.4 and its freegb library installed::

    python examples/four_strand_square_zero.py --singular Singular

The full presentation contains 31 relations. An explicitly selected 30-relation
subset (omitting T_4) supplies the completion and local ideal-membership proofs.
Verification must check all 31 relations, thereby also establishing that T_4 is
redundant. Neither a precomputed basis nor specialization data is used.

This is a bounded research computation, excluded from CI. The producer emits a
candidate, not a trusted dimension assertion. Replay proof.txt with the
independent GAP checker before interpreting the expected dimension 186.
"""
from __future__ import annotations

import argparse
import ast
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction
import json
import os
from pathlib import Path
import re
import resource
import shutil
import signal
import subprocess
import time

from zipper_algebra import Presentation
from zipper_algebra.algebra import AlgebraPresentation, CofactorIdentity, FreeAlgebra
from zipper_algebra.coefficients import SymbolicRing
from zipper_algebra.execution.artifacts import prepare_output_root
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.provenance import capture_provenance, file_identity, write_bundle_provenance
from zipper_algebra.verification.groebner import GroebnerCertificate, ProofTerm


FIELD = RationalFunctionField(('q', 't'))
NAMES = ('x1', 'x2', 'x3', 'y1', 'y2', 'y3')
ENGINE = Path(__file__).with_suffix('.sing')


def presentation():
    """Literal full legacy square-zero presentation, with exact field coefficients."""
    source = Presentation(4, 'square-zero')
    _, relations, _, tags = source.build(SymbolicRing())
    algebra = FreeAlgebra(FIELD, NAMES)
    q, t = (FIELD.parameter(name) for name in FIELD.parameters)
    converted = tuple(algebra.element({
        word: sum((c*q**i*t**j for (i, j), c in coefficient.v.items()), FIELD.zero)
        for word, coefficient in relation.terms.items()}) for relation in relations)
    return AlgebraPresentation(algebra, converted, tuple(tag.label(4) for tag in tags))


@lru_cache(maxsize=32768)
def parse_coefficient(text):
    """Parse exact arithmetic without eval, attributes, calls, or implicit names."""
    if len(text) > 100_000 or '**' in text or not re.fullmatch(r'[qt0-9+*/^() -]+', text):
        raise ValueError('invalid coefficient syntax')
    try:
        tree = ast.parse(text.strip().replace('^', '**'), mode='eval')
    except (SyntaxError, RecursionError) as error:
        raise ValueError('invalid coefficient syntax') from error
    if sum(1 for _ in ast.walk(tree)) > 20_000:
        raise ValueError('coefficient expression exceeds the parser budget')

    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return FIELD.scalar(node.value)
        if isinstance(node, ast.Name) and node.id in FIELD.parameters:
            return FIELD.parameter(node.id)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = visit(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp):
            if isinstance(node.op, ast.Pow):
                if (not isinstance(node.right, ast.Constant) or type(node.right.value) is not int
                        or not 0 <= node.right.value <= 128):
                    raise ValueError('coefficient exponent must be an integer from zero to 128')
                return visit(node.left)**node.right.value
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
        raise ValueError('invalid coefficient expression')

    try:
        return visit(tree.body)
    except RecursionError as error:
        raise ValueError('coefficient expression exceeds the parser depth') from error


def parse_word(text):
    if not re.fullmatch(r'(?:[1-6],)*', text):
        raise ValueError('invalid six-generator word')
    word = tuple(int(letter) for letter in text.split(',') if letter)
    if len(word) > 48:
        raise ValueError('word exceeds the producer degree budget')
    return word


def parse_output(text, *, inputs=31):
    """Read framed sparse polynomials and acyclic local proofs from a solver log.

    Ordinary Singular progress text is ignored. Malformed framed records,
    duplicate terms, zero terms, noncontiguous nodes, and forward references are
    rejected. The runner separately requires the appropriate completion marker.
    """
    polynomials, derivations = defaultdict(dict), defaultdict(list)
    for line in text.splitlines():
        if not line.startswith(('G|', 'D|')):
            continue
        parts = line.split('|')
        if len(parts) != (6 if parts[0] == 'G' else 7):
            raise ValueError('malformed solver record')
        if not re.fullmatch(r'[1-9][0-9]*', parts[1]):
            raise ValueError('invalid proof node index')
        index, coefficient = int(parts[1]), parse_coefficient(parts[2])
        if coefficient.is_zero():
            raise ValueError('solver records must omit zero coefficients')
        left = parse_word(parts[3])
        if parts[0] == 'G':
            if parts[4:] != ['0', ''] or left in polynomials[index]:
                raise ValueError('invalid or duplicate polynomial term')
            polynomials[index][left] = coefficient
        else:
            if parts[4] not in ('0', '1') or not re.fullmatch(r'[1-9][0-9]*', parts[5]):
                raise ValueError('invalid derivation reference')
            reference, right = int(parts[5]), parse_word(parts[6])
            if (parts[4] == '0' and reference > inputs
                    or parts[4] == '1' and reference >= index):
                raise ValueError('out-of-range or forward derivation reference')
            derivations[index].append((coefficient, left, parts[4], reference, right))
    if not polynomials or set(polynomials) != set(range(1, len(polynomials)+1)):
        raise ValueError('solver polynomials must have contiguous positive indices')
    if set(derivations) - set(polynomials):
        raise ValueError('derivation has no polynomial')
    return dict(polynomials), dict(derivations)


def word_source(word):
    return '*'.join(NAMES[index-1] for index in word) if word else '1'


def coefficient_source(coefficient):
    """Serialize field arithmetic with explicit Singular number coefficients.

    Singular's bare integer division truncates: `(1/2)*x1` is zero. Construct
    each scalar in the coefficient field before any division takes place.
    """
    def polynomial(terms):
        result = []
        for monomial, value in terms:
            rational = Fraction(value)
            scalar = f'number({rational.numerator})'
            if rational.denominator != 1:
                scalar = f'({scalar}/number({rational.denominator}))'
            powers = [f'{name}^{power}' for name, power in zip(FIELD.parameters, monomial) if power]
            result.append('*'.join([scalar, *powers]))
        return '('+'+'.join(result)+')' if result else 'number(0)'
    return f'({polynomial(coefficient.numerator)}/{polynomial(coefficient.denominator)})'


def polynomial_source(terms):
    return '+'.join(f'{coefficient_source(coefficient)}*{word_source(word)}'
                    for word, coefficient in terms.items()) or '0'


def header(full, *, degree, library=None, library_version=None):
    subset = [index for index, label in enumerate(full.labels) if label != 'T_4']
    library_name = 'freegb.lib' if library is None else str(library)
    if any(character in library_name for character in ('\\', '"', '\n', '\r')):
        raise ValueError('unsupported Singular library path')
    source = (f'LIB "{library_name}";\n'
              'print("SOLVER_VERSION|"+string(system("version")));\n'
              + (f'print("FREEGB_VERSION|{library_version}");\n' if library_version is not None else '') +
              'ring r=(0,q,t),(y3,y2,y1,x3,x2,x1),(c,Dp);\n'
              f'ring R=freeAlgebra(r,{degree});\nshort=0;\nideal I=')
    source += ',\n'.join(polynomial_source(dict(full.relations[index].terms))
                           for index in subset)+';\n'
    return source, subset


def run_singular(binary, source, directory, stage, *, timeout, memory_mib):
    """Run one Linux-local stage with wall-clock and virtual-memory limits."""
    script, log = directory/(stage+'.sing'), directory/(stage+'.log')
    script.write_text(source, encoding='utf-8')

    def limit_memory():
        maximum = memory_mib*1024*1024
        resource.setrlimit(resource.RLIMIT_AS, (maximum, maximum))

    started = time.monotonic()
    with log.open('w', encoding='utf-8') as stream:
        process = subprocess.Popen([str(binary), '-q', str(script.resolve())],
                                   stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT,
                                   start_new_session=True, preexec_fn=limit_memory)
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired as error:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            raise RuntimeError(f'{stage} exceeded {timeout}s; inspect {log}') from error
    if code:
        raise RuntimeError(f'{stage} exited with status {code}; inspect {log}')
    text = log.read_text(encoding='utf-8')
    if any(line.lstrip().startswith('?') for line in text.splitlines()):
        raise RuntimeError(f'{stage} reported a Singular error; inspect {log}')
    return text, time.monotonic()-started


def parse_bridge(proof_text, bridge_text, *, proved_count, inputs):
    """A later stage may append nodes, never modify a previously proved node."""
    for line in bridge_text.splitlines():
        if line.startswith(('G|', 'D|')):
            parts = line.split('|')
            if len(parts) < 2 or not re.fullmatch(r'[1-9][0-9]*', parts[1]):
                raise ValueError('invalid bridge node index')
            if int(parts[1]) <= proved_count:
                raise ValueError('bridge attempted to modify an existing proof node')
    return parse_output(proof_text+'\n'+bridge_text, inputs=inputs)


def solver_environment(binary, output, *, timeout, memory_mib):
    """Resolve and identify the actual freegb source before expensive work."""
    source = ('print("SOLVER_VERSION|"+string(system("version")));\n'
              'print("LIBRARY_PATH|"+string(system("SingularLib")));\nquit;\n')
    text, seconds = run_singular(binary, source, output, 'environment',
                                 timeout=min(timeout, 15), memory_mib=memory_mib)
    versions = [line.split('|', 1)[1] for line in text.splitlines() if line.startswith('SOLVER_VERSION|')]
    paths = [line.split('|', 1)[1] for line in text.splitlines() if line.startswith('LIBRARY_PATH|')]
    if len(versions) != 1 or len(paths) != 1 or not re.fullmatch(r'[0-9]+', versions[0]):
        raise RuntimeError('Singular environment probe did not identify the runtime')
    directories = [Path.cwd(), *(Path(path) for path in paths[0].split(os.pathsep) if path)]
    library = next((directory/'freegb.lib' for directory in directories
                    if (directory/'freegb.lib').is_file()), None)
    if library is None:
        raise RuntimeError('freegb.lib was not found in the Singular library path')
    library = library.resolve()
    match = re.search(r'^version\s*=\s*"([^"\r\n]+)"', library.read_text(encoding='utf-8'), re.MULTILINE)
    if match is None or any(character in match[1] for character in ('\\', '\n', '\r')):
        raise RuntimeError('freegb.lib has no usable version declaration')
    metadata = {'singular_version_id': versions[0], 'freegb_version': match[1].strip(),
                'freegb_source': {'name': library.name, **file_identity(library)}}
    with (output/'environment.log').open('a', encoding='utf-8') as stream:
        stream.write('FREEGB_VERSION|'+metadata['freegb_version']+'\n')
        stream.write('FREEGB_SHA256|'+metadata['freegb_source']['sha256']+'\n')
    return library, metadata, seconds


def make_certificate(full, candidate, proved, derivations, *, subset):
    """Put the compact candidate first; preserve and remap all local proof nodes."""
    count, targets = len(proved)-len(candidate), len(candidate)
    if count <= 0 or set(derivations) != set(proved):
        raise ValueError('every proof node must have a derivation')
    order = list(range(count+1, count+targets+1))+list(range(1, count+1))
    remap = {old: new for new, old in enumerate(order)}
    if any(proved[count+i] != candidate[i] for i in candidate):
        raise ValueError('bridge polynomials disagree with the candidate')
    basis = tuple(full.algebra.element(proved[index]) for index in order)
    rows = tuple(tuple(ProofTerm(c, left, 'input' if kind == '0' else 'basis',
                                subset[reference-1] if kind == '0' else remap[reference], right)
                       for c, left, kind, reference, right in derivations[index]) for index in order)
    return GroebnerCertificate(full.content_hash(), basis, rows)


def square_zero_identity(full):
    """Literal three-term cofactor identity for w², without quotient reduction."""
    algebra = full.algebra
    q, t = (algebra.field.parameter(name) for name in ('q', 't'))
    x, y = (algebra.generator(name) for name in ('x1', 'y1'))
    slider = q*y+(1-q)-x
    shifted = x-t
    w = shifted*slider
    labels = full.labels
    identity = CofactorIdentity(full.content_hash(), w*w, (
        (algebra.one, labels.index('T_2'), slider),
        (q*shifted, labels.index('INV_1_m'), slider),
        (-q*shifted, labels.index('INV_1_p'), slider)))
    if not identity.verify(full):
        raise ArithmeticError('the literal square-zero cofactor identity failed')
    return identity


def framed_proof(full, certificate, candidate_count):
    """Plain data for an independent CAS; all references here are one-based."""
    lines = [f'META|{len(full.relations)}|{candidate_count}|{len(certificate.basis)}']
    def word(value):
        return ''.join(f'{index},' for index in value)
    for kind, polynomials in (('I', full.relations), ('G', certificate.basis)):
        for index, polynomial in enumerate(polynomials, 1):
            lines.extend(f'{kind}|{index}|{coefficient}|{word(monomial)}|0|'
                         for monomial, coefficient in polynomial.terms)
    for index, row in enumerate(certificate.derivations, 1):
        lines.extend(f'D|{index}|{term.coefficient}|{word(term.left)}|'
                     f'{int(term.kind == "basis")}|{term.index+1}|{word(term.right)}' for term in row)
    lines.append(f'END|{len(full.relations)}|{candidate_count}|{len(certificate.basis)}')
    return '\n'.join(lines)+'\n'


def produce(binary, output, *, timeout=180, memory_mib=6144, max_nodes=600, degree=24):
    """Build the full presentation, complete an auxiliary subset, and export proofs."""
    if timeout <= 0 or memory_mib < 256 or max_nodes < 31 or not 12 <= degree <= 48:
        raise ValueError('invalid resource limits')
    if output.exists() and any(output.iterdir()):
        raise ValueError('output directory must be empty; use a fresh run directory')
    prepare_output_root(output)
    producer = capture_provenance()
    checker = Path(__file__).with_name('verify_four_strand_square_zero.g')
    input_paths = (__file__, ENGINE, binary, checker)
    inputs = [{'name': Path(path).name, **file_identity(path)} for path in input_paths]
    full = presentation()
    identity = square_zero_identity(full)
    library, runtime, environment_seconds = solver_environment(
        binary, output, timeout=timeout, memory_mib=memory_mib)
    inputs.append(runtime['freegb_source'])
    (output/'presentation.json').write_text(json.dumps(full.to_json(), indent=2)+'\n')
    common, subset = header(full, degree=degree, library=library, library_version=runtime['freegb_version'])
    engine = ENGINE.read_text(encoding='utf-8')
    engine_prefix = f'int maxNodes={max_nodes};int maxDegree={degree};\n'
    candidate_source, _ = header(full, degree=12, library=library, library_version=runtime['freegb_version'])
    candidate_source += ('option(redSB);option(redTail);\n'+engine_prefix+engine+
                         '\nideal J=twostd(I);\n'
                         'int k;list terms;for(k=1;k<=size(J);k++){B[k]=J[k];emitNode(k,terms);}\n'
                         'print("CANDIDATE_DONE|"+string(size(J)));quit;\n')
    candidate_text, candidate_seconds = run_singular(
        binary, candidate_source, output, 'candidate', timeout=timeout, memory_mib=memory_mib)
    candidate, _ = parse_output(candidate_text)
    if f'CANDIDATE_DONE|{len(candidate)}' not in candidate_text.splitlines():
        raise RuntimeError('candidate completion marker missing')
    candidate_indices = sorted(candidate)
    heads = [max(candidate[index], key=lambda word: (len(word), word)) for index in candidate_indices]
    missing = ','.join('"'+''.join(str(7-index) for index in word)+'"' for word in heads)
    proof_source = (common+'option(noredSB);option(noredTail);\n'+engine_prefix+engine+
                    f'\nMissing=list({missing});\ncompleteToHeads();\nquit;\n')
    proof_text, proof_seconds = run_singular(
        binary, proof_source, output, 'derivation', timeout=timeout, memory_mib=memory_mib)
    proved, derivations = parse_output(proof_text, inputs=len(subset))
    if not any(re.fullmatch(r'FINISH\|[0-9]+\|0\|[0-9]+', line) for line in proof_text.splitlines()):
        raise RuntimeError('derivation did not cover every candidate leading word within its limits')
    bridge_source = common+engine_prefix+engine+'\n'
    for index, terms in sorted(proved.items()):
        bridge_source += f'B[{index}]={polynomial_source(terms)};Words[{index}]=wstr(B[{index}]);\n'
    bridge_source += (f'nodeCount={len(proved)};ideal Targets='+
                      ',\n'.join(polynomial_source(candidate[index]) for index in candidate_indices)+
                      ';\nbridgeTargets();quit;\n')
    bridge_text, bridge_seconds = run_singular(
        binary, bridge_source, output, 'bridge', timeout=timeout, memory_mib=memory_mib)
    # Contiguity is checked over the complete stream; the bridge extends the DAG.
    all_proved, all_derivations = parse_bridge(
        proof_text, bridge_text, proved_count=len(proved), inputs=len(subset))
    if f'BRIDGE_DONE|{len(candidate)}' not in bridge_text.splitlines():
        raise RuntimeError('candidate ideal-membership bridge failed')
    certificate = make_certificate(full, candidate, all_proved, all_derivations, subset=subset)
    (output/'certificate.json').write_text(json.dumps(certificate.to_json(), separators=(',', ':'))+'\n')
    (output/'proof.txt').write_text(framed_proof(full, certificate, len(candidate)))
    (output/'w-square-identity.json').write_text(json.dumps(identity.to_json(), indent=2)+'\n')
    report = {'schema': 'zipper-four-strand-square-zero-run-v1',
              'presentation_hash': full.content_hash(), 'field': FIELD.to_json(), 'runtime': runtime,
              'input_relations': len(full.relations), 'auxiliary_input_labels': [full.labels[i] for i in subset],
              'candidate_polynomials': len(candidate), 'proof_nodes': len(certificate.basis),
              'local_proof_terms': sum(map(len, certificate.derivations)),
              'verification': 'candidate produced; independent verification required',
              'limits': {'timeout_seconds_per_stage': timeout, 'memory_mib': memory_mib,
                         'maximum_proof_nodes': max_nodes, 'maximum_word_degree': degree},
              'seconds': {'environment': environment_seconds, 'candidate': candidate_seconds,
                          'derivation': proof_seconds, 'bridge': bridge_seconds}}
    (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    write_bundle_provenance(output, sorted(path.name for path in output.iterdir() if path.is_file()),
                            producer=producer, inputs=inputs)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--singular', default='Singular', help='path or command for Singular 4.4')
    parser.add_argument('--output', type=Path, default=Path('runs/four-strand-square-zero'))
    parser.add_argument('--timeout', type=int, default=180, help='wall-clock seconds per solver stage')
    parser.add_argument('--memory-mib', type=int, default=6144, help='solver virtual-memory cap')
    parser.add_argument('--max-nodes', type=int, default=600)
    parser.add_argument('--degree', type=int, default=24)
    args = parser.parse_args(argv)
    executable = shutil.which(args.singular)
    if executable is None:
        parser.error(f'Singular executable not found: {args.singular}')
    try:
        report = produce(Path(executable).resolve(), args.output, timeout=args.timeout,
                         memory_mib=args.memory_mib, max_nodes=args.max_nodes, degree=args.degree)
    except (ValueError, RuntimeError, ArithmeticError, OSError, UnicodeError) as error:
        parser.exit(1, f'{error}\n')
    print(f"Produced {report['candidate_polynomials']} candidate polynomials and "
          f"{report['proof_nodes']} proof nodes in {args.output}; run the independent checker.")
    return report


if __name__ == '__main__':
    main()
