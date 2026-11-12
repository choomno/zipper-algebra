"""Produce an explicit original-input cofactor from a finite-algebra certificate."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from zipper_algebra.algebra import AlgebraElement, AlgebraPresentation, CofactorIdentity
from zipper_algebra.cofactor.artifacts import write_json_atomic
from zipper_algebra.cofactor.general import cofactor_representation, verify_cofactor
from zipper_algebra.execution.artifacts import prepare_output_root, validate_source_snapshot
from zipper_algebra.execution.dependencies import reconstruction_provenance
from zipper_algebra.provenance import capture_provenance, file_identity, write_bundle_provenance
from zipper_algebra.verification.finite import WordBasisCertificate


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--presentation', type=Path, required=True)
    parser.add_argument('--word-basis', type=Path, required=True)
    parser.add_argument('--target', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--backend', choices=('literal', 'finiteflow'), default='literal')
    parser.add_argument('--max-operations', type=int, default=10_000_000)
    parser.add_argument('--max-contexts', type=int, default=100_000)
    parser.add_argument('--max-nodes', type=int, default=100_000)
    parser.add_argument('--max-dimension', type=int, default=10_000)
    parser.add_argument('--max-matrix-entries', type=int, default=100_000,
                        help='Python sparse matrix entry bound for the finiteflow backend')
    parser.add_argument('--max-primes', type=int, default=5, help='FiniteFlow prime budget')
    parser.add_argument('--max-degree', type=int, default=24, help='FiniteFlow degree search bound')
    parser.add_argument('--threads', type=int, default=1, help='FiniteFlow reconstruction workers')
    parser.add_argument('--start-prime', type=int, default=0, help='FiniteFlow starting prime index')
    parser.add_argument('--shift', type=int, nargs='*', help='FiniteFlow integer parameter shift')
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or
                            any(path.name != '.gitignore' for path in output.iterdir())):
        parser.error('output must be a new or empty directory')
    prepare_output_root(output)
    producer = capture_provenance()
    driver_identity = file_identity(Path(__file__))
    inputs = []
    report = {'schema': 'zipper-cofactor-reduction-run-v1', 'status': 'inconclusive',
              'verified': False, 'producer': producer, 'backend': args.backend,
              'limits': {'max_operations': args.max_operations, 'max_contexts': args.max_contexts,
                         'max_nodes': args.max_nodes, 'max_dimension': args.max_dimension},
              'computation_source': {'name': Path(__file__).name, **driver_identity}}
    documents = []
    native_options = {'max_matrix_entries': args.max_matrix_entries, 'max_primes': args.max_primes,
                      'max_degree': args.max_degree, 'threads': args.threads,
                      'start_prime': args.start_prime, 'shift': args.shift}
    try:
        decoded = []
        for name, path in (('presentation.json', args.presentation),
                           ('word-basis-certificate.json', args.word_basis), ('target.json', args.target)):
            content = path.read_bytes()
            inputs.append({'name': name, 'bytes': len(content),
                           'sha256': hashlib.sha256(content).hexdigest()})
            decoded.append(json.loads(content))
        presentation = AlgebraPresentation.from_json(decoded[0])
        certificate = WordBasisCertificate.from_json(decoded[1])
        target = AlgebraElement.from_json(decoded[2])
        del decoded
        dependencies = reconstruction_provenance() if args.backend == 'finiteflow' else None
        if dependencies is not None:
            report['reconstruction'] = {'dependencies': dependencies, 'options': native_options}
        result = cofactor_representation(presentation, certificate, target,
                    backend=args.backend, **native_options,
                    max_operations=args.max_operations, max_contexts=args.max_contexts,
                    max_nodes=args.max_nodes, max_dimension=args.max_dimension)
        del certificate
        expected_representative, statistics = result.remainder, result.statistics
        write_json_atomic(output / 'presentation.json', presentation.to_json())
        write_json_atomic(output / 'target.json', target.to_json())
        write_json_atomic(output / 'representative.json', result.remainder.to_json())
        write_json_atomic(output / 'cofactor.json', result.identity.to_json())
        del result
        # Bind the saved problem and normal form before checking the explicit proof.
        saved_presentation = AlgebraPresentation.from_json(json.loads((output / 'presentation.json').read_text()))
        saved_target = AlgebraElement.from_json(json.loads((output / 'target.json').read_text()))
        if saved_presentation.content_hash() != presentation.content_hash() or saved_target != target:
            raise ValueError('saved presentation or target differs from the consumed problem')
        representative = AlgebraElement.from_json(json.loads((output / 'representative.json').read_text()))
        saved = CofactorIdentity.from_json(json.loads((output / 'cofactor.json').read_text()))
        if representative != expected_representative:
            raise ValueError('saved representative differs from the verified normal form')
        if saved.target != target-representative:
            raise ValueError('saved cofactor target differs from target minus representative')
        if not verify_cofactor(saved_presentation, saved, max_operations=args.max_operations):
            raise ValueError('saved cofactor fails its exact original-input identity')
        if dependencies is not None and reconstruction_provenance() != dependencies:
            raise ValueError('reconstruction dependencies changed during cofactor production')
        validate_source_snapshot(producer)
        if file_identity(Path(__file__)) != driver_identity:
            raise ValueError('computation source changed during cofactor production')
        documents = ['presentation.json', 'target.json', 'representative.json', 'cofactor.json']
        report.update(status='verified', verified=True,
                      presentation_hash=presentation.content_hash(),
                      in_ideal=representative.is_zero(),
                      cofactor_terms=len(saved.terms), statistics=statistics,
                      artifacts={name: str(output / name) for name in documents})
    except (OSError, ValueError, TypeError, ImportError, ArithmeticError, LookupError, RuntimeError) as error:
        for name in ('presentation.json', 'target.json', 'representative.json', 'cofactor.json'):
            (output / name).unlink(missing_ok=True)
        report.update(reason=str(error))
    report['inputs'] = inputs
    write_json_atomic(output / 'report.json', report)
    write_bundle_provenance(output, documents + ['report.json'], producer=producer,
                            inputs=inputs + [{'name': Path(__file__).name, **driver_identity}])
    print(json.dumps(report, sort_keys=True))
    return 0 if report['verified'] else 3


if __name__ == '__main__':
    raise SystemExit(main())
