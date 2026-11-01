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
from zipper_algebra.provenance import capture_provenance, file_identity, write_bundle_provenance
from zipper_algebra.verification.derivations import VerificationLimit
from zipper_algebra.verification.finite import WordBasisCertificate


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--presentation', type=Path, required=True)
    parser.add_argument('--word-basis', type=Path, required=True)
    parser.add_argument('--target', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-operations', type=int, default=10_000_000)
    parser.add_argument('--max-contexts', type=int, default=100_000)
    parser.add_argument('--max-nodes', type=int, default=100_000)
    parser.add_argument('--max-dimension', type=int, default=10_000)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    prepare_output_root(output)
    if any(path.name != '.gitignore' for path in output.iterdir()):
        parser.error('output must be a new or empty directory')
    producer = capture_provenance()
    driver_identity = file_identity(Path(__file__))
    inputs = []
    report = {'schema': 'zipper-cofactor-reduction-run-v1', 'status': 'inconclusive',
              'verified': False, 'producer': producer,
              'limits': {'max_operations': args.max_operations, 'max_contexts': args.max_contexts,
                         'max_nodes': args.max_nodes, 'max_dimension': args.max_dimension},
              'computation_source': {'name': Path(__file__).name, **driver_identity}}
    documents = []
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
        result = cofactor_representation(presentation, certificate, target,
                    max_operations=args.max_operations, max_contexts=args.max_contexts,
                    max_nodes=args.max_nodes, max_dimension=args.max_dimension)
        write_json_atomic(output / 'representative.json', result.remainder.to_json())
        write_json_atomic(output / 'cofactor.json', result.identity.to_json())
        # Check the actual serialized identity against the original input.
        saved = CofactorIdentity.from_json(json.loads((output / 'cofactor.json').read_text()))
        if not verify_cofactor(presentation, saved, max_operations=args.max_operations):
            raise ValueError('saved cofactor fails its exact original-input identity')
        validate_source_snapshot(producer)
        if file_identity(Path(__file__)) != driver_identity:
            raise ValueError('computation source changed during cofactor production')
        documents = ['representative.json', 'cofactor.json']
        report.update(status='verified', verified=True,
                      presentation_hash=presentation.content_hash(),
                      in_ideal=result.remainder.is_zero(),
                      cofactor_terms=len(result.identity.terms), statistics=result.statistics,
                      artifacts={name: str(output / name) for name in documents})
    except (OSError, ValueError, TypeError, ImportError, ArithmeticError, VerificationLimit) as error:
        report.update(reason=str(error))
    report['inputs'] = inputs
    write_json_atomic(output / 'report.json', report)
    write_bundle_provenance(output, documents + ['report.json'], producer=producer,
                            inputs=inputs + [{'name': Path(__file__).name, **driver_identity}])
    print(json.dumps(report, sort_keys=True))
    return 0 if report['verified'] else 3


if __name__ == '__main__':
    raise SystemExit(main())
