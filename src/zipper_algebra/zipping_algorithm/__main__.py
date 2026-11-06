"""Bounded zipper computations and standalone JSON certificate checking."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
from typing import Any, Sequence

from .catalog import Catalog
from .algorithm import AlgorithmSession, check_certificate as check_algorithm
from .enumeration import exhaust
from .general_enumeration import exhaust_general
from .interpreter import Limits, Request, RunState, check_certificate, execute
from .transport import (
    BandTransport,
    TransportError,
    TransportStep,
    check_transport,
    transport_band,
)
from .two_strand import check_two_strand, reduce_two_strand
from .verification import (
    LocalVerificationLimit,
    VerificationLimits,
    certify_rule,
    check_local_certificate,
)
from .upward import check_certificate as check_upward
from .upward import initialize, run as run_upward


def _object(value: Any) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError("A JSON object is required")
    return value


def _read(path: str) -> dict[str, Any]:
    return _object(json.loads(Path(path).read_text(encoding="utf-8")))


def _check(certificate: dict[str, Any]) -> dict[str, Any]:
    schema = certificate.get("schema")
    if schema == "zipper-collected-run-v1":
        from .collected import check_certificate as check_collected

        return {"pass": check_collected(certificate), "scope": certificate["scope"]}
    if schema == "zipper-algorithm-run-v1":
        return {"pass": check_algorithm(certificate), "scope": certificate["scope"]}
    if schema == "zipper-contextual-replay-v1":
        return {"pass": check_certificate(certificate), "scope": certificate["scope"]}
    if schema == "zipper-normal-upward-v1":
        return {"pass": check_upward(certificate), "scope": certificate["scope"]}
    if schema == "zipper-local-original-relator-v1":
        return {
            "pass": check_local_certificate(certificate),
            "scope": certificate["scope"],
        }
    if schema == "zipper-two-strand-run-v2":
        return check_two_strand(certificate)
    if schema == "zipper-band-transport-v1":
        if set(certificate) != {"schema", "witness"}:
            raise ValueError("Unexpected transport certificate fields")
        fields = dict(_object(certificate["witness"]))
        fields["word"] = tuple(fields["word"])
        fields["steps"] = tuple(
            TransportStep(**dict(_object(step), block=tuple(step["block"])))
            for step in fields["steps"]
        )
        check_transport(BandTransport(**fields))
        return {"pass": True, "scope": "whole_band_transport"}
    raise ValueError("Unsupported certificate schema")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    instantiate = commands.add_parser(
        "instantiate", help="Instantiate one exact local identity"
    )
    instantiate.add_argument("rule")
    instantiate.add_argument(
        "--parameters", required=True, help="JSON original upward parameter bindings"
    )
    instantiate.add_argument("--words", default="{}", help="JSON named braid words")
    instantiate.add_argument("--expand-sliders", action="store_true")
    instantiate.add_argument("--max-terms", type=int, default=100000)
    verify = commands.add_parser(
        "verify-rule", help="Certify a local identity against original zipper relators"
    )
    verify.add_argument("rule")
    verify.add_argument(
        "--parameters", required=True, help="JSON original upward parameter bindings"
    )
    verify.add_argument("--words", default="{}", help="JSON named braid words")
    verify.add_argument("--limits", default="{}", help="JSON exact verification limits")
    run = commands.add_parser("run", help="Run the literal-word zipper procedure")
    run.add_argument("--n", type=int, required=True)
    run.add_argument("--word", required=True, help="JSON literal signed braid word")
    run.add_argument("--limits", default="{}", help="JSON exact computation limits")
    run.add_argument("--trace", action="store_true")
    collected = commands.add_parser(
        "run-collected", help="Collect exact formal sums and rewrite terms in parallel"
    )
    collected.add_argument("--n", type=int, required=True)
    collected.add_argument("--word", required=True, help="JSON literal signed braid word")
    collected.add_argument("--limits", default="{}", help="JSON exact computation limits")
    collected.add_argument("--trace", action="store_true")
    collected.add_argument("--workers", type=int, default=1)
    collected.add_argument("--batch-size", type=int, default=64)
    general = commands.add_parser(
        "exhaust-general",
        help="Enumerate literal words with the general zipper procedure",
    )
    general.add_argument("output", help="Directory for SQLite checkpoints and summary")
    general.add_argument("--ranks", default="[3]", help="JSON ordered list of ranks")
    general.add_argument("--seconds", type=float, default=60)
    general.add_argument("--max-length", type=int)
    general.add_argument(
        "--limits",
        default='{"max_steps":5000,"max_terms":1024}',
        help="JSON exact per-word limits",
    )
    general.add_argument("--trace", action="store_true")
    enumeration = commands.add_parser(
        "exhaust", help="Enumerate and exactly check literal braid words"
    )
    enumeration.add_argument(
        "output", help="Directory for SQLite checkpoints and summary"
    )
    enumeration.add_argument(
        "--ranks", default="[1,2]", help="JSON ordered list of ranks"
    )
    enumeration.add_argument("--seconds", type=float, default=60)
    enumeration.add_argument("--max-length", type=int)
    enumeration.add_argument("--max-steps", type=int, default=5000)
    enumeration.add_argument("--max-terms", type=int, default=1024)
    enumeration.add_argument("--trace", action="store_true")
    replay = commands.add_parser(
        "replay", help="Execute an explicit contextual schedule"
    )
    replay.add_argument(
        "input", help="JSON object with state, requests, optional limits"
    )
    check = commands.add_parser("check", help="Verify a saved JSON certificate")
    check.add_argument("certificate")
    upper = commands.add_parser("upper", help="Run a normal upward pass")
    upper.add_argument("--left", required=True, help="JSON signed braid word")
    upper.add_argument("--friends", type=int, default=0)
    transport = commands.add_parser(
        "transport", help="Transport a slider through whole-band crossings"
    )
    transport.add_argument("--word", required=True, help="JSON signed braid word")
    transport.add_argument("--mirrored", action="store_true")
    for command in (upper, transport):
        command.add_argument("--n", type=int, required=True)
        command.add_argument("--size", type=int, required=True)
        command.add_argument("--start", type=int, required=True)
    two = commands.add_parser(
        "two-strand", help="Run the complete two-strand procedure"
    )
    two.add_argument("--word", required=True, help="JSON signed braid word")
    for command in (upper, two):
        command.add_argument("--max-steps", type=int, default=100000)
        command.add_argument("--max-terms", type=int, default=100000)
        command.add_argument(
            "--trace",
            action="store_true",
            help="Include every transition in the JSON output",
        )
    return parser


def _execute(args: argparse.Namespace) -> tuple[dict[str, Any], int]:
    if args.command == "run-collected":
        from .collected import reduce_word as reduce_collected

        result = reduce_collected(
            json.loads(args.word), n=args.n,
            limits=Limits(**_object(json.loads(args.limits))), trace=args.trace,
            workers=args.workers, batch_size=args.batch_size,
        )
        return result, 0 if result["status"] == "complete" else 3
    if args.command == "instantiate":
        return Catalog.bundled().instantiate(
            args.rule,
            _object(json.loads(args.parameters)),
            words=_object(json.loads(args.words)),
            expand_sliders=args.expand_sliders,
            max_terms=args.max_terms,
        ), 0
    if args.command == "run":
        result = AlgorithmSession().reduce(
            json.loads(args.word),
            n=args.n,
            limits=Limits(**_object(json.loads(args.limits))),
            trace=args.trace,
        )
        return result, 0 if result["status"] == "complete" else 3
    if args.command == "exhaust-general":
        result = exhaust_general(
            args.output,
            ranks=json.loads(args.ranks),
            wall_seconds=args.seconds,
            max_length=args.max_length,
            limits=Limits(**_object(json.loads(args.limits))),
            trace=args.trace,
        )
        return result, 0 if result["stop_reason"] in (
            "max_length",
            "all_nonempty_layers_exhausted",
        ) else 3
    if args.command == "verify-rule":
        return certify_rule(
            args.rule,
            _object(json.loads(args.parameters)),
            words=_object(json.loads(args.words)),
            limits=VerificationLimits(**_object(json.loads(args.limits))),
        ), 0
    if args.command == "exhaust":
        result = exhaust(
            args.output,
            ranks=json.loads(args.ranks),
            wall_seconds=args.seconds,
            max_length=args.max_length,
            max_steps=args.max_steps,
            max_terms=args.max_terms,
            trace=args.trace,
        )
        complete = (
            result["stop_reason"] in ("max_length", "all_nonempty_layers_exhausted")
            and not result["unsupported_ranks"]
        )
        return result, 0 if complete else 3
    if args.command == "check":
        return _check(_read(args.certificate)), 0
    if args.command == "replay":
        specification = _read(args.input)
        if not {"state", "requests"} <= set(specification) or set(specification) - {
            "state",
            "requests",
            "limits",
        }:
            raise ValueError(
                "Replay input requires state, requests, and optional limits"
            )
        result = execute(
            RunState.from_json(specification["state"]),
            [Request(**_object(x)) for x in specification["requests"]],
            limits=Limits(**_object(specification.get("limits", {}))),
        )
        return result, 3 if result["status"] == "incomplete" else 0
    if args.command == "upper":
        result = run_upward(
            initialize(
                args.n,
                json.loads(args.left),
                size=args.size,
                start=args.start,
                friends=args.friends,
            ),
            limits=Limits(max_steps=args.max_steps, max_terms=args.max_terms),
            trace=args.trace,
        )
        return result, 0 if result["status"] == "upper_pass_complete" else 3
    if args.command == "two-strand":
        result = reduce_two_strand(
            json.loads(args.word),
            max_steps=args.max_steps,
            max_terms=args.max_terms,
            trace=args.trace,
        )
        return result, 0 if result["status"] == "complete" else 3
    witness = transport_band(
        json.loads(args.word),
        n=args.n,
        size=args.size,
        start=args.start,
        mirrored=args.mirrored,
    )
    check_transport(witness)
    return {"schema": "zipper-band-transport-v1", "witness": asdict(witness)}, 0


def main(argv: Sequence[str] | None = None) -> int:
    """Print exact JSON results; use exits 2 for invalid and 3 for incomplete input."""
    args = _parser().parse_args(argv)
    try:
        result, status = _execute(args)
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    except LocalVerificationLimit as exc:
        print(
            json.dumps({"status": "incomplete", "diagnostic": str(exc)}),
            file=sys.stderr,
        )
        return 3
    except TransportError as exc:
        print(
            json.dumps({"status": "unsupported", "diagnostic": str(exc)}),
            file=sys.stderr,
        )
        return 3
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as exc:
        print(
            json.dumps({"status": "invalid", "diagnostic": str(exc)}), file=sys.stderr
        )
        return 2
    return status


if __name__ == "__main__":
    raise SystemExit(main())
