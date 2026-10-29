"""Installed commands for presentations, exact checking and solver-backed discovery."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

from .presentation import Presentation, EXTRA_CATALOGUE
from .cofactor.cofactor_io import Cofactor
from .cofactor.verify import check, pit_check
from .provenance import capture_provenance, write_output_provenance


def _integer_at_least(minimum):
    def parse(value):
        try:
            integer = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError("expected an integer") from None
        if integer < minimum:
            raise argparse.ArgumentTypeError(f"expected an integer at least {minimum}")
        return integer
    return parse


def _presentation_arguments(parser, *, require_n):
    parser.add_argument("--n", type=_integer_at_least(2), required=require_n,
                        default=None if require_n else 4, help="number of strands")
    parser.add_argument("--untwist", default="bigelow",
                        choices=("bigelow", "square-zero", "silly_question", "none"))
    parser.add_argument("--framing-power", type=_integer_at_least(1), default=1)
    parser.add_argument("--theta-power", type=_integer_at_least(1), default=1)
    parser.add_argument("--max-zipper-level", type=_integer_at_least(2))
    parser.add_argument("--max-untwist-level", type=_integer_at_least(2))
    parser.add_argument("--zero-from", type=_integer_at_least(2))
    parser.add_argument("--extra", action="append", default=[], choices=tuple(EXTRA_CATALOGUE))
    parser.add_argument("--json", action="store_true", help="write a machine-readable report")


def _presentation(args):
    for name in ("max_zipper_level", "max_untwist_level", "zero_from"):
        value = getattr(args, name)
        if value is not None and value > args.n:
            raise ValueError(f"--{name.replace('_', '-')} cannot exceed --n")
    return Presentation(args.n, args.untwist, args.framing_power, args.theta_power,
                        args.max_zipper_level, args.max_untwist_level, args.zero_from,
                        tuple(args.extra))


def _parameter_pair(value):
    try:
        q, t = value.split(",")
        return int(q), int(t)
    except ValueError:
        raise argparse.ArgumentTypeError("expected an integer pair q,t") from None


def parser():
    command = argparse.ArgumentParser(
        prog="zipper-algebra", description="Exact tools for zipper algebra presentations")
    command.add_argument("--version", action="version", version="%(prog)s 0.1.0.dev0")
    subcommands = command.add_subparsers(dest="command", required=True)
    describe = subcommands.add_parser("describe", help="describe an algebra presentation")
    _presentation_arguments(describe, require_n=False)
    verify = subcommands.add_parser("verify", help="verify a free-algebra membership certificate")
    verify.add_argument("certificate", type=Path)
    _presentation_arguments(verify, require_n=True)
    verify.add_argument("--target", help="target slider or additional relation name")
    verify.add_argument("--pit", action="store_true", help="also check the exact finite product grid")
    emit = subcommands.add_parser("emit", help="emit only the specified defining relations")
    _presentation_arguments(emit, require_n=True)
    emit.add_argument("--q", type=int, default=5)
    emit.add_argument("--t", type=int, default=7)
    emit.add_argument("--char", type=_integer_at_least(0), default=0)
    emit.add_argument("--out", type=Path, required=True)
    compose = subcommands.add_parser("compose", help="compose and exactly verify symbolic proof data")
    from .cofactor.compose import add_arguments
    add_arguments(compose)
    discover = subcommands.add_parser("discover", help="produce a certificate from solver proofs")
    _presentation_arguments(discover, require_n=True)
    discover.add_argument("--target", required=True)
    discover.add_argument("--solver", type=Path, required=True)
    discover.add_argument("--output", type=Path, required=True)
    discover.add_argument("--cache", type=Path)
    discover.add_argument("--anchor", type=_parameter_pair, default=(5, 7), metavar="Q,T")
    discover.add_argument("--start", type=_integer_at_least(1), default=100)
    discover.add_argument("--step", type=_integer_at_least(1), default=40)
    discover.add_argument("--max", dest="maxpts", type=_integer_at_least(1), default=260)
    discover.add_argument("--jobs", type=_integer_at_least(1), default=1)
    discover.add_argument("--atoms", type=Path)
    discover.add_argument("--timeout", type=float, default=1800)
    discover.add_argument("--max-iterations", type=_integer_at_least(1), default=100000)
    discover.add_argument("--max-degree", type=_integer_at_least(1))
    discover.add_argument("--solver-revision", default="unrecorded")
    discover.add_argument("--solver-io", choices=("files", "sqlite"), default="files",
                          help="native solver I/O route; sqlite requires a SQLite-enabled f4ncgb build")
    discover.add_argument("--refresh", action="store_true")
    export = subcommands.add_parser("export-sample", help="export one stored sample for inspection")
    export.add_argument("--store", type=Path, required=True)
    export.add_argument("--key", required=True)
    export.add_argument("--generation")
    export.add_argument("--output", type=Path, required=True)
    export.add_argument("--json", action="store_true")
    return command


def main(argv=None):
    args = parser().parse_args(argv)
    if args.command == "compose":
        from .cofactor.compose import run as run_composition
        return run_composition(args)
    try:
        if args.command == "export-sample":
            from .execution.store import export
            producer = capture_provenance()
            record = export(args.store, args.key, args.output, generation=args.generation)
            report = {"store": str(args.store), "key": args.key,
                      "generation": record["generation"], "output": str(args.output),
                      "original_producer": record["request"]["producer"], "producer": producer}
            print(json.dumps(report, sort_keys=True) if args.json
                  else f"Exported sample {args.key} to {args.output}")
            return 0
        presentation = _presentation(args)
        if args.command == "describe":
            report = {"presentation": presentation.slug(),
                      "presentation_hash": presentation.content_hash(),
                      "n": presentation.n, "variables": presentation.var_names(),
                      "relations": presentation.labels()}
            if args.json:
                report["producer"] = capture_provenance()
                print(json.dumps(report, sort_keys=True))
            else:
                print(report["presentation"], report["presentation_hash"],
                      len(report["relations"]), "relations")
                for index, label in enumerate(report["relations"]):
                    print(f"  R_{index:<3} {label}")
            return 0

        if args.command == "emit":
            producer = capture_provenance()
            args.out.parent.mkdir(parents=True, exist_ok=True)
            count = presentation.emit_ms(str(args.out), args.q, args.t, args.char, provenance=False)
            write_output_provenance(args.out, producer=producer)
            report = {"input": str(args.out), "relations": count, "producer": producer,
                      "presentation_hash": presentation.content_hash()}
            print(json.dumps(report, sort_keys=True) if args.json
                  else f"Wrote {count} defining relations to {args.out}")
            return 0
        if args.command == "discover":
            from .cofactor.atoms import AtomSet
            from .execution.config import RunConfig
            from .execution.pipeline import run as run_discovery
            config = RunConfig(args.solver, args.output, args.cache, args.timeout,
                               args.max_iterations, args.max_degree,
                               solver_revision=args.solver_revision, solver_io=args.solver_io)
            atoms = None if args.atoms is None else AtomSet.load(args.atoms)
            report = run_discovery(presentation, args.target, (*args.anchor, 0), config,
                         start=args.start, step=args.step, maxpts=args.maxpts,
                         jobs=args.jobs, atoms=atoms, refresh=args.refresh, verbose=not args.json)
            print(json.dumps(report, sort_keys=True) if args.json else report["status"])
            return 0 if report["verified"] else 3

        producer = capture_provenance()
        content = args.certificate.read_bytes()
        inputs = [{"name": args.certificate.name, "bytes": len(content),
                   "sha256": hashlib.sha256(content).hexdigest()}]
        certificate = Cofactor.from_text(content.decode('utf-8'), pres=presentation,
                                         target=args.target, name=str(args.certificate))
        target = args.target if args.target is not None else certificate.target
        if not target:
            raise ValueError("this certificate has no target; supply --target")
        exact = check(certificate, presentation, target, verbose=not args.json)
        report = {"pass": exact["pass"], "exact": exact, "producer": producer, "inputs": inputs}
        if args.pit:
            grid = pit_check(certificate, presentation, target, verbose=not args.json)
            report["pit"] = grid
            report["pass"] = report["pass"] and grid["pass"]
        if args.json:
            print(json.dumps(report, sort_keys=True))
        return 0 if report["pass"] else 1
    except (OSError, ValueError, TypeError, KeyError, ZeroDivisionError) as error:
        if args.json:
            print(json.dumps({"pass": False, "error": str(error), "status": "invalid-input"},
                             sort_keys=True))
        else:
            print(f"zipper-algebra: {error}", file=sys.stderr)
        return 2
