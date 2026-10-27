"""Command line presentation descriptions and exact certificate verification."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .presentation import Presentation, EXTRA_CATALOGUE
from .cofactor.cofactor_io import Cofactor
from .cofactor.verify import check, pit_check


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
    return command


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        presentation = _presentation(args)
        if args.command == "describe":
            report = {"presentation": presentation.slug(),
                      "presentation_hash": presentation.content_hash(),
                      "n": presentation.n, "variables": presentation.var_names(),
                      "relations": presentation.labels()}
            if args.json:
                print(json.dumps(report, sort_keys=True))
            else:
                print(report["presentation"], report["presentation_hash"],
                      len(report["relations"]), "relations")
                for index, label in enumerate(report["relations"]):
                    print(f"  R_{index:<3} {label}")
            return 0

        certificate = Cofactor.load(str(args.certificate), pres=presentation, target=args.target)
        target = args.target if args.target is not None else certificate.target
        if not target:
            raise ValueError("this certificate has no target; supply --target")
        exact = check(certificate, presentation, target, verbose=not args.json)
        report = {"pass": exact["pass"], "exact": exact}
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
