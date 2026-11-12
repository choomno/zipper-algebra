"""Produce a certified generic three-strand quotient from fresh f4ncgb proofs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from zipper_algebra.algebra import AlgebraPresentation
from zipper_algebra.execution.config import RunConfig
from zipper_algebra.execution.production import produce_algebra
from zipper_algebra.fields import RationalFunctionField


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--degree', type=int, choices=(1, 2), default=1)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--solver-io', choices=('files', 'sqlite'), default='files')
    parser.add_argument('--timeout', type=float, default=60)
    args = parser.parse_args(argv)
    coefficient_names = tuple(f'c{i}' for i in range(args.degree))
    names = ('q',) + coefficient_names + (('r',) if args.degree == 2 else ())
    field = RationalFunctionField(names)
    coefficients = tuple(field.parameter(name) for name in coefficient_names) + (field.one,)
    untwistings = {2: coefficients}
    if args.degree == 2:
        untwistings[3] = (-field.parameter('r'), field.one)
    presentation = AlgebraPresentation.zipper(3, field=field, q=field.parameter('q'),
                                              untwistings=untwistings)
    # A rational anchor supplies a proof proposal. The returned statement and
    # all certificate coefficients concern Q(q,c0,...) itself.
    point = dict(zip(names, (5, 7, 11, 13)))
    config = RunConfig(args.solver, args.output, timeout_seconds=args.timeout,
                       solver_io=args.solver_io)
    report = produce_algebra(presentation, point, config, threads=args.threads)
    print(json.dumps(report, sort_keys=True))
    return 0 if report['verified'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
