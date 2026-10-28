"""Enforce separate line and branch coverage targets for the maintained package."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def check(report, root=ROOT, minimum=99.5):
    root = Path(root).resolve()
    source = root / 'src/zipper_algebra'
    expected = {path.resolve() for path in source.rglob('*.py')}
    files = report['files']
    measured = {Path(name).resolve() if Path(name).is_absolute()
                else (root / name).resolve(): data for name, data in files.items()}
    if set(measured) != expected:
        raise ValueError('coverage must report every maintained Python file, and only those files')
    if any(data.get('excluded_lines') for data in measured.values()):
        raise ValueError('maintained source must have no excluded lines')
    sums = {key: sum(data['summary'][key] for data in measured.values())
            for key in ('num_statements', 'covered_lines', 'num_branches', 'covered_branches')}
    if not sums['num_statements'] or not sums['num_branches']:
        raise ValueError('line and branch measurements are both required')
    lines = 100 * sums['covered_lines'] / sums['num_statements']
    branches = 100 * sums['covered_branches'] / sums['num_branches']
    print(f"Maintained package: {lines:.2f}% lines "
          f"({sums['covered_lines']}/{sums['num_statements']}), "
          f"{branches:.2f}% branches ({sums['covered_branches']}/{sums['num_branches']})")
    if lines < minimum or branches < minimum:
        raise ValueError(f'both line and branch coverage must be at least {minimum}% before rounding')
    return sums


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', nargs='?', default='coverage.json')
    args = parser.parse_args()
    try:
        check(json.loads(Path(args.report).read_text(encoding='utf-8')))
    except (KeyError, OSError, TypeError, ValueError) as error:
        raise SystemExit(f'Coverage check failed: {error}')
