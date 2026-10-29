"""Run the current unit and integration tests."""
from __future__ import annotations

import argparse
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def cases(suite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from cases(item)
        else:
            yield item


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list discovered test identities")
    args = parser.parse_args(argv)
    loader = unittest.TestLoader()
    suite = loader.discover(str(ROOT / "tests"), top_level_dir=str(ROOT))
    if loader.errors:
        raise RuntimeError("Test discovery failed:\n" + "\n".join(loader.errors))
    if not suite.countTestCases():
        raise RuntimeError("No tests discovered")
    if args.list:
        for case in cases(suite):
            print(case.id())
        return 0
    return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
