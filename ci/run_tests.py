"""Run the complete suite or its complementary platform/reproduction selections."""
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


def select(name):
    loader = unittest.TestLoader()
    discovered = loader.discover(str(ROOT / "tests"), top_level_dir=str(ROOT))
    if loader.errors:
        raise RuntimeError("Test discovery failed:\n" + "\n".join(loader.errors))
    selected = []
    for case in cases(discovered):
        historical = case.id().startswith("tests.reproduction.")
        if name == "all" or (name == "reproduction") == historical:
            selected.append(case)
    if not selected:
        raise RuntimeError(f"The {name} selection contains no tests")
    return unittest.TestSuite(selected)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=("all", "fast", "reproduction"), default="all")
    parser.add_argument("--list", action="store_true", help="list selected case identities")
    args = parser.parse_args(argv)
    suite = select(args.suite)
    if args.list:
        for case in cases(suite):
            print(case.id())
        return 0
    return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
