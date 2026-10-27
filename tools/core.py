#!/usr/bin/env python3
"""Compatibility entry point for :mod:`zipper_algebra.core`.

The maintained implementation is installed from ``src/zipper_algebra``.
"""
from importlib import import_module
import runpy
import sys

if __name__ == "__main__":
    runpy.run_module("zipper_algebra.core", run_name="__main__")
else:
    sys.modules[__name__] = import_module("zipper_algebra.core")
