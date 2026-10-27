"""Exact algebraic infrastructure for zipper relation research.

The reference documentation specifies coefficient domains, ambient algebras and
the conclusions established by the individual certificate checks.
"""

from importlib import import_module

__version__ = "0.1.0.dev0"
__all__ = ["Atom", "AtomSet", "Cofactor", "Presentation", "__version__"]


def __getattr__(name):
    modules = {"Presentation": ".presentation", "Atom": ".cofactor.atoms",
               "AtomSet": ".cofactor.atoms", "Cofactor": ".cofactor.cofactor_io"}
    if name not in modules:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(modules[name], __name__), name)
    globals()[name] = value
    return value
