"""Exact zipper rewrite data, contextual replay, and checked control operations.

The complete two-strand procedure and the upward-pass controller have distinct
output contracts; see their modules before interpreting a completion status.
"""

from .catalog import Catalog, CatalogError

__all__ = ["Catalog", "CatalogError"]
