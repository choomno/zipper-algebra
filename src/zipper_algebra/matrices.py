"""Immutable sparse matrices over explicitly declared exact coefficient fields."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json

from .fields import FieldElement, RationalFunctionField


@dataclass(frozen=True, slots=True)
class ExactMatrix:
    """A matrix acting on column vectors, with explicit dimensions even when empty.

    Entries are triples ``(row, column, coefficient)``. Coordinates are zero based;
    repeated coordinates are rejected before zero coefficients are removed.
    """

    field: RationalFunctionField
    nrows: int
    ncols: int
    entries: tuple[tuple[int, int, FieldElement], ...] = ()

    def __post_init__(self):
        if not isinstance(self.field, RationalFunctionField):
            raise TypeError("a matrix requires a declared coefficient field")
        if any(type(n) is not int or n < 0 for n in (self.nrows, self.ncols)):
            raise ValueError("matrix dimensions must be nonnegative integers")
        seen = set()
        entries = []
        for row, column, coefficient in self.entries:
            if (type(row) is not int or type(column) is not int
                    or not 0 <= row < self.nrows or not 0 <= column < self.ncols):
                raise ValueError("matrix coordinate is outside its declared shape")
            if (row, column) in seen:
                raise ValueError("matrix coordinates must be distinct")
            seen.add((row, column))
            value = self.field.scalar(coefficient)
            if not value.is_zero():
                entries.append((row, column, value))
        object.__setattr__(self, "entries", tuple(sorted(entries)))

    @classmethod
    def from_rows(cls, field, rows, *, ncols=None):
        """Snapshot dense rows; specify ``ncols`` to represent a zero-row matrix."""
        rows = tuple(tuple(row) for row in rows)
        width = (len(rows[0]) if rows else 0) if ncols is None else ncols
        if any(len(row) != width for row in rows):
            raise ValueError("matrix rows must have the declared common length")
        return cls(field, len(rows), width,
                   tuple((i, j, value) for i, row in enumerate(rows)
                         for j, value in enumerate(row)))

    def to_json(self):
        return {"schema": "zipper-matrix-v1", "field": self.field.to_json(),
                "shape": [self.nrows, self.ncols], "orientation": "column-vectors",
                "entries": [[i, j, c.to_json()] for i, j, c in self.entries]}

    @classmethod
    def from_json(cls, value):
        if (not isinstance(value, dict)
                or set(value) != {"schema", "field", "shape", "orientation", "entries"}
                or value["schema"] != "zipper-matrix-v1"
                or value["orientation"] != "column-vectors"
                or not isinstance(value["shape"], list) or len(value["shape"]) != 2
                or not isinstance(value["entries"], list)
                or any(not isinstance(entry, list) or len(entry) != 3
                       for entry in value["entries"])):
            raise ValueError("invalid exact matrix schema")
        field = RationalFunctionField.from_json(value["field"])
        return cls(field, value["shape"][0], value["shape"][1],
                   tuple((i, j, FieldElement.from_json(c)) for i, j, c in value["entries"]))

    def content_hash(self):
        """SHA-256 of canonical exact entries, dimensions, and complete field."""
        payload = json.dumps(self.to_json(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()
