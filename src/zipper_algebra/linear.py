"""Exact matrix rank over the rationals or an explicitly checked prime field."""
from __future__ import annotations

from fractions import Fraction
from math import isqrt


def validate_field_characteristic(characteristic: int) -> int:
    """Accept characteristic zero or a prime, using exact trial division."""
    if type(characteristic) is not int or characteristic < 0 or characteristic == 1:
        raise ValueError("field characteristic must be zero or a prime integer")
    if characteristic and any(characteristic % d == 0
                              for d in range(2, isqrt(characteristic) + 1)):
        raise ValueError("field characteristic must be zero or a prime integer")
    return characteristic


def _scalar(value, characteristic):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("matrix entries must be exact integers or Fractions")
    value = Fraction(value)
    if characteristic == 0:
        return value
    try:
        inverse = pow(value.denominator, -1, characteristic)
    except ValueError:
        raise ZeroDivisionError("matrix denominator vanishes in the prime field") from None
    return value.numerator * inverse % characteristic


def matrix_rank(rows, *, characteristic: int = 0) -> int:
    """Return the exact rank of a rectangular matrix, without mutating its rows.

    A rank over a prime field is the rank of that supplied matrix. It is not
    asserted to equal the rational rank or a generic rank over parameters.
    """
    validate_field_characteristic(characteristic)
    matrix = [[_scalar(value, characteristic) for value in row] for row in rows]
    ncols = len(matrix[0]) if matrix else 0
    if any(len(row) != ncols for row in matrix):
        raise ValueError("matrix rows must have equal lengths")
    rank = 0
    for column in range(ncols):
        pivot = next((r for r in range(rank, len(matrix)) if matrix[r][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = (pow(matrix[rank][column], -1, characteristic) if characteristic
                   else 1 / matrix[rank][column])
        for row in range(rank + 1, len(matrix)):
            if not matrix[row][column]:
                continue
            factor = matrix[row][column] * inverse
            for col in range(column, ncols):
                value = matrix[row][col] - factor * matrix[rank][col]
                matrix[row][col] = value % characteristic if characteristic else value
        rank += 1
        if rank == len(matrix):
            break
    return rank
