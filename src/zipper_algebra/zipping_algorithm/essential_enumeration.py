"""Enumerate full-support signed words, optionally modulo far commutations.

At rank n every index 1,...,n-1 must occur. No inverse cancellation or
three-letter braid relation is imposed. The optional representatives are the
lexicographically least words in their far-commutation classes, using the same
signed alphabet order as literal enumeration. This is a combinatorial domain;
it does not itself certify reuse of a zipper calculation for other inputs.

A word fails lexicographic minimality precisely when it contains b u a with
b > a and a independent of every letter of b u. Such an a can move left;
conversely, a smaller representative's first differing letter must move left
past exactly such an independent suffix. After a prefix, F records the letters
that would create this obstruction. Appending a is allowed iff a is not in F,
and updates F to (F intersect I(a)) union {b < a : b in I(a)}. Here I(a) is
the signed alphabet independent of a. This finite automaton, together with the
missing-support mask, provides the completion counts used below.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Iterator, Sequence

from .enumeration import _integer, alphabet


class EssentialWords:
    """Random-access, lexicographic enumeration of one finite essential layer.

    A prefix automaton rejects a letter if it could move left past commuting
    letters to produce a smaller word. Dynamic programming counts admissible
    completions with all generator indices present, so seeking does not scan
    rejected literal words or earlier representatives. Counts use Python ints.
    An instance owns its cache; release it when its rank/length layer is done.
    Rank one contains the empty word only.
    """

    def __init__(self, n: int, length: int, *, commute: bool = True) -> None:
        self.letters = alphabet(n)
        _integer(length, "length", 0)
        if type(commute) is not bool:
            raise ValueError("commute must be a boolean")
        self.n = n
        self.length = length
        self.commute = commute
        self._full = (1 << (n - 1)) - 1
        self._support = tuple(1 << (abs(a) - 1) for a in self.letters)
        self._independent = tuple(sum(1 << j for j, b in enumerate(self.letters)
                                     if abs(abs(a) - abs(b)) >= 2)
                                  for a in self.letters)
        self._smaller = tuple(mask & ((1 << j) - 1)
                              for j, mask in enumerate(self._independent))

        @lru_cache(maxsize=None)
        def count_literal(remaining: int, missing_count: int) -> int:
            if missing_count > remaining:
                return 0
            if remaining == 0:
                return int(missing_count == 0)
            value = 2 * (n - 1 - missing_count) * count_literal(remaining - 1, missing_count)
            if missing_count:
                value += 2 * missing_count * count_literal(remaining - 1, missing_count - 1)
            return value

        # Capture immutable transition tables, not the owning layer instance.
        support, independent, smaller = self._support, self._independent, self._smaller
        letter_count = len(self.letters)

        @lru_cache(maxsize=None)
        def count_trace(remaining: int, missing: int, forbidden: int) -> int:
            if missing.bit_count() > remaining:
                return 0
            if remaining == 0:
                return int(missing == 0)
            return sum(count_trace(remaining - 1, missing & ~support[j],
                                   (forbidden & independent[j]) | smaller[j])
                       for j in range(letter_count) if not forbidden & (1 << j))

        self._literal_count = count_literal
        self._trace_count = count_trace

    def _count(self, remaining: int, missing: int, forbidden: int) -> int:
        if self.commute:
            return self._trace_count(remaining, missing, forbidden)
        return self._literal_count(remaining, missing.bit_count())

    def _next(self, missing: int, forbidden: int, j: int) -> tuple[int, int]:
        return (missing & ~self._support[j],
                (forbidden & self._independent[j]) | self._smaller[j] if self.commute else 0)

    @property
    def count(self) -> int:
        """Exact number of representatives in this layer."""
        return self._count(self.length, self._full, 0)

    def word_at(self, ordinal: int) -> tuple[int, ...]:
        """Return a zero-based representative without enumerating its predecessors."""
        _integer(ordinal, "ordinal", 0)
        if ordinal >= self.count:
            raise IndexError("ordinal is outside the essential layer")
        missing, forbidden = self._full, 0
        word = []
        for remaining in range(self.length - 1, -1, -1):
            for j, a in enumerate(self.letters):
                if forbidden & (1 << j):
                    continue
                next_missing, next_forbidden = self._next(missing, forbidden, j)
                size = self._count(remaining, next_missing, next_forbidden)
                if ordinal >= size:
                    ordinal -= size
                else:
                    word.append(a)
                    missing, forbidden = next_missing, next_forbidden
                    break
        return tuple(word)

    def ordinal(self, word: Sequence[int]) -> int:
        """Rank a representative, rejecting noncanonical or incomplete support."""
        if len(word) != self.length:
            raise ValueError("word has the wrong length")
        missing, forbidden, ordinal = self._full, 0, 0
        for offset, a in enumerate(word):
            if type(a) is not int or a not in self.letters:
                raise ValueError("word contains an invalid generator")
            chosen = self.letters.index(a)
            if forbidden & (1 << chosen):
                raise ValueError("word is not a far-commutation representative")
            for j in range(chosen):
                if not forbidden & (1 << j):
                    next_missing, next_forbidden = self._next(missing, forbidden, j)
                    ordinal += self._count(self.length - offset - 1, next_missing, next_forbidden)
            missing, forbidden = self._next(missing, forbidden, chosen)
        if missing:
            raise ValueError("word does not have full generator support")
        return ordinal

    def iter_words(self, start: int = 0, stop: int | None = None) -> Iterator[tuple[int, ...]]:
        """Yield an ordinal interval, pruning entire excluded prefix subtrees."""
        _integer(start, "start", 0)
        if stop is None:
            stop = self.count
        _integer(stop, "stop", 0)
        if not start <= stop <= self.count:
            raise ValueError("invalid ordinal interval")
        prefix: list[int] = []

        def visit(remaining: int, missing: int, forbidden: int,
                  first: int, last: int) -> Iterator[tuple[int, ...]]:
            if first == last:
                return
            if remaining == 0:
                yield tuple(prefix)
                return
            offset = 0
            for j, a in enumerate(self.letters):
                if forbidden & (1 << j):
                    continue
                next_missing, next_forbidden = self._next(missing, forbidden, j)
                size = self._count(remaining - 1, next_missing, next_forbidden)
                left, right = max(first, offset), min(last, offset + size)
                if left < right:
                    prefix.append(a)
                    yield from visit(remaining - 1, next_missing, next_forbidden,
                                     left - offset, right - offset)
                    prefix.pop()
                offset += size
                if offset >= last:
                    break

        yield from visit(self.length, self._full, 0, start, stop)
