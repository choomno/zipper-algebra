"""Full-support literal inputs remaining after initial singleton peeling.

The distinguished original bottom strand n starts at position n. Reading the
input RIGHT TO LEFT, a crossing adjacent to its current position moves it to
the other position; it is over exactly when (position == abs(letter)) equals
(letter > 0). If it is never under, the descending-prefix theorem reduces the
calculation to a smaller-rank literal core plus deterministic extraction data.

This module enumerates the complementary inputs: full generator support and
at least one undercrossing of strand n. No commutation, cancellation, or braid
relation is imposed. Ordering is colexicographic: reversed words are ordered
lexicographically in the usual signed alphabet. It is deliberately distinct
from EssentialWords' left-to-right ordinals. Rank one has no such inputs; the
global identity must be handled separately by any coverage construction.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Iterator, Sequence

from .enumeration import _integer, alphabet


def has_bottom_n_undercrossing(word: Sequence[int], *, n: int) -> bool:
    """Whether original bottom strand n passes under another input strand.

    Validate the entire word even after an undercrossing has been encountered.
    """
    letters = alphabet(n)
    position, seen = n, False
    for letter in reversed(word):
        if type(letter) is not int or letter not in letters:
            raise ValueError("word contains an invalid generator")
        index = abs(letter)
        if position in (index, index + 1):
            seen |= (position == index) != (letter > 0)
            position = 2 * index + 1 - position
    return seen


class NonPeelableWords:
    """Exact counts and random access for one full-support unpeeled layer.

    Completion state is remaining length, missing support, distinguished-strand
    position, and whether an undercrossing has occurred. After it occurs only
    the number of missing generators matters. Before it occurs the desired
    count is all full-support completions minus those avoiding undercrossings.
    The latter automaton uses two signed choices for a far generator and one
    overcrossing choice for an adjacent generator. Caches belong to this
    instance, contain immutable states, and do not retain the instance itself.
    """

    def __init__(self, n: int, length: int) -> None:
        self.letters = alphabet(n)
        _integer(length, "length", 0)
        self.n, self.length = n, length
        self._full = (1 << (n - 1)) - 1
        self._support = tuple(1 << (abs(a) - 1) for a in self.letters)

        @lru_cache(maxsize=None)
        def literal(remaining: int, missing_count: int) -> int:
            if missing_count > remaining:
                return 0
            if remaining == 0:
                return int(missing_count == 0)
            value = 2 * (n - 1 - missing_count) * literal(remaining - 1, missing_count)
            if missing_count:
                value += 2 * missing_count * literal(remaining - 1, missing_count - 1)
            return value

        @lru_cache(maxsize=None)
        def avoiding(remaining: int, missing: int, position: int) -> int:
            if missing.bit_count() > remaining:
                return 0
            if remaining == 0:
                return int(missing == 0)
            result = 0
            for index in range(1, n):
                adjacent = position in (index, index + 1)
                following = 2 * index + 1 - position if adjacent else position
                result += (1 if adjacent else 2) * avoiding(
                    remaining - 1, missing & ~(1 << (index - 1)), following)
            return result

        self._literal = literal
        self._avoiding = avoiding

    def _count(self, remaining: int, missing: int, position: int, seen: bool) -> int:
        total = self._literal(remaining, missing.bit_count())
        return total if seen else total - self._avoiding(remaining, missing, position)

    def _next(self, missing: int, position: int, seen: bool, j: int) -> tuple[int, int, bool]:
        letter = self.letters[j]
        index = abs(letter)
        if position in (index, index + 1):
            seen |= (position == index) != (letter > 0)
            position = 2 * index + 1 - position
        return missing & ~self._support[j], position, seen

    @property
    def count(self) -> int:
        """Exact number of retained literal words."""
        return self._count(self.length, self._full, self.n, False)

    def word_at(self, ordinal: int) -> tuple[int, ...]:
        """Unrank a word in colex order, without scanning earlier words."""
        _integer(ordinal, "ordinal", 0)
        if ordinal >= self.count:
            raise IndexError("ordinal is outside the nonpeelable layer")
        state = (self._full, self.n, False)
        backward: list[int] = []
        for remaining in range(self.length - 1, -1, -1):
            # The validated ordinal belongs to a child interval.
            j = 0
            while True:
                letter = self.letters[j]
                following = self._next(*state, j)
                size = self._count(remaining, *following)
                if ordinal >= size:
                    ordinal -= size
                    j += 1
                else:
                    backward.append(letter)
                    state = following
                    break
        return tuple(reversed(backward))

    def ordinal(self, word: Sequence[int]) -> int:
        """Rank a valid word; reject missing support or a peelable input."""
        if len(word) != self.length:
            raise ValueError("word has the wrong length")
        state = (self._full, self.n, False)
        ordinal = 0
        for offset, letter in enumerate(reversed(word)):
            if type(letter) is not int or letter not in self.letters:
                raise ValueError("word contains an invalid generator")
            chosen = self.letters.index(letter)
            for j in range(chosen):
                ordinal += self._count(self.length - offset - 1, *self._next(*state, j))
            state = self._next(*state, chosen)
        if state[0]:
            raise ValueError("word does not have full generator support")
        if not state[2]:
            raise ValueError("word admits initial singleton peeling")
        return ordinal

    def iter_words(self, start: int = 0, stop: int | None = None) -> Iterator[tuple[int, ...]]:
        """Yield a colex ordinal interval, pruning whole excluded subtrees."""
        _integer(start, "start", 0)
        if stop is None:
            stop = self.count
        _integer(stop, "stop", 0)
        if not start <= stop <= self.count:
            raise ValueError("invalid ordinal interval")
        backward: list[int] = []

        def visit(remaining: int, missing: int, position: int, seen: bool,
                  first: int, last: int) -> Iterator[tuple[int, ...]]:
            if first == last:
                return
            if remaining == 0:
                yield tuple(reversed(backward))
                return
            offset = 0
            j = 0
            while offset < last:
                letter = self.letters[j]
                following = self._next(missing, position, seen, j)
                size = self._count(remaining - 1, *following)
                left, right = max(first, offset), min(last, offset + size)
                if left < right:
                    backward.append(letter)
                    yield from visit(remaining - 1, *following, left - offset, right - offset)
                    backward.pop()
                offset += size
                j += 1

        yield from visit(self.length, self._full, self.n, False, start, stop)
