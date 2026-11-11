"""Exact input witnesses for parabolic blocks and distant-generator commutations.

These certificates concern words and Artin relations only. They do not certify
that two executions of the zipper scheduler have identical outputs or costs.
"""
from __future__ import annotations

from dataclasses import dataclass
import heapq
from typing import Sequence


def _validated(word: Sequence[int], n: int) -> tuple[int, ...]:
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    result = tuple(word)
    if any(type(x) is not int or not 0 < abs(x) < n for x in result):
        raise ValueError("word contains an invalid signed generator")
    return result


def independent(a: int, b: int) -> bool:
    """The two signed generators commute by the distant Artin relation."""
    return abs(abs(a) - abs(b)) >= 2


def trace_permutation(word: Sequence[int]) -> tuple[int, ...]:
    """Indices giving the lexicographically least distant-commutation word.

    Entries must be nonzero signed integers. Occurrence indices break ties,
    making the witness deterministic even for repeated letters.
    """
    letters = tuple(word)
    if any(type(x) is not int or x == 0 for x in letters):
        raise ValueError("word contains an invalid signed generator")
    successors: list[list[int]] = [[] for _ in letters]
    indegrees = [0] * len(letters)
    for j, right in enumerate(letters):
        for i in range(j):
            if not independent(letters[i], right):
                successors[i].append(j)
                indegrees[j] += 1
    ready = [(letters[i], i) for i, degree in enumerate(indegrees) if degree == 0]
    heapq.heapify(ready)
    result = []
    while ready:
        _, i = heapq.heappop(ready)
        result.append(i)
        for j in successors[i]:
            indegrees[j] -= 1
            if indegrees[j] == 0:
                heapq.heappush(ready, (letters[j], j))
    return tuple(result)


def check_commutation(word: Sequence[int], permutation: Sequence[int]) -> bool:
    """Independently verify an occurrence permutation using all inversions."""
    letters, order = tuple(word), tuple(permutation)
    if any(type(x) is not int or x == 0 for x in letters):
        return False
    if any(type(i) is not int for i in order) or sorted(order) != list(range(len(letters))):
        return False
    positions = [0] * len(order)
    for position, occurrence in enumerate(order):
        positions[occurrence] = position
    return all(positions[i] < positions[j] or independent(letters[i], letters[j])
               for j in range(len(letters)) for i in range(j))


@dataclass(frozen=True)
class Block:
    """A connected support interval, embedded by generator index + offset."""

    offset: int
    n: int
    positions: tuple[int, ...]
    word: tuple[int, ...]


@dataclass(frozen=True)
class Coverage:
    n: int
    word: tuple[int, ...]
    commute: bool
    blocks: tuple[Block, ...]

    @property
    def permutation(self) -> tuple[int, ...]:
        return tuple(i for block in self.blocks for i in block.positions)


def decompose(word: Sequence[int], *, n: int, commute: bool = True) -> Coverage:
    """Factor across missing generators, translate blocks, optionally commute.

    Positions retain the complete literal input, including all inverse pairs.
    The empty word has no blocks; ``n`` still records its ambient identity.
    """
    letters = _validated(word, n)
    if type(commute) is not bool:
        raise ValueError("commute must be a boolean")
    intervals: list[list[int]] = []
    for generator in sorted({abs(x) for x in letters}):
        if not intervals or generator != intervals[-1][-1] + 1:
            intervals.append([])
        intervals[-1].append(generator)
    blocks = []
    for interval in intervals:
        start, end = interval[0], interval[-1]
        positions = tuple(i for i, x in enumerate(letters) if start <= abs(x) <= end)
        normalized = tuple((1 if letters[i] > 0 else -1) * (abs(letters[i]) - start + 1)
                           for i in positions)
        order = trace_permutation(normalized) if commute else tuple(range(len(positions)))
        blocks.append(Block(start - 1, end - start + 2,
                            tuple(positions[i] for i in order), tuple(normalized[i] for i in order)))
    result = Coverage(n, letters, commute, tuple(blocks))
    if not check_coverage(result):
        raise AssertionError("constructed an invalid input coverage witness")
    return result


def check_coverage(certificate: Coverage) -> bool:
    """Check factorization without calling the decomposition constructor.

    Lexicographic minimality is checked by a separate greedy criterion: no
    smaller letter can cross all intervening letters to an earlier position.
    """
    if not isinstance(certificate, Coverage) or type(certificate.commute) is not bool:
        return False
    try:
        letters = _validated(certificate.word, certificate.n)
    except (TypeError, ValueError):
        return False
    if not isinstance(certificate.blocks, tuple) or any(
            not isinstance(block, Block) or not isinstance(block.positions, tuple)
            or not isinstance(block.word, tuple) for block in certificate.blocks):
        return False
    if not check_commutation(letters, certificate.permutation):
        return False
    last_end = -1
    for block in certificate.blocks:
        if (not isinstance(block, Block) or type(block.offset) is not int
                or type(block.n) is not int or block.offset < 0 or block.n < 2
                or block.offset + block.n > certificate.n or block.offset <= last_end):
            return False
        # A missing generator must separate consecutive blocks.
        last_end = block.offset + block.n - 1
        if len(block.word) != len(block.positions):
            return False
        try:
            normalized = _validated(block.word, block.n)
        except (TypeError, ValueError):
            return False
        if {abs(x) for x in normalized} != set(range(1, block.n)):
            return False
        for i, x in zip(block.positions, normalized):
            if letters[i] != (1 if x > 0 else -1) * (abs(x) + block.offset):
                return False
        if not certificate.commute:
            if tuple(sorted(block.positions)) != block.positions:
                return False
        else:
            for j, letter in enumerate(normalized):
                for i in range(j - 1, -1, -1):
                    if not independent(normalized[i], letter):
                        break
                    if letter < normalized[i]:
                        return False
    return True


def adjacent_swaps(certificate: Coverage) -> tuple[int, ...]:
    """Expand the witness into adjacent swaps, each index swapping i and i+1."""
    if not check_coverage(certificate):
        raise ValueError("invalid coverage certificate")
    order = list(range(len(certificate.word)))
    swaps = []
    for target, occurrence in enumerate(certificate.permutation):
        position = order.index(occurrence, target)
        while position > target:
            swaps.append(position - 1)
            order[position - 1], order[position] = order[position], order[position - 1]
            position -= 1
    return tuple(swaps)
