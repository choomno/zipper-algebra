"""Explicit transport through whole-band crossings in the braid group algebra.

If ``P`` passes this parser, its witness proves ``X_(i,j) P = P X_(i,j')``
over Z[q,q^-1]. Products and offsets use algebraic left-to-right order.
The same identities hold for the mirrored slider: each whole-band crossing
transports every braid generator in the band by shifting its index, hence
transports every polynomial in those generators. Disjoint crossings commute.
For i=2 the internal braid group is cyclic, so internal twists also commute.

This is a sufficient, explicitly checked class of transport contexts. A
partial tape crossing arising from a flype requires additional half twists;
it is deliberately rejected here. No completed-phase invariant is inferred.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


class TransportError(ValueError):
    """An unsupported context, with its exact unconsumed word position."""

    def __init__(self, reason: str, *, offset: int, start: int,
                 remaining: tuple[int, ...]):
        self.offset = offset
        self.start = start
        self.remaining = remaining
        super().__init__(f"{reason} at offset {offset}, band starts at {start}; "
                         f"remaining={remaining!r}")


@dataclass(frozen=True)
class TransportStep:
    kind: str
    offset: int
    start_before: int
    start_after: int
    block: tuple[int, ...]


@dataclass(frozen=True)
class BandTransport:
    n: int
    size: int
    start: int
    end: int
    mirrored: bool
    word: tuple[int, ...]
    steps: tuple[TransportStep, ...]

    @property
    def before(self) -> tuple[str, int, int, bool, bool, bool]:
        return ("slider", self.size, self.start, self.mirrored, False, False)

    @property
    def after(self) -> tuple[str, int, int, bool, bool, bool]:
        return ("slider", self.size, self.end, self.mirrored, False, False)


def _validate(word: Sequence[int], n: int, size: int, start: int,
              mirrored: bool) -> tuple[int, ...]:
    if (type(n) is not int or type(size) is not int or type(start) is not int
            or not 1 <= size <= n or not 1 <= start <= n - size + 1
            or type(mirrored) is not bool):
        raise ValueError("valid rank, band interval and boolean orientation required")
    result = tuple(word)
    if any(type(x) is not int or not 1 <= abs(x) < n for x in result):
        raise ValueError("word must contain signed braid generators in B_n")
    return result


def transport_band(word: Sequence[int], *, n: int, size: int, start: int,
                   mirrored: bool = False) -> BandTransport:
    """Parse a sufficient band context and retain every letter unchanged.

    A rightward move consumes descending generators from j+i-1 to j, all
    with the same sign; a leftward move consumes ascending generators from
    j-1 to j+i-2. Each move changes the slider anchor by one. Failure never
    commutes a partial block or an arbitrary internal crossing.
    """
    braid = _validate(word, n, size, start, mirrored)
    position = start
    offset = 0
    steps: list[TransportStep] = []
    while offset < len(braid):
        letter = braid[offset]
        index = abs(letter)
        sign = 1 if letter > 0 else -1
        after = position
        block: tuple[int, ...]
        if index < position - 1 or index >= position + size:
            kind, block = "far", (letter,)
        elif size == 1:
            kind, block = "singleton", (letter,)
            after += 1 if index == position else -1
        elif size == 2 and index == position:
            kind, block = "internal_i2", (letter,)
        elif index == position + size - 1:
            kind = "right"
            block = tuple(sign * x for x in range(position + size - 1, position - 1, -1))
            after += 1
        elif index == position - 1:
            kind = "left"
            block = tuple(sign * x for x in range(position - 1, position + size - 1))
            after -= 1
        else:
            raise TransportError("internal crossing requires a flype or another explicit witness",
                                 offset=offset, start=position, remaining=braid[offset:])
        if braid[offset:offset + len(block)] != block:
            raise TransportError("whole-band crossing required",
                                 offset=offset, start=position, remaining=braid[offset:])
        steps.append(TransportStep(kind, offset, position, after, block))
        position = after
        offset += len(block)
    return BandTransport(n, size, start, position, mirrored, braid, tuple(steps))


def check_transport(witness: BandTransport) -> None:
    """Check a saved structural witness without rerunning the producer parser.

    This checks the local naturality identities and their concatenation. It
    does not certify that an algorithm's processed word has this structure.
    """
    if (not isinstance(witness, BandTransport) or not isinstance(witness.word, tuple)
            or not isinstance(witness.steps, tuple)):
        raise ValueError("immutable band transport witness required")
    braid = _validate(witness.word, witness.n, witness.size, witness.start, witness.mirrored)
    if type(witness.end) is not int or not 1 <= witness.end <= witness.n - witness.size + 1:
        raise ValueError("invalid final band interval")
    position, offset = witness.start, 0
    for step in witness.steps:
        if not isinstance(step, TransportStep):
            raise ValueError("transport step record required")
        if (type(step.offset) is not int or type(step.start_before) is not int
                or type(step.start_after) is not int or step.offset != offset
                or step.start_before != position or not step.block
                or not isinstance(step.block, tuple)
                or any(type(x) is not int or not 1 <= abs(x) < witness.n for x in step.block)):
            raise ValueError("invalid transport step location or block")
        sign = 1 if step.block[0] > 0 else -1
        index = abs(step.block[0])
        i = witness.size
        if step.kind == "far":
            valid = (len(step.block) == 1 and step.start_after == position
                     and (index < position - 1 or index >= position + i))
        elif step.kind == "internal_i2":
            valid = (i == 2 and len(step.block) == 1 and index == position
                     and step.start_after == position)
        elif step.kind == "singleton":
            valid = (i == 1 and len(step.block) == 1
                     and index in (position - 1, position)
                     and step.start_after == position + (1 if index == position else -1))
        elif step.kind == "right":
            valid = (step.block == tuple(sign * x for x in range(position + i - 1, position - 1, -1))
                     and step.start_after == position + 1)
        elif step.kind == "left":
            valid = (step.block == tuple(sign * x for x in range(position - 1, position + i - 1))
                     and step.start_after == position - 1)
        else:
            valid = False
        if (not valid or braid[offset:offset + len(step.block)] != step.block
                or not 1 <= step.start_after <= witness.n - i + 1):
            raise ValueError("invalid local band transport identity")
        position = step.start_after
        offset += len(step.block)
    if offset != len(braid) or position != witness.end:
        raise ValueError("transport witness does not cover the entire word")
