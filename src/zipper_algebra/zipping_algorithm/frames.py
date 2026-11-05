"""Exact transport of a marked tape carrying a half-twist power.

For ``side='left'`` the witness proves U H = H' U'; for ``side='right'``
it proves H U = U' H'. Here H is Delta_(size,start)^power. Products and
offsets are in algebraic order. Internal tape crossings are reflected by
conjugation when power is odd; whole-tape crossings shift the anchor.

These are braid identities over any coefficient ring. They provide a
checked primitive for retained flype frames, without asserting that every
reachable zipper context admits the required marked-tape decomposition.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


class FrameError(ValueError):
    """The next crossing does not contain a complete marked tape."""


@dataclass(frozen=True)
class FrameStep:
    kind: str
    offset: int
    start_before: int
    start_after: int
    block: tuple[int, ...]
    emitted: tuple[int, ...]


@dataclass(frozen=True)
class FrameTransport:
    n: int
    size: int
    start: int
    end: int
    power: int
    side: str
    word: tuple[int, ...]
    transformed: tuple[int, ...]
    steps: tuple[FrameStep, ...]


def _validate(
    word: Sequence[int], n: int, size: int, start: int, power: int, side: str
) -> tuple[int, ...]:
    if (
        any(type(x) is not int for x in (n, size, start, power))
        or not 1 <= size <= n
        or not 1 <= start <= n - size + 1
        or side not in ("left", "right")
    ):
        raise ValueError("valid rank, tape interval, integer power and side required")
    result = tuple(word)
    if any(type(x) is not int or not 1 <= abs(x) < n for x in result):
        raise ValueError("word must contain signed braid generators in B_n")
    return result


def _local(
    letter: int, size: int, start: int, power: int, side: str
) -> tuple[str, int, tuple[int, ...], tuple[int, ...]]:
    index, sign = abs(letter), 1 if letter > 0 else -1
    if index < start - 1 or index >= start + size:
        return "far", start, (letter,), (letter,)
    if start <= index < start + size - 1:
        # Delta sigma_r Delta^-1 = sigma_(2j+t-2-r), and the
        # conjugation has order two; this works for every integer power.
        reflected = 2 * start + size - 2 - index if power % 2 else index
        return "internal", start, (letter,), (sign * reflected,)
    if size == 1:
        after = start + (1 if index == start else -1)
        return "singleton", after, (letter,), (letter,)
    if index == (start + size - 1 if side == "left" else start - 1):
        after = start + (1 if side == "left" else -1)
        low = start if side == "left" else start - 1
        block = tuple(sign * x for x in range(low, low + size))
    else:
        after = start + (-1 if side == "left" else 1)
        low = start - 1 if side == "left" else start
        block = tuple(sign * x for x in range(low + size - 1, low - 1, -1))
    # Whole-band naturality transports every internal braid generator by
    # an index shift, hence every power of the internal half twist.
    return "right" if after > start else "left", after, block, block


def transport_frame(
    word: Sequence[int],
    *,
    n: int,
    size: int,
    start: int,
    power: int = 1,
    side: str = "left",
) -> FrameTransport:
    """Transport a half-twist power, retaining a complete structural witness.

    Left transport scans U from its algebraic right end; right transport
    scans from its algebraic left end. No half-twist expansion is needed.
    The input and output crossing counts agree exactly, including power zero.
    """
    braid = _validate(word, n, size, start, power, side)
    position = start
    cursor = len(braid) if side == "left" else 0
    pieces: list[tuple[int, ...]] = []
    steps: list[FrameStep] = []
    while cursor > 0 if side == "left" else cursor < len(braid):
        letter = braid[cursor - 1] if side == "left" else braid[cursor]
        kind, after, block, emitted = _local(letter, size, position, power, side)
        offset = cursor - len(block) if side == "left" else cursor
        if (
            offset < 0
            or braid[offset : offset + len(block)] != block
            or not 1 <= after <= n - size + 1
        ):
            raise FrameError(
                f"whole marked tape required at offset {offset}; "
                f"tape starts at {position}, side={side}"
            )
        steps.append(FrameStep(kind, offset, position, after, block, emitted))
        pieces.append(emitted)
        cursor = offset if side == "left" else cursor + len(block)
        position = after
    if side == "left":
        pieces.reverse()
    transformed = tuple(x for piece in pieces for x in piece)
    return FrameTransport(
        n, size, start, position, power, side, braid, transformed, tuple(steps)
    )


def check_frame(witness: FrameTransport) -> None:
    """Check local Garside/naturality identities and full word coverage.

    This does not rerun the producer parser. It checks the supplied moves,
    their actual input intervals and the entire emitted word.
    """
    if (
        not isinstance(witness, FrameTransport)
        or not isinstance(witness.word, tuple)
        or not isinstance(witness.transformed, tuple)
        or not isinstance(witness.steps, tuple)
    ):
        raise ValueError("immutable frame witness required")
    braid = _validate(
        witness.word,
        witness.n,
        witness.size,
        witness.start,
        witness.power,
        witness.side,
    )
    if (
        type(witness.end) is not int
        or not 1 <= witness.end <= witness.n - witness.size + 1
    ):
        raise ValueError("invalid final tape interval")
    position = witness.start
    cursor = len(braid) if witness.side == "left" else 0
    pieces: list[tuple[int, ...]] = []
    for step in witness.steps:
        if (
            not isinstance(step, FrameStep)
            or any(
                type(x) is not int
                for x in (step.offset, step.start_before, step.start_after)
            )
            or not isinstance(step.block, tuple)
            or not step.block
            or not isinstance(step.emitted, tuple)
            or any(
                type(x) is not int or not 1 <= abs(x) < witness.n
                for x in step.block + step.emitted
            )
            or step.start_before != position
        ):
            raise ValueError("invalid frame step")
        # The step's declared local relation is reconstructed without
        # searching for another parse or invoking transport_frame.
        letter = step.block[-1] if witness.side == "left" else step.block[0]
        expected = _local(letter, witness.size, position, witness.power, witness.side)
        offset = cursor - len(step.block) if witness.side == "left" else cursor
        if (
            (step.kind, step.start_after, step.block, step.emitted) != expected
            or step.offset != offset
            or offset < 0
            or braid[offset : offset + len(step.block)] != step.block
            or not 1 <= step.start_after <= witness.n - witness.size + 1
        ):
            raise ValueError("invalid local frame identity")
        pieces.append(step.emitted)
        cursor = offset if witness.side == "left" else cursor + len(step.block)
        position = step.start_after
    if witness.side == "left":
        pieces.reverse()
    transformed = tuple(x for piece in pieces for x in piece)
    boundary = 0 if witness.side == "left" else len(braid)
    if (
        cursor != boundary
        or position != witness.end
        or transformed != witness.transformed
    ):
        raise ValueError("frame witness does not cover the entire identity")
