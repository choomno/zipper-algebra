"""Independent braid identities and persisted-witness checks for flype frames."""

from dataclasses import replace

import pytest

from zipper_algebra.zipping_algorithm.frames import (
    FrameError,
    FrameStep,
    check_frame,
    transport_frame,
)


def _inverse(word):
    return tuple(-x for x in reversed(word))


def _substitute(word, images):
    stack = []
    for x in word:
        for y in images[x - 1] if x > 0 else _inverse(images[-x - 1]):
            if stack and stack[-1] == -y:
                stack.pop()
            else:
                stack.append(y)
    return tuple(stack)


def _key(word, n):
    """Faithful Artin action independent of the frame relation table."""
    images = tuple((i,) for i in range(1, n + 1))
    for letter in word:
        j = abs(letter)
        elementary = [(i,) for i in range(1, n + 1)]
        elementary[j - 1 : j + 1] = (
            [(j, j + 1, -j), (j,)] if letter > 0 else [(j + 1,), (-j - 1, j, j + 1)]
        )
        images = tuple(_substitute(w, images) for w in elementary)
    return images


def _half(size, start, power):
    half = tuple(
        x for length in range(size - 1, 0, -1) for x in range(start, start + length)
    )
    return (half if power >= 0 else _inverse(half)) * abs(power)


def _check_identity(witness):
    check_frame(witness)
    old = _half(witness.size, witness.start, witness.power)
    new = _half(witness.size, witness.end, witness.power)
    if witness.side == "left":
        assert _key(witness.word + old, witness.n) == _key(
            new + witness.transformed, witness.n
        )
    else:
        assert _key(old + witness.word, witness.n) == _key(
            witness.transformed + new, witness.n
        )


def test_frame_identities_all_small_intervals_signs_and_both_sides():
    for n in range(2, 7):
        for size in range(1, n + 1):
            for start in range(1, n - size + 2):
                for side in ("left", "right"):
                    for sign in (-1, 1):
                        contexts = [(sign * r,) for r in range(start, start + size - 1)]
                        contexts += [
                            (sign * r,)
                            for r in range(1, n)
                            if r < start - 1 or r >= start + size
                        ]
                        if start + size <= n:
                            block = tuple(sign * x for x in range(start, start + size))
                            contexts.append(block if side == "left" else block[::-1])
                        if start > 1:
                            block = tuple(
                                sign * x for x in range(start - 1, start + size - 1)
                            )
                            contexts.append(block[::-1] if side == "left" else block)
                        for word in contexts:
                            for power in (-2, -1, 0, 1, 2):
                                _check_identity(
                                    transport_frame(
                                        word,
                                        n=n,
                                        size=size,
                                        start=start,
                                        power=power,
                                        side=side,
                                    )
                                )


def test_composed_left_and_right_transport_preserves_crossing_count():
    left = transport_frame((2, 1, 2, 6, 1), n=7, size=2, start=1, power=-3)
    assert [step.kind for step in left.steps] == [
        "internal",
        "far",
        "right",
        "internal",
    ]
    assert len(left.word) == len(left.transformed)
    _check_identity(left)
    right = transport_frame(
        tuple(reversed(left.word)), n=7, size=2, start=1, power=-3, side="right"
    )
    assert right.transformed == tuple(reversed(left.transformed))
    _check_identity(right)
    _check_identity(transport_frame((), n=1, size=1, start=1, power=1000000))


@pytest.mark.parametrize(
    "side,word,start",
    [("left", (3,), 1), ("right", (3,), 1), ("left", (2, 1), 2), ("right", (1, 2), 2)],
)
def test_partial_tape_refused(side, word, start):
    with pytest.raises(FrameError, match="whole marked tape"):
        transport_frame(word, n=5, size=3, start=start, side=side)


@pytest.mark.parametrize(
    "changed",
    [
        {"n": True},
        {"size": False},
        {"start": 1.0},
        {"power": True},
        {"n": 0},
        {"size": 0},
        {"size": 4},
        {"start": 0},
        {"start": 3},
        {"side": "up"},
        {"word": (0,)},
        {"word": (3,)},
        {"word": (True,)},
    ],
)
def test_invalid_inputs(changed):
    data = dict(word=(), n=3, size=2, start=1)
    data.update(changed)
    with pytest.raises(ValueError):
        transport_frame(**data)


def test_checker_rejects_malformed_and_falsified_complete_witnesses():
    witness = transport_frame((1, 2), n=3, size=2, start=1)
    step = witness.steps[0]
    bad = [
        None,
        replace(witness, word=[1, 2]),
        replace(witness, transformed=[]),
        replace(witness, steps=[]),
        replace(witness, end=True),
        replace(witness, end=4),
        replace(witness, steps=(None,)),
        replace(witness, transformed=(2, 1)),
        replace(witness, steps=()),
        replace(witness, word=(-1, -2)),
    ]
    for change in (
        {"kind": "unknown"},
        {"offset": 1},
        {"offset": True},
        {"start_before": 2},
        {"start_after": 1},
        {"block": ()},
        {"block": [1, 2]},
        {"emitted": []},
        {"block": (True,)},
        {"emitted": (0,)},
        {"block": (2,)},
        {"emitted": (2, 1)},
    ):
        bad.append(replace(witness, steps=(replace(step, **change),)))
    for invalid in bad:
        with pytest.raises(ValueError):
            check_frame(invalid)
    assert isinstance(step, FrameStep)
