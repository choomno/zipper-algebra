"""Exact generic-q braid-group-algebra checks for band transport."""
from collections import defaultdict
from dataclasses import replace

import pytest

from zipper_algebra.zipping_algorithm.transport import (
    BandTransport,
    TransportError,
    TransportStep,
    check_transport,
    transport_band,
)


def _inverse(word):
    return tuple(-x for x in reversed(word))


def _substitute(word, images):
    stack = []
    for x in word:
        image = images[x - 1] if x > 0 else _inverse(images[-x - 1])
        for y in image:
            if stack and stack[-1] == -y:
                stack.pop()
            else:
                stack.append(y)
    return tuple(stack)


def _artin_key(word, n):
    """Faithful Artin action, independent of the transport code."""
    images = tuple((i,) for i in range(1, n + 1))
    for letter in word:
        j = abs(letter)
        generator = [(i,) for i in range(1, n + 1)]
        if letter > 0:
            generator[j - 1:j + 1] = [(j, j + 1, -j), (j,)]
        else:
            generator[j - 1:j + 1] = [(j + 1,), (-j - 1, j, j + 1)]
        images = tuple(_substitute(w, images) for w in generator)
    return images


def _slider(size, start, mirrored):
    """Independent direct defining-polynomial expansion over Z[q,q^-1]."""
    if size == 1:
        return {(0, ()): 1}
    out = {(1, (-start,)): 1, (0, ()): 1, (1, ()): -1, (0, (start,)): -1}
    for level in range(3, size + 1):
        nxt = defaultdict(int)
        positive = tuple(range(start, start + level - 1))
        for (q, word), coefficient in out.items():
            nxt[(q + level - 1, tuple(-x for x in positive) + word)] += coefficient
            nxt[(q, positive + word)] -= coefficient
        out = dict(nxt)
    if mirrored:
        out = {(q, tuple((1 if x > 0 else -1) * (2 * start + size - 2 - abs(x))
                         for x in word)): c for (q, word), c in out.items()}
    return out


def _expanded_side(size, start, mirrored, prefix, suffix, n):
    out = defaultdict(int)
    for (q, word), coefficient in _slider(size, start, mirrored).items():
        out[(q, _artin_key(prefix + word + suffix, n))] += coefficient
    return {key: value for key, value in out.items() if value}


def _assert_identity(witness):
    check_transport(witness)
    assert _expanded_side(witness.size, witness.start, witness.mirrored,
                          (), witness.word, witness.n) == _expanded_side(
        witness.size, witness.end, witness.mirrored, witness.word, (), witness.n)


def test_whole_band_moves_all_small_ranks_generic_q():
    for n in range(2, 7):
        for size in range(1, min(n, 4) + 1):
            for start in range(1, n - size + 2):
                for sign in (-1, 1):
                    contexts = []
                    if start + size <= n:
                        contexts.append(tuple(sign * x for x in range(start + size - 1, start - 1, -1)))
                    if start > 1:
                        contexts.append(tuple(sign * x for x in range(start - 1, start + size - 1)))
                    for word in contexts:
                        for mirrored in (False, True):
                            _assert_identity(transport_band(word, n=n, size=size,
                                                            start=start, mirrored=mirrored))


def test_disjoint_and_two_strand_internal_crossings():
    witness = transport_band((1, -6, 3, -3, 5, 4, 3, 3, 4), n=7, size=2, start=3)
    assert [s.kind for s in witness.steps] == ["far", "far", "internal_i2", "internal_i2",
                                              "far", "right", "left"]
    _assert_identity(witness)


def test_two_strand_internal_twists_without_inverse_cancellation():
    for n in range(2, 8):
        for start in range(1, n):
            for sign in (-1, 1):
                for power in (1, 2, 5):
                    for mirrored in (False, True):
                        witness = transport_band((sign * start,) * power,
                                                 n=n, size=2, start=start, mirrored=mirrored)
                        assert witness.end == start
                        assert all(step.kind == "internal_i2" for step in witness.steps)
                        _assert_identity(witness)


def test_corrected_anchor_retains_complete_processed_word():
    witness = transport_band((-2, -1), n=3, size=2, start=1)
    assert witness.word == (-2, -1)
    assert witness.before == ("slider", 2, 1, False, False, False)
    assert witness.after == ("slider", 2, 2, False, False, False)
    assert witness.steps == (TransportStep("right", 0, 1, 2, (-2, -1)),)
    _assert_identity(witness)
    assert _expanded_side(2, 1, False, (), witness.word, 3) != _expanded_side(
        2, 1, False, witness.word, (), 3)


def test_composed_transport_and_empty_context():
    witness = transport_band((-4, -3, -2, 6, 2, 3, 4), n=7, size=3, start=2, mirrored=True)
    assert witness.end == 2
    assert [s.offset for s in witness.steps] == [0, 3, 4]
    _assert_identity(witness)
    _assert_identity(transport_band((), n=4, size=4, start=1))


def test_reflection_and_rotation_transport_anchors():
    original = transport_band((-4, -3, -2, 6), n=7, size=3, start=2)
    reflected_word = tuple((1 if x > 0 else -1) * (7 - abs(x)) for x in original.word)
    reflected = transport_band(reflected_word, n=7, size=3, start=4, mirrored=True)
    assert reflected.end == 3
    _assert_identity(reflected)
    # phi is reflection followed by reversal, so the two endpoints exchange.
    rotated = transport_band(tuple(reversed(reflected_word)), n=7, size=3, start=3)
    assert rotated.end == 4
    _assert_identity(rotated)


@pytest.mark.parametrize("word,size,start,offset", [
    ((3,), 3, 1, 0),
    ((3, -2, 1), 3, 1, 0),
    ((1, 2), 3, 2, 0),
    ((5, 2), 3, 1, 1),
])
def test_partial_tape_and_internal_crossings_refused(word, size, start, offset):
    with pytest.raises(TransportError) as caught:
        transport_band(word, n=6, size=size, start=start)
    assert caught.value.offset == offset
    assert caught.value.remaining == word[offset:]
    assert "offset" in str(caught.value)


@pytest.mark.parametrize("kwargs", [
    {"n": True}, {"size": False}, {"start": True}, {"n": 0}, {"size": 0},
    {"size": 4}, {"start": 0}, {"start": 3}, {"mirrored": 1},
    {"word": (0,)}, {"word": (3,)}, {"word": (True,)}, {"word": (1.0,)},
])
def test_invalid_input(kwargs):
    data = dict(word=(), n=3, size=2, start=1)
    data.update(kwargs)
    with pytest.raises(ValueError):
        transport_band(**data)


def test_checker_rejects_incomplete_or_forged_witness():
    witness = transport_band((-2, -1), n=3, size=2, start=1)
    step, = witness.steps
    bad_steps = [replace(step, offset=1), replace(step, offset=False),
                 replace(step, start_before=2), replace(step, start_before=True),
                 replace(step, start_after=True), replace(step, block=()),
                 replace(step, block=[-2, -1]), replace(step, block=(True,)),
                 replace(step, block=(3,)), replace(step, kind="unknown"),
                 replace(step, start_after=1), replace(step, block=(-2, 1)),
                 replace(step, kind="far"), replace(step, kind="internal_i2"),
                 replace(step, kind="singleton"), replace(step, kind="left")]
    for bad_step in bad_steps:
        with pytest.raises(ValueError):
            check_transport(replace(witness, steps=(bad_step,)))
    for bad in (replace(witness, end=True), replace(witness, end=0),
                replace(witness, end=1), replace(witness, steps=()),
                replace(witness, word=(-2, -1, 2))):
        with pytest.raises(ValueError):
            check_transport(bad)
    with pytest.raises(ValueError):
        check_transport(BandTransport(2, 1, 1, 1, False, (1,),
                                      (TransportStep("right", 0, 1, 0, (1,)),)))


def test_checker_rejects_malformed_record_containers():
    witness = transport_band((-2, -1), n=3, size=2, start=1)
    for malformed in ({}, replace(witness, word=[-2, -1]),
                      replace(witness, steps=list(witness.steps)),
                      replace(witness, steps=({"kind": "right"},))):
        with pytest.raises(ValueError):
            check_transport(malformed)
