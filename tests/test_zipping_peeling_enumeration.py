"""Literal label tracking independently checks the peeling automaton."""
from __future__ import annotations

import itertools
import random
import weakref

import pytest

from zipper_algebra.zipping_algorithm.enumeration import alphabet
from zipper_algebra.zipping_algorithm.essential_enumeration import EssentialWords
from zipper_algebra.zipping_algorithm.peeling_enumeration import (
    NonPeelableWords, has_bottom_n_undercrossing,
)


def label_scan(word, n):
    labels = list(range(1, n + 1))
    under = False
    for letter in reversed(word):
        j = abs(letter) - 1
        left, right = labels[j:j + 2]
        lower = right if letter > 0 else left
        under |= lower == n
        labels[j:j + 2] = [right, left]
    return under


@pytest.mark.parametrize("n,length", [(1, 0), (1, 2), (2, 0), (2, 5), (3, 4), (4, 4), (5, 5)])
def test_exhaustive_independent_labels_and_random_access(n, length):
    expected = []
    for word in itertools.product(alphabet(n), repeat=length):
        under = label_scan(word, n)
        assert has_bottom_n_undercrossing(word, n=n) == under
        if under and {abs(a) for a in word} == set(range(1, n)):
            expected.append(word)
    expected.sort(key=lambda w: w[::-1])
    layer = NonPeelableWords(n, length)
    assert layer.count == len(expected)
    assert list(layer.iter_words()) == expected
    for i, word in enumerate(expected):
        assert layer.ordinal(word) == i
        assert layer.word_at(i) == word
    for start in range(0, len(expected) + 1, max(1, len(expected) // 7)):
        stop = min(len(expected), start + 13)
        assert list(layer.iter_words(start, stop)) == expected[start:stop]


def test_large_counts_seek_and_inverse_pairs():
    layer = NonPeelableWords(8, 10)
    assert layer.count == 22273020209
    assert EssentialWords(8, 10, commute=False).count - layer.count == 8073424591
    randomizer = random.Random(294)
    for ordinal in [0, layer.count - 1, *(randomizer.randrange(layer.count) for _ in range(20))]:
        word = layer.word_at(ordinal)
        assert label_scan(word, 8)
        assert layer.ordinal(word) == ordinal
        assert list(layer.iter_words(ordinal, ordinal + 1)) == [word]
    # The domain retains inverse pairs whenever the distinguished strand goes under.
    assert (-1, 1) in list(NonPeelableWords(2, 2).iter_words())
    for length in range(1, 11):
        assert NonPeelableWords(2, length).count == 2 ** length - 1


@pytest.mark.parametrize("n,length", [(0, 2), (True, 2), (2, -1), (2, True)])
def test_invalid_layer(n, length):
    with pytest.raises(ValueError):
        NonPeelableWords(n, length)


def test_invalid_ordinals_words_and_intervals():
    layer = NonPeelableWords(3, 2)
    for ordinal in (-1, True, 1.5):
        with pytest.raises(ValueError):
            layer.word_at(ordinal)
    with pytest.raises(IndexError):
        layer.word_at(layer.count)
    for word in [(1,), (1, 1), (0, 1), (True, 1), (-1, -2)]:
        with pytest.raises(ValueError):
            layer.ordinal(word)
    for start, stop in [(-1, 0), (2, 1), (0, layer.count + 1), (True, 2)]:
        with pytest.raises(ValueError):
            list(layer.iter_words(start, stop))
    assert list(layer.iter_words(layer.count, layer.count)) == []
    with pytest.raises(ValueError):
        NonPeelableWords(1, 0).ordinal(())
    for word in [(0,), (True,), (2,), (0, 1)]:
        with pytest.raises(ValueError):
            has_bottom_n_undercrossing(word, n=2)
    with pytest.raises(ValueError):
        has_bottom_n_undercrossing((), n=0)


def test_caches_do_not_retain_layer():
    layer = NonPeelableWords(5, 7)
    assert layer.count > 0
    ref = weakref.ref(layer)
    del layer
    assert ref() is None
