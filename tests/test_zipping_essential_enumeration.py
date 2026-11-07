"""Independent orbit checks and random access for essential enumeration."""
from __future__ import annotations

import itertools
import random
import weakref

import pytest

from zipper_algebra.zipping_algorithm.enumeration import alphabet
from zipper_algebra.zipping_algorithm.essential_enumeration import EssentialWords


def orbit(word):
    seen, queue = {word}, [word]
    while queue:
        current = queue.pop()
        for i in range(len(current) - 1):
            if abs(abs(current[i]) - abs(current[i + 1])) >= 2:
                other = current[:i] + (current[i + 1], current[i]) + current[i + 2:]
                if other not in seen:
                    seen.add(other)
                    queue.append(other)
    return seen


@pytest.mark.parametrize("n,length", [(1, 0), (1, 2), (2, 0), (2, 4), (3, 4), (4, 4), (4, 6), (5, 4)])
@pytest.mark.parametrize("commute", [False, True])
def test_against_literal_words_and_independent_swap_orbits(n, length, commute):
    raw = {w for w in itertools.product(alphabet(n), repeat=length)
           if {abs(a) for a in w} == set(range(1, n))}
    expected = sorted({min(orbit(w)) for w in raw} if commute else raw)
    layer = EssentialWords(n, length, commute=commute)
    assert layer.count == len(expected)
    assert list(layer.iter_words()) == expected
    for i, word in enumerate(expected):
        assert layer.word_at(i) == word
        assert layer.ordinal(word) == i
    for start in range(0, len(expected) + 1, max(1, len(expected) // 7)):
        stop = min(len(expected), start + 9)
        assert list(layer.iter_words(start, stop)) == expected[start:stop]


def test_known_full_support_totals_and_signed_inverse_retention():
    assert sum(EssentialWords(4, length, commute=True).count for length in range(11)) == 19606864
    assert sum(EssentialWords(4, length, commute=False).count for length in range(11)) == 68371248
    for commute in (False, True):
        layer = EssentialWords(2, 2, commute=commute)
        assert layer.count == 4
        assert (1, -1) in list(layer.iter_words())
        assert (-1, 1) in list(layer.iter_words())
    # Adjacent braid relations are not imposed.
    layer = EssentialWords(3, 3, commute=True)
    assert layer.ordinal((1, 2, 1)) != layer.ordinal((2, 1, 2))


def test_nonlocal_commutation_and_sign_order():
    layer = EssentialWords(4, 4, commute=True)
    with pytest.raises(ValueError, match="representative"):
        layer.ordinal((3, 3, 1, 2))
    assert layer.word_at(layer.ordinal((-3, 1, 3, 2))) == (-3, 1, 3, 2)


@pytest.mark.parametrize("args,kwargs", [((0, 2), {}), ((True, 2), {}), ((2, -1), {}),
                                         ((2, True), {}), ((2, 2), {"commute": 1})])
def test_invalid_layer(args, kwargs):
    with pytest.raises(ValueError):
        EssentialWords(*args, **kwargs)


def test_invalid_ordinals_and_words():
    layer = EssentialWords(4, 4, commute=True)
    for ordinal in (-1, True, 1.5):
        with pytest.raises(ValueError):
            layer.word_at(ordinal)
    with pytest.raises(IndexError):
        layer.word_at(layer.count)
    for word in ((1,), (1, 1, 1, 1), (0, 1, 2, 3), (True, 1, 2, 3)):
        with pytest.raises(ValueError):
            layer.ordinal(word)
    for start, stop in ((-1, 0), (2, 1), (0, layer.count + 1), (True, 2)):
        with pytest.raises(ValueError):
            list(layer.iter_words(start, stop))
    assert list(layer.iter_words(layer.count, layer.count)) == []


def independent_full_support_trace_count(n, length):
    """Inclusion-exclusion of clique-polynomial inverses, no prefix DFA."""
    answer = 0
    for mask in range(1 << (n - 1)):
        coefficients = [1] + [0] * length
        for clique in range(1, 1 << (n - 1)):
            if clique & ~mask or clique & (clique << 1):
                continue
            size = clique.bit_count()
            if size <= length:
                coefficients[size] += (-2) ** size
        growth = [1]
        for degree in range(1, length + 1):
            growth.append(-sum(coefficients[k] * growth[degree - k]
                               for k in range(1, degree + 1)))
        answer += (-1) ** (n - 1 - mask.bit_count()) * growth[length]
    return answer


@pytest.mark.parametrize("n", range(2, 9))
def test_large_counts_against_independent_clique_inversion(n):
    for length in (n - 2, n - 1, 10):
        layer = EssentialWords(n, length)
        assert layer.count == independent_full_support_trace_count(n, length)
        if layer.count:
            for ordinal in (0, layer.count // 2, layer.count - 1):
                assert layer.ordinal(layer.word_at(ordinal)) == ordinal


def test_large_samples_against_dependency_poset_and_intervals():
    from zipper_algebra.zipping_algorithm.word_coverage import trace_permutation

    layer = EssentialWords(8, 10)
    rng = random.Random(1707)
    ordinals = {0, layer.count - 2, layer.count // 2}
    ordinals.update(rng.randrange(layer.count - 1) for _ in range(300))
    for ordinal in sorted(ordinals):
        word, following = layer.word_at(ordinal), layer.word_at(ordinal + 1)
        assert word < following
        for representative in (word, following):
            permutation = trace_permutation(representative)
            assert tuple(representative[i] for i in permutation) == representative
        assert list(layer.iter_words(ordinal, ordinal + 2)) == [word, following]
    for start in (layer.count // 3, layer.count // 2, layer.count - 73):
        words = list(layer.iter_words(start, start + 73))
        assert words == [layer.word_at(i) for i in range(start, start + 73)]
        assert words == sorted(set(words))


def test_count_closures_do_not_keep_layer_instance_alive():
    layer = EssentialWords(8, 10)
    assert layer.count > 0
    reference = weakref.ref(layer)
    # Holding a counter must not retain the layer and its other counter.
    counter = layer._trace_count
    del layer
    assert reference() is None
    assert counter(0, 0, 0) == 1
