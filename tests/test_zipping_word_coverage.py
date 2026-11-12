"""Input coverage is checked independently of the zipper transition engine."""
from dataclasses import replace
import itertools

import pytest

from zipper_algebra.zipping_algorithm.word_coverage import (
    Block, Coverage, adjacent_swaps, check_commutation, check_coverage, decompose,
    independent, trace_permutation,
)


def orbit(word):
    reached = {word}
    todo = [word]
    while todo:
        current = todo.pop()
        for i in range(len(current) - 1):
            if independent(current[i], current[i + 1]):
                other = current[:i] + (current[i + 1], current[i]) + current[i + 2:]
                if other not in reached:
                    reached.add(other)
                    todo.append(other)
    return reached


def test_all_small_words_against_explicit_commutation_orbits():
    for length in range(5):
        for word in itertools.product((-3, -2, -1, 1, 2, 3), repeat=length):
            permutation = trace_permutation(word)
            assert check_commutation(word, permutation)
            assert tuple(word[i] for i in permutation) == min(orbit(word))
            for commute in (False, True):
                certificate = decompose(word, n=4, commute=commute)
                assert check_coverage(certificate)
                expanded = list(word)
                for i in adjacent_swaps(certificate):
                    assert independent(expanded[i], expanded[i + 1])
                    expanded[i], expanded[i + 1] = expanded[i + 1], expanded[i]
                assert expanded == [word[i] for i in certificate.permutation]
                assert sum(map(lambda b: len(b.word), certificate.blocks)) == length


def occurrence_permutation(source, target):
    """Match repeated occurrences stably, without using the production DAG."""
    remaining = list(enumerate(source))
    result = []
    for letter in target:
        at = next(i for i, (_, value) in enumerate(remaining) if value == letter)
        occurrence, _ = remaining.pop(at)
        result.append(occurrence)
    return tuple(result)


def assert_orbit_checker_selects_only_minimum(word, n):
    candidates = orbit(word)
    least = min(candidates)
    for candidate in candidates:
        permutation = occurrence_permutation(word, candidate)
        assert check_commutation(word, permutation)
        certificate = Coverage(n, word, True, (Block(0, n, permutation, candidate),))
        assert check_coverage(certificate) == (candidate == least)


def test_checker_rejects_every_nonminimal_small_full_support_orbit_word():
    for length in (3, 4):
        for word in itertools.product((-3, -2, -1, 1, 2, 3), repeat=length):
            if {abs(x) for x in word} == {1, 2, 3}:
                assert_orbit_checker_selects_only_minimum(word, 4)


@pytest.mark.parametrize("word", [
    (4, 3, 1, 2), (-1, -2, -4, -3), (4, -3, -1, 2),
    (5, 4, 2, 3, 1), (-1, -2, -4, -3, -5), (5, -4, 2, -3, 1),
])
def test_checker_dependency_patterns_at_larger_ranks(word):
    # The middle of the first three letters depends on the outer letter,
    # while the final letter commutes with both; signs change lexical order.
    assert not independent(word[0], word[1])
    assert independent(word[0], word[2])
    assert independent(word[1], word[2])
    assert_orbit_checker_selects_only_minimum(word, max(map(abs, word)) + 1)


def test_disconnected_blocks_translation_and_inverse_pairs():
    certificate = decompose((6, -2, 5, 2, -6, 9), n=11, commute=True)
    assert certificate.blocks == (
        Block(1, 2, (1, 3), (-1, 1)),
        Block(4, 3, (0, 2, 4), (2, 1, -2)),
        Block(8, 2, (5,), (1,)),
    )
    assert check_coverage(certificate)
    assert decompose((), n=8).blocks == ()


def test_signed_letters_on_same_generator_do_not_commute():
    assert not check_commutation((1, -1), (1, 0))
    assert not check_commutation((1, 2), (1, 0))
    assert check_commutation((1, -3), (1, 0))


def test_tampered_certificates():
    certificate = decompose((3, 1, 2), n=4, commute=True)
    block = certificate.blocks[0]
    assert not check_coverage(replace(certificate, blocks=(replace(block, word=(-1, 3, 2)),)))
    assert not check_coverage(replace(certificate, blocks=(replace(block, positions=(0, 1, 1)),)))
    assert not check_coverage(replace(certificate, blocks=(replace(block, offset=1),)))
    # Valid commuting order, but not the required least representative.
    assert not check_coverage(Coverage(4, (3, 1, 2), True, (Block(0, 4, (0, 1, 2), (3, 1, 2)),)))
    assert not check_coverage(replace(certificate, commute=False))
    # Splitting adjacent supports would incorrectly factor a connected block.
    assert not check_coverage(Coverage(3, (1, 2), False,
                                     (Block(0, 2, (0,), (1,)), Block(1, 2, (1,), (1,)))))
    with pytest.raises(ValueError, match="invalid coverage"):
        adjacent_swaps(replace(certificate, n=1))


@pytest.mark.parametrize("word,n", [((0,), 3), ((3,), 3), ((True,), 3), ((1,), True), ((), 0)])
def test_invalid_inputs(word, n):
    with pytest.raises(ValueError):
        decompose(word, n=n)


def test_invalid_commutation_inputs():
    assert not check_commutation((1,), (True,))
    assert not check_commutation((1,), (1,))
    assert not check_commutation((0,), (0,))
    with pytest.raises(ValueError):
        trace_permutation((0,))
    with pytest.raises(ValueError):
        decompose((), n=1, commute=1)


def test_checker_rejects_malformed_shapes():
    certificate = decompose((1,), n=2)
    assert not check_coverage(replace(certificate, blocks=(None,)))
    assert not check_coverage(replace(certificate, blocks=(replace(certificate.blocks[0], positions=None),)))
    assert not check_coverage(replace(certificate, word=None))
    assert not check_coverage(replace(certificate, commute=1))
