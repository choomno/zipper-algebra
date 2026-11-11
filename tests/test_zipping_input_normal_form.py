"""Current-procedure request witnesses and independent output comparisons."""
from copy import deepcopy
from itertools import product
import json

import pytest

from zipper_algebra.zipping_algorithm.input_normal_form import (
    POLICY, canonical_input_json, check_input_normal_form, input_normal_form,
    leaf_requests, reconstruction_key, reconstruction_plan,
)
from zipper_algebra.zipping_algorithm.collected import CollectedSession


def test_split_translation_and_idle_maps():
    form = input_normal_form((2, -5, 2), n=7)
    assert check_input_normal_form(form)
    assert form['root']['kind'] == 'split'
    assert [edge['offset'] for edge in form['root']['children']] == [6, 4, 3, 1, 0]
    assert [(r['n'], r['word']) for r in leaf_requests(form)] == [(2, [1, 1])]
    leaf = next(a for a in reconstruction_plan(form)['atoms'] if a['kind'] == 'leaf')
    assert leaf['bottom_labels'] == [2, 3]
    assert leaf['top_labels'] == [2, 3]


def test_nonidle_peel_and_same_reconstruction():
    a = input_normal_form((1, -3, 1, 3, 1, -3, 2, 1), n=4)
    b = input_normal_form((-3, 1, 1, 3, 1, -3, 2, 1), n=4)
    assert a['root']['kind'] == 'peel'
    assert a['root']['top_port'] == 3
    assert a['root']['child']['top_labels'] == [1, 2, 4]
    assert a['root']['child']['word'] == [1, 1, 1, 2, 1]
    assert reconstruction_key(a) == reconstruction_key(b)
    assert canonical_input_json(a) != canonical_input_json(b)
    assert reconstruction_plan(a)['atoms'][0] == {'kind': 'singleton', 'bottom': 4, 'top': 3}


def test_identity_and_leaf_deduplication():
    empty = input_normal_form([], n=2)
    assert reconstruction_key(empty) == reconstruction_key(input_normal_form([1, -1], n=2))
    assert reconstruction_key(empty) == reconstruction_key(input_normal_form([-1, 1], n=2))
    form = input_normal_form([1, 3], n=4)
    assert len(leaf_requests(form)) == 1
    assert len(reconstruction_plan(form)['atoms']) == 2
    assert leaf_requests(empty) == []
    assert reconstruction_plan(input_normal_form([], n=3))['atoms'] == [
        {'kind': 'singleton', 'bottom': 3, 'top': 3},
        {'kind': 'singleton', 'bottom': 2, 'top': 2},
        {'kind': 'singleton', 'bottom': 1, 'top': 1}]


def test_no_unproved_higher_rank_identification():
    a = input_normal_form([1, -3, 2, -3], n=4)
    b = input_normal_form([-3, 1, 2, -3], n=4)
    assert a['root']['kind'] == b['root']['kind'] == 'leaf'
    assert reconstruction_key(a) != reconstruction_key(b)
    assert input_normal_form([2, 1, -1], n=3)['root']['word'] == [2, 1, -1]


def test_rank_two_exponent_and_full_output():
    session = CollectedSession()
    grouped = {}
    for length in range(6):
        for word in product((-1, 1), repeat=length):
            form = input_normal_form(word, n=2)
            assert check_input_normal_form(form)
            key = reconstruction_key(form)
            output = sorted(json.dumps(x, sort_keys=True) for x in session.reduce(word, n=2)['output'])
            if key in grouped:
                assert grouped[key] == output
            grouped[key] = output
            canonical = (1 if sum(word) > 0 else -1,) * abs(sum(word))
            assert key == reconstruction_key(input_normal_form(canonical, n=2))
    assert len(grouped) == 11


@pytest.mark.parametrize('n,cutoff', [(1, 0), (2, 4), (3, 3), (4, 3)])
def test_exhaustive_witnesses_and_roundtrip(n, cutoff):
    alphabet = tuple(range(1-n, 0)) + tuple(range(1, n))
    for length in range(cutoff + 1):
        for word in product(alphabet, repeat=length):
            form = input_normal_form(word, n=n)
            assert check_input_normal_form(form)
            assert json.loads(canonical_input_json(form)) == form
            assert form['root']['word'] == list(word)
            for leaf in leaf_requests(form):
                assert {abs(a) for a in leaf['word']} == set(range(1, leaf['n']))
                assert input_normal_form(leaf['word'], n=leaf['n'])['root']['kind'] == 'leaf'


def test_constructor_and_plan_do_not_alias_inputs():
    word = [1, 1]
    form = input_normal_form(word, n=2)
    word.append(1)
    request = leaf_requests(form)[0]
    request['word'].append(1)
    plan = reconstruction_plan(form)
    plan['atoms'][0]['top_labels'].append(3)
    assert form['root']['word'] == [1, 1]
    assert check_input_normal_form(form)


@pytest.mark.parametrize('word,n', [([True], 2), ([0], 2), ([2], 2), ([1], True), ([], 0), ('1', 2)])
def test_invalid_inputs(word, n):
    with pytest.raises(ValueError):
        input_normal_form(word, n=n)


def test_tampered_witnesses_rejected():
    forms = [input_normal_form([1, 3], n=4), input_normal_form([-1], n=2),
             input_normal_form([1], n=2), input_normal_form([1, -1], n=2)]
    mutations = [
        lambda f: f.update(policy=POLICY + '-other'),
        lambda f: f.update(schema='other'),
        lambda f: f.update(extra=1),
        lambda f: f['root'].update(rank=True),
        lambda f: f['root'].update(top_labels=[True] + f['root']['top_labels'][1:]),
        lambda f: f['root'].update(word=[0]),
        lambda f: f['root'].update(extra='unknown'),
    ]
    for original in forms:
        for mutate in mutations:
            form = deepcopy(original)
            mutate(form)
            assert not check_input_normal_form(form)
            with pytest.raises(ValueError):
                canonical_input_json(form)
    form = deepcopy(forms[0]); form['root']['children'].reverse()
    assert not check_input_normal_form(form)
    form = deepcopy(forms[0]); form['root']['children'][0]['positions'] = [0]
    assert not check_input_normal_form(form)
    form = deepcopy(forms[1]); form['root']['top_port'] = 2
    assert not check_input_normal_form(form)
    form = deepcopy(forms[3]); form['root']['exponent'] = True
    assert not check_input_normal_form(form)
    assert not check_input_normal_form(None)
    assert not check_input_normal_form({'schema': 'other'})


def test_strict_json_containers_and_complete_split_witness():
    from collections import UserDict

    base = input_normal_form([1, 3], n=4)
    malformed = []
    form = deepcopy(base)
    form['root'] = UserDict(form['root'])
    malformed.append(form)
    form = deepcopy(base)
    form['root']['word'] = tuple(form['root']['word'])
    malformed.append(form)
    form = deepcopy(base)
    form['root']['children'] = tuple(form['root']['children'])
    malformed.append(form)
    form = deepcopy(base)
    form['root']['children'].pop()
    malformed.append(form)
    form = deepcopy(base)
    form['root']['children'][0]['node']['word'] = [1, 1]
    malformed.append(form)
    form = deepcopy(base)
    form['root']['children'][0]['node']['top_labels'] = [1, 2]
    malformed.append(form)
    for form in malformed:
        assert not check_input_normal_form(form)
        with pytest.raises(ValueError):
            reconstruction_plan(form)


def test_constructor_fails_closed_on_incorrect_peeling_witness(monkeypatch):
    import importlib

    module = importlib.import_module('zipper_algebra.zipping_algorithm.input_normal_form')
    # The real negative crossing sends bottom2 to top1, not top2. The verifier
    # has its own single-position scan, so a corrupted producer is detected.
    monkeypatch.setattr(module, '_peel', lambda word, n: (True, [], 2))
    with pytest.raises(AssertionError, match='invalid generated'):
        input_normal_form([-1], n=2)
