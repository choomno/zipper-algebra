"""Exact interchange identities and rejection boundaries for terminal outputs."""
from copy import deepcopy
import json

import pytest

from zipper_algebra.zipping_algorithm.collected import CollectedSession
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.output_normal_form import (
    canonical_bytes, canonical_output, engine_output, from_certificate, output_fingerprint,
)


@pytest.fixture(scope="module")
def runs():
    session = CollectedSession()
    return [session.reduce(word, n=n) for n, word in ((1, []), (2, [1]), (3, [1, 2, 1]))]


def test_current_outputs_and_json_roundtrip(runs):
    for run in runs:
        assert run['status'] == 'complete'
        form = from_certificate(run)
        assert form['completion'] == 'complete'
        assert form == json.loads(canonical_bytes(form))
        assert canonical_bytes(form) == canonical_bytes(json.loads(canonical_bytes(form)))
        assert len(output_fingerprint(form)) == 64
        assert all(type(x) is str for term in form['terms'] for pair in term['coefficient'] for x in pair)


def test_term_order_key_order_and_exact_collection(runs):
    run = runs[1]
    original = deepcopy(run['output'])
    messy = list(reversed(deepcopy(original)))
    term = deepcopy(messy[0])
    term['coefficient'] = [[2, 3], [2, -3], [0, 5]]
    opposite = deepcopy(term)
    opposite['coefficient'] = [[0, -5]]
    messy.extend([term, opposite])
    messy = [dict(reversed(list(t.items()))) for t in messy]
    assert canonical_output(messy, n=2, complete=True) == from_certificate(run)
    assert run['output'] == original


def test_wide_coefficients_without_global_integer_conversion_changes(runs):
    import sys
    before = sys.get_int_max_str_digits()
    term = deepcopy(runs[0]['output'][0])
    huge = 10**5000 + 7
    term['coefficient'] = [[huge, huge], [-huge, -huge]]
    form = canonical_output([term], n=1, complete=True)
    encoded = canonical_bytes(form)
    assert len(encoded) > 20000
    assert canonical_bytes(json.loads(encoded)) == encoded
    assert canonical_output(engine_output(form), n=1, complete=True) == form
    assert sys.get_int_max_str_digits() == before


def test_complete_partial_rank_and_extension_identity(runs):
    full = from_certificate(runs[0])
    partial = canonical_output(runs[0]['output'], n=1, complete=False)
    assert output_fingerprint(full) != output_fingerprint(partial)
    term = deepcopy(runs[0]['output'][0])
    term['mathematical_extension'] = {'flag': True, 'data': (1, 'x', None)}
    extended = canonical_output([term], n=1, complete=True)
    assert extended['terms'][0]['mathematical_extension']['data'] == [1, 'x', None]
    assert output_fingerprint(extended) != output_fingerprint(full)
    term['mathematical_extension']['flag'] = 1
    assert output_fingerprint(canonical_output([term], n=1, complete=True)) != output_fingerprint(extended)
    assert output_fingerprint(canonical_output([], n=1, complete=True)) != output_fingerprint(canonical_output([], n=2, complete=True))


def test_partial_certificate_is_not_completion_or_checkpoint(runs):
    run = CollectedSession().reduce([1], n=2, limits=Limits(max_steps=1))
    assert from_certificate(run)['completion'] == 'partial'
    run['status'] = 'control_obligation'
    assert from_certificate(run)['completion'] == 'partial'
    run['status'] = 'complete'
    with pytest.raises(ValueError, match='pending'):
        from_certificate(run)
    for change in ({'schema': 'unknown'}, {'coefficient_ring': 'Q'}, {'status': 'failed'}, {'pending': None}):
        bad = dict(runs[0], **change)
        with pytest.raises(ValueError):
            from_certificate(bad)


@pytest.mark.parametrize('coefficient', [None, [1], [[1]], [[True, 1]], [[0, False]],
                                        [[0, 1.5]], [['01', '2']], [['-0', '1']],
                                        [['', '2']], [['+', '1']], [['١', '2']]])
def test_bad_coefficients_rejected(runs, coefficient):
    term = deepcopy(runs[0]['output'][0])
    term['coefficient'] = coefficient
    with pytest.raises(ValueError):
        canonical_output([term], n=1, complete=True)


@pytest.mark.parametrize('mutation', [
    lambda t: t.pop('components'),
    lambda t: t.update(word=None),
    lambda t: t.update(word=[0]),
    lambda t: t.update(word=[2]),
    lambda t: t.update(word=[True]),
    lambda t: t.update(word=[['slider', 2, 1, 0, False, False]]),
    lambda t: t.update(components=None),
    lambda t: t.update(components=[{}]),
    lambda t: t['components'][0].update(rank=1),
    lambda t: t['components'][0].update(size=True),
    lambda t: t['components'][0].update(component_word=[]),
    lambda t: t['components'][0].update(top_route=[['slider', 1, 1, False, False, False]]),
    lambda t: t['components'][0]['record'].pop('orientation'),
    lambda t: t['components'][0]['record'].update(bottom_ports=[2, 1]),
    lambda t: t['components'][0]['record'].update(original_top_ports=[1]),
    lambda t: t['components'][0]['record'].update(twist_powers=[]),
    lambda t: t['components'][0]['record'].update(orientation='barred'),
    lambda t: t['components'][0]['record'].update(orientation='mirrored'),
    lambda t: t['components'][0]['record'].update(canonical_parities=None),
    lambda t: t['components'][0]['record'].update(canonical_parities={'upper': [True], 'lower_inverse': []}),
    lambda t: t.update(extra=1.2),
    lambda t: t.update(extra={1: 'x'}),
])
def test_bad_full_term_structures_rejected(runs, mutation):
    term = deepcopy(next(t for t in runs[1]['output'] if t['components'][0]['size'] == 2))
    mutation(term)
    with pytest.raises(ValueError):
        canonical_output([term], n=2, complete=True)


def test_partition_reconstruction_hierarchy_and_flag_preservation(runs):
    singleton_term = deepcopy(next(t for t in runs[1]['output'] if len(t['components']) == 2))
    wrong_ports = deepcopy(singleton_term)
    record = wrong_ports['components'][0]['record']
    record['original_bottom_ports'] = [3 - record['original_bottom_ports'][0]]
    with pytest.raises(ValueError, match='inherited endpoint'):
        canonical_output([wrong_ports], n=2, complete=True)
    singleton_term['components'].pop()
    with pytest.raises(ValueError, match='partition'):
        canonical_output([singleton_term], n=2, complete=True)
    term = deepcopy(next(t for t in runs[2]['output'] if t['components'][0]['size'] == 3))
    term['components'][0]['record']['input_hierarchy'] = [1, 1, 3]
    with pytest.raises(ValueError, match='hierarchy'):
        canonical_output([term], n=3, complete=True)
    term = deepcopy(next(t for t in runs[1]['output'] if t['components'][0]['size'] == 2))
    baseline = canonical_output([term], n=2, complete=True)
    term['word'][0][3] = True
    with pytest.raises(ValueError, match='reconstruction'):
        canonical_output([term], n=2, complete=True)
    term['components'][0]['component_word'][0][3] = True
    term['components'][0]['record']['orientation'] = 'mirrored'
    assert output_fingerprint(canonical_output([term], n=2, complete=True)) != output_fingerprint(baseline)


def test_normal_form_validation_and_empty_sum(runs):
    form = canonical_output([], n=1, complete=True)
    assert form['terms'] == []
    for change in ({'schema': 'wrong'}, {'completion': 'unknown'}, {'extra': True}, {'n': True}, {'coefficient_ring': 'Q'}):
        with pytest.raises(ValueError):
            canonical_bytes(dict(form, **change))
    for n, complete, output in ((0, True, []), (True, True, []), (1, 1, []), (1, True, {})):
        with pytest.raises(ValueError):
            canonical_output(output, n=n, complete=complete)
    with pytest.raises(ValueError):
        canonical_output([None], n=1, complete=True)
    zero = deepcopy(runs[0]['output'][0])
    zero['coefficient'] = [[0, 1], [0, -1]]
    assert canonical_output([zero], n=1, complete=True) == form


@pytest.mark.parametrize("side", ["top", "bottom"])
def test_noncanonical_ambient_route_rejected_even_with_matching_word(runs, side):
    term = deepcopy(next(t for t in runs[1]['output'] if len(t['components']) == 2))
    # This inverse pair changes no braid and the full literal word still agrees
    # with its contexts, but it is not the canonical stable-partition route.
    if side == 'top':
        term['components'][0]['top_route'] = [1, -1] + term['components'][0]['top_route']
        term['word'] = [1, -1] + term['word']
    else:
        term['components'][0]['bottom_route'] += [1, -1]
        term['word'] += [1, -1]
    with pytest.raises(ValueError, match='canonical port routes'):
        canonical_output([term], n=2, complete=True)


def test_native_routed_examples_have_canonical_ambient_routes():
    from zipper_algebra.zipping_algorithm.native import NativeSession
    with NativeSession() as session:
        for n, word in ((2, [1]), (3, [2, 1, -2]), (3, [1, 2, 1]),
                        (4, [1, -3, 1, 3, 1, -3, 2, 1])):
            run = session.reduce(word, n=n)
            assert run['status'] == 'complete'
            form = from_certificate(run)
            assert canonical_bytes(form) == canonical_bytes(
                canonical_output(engine_output(form), n=n, complete=True))
