"""Independent braid identities and failures for explicitly marked ribbons."""
from dataclasses import replace
import random

import pytest

from zipper_algebra.zipping_algorithm.extraction import (
    NestedRibbon, RibbonBlock, RibbonError, check_extraction, check_nested,
    extract_slider, recognize_nested, reconstruct_nested, split_layered, check_split,
    normalize_canonical_slider, check_canonical_slider, BoundaryRoute,
    RibbonLimit, RibbonGathering, recognize_nested_strict,
)


def _inverse(word):
    return tuple(-x for x in reversed(word))


def _key(word, n):
    images = tuple((x,) for x in range(1, n + 1))
    for letter in word:
        at = abs(letter)
        generator = [(x,) for x in range(1, n + 1)]
        generator[at - 1:at + 1] = ([(at, at + 1, -at), (at,)] if letter > 0 else
                                          [(at + 1,), (-at - 1, at, at + 1)])
        new = []
        for image in generator:
            stack = []
            for x in image:
                for y in images[x - 1] if x > 0 else _inverse(images[-x - 1]):
                    if stack and stack[-1] == -y:
                        stack.pop()
                    else:
                        stack.append(y)
            new.append(tuple(stack))
        images = tuple(new)
    return images


def _cross(ports, labels, size, sign):
    """Independent geometric shuffle: move the singleton across its old child."""
    singleton = labels[size - 1]
    old = set(labels[:size - 1])
    at = ports.index(singleton)
    direction = -1 if at and ports[at - 1] in old else 1
    letters = []
    for _ in range(size - 1):
        nxt = at + direction
        letter = min(at, nxt) + 1
        letters.append(sign * letter)
        ports[at], ports[nxt] = ports[nxt], ports[at]
        at = nxt
    return letters


def test_interleaved_nested_crossings_reconstruct_exact_braids():
    rng = random.Random(5633)
    for n in range(2, 6):
        labels = tuple(range(1, n + 1))
        for _ in range(30):
            ports = [1]
            for leaf in labels[1:]:
                ports.insert(0 if rng.randrange(2) else len(ports), leaf)
            initial = tuple(ports)
            chron = []
            for _ in range(9):
                chron.extend(_cross(ports, labels, rng.randint(2, n), rng.choice((-1, 1))))
            word = tuple(reversed(chron))
            witness = recognize_nested(word, labels=labels, initial_ports=initial)
            check_nested(witness)
            assert witness.final_ports == tuple(ports)
            assert _key(word, n) == _key(witness.reconstructed_word, n)
            assert len(witness.reconstructed_word) <= len(word)


def test_thesis_lasso_and_half_twist_have_explicit_tree_witnesses():
    # Theta_3 = sigma_1 sigma_2^2 sigma_1 wraps the singleton 1 around child {2,3}.
    lasso = recognize_nested((1, 2, 2, 1), labels=(2, 3, 1), initial_ports=(1, 2, 3))
    assert lasso.crossing_counts == (0, 2)
    assert lasso.twist_powers == (0, 1)
    assert lasso.boundary_parities == (0,)
    # Delta_3 has one child half-turn and one crossing of that child with leaf 3.
    half = recognize_nested((1, 2, 1), labels=(1, 2, 3))
    assert half.crossing_counts == (1, 1)
    assert half.twist_powers == (1, 0)
    assert half.boundary_parities == (1,)
    for word in ((1, 2, 1), (-1, -2, -1), (1, 2, 1, 1, 2, 1)):
        witness = recognize_nested(word, labels=(1, 2, 3))
        check_nested(witness)
        assert _key(word, 3) == _key(witness.reconstructed_word, 3)


def test_negative_odd_counts_retain_the_boundary_half_crossing():
    witness = recognize_nested((-1, -2), labels=(1, 2, 3))
    assert witness.crossing_counts == (0, -1)
    assert witness.twist_powers == (0, -1)
    assert witness.boundary_parities == (1,)
    assert 2 * witness.twist_powers[1] + witness.boundary_parities[0] == -1


def test_empty_ribbon_and_arbitrary_original_integer_labels():
    witness = recognize_nested((), labels=(23,))
    check_nested(witness)
    assert witness.twist_powers == witness.boundary_parities == ()
    assert reconstruct_nested((23,), ()) == ()
    pair = recognize_nested((1, -1, 1), labels=(23, 8))
    assert pair.crossing_counts == (1,)
    assert pair.final_ports == (8, 23)


@pytest.mark.parametrize("word,labels,ports", [
    ((), (), None), ((), (1, 1), None), ((), (True,), None),
    ((), (1, 2), (1,)), ((), (1, 2), (1, 3)), ((), (1, 2), (1, True)),
    ((), (1, 2, 3), (1, 3, 2)), ((0,), (1, 2), None),
    ((2,), (1, 2), None), ((True,), (1, 2), None),
    ((2,), (1, 2, 3), None), ((2, -2), (1, 2, 3), None),
    ((-1, 2), (1, 2, 3), None),
])
def test_invalid_or_partial_nested_witnesses_are_refused(word, labels, ports):
    with pytest.raises(RibbonError):
        recognize_nested_strict(word, labels=labels, initial_ports=ports)
    if word != (2, -2):
        with pytest.raises(RibbonError):
            recognize_nested(word, labels=labels, initial_ports=ports)


def test_moving_nonadjacent_children_is_refused():
    # Tree is valid initially, but the explicit private helper also guards an
    # actual boundary obtained when an unrecognized external strand intervenes.
    from zipper_algebra.zipping_algorithm.extraction import _children
    with pytest.raises(RibbonError, match="not contiguous"):
        _children([1, 3, 2], (1, 2, 3), 3)
    with pytest.raises(RibbonError, match="not adjacent"):
        _children([1, 2, 4, 3], (1, 2, 3, 4), 3)


@pytest.mark.parametrize("counts", [(), (1,), (1, 2, 3), (True, 0)])
def test_invalid_node_counts(counts):
    with pytest.raises(RibbonError):
        reconstruct_nested((1, 2, 3), counts)


def test_checker_rejects_changes_to_all_mathematical_witness_parts():
    witness = recognize_nested((1, 2, 1), labels=(1, 2, 3))
    for changed in [
        replace(witness, crossing_counts=(2, 1)),
        replace(witness, final_ports=(1, 2, 3)),
        replace(witness, reconstructed_word=()),
        replace(witness, blocks=(RibbonBlock(3, 1, (1, 2)), RibbonBlock(2, 1, (1,)))),
        replace(witness, blocks=(RibbonBlock(2, -1, (1,)),) + witness.blocks[1:]),
        replace(witness, blocks=(RibbonBlock(2, 1, (0,)),) + witness.blocks[1:]),
        replace(witness, blocks=(RibbonBlock(1, 1, ()),)),
        replace(witness, blocks=(RibbonBlock(2, True, (1,)),)),
        replace(witness, blocks=(None,)),
        replace(witness, blocks=list(witness.blocks)),
    ]:
        with pytest.raises(RibbonError):
            check_nested(changed)
    with pytest.raises(RibbonError):
        check_nested(None)
    assert isinstance(witness, NestedRibbon)


def test_coupon_extraction_preserves_both_sides_and_orientation_flags():
    token = ("slider", 3, 1, True, False, False)
    original = (1, 2, 1, token, -1, -2)
    witness = extract_slider(original, upper_labels=(1, 2, 3), lower_labels=(1, 2, 3))
    check_extraction(witness)
    assert witness.slider == token
    assert witness.reconstructed_word.count(token) == 1
    assert _key(original[:3], 3) == _key(witness.upper.reconstructed_word, 3)
    assert _key(original[4:], 3) == _key(witness.lower.reconstructed_word, 3)
    with pytest.raises(RibbonError):
        check_extraction(replace(witness, reconstructed_word=(token,)))
    with pytest.raises(RibbonError):
        check_extraction(None)


@pytest.mark.parametrize("word,upper,lower", [
    ((), (1,), (1,)), ((1,), (1, 2), (1, 2)),
    ((("slider", 1, 1, False, False, False),) * 2, (1,), (1,)),
    ((("slider", 2, 2, False, False, False),), (1, 2), (1, 2)),
    ((("slider", True, 1, False, False, False),), (1,), (1,)),
    ((("slider", 1, True, False, False, False),), (1,), (1,)),
    ((("slider", 0, 1, False, False, False),), (1,), (1,)),
    ((("slider", 1, 1, 1, False, False),), (1,), (1,)),
    ((("slider", 1),), (1,), (1,)),
    (("slider",), (1,), (1,)),
    ((("slider", 2, 1, False, False, False),), (1,), (1, 2)),
])
def test_invalid_coupon_or_missing_marked_tree_refused(word, upper, lower):
    with pytest.raises(RibbonError):
        extract_slider(word, upper_labels=upper, lower_labels=lower)


def _slider_polynomial(size, flags):
    """Independent exact definition; no catalog or extractor code is used."""
    from collections import defaultdict
    if size == 1:
        terms = {(0, ()): 1}
    else:
        terms = {(1, (-1,)): 1, (0, ()): 1, (1, ()): -1, (0, (1,)): -1}
    for rank in range(3, size + 1):
        out = defaultdict(int)
        fat = tuple(range(1, rank))
        for (power, word), coefficient in terms.items():
            out[(power + rank - 1, tuple(-x for x in fat) + word)] += coefficient
            out[(power, fat + word)] -= coefficient
        terms = dict(out)
    mirrored, reversed_product, barred = flags
    transformed = defaultdict(int)
    for (power, word), coefficient in terms.items():
        if mirrored:
            word = tuple((1 if x > 0 else -1) * (size - abs(x)) for x in word)
        if reversed_product:
            word = tuple(reversed(word))
        if barred:
            word = tuple(-x for x in word)
        transformed[(power, word)] += coefficient
    return dict(transformed)


def _expanded_coupon_diagram(word, n):
    from collections import defaultdict
    at = next(k for k, value in enumerate(word) if type(value) is not int)
    token = word[at]
    result = defaultdict(int)
    for (power, braid), coefficient in _slider_polynomial(token[1], token[3:]).items():
        shifted = tuple((1 if x > 0 else -1) * (abs(x) + token[2] - 1) for x in braid)
        result[(power, _key(word[:at] + shifted + word[at + 1:], n))] += coefficient
    return {key: value for key, value in result.items() if value}


def test_full_slider_extraction_in_the_generic_braid_group_algebra():
    import itertools
    for n in range(1, 5):
        tree = tuple(range(1, n + 1))
        ports = list(tree)
        chronological = []
        if n > 1:
            for size, sign in [(n, 1), (2, 1), (n, -1), (2, -1), (n, 1)]:
                chronological.extend(_cross(ports, tree, size, sign))
        upper = tuple(reversed(chronological))
        lower = _inverse(upper)
        for flags in itertools.product((False, True), repeat=3):
            token = ("slider", n, 1, *flags)
            original = upper + (token,) + lower
            witness = extract_slider(original, upper_labels=tree, lower_labels=tree,
                                     lower_ports=tuple(ports))
            check_extraction(witness)
            assert _expanded_coupon_diagram(original, n) == _expanded_coupon_diagram(
                witness.reconstructed_word, n)


def test_layered_component_deletion_reconstructs_exact_generic_coupon_diagrams():
    from zipper_algebra.zipping_algorithm.extraction import _over_route
    rng = random.Random(8541)
    for n in range(2, 6):
        for size in range(2, min(n, 4) + 1):
            selected = tuple(sorted(rng.sample(range(1, n + 1), size)))
            bottom_route = _over_route(n, selected)
            for mirrored in (False, True):
                chronological = list(reversed(bottom_route))
                chronological.append(("slider", size, 1, mirrored, False, False))
                active = [True] * size + [False] * (n - size)
                for _ in range(10):
                    at = rng.randrange(n - 1)
                    if active[at] != active[at + 1]:
                        sign = 1 if active[at] else -1
                    else:
                        sign = rng.choice((-1, 1))
                    chronological.append(sign * (at + 1))
                    active[at], active[at + 1] = active[at + 1], active[at]
                word = tuple(reversed(chronological))
                witness = split_layered(word, n=n, bottom_ports=selected)
                check_split(witness)
                assert witness.top_ports == tuple(i for i, value in enumerate(active, 1) if value)
                assert _expanded_coupon_diagram(word, n) == _expanded_coupon_diagram(
                    witness.reconstructed_word, n)
                assert sum(type(x) is not int for x in witness.component_word) == 1
                assert all(1 <= abs(x) < n - size for x in witness.complement_word)


def test_singleton_split_and_reindexing():
    for word, selected in [((1,), (1,)), ((-1,), (2,)), ((), (1,)),
                           ((-1, -2, 2, 1), (1,))]:
        witness = split_layered(word, n=3, bottom_ports=selected)
        check_split(witness)
        assert _key(word, 3) == _key(witness.reconstructed_word, 3)
        assert not witness.component_word
    pure = split_layered((3, 1, -3), n=4, bottom_ports=(1,))
    assert pure.complement_word == (2, -2)
    with_coupon = split_layered((("slider", 1, 2, False, False, False),),
                                n=3, bottom_ports=(2,))
    check_split(with_coupon)
    assert with_coupon.component_word == (("slider", 1, 1, False, False, False),)


@pytest.mark.parametrize("word,n,selected", [
    ((), True, (1,)), ((), 0, (1,)), ((), 2, ()), ((), 2, (True,)),
    ((), 2, (0,)), ((), 2, (3,)), ((), 3, (2, 1)), ((), 3, (1, 1)),
    ((0,), 2, (1,)), ((2,), 2, (1,)), ((1,), 2, (2,)), ((-1,), 2, (1,)),
    ((("slider", 2, 1, False, False, False),), 3, (1, 3)),
    ((("slider", 2, 2, False, False, False),), 3, (1, 2)),
    ((("slider", 1, True, False, False, False),), 1, (1,)),
    ((("slider", 1, 1, 0, False, False),), 1, (1,)),
    ((("slider", 1, 1, False, False, False),) * 2, 1, (1,)),
    ((), 2, (1, 2)), (("bad",), 1, (1,)),
])
def test_split_refuses_invalid_or_nonlayered_components(word, n, selected):
    with pytest.raises(RibbonError):
        split_layered(word, n=n, bottom_ports=selected)


def test_split_checker_binds_every_output_and_source_field():
    witness = split_layered((1,), n=3, bottom_ports=(1,))
    for changed in [replace(witness, top_ports=(1,)), replace(witness, component_word=(1,)),
                    replace(witness, complement_word=(1,)), replace(witness, bottom_route=(1,)),
                    replace(witness, top_route=()), replace(witness, reconstructed_word=())]:
        with pytest.raises(RibbonError):
            check_split(changed)
    with pytest.raises(RibbonError):
        check_split(None)


def test_noninitial_slider_anchor_has_exact_over_route_reconstruction():
    for mirrored in (False, True):
        token = ("slider", 3, 2, mirrored, False, False)
        original = (2, token, -3)
        witness = split_layered(original, n=5, bottom_ports=(2, 3, 4))
        check_split(witness)
        assert witness.component_word == (1, ("slider", 3, 1, mirrored, False, False), -2)
        assert witness.complement_word == ()
        assert _expanded_coupon_diagram(original, 5) == _expanded_coupon_diagram(
            witness.reconstructed_word, 5)


def test_two_sided_twists_collect_through_the_generic_coupon_with_boundary_flags():
    rng = random.Random(59113)
    for n in range(1, 5):
        ports = tuple(range(1, n + 1))
        for mirrored in (False, True):
            input_tree = tuple(reversed(ports)) if mirrored else ports
            output_tree = ports if mirrored else tuple(reversed(ports))
            for _ in range(8):
                sides = []
                for tree in (output_tree, input_tree):
                    current = list(ports)
                    chronological = []
                    if n > 1:
                        for _ in range(7):
                            chronological.extend(_cross(current, tree, rng.randint(2, n),
                                                        rng.choice((-1, 1))))
                    sides.append(tuple(reversed(chronological)))
                original = sides[0] + (("slider", n, 1, mirrored, False, False),) + _inverse(sides[1])
                witness = normalize_canonical_slider(original)
                check_canonical_slider(witness)
                assert witness.word == original
                assert witness.slider == ("slider", n, 1, False, False, False)
                assert witness.upper.labels == tuple(reversed(ports))
                assert witness.lower_inverse.labels == ports
                assert _expanded_coupon_diagram(original, n) == _expanded_coupon_diagram(
                    witness.normalized_word, n)
                assert len(witness.twist_powers) == n - 1
                route_bound = max(0, n * (n - 1) // 2 - 1)
                assert len(witness.top_route) <= route_bound
                assert len(witness.bottom_route) <= route_bound


def test_base_and_full_lasso_twists_have_independently_specified_collected_powers():
    normal = normalize_canonical_slider((("slider", 3, 1, False, False, False), 1))
    assert normal.twist_powers == (1, 0)
    assert normal.normalized_word == (2, ("slider", 3, 1, False, False, False))
    mirrored = normalize_canonical_slider((("slider", 3, 1, True, False, False), 2))
    assert mirrored.twist_powers == (1, 0)
    assert mirrored.normalized_word == (
        2, 1, 2, ("slider", 3, 1, False, False, False), -2, -1)
    full = normalize_canonical_slider((("slider", 3, 1, False, False, False), 2, 1, 1, 2))
    assert full.twist_powers == (0, 1)
    assert full.normalized_word == (1, 2, 2, 1, ("slider", 3, 1, False, False, False))


@pytest.mark.parametrize("word", [(), (1,), (("slider", 2, 1, False, False, False),) * 2,
    (("slider", 2, 1, False, True, False),), (("slider", 2, 1, False, False, True),),
    (("slider", 2, 2, False, False, False),), (("slider", True, 1, False, False, False),),
    (("slider", 1, True, False, False, False),), (("bad",),),
    (("slider", 3, 1, False, False, False), 2),
])
def test_canonical_slider_refuses_nonstandard_coupon_or_missing_ribbon(word):
    with pytest.raises(RibbonError):
        normalize_canonical_slider(word)


def test_canonical_slider_checker_binds_collected_vector_and_routes():
    witness = normalize_canonical_slider((("slider", 3, 1, False, False, False), 1))
    for changed in [replace(witness, twist_powers=(0, 0)), replace(witness, top_route=(1,)),
                    replace(witness, bottom_route=(1,)), replace(witness, normalized_word=())]:
        with pytest.raises(RibbonError):
            check_canonical_slider(changed)
    with pytest.raises(RibbonError):
        check_canonical_slider(None)


def test_inherited_terminal_friend_route_is_retained_outside_the_body():
    # Actual growth term from the B3 input sigma_1 sigma_2.
    token = ("slider", 3, 1, True, False, False)
    original = (2, -1, token, 2)
    route = BoundaryRoute(3, (2,))
    witness = normalize_canonical_slider(original, top_routes=(route,))
    check_canonical_slider(witness)
    assert witness.twist_powers == (0, 0)
    assert witness.normalized_word == (
        2, 2, 1, ("slider", 3, 1, False, False, False), -2, -1)
    assert witness.inherited_top == (route,)
    assert _expanded_coupon_diagram(original, 3) == _expanded_coupon_diagram(witness.normalized_word, 3)
    with pytest.raises(RibbonError):
        normalize_canonical_slider(original)


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_final_unflype_transports_every_descendant_and_keeps_default_coupon(n):
    original = (("slider", n, 1, True, False, False),)
    witness = normalize_canonical_slider(original)
    check_canonical_slider(witness)
    assert witness.slider == ("slider", n, 1, False, False, False)
    assert witness.upper.crossing_counts == (1,) * (n - 1)
    assert witness.lower_inverse.crossing_counts == (1,) * (n - 1)
    assert witness.twist_powers == (0,) * (n - 1)
    assert _expanded_coupon_diagram(original, n) == _expanded_coupon_diagram(
        witness.normalized_word, n)
    # The base crossing is absorbed into the twist vector; only higher-node
    # parity routes remain, even when every level changes its boundary order.
    assert len(witness.top_route) == len(witness.bottom_route) == n * (n - 1) // 2 - 1
    with pytest.raises(RibbonError):
        check_canonical_slider(replace(witness, slider=original[0]))


def test_finite_inherited_routes_on_both_sides_need_no_commuting_identity():
    token = ("slider", 4, 1, False, False, False)
    top = (BoundaryRoute(3, (-2, -1)), BoundaryRoute(4, (3,)))
    bottom = (BoundaryRoute(4, (-3, -2)), BoundaryRoute(2, (1,)))
    top_word = tuple(x for factor in top for x in factor.word)
    bottom_word = tuple(x for factor in bottom for x in factor.word)
    original = top_word + (token, 1) + bottom_word
    witness = normalize_canonical_slider(original, top_routes=top, bottom_routes=bottom)
    check_canonical_slider(witness)
    assert witness.normalized_word == top_word + (3, token) + bottom_word
    assert _expanded_coupon_diagram(original, 4) == _expanded_coupon_diagram(witness.normalized_word, 4)
    with pytest.raises(RibbonError):
        check_canonical_slider(replace(witness, inherited_top=()))
    empty = normalize_canonical_slider((token,), top_routes=(BoundaryRoute(4, ()),))
    check_canonical_slider(empty)
    assert empty.normalized_word == (token,)


@pytest.mark.parametrize("top,bottom", [
    ((None,), ()), ((BoundaryRoute(True, ()),), ()), ((BoundaryRoute(1, ()),), ()),
    ((BoundaryRoute(4, ()),), ()), ((BoundaryRoute(2, []),), ()),
    ((BoundaryRoute(2, ()), BoundaryRoute(2, ())), ()),
    ((BoundaryRoute(2, (1, 1)),), ()), ((BoundaryRoute(2, (0,)),), ()),
    ((BoundaryRoute(2, (2,)),), ()), ((), (BoundaryRoute(2, (1,)),)),
    ((BoundaryRoute(2, (1,)),), ()),
])
def test_invalid_or_nonliteral_inherited_boundaries_are_refused(top, bottom):
    token = ("slider", 3, 1, False, False, False)
    with pytest.raises(RibbonError):
        normalize_canonical_slider((token,), top_routes=top, bottom_routes=bottom)


def test_b5_generated_word_has_exact_semantic_nested_gathering():
    a = (1, -4, -3, -2, -1, -2, -3, -2, -4, -3, -2)
    b = (1, 2, 3, 1, 2, 1, 1, 2, 3, 3, 2, 1, -3, 2, 3)
    def phi(word):
        return tuple((1 if x > 0 else -1) * (5 - abs(x)) for x in reversed(word))
    word = phi(b) + (("slider", 5, 1, False, False, False),) + phi(a)
    witness = normalize_canonical_slider(word)
    check_canonical_slider(witness)
    assert witness.upper.crossing_counts == (0, 2, 3, 0)
    assert witness.lower_inverse.crossing_counts == (0, 1, 1, 1)
    assert witness.twist_powers == (0, 1, 1, 0)
    assert isinstance(witness.upper.gathering, RibbonGathering)
    assert _key(witness.upper.word, 5) == _key(witness.upper.reconstructed_word, 5)
    assert _expanded_coupon_diagram(word, 5) == _expanded_coupon_diagram(witness.normalized_word, 5)
    with pytest.raises(RibbonError):
        recognize_nested_strict(phi(b), labels=(5, 4, 3, 2, 1), initial_ports=(1, 2, 3, 4, 5))


def test_inverse_cancellation_is_recognized_but_writhe_zero_commutator_is_refused():
    inverse_pair = recognize_nested((2, -2), labels=(1, 2, 3))
    check_nested(inverse_pair)
    assert inverse_pair.crossing_counts == (0, 0)
    assert inverse_pair.reconstructed_word == ()
    assert inverse_pair.gathering is not None
    commutator = (1, 1, 2, 2, -1, -1, -2, -2)
    with pytest.raises(RibbonError, match="not equal"):
        recognize_nested(commutator, labels=(1, 2, 3))
    with pytest.raises(RibbonError, match="not equal"):
        check_nested(replace(inverse_pair, word=commutator))
    with pytest.raises(RibbonError, match="pair crossing counts"):
        recognize_nested((2, 2), labels=(1, 2, 3))


def test_saved_exact_action_witness_is_checked_and_immutable():
    witness = recognize_nested((2, -2), labels=(1, 2, 3))
    for changed in [replace(witness, gathering={}),
                    replace(witness, gathering=replace(witness.gathering, action_sha256="wrong")),
                    replace(witness, gathering=replace(witness.gathering, max_operations=True))]:
        with pytest.raises(RibbonError):
            check_nested(changed)
    canonical = recognize_nested((), labels=(1, 2, 3))
    assert canonical.gathering is None


@pytest.mark.parametrize("kwargs", [{"max_image_letters": 0}, {"max_image_letters": True},
    {"max_operations": 0}, {"max_operations": False}])
def test_invalid_gathering_limits_are_refused_even_on_empty_input(kwargs):
    with pytest.raises(RibbonError):
        recognize_nested((), labels=(1, 2, 3), **kwargs)
    with pytest.raises(RibbonError):
        normalize_canonical_slider((("slider", 3, 1, False, False, False),), **kwargs)


@pytest.mark.parametrize("kwargs", [{"max_image_letters": 2}, {"max_image_letters": 3},
                                   {"max_operations": 1}, {"max_operations": 4}])
def test_exact_gathering_limits_are_incomplete_computations(kwargs):
    with pytest.raises(RibbonLimit):
        recognize_nested((2, -2), labels=(1, 2, 3), **kwargs)
    token = ("slider", 3, 1, False, False, False)
    # The upper suffix hierarchy fails strict recognition for the first pair.
    with pytest.raises(RibbonLimit):
        normalize_canonical_slider((1, -1, token), **kwargs)


def test_action_image_growth_has_a_bound_before_intermediate_allocation():
    from zipper_algebra.zipping_algorithm.extraction import _artin_action
    with pytest.raises(RibbonLimit, match="image"):
        _artin_action((1,), 2, max_image_letters=2, max_operations=100, operations=[0])
    expected = ((1, 2, -1), (1,))
    assert _artin_action((1,), 2, max_image_letters=4,
                         max_operations=100, operations=[0]) == expected
