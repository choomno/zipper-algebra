"""Deep derivation graphs retain exact contexts and formal denominators."""
from fractions import Fraction

import pytest

from zipper_algebra import Presentation
from zipper_algebra.coefficients import NumericRing
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.compose import Rat, _order, compose, expand_over_inputs
from zipper_algebra.cofactor.lift import Lift
from zipper_algebra.cofactor.skeleton import Skeleton, evaluate


def term(coefficient=1, *, denominator=(0, 0), left=(), kind="basis", ref="", right=()):
    numerator = {(0, 0): Fraction(coefficient)} if coefficient else {}
    return denominator, numerator, left, kind, ref, right


def atom_family():
    return AtomSet([Atom("q", {(1, 0): 1}), Atom("t", {(0, 1): 1})])


def deep_entries(depth=1501):
    entries = {f"g{index}": [term(ref=f"g{index + 1}")]
               for index in range(depth - 1)}
    entries["g0"] = [term(2, denominator=(1, 0), left=(1,), ref="g1", right=(4,))]
    entries[f"g{depth - 1}"] = [term(3, denominator=(0, 2), left=(2,),
                                         kind="input", ref=0, right=(3,))]
    return entries


def test_forward_chain_beyond_python_recursion_depth_preserves_the_identity():
    presentation = Presentation(3)
    lift = Lift(atom_family(), deep_entries(), 0, 0, (5, 7, 0), {})
    expanded = expand_over_inputs(lift, pres=presentation)
    assert len(expanded) == 1501
    assert expanded["g0"] == {
        ((1, 2), 0, (3, 4)): ({(0, 0): Fraction(6)}, (1, 2)),
    }
    arithmetic = Rat(lift.atoms)
    terms = [(arithmetic.eval(coefficient, 5, 7), left, index, right)
             for (left, index, right), coefficient in expanded["g0"].items()]
    actual = evaluate(terms, presentation.relations(NumericRing(5, 7, 0)), 0)
    # 6/(q*t^2) * x1*x2 * (x1*y1-1) * y1*y2, with q=5,t=7.
    assert {word: coefficient.v for word, coefficient in actual.terms.items()} == {
        (1, 2, 1, 3, 3, 4): Fraction(6, 245),
        (1, 2, 3, 4): Fraction(-6, 245),
    }


def test_deep_zero_coefficient_cycle_is_rejected_before_cancellation():
    entries = deep_entries()
    entries["g1500"] = [term(0, ref="g0")]
    with pytest.raises(ValueError, match="cycle"):
        expand_over_inputs(Lift(atom_family(), entries, 0, 0, (5, 7, 0), {}))


def test_independent_node_order_preserves_certificate_body_and_denominator():
    presentation = Presentation(2)
    atoms = atom_family()
    # Each leaf is the first inverse relation, stored as q/q and t/t.
    entries = {
        "outer": [term(2, ref="b"), term(-1, ref="a")],
        "a": [((0, 1), {(0, 1): Fraction(1)}, (), "input", 0, ())],
        "unused": [term(kind="input", ref=1)],
        "b": [((1, 0), {(1, 0): Fraction(1)}, (), "input", 0, ())],
    }
    skeleton = Skeleton(presentation.slug(), presentation.content_hash(), "input relation",
                        (5, 7, 0), [{"basis_key": "outer", "left": (), "right": (),
                                     "mult": Fraction(1)}])
    expected = ["() * INV_1_p * ()  ::  (1*q^1*t^1) / (q*t)"]
    orderings = (entries, {key: entries[key] for key in ("a", "b", "unused", "outer")})
    first_order, second_order = (_order(ordered) for ordered in orderings)
    assert first_order.index("b") < first_order.index("a")
    assert second_order.index("a") < second_order.index("b")
    for ordered in orderings:
        lift = Lift(atoms, ordered, 0, 0, (5, 7, 0), {})
        certificate = compose(lift, skeleton, presentation,
                              multipliers=[({(0, 0): Fraction(1)}, (0, 0))], verbose=False)
        assert certificate.body_lines() == expected
        assert certificate.evaluate(5, 7, presentation).terms == (
            presentation.relations(NumericRing(5, 7, 0))[0].terms)
