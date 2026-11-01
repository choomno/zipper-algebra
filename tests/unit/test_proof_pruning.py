"""Pruning preserves certified actions and every raw dependency of their proofs."""
from dataclasses import replace
from fractions import Fraction

import pytest

from zipper_algebra.action_production import _prune_certificate, word_basis_from_rewrites
from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.cofactor.general import expand_derivation, verify_cofactor
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.derivations import DerivationTerm, _Budget
from zipper_algebra.verification.finite import verify_word_basis
from zipper_algebra.verification.groebner import GroebnerCertificate, ProofTerm, VerificationLimit


def square_zero():
    algebra = FreeAlgebra(RationalFunctionField(), ('x',))
    x = algebra.generator('x')
    presentation = AlgebraPresentation(algebra, (x*x,), ('square',))
    one = algebra.field.one
    proposal = GroebnerCertificate(presentation.content_hash(), (x*x,),
                                   ((ProofTerm(one, (), 'input', 0, ()),),))
    return presentation, proposal, x


def test_unused_false_proposal_is_not_a_claim_in_the_returned_certificate():
    p, g, x = square_zero()
    # x^3 is in the ideal, but this empty derivation does not prove it. The
    # shortest leading word x^2 makes that redundant rule unnecessary.
    g = replace(g, basis=(x*x, x**3), derivations=(g.derivations[0], ()))
    result = word_basis_from_rewrites(p, g)
    assert result.basis_words == ((), (1,))
    assert x**3 not in result.proof.nodes
    assert verify_word_basis(p, result).dimension == 2
    # The original proposal is neither mutated nor asserted to be checked.
    assert g.derivations[1] == () and g.basis[1] == x**3


def test_a_false_identity_retained_by_a_zero_edge_still_fails():
    p, g, x = square_zero()
    zero = x.algebra.field.zero
    g = replace(g, basis=(x*x, x**3), derivations=(
        g.derivations[0] + (ProofTerm(zero, (), 'basis', 1, ()),), ()))
    with pytest.raises(ValueError, match='local input-ideal identity'):
        word_basis_from_rewrites(p, g)


def test_zero_edges_are_retained_and_remapped_even_when_arithmetically_unnecessary():
    p, g, x = square_zero()
    field = x.algebra.field
    # The leading x^4 rule is unused. Its removal shifts the x^2 node and its
    # forward zero reference to the honest contextual x^3 derivation.
    g = replace(g, basis=(x**4, x*x, x**3), derivations=(
        (ProofTerm(field.one, (1,1), 'input', 0, ()),),
        (ProofTerm(field.one, (), 'input', 0, ()),
         ProofTerm(field.zero, (), 'basis', 2, ())),
        (ProofTerm(field.one, (1,), 'input', 0, ()),)))
    result = word_basis_from_rewrites(p, g)
    assert x**4 not in result.proof.nodes
    assert result.proof.nodes[:2] == (x*x, x**3)
    assert result.proof.derivations[0][1] == DerivationTerm(field.zero, (), 'node', 1, ())
    assert verify_word_basis(p, result).dimension == 2


def test_forward_shared_noncommuting_proofs_survive_removed_prefix_and_gaps():
    algebra = FreeAlgebra(RationalFunctionField(), ('x', 'y'))
    x, y = map(algebra.generator, algebra.generators)
    f = algebra.field
    p = AlgebraPresentation(algebra, (x*x, y*y, y*x+x*y), ('xx', 'yy', 'mixed'))
    basis = (x**3, x*x, 2*x**3, y*y, y*x+x*y, 2*x*x, x*y*x*x*y*x)
    rows = ((ProofTerm(f.one, (1,), 'input', 0, ()),),
            (ProofTerm(f.scalar(Fraction(1,2)), (), 'basis', 5, ()),),
            (ProofTerm(f.scalar(2), (), 'basis', 0, ()),),
            (ProofTerm(f.one, (), 'input', 1, ()),),
            (ProofTerm(f.one, (), 'input', 2, ()),),
            (ProofTerm(f.scalar(2), (), 'input', 0, ()),),
            (ProofTerm(f.one, (1,2), 'input', 0, (2,1)),))
    g = GroebnerCertificate(p.content_hash(), basis, rows)
    result = word_basis_from_rewrites(p, g)
    assert result.basis_words == ((), (1,), (2,), (1,2))
    assert result.proof.nodes[:4] == (x*x, y*y, y*x+x*y, 2*x*x)
    assert result.proof.derivations[0][0].index == 3  # old forward reference 5
    assert all(node not in result.proof.nodes for node in (x**3, 2*x**3, basis[-1]))
    finite = verify_word_basis(p, result)
    assert finite.dimension == 4
    for node in {index for row in result.action_nodes for index in row}:
        identity = expand_derivation(p, result.proof, node)
        assert verify_cofactor(p, identity)
    # yxy = -xyy, so a nonempty right context survives the pruning/remapping.
    assert any(term.right for row in result.proof.derivations for term in row)


@pytest.mark.parametrize('coefficient', [0, 1])
@pytest.mark.parametrize('kind', ['input', 'basis'])
def test_original_proposal_bounds_are_checked_before_appending_generated_nodes(coefficient, kind):
    p, g, x = square_zero()
    # In the augmented graph, index 2 would name the generated zero node. It
    # must not turn an invalid original basis reference into a valid one.
    index = 1 if kind == 'input' else 2
    row = (ProofTerm(x.algebra.field.scalar(coefficient), (), kind, index, ()),)
    g = replace(g, basis=(x*x, x**3), derivations=(g.derivations[0], row))
    with pytest.raises(ValueError, match='raw rewrite reference'):
        word_basis_from_rewrites(p, g)


@pytest.mark.parametrize('coefficient', [0, 1])
def test_unused_raw_cycles_are_rejected_before_pruning(coefficient):
    p, g, x = square_zero()
    value = x.algebra.field.scalar(coefficient)
    g = replace(g, basis=(x*x, x**3, 2*x**3), derivations=(g.derivations[0],
        (ProofTerm(value, (), 'basis', 2, ()),),
        (ProofTerm(value, (), 'basis', 1, ()),)))
    with pytest.raises(ValueError, match='cyclic derivation graph'):
        word_basis_from_rewrites(p, g)


def test_unit_root_keeps_only_its_proof_and_remaps_zero_algebra():
    algebra = FreeAlgebra(RationalFunctionField(), ('x', 'y'))
    x, y = map(algebra.generator, algebra.generators)
    f = algebra.field
    p = AlgebraPresentation(algebra, (algebra.scalar(3),), ('unit',))
    # The first proposal has no honest proof, but the scalar rule supersedes it.
    g = GroebnerCertificate(p.content_hash(), (x*y, algebra.scalar(3)),
        ((), (ProofTerm(f.one, (), 'input', 0, ()),)))
    result = word_basis_from_rewrites(p, g)
    assert result.basis_words == () and all(not row for row in result.action_nodes)
    assert result.unit_node == 1
    assert result.proof.nodes == (algebra.scalar(3), algebra.one)
    assert result.proof.derivations[1][0].index == 0
    assert verify_word_basis(p, result).dimension == 0


def test_empty_alphabet_needs_no_proof_nodes():
    algebra = FreeAlgebra(RationalFunctionField(), ())
    p = AlgebraPresentation(algebra, (), ())
    g = GroebnerCertificate(p.content_hash(), (), ())
    result = word_basis_from_rewrites(p, g)
    assert result.proof.nodes == result.proof.derivations == ()
    assert verify_word_basis(p, result).dimension == 1


def test_deep_forward_chain_is_validated_without_python_recursion():
    p, g, x = square_zero()
    count = 1500
    one = x.algebra.field.one
    rows = tuple((ProofTerm(one, (), 'basis', index+1, ()),) for index in range(count-1))
    g = replace(g, basis=(x*x,)*count, derivations=rows+g.derivations)
    result = word_basis_from_rewrites(p, g)
    assert len(result.proof.nodes) == count+2
    assert verify_word_basis(p, result).dimension == 2


def test_pruning_an_already_closed_graph_preserves_the_existing_objects():
    p, g, _ = square_zero()
    result = word_basis_from_rewrites(p, g)
    assert _prune_certificate(p, result, _Budget(1000)) is result


@pytest.mark.parametrize('kind', ['input', 'node'])
def test_pruner_validates_unused_raw_edges_in_a_loaded_finite_certificate(kind):
    p, g, x = square_zero()
    c = word_basis_from_rewrites(p, g)
    f = x.algebra.field
    wrong = DerivationTerm(f.zero, (), kind, 1000, ())
    graph = replace(c.proof, nodes=c.proof.nodes+(x**3,),
                    derivations=c.proof.derivations+((wrong,),))
    c = replace(c, proof=graph)
    with pytest.raises(ValueError, match='raw derivation reference'):
        _prune_certificate(p, c, _Budget(1000))


@pytest.mark.parametrize('zero_algebra', [False, True])
def test_pruner_rejects_missing_action_or_unit_roots(zero_algebra):
    p, g, x = square_zero()
    if zero_algebra:
        p = AlgebraPresentation(x.algebra, (x.algebra.one,), ('unit',))
        g = GroebnerCertificate(p.content_hash(), p.relations,
                               ((ProofTerm(x.algebra.field.one, (), 'input', 0, ()),),))
    c = word_basis_from_rewrites(p, g)
    if zero_algebra:
        c = replace(c, unit_node=len(c.proof.nodes))
    else:
        c = replace(c, action_nodes=((len(c.proof.nodes),)*2,))
    with pytest.raises(ValueError, match='finite proof root'):
        _prune_certificate(p, c, _Budget(1000))


def test_raw_structure_validation_consumes_a_bounded_work_budget():
    p, g, _ = square_zero()
    c = word_basis_from_rewrites(p, g)
    with pytest.raises(VerificationLimit, match='operation budget'):
        _prune_certificate(p, c, _Budget(1))
