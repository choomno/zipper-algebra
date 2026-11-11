"""Proof validation uses local DAG equations, including cancelled contexts."""
from fractions import Fraction

import pytest

from zipper_algebra.coefficients import Coef
from zipper_algebra.execution.artifacts import ArtifactError, validate_proofs
from zipper_algebra.polynomial import Poly


def artifacts(basis, proofs):
    return {'basis.txt': basis.encode(), 'proof.txt': proofs.encode()}


@pytest.mark.timeout(2)
def test_shared_derivations_do_not_require_exponential_cofactor_expansion():
    # Every named polynomial is 1.  The local equation is a*1*b-ab*1+1=1.
    # Expanding the contextual proof instead produces exponentially many
    # distinct (left,input,right) triples that cancel only as free polynomials.
    count = 40
    proofs = ['[i0]'] + [f'x1*[{index-1}]*x2 -x1*x2*[{index-1}] +[i0]'
                         for index in range(1, count)]
    inputs = [Poly.scalar(Coef.one(0))]
    assert validate_proofs(artifacts('1\n' * count, '\n'.join(proofs) + '\n'),
                           inputs, ['x1', 'x2']) == count


def test_forward_references_are_checked_locally_but_cycles_are_rejected():
    one = Poly.scalar(Coef.one(0))
    assert validate_proofs(artifacts('1\n1\n', '[1]\n[i0]\n'), [one], ['x']) == 2
    with pytest.raises(ArtifactError, match='cyclic proof'):
        validate_proofs(artifacts('1\n1\n', '[1]\n[0]\n'), [one], ['x'])
    # A formally cancelled dependency must still be checked for circularity.
    with pytest.raises(ArtifactError, match='cyclic proof'):
        validate_proofs(artifacts('1\n', '[i0] +0*[0]\n'), [one], ['x'])


def test_each_local_equation_is_required_even_when_another_one_is_correct():
    one = Poly.scalar(Coef.one(0))
    with pytest.raises(ArtifactError, match='basis entry 1 fails its input-ideal identity'):
        validate_proofs(artifacts('1\n1\n', '[i0]\n2*[0]\n'), [one], ['x'])
    with pytest.raises(ArtifactError, match='raw proof reference'):
        validate_proofs(artifacts('1\n', '[i0] +0*[i1]\n'), [one], ['x'])


def test_noncommuting_contexts_and_rational_scalars_are_preserved():
    relation = Poly(0)
    relation.add_term((1, 2), Coef(1, 0))
    relation.add_term((2, 1), Coef(-1, 0))
    source = artifacts('1/2*x2*x1*x2*x1 -1/2*x2*x2*x1*x1\n',
                       '1/2*x2*[i0]*x1\n')
    assert validate_proofs(source, [relation], ['x1', 'x2']) == 1
    wrong = dict(source, **{'proof.txt': b'1/2*x1*[i0]*x2\n'})
    with pytest.raises(ArtifactError, match='input-ideal identity'):
        validate_proofs(wrong, [relation], ['x1', 'x2'])
    assert relation.terms[(1, 2)].v == Fraction(1)
