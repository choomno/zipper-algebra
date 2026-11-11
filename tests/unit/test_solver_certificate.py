"""Numeric solver decoding cannot substitute for completion verification."""
from fractions import Fraction

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.f4ncgb import certificate_from_solver
from zipper_algebra.verification.groebner import verify_groebner


def dual_numbers(field=None):
    algebra = FreeAlgebra(field or RationalFunctionField(), ('x',))
    return AlgebraPresentation(algebra, (algebra.generator('x')**2,), ('square',))


def test_numeric_forward_derivations_and_nonmonic_scaling_are_exact():
    presentation = dual_numbers()
    certificate = certificate_from_solver(presentation, '1/2*x*x\nx*x\n', '1/2*[1]\n[i0]\n')
    verified = verify_groebner(presentation, certificate)
    assert verified.dimension == 2
    assert certificate.basis[0] == Fraction(1, 2)*presentation.relations[0]
    assert verified.normal_form(presentation.algebra.generator('x')**3).is_zero()


def test_prime_field_interpretation_and_undefined_rational_rejection():
    presentation = dual_numbers(RationalFunctionField((), 5))
    certificate = certificate_from_solver(presentation, '3*x*x\n', '3*[i0]\n')
    assert verify_groebner(presentation, certificate).dimension == 2
    with pytest.raises(ZeroDivisionError):
        certificate_from_solver(presentation, '1/5*x*x\n', '1/5*[i0]\n')


def test_decoder_does_not_assert_membership_or_completion():
    presentation = dual_numbers()
    false_proof = certificate_from_solver(presentation, 'x*x\n', '2*[i0]\n')
    with pytest.raises(ValueError, match='derivation'):
        verify_groebner(presentation, false_proof)
    with pytest.raises(ValueError):
        certificate_from_solver(presentation, '0\n', '[i0]\n')
    with pytest.raises(ValueError):
        certificate_from_solver(presentation, 'x*x\nx\n', '[i0]\n')


@pytest.mark.parametrize('basis,proof', [
    ('x*x trailing\n', '[i0]\n'), ('x*x\n', '[i0] junk\n'),
    ('x*x\n', '[i0] +\n'), ('x*x\n', '0*[i999]\n'),
    ('x*x\n', '[0]\n'),
])
def test_invalid_or_false_raw_data_are_rejected_at_the_correct_boundary(basis, proof):
    presentation = dual_numbers()
    with pytest.raises(ValueError):
        verify_groebner(presentation, certificate_from_solver(presentation, basis, proof))


@pytest.mark.parametrize('presentation,basis,proof', [
    (None, '', ''), (dual_numbers(), b'x*x', '[i0]'), (dual_numbers(), 'x*x', b'[i0]'),
])
def test_decoder_requires_typed_presentation_and_text(presentation, basis, proof):
    with pytest.raises(TypeError):
        certificate_from_solver(presentation, basis, proof)
