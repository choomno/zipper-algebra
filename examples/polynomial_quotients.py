"""Small exact polynomial-quotient computations, without a Gröbner-basis solver.

Run after installing the package::

    python examples/polynomial_quotients.py
    python examples/polynomial_quotients.py --output /tmp/quotients.json

Each example verifies its literal defining relations in a companion
representation. A cofactor identity supplies the upper dimension bound; the
cyclic coordinate basis supplies the matching lower bound. No numerical
sampling is used. Optional JSON output includes the identities and a provenance
sidecar; the script itself is recorded as an input file.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from zipper_algebra.algebra import AlgebraPresentation, CofactorIdentity
from zipper_algebra.execution.artifacts import prepare_output_root
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.provenance import capture_provenance, write_output_provenance


def _require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def _identity(field, size):
    return [[field.scalar(int(i == j)) for j in range(size)] for i in range(size)]


def _zero(field, size):
    return [[field.zero for _ in range(size)] for _ in range(size)]


def _multiply(field, left, right):
    size = len(left)
    out = _zero(field, size)
    for row in range(size):
        for mid in range(size):
            if left[row][mid].is_zero():
                continue
            for col in range(size):
                if not right[mid][col].is_zero():
                    out[row][col] = out[row][col] + left[row][mid] * right[mid][col]
    return out


def _matrix_polynomial(field, coefficients, matrix):
    """Independent power sum, as opposed to free-algebra Horner evaluation."""
    size = len(matrix)
    out = _zero(field, size)
    power = _identity(field, size)
    for coefficient in coefficients:
        for row in range(size):
            for col in range(size):
                out[row][col] = out[row][col] + coefficient * power[row][col]
        power = _multiply(field, power, matrix)
    return out


def _evaluate(element, matrices):
    field = element.algebra.field
    size = len(matrices[0])
    out = _zero(field, size)
    for word, coefficient in element.terms:
        term = _identity(field, size)
        for letter in word:
            term = _multiply(field, term, matrices[letter - 1])
        for row in range(size):
            for col in range(size):
                out[row][col] = out[row][col] + coefficient * term[row][col]
    return out


def _polynomial(algebra, coefficients):
    x = algebra.generator('x1')
    out = algebra.zero
    power = algebra.one
    for coefficient in coefficients:
        out = out + algebra.scalar(coefficient) * power
        power = power * x
    return out


def certify_two_strand(presentation, q, p, expanded_product):
    """Certify a two-strand quotient using a supplied, expanded product.

    ``expanded_product`` is independently specified as p(x)(x-1)(x+q).
    Its low zero coefficients are removed using the invertibility of x. The
    returned cofactor identity verifies that removal in the *free* algebra,
    rather than assuming that Laurent cancellation holds there.
    """
    algebra = presentation.algebra
    field = algebra.field
    q = field.scalar(q)
    if algebra.generators != ('x1', 'y1'):
        raise ValueError('the companion example requires the two-strand x1,y1 alphabet')
    if not {'INV_1_p', 'INV_1_m', 'T_2'} <= set(presentation.labels):
        raise ValueError('the two inverse relations and level-two untwisting must be present')
    x, y = algebra.generator('x1'), algebra.generator('y1')
    for label, expected in (('INV_1_p', x*y-algebra.one), ('INV_1_m', y*x-algebra.one)):
        if presentation.relations[presentation.labels.index(label)] != expected:
            raise ValueError(f'{label} must be the literal inverse relation')
    coefficients = [field.scalar(value) for value in expanded_product]
    if not any(not coefficient.is_zero() for coefficient in coefficients):
        raise ValueError('the expanded quotient polynomial must be nonzero')
    while coefficients[-1].is_zero():
        coefficients.pop()
    removed = 0
    while coefficients[0].is_zero():
        coefficients.pop(0)
        removed += 1
    size = len(coefficients) - 1
    if size < 1:
        raise ValueError('these examples require a nonzero quotient of positive dimension')
    leading = coefficients[-1]
    monic = [coefficient / leading for coefficient in coefficients]
    matrix_x = _zero(field, size)
    for column in range(size - 1):
        matrix_x[column + 1][column] = field.one
    for row in range(size):
        matrix_x[row][-1] = -monic[row]
    # Invert the companion matrix using its nonzero constant coefficient.
    matrix_y = _matrix_polynomial(field, [-c / monic[0] for c in monic[1:]], matrix_x)
    unit = _identity(field, size)
    _require(_multiply(field, matrix_x, matrix_y) == unit, 'x*y is not the identity matrix')
    _require(_multiply(field, matrix_y, matrix_x) == unit, 'y*x is not the identity matrix')
    _require(all(_evaluate(relation, [matrix_x, matrix_y]) == _zero(field, size)
                 for relation in presentation.relations),
             'a defining relation does not vanish in the companion representation')
    # x^j e_0=e_j, independently proving that the proposed words are independent.
    power = unit
    for column in range(size):
        _require([row[0] for row in power] == [field.scalar(int(row == column))
                                              for row in range(size)],
                 'the claimed cyclic words do not give the coordinate basis')
        power = _multiply(field, power, matrix_x)

    polynomial = _polynomial(algebra, p)
    reduced_product = _polynomial(algebra, coefficients)
    labels = presentation.labels
    inverse_forward = labels.index('INV_1_p')
    inverse_reverse = labels.index('INV_1_m')
    untwisting = labels.index('T_2')
    # F=q*p*(xy-1)-x*T_2. For F=x^k*G, remove x^k with:
    # y^k*x^k-1 = sum_{j=0}^{k-1} y^j*(yx-1)*x^j.
    terms = [(y**removed * algebra.scalar(q) * polynomial, inverse_forward, algebra.one),
             (-(y**removed * x), untwisting, algebra.one)]
    terms.extend((-(y**j), inverse_reverse, x**j * reduced_product)
                 for j in range(removed))
    certificate = CofactorIdentity(presentation.content_hash(), reduced_product, tuple(terms))
    _require(certificate.verify(presentation), 'the proposed upper-bound cofactor identity is false')
    matrix_witness = {
        'schema': 'zipper-two-strand-matrix-example-v1',
        'dimension': size,
        # Sparse matrix coordinates are zero-based; free-word indices are one-based.
        'generator_actions': {
            name: [[row, col, matrix[row][col].to_json()]
                   for row in range(size) for col in range(size)
                   if not matrix[row][col].is_zero()]
            for name, matrix in zip(algebra.generators, (matrix_x, matrix_y))
        },
        'basis_words': [[1]*power for power in range(size)],
        'cyclic_vector': [[0, field.one.to_json()]],
    }
    return {
        'dimension': size,
        'removed_invertible_x_factors': removed,
        'literal_relations_verified': len(presentation.relations),
        'cyclic_basis_verified': True,
        'cofactor_verified': True,
        'presentation_hash': presentation.content_hash(),
        'presentation': presentation.to_json(),
        'field': field.to_json(),
        'cofactor': certificate.to_json(),
        'matrix_witness': matrix_witness,
    }, (matrix_x, matrix_y)


def run_examples():
    """Return small exact controls, including distinct exceptional specializations."""
    examples = {}
    field = RationalFunctionField(('q', 'a', 'b'))
    q, a, b = (field.parameter(name) for name in ('q', 'a', 'b'))
    p = [b, a, field.one]
    # Explicit quartic coefficients are the independent oracle for the constructor.
    quartic = [-b*q, b*(q-1)-a*q, b+a*(q-1)-q, a+q-1, field.one]
    presentation = AlgebraPresentation.zipper(2, field=field, q=q, untwistings={2: p})
    examples['three_parameters'], _ = certify_two_strand(presentation, q, p, quartic)
    assignments = {'q': 2, 'a': 1, 'b': 3}
    prime_presentation = presentation.specialize(assignments, characteristic=5)
    examples['prime_field_specialization'], _ = certify_two_strand(
        prime_presentation, q.specialize(assignments, characteristic=5),
        [c.specialize(assignments, characteristic=5) for c in p],
        [c.specialize(assignments, characteristic=5) for c in quartic])

    repeated_field = RationalFunctionField(('u',))
    u = repeated_field.parameter('u')
    q2 = repeated_field.scalar(2)
    repeated_p = [u*u, -2*u, repeated_field.one]
    repeated_quartic = [-2*u*u, u*u+4*u, u*u-2*u-2, 1-2*u, repeated_field.one]
    repeated_presentation = AlgebraPresentation.zipper(
        2, field=repeated_field, q=q2, untwistings={2: repeated_p})
    repeated, matrices = certify_two_strand(
        repeated_presentation, q2, repeated_p, repeated_quartic)
    algebra = repeated_presentation.algebra
    x = algebra.generator('x1')
    h = (x-algebra.one)*(x+algebra.scalar(2))
    nilpotent = (x-algebra.scalar(u))*h
    nilpotent_matrix = _evaluate(nilpotent, matrices)
    _require(nilpotent_matrix != _zero(repeated_field, 4), 'the square-zero element is zero')
    _require(_multiply(repeated_field, nilpotent_matrix, nilpotent_matrix) == _zero(repeated_field, 4),
             'the proposed nilpotent matrix does not have square zero')
    polynomial = _polynomial(algebra, repeated_p)
    nilpotent_certificate = CofactorIdentity(
        repeated_presentation.content_hash(), nilpotent*nilpotent,
        ((h*algebra.scalar(2)*polynomial, repeated_presentation.labels.index('INV_1_p'), algebra.one),
         (-(h*x), repeated_presentation.labels.index('T_2'), algebra.one)))
    _require(nilpotent_certificate.verify(repeated_presentation),
             'the proposed square-zero cofactor identity is false')
    repeated['nonzero_square_zero_element_verified'] = True
    repeated['square_zero_element'] = nilpotent.to_json()
    repeated['square_zero_cofactor'] = nilpotent_certificate.to_json()
    examples['repeated_root'] = repeated

    parameter_field = RationalFunctionField(('s',))
    s = parameter_field.parameter('s')
    one = parameter_field.one
    q2 = parameter_field.scalar(2)
    cases = {
        # p=s*x²+x+1: top degree drops at s=0.
        'leading_coefficient': ([one, one, s], [-2*one, -one, 2*one-2*s, one+s, s]),
        # p=x²-s: at s=0 its x² factor is a Laurent unit.
        'constant_coefficient': ([-s, parameter_field.zero, one], [2*s, -s, -s-2, one, one]),
    }
    for name, (polynomial, product) in cases.items():
        family = AlgebraPresentation.zipper(2, field=parameter_field, q=q2,
                                           untwistings={2: polynomial})
        examples[name+'_generic'], _ = certify_two_strand(family, q2, polynomial, product)
        specialized = family.specialize({'s': 0})
        examples[name+'_zero'], _ = certify_two_strand(
            specialized, q2.specialize({'s': 0}),
            [c.specialize({'s': 0}) for c in polynomial],
            [c.specialize({'s': 0}) for c in product])
    return examples


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write exact identities and provenance to JSON')
    args = parser.parse_args(argv)
    producer = capture_provenance() if args.output is not None else None
    report = run_examples()
    for name, result in report.items():
        print(f"{name}: dimension {result['dimension']}; exact cofactor and cyclic representation verified")
    if args.output is not None:
        if not args.output.parent.is_dir():
            prepare_output_root(args.output.parent)
        args.output.write_text(json.dumps(report, sort_keys=True, indent=2)+'\n', encoding='utf-8')
        write_output_provenance(args.output, producer=producer, inputs=[__file__])
    return report


if __name__ == '__main__':
    main()
