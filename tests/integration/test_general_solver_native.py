"""Tiny installed native acquisition/reconstruction/action integration controls."""
import json
import os
from pathlib import Path

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.execution.config import RunConfig
from zipper_algebra.execution.production import produce_algebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.finite import WordBasisCertificate, verify_word_basis
from zipper_algebra.verification.groebner import GroebnerCertificate, verify_groebner


@pytest.mark.parametrize('io_mode', ['files', 'sqlite'])
def test_native_literal_generic_family_reloads_both_certificates(tmp_path, io_mode):
    binary = os.environ.get('F4NCGB_BINARY')
    if not binary:
        pytest.skip('set F4NCGB_BINARY to the pinned SQLite-enabled solver')
    pytest.importorskip('fflow')
    pytest.importorskip('sympy')
    field = RationalFunctionField(('alpha', 'beta', 'gamma'))
    a, b, c = (field.parameter(name) for name in field.parameters)
    parent = FreeAlgebra(field, ('u', 'v'))
    u, v = parent.generator('u'), parent.generator('v')
    # This is a tiny noncommutative integration control, not a census search.
    pres = AlgebraPresentation(parent, (parent.zero, u*u-a, v*v-b, v*u+u*v-c),
                               ('zero', 'u2', 'v2', 'mixed'))
    config = RunConfig(Path(binary), tmp_path / io_mode, timeout_seconds=15, solver_io=io_mode)
    report = produce_algebra(pres, {'gamma': 5, 'beta': 3, 'alpha': 2}, config)
    assert report['verified'], report
    assert report['dimension'] == 4
    assert report['word_basis_verified'] is True
    assert report['groebner_complete'] is None
    assert 'basis_complete' not in report
    assert report['solver']['request']['input_relation_indices'] == [1, 2, 3]
    artifacts = report['artifacts']
    saved = AlgebraPresentation.from_json(json.loads(Path(artifacts['presentation.json']).read_text()))
    generic = GroebnerCertificate.from_json(json.loads(Path(artifacts['rewrite-candidate.json']).read_text()))
    certificate = WordBasisCertificate.from_json(json.loads(Path(artifacts['word-basis-certificate.json']).read_text()))
    assert saved.content_hash() == pres.content_hash()
    assert verify_groebner(saved, generic).dimension == 4
    algebra = verify_word_basis(saved, certificate)
    assert algebra.dimension == 4
    vu = algebra.coordinates(v*u)
    assert vu == algebra.coordinates(-u*v+c)
    assert vu != algebra.coordinates(u*v)
    assert config.output_root.joinpath('.gitignore').read_text() == '*\n!.gitignore\n'
