"""Exact presentation acquisition with independent literal proof expectations."""
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
import sys
from unittest.mock import patch

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.execution.artifacts import ArtifactError, canonical_json, digest_bytes
from zipper_algebra.execution.config import RunConfig
from zipper_algebra.execution.general import run_algebra_point, solver_certificate
from zipper_algebra.execution import store
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.derivations import VerificationLimit


PRODUCER = '''
import argparse
from pathlib import Path
p=argparse.ArgumentParser()
for flag in ("-i", "-o", "-p"): p.add_argument(flag, required=True)
p.add_argument("--mode", default="normal")
a,_=p.parse_known_args()
lines=Path(a.i).read_text().splitlines()[2:]
basis="\\n".join(line.rstrip(",").replace("-", " -").replace("+", " +") for line in lines)+"\\n"
proof="\\n".join("[i%d]"%i for i in range(len(lines)))+"\\n"
if a.mode=="wrong": proof="[i99]\\n"*len(lines)
if a.mode=="identity": proof="2*[i0]\\n"*len(lines)
if a.mode=="empty": basis=proof=""
if a.mode=="cycle": proof="[0]\\n"*len(lines)
if a.mode=="zero": basis="0\\n";proof="0*[i0]\\n"
if a.mode=="forward": basis=basis+basis;proof="[1]\\n[i0]\\n"
Path(a.o).write_text(basis)
if a.mode=="badbytes": Path(a.p).write_bytes(b"\\xff")
elif a.mode!="missing": Path(a.p).write_text(proof)
'''


@pytest.fixture
def setup(tmp_path):
    field = RationalFunctionField(("c0", "c1"))
    algebra = FreeAlgebra(field, ("a", "b"))
    a, b = algebra.generator("a"), algebra.generator("b")
    c0, c1 = field.parameter("c0"), field.parameter("c1")
    presentation = AlgebraPresentation(algebra, (algebra.zero, a*a + c0*a + c1, b),
                                       ("zero", "quadratic", "kill_b"))
    source = tmp_path / "producer.py"
    source.write_text(PRODUCER)
    config = RunConfig(Path(sys.executable), tmp_path / "out", command_prefix=(source,),
                       flags=(), timeout_seconds=3)
    return presentation, config


def test_exact_rationals_and_zero_input_index_mapping(setup):
    presentation, config = setup
    values = {"c1": Fraction(-3, 2), "c0": Fraction(2, 3)}
    run = run_algebra_point(presentation, values, config)
    assert run.files["input.ms"] == b"a,b\n0\n-3/2+2/3*a+1*a*a,\n1*b\n"
    assert run.record["request"]["point"] == [["c0", "2/3"], ["c1", "-3/2"]]
    assert run.record["request"]["input_relation_indices"] == [1, 2]
    assert run.record["status"] == "proof-validated"
    assert not run.record["basis_complete"]
    certificate = solver_certificate(run, presentation)
    assert certificate.presentation_hash == presentation.specialize(values).content_hash()
    assert tuple(row[0].index for row in certificate.derivations) == (1, 2)
    assert certificate.basis == presentation.specialize(values).relations[1:]
    assert run.files["proof.txt"] == b"[i0]\n[i1]\n"


def test_cache_refresh_copy_and_proof_rechecking(setup, tmp_path):
    presentation, config = setup
    config = replace(config, cache_root=tmp_path / "cache")
    values = {"c0": 3, "c1": 2}
    first = run_algebra_point(presentation, values, config)
    assert not first.hit
    reused = run_algebra_point(presentation, dict(reversed(list(values.items()))), config)
    assert reused.hit and first.key == reused.key
    copied = run_algebra_point(presentation, values, replace(config, output_root=tmp_path / "other"))
    assert copied.hit and copied.files == first.files
    assert not run_algebra_point(presentation, values, config, refresh=True).hit
    # Coherently rehashed external bytes must still pass the mathematical check.
    files = {**first.files, "proof.txt": b"2*[i0]\n[i1]\n"}
    store.publish(first.store, first.key, files, first.record["request"], 2)
    with pytest.raises(ArtifactError, match="local input-ideal identity"):
        run_algebra_point(presentation, values, config)
    # Cached count is independently checked against decoded rows.
    store.publish(first.store, first.key, first.files, first.record["request"], 99)
    with pytest.raises(ArtifactError, match="entry count"):
        run_algebra_point(presentation, values, config)


@pytest.mark.parametrize("mode", ["wrong", "identity", "empty", "cycle", "zero", "badbytes", "missing"])
def test_false_incomplete_and_malformed_outputs_do_not_publish(setup, mode):
    presentation, config = setup
    with pytest.raises(ArtifactError):
        run_algebra_point(presentation, {"c0": 3, "c1": 2},
                          replace(config, flags=("--mode", mode)))
    assert not (config.output_root / "samples.sqlite3").exists()


def test_forward_shared_references_are_retained(setup):
    presentation, config = setup
    presentation = AlgebraPresentation(presentation.algebra, presentation.relations[1:2], ("quad",))
    run = run_algebra_point(presentation, {"c0": 3, "c1": 2},
                            replace(config, flags=("--mode", "forward")))
    certificate = solver_certificate(run, presentation)
    assert certificate.derivations[0][0].kind == "basis"
    assert certificate.derivations[0][0].index == 1


def test_full_ascii_alphabet_and_constant_field_are_preserved(setup):
    _, config = setup
    algebra = FreeAlgebra(RationalFunctionField(()), ("a", "a10", "_z"))
    elements = tuple(algebra.generator(name) for name in algebra.generators)
    presentation = AlgebraPresentation(algebra, elements, ("first", "prefix", "underscore"))
    run = run_algebra_point(presentation, {}, config)
    assert run.files["input.ms"] == b"a,a10,_z\n0\n1*a,\n1*a10,\n1*_z\n"
    assert solver_certificate(run, presentation).basis == elements


def test_replay_accepts_older_literal_single_generator_input(setup):
    _, config = setup
    algebra = FreeAlgebra(RationalFunctionField(()), ("x1",))
    relation = algebra.generator("x1") - 1
    presentation = AlgebraPresentation(algebra, (relation,), ("linear",))
    with patch("zipper_algebra.execution.general._solver_names", return_value=["x1"]):
        old = run_algebra_point(presentation, {}, config)
    new = run_algebra_point(presentation, {}, config)
    assert old.files["input.ms"] == b"x1\n0\n-1+1*x1\n"
    assert new.files["input.ms"] == b"z01\n0\n-1+1*z01\n"
    assert old.key != new.key
    assert solver_certificate(old, presentation) == solver_certificate(new, presentation)


@pytest.mark.parametrize("values", [{}, {"c0": 3}, {"c0": 3, "c1": 4, "x": 0},
                                     {"c0": True, "c1": 4}, {"c0": 3.0, "c1": 4}, []])
def test_requires_complete_exact_rational_point(setup, values):
    presentation, config = setup
    with pytest.raises(ValueError, match="exact rational"):
        run_algebra_point(presentation, values, config)
    assert not config.output_root.exists()


def test_wrong_parent_characteristic_poles_and_trivial_inputs(setup):
    presentation, config = setup
    with pytest.raises(TypeError, match="AlgebraPresentation"):
        run_algebra_point(None, {}, config)
    prime = FreeAlgebra(RationalFunctionField((), 7), ("a",))
    with pytest.raises(ValueError, match="characteristic zero"):
        run_algebra_point(AlgebraPresentation(prime, (prime.generator("a"),), ("a",)), {}, config)
    algebra = presentation.algebra
    # Build a literal coefficient with a genuine pole at c0=0.
    pole = AlgebraPresentation(algebra, (algebra.element({(1,): 1/algebra.field.parameter("c0")}),), ("pole",))
    with pytest.raises(ZeroDivisionError):
        run_algebra_point(pole, {"c0": 0, "c1": 1}, config)
    for trivial in (AlgebraPresentation(algebra, (), ()),
                    AlgebraPresentation(algebra, (algebra.zero,), ("zero",))):
        with pytest.raises(ValueError, match="no numerical proof acquisition"):
            run_algebra_point(trivial, {"c0": 1, "c1": 1}, config)
    scalars = FreeAlgebra(RationalFunctionField(()), ())
    with pytest.raises(ValueError, match="one generator"):
        run_algebra_point(AlgebraPresentation(scalars, (scalars.one,), ("one",)), {}, config)


@pytest.mark.parametrize("limit", [0, -1, True, 1.0])
def test_positive_validation_limits(setup, limit):
    presentation, config = setup
    with pytest.raises(ValueError, match="positive integers"):
        run_algebra_point(presentation, {"c0": 3, "c1": 2}, config, max_nodes=limit)


def test_validation_limits_are_enforced_before_publication(setup):
    presentation, config = setup
    for limits in ({"max_operations": 1}, {"max_nodes": 1}):
        with pytest.raises(VerificationLimit):
            run_algebra_point(presentation, {"c0": 3, "c1": 2}, config, **limits)
    assert not (config.output_root / "samples.sqlite3").exists()


def test_v4_schema_and_certificate_binding(setup):
    presentation, config = setup
    run = run_algebra_point(presentation, {"c0": 3, "c1": 2}, config)
    request = run.record["request"]
    for malformed in ({**request, "target": "private"}, {key: value for key, value in request.items()
                                                         if key != "point"}):
        with pytest.raises(ArtifactError, match="request-v4"):
            store.publish(run.store, digest_bytes(canonical_json(malformed)), run.files, malformed, 2)
    changed = AlgebraPresentation(presentation.algebra, presentation.relations, ("changed", "quad", "b"))
    with pytest.raises(ArtifactError, match="literal algebra presentation"):
        solver_certificate(run, changed)
    for field, value in (("point", [["c0", "03"], ["c1", "2"]]),
                         ("input_relation_indices", [0, 1]),
                         ("specialized_presentation_hash", "0"*64),
                         ("input_sha256", "0"*64),
                         ("point", [["missing", "1"]]),
                         ("point", [["c0", "1/0"], ["c1", "2"]])):
        bad = replace(run, record={**run.record, "request": {**request, field: value}})
        with pytest.raises(ArtifactError):
            solver_certificate(bad, presentation)
    bad = replace(run, files={**run.files, "input.ms": b"changed"})
    with pytest.raises(ArtifactError, match="specialized input"):
        solver_certificate(bad, presentation)


def test_native_sqlite_transport_uses_same_validation(setup):
    presentation, config = setup
    reference = run_algebra_point(presentation, {"c0": 3, "c1": 2}, config)
    sqlite_config = replace(config, output_root=config.output_root / "sqlite", solver_io="sqlite")
    with patch("zipper_algebra.execution.native_io.execute", return_value=(reference.files, 0.25)):
        native = run_algebra_point(presentation, {"c0": 3, "c1": 2}, sqlite_config)
    assert native.files == reference.files
    assert native.seconds == 0.25 and native.key != reference.key
