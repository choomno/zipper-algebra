"""Independent original-Z cofactor checks and malformed witness regressions."""
from copy import deepcopy
import json

import pytest

from zipper_algebra.zipping_algorithm.catalog import Catalog, CatalogError, condition
from zipper_algebra.zipping_algorithm.verification import (
    LocalVerificationError, LocalVerificationLimit, VerificationLimits,
    certify_rule, check_local_certificate,
)


def _small_binding(rule, catalog):
    for n in range(2, 5):
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                for k in range(n + 1):
                    for l in range(1, n):
                        env = dict(n=n, i=i, j=j, k=k, l=l, eps=1)
                        if not condition(rule["condition"], env):
                            continue
                        try:
                            catalog.instantiate(rule["id"], env, expand_sliders=True)
                        except CatalogError:
                            continue
                        return env
    raise AssertionError("missing bounded witness for " + rule["id"])


def _nearest(mode="I", **kwargs):
    return certify_rule("zip.nearest.negative.i_ge_3." + mode,
                        dict(n=5, i=3, j=1, k=2, l=1), **kwargs)


def test_every_rule_family_and_orientation_against_original_definitions():
    catalog = Catalog.bundled()
    checks = 0
    for base in catalog.rules:
        if not base["id"].endswith(".I"):
            continue
        parameters = _small_binding(base, catalog)
        assert parameters["n"] <= 4
        for mode in ("I", "M", "phi", "Mphi"):
            certificate = certify_rule(base["base_id"] + "." + mode, parameters, catalog=catalog)
            assert check_local_certificate(json.loads(json.dumps(certificate)), catalog=catalog)
            assert certificate["coefficient_ring"] == "Z[q,q^-1]"
            if "original_relator_cofactor" in base:
                assert certificate["cofactors"]
                assert certificate["cofactors"][0]["relator"]["symmetry"] == mode
            else:
                assert certificate["cofactors"] == []
            checks += 1
    assert checks == 128


def test_anti_symmetry_transports_both_contexts_in_reversed_order():
    base = _nearest()
    phi = _nearest("phi")
    cofactor, = base["cofactors"]
    transformed, = phi["cofactors"]
    assert cofactor["left"] == [-1, 3, 4]
    assert transformed["left"] == []
    assert transformed["right"] == [1, 2, -4]
    assert transformed["q_power"] == cofactor["q_power"]
    assert transformed["integer"] == cofactor["integer"]
    assert check_local_certificate(phi)


@pytest.mark.parametrize("mutation", [
    lambda c: c.update(extra=True),
    lambda c: c.update(schema="other"),
    lambda c: c.update(scope="free_algebra_certificate"),
    lambda c: c.update(coefficient_ring="Q"),
    lambda c: c.update(source_sha256="0" * 64),
    lambda c: c.update(catalog_sha256="0" * 64),
    lambda c: c.update(rule_sha256="0" * 64),
    lambda c: c["lhs"][0].update(integer=999),
    lambda c: c["cofactors"][0].update(integer=1),
    lambda c: c["cofactors"][0].update(left=[0]),
    lambda c: c["cofactors"][0]["relator"].update(family="T"),
    lambda c: c["cofactors"][0]["relator"].update(level=True),
    lambda c: c["cofactors"][0]["relator"].update(start=0),
    lambda c: c["cofactors"][0]["relator"].update(symmetry="unknown"),
    lambda c: c["cofactors"][0]["relator"].update(definition_sha256="0" * 64),
    lambda c: c.update(cofactors=None),
    lambda c: c.update(parameters={"n": True}),
    lambda c: c.update(parameters={"n": 0}),
    lambda c: c.update(parameters={}),
    lambda c: c.update(words={"P": [False]}),
    lambda c: c.update(words={"P": None}),
])
def test_saved_cofactor_tampering_refused(mutation):
    certificate = _nearest()
    mutation(certificate)
    with pytest.raises((LocalVerificationError, CatalogError)):
        check_local_certificate(certificate)


def test_missing_cofactor_does_not_promote_nonzero_residual_to_identity():
    certificate = _nearest()
    certificate["cofactors"] = []
    with pytest.raises(LocalVerificationError, match="nonzero"):
        check_local_certificate(certificate)
    catalog = Catalog.bundled()
    rule = catalog.get("zip.nearest.negative.i_ge_3.I")
    del rule["original_relator_cofactor"]
    altered = Catalog({"schema_version": "1.0", "rules": [rule]})
    with pytest.raises(LocalVerificationError, match="nonzero"):
        _nearest(catalog=altered)


@pytest.mark.parametrize("limits", [
    VerificationLimits(max_operations=1),
    VerificationLimits(max_image_letters=4),
    VerificationLimits(max_integer_bits=1),
    VerificationLimits(max_word_size=1),
    VerificationLimits(max_terms=1),
])
def test_explicit_verification_limits_refuse_incomplete_checks(limits):
    with pytest.raises((LocalVerificationLimit, CatalogError)):
        _nearest(limits=limits)


def test_malformed_limits_and_parameters():
    with pytest.raises(LocalVerificationError):
        VerificationLimits(max_operations=True)
    for parameters, words in [(None, {}), ({"n": 3}, []), ({3: 2, "n": 3}, {}),
                              ({"n": 3}, {3: []})]:
        with pytest.raises(LocalVerificationError):
            certify_rule("single.far.positive.I", parameters, words=words)


def test_valid_different_cofactor_representation_is_checked_as_an_identity():
    certificate = _nearest()
    first = certificate["cofactors"][0]
    # Split the same coefficient into two integer summands; no template equality
    # is required, since the checker verifies the actual contextual identity.
    doubled = deepcopy(first)
    doubled["integer"] *= 2
    opposite = deepcopy(first)
    opposite["integer"] *= -1
    certificate["cofactors"] = [doubled, opposite]
    assert check_local_certificate(certificate)


def test_malformed_expanded_sides_and_large_saved_lists():
    for value in (None, [{"integer": 1, "q_power": 0, "word": []}] * 2):
        certificate = _nearest()
        certificate["lhs"] = value
        with pytest.raises(LocalVerificationError):
            check_local_certificate(certificate)
    certificate = _nearest()
    certificate["lhs"] = certificate["lhs"] * 101
    with pytest.raises(LocalVerificationLimit, match="term"):
        check_local_certificate(certificate, limits=VerificationLimits(max_terms=100))
    certificate = _nearest()
    certificate["cofactors"] *= 101
    with pytest.raises(LocalVerificationLimit, match="cofactor"):
        check_local_certificate(certificate, limits=VerificationLimits(max_terms=100))


def test_guard_and_word_encoding_limits_are_checked_before_algebra():
    with pytest.raises(LocalVerificationError, match="catalog"):
        certify_rule("single.far.positive.I", dict(n=3, i=1, j=1, k=0, l=1))
    with pytest.raises(LocalVerificationLimit, match="word"):
        certify_rule("single.far.positive.I", dict(n=3, i=1, j=1, k=0, l=2),
                     words={"P": [1, 1]}, limits=VerificationLimits(max_word_size=1))


@pytest.mark.parametrize("cap", [2, 3])
def test_artin_intermediate_and_complete_images_are_bounded(cap):
    with pytest.raises(LocalVerificationLimit, match="Artin image"):
        certify_rule("zip.internal.i2.positive.I", dict(n=2, i=2, j=1, k=0, l=1),
                     limits=VerificationLimits(max_image_letters=cap))


def test_relator_expansion_word_bound_applies_even_to_zero_cofactors():
    certificate = certify_rule("single.far.positive.I", dict(n=4, i=1, j=1, k=0, l=3))
    relator_source = certify_rule("zip.nearest.negative.i_ge_3.I", dict(n=4, i=3, j=1, k=1, l=1))
    term = deepcopy(relator_source["cofactors"][0])
    term.update(integer=0, left=[], right=[])
    certificate["cofactors"] = [term]
    with pytest.raises(LocalVerificationLimit, match="word"):
        check_local_certificate(certificate, limits=VerificationLimits(max_word_size=3))


def test_invalid_external_cofactor_relator_family_refused():
    rule = Catalog.bundled().get("zip.nearest.negative.i_ge_3.I")
    rule["original_relator_cofactor"]["relator"]["family"] = "T"
    altered = Catalog({"schema_version": "1.0", "rules": [rule]})
    with pytest.raises(LocalVerificationError, match="family"):
        _nearest(catalog=altered)


def test_intermediate_uncollected_difference_term_limit():
    # A syntactically valid external candidate exceeds the difference budget
    # before its (false) braid-group identity is considered.
    rule = Catalog.bundled().get("single.far.positive.I")
    rule["lhs"] = [{"coefficient": [{"integer": 1, "q_power": 0}],
                    "word": {"type": "sigma", "index": 1}}]
    rule["rhs"] = [{"coefficient": [{"integer": 1, "q_power": 0}],
                    "word": {"type": "sigma", "index": 2}},
                   {"coefficient": [{"integer": 1, "q_power": 0}],
                    "word": {"type": "identity"}}]
    altered = Catalog({"schema_version": "1.0", "rules": [rule]})
    with pytest.raises(LocalVerificationLimit, match="term"):
        certify_rule(rule["id"], dict(n=3, i=1, j=1, k=0, l=2), catalog=altered,
                     limits=VerificationLimits(max_terms=2))


def test_every_applicable_zipper_family_has_independently_verified_BM_instances():
    catalog = Catalog.bundled()
    checks = 0
    covered = set()
    for base in catalog.rules:
        if not base["id"].startswith("zip.") or not base["id"].endswith(".I"):
            continue
        for size in (2, 3, 4):
            binding = None
            for n in range(size, size + 4):
                for j in range(1, n - size + 2):
                    for k in range(n - size + 1):
                        for l in range(1, n):
                            env = dict(n=n, i=size, j=j, k=k, l=l)
                            if condition(base["condition"], env):
                                binding = env
                                break
                        if binding:
                            break
                    if binding:
                        break
                if binding:
                    break
            if binding is None:
                continue  # Rank-specific exceptional schemata have disjoint guards.
            certificate = certify_rule(base["base_id"] + ".BM", binding, catalog=catalog)
            assert check_local_certificate(json.loads(json.dumps(certificate)), catalog=catalog)
            for term in certificate["cofactors"]:
                assert term["relator"]["symmetry"] == "M"
            checks += 1
            covered.add(base["base_id"])
    assert len(covered) == 20
    assert checks == 44


def test_BM_cofactor_unit_and_contexts_preserve_original_relator():
    for size in (2, 3, 4):
        family = "zip.nearest.negative.i2" if size == 2 else "zip.nearest.negative.i_ge_3"
        parameters = dict(n=size + 2, i=size, j=1, k=2, l=1)
        certificate = certify_rule(family + ".BM", parameters)
        term, = certificate["cofactors"]
        assert term["integer"] == 1
        assert term["q_power"] == 1 - size
        assert term["relator"]["symmetry"] == "M"
        assert term["right"] == []
        assert all(x < 0 for x in term["left"][1:])
        altered = deepcopy(certificate)
        altered["cofactors"][0]["q_power"] += 1
        with pytest.raises(LocalVerificationError, match="nonzero"):
            check_local_certificate(altered)


def test_BM_source_hash_and_derived_rule_identity_are_bound():
    catalog = Catalog.bundled()
    source = catalog.get("zip.nearest.negative.i_ge_3.M")
    source["notes"] += " Updated source record."
    custom = Catalog({"schema_version": "1.0", "rules": [
        catalog.get("zip.nearest.negative.i_ge_3.I"), source]})
    certificate = certify_rule("zip.nearest.negative.i_ge_3.BM", dict(n=4, i=3, j=1, k=1, l=1), catalog=custom)
    certificate["catalog_sha256"] = catalog.sha256
    with pytest.raises(LocalVerificationError, match="rule hash"):
        check_local_certificate(certificate, catalog=catalog)
    with pytest.raises(LocalVerificationError):
        certify_rule("single.far.positive.BM", dict(n=3, i=1, j=1, k=0, l=2))
    with pytest.raises(LocalVerificationError):
        certify_rule(False, dict(n=3, i=1, j=1, k=0, l=2))
