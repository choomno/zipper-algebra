"""Independent caches are derived from Python execution, never candidate output."""
from copy import deepcopy

import pytest

from zipper_algebra.zipping_algorithm import collected, independent_replay, native
from zipper_algebra.zipping_algorithm.interpreter import Limits


def fake_session():
    session = object.__new__(native.NativeSession)
    session.metadata = {"source_sha256": {"engine.cpp": "native"}, "abi": native.ABI,
                        "binary_sha256": "binary", "compiler": "g++", "compiler_version": "test",
                        "flags": [], "wrapper_source_sha256": "wrapper", "link_libraries": []}
    session._catalog_sha256 = collected.Catalog.bundled().sha256
    return session


def certificate(word=(1, -2), n=3, trace=False):
    session = fake_session()
    result = native._json_arrays(collected.CollectedSession().reduce(word, n=n, trace=trace))
    result.update(schema="zipper-native-run-v1", scope="literal_word_native_collected_spanning_procedure",
                  python_reference_source_sha256=result["source_sha256"],
                  source_sha256=session.metadata["source_sha256"],
                  native_build={key: value for key, value in session.metadata.items() if key != "source_sha256"})
    return session, result


@pytest.mark.parametrize("trace", [False, True])
def test_independent_full_replay_retains_no_native_outputs(monkeypatch, trace):
    session, result = certificate(trace=trace)
    verifier = independent_replay.IndependentVerifier()
    assert verifier.check_native(session, result)
    before = verifier.cache_info()
    assert before["evaluations"] > 0
    monkeypatch.setattr(collected, "_transition", lambda *args: pytest.fail("warm replay evaluated"))
    assert verifier.check_native(session, result)
    after = verifier.cache_info()
    assert after["evaluations"] == before["evaluations"]
    assert after["hits"] > before["hits"] and after["cross_input_hits"] > 0
    changed = deepcopy(result)
    changed["output"][0]["coefficient"][0][1] += 1
    with pytest.raises(native.NativeError, match="differs from independent"):
        verifier.check_native(session, changed)
    assert verifier.check_native(session, result)


@pytest.mark.parametrize("field,value", [("n", True), ("trace_enabled", 0), ("input_word", [1.0]),
                                         ("status", "invented"), ("scope", "invented"), ("unexpected", True)])
def test_strict_candidate_fields(field, value):
    session, result = certificate(word=[1], n=2)
    verifier = independent_replay.IndependentVerifier()
    assert verifier.check_native(session, result)
    result[field] = value
    with pytest.raises(native.NativeError):
        verifier.check_native(session, result)


def test_source_drift_rejects_before_warm_cache(monkeypatch):
    session, result = certificate()
    verifier = independent_replay.IndependentVerifier()
    assert verifier.check_native(session, result)
    monkeypatch.setattr(independent_replay, "source_identity", lambda: {"changed": "source"})
    with pytest.raises(native.NativeError, match="sources changed"):
        verifier.check_native(session, result)


@pytest.mark.parametrize("kwargs", [{"cache_bytes": -1}, {"cache_entries": True},
                                   {"rule_cache_bytes": 1.0}, {"rule_cache_entries": -2}])
def test_invalid_budgets(kwargs):
    with pytest.raises(collected.CollectedError):
        independent_replay.IndependentVerifier(**kwargs)


@pytest.mark.parametrize("kwargs", [{"cache_bytes": 0}, {"cache_entries": 0},
                                   {"cache_bytes": 1}, {"cache_entries": 1}])
def test_disabled_oversized_and_eviction(kwargs):
    session, result = certificate()
    verifier = independent_replay.IndependentVerifier(**kwargs)
    assert verifier.check_native(session, result)
    info = verifier.cache_info()
    if kwargs.get("cache_entries") == 1:
        assert info["evictions"] > 0 and info["entries"] == 1
    else:
        assert info["entries"] == 0
    if kwargs.get("cache_bytes") == 1:
        assert info["oversized"] > 0
    verifier.clear_cache()
    assert verifier.cache_info()["evaluations"] == 0


def test_key_retains_context_labels_and_limits(monkeypatch):
    verifier = independent_replay.IndependentVerifier()
    monkeypatch.setattr(collected, "_transition", lambda work, limits, trace, catalog:
                        {"children": [deepcopy(work)], "output": [], "unit_evidence": {}})
    base = {"phase": "extract", "depth": 0, "coefficient": [[0, 1]], "contexts": [], "labels": [1]}
    variants = [base, {**base, "contexts": [{"x": 1}]}, {**base, "labels": [2]}]
    for work in variants:
        results, pool = verifier._execute_batch([("", work)], Limits(), False, 1, 64, None)
        assert results[0]["children"][0] == {**work, "depth": 1} and pool is None
    verifier._execute_batch([("", base)], Limits(max_steps=1), False, 1, 64, None)
    verifier._execute_batch([("", base)], Limits(), True, 1, 64, None)
    assert verifier.cache_info()["misses"] == len(variants) + 2
    changed = {**base, "depth": 19, "coefficient": [[-7, -(2**32771)], [2, 3]]}
    result, _ = verifier._execute_batch([("", changed)], Limits(), False, 1, 64, None)
    assert result[0]["children"][0]["coefficient"] == changed["coefficient"]
    assert result[0]["children"][0]["depth"] == 20
    result[0]["children"][0]["labels"].append(99)
    again, _ = verifier._execute_batch([("", base)], Limits(), False, 1, 64, None)
    assert again[0]["children"][0]["labels"] == [1]
    assert verifier.cache_info()["hits"] == 2


def test_native_default_explicit_and_cold_paths():
    session, result = certificate()
    assert session.check(result)
    assert session.check(result)
    assert session._verifier.cache_info()["hits"] > 0
    assert session.check(result, verifier=independent_replay.IndependentVerifier())
    assert session.check(result, cold=True)
    for kwargs in ({"cold": 1}, {"cold": True, "verifier": independent_replay.IndependentVerifier()}):
        with pytest.raises(native.NativeError):
            session.check(result, **kwargs)


def test_frozen_rules_match_and_cannot_be_mutated():
    plain = collected.Catalog.bundled()
    frozen = independent_replay._ReplayCatalog.bundled()
    for rule in plain.rules:
        key = rule["id"]
        candidate = frozen.get(key)
        assert native._strict_equal(candidate, plain.get(key))
        candidate["lhs"][0]["coefficient"][0]["integer"] = 999999
        assert native._strict_equal(frozen.get(key), plain.get(key))
    env = {"n": 4, "i": 2, "j": 1, "k": 0, "l": 2, "eps": 1}
    assert native._strict_equal(frozen.select_rules(env), plain.select_rules(env))
    with pytest.raises(independent_replay.CatalogError, match="Unknown rule"):
        frozen.get("not-a-rule")


def test_malformed_and_refused_replay():
    session, result = certificate()
    verifier = independent_replay.IndependentVerifier()
    del result["limits"]
    with pytest.raises(native.NativeError, match="Malformed"):
        verifier.check_native(session, result)
    initial = collected.reference._single({
        "rank": 3, "coefficient": [[0, 1]], "contexts": [], "world_phi": False,
        "boundary_left": [], "boundary_right": [], "boundary_left_factors": [],
        "boundary_right_factors": [], "depth": 4, "bottom_labels": [1, 2, 3],
        "top_labels": [1, 2, 3]}, [1, 1])
    for limits in (Limits(expansion_budget=1), Limits(max_integer_bits=1), Limits()):
        for trace in (True, False):
            for work in (initial, {**initial, "coefficient": [[0, 2**1000]], "depth": 17}):
                expected = collected._transition(work, limits, trace, verifier._catalog)
                for _ in range(2):
                    result, _ = verifier._execute_batch([("key", work)], limits, trace, 1, 64, None)
                    assert native._strict_equal(result[0], expected)


def test_type_checking_imports_resolve_without_building(monkeypatch):
    import runpy
    import typing

    monkeypatch.setattr(typing, "TYPE_CHECKING", True)
    for module in (independent_replay, native):
        result = runpy.run_path(module.__file__, run_name=module.__name__)
        assert "IndependentVerifier" in result


def test_unit_cache_all_real_phases_exactly_matches_unfactored_python():
    seen = {}

    class Capture(collected.CollectedSession):
        def _execute_batch(self, selected, limits, trace, workers, batch_size, pool):
            for _, work in selected:
                seen.setdefault(work["phase"], deepcopy(work))
            return super()._execute_batch(selected, limits, trace, workers, batch_size, pool)

    Capture().reduce([1, -2, 1], n=3, trace=True, batch_size=4)
    assert set(seen) == {"single", "zip_up", "zip_down", "layer", "extract"}
    verifier = independent_replay.IndependentVerifier()
    for base in seen.values():
        for trace in (False, True):
            for limits in (Limits(), Limits(max_integer_bits=1), Limits(expansion_budget=1)):
                for coefficient, depth in (([[0, 1]], 0), ([[-7, -(2**32771)], [2, 3]], 17)):
                    work = deepcopy(base)
                    collected._set_coefficient(work, coefficient)
                    work["depth"] = depth
                    expected = collected._transition(work, limits, trace, verifier._catalog)
                    for _ in range(2):
                        actual, _ = verifier._execute_batch([("", work)], limits, trace, 1, 64, None)
                        assert native._strict_equal(actual[0], expected)


def test_cold_audit_rejects_missing_and_changed_math():
    session, result = certificate()
    missing = deepcopy(result)
    del missing["limits"]
    with pytest.raises(native.NativeError, match="Malformed"):
        session.check(missing, cold=True)
    result["output"][0]["coefficient"][0][1] += 1
    with pytest.raises(native.NativeError, match="differs from independent"):
        session.check(result, cold=True)


@pytest.mark.parametrize("operation", ["clear_cache", "cache_info"])
def test_cache_operations_wait_for_active_verification(operation):
    from concurrent.futures import ThreadPoolExecutor
    from threading import Event, Lock

    verifier = independent_replay.IndependentVerifier()
    entered = Event()
    lock = Lock()

    class ObservedLock:
        def __enter__(self):
            entered.set()
            lock.acquire()

        def __exit__(self, *args):
            lock.release()

    verifier._lock = ObservedLock()
    with ThreadPoolExecutor(max_workers=1) as pool:
        with lock:
            future = pool.submit(getattr(verifier, operation))
            assert entered.wait(5)
            assert not future.done()
        future.result(timeout=5)


@pytest.mark.parametrize("value", [None, [], "certificate", 1])
def test_non_object_certificate_rejected(value):
    with pytest.raises(native.NativeError):
        fake_session().check(value)


def test_wrapper_comparator_is_bound_in_verifier_identity():
    import hashlib
    from pathlib import Path

    assert independent_replay.source_identity()["native.py"] == hashlib.sha256(
        Path(native.__file__).read_bytes()).hexdigest()
