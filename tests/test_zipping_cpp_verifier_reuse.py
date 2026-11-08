"""Independent completion replay, context restoration, and tamper resistance."""
from copy import deepcopy
from dataclasses import asdict
from threading import Lock

import pytest

from zipper_algebra.zipping_algorithm import collected, native
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier, VerifierError, VerifierIncomplete

pytestmark = pytest.mark.timeout(300)

EXECUTION = {"schedule": "intrinsic-completion-dag-v1",
             "collection": "full-output-except-scalar-v1",
             "local_term_floor": 4, "local_support_floor": 2}


def reused_envelope(certificate):
    result = deepcopy(certificate)
    result.pop("transitions")
    result["schema"] = "zipper-native-reuse-v1"
    result["scope"] = "literal_word_native_reused_spanning_procedure"
    result["execution"] = EXECUTION.copy()
    return result


@pytest.fixture(scope="module")
def material():
    # Standalone expected outputs use the original stable-prefix scheduler.
    with native.NativeSession() as producer:
        cases = []
        for n, word in ((3, [1, -2, 1]), (4, [1, 3, -2]),
                        (5, [4, -2, 1, -4]), (8, [-7, 1, 7, 2])):
            cert = producer.reduce(word, n=n, limits=Limits(max_steps=100000))
            assert cert["status"] == "complete"
            cases.append(cert)
        yield producer.metadata, cases


@pytest.mark.parametrize("cache_entries,cache_bytes", [(0, 0), (1, 10000000), (10000, 100000000)])
def test_completed_replay_matches_independent_standalone_outputs(material, cache_entries, cache_bytes):
    metadata, cases = material
    with CppVerifier(producer_metadata=metadata, cache_entries=cache_entries,
                     cache_bytes=cache_bytes) as verifier:
        for _ in range(2):
            for cert in cases:
                assert verifier.check_reused(reused_envelope(cert))
                assert verifier.last_result["verification_mode"] == "independent-intrinsic-completion-v1"
        info = verifier.cache_info()["completed"]
        assert info["entries"] <= cache_entries
        assert info["bytes"] <= cache_bytes
        if cache_entries == 0:
            assert info["entries"] == info["hits"] == 0
        elif cache_entries == 1:
            assert info["evictions"] > 0
        else:
            assert info["cross_input_hits"] > 0
            assert verifier.last_result["metrics"]["expanded_nodes"] == 0
        verifier.clear_cache()
        assert verifier.cache_info()["completed"]["entries"] == 0
        assert verifier.check_reused(reused_envelope(cases[0]))


@pytest.mark.parametrize("field", ["coefficient", "ports", "word", "status", "pending", "trace",
                                   "schema", "scope", "source", "limits", "execution", "extra", "type"])
def test_reuse_rejects_changed_math_contract_and_provenance(material, field):
    metadata, cases = material
    good = reused_envelope(cases[0])
    with CppVerifier(producer_metadata=metadata) as verifier:
        assert verifier.check_reused(good)
        bad = deepcopy(good)
        if field == "coefficient":
            bad["output"][0]["coefficient"][0][1] += 1
        elif field == "ports":
            bad["output"][0]["components"][0]["record"]["original_bottom_ports"][0] += 1
        elif field == "word":
            bad["output"][0]["word"].append(1)
        elif field == "status":
            bad["status"] = "incomplete"
        elif field == "pending":
            bad["pending"] = [{}]
        elif field == "trace":
            bad["trace_enabled"] = True
        elif field == "schema":
            bad["schema"] = "zipper-native-run-v1"
        elif field == "scope":
            bad["scope"] = "unproved_quotient"
        elif field == "source":
            bad["catalog_sha256"] = "0" * 64
        elif field == "limits":
            bad["limits"]["max_steps"] = 1
        elif field == "execution":
            bad["execution"]["schedule"] = "stable-prefix-atomic-v1"
        elif field == "extra":
            bad["extra"] = True
        else:
            bad["n"] = float(bad["n"])
        with pytest.raises(VerifierError) as refused:
            verifier.check_reused(bad)
        assert isinstance(refused.value, VerifierIncomplete) == (field == "limits")
        assert verifier.check_reused(good)


def test_output_order_and_exact_collection_not_producer_telemetry(material):
    metadata, cases = material
    cert = reused_envelope(cases[0])
    cert["output"].reverse()
    term = cert["output"][0]
    duplicate = deepcopy(term)
    duplicate["coefficient"] = [[q, -c] for q, c in term["coefficient"]]
    term["coefficient"] = [[q, 2 * c] for q, c in term["coefficient"]]
    cert["output"].append(duplicate)
    cert["metrics"] = {"untrusted_telemetry": True}
    with CppVerifier(producer_metadata=metadata) as verifier:
        assert verifier.check_reused(cert)


def test_reuse_does_not_call_producer_or_python_math(material, monkeypatch):
    metadata, cases = material

    def forbidden(*args, **kwargs):
        raise AssertionError("Independent completion called another mathematical implementation")

    with CppVerifier(producer_metadata=metadata) as verifier:
        monkeypatch.setattr(native.NativeSession, "reduce", forbidden)
        monkeypatch.setattr(native.NativeSession, "_request", forbidden)
        monkeypatch.setattr(collected.CollectedSession, "reduce", forbidden)
        monkeypatch.setattr(collected, "_transition", forbidden)
        assert verifier.check_reused(reused_envelope(cases[1]))


def test_incomplete_root_is_never_retained_as_a_complete_result(material):
    metadata, cases = material
    with CppVerifier(producer_metadata=metadata) as verifier:
        for _ in range(2):
            result = verifier.diagnostic("replay_reused", n=3, input_word=[1, -2, 1],
                                         limits=asdict(Limits(max_steps=1)))
            assert result["status"] == "incomplete"
            assert result["output"] == []
            assert result["metrics"]["expanded_nodes"] == 1
        assert verifier.check_reused(reused_envelope(cases[0]))


def test_full_replay_context_normalization_keeps_certificate_exact(material):
    metadata, cases = material
    with CppVerifier(producer_metadata=metadata) as verifier:
        for cert in cases:
            assert verifier.check(cert)
        assert verifier.cache_info()["cross_input_hits"] > 0


def test_cold_verifier_resource_refusal_is_distinct_from_wrong_output():
    limits = Limits(max_steps=2)
    with native.NativeSession() as producer:
        seed = producer.reduce_reused([], n=1, limits=limits)
        warm = producer.reduce_reused([], n=2, limits=limits)
        assert seed["status"] == warm["status"] == "complete"
        with CppVerifier(producer_metadata=producer.metadata) as verifier:
            with pytest.raises(VerifierIncomplete) as refusal:
                verifier.check_reused(warm)
            assert refusal.value.replay_status == "incomplete"
            assert "limit" in refusal.value.replay_diagnostic["message"]
            assert verifier.last_result["replay_status"] == "incomplete"
            assert verifier.check_reused(seed)
            assert verifier.check_reused(warm)
            changed = deepcopy(warm)
            changed["output"][0]["coefficient"][0][1] += 1
            with pytest.raises(VerifierError) as rejected:
                verifier.check_reused(changed)
            assert not isinstance(rejected.value, VerifierIncomplete)
            assert verifier.last_result["replay_status"] == "complete"


@pytest.mark.parametrize("n,word", [(0, []), (3, [0]), (3, [3]), (3, [1.0]), (3, "word")])
def test_reused_replay_validates_literal_input(material, n, word):
    metadata, _ = material
    with CppVerifier(producer_metadata=metadata) as verifier:
        with pytest.raises(VerifierError):
            verifier.diagnostic("replay_reused", n=n, input_word=word, limits=asdict(Limits()))


def bare_verifier():
    verifier = object.__new__(CppVerifier)
    verifier._lock = Lock()
    verifier._process = None
    verifier._stderr = None
    verifier.last_result = {"ok": True, "replay_status": "complete"}
    verifier.last_call_seconds = 9.0
    return verifier


@pytest.mark.parametrize("method", ["check", "check_reused"])
@pytest.mark.parametrize("failure", ["transport", "validation", "source"])
def test_failed_check_clears_previous_success(monkeypatch, method, failure):
    verifier = bare_verifier()

    def fail(*args, **kwargs):
        raise VerifierError("failed request")

    monkeypatch.setattr(verifier, "_assert_source", fail if failure == "source" else lambda: None)
    monkeypatch.setattr(verifier, "_request", fail)
    with pytest.raises(VerifierError):
        getattr(verifier, method)(None if failure == "validation" else {})
    assert verifier.last_result == {}
    assert verifier.last_call_seconds == 0.0


@pytest.mark.parametrize("status", ["incomplete", "control_obligation", "unsupported", "complete", "unknown"])
def test_only_unresolved_replay_statuses_are_retryable(monkeypatch, status):
    verifier = bare_verifier()
    monkeypatch.setattr(verifier, "_assert_source", lambda: None)
    monkeypatch.setattr(verifier, "_request", lambda _: {
        "ok": False, "error": "refused", "replay_status": status,
        "replay_diagnostic": {"message": "reason"}})
    with pytest.raises(VerifierError) as error:
        verifier.check_reused({})
    assert isinstance(error.value, VerifierIncomplete) == (status in {
        "incomplete", "control_obligation", "unsupported"})
    if isinstance(error.value, VerifierIncomplete):
        assert error.value.replay_status == status
        assert error.value.replay_diagnostic == {"message": "reason"}


def test_full_verifier_resource_refusal_preserves_diagnostic(material, monkeypatch):
    """Cold full replay can exhaust resources independently of its producer."""
    verifier = bare_verifier()
    monkeypatch.setattr(verifier, '_assert_source', lambda: None)
    monkeypatch.setattr(verifier, '_request', lambda _: {
        'ok': False, 'replay_status': 'incomplete', 'error': 'test replay budget',
        'replay_diagnostic': {'message': 'coefficient limit'}})
    with pytest.raises(VerifierIncomplete, match='test replay budget') as error:
        verifier.check(material[1][0])
    assert error.value.replay_diagnostic == {'message': 'coefficient limit'}
    assert verifier.last_result['ok'] is False
