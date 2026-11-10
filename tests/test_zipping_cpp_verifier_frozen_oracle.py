"""Differential contract tests against the frozen independent JSON verifier.

The producer supplies authenticated build metadata only. All expected transitions,
refusals, traces, outputs and cache accounting come from revision 9c477dd.
"""
from collections import deque
from contextlib import ExitStack
from copy import deepcopy
from dataclasses import asdict
import hashlib
import io
import json
from pathlib import Path
import tarfile
import subprocess
from types import ModuleType

import pytest

from zipper_algebra.zipping_algorithm import native, native_verifier
from zipper_algebra.zipping_algorithm.interpreter import Limits

pytestmark = pytest.mark.timeout(600)
_FIXTURE = Path(__file__).parent / "native" / "frozen_verifier_9c477dd"
_PHASES = {"single", "zip_up", "zip_down", "layer", "extract"}


def _frozen_module(destination):
    """Authenticate every frozen input before compiling a separate executable."""
    manifest = json.loads(_FIXTURE.with_suffix(".json").read_text())
    payload = _FIXTURE.with_suffix(".tar.gz").read_bytes()
    assert manifest["archive_sha256"] == (
        "7c41dd88004a8f39b595e32b738faa5c54fd50b04b1544620031d648274266ba")
    assert hashlib.sha256(payload).hexdigest() == manifest["archive_sha256"]
    assert manifest["commit"] == "9c477dd404738be0a7133bef75b2d4adc37c2616"
    members = {}
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                assert member.isdir()
                continue
            assert member.name in manifest["files"]
            data = archive.extractfile(member).read()
            assert hashlib.sha256(data).hexdigest() == manifest["files"][member.name]
            target = destination / member.name
            assert target.resolve().is_relative_to(destination.resolve())
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
            members[member.name] = True
    assert members.keys() == manifest["files"].keys()
    # Execute the transport wrapper in a separate global namespace. Its source,
    # producer and reference bindings remain authenticated; only the independent
    # verifier source inventory points to the frozen snapshot.
    wrapper = ModuleType(native_verifier.__package__ + ".frozen_oracle")
    wrapper.__package__ = native_verifier.__package__
    wrapper.__file__ = native_verifier.__file__
    exec(compile(Path(wrapper.__file__).read_text(), wrapper.__file__, "exec"), wrapper.__dict__)
    wrapper.SOURCE_DIRECTORY = destination / "native" / "zipper_verifier"
    return wrapper


@pytest.fixture(scope="module")
def oracle_material(tmp_path_factory):
    directory = tmp_path_factory.mktemp("frozen-json-verifier")
    wrapper = _frozen_module(directory)
    with native.NativeSession() as producer:
        yield wrapper, producer.metadata, directory / "build"


@pytest.fixture(scope="module")
def oracle(oracle_material):
    wrapper, metadata, build = oracle_material
    with wrapper.CppVerifier(producer_metadata=metadata, build_directory=build) as session:
        yield session


@pytest.fixture(scope="module")
def candidate(oracle_material):
    _, metadata, _ = oracle_material
    with native_verifier.CppVerifier(producer_metadata=metadata) as session:
        yield session


def _single(rank, word):
    return {"rank": rank, "coefficient": [[0, 1]], "contexts": [],
            "world_phi": False, "boundary_left": [], "boundary_right": [],
            "boundary_left_factors": [], "boundary_right_factors": [],
            "depth": 0, "bottom_labels": list(range(1, rank + 1)),
            "top_labels": list(range(1, rank + 1)), "phase": "single",
            "left": list(word), "processed": [], "position": rank}


def _request(session, op, **parameters):
    # Preserve the entire native refusal object, including diagnostic precedence.
    return session._request(dict(op=op, **native_verifier._wire(parameters)))


def _first_difference(actual, expected, path=""):
    if type(actual) is not type(expected):
        return path, type(actual).__name__, type(expected).__name__
    if isinstance(actual, dict):
        if actual.keys() != expected.keys():
            return path, sorted(actual.keys() - expected.keys()), sorted(expected.keys() - actual.keys())
        for key in actual:
            difference = _first_difference(actual[key], expected[key], path + "/" + str(key))
            if difference is not None:
                return difference
    elif isinstance(actual, list):
        if len(actual) != len(expected):
            return path + "/length", len(actual), len(expected)
        for index, (left, right) in enumerate(zip(actual, expected)):
            difference = _first_difference(left, right, path + "/" + str(index))
            if difference is not None:
                return difference
    elif actual != expected:
        return path, actual, expected
    return None


def _equal(actual, expected, label):
    difference = _first_difference(actual, expected)
    assert difference is None, (label, difference)


_SEED_DRIVER = r"""
#include "algorithm.hpp"
#include "controllers.hpp"
#include <iostream>
using namespace zipper_verify;
int main() {
    J result=J::array();
    for(I n:{5,9}) for(I size:{3,4}) for(bool mirror:{false,true})
        for(const auto &[out,in]:std::vector<std::pair<I,I>>{{0,0},{1,-1},{-2,2}}) {
            I start=2;
            J state=framed_initialize(n,J::array({1,-(n-1)}),size,start,0,mirror,
                         out,in,J::array({1}),J::array({-1}),J::array(),"bar_mirror",Limits());
            J work=initial_work(n,J::array());
            work.erase("left");work.erase("processed");work.erase("position");
            work["phase"]="zip_up";work["branch"]=state.at("branches").at(0);
            result.push_back(std::move(work));
        }
    std::cout<<result<<'\n';
}
"""


@pytest.fixture(scope="module")
def framed_seeds(oracle_material):
    wrapper, _, build = oracle_material
    # All objects are authenticated by the frozen wrapper build. The seed
    # generator links that complete legacy implementation, never current helpers.
    _, metadata, _ = wrapper.build_verifier(build_directory=build)
    source = build.parent / "frozen-seeds.cpp"
    source.write_text(_SEED_DRIVER)
    binary = build.parent / "frozen-seeds"
    objects = []
    for name in sorted(wrapper.SOURCE_DIRECTORY.glob("*.cpp")):
        if name.name != "main.cpp":
            matches = list((build / "objects").rglob(name.stem + ".o"))
            assert len(matches) == 1
            objects.extend(matches)
    subprocess.run([metadata["compiler"], *wrapper.FLAGS, "-I", str(wrapper.SOURCE_DIRECTORY),
                    str(source), *map(str, objects),
                    *("-l" + name for name in wrapper.LIBRARIES), "-o", str(binary)],
                   check=True, capture_output=True, text=True)
    seeds = json.loads(subprocess.run([str(binary)], check=True, capture_output=True, text=True).stdout)
    assert len(seeds) == 24
    return seeds


def _geometry(work):
    control = work.get("branch", {}).get("control", {})
    return tuple(control.get(key) for key in ("i", "mirrored", "output_power", "input_power", "phase")) + (
        len(work["contexts"]), work["world_phi"])


@pytest.fixture(scope="module")
def reached_states(oracle, framed_seeds):
    """Explore only the frozen controller; retain rule and geometry diversity."""
    starts = [(3, [1, -2, 1]), (4, [-3, -2, -1]),
              (5, [4, -3, 2, -1]), (5, [-4, -3, 1, 2, 3, 4])]
    queue = deque([*(_single(n, word) for n, word in starts), *framed_seeds])
    seen = set()
    selected = {}
    for _ in range(1800):
        if not queue:
            break
        work = queue.popleft()
        identity = json.dumps(work, sort_keys=True)
        if identity in seen:
            continue
        seen.add(identity)
        result = _request(oracle, "transition", work=work, limits=asdict(Limits()), trace=True)
        assert result.get("ok") is not False, result
        key = (work["phase"], result.get("unit_evidence", {}).get("rule_id"), _geometry(work))
        selected.setdefault(key, work)
        if "children" in result:
            queue.extend(result["children"])
    states = list(selected.values())
    assert {work["phase"] for work in states} == _PHASES
    controls = [work["branch"]["control"] for work in states if "branch" in work]
    assert {control["i"] for control in controls} >= {2, 3, 4}
    assert {control.get("mirrored") for control in controls} >= {False, True}
    assert {control.get("output_power") for control in controls} >= {-2, 0, 1}
    assert {control.get("input_power") for control in controls} >= {-1, 0, 2}
    return states


@pytest.mark.parametrize("trace", [False, True])
def test_reached_phase_transitions_equal_frozen_json(oracle, candidate, reached_states, trace):
    for original in reached_states:
        work = deepcopy(original)
        work["depth"] = 7
        coefficient = [[-2, 2**256 + 3], [1, -1]]
        work["coefficient"] = coefficient
        if "branch" in work:
            work["branch"]["coefficient"] = coefficient
        parameters = dict(work=work, limits=asdict(Limits()), trace=trace)
        _equal(_request(candidate, "typed_transition", **parameters),
               _request(oracle, "transition", **parameters), (work["phase"], trace))


def test_legacy_phase_malformed_fields_and_resource_refusals(oracle, candidate, reached_states):
    by_phase = {work["phase"]: work for work in reached_states}
    for phase, original in by_phase.items():
        mutations = []
        for field in ("rank", "coefficient", "contexts", "depth", "bottom_labels", "top_labels"):
            absent = deepcopy(original)
            absent.pop(field)
            mutations.append(("missing " + field, absent))
        for field, value in (("rank", True), ("depth", 2**63), ("contexts", {}),
                             ("coefficient", [[0, 0]]), ("coefficient", [[0, 1], [0, 2]]),
                             ("world_phi", "false"), ("bottom_labels", [True]),
                             ("boundary_left_factors", [None]), ("coefficient", []),
                             ("extension", {"unicode": "λ", "value": [True, None]})):
            bad = deepcopy(original)
            bad[field] = value
            mutations.append((field, bad))
        # JSON control records can retain extension data and skip fields on a
        # zero branch. These cases check the explicit legacy compatibility API.
        zero = deepcopy(original)
        zero["coefficient"] = []
        zero["contexts"] = [None]
        if "branch" in zero:
            zero["branch"]["coefficient"] = []
            zero["branch"]["extension"] = {"unused": [None, True]}
        mutations.append(("zero scalar with dead malformed context", zero))
        for key in ("branch", "control"):
            if isinstance(original.get(key), dict):
                extra = deepcopy(original)
                extra[key]["extension"] = {"unicode": "λ"}
                mutations.append((key + " extension", extra))
        if isinstance(original.get("branch"), dict):
            for key in ("control", "components", "record"):
                if isinstance(original["branch"].get(key), dict):
                    extra = deepcopy(original)
                    extra["branch"][key]["extension"] = [True, None]
                    mutations.append(("branch " + key + " extension", extra))
        for label, work in mutations:
            parameters = dict(work=work, limits=asdict(Limits()), trace=True)
            _equal(_request(candidate, "transition", **parameters),
                   _request(oracle, "transition", **parameters), (phase, label))
        for bound in ("max_word_size", "max_terms", "max_coefficient_support",
                      "max_integer_bits", "expansion_budget"):
            limits = asdict(Limits())
            limits[bound] = 0
            parameters = dict(work=original, limits=limits, trace=False)
            _equal(_request(candidate, "typed_transition", **parameters),
                   _request(oracle, "transition", **parameters), (phase, bound))


@pytest.mark.parametrize("n,word", [(1, []), (2, [1, -1]), (3, [1, -2, 1]),
                                    (4, [1, 3, -2]), (8, [-7, 1, 7, 2]),
                                    (16, [15, -1]), (33, [-32, 1])])
@pytest.mark.parametrize("trace,batch", [(False, 1), (True, 1), (True, 4)])
def test_standalone_output_trace_atomic_batches_and_arbitrary_rank(oracle, candidate, n, word, trace, batch):
    parameters = dict(n=n, input_word=word, limits=asdict(Limits()),
                      batch_size=batch, trace_enabled=trace)
    _equal(_request(candidate, "replay", **parameters),
           _request(oracle, "replay", **parameters), (n, word, trace, batch))


def _all_differences(actual, expected, path=""):
    if type(actual) is not type(expected):
        yield path, actual, expected
    elif isinstance(actual, dict):
        for key in actual.keys() | expected.keys():
            if key not in actual or key not in expected:
                yield path + "/" + str(key), actual.get(key), expected.get(key)
            else:
                yield from _all_differences(actual[key], expected[key], path + "/" + str(key))
    elif isinstance(actual, list):
        if len(actual) != len(expected):
            yield path + "/length", len(actual), len(expected)
        for index, (left, right) in enumerate(zip(actual, expected)):
            yield from _all_differences(left, right, path + "/" + str(index))
    elif actual != expected:
        yield path, actual, expected


def _equal_completion(actual, expected, label, raw_differences, record_path, prefix="/cache"):
    differences = list(_all_differences(actual, expected))
    if differences:
        raw_differences.append({"case": label, "differences": differences})
        record_path.write_text(json.dumps(raw_differences, indent=2) + "\n")
    # The typed chosen-rule path holds its independently checked bound handle.
    # The JSON core.apply_rule path performs a second lookup. Omitting that
    # redundant lookup changes only these internal bound-cache call counters in
    # this corpus. All selection, residency, eviction, and completion fields must
    # still agree exactly with the frozen implementation.
    allowed = {prefix + "/bound_rules/" + field for field in ("hits", "misses", "oversized")}
    unexpected = [difference for difference in differences if difference[0] not in allowed]
    assert not unexpected, (label, unexpected[:3])


@pytest.mark.parametrize("cache_entries,cache_bytes", [(0, 0), (1, 10000000),
                                                        (10000, 4096), (10000, 10000000)])
def test_fused_completion_cache_admission_eviction_and_limits(oracle_material, cache_entries, cache_bytes, tmp_path):
    wrapper, metadata, build = oracle_material
    raw_differences = []
    record_path = tmp_path / "frozen-native-differences.json"
    with ExitStack() as stack:
        old = stack.enter_context(wrapper.CppVerifier(producer_metadata=metadata, build_directory=build,
                                  cache_entries=cache_entries, cache_bytes=cache_bytes))
        new = stack.enter_context(native_verifier.CppVerifier(producer_metadata=metadata,
                                  cache_entries=cache_entries, cache_bytes=cache_bytes))
        cases = [(3, [1, -2, 1]), (4, [1, 3, -2]), (3, [-1, 2, -1]), (3, [1, -2, 1])]
        settings = [asdict(Limits())]
        for field, value in (("max_steps", 0), ("max_steps", 1), ("max_steps", 8),
                             ("max_terms", 1), ("max_word_size", 2),
                             ("max_coefficient_support", 1), ("max_integer_bits", 1)):
            limits = asdict(Limits())
            limits[field] = value
            settings.append(limits)
        for limits in settings:
            for n, word in cases:
                parameters = dict(n=n, input_word=word, limits=limits)
                _equal_completion(_request(new, "replay_reused", **parameters),
                       _request(old, "replay_reused", **parameters), (cache_entries, limits, n, word),
                       raw_differences, record_path)
                _equal_completion(new.cache_info(), old.cache_info(),
                       (cache_entries, "cache accounting"), raw_differences, record_path, prefix="")
        new.clear_cache()
        old.clear_cache()
        _equal_completion(new.cache_info(), old.cache_info(), "clear cache",
                          raw_differences, record_path, prefix="")


@pytest.mark.parametrize("n,word", [(0, []), (True, []), (3, [True]), (3, [0]),
                                    (3, [3]), (3, [-(2**63)]),
                                    (3, [["S", 1, 2]]), (3, [1.0]),
                                    (0, [0]), (0, [["S", 1, 2]]),
                                    (3, {}), (3, None), (3, True)])
def test_public_literal_input_refusals_equal_frozen(oracle, candidate, n, word):
    for op in ("replay", "replay_reused"):
        parameters = dict(n=n, input_word=word, limits=asdict(Limits()),
                          batch_size=4, trace_enabled=True)
        _equal(_request(candidate, op, **parameters),
               _request(oracle, op, **parameters), (op, n, word))


def test_typed_diagnostic_refuses_missing_mandatory_state_fields(candidate, reached_states):
    # The legacy diagnostic has a wider domain. These checks concern the typed
    # diagnostic's own mandatory state schema, without requiring legacy errors.
    by_phase = {work["phase"]: work for work in reached_states}
    for phase, original in by_phase.items():
        fields = {"rank", "coefficient", "contexts", "depth", "bottom_labels", "top_labels"}
        fields |= {"single": {"left", "processed", "position"},
                   "zip_up": {"branch"}, "zip_down": {"branch"},
                   "layer": {"word", "payload_length", "completed_sides"},
                   "extract": {"word"}}[phase]
        for field in fields:
            work = deepcopy(original)
            del work[field]
            result = _request(candidate, "typed_transition", work=work,
                              limits=asdict(Limits()), trace=True)
            assert result.get("ok") is False, (phase, field, result)


def test_typed_diagnostic_rejects_unknown_owned_record_keys(candidate, reached_states):
    cases = []
    for phase in _PHASES:
        work = deepcopy(next(work for work in reached_states if work["phase"] == phase))
        work["extension"] = {"unicode": "λ", "value": [True, None]}
        cases.append((phase + " work", work))
    branch_work = next(work for work in reached_states if "branch" in work)
    for path in (("branch",), ("branch", "control")):
        work = deepcopy(branch_work)
        record = work
        for key in path:
            record = record[key]
        record["extension"] = {"unused": [None, True]}
        cases.append(("/".join(path), work))
    context_work = next(work for work in reached_states if work["contexts"])
    for path in (("contexts", 0), ("contexts", 0, "record"),
                 ("contexts", 0, "record", "canonical_parities")):
        work = deepcopy(context_work)
        record = work
        for key in path:
            record = record[key]
        record["extension"] = []
        cases.append(("/".join(map(str, path)), work))
    for label, work in cases:
        result = _request(candidate, "typed_transition", work=work,
                          limits=asdict(Limits()), trace=False)
        assert result.get("ok") is False, (label, result)



def test_public_typed_diagnostic_uses_real_independent_binary(candidate, oracle):
    parameters = dict(work=_single(3, [1, -2, 1]), limits=asdict(Limits()),
                      trace_enabled=True)
    _equal(candidate.diagnostic("typed_transition", **parameters),
           oracle.diagnostic("transition", **parameters), "public typed diagnostic")
