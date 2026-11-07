"""Compare independent phase replay with the Python specification on small states."""

from collections import deque
from copy import deepcopy
from dataclasses import asdict

import pytest

from zipper_algebra.zipping_algorithm import algorithm, collected, native
from zipper_algebra.zipping_algorithm.catalog import Catalog
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier, _wire

pytestmark = pytest.mark.timeout(300)
_PHASES = {"single", "zip_up", "zip_down", "layer", "extract"}


def _common(rank):
    return {
        "rank": rank,
        "coefficient": [[0, 1]],
        "contexts": [],
        "world_phi": False,
        "boundary_left": [],
        "boundary_right": [],
        "boundary_left_factors": [],
        "boundary_right_factors": [],
        "depth": 0,
        "bottom_labels": list(range(1, rank + 1)),
        "top_labels": list(range(1, rank + 1)),
    }


@pytest.fixture(scope="module")
def phase_cases():
    """Keep one reached state per phase and starting geometry, without saved data."""
    limits = Limits()
    catalog = Catalog.bundled()
    starts = [("literal", algorithm._single(_common(3), [1, -2, 1]))]
    for mirrored in (False, True):
        starts.append((
            "mirrored" if mirrored else "ordinary",
            algorithm._new_zip(
                _common(5), left=[], right=[], size=3, start=2,
                mirrored=mirrored, limits=limits,
            ),
        ))
    cases = []
    for geometry, initial in starts:
        queue = deque([initial])
        phases = set()
        for _ in range(48):
            if not queue:
                break
            work = queue.popleft()
            if work["phase"] not in phases:
                cases.append((geometry, work))
                phases.add(work["phase"])
            result = collected._transition(work, limits, False, catalog)
            assert "error" not in result
            queue.extend(result["children"])
        if geometry != "literal":
            assert {"zip_up", "zip_down", "layer", "extract"} <= phases
    assert {work["phase"] for _, work in cases} == _PHASES
    return cases, limits, catalog


@pytest.fixture(scope="module")
def verifier():
    with native.NativeSession() as producer:
        with CppVerifier(producer_metadata=producer.metadata) as session:
            yield session


@pytest.mark.parametrize("trace", [False, True], ids=["compact", "full-trace"])
@pytest.mark.parametrize("wide", [False, True], ids=["unit", "wide-laurent"])
def test_all_phase_transitions_match_python(verifier, phase_cases, trace, wide):
    cases, limits, catalog = phase_cases
    for geometry, original in cases:
        work = deepcopy(original)
        work["depth"] = 7
        coefficient = [[-2, 2**2048 + 3], [1, -1]] if wide else [[0, 1]]
        collected._set_coefficient(work, coefficient)
        expected = collected._transition(work, limits, trace, catalog)
        actual = verifier.diagnostic(
            "transition", work=work, limits=asdict(limits), trace=trace,
        )
        assert native._strict_equal(_wire(expected), actual), (geometry, work["phase"])
        assert all(child["depth"] == 8 for child in actual["children"])
