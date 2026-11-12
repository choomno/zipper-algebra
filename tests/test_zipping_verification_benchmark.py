"""Storage and exact-input contracts of independent verification benchmarking."""
import hashlib
import json
import sqlite3
import zlib

import pytest

from experiments import zipper_verification_benchmark as benchmark
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import _decode_coefficients, NativeSession


@pytest.mark.parametrize("cases", [[], {}, [{"n": True, "word": []}], [{"n": 3, "word": [3]}],
                                   [{"n": 3, "word": [True]}], [{"n": 3, "word": [1.0]}],
                                   [{"n": 2, "word": []}] * 2])
def test_reject_invalid_or_repeated_inputs(cases):
    with pytest.raises(ValueError):
        benchmark.validate_cases(cases)


def test_sqlite_retains_exact_certificate_and_separate_timings(tmp_path):
    database = tmp_path / "measurements.sqlite3"
    campaign = benchmark.run(database, [{"n": 2, "word": [1, -1]}], limits=Limits(max_steps=64))
    with sqlite3.connect(database) as db:
        assert db.execute("SELECT status FROM verification_campaigns").fetchone() == ("complete",)
        row = db.execute("SELECT certificate,certificate_sha256,producer_seconds,cold_seconds,retained_seconds "
                         "FROM verification_measurements WHERE campaign=?", (campaign,)).fetchone()
        assert db.execute("PRAGMA foreign_key_check").fetchall() == []
    payload = zlib.decompress(row[0])
    assert hashlib.sha256(payload).hexdigest() == row[1]
    assert all(value >= 0 for value in row[2:])
    certificate = _decode_coefficients(json.loads(payload))
    assert certificate["input_word"] == [1, -1]
    with NativeSession() as session:
        assert session.check(certificate, cold=True)
