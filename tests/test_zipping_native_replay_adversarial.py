"""Independent acceptance tests that deny the verifier access to producer work."""
from copy import deepcopy

import pytest

from zipper_algebra.zipping_algorithm.independent_replay import IndependentVerifier
from zipper_algebra.zipping_algorithm.native import NativeError, NativeSession

pytestmark = pytest.mark.timeout(300)


def test_verification_uses_no_native_callbacks_and_survives_candidate_corruption(monkeypatch):
    with NativeSession() as producer:
        certificate = producer.reduce([1, -2, 1], n=3, trace=True, batch_size=4)
        saved = deepcopy(certificate)
        verifier = IndependentVerifier()

        def forbidden(*args, **kwargs):
            raise AssertionError("independent verification called its producer")

        # Deny both full native reductions and direct native protocol operations.
        monkeypatch.setattr(producer, "reduce", forbidden)
        monkeypatch.setattr(producer, "_request", forbidden)
        assert producer.check(certificate, verifier=verifier)
        before = verifier.cache_info()
        assert before["evaluations"] > 0

        # A changed terminal formal sum, scheduler metric, or rewrite witness
        # must not be accepted merely because this input has been checked.
        for key, replacement in (("output", []), ("transitions", []),
                                 ("metrics", {**certificate["metrics"], "steps": 0})):
            forged = deepcopy(certificate)
            assert forged[key] != replacement
            forged[key] = replacement
            with pytest.raises(NativeError, match="differs from independent"):
                producer.check(forged, verifier=verifier)
            assert producer.check(saved, verifier=verifier)
        assert verifier.cache_info()["hits"] > before["hits"]
        assert certificate == saved


def test_interrupt_while_capturing_stderr_still_reaps_native_child():
    class InterruptedClose(NativeSession):
        interrupt_stderr_read = False

        @property
        def _stderr(self):
            if self.interrupt_stderr_read:
                self.interrupt_stderr_read = False
                raise KeyboardInterrupt("during close detachment")
            return self.__dict__.get("saved_stderr")

        @_stderr.setter
        def _stderr(self, value):
            self.__dict__["saved_stderr"] = value

    session = InterruptedClose()
    process, stderr = session._process, session._stderr
    assert process is not None and stderr is not None
    session.interrupt_stderr_read = True
    with pytest.raises(KeyboardInterrupt, match="during close detachment"):
        session.close()
    assert process.poll() is not None
    assert session.pid is None
    assert process.stdin.closed and process.stdout.closed and stderr.closed
    session.close()
