"""Exercise frame expansion boundaries directly, without extending the wire API."""

from pathlib import Path
import shutil
import subprocess

import pytest


_DRIVER = r"""
#include "controllers.hpp"
#include <cassert>
#include <limits>
using namespace zipper_verify;
int main() {
    const I hi = std::numeric_limits<I>::max();
    const I lo = std::numeric_limits<I>::min();
    assert(frame_power(hi, 1, 0, 0) == J::array());
    assert(frame_power(1, 1, hi, 0) == J::array());
    assert(frame_power(1, 1, -hi, 0) == J::array());
    for (I size : {I(1), hi}) {
        bool caught = false;
        try { frame_power(size, 1, lo, 0); }
        catch (const std::exception&) { caught = true; }
        assert(caught);
    }
    bool caught = false;
    try { frame_power(1, 1, 0, -1); }
    catch (const LimitExceeded&) { caught = true; }
    assert(caught);
    assert(frame_power(3, 2, 1, 3) == J::array({2, 3, 2}));
    assert(frame_power(3, 2, -1, 3) == J::array({-2, -3, -2}));
    caught = false;
    try { frame_power(3, 2, 1, 2); }
    catch (const LimitExceeded&) { caught = true; }
    assert(caught);
}
"""


def test_frame_power_zero_output_and_checked_boundaries(tmp_path: Path) -> None:
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("C++ compiler is unavailable")
    source_dir = Path(__file__).resolve().parents[1] / "native" / "zipper_verifier"
    driver = tmp_path / "frame_boundary.cpp"
    driver.write_text(_DRIVER, encoding="utf-8")
    executable = tmp_path / "frame_boundary"
    # Garbage collection keeps this small driver independent of controller callers.
    subprocess.run(
        [compiler, "-std=c++17", "-O0", "-ffunction-sections", "-fdata-sections",
         "-I", str(source_dir), str(driver), str(source_dir / "core.cpp"),
         str(source_dir / "controllers.cpp"), "-Wl,--gc-sections",
         "-lgmpxx", "-lgmp", "-lcrypto", "-o", str(executable)],
        check=True, capture_output=True, text=True, timeout=120,
    )
    # A huge zero-output frame must finish without constructing or looping over it.
    subprocess.run([str(executable)], check=True, capture_output=True, timeout=5)
