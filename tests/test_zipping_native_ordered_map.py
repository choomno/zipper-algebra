"""Allocation failure cannot leave an unindexed native scheduler/cache entry."""
from pathlib import Path
import shutil
import subprocess

import pytest

from zipper_algebra.zipping_algorithm.native import SOURCE_DIRECTORY


def test_ordered_map_insert_rolls_back_after_allocation_failure(tmp_path: Path):
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("native C++ compiler unavailable")
    source = tmp_path / "failure.cpp"
    source.write_text(r'''
#include <cassert>
#include <cstdlib>
#include <new>
#include <string>
#include "ordered.hpp"

static long remaining = -1;
void *operator new(std::size_t size) {
    if (remaining == 0) throw std::bad_alloc();
    if (remaining > 0) --remaining;
    if (auto *result = std::malloc(size ? size : 1)) return result;
    throw std::bad_alloc();
}
void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept { std::free(pointer); }
int main() {
    int failures = 0, successes = 0;
    for (int populated : {1, 13}) {
        for (int fail_at = 0; fail_at < 12; ++fail_at) {
            zipper::OrderedMap<int> map;
            for (int i = 0; i < populated; ++i)
                map.set("existing" + std::to_string(i), i);
            std::string candidate(80, 'x');
            const auto expected_size = map.size();
            remaining = fail_at;
            bool threw = false;
            try { map.set(std::move(candidate), 99); }
            catch (const std::bad_alloc &) { threw = true; }
            remaining = -1;
            if (threw) {
                ++failures;
                assert(map.size() == expected_size);
                assert(!map.contains(std::string(80, 'x')));
            } else {
                ++successes;
                assert(map.size() == expected_size + 1);
                assert(*map.get(std::string(80, 'x')) == 99);
            }
            std::size_t count = 0;
            for (const auto &[key, value] : map) {
                assert(map.contains(key));
                assert(*map.get(key) == value);
                ++count;
            }
            assert(count == map.size());
            map.set("after", 100);
            assert(*map.get("after") == 100);
            map.erase("after");
            assert(!map.contains("after"));
            int expected = 0;
            for (const auto &[key, value] : map) {
                if (key.rfind("existing", 0) == 0) assert(value == expected++);
            }
            assert(expected == populated);
        }
    }
    assert(failures > 0 && successes > 0);
    zipper::OrderedMap<int> touched;
    touched.set("a", 1); touched.set("b", 2); touched.set("c", 3);
    const std::string b = "b", missing = "missing";
    auto *before = touched.get(b);
    remaining = 0;
    assert(touched.touch(b) == before);
    assert(touched.touch(missing) == nullptr);
    remaining = -1;
    auto it = touched.begin();
    assert(it++->first == "a");
    assert(it++->first == "c");
    assert(it++->first == "b");
    assert(it == touched.end());
}
''', encoding="utf-8")
    binary = tmp_path / "failure"
    subprocess.run([compiler, "-std=c++17", "-O0", "-Wall", "-Wextra",
                    "-I", str(SOURCE_DIRECTORY), str(source), "-o", str(binary)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True, capture_output=True, text=True)
