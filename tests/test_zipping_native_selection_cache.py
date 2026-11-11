"""Exact selector-key, bounded retention, and warmed Catalog ownership checks."""
from pathlib import Path
import shutil
import subprocess

import pytest

from zipper_algebra.zipping_algorithm.native import SOURCE_DIRECTORY, NativeSession

pytestmark = pytest.mark.timeout(300)


def test_native_selection_cache_filters_copy_move_and_budgets(tmp_path: Path):
    compiler = shutil.which('g++')
    if compiler is None:
        pytest.skip('native C++ compiler unavailable')
    source = tmp_path / 'selection.cpp'
    source.write_text(r'''
#include "engine.hpp"
#include <cassert>
#include <optional>
#include <string>
#include <new>
#include <cstdlib>
static long remaining = -1;
void *operator new(std::size_t size) {
    // Inject one failure, allowing JSON's exception cleanup to allocate.
    if (remaining == 0) { remaining = -1; throw std::bad_alloc(); }
    if (remaining > 0) --remaining;
    if (auto *result = std::malloc(size ? size : 1)) return result;
    throw std::bad_alloc();
}
void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept { std::free(pointer); }
namespace zipper {
// This isolated fixture has no .M rules requiring source derivation. The stub
// satisfies the unused derived-rule branch while testing the real selector.
std::string sha256_json(const Json &) { return "fixture-unused"; }
}
int main() {
    using namespace zipper;
    auto rule = [](const char *id, Int i) {
        return Json{{"id", id}, {"execution", {{"algebra_verified", true}}},
                    {"lhs", Json::array()}, {"rhs", Json::array()},
                    {"condition", {{"op", "eq"}, {"left", "i"}, {"right", i}}}};
    };
    Json source = {{"rules", Json::array({rule("zip.first.I", 2),
                   rule("zip.flype.test.I", 2), rule("other.I", 3), rule("zip.other.BM", 2)})}};
    auto ids = [](const auto &selected) {
        std::vector<std::string> out;
        for (const auto *item : selected) out.push_back(item->id);
        return out;
    };
    Catalog plain(source), cached(source);
    cached.configure_selection_cache(64, 100000);
    for (Int i : {2, 3, -2}) for (bool exclude : {false, true})
        for (const std::string &prefix : {std::string("zip."), std::string("other")})
            for (const std::string &suffix : {std::string(".I"), std::string(".BM")}) {
                Env env{{"i", i}};
                const auto expected = ids(plain.select(env, prefix, suffix, exclude));
                assert(ids(cached.select(env, prefix, suffix, exclude)) == expected);
                assert(ids(cached.select(env, prefix, suffix, exclude)) == expected);
            }
    assert(cached.selection_cache_info().at("hits") == 24);
    const auto misses = cached.selection_cache_info().at("misses").get<int>();
    cached.select({{"i", 2}, {"unused", -9}}, "zip.", ".I", false);
    assert(cached.selection_cache_info().at("misses") == misses + 1);
    std::optional<Catalog> copy;
    {
        Catalog original(source);
        original.configure_selection_cache(8, 100000);
        original.select({{"i", 2}}, "zip.", ".I", false);
        copy.emplace(original);
    }
    const auto selection = copy->select({{"i", 2}}, "zip.", ".I", false);
    assert(selection.size() == 2);
    for (const auto *item : selection) assert(item == &copy->get(item->id));
    Catalog assigned(source);
    assigned = *copy;
    copy.reset();
    Catalog moved(std::move(assigned));
    assert(moved.selection_cache_info().at("entries") == 0);
    assert(assigned.selection_cache_info().at("bytes") == 0);
    assert(assigned.selection_cache_info().at("entries") == 0);
    assert(assigned.selection_cache_info().at("hits") == 0);
    assert(assigned.select({{"i", 2}}, "zip.", ".I", false).empty());
    for (const auto *item : moved.select({{"i", 2}}, "zip.", ".I", false))
        assert(item == &moved.get(item->id));
    moved.clear_selection_cache();
    assert(moved.selection_cache_info().at("entries") == 0);
    assert(moved.selection_cache_info().at("bytes") == 0);
    moved.configure_selection_cache(1, 100000);
    for (Int i : {2, 3, 2}) moved.select({{"i", i}}, "zip.", ".I", false);
    assert(moved.selection_cache_info().at("entries") == 1);
    assert(moved.selection_cache_info().at("evictions") == 2);
    moved.configure_selection_cache(4, 1);
    assert(moved.select({{"i", 2}}, "zip.", ".I", false).size() == 2);
    assert(moved.selection_cache_info().at("entries") == 0);
    assert(moved.selection_cache_info().at("oversized") == 1);
    moved.configure_selection_cache(0, 10000);
    moved.select({{"i", 2}}, "zip.", ".I", false);
    moved.select({{"i", 2}}, "zip.", ".I", false);
    assert(moved.selection_cache_info().at("entries") == 0);
    assert(moved.selection_cache_info().at("hits") == 0);
    moved.configure_selection_cache(4, 10000);
    bool refused = false;
    try { moved.select({}, "zip.", ".I", false); }
    catch (const CatalogError &) { refused = true; }
    assert(refused && moved.selection_cache_info().at("entries") == 0);
    int failures = 0;
    for (long fail_at : {0, 1, 2, 3, 4, 8, 16, 32, 64, 128, 256, 1024}) {
        Json old_source = {{"rules", Json::array({rule("old.I", 7)})}};
        Catalog target(old_source);
        target.configure_selection_cache(9, 4000);
        target.select({{"i", 7}}, "old", ".I", false);
        const auto old_info = target.selection_cache_info();
        remaining = fail_at;
        bool failed = false;
        try { target = cached; }
        catch (const std::bad_alloc &) { failed = true; }
        remaining = -1;
        if (failed) {
            ++failures;
            assert(target.selection_cache_info() == old_info);
            assert(target.get("old.I").id == "old.I");
            assert(target.select({{"i", 7}}, "old", ".I", false).size() == 1);
        } else {
            assert(target.selection_cache_info().at("entries") == 0);
            assert(target.selection_cache_info().at("bytes") == 0);
            assert(target.select({{"i", 2}}, "zip.", ".I", false).size() == 2);
        }
    }
    assert(failures > 0);
    Catalog move_source(source), move_target(source);
    move_source.configure_selection_cache(8, 4000);
    move_source.select({{"i", 2}}, "zip.", ".I", false);
    move_target = std::move(move_source);
    assert(move_target.selection_cache_info().at("entries") == 0);
    assert(move_source.selection_cache_info().at("entries") == 0);
    assert(move_source.selection_cache_info().at("bytes") == 0);
    assert(move_source.selection_cache_info().at("misses") == 0);
    assert(move_source.select({{"i", 2}}, "zip.", ".I", false).empty());
    assert(move_target.select({{"i", 2}}, "zip.", ".I", false).size() == 2);
}
''', encoding='utf-8')
    binary = tmp_path / 'selection'
    subprocess.run([compiler, '-std=c++17', '-O0', '-ffunction-sections', '-fdata-sections',
                    '-I', str(SOURCE_DIRECTORY), str(source), str(SOURCE_DIRECTORY/'core.cpp'),
                    '-Wl,--gc-sections', '-o', str(binary)], check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True, capture_output=True, text=True)


def test_session_selection_budget_and_clear_are_explicit():
    with NativeSession(cache_entries=17, cache_bytes=4096) as session:
        session.reduce_reused([1, 2, -1], n=3)
        info = session.cache_info()['rule_selection']
        assert info['max_entries'] == 17 and info['max_bytes'] == 4096
        assert info['misses'] > 0 and info['hits'] > 0
        assert info['entries'] <= 17 and info['bytes'] <= 4096
        session.clear_cache()
        info = session.cache_info()['rule_selection']
        assert info['entries'] == info['hits'] == info['misses'] == info['bytes'] == 0


def test_single_alias_admission_recharges_its_payload():
    # Exactly two intrinsic states share one immutable terminal sum. With a
    # one-entry budget the child's final alias must be evicted for its parent.
    for entries in (1, 2):
        with NativeSession(cache_entries=entries, cache_bytes=100000) as session:
            result = session.reduce_reused([], n=1)
            assert result['status'] == 'complete'
            assert result['metrics']['transition_evaluations'] == 2
            info = session.cache_info()['completed']
            assert info['entries'] == entries
            assert info['unique_payloads'] == 1
            assert info['evictions'] == 2 - entries
            assert info['unique_payload_bytes'] > 0
            assert info['bytes'] == info['key_bytes'] + info['unique_payload_bytes']
