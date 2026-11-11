"""Exact independent guard-cache keys, refusal boundaries, and borrowed records."""
from pathlib import Path
import shutil
import subprocess

import pytest

_DRIVER = r'''
#include "catalog.hpp"
#include <cassert>
#include <type_traits>
using namespace zipper_verify;
int main() {
    auto rule=[](const char *id,J guard) {
        return J{{"id",id},{"lhs",J::array()},{"rhs",J::array()},{"condition",guard}};
    };
    J source={{"rules",J::array({
        rule("zip.a.I",{{"op","ge"},{"left","n"},{"right",3}}),
        rule("zip.b.I",{{"op","eq"},{"left","i"},{"right",2}}),
        rule("zip.c.M",{{"op","le"},{"left","k"},{"right",1}})})}};
    Catalog cached(source), reference(source);
    static_assert(!std::is_copy_constructible_v<Catalog> && !std::is_copy_assignable_v<Catalog>);
    static_assert(!std::is_move_constructible_v<Catalog> && !std::is_move_assignable_v<Catalog>);
    cached.configure_selection_cache(20,100000);
    J env={{"n",4},{"i",2},{"k",0}};
    static_assert(std::is_same_v<decltype(cached.get("zip.a.I")),const J &>);
    const J *borrowed=&cached.get("zip.a.I");
    J copy=*borrowed; copy["id"]="modified";
    assert(cached.get("zip.a.I").at("id")=="zip.a.I");
    auto check=[&](const J &e,const std::string &prefix,const std::string &suffix,
                   const std::vector<std::string> &exclude) {
        auto wanted=reference.select_rule_ids(e,prefix,suffix,exclude);
        assert(cached.select_rule_ids(e,prefix,suffix,exclude)==wanted);
        assert(cached.select_rule_ids(e,prefix,suffix,exclude)==wanted);
    };
    check(env,"zip.",".I",{});
    check(J{{"n",2},{"i",3},{"k",0}},"zip.",".I",{}); // cached empty result
    check(env,"zip.",".M",{});
    check(env,"missing.",".I",{});
    check(env,"zip.",".I",{"zip.a"});
    check(env,"zip.",".I",{"zip.b"});
    check(env,"zip.",".I",{"zip.a","zip.b"});
    auto misses=as_i(cached.selection_cache_info().at("misses"));
    check(env,"zip.",".I",{"zip.b","zip.a"});
    assert(as_i(cached.selection_cache_info().at("misses"))==misses+1);
    J extra=env;extra["extra"]=17;
    check(extra,"zip.",".I",{});
    extra["coefficient"]="unused malformed coefficient extension";
    check(extra,"zip.",".I",{});
    extra["unused"]=std::string(1,static_cast<char>(0xff));
    auto hits=as_i(cached.selection_cache_info().at("hits"));
    check(extra,"zip.",".I",{});
    assert(as_i(cached.selection_cache_info().at("hits"))==hits);
    auto returned=cached.select_rule_ids(env,"zip.",".I",{});
    returned.clear();
    assert(cached.select_rule_ids(env,"zip.",".I",{}).size()==2);
    for (J malformed : {J{{"n",4}},J{{"n",4},{"i",2.0}},J{{"n",true},{"i",2}}}) {
        for (int repeat=0;repeat<2;++repeat) {
            bool failed=false;
            try { cached.select_rule_ids(malformed,"zip.",".I",{}); }
            catch (const CatalogError &) { failed=true; }
            assert(failed);
        }
    }
    cached.configure_selection_cache(1,100000);
    for (I n=2;n<8;++n) {env["n"]=n;check(env,"zip.",".I",{});}
    assert(cached.selection_cache_info().at("entries")==1);
    assert(as_i(cached.selection_cache_info().at("evictions"))>0);
    assert(&cached.get("zip.a.I")==borrowed);
    for (auto budget: {std::pair<size_t,size_t>{0,100000},{10,0},{10,1}}) {
        cached.configure_selection_cache(budget.first,budget.second);
        check(env,"zip.",".I",{});
        assert(cached.selection_cache_info().at("entries")==0);
        assert(cached.selection_cache_info().at("hits")==0);
    }
    cached.clear_selection_cache();
    assert(cached.selection_cache_info().at("bytes")==0);
    assert(cached.selection_cache_info().at("misses")==0);
}
'''


def test_independent_selection_cache_and_immutable_borrows(tmp_path: Path) -> None:
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("C++ compiler is unavailable")
    source = Path(__file__).resolve().parents[1] / "native" / "zipper_verifier"
    driver = tmp_path / "selection.cpp"
    driver.write_text(_DRIVER)
    binary = tmp_path / "selection"
    subprocess.run([compiler, "-std=c++17", "-O0", "-ffunction-sections", "-fdata-sections",
                    "-I", str(source), str(driver), str(source / "core.cpp"), str(source / "catalog.cpp"),
                    "-Wl,--gc-sections", "-lgmpxx", "-lgmp", "-lcrypto", "-o", str(binary)],
                   check=True, capture_output=True, text=True, timeout=120)
    subprocess.run([str(binary)], check=True, capture_output=True, text=True, timeout=10)
