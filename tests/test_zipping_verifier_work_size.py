"""Count-only verifier observation retains the materialized coordinate checks."""
from pathlib import Path
import shutil
import subprocess

import pytest

from zipper_algebra.zipping_algorithm.native_verifier import SOURCE_DIRECTORY


@pytest.mark.timeout(300)
def test_work_size_matches_materialized_expression_and_refusals(tmp_path: Path):
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("native C++ compiler unavailable")
    source = tmp_path / "work_size.cpp"
    source.write_text(r'''
#include "algorithm.hpp"
#include "work_size.hpp"
#include <cassert>
#include <functional>
#include <random>
#include <string>
using namespace zipper_verify;
std::string observe(const std::function<I()> &fn) {
    try { return "ok:"+std::to_string(fn()); }
    catch (const LimitExceeded &) { return "overflow"; }
    catch (const Error &) { return "invalid"; }
    catch (const J::exception &) { return "invalid"; }
}
void compare(const J &work) {
    assert(observe([&] {return length(work_expression(work));})==
           observe([&] {return work_expression_size(work);}));
}
int main() {
    std::mt19937 rng(803196);
    const I lo=std::numeric_limits<I>::min(), hi=std::numeric_limits<I>::max();
    J tokens=J::array({1,-1,0,lo,hi,slider(1,1),slider(hi,lo),slider(lo,hi),
                      nullptr,true,1.5,"bad",J::array({"slider",1,1,false,false,0})});
    for (int trial=0;trial<12000;++trial) {
        J work=initial_work(4,J::array({1,-2,3}));
        work["phase"]=J::array({"single","zip_up","zip_down","extract","other"})[rng()%5];
        work["position"]=tokens[rng()%(trial<6000?3:5)];
        auto word=[&] {
            J out=J::array();
            for (int i=0,n=rng()%5;i<n;++i) out.push_back(tokens[rng()%(trial<6000?3:tokens.size())]);
            return out;
        };
        work["left"]=word();work["processed"]=word();work["word"]=word();
        work["branch"]={{"word",word()}};
        work["boundary_left"]=word();work["boundary_right"]=word();
        work["world_phi"]=bool(rng()%2);
        for (int i=0,n=rng()%4;i<n;++i)
            work["contexts"].push_back({{"size",tokens[rng()%(trial<6000?3:5)]},
                {"top_route",word()},{"component_word",word()},{"bottom_route",word()}});
        compare(work);
    }
    J work=initial_work(3,J::array({1,-2}));
    for (const char *field:{"left","processed","boundary_left","boundary_right"}) {
        J bad=work;bad[field]=J::object();compare(bad);
        bad=work;bad.erase(field);compare(bad);
    }
    for (const char *field:{"rank","phase","position","world_phi","contexts"}) {
        J bad=work;bad.erase(field);compare(bad);
    }
    // Reflection twice cancels algebraically, but the first overflow remains.
    work["phase"]="zip_down";work["world_phi"]=true;
    work["branch"]={{"word",J::array({slider(hi,hi)})}};compare(work);
    assert(observe([&] {return work_expression_size(work);})=="overflow");
    // Inner then outer shifts also preserve intermediate overflow.
    work=initial_work(3,J::array());work["left"]=J::array({hi});
    for (I size:{I(1),I(-1)}) work["contexts"].push_back({{"size",size},
        {"top_route",J::array()},{"component_word",J::array()},{"bottom_route",J::array()}});
    compare(work);
}
''', encoding="utf-8")
    binary = tmp_path / "work_size"
    subprocess.run([compiler, "-std=c++17", "-O1", "-ffunction-sections", "-fdata-sections",
                    "-I", str(SOURCE_DIRECTORY), str(source),
                    *[str(SOURCE_DIRECTORY / name) for name in
                      ("core.cpp", "algorithm.cpp", "algorithm_extract.cpp")],
                    "-Wl,--gc-sections", "-lgmpxx", "-lgmp", "-lcrypto", "-o", str(binary)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True, capture_output=True, text=True)
