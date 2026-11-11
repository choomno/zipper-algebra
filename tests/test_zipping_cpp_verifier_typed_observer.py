"""Typed count-only observation preserves every frozen coordinate-map refusal."""
from pathlib import Path
import hashlib
import io
import json
import shutil
import subprocess
import tarfile

import pytest

from zipper_algebra.zipping_algorithm.native_verifier import SOURCE_DIRECTORY

_DRIVER = r'''
#ifdef TYPED_OBSERVER
#include "typed_state.hpp"
#include "typed_work_size.hpp"
#else
#include "algorithm.hpp"
#include "work_size.hpp"
#endif
#include <functional>
#include <iostream>
#include <random>
#include <typeinfo>
using namespace zipper_verify;
std::string observe(const std::function<I()> &fn) {
    try { return "ok:"+std::to_string(fn()); }
    catch (const Error &e) { return std::string(typeid(e).name())+":"+e.what(); }
    catch (const J::exception &e) { return std::string(typeid(e).name())+":"+e.what(); }
}
void compare(const J &j) {
#ifdef TYPED_OBSERVER
    // Only expression fields are relevant here. Construct those fields directly
    // so adversarial coordinates need not satisfy controller geometry schemas.
    typed::Work w;
    w.rank=as_i(j.at("rank"));w.world_phi=j.at("world_phi");
    w.boundary_left=typed::decode_word(j.at("boundary_left"));
    w.boundary_right=typed::decode_word(j.at("boundary_right"));
    for(const auto &c:j.at("contexts")) {
        typed::Context context;context.size=as_i(c.at("size"));
        context.top_route=typed::decode_word(c.at("top_route"));
        context.component_word=typed::decode_word(c.at("component_word"));
        context.bottom_route=typed::decode_word(c.at("bottom_route"));
        w.contexts.push_back(std::move(context));
    }
    const std::string phase=j.at("phase");
    if(phase=="single")w.phase=typed::Single{typed::decode_word(j.at("left")),typed::decode_word(j.at("processed")),as_i(j.at("position"))};
    else if(phase=="zip_up" || phase=="zip_down") {
        typed::Branch branch;branch.word=typed::decode_word(j.at("branch").at("word"));
        if(phase=="zip_up")w.phase=typed::ZipUp{std::move(branch)};
        else w.phase=typed::ZipDown{std::move(branch)};
    } else if(phase=="layer")w.phase=typed::Layer{typed::decode_word(j.at("word")),0,0};
    else w.phase=typed::Extract{typed::decode_word(j.at("word"))};
    auto materialized=observe([&]{return typed::length(typed::work_expression(w));});
    auto counted=observe([&]{return typed::work_expression_size(w);});
#else
    auto materialized=observe([&]{return length(work_expression(j));});
    auto counted=observe([&]{return work_expression_size(j);});
#endif
    if(materialized!=counted) {
        std::cerr<<"observer differs: "<<j<<"\n"<<materialized<<"\n"<<counted<<'\n';
        std::exit(3);
    }
    std::cout<<materialized<<'\n';
}
J base(I rank) {
    return {{"rank",rank},{"phase","single"},{"left",J::array()},
        {"processed",J::array()},{"position",rank},{"contexts",J::array()},
        {"world_phi",false},{"boundary_left",J::array()},{"boundary_right",J::array()}};
}
int main() {
    std::mt19937 rng(803196);
    const I lo=std::numeric_limits<I>::min(),hi=std::numeric_limits<I>::max();
    J tokens=J::array({1,-1,0,slider(1,1),lo,hi,slider(hi,lo),slider(lo,hi)});
    J coordinates=J::array({-1,0,1,lo,hi});
    for(int trial=0;trial<12000;++trial) {
        J w=base(1+rng()%16);
        w["phase"]=J::array({"single","zip_up","zip_down","layer","extract"})[rng()%5];
        w["position"]=coordinates[rng()%(trial<6000?3:5)];
        auto word=[&] {
            J out=J::array();
            for(int k=0,n=rng()%5;k<n;++k)out.push_back(tokens[rng()%(trial<6000?4:8)]);
            return out;
        };
        w["left"]=word();w["processed"]=word();w["word"]=word();
        w["branch"]={{"word",word()}};w["world_phi"]=bool(rng()%2);
        w["boundary_left"]=word();w["boundary_right"]=word();
        for(int k=0,n=rng()%4;k<n;++k)w["contexts"].push_back({
            {"size",coordinates[rng()%(trial<6000?3:5)]},{"top_route",word()},
            {"component_word",word()},{"bottom_route",word()}});
        compare(w);
    }
    // Two reflections cancel formally; the first signed-coordinate overflow
    // remains observable before the second reflection can undo it.
    J w=base(3);w["phase"]="zip_down";w["world_phi"]=true;
    w["branch"]={{"word",J::array({slider(hi,hi)})}};compare(w);
    // Inner +1 then outer -1 similarly cannot conceal its first overflow.
    w=base(3);w["left"]=J::array({hi});
    for(I size:{I(-1),I(1)})w["contexts"].push_back({{"size",size},
        {"top_route",J::array()},{"component_word",J::array()},{"bottom_route",J::array()}});
    compare(w);
}
'''


def _frozen_sources(destination: Path) -> Path:
    fixture = Path(__file__).parent / "native/frozen_verifier_9c477dd"
    manifest = json.loads(fixture.with_suffix(".json").read_text())
    payload = fixture.with_suffix(".tar.gz").read_bytes()
    assert manifest["commit"] == "9c477dd404738be0a7133bef75b2d4adc37c2616"
    assert hashlib.sha256(payload).hexdigest() == manifest["archive_sha256"]
    names = set()
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
            names.add(member.name)
    assert names == manifest["files"].keys()
    return destination / "native/zipper_verifier"


@pytest.mark.timeout(300)
def test_typed_observer_matches_frozen_materializer_and_exact_refusals(tmp_path: Path):
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("native C++ compiler unavailable")
    frozen = _frozen_sources(tmp_path / "frozen")
    driver = tmp_path / "typed_observer.cpp"
    driver.write_text(_DRIVER)
    outputs = []
    for typed, directory, units in [
        (False, frozen, ("core.cpp", "algorithm.cpp", "algorithm_extract.cpp")),
        (True, SOURCE_DIRECTORY, ("core.cpp", "typed_word.cpp", "typed_algorithm.cpp",
                                 "typed_extract.cpp", "typed_state.cpp")),
    ]:
        binary = tmp_path / ("typed" if typed else "frozen-json")
        subprocess.run(
            [compiler, "-std=c++17", "-O1", "-ffunction-sections", "-fdata-sections",
             *(["-DTYPED_OBSERVER"] if typed else []), "-I", str(directory), str(driver),
             *[str(directory / unit) for unit in units], "-Wl,--gc-sections",
             "-lgmpxx", "-lgmp", "-lcrypto", "-o", str(binary)],
            check=True, capture_output=True, text=True, timeout=180,
        )
        result = subprocess.run([str(binary)], check=True, capture_output=True, timeout=30)
        outputs.append(result.stdout.splitlines())
    assert len(outputs[0]) == len(outputs[1]) == 12002
    assert outputs[0] == outputs[1]
    for result in outputs[0][-2:]:
        assert b"signed integer coordinate overflow" in result
