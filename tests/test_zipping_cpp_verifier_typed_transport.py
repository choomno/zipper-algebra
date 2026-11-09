"""Typed verifier transport preserves exact legacy witnesses and rejection checks."""
from pathlib import Path
import shutil
import subprocess

import pytest

_DRIVER = r'''
#include "transport.hpp"
#include <cassert>
#include <functional>
#include <typeinfo>
using namespace zipper_verify;
std::pair<std::string,std::string> failure(const std::function<void()> &f) {
    try {f();} catch(const Error &e) {return {typeid(e).name(),e.what()};}
    assert(false);return {};
}
void rejects(const std::function<void()> &f) {failure(f);}
void reject_frame(const FrameWitness &w) {
    assert(failure([&]{check_frame(w);})==failure([&]{check_frame(w.json());}));
}
void reject_band(const BandWitness &w) {
    assert(failure([&]{check_transport(w);})==failure([&]{check_transport(w.json());}));
}
void frame_match(const J &word,I n,I size,I start,I power,const std::string &side) {
    const J legacy=transport_frame(word,n,size,start,power,side);
    auto typed=typed_transport_frame(word,n,size,start,power,side);
    assert(canonical_json(typed.json())==canonical_json(legacy));
    check_frame(legacy);check_frame(typed);
    assert(typed.steps.size()==4);
    auto changed=typed;changed.n=0;reject_frame(changed);
    changed=typed;changed.size=0;reject_frame(changed);
    changed=typed;changed.start=0;reject_frame(changed);
    changed=typed;changed.side="invalid";reject_frame(changed);
    for(size_t k=0;k<typed.steps.size();++k) {
        assert(length(typed.steps[k].block)==size);
        auto changed=typed;changed.steps[k].start_before+=1;reject_frame(changed);
        changed=typed;changed.steps[k].start_after+=1;reject_frame(changed);
        changed=typed;changed.steps[k].block[0]=neg(as_i(changed.steps[k].block[0]));reject_frame(changed);
        changed=typed;changed.steps[k].emitted[0]=neg(as_i(changed.steps[k].emitted[0]));reject_frame(changed);
    }
}
void band_match(const J &word,I n,I size,I start,bool mirror) {
    const J legacy=transport_band(word,n,size,start,mirror);
    auto typed=typed_transport_band(word,n,size,start,mirror);
    assert(canonical_json(typed.json())==canonical_json(legacy));
    check_transport(legacy);check_transport(typed);
    assert(typed.steps.size()==4);
    auto changed=typed;changed.n=0;reject_band(changed);
    changed=typed;changed.size=0;reject_band(changed);
    changed=typed;changed.start=0;reject_band(changed);
    for(size_t k=0;k<typed.steps.size();++k) {
        assert(length(typed.steps[k].block)==size);
        auto changed=typed;changed.steps[k].start_before+=1;reject_band(changed);
        changed=typed;changed.steps[k].start_after+=1;reject_band(changed);
        changed=typed;changed.steps[k].block[0]=neg(as_i(changed.steps[k].block[0]));reject_band(changed);
    }
}
J moving_word(I size,I start,I direction,I sign,bool frame,bool left_side) {
    J word=J::array();I position=start;
    for(I hop=0;hop<4;++hop) {
        const I move=hop%2?-direction:direction;
        const I low=position+(move>0?0:-1);
        const bool descending=frame&&left_side?move<0:move>0;
        J block=J::array();
        for(I j=0;j<size;++j) block.push_back(sign*(low+(descending?size-1-j:j)));
        if(frame&&left_side) word=concat({block,word});else append(word,block);
        position+=move;
    }
    return word;
}
int main() {
    size_t frames=0,bands=0;
    for(I n=1;n<=5;++n) for(I size=1;size<=n;++size) for(I start=1;start<=n-size+1;++start) {
        std::vector<J> words={J::array()};
        for(I a=1-n;a<n;++a) if(a) {
            words.push_back(J::array({a}));
            for(I b=1-n;b<n;++b) if(b) words.push_back(J::array({a,b}));
        }
        for(const auto &word:words) {
            for(const auto &side:{"left","right"}) for(I power:{-3,-2,-1,0,1,2,3}) {
                J legacy;bool valid=true;
                try {legacy=transport_frame(word,n,size,start,power,side);}
                catch(const Error &) {valid=false;}
                if(!valid) {
                    assert(failure([&]{typed_transport_frame(word,n,size,start,power,side);})==
                           failure([&]{transport_frame(word,n,size,start,power,side);}));continue;
                }
                auto typed=typed_transport_frame(word,n,size,start,power,side);
                assert(canonical_json(typed.json())==canonical_json(legacy));check_frame(legacy);check_frame(typed);++frames;
                auto changed=typed;changed.end=0;reject_frame(changed);
                changed=typed;changed.transformed.push_back(1);reject_frame(changed);
                if(!typed.steps.empty()) {
                    changed=typed;changed.steps[0].offset+=1;reject_frame(changed);
                    changed=typed;changed.steps[0].kind="invalid";reject_frame(changed);
                    changed=typed;changed.steps[0].block=J::array();reject_frame(changed);
                    changed=typed;changed.steps[0].emitted=J::array();reject_frame(changed);
                    changed=typed;changed.steps.pop_back();reject_frame(changed);
                }
            }
            for(bool mirrored:{false,true}) {
                J legacy;bool valid=true;
                try {legacy=transport_band(word,n,size,start,mirrored);}
                catch(const Error &) {valid=false;}
                if(!valid) {
                    assert(failure([&]{typed_transport_band(word,n,size,start,mirrored);})==
                           failure([&]{transport_band(word,n,size,start,mirrored);}));continue;
                }
                auto typed=typed_transport_band(word,n,size,start,mirrored);
                assert(canonical_json(typed.json())==canonical_json(legacy));check_transport(legacy);check_transport(typed);++bands;
                auto changed=typed;changed.end=0;reject_band(changed);
                changed=typed;changed.word.push_back(1);reject_band(changed);
                if(!typed.steps.empty()) {
                    changed=typed;changed.steps[0].offset+=1;reject_band(changed);
                    changed=typed;changed.steps[0].kind="invalid";reject_band(changed);
                    changed=typed;changed.steps[0].block=J::array();reject_band(changed);
                    changed=typed;changed.steps[0].start_after+=1;reject_band(changed);
                    changed=typed;changed.steps.pop_back();reject_band(changed);
                }
            }
        }
    }
    assert(frames>5000&&bands>1000);
    // Whole bands of size >=3 require words longer than the exhaustive small corpus.
    // Exercise four successive crossings, both movement directions, signs, mirrors,
    // and frame traversal sides at every available interval through rank seven.
    size_t whole_frames=0,whole_bands=0;
    for(I n=4;n<=7;++n) for(I size=3;size<n;++size)
        for(I start=1;start<=n-size+1;++start) for(I direction:{-1,1}) {
            if(start+direction<1||start+direction>n-size+1) continue;
            for(I sign:{-1,1}) {
                for(bool left_side:{false,true}) for(I power:{-3,-2,-1,0,1,2,3}) {
                    J word=moving_word(size,start,direction,sign,true,left_side);
                    frame_match(word,n,size,start,power,left_side?"left":"right");++whole_frames;
                }
                for(bool mirror:{false,true}) {
                    J word=moving_word(size,start,direction,sign,false,false);
                    band_match(word,n,size,start,mirror);++whole_bands;
                }
            }
        }
    assert(whole_frames==1120&&whole_bands==160);
    auto parity=typed_transport_frame(J::array({1}),4,3,1,0);
    parity.power=1;reject_frame(parity);
    for(J bad:{J::array({true}),J::array({1.0}),J::array({0}),J::array({5}),J::object()}) {
        assert(failure([&]{typed_transport_frame(bad,5,2,1);})==failure([&]{transport_frame(bad,5,2,1);}));
        assert(failure([&]{typed_transport_band(bad,5,2,1);})==failure([&]{transport_band(bad,5,2,1);}));
    }
    rejects([]{typed_transport_frame(J::array(),3,0,1);});
    rejects([]{typed_transport_frame(J::array(),3,2,1,1,"invalid");});
    rejects([]{typed_transport_band(J::array(),3,2,3);});
}
'''


def test_typed_transport_witnesses_and_adversarial_checks(tmp_path: Path) -> None:
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("C++ compiler is unavailable")
    source = Path(__file__).resolve().parents[1] / "native" / "zipper_verifier"
    driver = tmp_path / "typed_transport.cpp"
    driver.write_text(_DRIVER)
    binary = tmp_path / "typed_transport"
    subprocess.run(
        [compiler, "-std=c++17", "-O1", "-ffunction-sections", "-fdata-sections",
         "-I", str(source), str(driver), str(source / "core.cpp"), str(source / "bands.cpp"),
         "-Wl,--gc-sections", "-lgmpxx", "-lgmp", "-lcrypto", "-o", str(binary)],
        check=True, capture_output=True, text=True, timeout=120,
    )
    subprocess.run([str(binary)], check=True, capture_output=True, text=True, timeout=30)
