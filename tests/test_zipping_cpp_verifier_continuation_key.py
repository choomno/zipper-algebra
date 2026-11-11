"""Continuation keys preserve the prior exact encoding without copying states."""
from pathlib import Path
import shutil
import subprocess

import pytest

_DRIVER = r'''
#include "core.hpp"
#include <cassert>
#include <random>
#include <typeinfo>
#include <limits>
using namespace zipper_verify;
std::string old_key(J value, bool pending) {
    value.erase("coefficient");
    if (pending) {
        value.erase("depth");
        if (value.contains("branch")) value["branch"].erase("coefficient");
    }
    return canonical_json(value);
}
int main() {
    J left=J::array({1,-2}), right=J::array({3,4});
    assert(concat({left,right})==J::array({1,-2,3,4}));
    assert(concat({concat({J::array({1}),J::array({2})}),right,J::array({-3})})==J::array({1,2,3,4,-3}));
    left=concat({left,left});
    assert(left==J::array({1,-2,1,-2}));
    assert(concat({}).empty());
    bool rejected=false;
    try {concat({J::array(),J(1)});} catch(const ReplayError &) {rejected=true;}
    assert(rejected);
    J coef=J::array({J::array({-3,"1234567890123456789012345678901234567890"})});
    J branch={{"coefficient",coef},{"depth",4},{"word",J::array({-1,2})},
              {"control",{{"coefficient",coef},{"depth",17},{"branch",{{"coefficient",coef}}}}}};
    J value={{"branch",branch},{"coefficient",coef},{"depth",5},{"rank",4},
             {"word",J::array({1,2})},{"unicode",u8"Θ\n\"\\"},
             {u8"é",{{"coefficient",coef},{"depth",19},{"branch",branch}}}};
    std::mt19937 random(617);
    for (int i=0;i<3000;++i) {
        value["depth"]=random(); value["branch"]["depth"]=random();
        value["word"]=J::array({I(random()%10)-5,I(random()%10)-5});
        value["branch"]["word"]=value["word"];
        for(bool pending:{false,true}) {
            J before=value;
            assert(continuation_json(value,pending)==old_key(value,pending));
            assert(value==before);
        }
    }
    J bad_coef=J::array({J::array({0,"01"})});
    J huge_exp=J::array({J::array({std::numeric_limits<std::uint64_t>::max(),"1"})});
    J invalid_text=std::string(1,static_cast<char>(0xff));
    for (J special:{J{{"coefficient",bad_coef}}, J{{"coefficient",huge_exp}},
                    J{{"coefficient",invalid_text}}, J{{"depth",invalid_text}},
                    J{{"coefficient",J::binary({1,2})}}, J{{"depth",J::binary({1,2})}},
                    J{{std::string("coefficient\0x",13),bad_coef}},
                    J{{"nested",{{"coefficient",bad_coef}}}},
                    J{{"nested",{{"coefficient",huge_exp}}}},
                    J{{"nested",invalid_text}},
                    J{{"branch",{{"coefficient",bad_coef}}}},
                    J{{"branch",{{"coefficient",huge_exp}}}},
                    J{{"branch",{{"coefficient",invalid_text}}}},
                    J{{"branch",{{"nested",{{"coefficient",bad_coef}}}}}},
                    J{{"binary",J::binary({1,2})}},
                    J{{std::string(1,static_cast<char>(0xff)),7}},
                    J::object(), J{{"branch",J::object()}},J{{"coefficient",nullptr}},
                    J{{"branch",nullptr}},J{{"branch",J::array()}},J(nullptr),J::array()}) {
        for (bool pending:{false,true}) {
            bool old_failed=false,new_failed=false;
            std::string old_value,new_value,old_error,new_error;
            try {old_value=old_key(special,pending);} catch(const std::exception &e) {old_failed=true;old_error=std::string(typeid(e).name())+e.what();}
            try {new_value=continuation_json(special,pending);} catch(const std::exception &e) {new_failed=true;new_error=std::string(typeid(e).name())+e.what();}
            assert(old_failed==new_failed);
            assert(old_error==new_error);
            if(!old_failed) assert(old_value==new_value);
        }
    }
}
'''


def test_exact_continuation_encoding_without_state_copy(tmp_path: Path) -> None:
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("C++ compiler is unavailable")
    source = Path(__file__).resolve().parents[1] / "native" / "zipper_verifier"
    driver = tmp_path / "key.cpp"
    driver.write_text(_DRIVER)
    binary = tmp_path / "key"
    subprocess.run([compiler, "-std=c++17", "-O0", "-ffunction-sections", "-fdata-sections",
                    "-I", str(source), str(driver), str(source / "core.cpp"),
                    "-Wl,--gc-sections", "-lgmpxx", "-lgmp", "-lcrypto", "-o", str(binary)],
                   check=True, capture_output=True, text=True, timeout=120)
    subprocess.run([str(binary)], check=True, capture_output=True, text=True, timeout=30)
