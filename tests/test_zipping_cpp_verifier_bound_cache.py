"""Bound semilinear records keep independent validation and owned lifetimes."""
from pathlib import Path
import shutil
import subprocess

import pytest

_DRIVER = r'''
#include "catalog.hpp"
#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>
using namespace zipper_verify;
int main() {
    J term={{"coefficient",J::array({J{{"integer",1},{"q_power",0}}})},
            {"word",{{"type","slider"},{"size","i"},{"start",1}}}};
    J source={{"id","zip.bound.M"},{"lhs",J::array({term})},{"rhs",J::array({term})},
        {"condition",{{"op","ge"},{"left","i"},{"right",2}}},
        {"execution",{{"algebra_verified",true}}},{"symmetry",J::object()}};
    J second=source;second["id"]="zip.other.M";
    Catalog catalog(J{{"rules",J::array({source,second})}});
    catalog.configure_bound_cache(2,1000000);
    J env={{"n",4},{"i",3}};
    auto first=catalog.semilinear_rule_handle("zip.bound.M",env);
    static_assert(std::is_same_v<decltype(first),std::shared_ptr<const J>>);
    assert(*first==semilinear_zip_rule(source,env));
    assert(catalog.semilinear_rule_handle("zip.bound.M",env)==first);
    assert(catalog.bound_cache_info().at("hits")==1);
    J copy=catalog.semilinear_rule("zip.bound.M",env);copy["id"]="changed";
    assert(first->at("id")=="zip.bound.BM");
    auto other=catalog.semilinear_rule_handle("zip.other.M",env);
    assert(other->at("id")=="zip.other.BM");
    J changed=env;changed["unused"]=17;
    auto changed_record=catalog.semilinear_rule_handle("zip.bound.M",changed);
    assert(changed_record!=first && *changed_record==semilinear_zip_rule(source,changed));
    assert(*first==semilinear_zip_rule(source,env)); // old handle survives nested eviction
    assert(as_i(catalog.bound_cache_info().at("evictions"))>0);
    catalog.configure_bound_cache(20,1000000);
    auto signed_record=catalog.semilinear_rule_handle("zip.bound.M",env);
    J unsigned_env=env;unsigned_env["i"]=J::number_unsigned_t(3);
    auto unsigned_record=catalog.semilinear_rule_handle("zip.bound.M",unsigned_env);
    assert(signed_record!=unsigned_record);
    assert(signed_record->at("semilinear").at("parameters").at("i").is_number_integer());
    assert(unsigned_record->at("semilinear").at("parameters").at("i").is_number_unsigned());
    for (J extension : {J(2.0),J(nullptr),J("unused"),J::array({1}),J::binary({1,2}),
                        J(std::numeric_limits<double>::quiet_NaN())}) {
        J generic=env;generic["extra"]=extension;
        auto a=catalog.semilinear_rule_handle("zip.bound.M",generic);
        auto b=catalog.semilinear_rule_handle("zip.bound.M",generic);
        assert(a!=b); // noninteger extensions use the uncached original path
        assert(a->at("semilinear").at("parameters").at("extra").type()==extension.type());
        if (extension.is_number_float() && std::isnan(extension.get<double>()))
            assert(std::isnan(a->at("semilinear").at("parameters").at("extra").get<double>()));
        else assert(*a==semilinear_zip_rule(source,generic));
    }
    J invalid_utf8=env;invalid_utf8[std::string(1,static_cast<char>(0xff))]=1;
    auto uncached=catalog.semilinear_rule_handle("zip.bound.M",invalid_utf8);
    assert(*uncached==semilinear_zip_rule(source,invalid_utf8));
    assert(catalog.semilinear_rule_handle("zip.bound.M",invalid_utf8)!=uncached);
    auto entries=catalog.bound_cache_info().at("entries");
    for (J invalid : {J::object(),J{{"i",1}},J{{"i",2.0}},J{{"i",true}},
                       J{{"i",std::numeric_limits<I>::max()}}}) {
        for(int repeat=0;repeat<2;++repeat) {
            bool failed=false;
            try {catalog.semilinear_rule_handle("zip.bound.M",invalid);}
            catch(const Error &) {failed=true;}
            assert(failed && catalog.bound_cache_info().at("entries")==entries);
        }
    }
    catalog.clear_bound_cache();
    assert(catalog.bound_cache_info().at("entries")==0 && catalog.bound_cache_info().at("bytes")==0);
    assert(*first==semilinear_zip_rule(source,env)); // clear does not invalidate held records
    catalog.configure_bound_cache(1,1000000);
    auto old=catalog.semilinear_rule_handle("zip.bound.M",env);
    catalog.semilinear_rule_handle("zip.other.M",env);
    assert(*old==semilinear_zip_rule(source,env));
    const auto one_size=as_i(catalog.bound_cache_info().at("bytes"));
    catalog.configure_bound_cache(100,static_cast<size_t>(one_size)+10);
    catalog.semilinear_rule_handle("zip.bound.M",env);
    catalog.semilinear_rule_handle("zip.other.M",env);
    assert(catalog.bound_cache_info().at("entries")==1 && catalog.bound_cache_info().at("evictions")==1);
    for(auto budget : {std::pair<size_t,size_t>{0,1000000},{100,0},{100,1}}) {
        catalog.configure_bound_cache(budget.first,budget.second);
        auto a=catalog.semilinear_rule_handle("zip.bound.M",env);
        assert(*a==semilinear_zip_rule(source,env));
        assert(catalog.semilinear_rule_handle("zip.bound.M",env)!=a);
        assert(catalog.bound_cache_info().at("entries")==0);
    }
}
'''


def test_independent_bound_rule_cache_and_owned_lifetimes(tmp_path: Path) -> None:
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("C++ compiler is unavailable")
    source = Path(__file__).resolve().parents[1] / "native" / "zipper_verifier"
    driver = tmp_path / "bound.cpp"
    driver.write_text(_DRIVER)
    binary = tmp_path / "bound"
    subprocess.run([compiler, "-std=c++17", "-O0", "-ffunction-sections", "-fdata-sections",
                    "-I", str(source), str(driver), str(source / "core.cpp"), str(source / "catalog.cpp"),
                    "-Wl,--gc-sections", "-lgmpxx", "-lgmp", "-lcrypto", "-o", str(binary)],
                   check=True, capture_output=True, text=True, timeout=120)
    subprocess.run([str(binary)], check=True, capture_output=True, text=True, timeout=10)
