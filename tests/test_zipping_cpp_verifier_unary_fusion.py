"""Compare fused DAG execution with its retained unfused schedule in one build."""
from importlib.resources import as_file, files
from pathlib import Path
import shutil
import subprocess

import pytest

from zipper_algebra.zipping_algorithm import native_verifier as verifier

_DRIVER = r'''
#include "scheduler.hpp"
#include "catalog.hpp"
#include <fstream>
#include <iostream>
#include <typeinfo>
using namespace zipper_verify;

J capture(ReplaySession &session, I n, const Word &word, const Limits &limits) {
    try { return session.replay_reused(n,word,limits); }
    catch (const Error &error) {
        return {{"exception",typeid(error).name()},{"message",error.what()}};
    }
}
struct Pair {
    Catalog plain_catalog, fused_catalog;
    ReplaySession plain, fused;
    I calls=0;
    Pair(const J &catalog, std::size_t entries, std::size_t bytes):
        plain_catalog(catalog),fused_catalog(catalog),
        plain(plain_catalog,entries,bytes,false),fused(fused_catalog,entries,bytes,true) {}
    J run(I rank, const Word &word, const Limits &limits) {
        ++calls;
        J left=capture(plain,rank,word,limits),right=capture(fused,rank,word,limits);
        if (!strict_equal(left,right) || !strict_equal(plain.cache_info(),fused.cache_info())) {
            std::cerr << "call " << calls << " rank " << rank << " word " << word
                      << " limits " << limits.json() << "\nUNFUSED " << left
                      << "\nFUSED " << right << '\n';
            throw ReplayError("fused and unfused schedules differ");
        }
        return left;
    }
};
int main(int argc,char **argv) {
    if(argc!=2) return 2;
    std::ifstream stream(argv[1]); J catalog;stream>>catalog;
    I calls=0,completed=0,incomplete=0;
    bool saw_alias=false,saw_hit=false,saw_eviction=false,saw_branch=false;
    const std::vector<std::pair<I,Word>> cases={
        {1,J::array()},{2,J::array({1})},{2,J::array({-1,-1})},
        {3,J::array({1,-2,1})},{3,J::array({-1,-2,-1,-2})},
        {4,J::array({1,2,-3})},{4,J::array({-3,1,2})},
        {5,J::array({1,-4})},{8,J::array({7,-1})},{12,J::array({11})}
    };
    for(const auto &budget:std::vector<std::pair<std::size_t,std::size_t>>{
            {0,0},{1,4096},{8,8192},{64,65536},{4096,8*1024*1024}}) {
        Pair pair(catalog,budget.first,budget.second);
        Limits standard;standard.max_steps=3000;standard.max_terms=1000;
        for(const auto &[rank,word]:cases) {
            J full=pair.run(rank,word,standard);
            if(full.value("status","")=="complete") ++completed;
            if(full.contains("metrics")) {
                const auto &metrics=full.at("metrics");
                I expanded=as_i(metrics.at("expanded_nodes"));
                saw_branch=saw_branch || as_i(metrics.at("peak_live_terms"))>1;
                // Fresh limit keys prohibit a completion cached at another step
                // budget from concealing either side of the refusal boundary.
                for(I steps:std::vector<I>{0,1,2,3,4,5,8,13,
                        std::max<I>(0,expanded-1),expanded,checked_add(expanded,1)}) {
                    Limits limited=standard;limited.max_steps=steps;
                    J result=pair.run(rank,word,limited);
                    if(result.value("status","")=="incomplete") ++incomplete;
                }
            }
            for(I terms:{0,1,2,3,4,8,16}) {
                Limits limited=standard;limited.max_terms=terms;
                pair.run(rank,word,limited);
            }
            for(I size:{0,1,2,3,4,8}) {
                Limits limited=standard;limited.max_word_size=size;
                pair.run(rank,word,limited);
            }
            for(I bound:{0,1,2}) {
                Limits limited=standard;limited.max_coefficient_support=bound;
                pair.run(rank,word,limited);
                limited=standard;limited.max_integer_bits=bound;
                pair.run(rank,word,limited);
            }
            for(I expansion:{0,1,2,4}) {
                Limits limited=standard;limited.expansion_budget=expansion;
                pair.run(rank,word,limited);
            }
            pair.run(rank,word,standard); // repeat root with retained cache
        }
        // Interleaving roots exercises cross-input reuse and shared-payload
        // accounting after eviction, rather than only repeating one root.
        for(auto it=cases.rbegin();it!=cases.rend();++it)
            pair.run(it->first,it->second,standard);
        J cache=pair.plain.cache_info().at("completed");
        saw_hit=saw_hit || as_i(cache.at("hits"))>0;
        saw_eviction=saw_eviction || as_i(cache.at("evictions"))>0;
        saw_alias=saw_alias || as_i(cache.at("entries"))>as_i(cache.at("unique_outputs"));
        pair.plain.clear_cache();pair.fused.clear_cache();
        if(!strict_equal(pair.plain.cache_info(),pair.fused.cache_info())) return 3;
        pair.run(3,J::array({1,-2,1}),standard);
        calls+=pair.calls;
    }
    if(calls<1000 || !completed || !incomplete || !saw_hit || !saw_eviction || !saw_alias || !saw_branch) {
        std::cerr<<"coverage "<<calls<<' '<<completed<<' '<<incomplete<<' '
                 <<saw_hit<<saw_eviction<<saw_alias<<saw_branch<<'\n';return 4;
    }
    std::cout<<calls<<" exact schedule/cache comparisons\n";
}
'''


@pytest.mark.timeout(300)
def test_unary_fusion_preserves_microsteps_limits_and_cache(tmp_path: Path) -> None:
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("C++ compiler is unavailable")
    source = verifier.SOURCE_DIRECTORY
    # Reuse authenticated objects: the harness adds no second compilation of
    # the verifier translation units and links no producer implementation.
    executable, metadata, _ = verifier.build_verifier()
    cache = executable.parent.parent
    objects = []
    for unit in sorted(source.glob("*.cpp")):
        if unit.name == "main.cpp":
            continue
        inputs = {name: digest for name, digest in metadata["source_sha256"].items()
                  if name == unit.name or Path(name).suffix in (".hpp", ".h")}
        identity = {key: metadata[key] for key in
                    ("compiler", "compiler_version", "flags", "libraries")}
        identity.update(source=unit.name, inputs=inputs)
        objects.append(cache / "objects" / verifier._key(identity) / (unit.stem + ".o"))
    driver = tmp_path / "unary_fusion.cpp"
    driver.write_text(_DRIVER)
    binary = tmp_path / "unary_fusion"
    subprocess.run(
        [compiler, "-std=c++17", "-O2", "-I", str(source), str(driver),
         *map(str, objects), "-pthread", "-lgmpxx", "-lgmp", "-lcrypto", "-o", str(binary)],
        check=True, capture_output=True, text=True, timeout=120,
    )
    catalog = files(verifier.__package__).joinpath("data/full-catalog.json")
    with as_file(catalog) as path:
        result = subprocess.run([str(binary), str(path)], capture_output=True,
                                text=True, timeout=120)
    assert result.returncode == 0, result.stderr
