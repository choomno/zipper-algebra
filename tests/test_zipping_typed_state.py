"""Differential checks for independently typed extraction and ribbon witnesses."""

from pathlib import Path
import shutil
import subprocess

import pytest


_DRIVER = r"""
#include "typed_extraction.hpp"
#include "extraction.hpp"
#include "algorithm.hpp"
#include <cassert>
#include <iostream>
using namespace zipper_verify;
namespace t=zipper_verify::typed;
int main(){
 int cases=0;
 for(I n=1;n<=10;n++)for(I size=1;size<=n;size++)for(I start=1;start<=n-size+1;start++)for(bool mirror:{false,true}){
  t::Labels selected;for(I x=start;x<start+size;x++)selected.push_back(x);
  t::Word w{t::slider(size,start,mirror)};
  if(size>1){w.insert(w.begin(),I(start));w.push_back(I(-start));}
  auto a=zipper_verify::split_layered(t::encode_word(w),n,J(selected));auto b=t::encode_split(t::split_layered(w,n,selected));if(a!=b){std::cerr<<a<<"\n"<<b<<"\n";return 1;}++cases;
 }
 for(I n=1;n<=9;n++)for(I power=-3;power<=3;power++)for(bool mirror:{false,true}){
  t::Labels labels,ports,counts;for(I x=1;x<=n;x++)ports.push_back(x);labels=ports;std::reverse(labels.begin(),labels.end());for(I x=1;x<n;x++)counts.push_back(power+x%2);
  auto upper=t::reconstruct_nested(labels,counts,ports);auto lower=t::reconstruct_nested(ports,counts,ports);t::Word w;for(I x:upper)w.push_back(x);w.push_back(t::slider(n,1,mirror));for(auto x=lower.rbegin();x!=lower.rend();++x)w.push_back(neg(*x));
  std::string ea,eb;J a,b;try{a=zipper_verify::normalize_canonical_slider(t::encode_word(w));}catch(const Error &e){ea=e.what();}try{b=t::encode_canonical(t::normalize_canonical_slider(w,{},{ }));}catch(const Error &e){eb=e.what();}
  if(a!=b || ea!=eb){std::cerr<<"canonical "<<n<<" "<<power<<" "<<mirror<<"\n"<<ea<<" / "<<eb<<"\n"<<a<<"\n"<<b<<"\n";return 2;}++cases;
 }
 for(I n=1;n<=8;n++)for(I size=1;size<=n;size++)for(I start=1;start<=n-size+1;start++)for(bool mirror:{false,true})for(bool reflected:{false,true}){
  t::Work work;work.rank=n;work.coefficient[0]=1;work.world_phi=reflected;for(I p=1;p<=n;p++){work.bottom_labels.push_back(p);work.top_labels.push_back(p);}
  t::Word w{t::slider(size,start,mirror)};if(size>1){w.insert(w.begin(),I(start));w.push_back(I(-start));}work.phase=t::Extract{w};
  for(bool trace:{false,true}){J a,b;std::string ea,eb;try{a=zipper_verify::extraction_transition(t::encode_work(work),Limits{},trace);}catch(const Error &e){ea=e.what();}try{b=t::encode_transition(t::extraction_transition(work,Limits{},trace));}catch(const Error &e){eb=e.what();}if(a!=b || ea!=eb){std::cerr<<"extraction "<<n<<" "<<size<<" "<<start<<" "<<mirror<<" "<<reflected<<" "<<trace<<"\n"<<ea<<" / "<<eb<<"\n"<<a<<"\n"<<b<<"\n";return 3;}++cases;}
 }
 std::cout<<cases<<" exact helper parity cases\n";
}

"""


def test_typed_extraction_full_witness_and_refusal_parity(tmp_path: Path) -> None:
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("C++ compiler is unavailable")
    source = Path(__file__).resolve().parents[1] / "native" / "zipper_verifier"
    driver = tmp_path / "typed_extraction.cpp"
    driver.write_text(_DRIVER, encoding="utf-8")
    executable = tmp_path / "typed_extraction"
    names = ("core", "typed_word", "typed_state", "typed_ribbons",
             "typed_extraction", "typed_extract", "typed_algorithm",
             "typed_controller", "algorithm", "algorithm_extract",
             "ribbons", "extraction")
    subprocess.run(
        [compiler, "-std=c++17", "-O0", "-ffunction-sections", "-fdata-sections",
         "-I", str(source), str(driver),
         *(str(source / f"{name}.cpp") for name in names),
         "-Wl,--gc-sections", "-lgmpxx", "-lgmp", "-lcrypto", "-o", str(executable)],
        check=True, capture_output=True, text=True, timeout=180,
    )
    subprocess.run([str(executable)], check=True, capture_output=True, text=True, timeout=60)
