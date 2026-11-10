#pragma once
#include "typed_word.hpp"
#include <optional>
namespace zipper_verify::typed {
using Labels=std::vector<I>;
struct RibbonBlock { I node_size=0,sign=0; Labels chronological_word; };
struct Gathering { I max_image_letters=0,max_operations=0; std::string action_sha256; };
struct RibbonWitness {
 Labels labels,initial_ports,word,crossing_counts,final_ports,reconstructed_word;
 std::vector<RibbonBlock> blocks;
 std::optional<Gathering> gathering;
};
Labels reconstruct_nested(const Labels &labels,const Labels &counts,const Labels &initial);
RibbonWitness recognize_nested_strict(const Labels &word,const Labels &labels,const Labels &initial);
RibbonWitness recognize_nested(const Labels &word,const Labels &labels,const Labels &initial,I max_image_letters=100000,I max_operations=10000000);
J encode_ribbon(const RibbonWitness &);
}
