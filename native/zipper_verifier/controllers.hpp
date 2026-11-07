#pragma once
#include "core.hpp"
#include "catalog.hpp"
#include "upward.hpp"
#include <tuple>
#include <utility>
namespace zipper_verify {
std::pair<I,I> anchors(I i,I j,I k,bool mirrored);
Word frame_power(I size,I start,I exponent,I budget);
Word framed_core(I n,I i,I j,I k,bool mirrored,const std::string &policy="M");
Word friend_word(I n,const J &control);
J rule_for_event(Catalog &catalog,const J &event);
J framed_control(I i,I j,I k,bool mirrored,I output_power,I input_power,I left,I processed,I tail);
Word framed_expression(I n,const J &control,const Word &left,const Word &processed,const Word &tail,I budget);
J framed_initialize(I n,const Word &left,I size=3,I start=1,I friends=0,bool mirrored=false,
                    I output_power=0,I input_power=0,const Word &processed=J::array(),
                    const Word &tail=J::array(),const Word &output_payload=J::array(),
                    const std::string &policy="M",const Limits &limits=Limits());
std::tuple<J,Word,Word,Word> framed_parts(const J &branch,I n,I budget);
std::tuple<I,I,bool,I,I> branch_shape(const J &term,const J &env,bool mirrored);
std::pair<J,J> framed_step(const J &state,I branch_index,Catalog &catalog,const Limits &limits,bool trace=true);
}
