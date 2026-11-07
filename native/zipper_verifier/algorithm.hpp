#pragma once
#include "core.hpp"
#include <tuple>
namespace zipper_verify {
J common_work(const J &work);
std::tuple<Word,J,Word> slider_parts(const Word &word);
Word lift_contexts(const Word &word,const J &contexts);
Word work_expression(const J &work);
J initial_work(I rank,const Word &word);
// Unfactored, independently computed transition with scalar/depth restored.
J algorithm_transition(const J &work,Catalog &catalog,const Limits &limits,bool trace);
// Internal phase implementation. Inputs already have unit scalar.
J extraction_transition(const J &work,const Limits &limits,bool trace);
} // namespace zipper_verify
