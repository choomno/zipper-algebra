#pragma once
#include "core.hpp"
#include <tuple>
#include <utility>

namespace zipper_verify {
Word upward_friends(I size, I start, I friends);
J upward_control(I size, I start, I friends, I left, I processed, I tail);
J upward_initialize(I rank, const Word &left, I size=2, I start=1, I friends=0,
                    const Word &processed=J::array(), const Word &tail=J::array());
std::tuple<J,Word,Word,Word> upward_parts(const J &branch, I rank);
std::vector<J> upward_flatten(const J &node);
Word upward_literal(const std::vector<J> &nodes, const J &parameters, I budget);
J upward_successor_control(const J &term, const J &parameters, I remaining,
                           I processed, I tail, I budget, const std::string &rule_id);
std::pair<J,J> upward_step(const J &state, I selected, Catalog &catalog,
                          const Limits &limits, bool trace=true);
} // namespace zipper_verify
