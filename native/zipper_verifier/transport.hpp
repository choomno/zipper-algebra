#pragma once
#include "core.hpp"
#include <tuple>
#include <utility>
namespace zipper_verify {
J transport_frame(const Word &word, I n, I size, I start, I power=1,
                  const std::string &side="left");
void check_frame(const J &witness);
J transport_band(const Word &word, I n, I size, I start, bool mirrored=false);
void check_transport(const J &witness);
using TurnParts = std::tuple<J,Word,Word,Word>;
TurnParts turnaround_parts(const J &branch, I n, I budget=100000);
std::pair<J,J> turnaround(const J &branch, I n, Catalog &catalog,
                         const Limits &limits, bool trace=false,
                         bool want_certificate=true);
}
