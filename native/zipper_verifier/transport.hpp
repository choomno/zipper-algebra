#pragma once
#include "core.hpp"
#include <tuple>
#include <utility>
namespace zipper_verify {
// Owned internal witnesses avoid constructing object trees on untraced replay.
struct FrameStep {
    std::string kind;
    I offset, start_before, start_after;
    Word block, emitted;
};
struct FrameWitness {
    I n, size, start, end, power;
    std::string side;
    Word word, transformed;
    std::vector<FrameStep> steps;
    J json() const;
};
struct BandStep {
    std::string kind;
    I offset, start_before, start_after;
    Word block;
};
struct BandWitness {
    I n, size, start, end;
    bool mirrored;
    Word word;
    std::vector<BandStep> steps;
    J json() const;
};
FrameWitness typed_transport_frame(const Word &word, I n, I size, I start, I power=1,
                                   const std::string &side="left");
void check_frame(const FrameWitness &witness);
BandWitness typed_transport_band(const Word &word, I n, I size, I start, bool mirrored=false);
void check_transport(const BandWitness &witness);
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
