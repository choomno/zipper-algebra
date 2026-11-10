#pragma once
#include "typed_controller.hpp"
namespace zipper_verify::typed {
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
FrameWitness transport_frame(const Word &, I n, I size, I start, I power = 1,
                             const std::string &side = "left");
void check_frame(const FrameWitness &);
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
BandWitness transport_band(const Word &, I n, I size, I start, bool mirrored = false);
void check_transport(const BandWitness &);
using TurnParts = std::tuple<Control, Word, Word, Word>;
TurnParts turnaround_parts(const Branch &, I n, I budget = 100000);
std::pair<Branch, J> turnaround(const Branch &, I n, Catalog &, const Limits &, bool trace = false,
                                bool want_certificate = true);
} // namespace zipper_verify::typed
