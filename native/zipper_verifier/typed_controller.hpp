#pragma once
#include "typed_word.hpp"
#include "catalog.hpp"
#include <optional>
#include <tuple>
namespace zipper_verify::typed {
enum class ControlPhase {
    Up,
    UpComplete,
    FramedUp,
    FramedComplete,
    UnsupportedFrame,
    TurnaroundComplete
};
struct InitialSlider {
    I j = 0, k = 0, output_power = 0, input_power = 0, processed_length = 0;
    bool mirrored = false;
    std::string mirror_policy = "M";
};
struct ProcessedEvent {
    std::string rule_id;
    std::map<std::string, I> parameters;
    I rhs_index = 0;
    std::string mirror_policy = "M";
};
struct Control {
    ControlPhase phase = ControlPhase::UpComplete;
    I i = 0, j = 0, k = 0, left_length = 0, processed_length = 0, tail_length = 0;
    bool mirrored = false;
    I output_power = 0, input_power = 0;
    std::string mirror_policy = "M";
    Word output_payload;
    std::optional<InitialSlider> initial_slider;
    std::vector<ProcessedEvent> processed_events;
    I exterior_length = 0, payload_length = 0, routing_length = 0, transported_length = 0;
    std::string rule_id, message;
};
struct Branch {
    Word word;
    Laurent coefficient;
    Control control;
};
J encode_control(const Control &);
Control decode_control(const J &);
J encode_branch(const Branch &);
Branch decode_branch(const J &);
J encode_state(I rank, const std::vector<Branch> &);
Laurent add_scalar(const Laurent &, const Laurent &);
Laurent multiply_scalar(const Laurent &, const Laurent &);
void check_branch(const Branch &, I rank, const Limits &);
void check_state(I rank, const std::vector<Branch> &, const Limits &);
std::pair<I, I> anchors(I i, I j, I k, bool mirrored);
Word upward_friends(I size, I start, I friends);
Control upward_control(I size, I start, I friends, I left, I processed, I tail);
Branch upward_initialize(I rank, const Word &left, I size = 2, I start = 1, I friends = 0,
                         const Word &processed = {}, const Word &tail = {});
std::tuple<Control, Word, Word, Word> upward_parts(const Branch &, I rank);
Control upward_successor_control(const J &term, const J &parameters, I remaining, I processed,
                                 I tail, I budget, const std::string &rule_id);
Word upward_literal(const std::vector<J> &, const J &parameters, I budget);
Word frame_power(I size, I start, I exponent, I budget);
Word framed_core(I n, I i, I j, I k, bool mirrored, const std::string &policy = "M");
Word friend_word(I n, const Control &);
Control framed_control(I i, I j, I k, bool mirrored, I output_power, I input_power, I left,
                       I processed, I tail);
Word framed_expression(I n, const Control &, const Word &left, const Word &processed,
                       const Word &tail, I budget);
Branch framed_initialize(I n, const Word &left, I size = 3, I start = 1, I friends = 0,
                         bool mirrored = false, I output_power = 0, I input_power = 0,
                         const Word &processed = {}, const Word &tail = {},
                         const Word &output_payload = {}, const std::string &policy = "M",
                         const Limits &limits = Limits());
std::tuple<Control, Word, Word, Word> framed_parts(const Branch &, I n, I budget);
std::pair<std::vector<Branch>, J> upward_step(I rank, const Branch &, Catalog &, const Limits &,
                                              bool trace = true);
std::pair<std::vector<Branch>, J> framed_step(I rank, const Branch &, Catalog &, const Limits &,
                                              bool trace = true);
} // namespace zipper_verify::typed
