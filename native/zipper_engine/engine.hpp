#pragma once
#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include "vendor/nlohmann/json.hpp"
namespace zipper {
using Int = std::int64_t;
using BigInt = boost::multiprecision::cpp_int;
using Json = nlohmann::json;
using Braid = std::vector<Int>;
struct Token {
    Int generator = 0, size = 0, start = 0;
    bool mirrored = false, reversed = false, barred = false;
    Token() = default;
    Token(Int x) : generator(x) {}
    static Token slider(Int size, Int start, bool mirrored = false, bool reversed = false,
                        bool barred = false) {
        Token t;
        t.size = size;
        t.start = start;
        t.mirrored = mirrored;
        t.reversed = reversed;
        t.barred = barred;
        return t;
    }
    bool is_slider() const {
        return size != 0;
    }
    auto key() const {
        return std::tie(generator, size, start, mirrored, reversed, barred);
    }
    bool operator==(const Token &x) const {
        return key() == x.key();
    }
    bool operator!=(const Token &x) const {
        return !(*this == x);
    }
    bool operator<(const Token &x) const {
        return key() < x.key();
    }
};
using Word = std::vector<Token>;
using Laurent = std::map<Int, BigInt>;
using Env = std::map<std::string, Int>;
using NamedWords = std::map<std::string, Braid>;
struct EngineError : std::runtime_error {
    using std::runtime_error::runtime_error;
};
struct LimitExceeded : EngineError {
    using EngineError::EngineError;
};
struct ControlError : EngineError {
    using EngineError::EngineError;
};
struct CatalogError : EngineError {
    using EngineError::EngineError;
};
struct FrameError : EngineError {
    using EngineError::EngineError;
};
struct TransportError : EngineError {
    using EngineError::EngineError;
};
struct TurnaroundError : EngineError {
    using EngineError::EngineError;
};
struct Limits {
    Int max_steps = 10000, max_terms = 10000, max_word_size = 100000,
        max_coefficient_support = 10000, max_integer_bits = 100000, expansion_budget = 100000;
};
struct InitialSlider {
    Int j = 1, k = 0, output_power = 0, input_power = 0, processed_length = 0;
    bool mirrored = false;
    std::string mirror_policy = "M";
};
struct ProcessedEvent {
    std::string rule_id;
    Env parameters;
    Int rhs_index = 0;
    std::string mirror_policy = "M";
};
struct Control {
    std::string phase;
    Int i = 2, j = 1, k = 0, left_length = 0, processed_length = 0, tail_length = 0;
    bool mirrored = false;
    Int output_power = 0, input_power = 0;
    Braid output_payload;
    std::string mirror_policy = "M";
    InitialSlider initial_slider;
    std::vector<ProcessedEvent> processed_events;
    Int exterior_length = 0, payload_length = 0, routing_length = 0, transported_length = 0;
    std::string rule_id, message;
};
struct Branch {
    Word word;
    Laurent coefficient{{0, 1}};
    Control control;
};
struct PolynomialTerm {
    Int q = 0;
    Word word;
    BigInt coefficient = 1;
};
using Polynomial = std::vector<PolynomialTerm>; // insertion order of surviving (q,word) keys
struct Rule {
    std::string id, base_id;
    Json lhs, rhs, guard, source;
};
class Catalog {
    std::map<std::string, Rule> rules_;
    std::vector<std::string> ordered_ids_;

  public:
    explicit Catalog(const Json &source);
    const Rule &get(const std::string &id) const;
    std::vector<const Rule *> select(const Env &, const std::string &prefix,
                                     const std::string &suffix, bool exclude_flype = false) const;
};
Int checked_add(Int a, Int b);
Int checked_sub(Int a, Int b);
Int checked_mul(Int a, Int b);
Int integer(const Json &, const Env &);
bool condition(const Json &, const Env &);
Word expand_word(const Json &, const Env &, Int budget, const NamedWords &words = {});
Polynomial expand_polynomial(const Json &, const Env &, Int budget, const NamedWords &words = {});
Polynomial serialize_order(Polynomial);
std::vector<Json> flatten(const Json &);
Word literal(const std::vector<Json> &, const Env &, Int budget);
Laurent multiply(const Laurent &, const Laurent &);
void add_to(Laurent &, const Laurent &);
std::size_t integer_bits(const BigInt &);
void check_branch(const Branch &, const Limits &);
void validate_word(const Word &, Int n);
Word word_from_braid(const Braid &);
Braid braid_from_word(const Word &);
Word slice(const Word &, std::size_t begin, std::size_t end);
Braid slice(const Braid &, std::size_t begin, std::size_t end);
template <class T> void append(std::vector<T> &a, const std::vector<T> &b) {
    a.insert(a.end(), b.begin(), b.end());
}
template <class T> std::vector<T> concat(std::vector<T> a, const std::vector<T> &b) {
    append(a, b);
    return a;
}
std::vector<Branch> apply_rule(const Branch &, const Rule &, const Env &, Int offset,
                               const Limits &, const std::vector<Control> *controls = nullptr,
                               const NamedWords &words = {});
Json word_json(const Word &);
Word word_from_json(const Json &);
Json coefficient_json(const Laurent &);
Laurent coefficient_from_json(const Json &);
Json control_json(const Control &);
Control control_from_json(const Json &);
Json branch_json(const Branch &);
Branch branch_from_json(const Json &);
Json limits_json(const Limits &);
Limits limits_from_json(const Json &);
Json event_json(const ProcessedEvent &);
// Trace-only canonical reference certificates. Numeric coefficient strings are
// rendered as JSON integers when hashing, matching Python arbitrary integers.
void set_reference_sources(const Json &, const std::string &catalog_sha);
std::string reference_source(const std::string &name);
std::string reference_catalog_sha();
std::string canonical_json(const Json &);
std::string sha256_json(const Json &);
Json bound_rule_json(const Rule &, const Env &);
Json contextual_transition(const Branch &before, const std::vector<Branch> &after, const Rule &,
                           const Env &, Int offset, Int n, const Json &evidence,
                           const std::vector<Control> *controls = nullptr,
                           const NamedWords &words = {});
Json run_state_json(Int n, const std::vector<Branch> &);
Json metric_delta(const Branch &, const std::vector<Branch> &);
// Pass controllers.
Braid friends(Int i, Int j, Int k);
std::pair<Int, Int> anchors(Int i, Int j, Int k, bool mirrored);
Word core(Int n, Int i, Int j, Int k, bool mirrored, const std::string &policy = "M");
Braid friend_word(Int n, const Control &);
Braid power(Int size, Int start, Int exponent, Int budget);
Control upward_control(Int i, Int j, Int k, Int left, Int processed, Int tail);
Control framed_control(Int i, Int j, Int k, bool mirrored, Int a, Int b, Int left, Int processed,
                       Int tail);
Word expression(Int n, const Control &, const Braid &u, const Braid &p, const Braid &tail,
                Int budget);
Branch initialize_upward(Int n, const Braid &left, Int size = 2, Int start = 1, Int friends = 0,
                         const Braid &processed = {}, const Braid &tail = {});
Branch initialize_framed(Int n, const Braid &left, Int size, Int start, Int friends, bool mirrored,
                         Int output_power, Int input_power, const Braid &processed,
                         const Braid &tail, const Braid &output_payload, const std::string &policy,
                         const Limits &);
struct Parts {
    Control control;
    Braid u, p, tail;
};
Parts upward_parts(const Branch &, Int n);
Parts framed_parts(const Branch &, Int n, Int budget);
std::tuple<Int, Int, bool, Int, Int> branch_shape(const Json &, const Env &, bool mirrored);
struct BranchStep {
    std::vector<Branch> branches;
    std::string rule_id;
    Json evidence = Json::object();
};
BranchStep upward_step(const Branch &, Int n, const Catalog &, const Limits &, bool trace = false,
                       bool want_evidence = true);
BranchStep framed_step(const Branch &, Int n, const Catalog &, const Limits &, bool trace = false,
                       bool want_evidence = true);
// Geometry transport witnesses, kept typed across transitions.
struct FrameStep {
    std::string kind;
    Int offset = 0, start_before = 0, start_after = 0;
    Braid block, transformed;
};
struct FrameTransport {
    Int n = 0, size = 0, start = 0, end = 0, power = 0;
    std::string side;
    Braid word, transformed;
    std::vector<FrameStep> steps;
};
struct TransportStep {
    std::string kind;
    Int offset = 0, start_before = 0, start_after = 0;
    Braid block;
};
struct BandTransport {
    Int n = 0, size = 0, start = 0, end = 0;
    bool mirrored = false;
    Braid word;
    std::vector<TransportStep> steps;
};
FrameTransport transport_frame(const Braid &, Int n, Int size, Int start, Int power,
                               const std::string &side);
BandTransport transport_band(const Braid &, Int n, Int size, Int start, bool mirrored = false);
void check_frame(const FrameTransport &);
void check_transport(const BandTransport &);
Json frame_json(const FrameTransport &);
Json band_json(const BandTransport &);
Parts turnaround_parts(const Branch &, Int n, Int budget);
std::pair<Branch, Json> turnaround(const Branch &, Int n, const Catalog &, const Limits &,
                                   bool trace = false, bool want_certificate = true);
} // namespace zipper
