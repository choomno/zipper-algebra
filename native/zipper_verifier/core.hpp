#pragma once
// Independent verifier value model. No producer engine headers or objects.
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <gmpxx.h>
#include "vendor/nlohmann/json.hpp"

namespace zipper_verify {
using J = nlohmann::json;
using I = std::int64_t;
using Word = J; // JSON arrays; sliders are ["slider",size,start,M,R,B].
using Big = mpz_class;
using Laurent = std::map<I, Big>;
struct Error : std::runtime_error { using std::runtime_error::runtime_error; };
struct ReplayError : Error { using Error::Error; };
struct LimitExceeded : Error { using Error::Error; };
struct ControlError : Error { using Error::Error; };
struct ControlObligation : Error { using Error::Error; };
struct CatalogError : Error { using Error::Error; };
struct FrameError : Error { using Error::Error; };
struct TransportError : Error { using Error::Error; };
struct TurnaroundError : Error { using Error::Error; };
struct RibbonError : Error { using Error::Error; };
struct RibbonLimit : LimitExceeded { using LimitExceeded::LimitExceeded; };

bool is_int(const J &v);
I as_i(const J &v);
I checked_add(I a, I b);
I checked_sub(I a, I b);
I checked_mul(I a, I b);
I neg(I a);
I abs_i(I a);
I floor_div(I a, I b);
I length(const J &v);
Big integer_value(const J &v); // canonical decimal string or JSON integer, never float/bool
I integer_bits(const Big &value);
Laurent coefficient(const J &value, bool allow_zero = true);
J encoded_coefficient(const Laurent &value);
J one(); // [[0,"1"]]
J add_coefficients(const J &a, const J &b);
J multiply(const J &a, const J &b);
void set_coefficient(J &work_or_branch, const J &value);
J normalize_coefficients(const J &value); // wire numeric integers -> decimal strings

Word slice(const Word &word, I begin, I end = std::numeric_limits<I>::max());
Word concat(std::initializer_list<Word> words);
void append(Word &target, const Word &word);
Word inverse(const Word &word); // pure braid words only
Word bar(const Word &word);
Word reverse_word(const Word &word);
Word phi(const Word &word, I rank);
Word shifted(const Word &word, I shift);
J slider(I size, I start, bool mirrored=false, bool reversed=false, bool barred=false);
bool is_slider(const J &token);
void validate_word(const Word &word, I rank);
void validate_branch(const J &branch, I rank);
J branch(const Word &word, const J &coef=one(), const J &control=J::object());
J run_state(I rank, const J &branches);

struct Limits {
    I max_steps=10000, max_terms=10000, max_word_size=100000;
    I max_coefficient_support=10000, max_integer_bits=100000, expansion_budget=100000;
    Limits()=default;
    explicit Limits(const J &value);
    J json() const;
    void check(const J &state) const;
    void check_branch(const J &value) const;
};
J metrics(const J &state);
std::string canonical_json(const J &value); // Python sorted ensure_ascii JSON; coefficient strings emitted numerically
std::string sha256(const std::string &text);
std::string digest(const J &value);
bool strict_equal(const J &a, const J &b);
void set_reference_sources(const J &sources);
const J &reference_sources();
std::string source_digest(std::initializer_list<std::string> names);
J source_subset(std::initializer_list<std::string> names);

class Catalog;
// Request is the exact Python Request JSON contract. Result is {state,transition}.
J apply_rule(const J &state, const J &request, Catalog &catalog, const Limits &limits, bool trace=true);
J request(I branch_index, I offset, const std::string &rule_id, const J &parameters,
          const J &words=J::object(), const J &evidence=J::object(), const J &controls=nullptr);
} // namespace zipper_verify
