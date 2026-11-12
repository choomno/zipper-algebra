#pragma once
#include "core.hpp"
#include <map>
#include <string>
#include <vector>

namespace zipper_verify {
struct Monomial {
    I q = 0;
    Word word = J::array();
    Big coefficient = 0;
};
using Polynomial = std::vector<Monomial>; // order of first insertion, with zero terms removed
I integer(const J &expression, const J &environment);
bool condition(const J &expression, const J &environment);
Polynomial expand_word(const J &node, const J &environment,
                       const J &words=J::object(), bool expand_sliders=false,
                       I max_terms=100000);
Polynomial expand_polynomial(const J &terms, const J &environment,
                             const J &words=J::object(), bool expand_sliders=false,
                             I max_terms=100000);
J serialize_polynomial(const Polynomial &poly);
J semilinear_word(const J &word);
void validate_semilinear_bindings(const J &rule, const J &environment);
J semilinear_zip_rule(const J &rule, const J &environment);
J instantiate_rule(const J &rule, const J &environment, const J &words=J::object(),
                   bool expand_sliders=false, I max_terms=100000,
                   bool require_verified=false);

class Catalog {
    J catalog_;
    std::map<std::string, J> rules_;
    std::vector<std::string> order_;
    std::map<std::string, J> semilinear_templates_;
public:
    std::string sha256;
    explicit Catalog(const J &catalog, std::string source_sha256="");
    J get(const std::string &rule_id) const;
    std::vector<J> rules() const;
    std::vector<std::string> select_rule_ids(const J &environment,
        const std::string &prefix="", const std::string &suffix="",
        const std::vector<std::string> &exclude_prefixes={}) const;
    std::vector<J> select_rules(const J &environment,
        const std::string &prefix="", const std::string &suffix="",
        const std::vector<std::string> &exclude_prefixes={}) const;
    J semilinear_rule(const std::string &source_id, const J &environment);
    J instantiate(const std::string &rule_id, const J &environment,
        const J &words=J::object(), bool expand_sliders=false,
        I max_terms=100000, bool require_verified=true) const;
    J instantiate(const std::string &rule_id, const J &environment, I max_terms) const {
        return instantiate(rule_id, environment, J::object(), false, max_terms, true);
    }
};
} // namespace zipper_verify
