#pragma once
#include "core.hpp"
#include <map>
#include <list>
#include <memory>
#include <unordered_map>
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
    struct Selection { std::vector<std::string> ids; std::size_t bytes; };
    using Selections = std::list<std::pair<std::string,Selection>>;
    mutable Selections selections_;
    mutable std::unordered_map<std::string,Selections::iterator> selection_index_;
    std::size_t selection_max_entries_=0, selection_max_bytes_=0;
    mutable std::size_t selection_bytes_=0;
    mutable I selection_hits_=0, selection_misses_=0, selection_evictions_=0;
    struct Bound { std::shared_ptr<const J> rule; std::size_t bytes; };
    using Bounds = std::list<std::pair<std::string,Bound>>;
    Bounds bounds_;
    std::unordered_map<std::string,Bounds::iterator> bound_index_;
    std::size_t bound_max_entries_=0, bound_max_bytes_=0, bound_bytes_=0;
    I bound_hits_=0, bound_misses_=0, bound_evictions_=0, bound_oversized_=0;
public:
    std::string sha256;
    explicit Catalog(const J &catalog, std::string source_sha256="");
    Catalog(const Catalog &)=delete;
    Catalog &operator=(const Catalog &)=delete;
    Catalog(Catalog &&)=delete;
    Catalog &operator=(Catalog &&)=delete;
    const J &get(const std::string &rule_id) const;
    // The protocol gives each catalog one ReplaySession/cache configuration.
    void configure_selection_cache(std::size_t entries,std::size_t bytes);
    void clear_selection_cache();
    J selection_cache_info() const;
    void configure_bound_cache(std::size_t entries,std::size_t bytes);
    void clear_bound_cache();
    J bound_cache_info() const;
    std::vector<J> rules() const;
    std::vector<std::string> select_rule_ids(const J &environment,
        const std::string &prefix="", const std::string &suffix="",
        const std::vector<std::string> &exclude_prefixes={}) const;
    std::vector<J> select_rules(const J &environment,
        const std::string &prefix="", const std::string &suffix="",
        const std::vector<std::string> &exclude_prefixes={}) const;
    J semilinear_rule(const std::string &source_id, const J &environment);
    std::shared_ptr<const J> semilinear_rule_handle(const std::string &source_id,
                                                 const J &environment);
    J instantiate(const std::string &rule_id, const J &environment,
        const J &words=J::object(), bool expand_sliders=false,
        I max_terms=100000, bool require_verified=true) const;
    J instantiate(const std::string &rule_id, const J &environment, I max_terms) const {
        return instantiate(rule_id, environment, J::object(), false, max_terms, true);
    }
};
} // namespace zipper_verify
