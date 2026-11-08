#include "upward.hpp"
#include "catalog.hpp"
#include <algorithm>
#include <map>
#include <set>

namespace zipper_verify {
namespace {
bool pure(const Word &word) {
    return word.is_array() && std::all_of(word.begin(), word.end(), [](const J &x) { return is_int(x); });
}
bool fields(const J &object, std::initializer_list<const char *> names) {
    if (!object.is_object() || object.size() != names.size()) return false;
    return std::all_of(names.begin(), names.end(), [&](const char *name) { return object.contains(name); });
}
I extent(I size, I start, I friends) {
    return checked_sub(checked_add(checked_add(size, start), friends), 1);
}
J delta_for(const J &old, const J &successors) {
    if (successors.empty()) throw ControlError("upward rule has no surviving branches");
    I word_tokens = -length(old.at("word"));
    I support = -length(old.at("coefficient"));
    I max_word = 0, bits = 0;
    I q_min = std::numeric_limits<I>::max(), q_max = std::numeric_limits<I>::min();
    for (const auto &item : successors) {
        I size = length(item.at("word"));
        word_tokens = checked_add(word_tokens, size);
        support = checked_add(support, length(item.at("coefficient")));
        max_word = std::max(max_word, size);
        for (const auto &[q, c] : coefficient(item.at("coefficient"), false)) {
            q_min = std::min(q_min, q); q_max = std::max(q_max, q);
            bits = std::max(bits, integer_bits(c));
        }
    }
    return {{"word_tokens", word_tokens}, {"coefficient_support", support},
            {"max_word_tokens", max_word}, {"q_exponent_min", q_min},
            {"q_exponent_max", q_max}, {"max_integer_bits", bits}};
}
} // namespace

Word upward_friends(I size, I start, I friends) {
    Word result = J::array();
    if (size < 1 || friends < 0) throw ControlError("invalid friend interval");
    for (I a = friends; a > 0; --a) {
        const I first = checked_add(start, a);
        const I stop = checked_add(first, checked_sub(size, 1));
        for (I x = first; x < stop; ++x) result.push_back(x);
    }
    return result;
}

J upward_control(I size, I start, I friends, I left, I processed, I tail) {
    return {{"phase", left ? "zip_up" : "upper_pass_complete"},
            {"i", size}, {"j", start}, {"k", friends}, {"left_length", left},
            {"processed_length", processed}, {"tail_length", tail}};
}

J upward_initialize(I rank, const Word &left, I size, I start, I friends,
                    const Word &processed, const Word &tail) {
    if (size < 2 || size > rank || start < 1 || friends < 0 || extent(size, start, friends) > rank)
        throw ControlError("valid upward slider and friend intervals required");
    if (!pure(left) || !pure(processed) || !pure(tail))
        throw ControlError("exterior words must contain signed braid integers");
    Word word = concat({left, upward_friends(size, start, friends),
                        J::array({slider(size, start)}), processed, tail});
    J control = upward_control(size, start, friends, length(left), length(processed), length(tail));
    return run_state(rank, J::array({branch(word, one(), control)}));
}

std::tuple<J,Word,Word,Word> upward_parts(const J &item, I rank) {
    const J &c = item.at("control");
    if (!fields(c, {"phase", "i", "j", "k", "left_length", "processed_length", "tail_length"}))
        throw ControlError("upward state boundaries are missing or unexpected");
    if ((c.at("phase") != "zip_up" && c.at("phase") != "upper_pass_complete") ||
        !is_int(c.at("i")) || !is_int(c.at("j")) || !is_int(c.at("k")) ||
        !is_int(c.at("left_length")) || !is_int(c.at("processed_length")) || !is_int(c.at("tail_length")))
        throw ControlError("invalid upward state boundaries");
    const I size = as_i(c.at("i")), start = as_i(c.at("j")), friends = as_i(c.at("k"));
    const I a = as_i(c.at("left_length")), b = as_i(c.at("processed_length")), d = as_i(c.at("tail_length"));
    if (size < 2 || size > rank || start < 1 || friends < 0 ||
        extent(size, start, friends) > rank || std::min({a,b,d}) < 0)
        throw ControlError("invalid upward state boundaries");
    const Word core = concat({upward_friends(size, start, friends), J::array({slider(size, start)})});
    const Word &word = item.at("word");
    const I end_core = checked_add(a, length(core));
    if (length(word) != checked_add(checked_add(end_core, b), d) ||
        slice(word, a, end_core) != core ||
        c.at("phase") != (a ? "zip_up" : "upper_pass_complete"))
        throw ControlError("upward boundaries do not describe the full expression");
    Word u = slice(word, 0, a), p = slice(word, end_core, checked_add(end_core, b));
    Word t = slice(word, checked_add(end_core, b));
    if (!pure(u) || !pure(p) || !pure(t))
        throw ControlError("upward exterior words must be pure braid words");
    return {c,u,p,t};
}

std::vector<J> upward_flatten(const J &node) {
    if (node.at("type") != "concat") return {node};
    std::vector<J> result;
    for (const J &child : node.at("items")) {
        auto children = upward_flatten(child);
        result.insert(result.end(), children.begin(), children.end());
    }
    return result;
}

Word upward_literal(const std::vector<J> &nodes, const J &parameters, I budget) {
    J expression = {{"type", "concat"}, {"items", nodes}};
    auto expanded = expand_word(expression, parameters, J::object(), false, budget);
    if (expanded.empty()) throw CatalogError("literal word expansion is empty");
    return expanded.front().word;
}

J upward_successor_control(const J &term, const J &env, I remaining,
                           I processed, I tail, I budget, const std::string &rule_id) {
    auto nodes = upward_flatten(term.at("word"));
    const auto found = std::find_if(nodes.begin(), nodes.end(), [](const J &node) { return node.at("type") == "slider"; });
    if (found == nodes.end()) throw ControlError("catalog branch has no slider");
    const auto at = static_cast<std::size_t>(std::distance(nodes.begin(), found));
    const J &token = nodes[at];
    const I size = integer(token.at("size"), env), start = integer(token.at("start"), env);
    if (token.value("orientation", std::string("normal")) == "mirrored" && size > 2)
        return {{"phase", "unsupported_frame"}, {"rule_id", rule_id},
                {"message", "higher mirrored slider requires transport of retained half twists"}};
    I friends = 0;
    std::size_t before = at;
    if (at && nodes[at-1].at("type") == "block") {
        const J &block = nodes[at-1];
        if (integer(block.at("length"), env) == checked_sub(size, 1) &&
            integer(block.at("start"), env) == checked_add(start, 1) &&
            integer(block.value("sign", J(1)), env) == 1 && block.at("order") == "descending") {
            friends = integer(block.at("count"), env); before = at-1;
        }
    }
    const I emitted_left = length(upward_literal(std::vector<J>(nodes.begin(), nodes.begin()+before), env, budget));
    const I emitted_processed = length(upward_literal(std::vector<J>(nodes.begin()+at+1, nodes.end()), env, budget));
    return upward_control(size, start, friends, checked_add(remaining, emitted_left),
                          checked_add(processed, emitted_processed), tail);
}

std::pair<J,J> upward_step(const J &state, I selected, Catalog &catalog,
                          const Limits &limits, bool trace) {
    try {
        const J &branches = state.at("branches");
        if (selected < 0 || selected >= length(branches))
            throw ControlError("selected upward branch does not exist");
        const J old = branches.at(static_cast<std::size_t>(selected));
        const I rank = as_i(state.at("n"));
        auto [c,u,p,t] = upward_parts(old, rank);
        if (u.empty()) throw ControlError("upward pass is already complete");
        const I size = as_i(c.at("i"));
        J env = {{"n", rank}, {"i", c.at("i")}, {"j", c.at("j")}, {"k", c.at("k")},
                 {"l", abs_i(as_i(u.back()))}, {"eps", as_i(u.back()) > 0 ? 1 : -1}};
        const I prefix_end = checked_add(checked_add(length(u), checked_mul(checked_sub(size, 1), as_i(c.at("k")))), 1);
        std::vector<std::pair<const J *,I>> candidates;
        const std::vector<std::string> excludes = size == 2 ? std::vector<std::string>{"zip.flype."} : std::vector<std::string>{};
        for (const auto &id : catalog.select_rule_ids(env, "zip.", ".I", excludes)) {
            const J &rule=catalog.get(id);
            auto lhs = expand_polynomial(rule.at("lhs"), env, J::object(), false, limits.expansion_budget);
            if (lhs.empty()) throw CatalogError("upward rule has empty left side");
            const Word &occurrence = lhs.front().word;
            const I offset = checked_sub(prefix_end, length(occurrence));
            if (offset >= 0 && offset < length(u) && slice(old.at("word"), offset, prefix_end) == occurrence)
                candidates.emplace_back(&rule, offset);
        }
        if (candidates.size() != 1) {
            if (candidates.empty()) throw ControlError("missing full-tape rule");
            std::string message = "competing upward rules: ";
            for (std::size_t i=0; i<candidates.size(); ++i) {
                if (i) message += ", ";
                message += candidates[i].first->at("id").get<std::string>();
            }
            throw ControlError(message);
        }
        const auto &[selected_rule, offset] = candidates.front();
        const J &rule=*selected_rule;
        std::map<std::string,std::size_t> positions;
        std::vector<J> terms;
        for (const auto &term : rule.at("rhs")) {
            const std::string key = canonical_json(upward_literal({term.at("word")}, env, limits.expansion_budget));
            auto known = positions.find(key);
            if (known == positions.end()) { positions.emplace(key, terms.size()); terms.push_back(term); }
            else terms[known->second] = term;
        }
        J controls = J::array();
        const std::string id = rule.at("id").get<std::string>();
        for (const auto &term : terms)
            controls.push_back(upward_successor_control(term, env, offset, length(p), length(t), limits.expansion_budget, id));
        J evidence = {{"controller", "normal_upward"}, {"full_tape_literal_match", true},
                      {"i2_internal_precedes_flype", size == 2}};
        J applied = apply_rule(state, request(selected, offset, id, env, J::object(), evidence, controls), catalog, limits, trace);
        J following = applied.at("state"), normalizations = J::array();
        for (I idx=selected; idx<checked_add(selected, length(controls)); ++idx) {
            J &item = following["branches"].at(static_cast<std::size_t>(idx));
            bool changed = false;
            for (J &token : item["word"]) {
                if (is_slider(token) && token.at(1) == 2 && token.at(3).get<bool>()) {
                    token = slider(2, as_i(token.at(2))); changed = true;
                }
            }
            if (changed) normalizations.push_back({{"branch", idx}, {"identity", "mirrored_X2_equals_X2"}});
        }
        J successors = slice(following.at("branches"), selected, checked_add(selected, length(controls)));
        for (const auto &successor : successors)
            if (successor.at("control").at("phase") != "unsupported_frame") upward_parts(successor, rank);
        J event = {{"rule_id", id}, {"metric_delta", delta_for(old, successors)}};
        if (trace) {
            event["catalog_transition"] = applied.at("transition");
            event["normalizations"] = normalizations;
            event["successor_state"] = following;
        }
        return {following,event};
    } catch (const CatalogError &error) {
        if (std::string(error.what()).find("budget") != std::string::npos) throw LimitExceeded(error.what());
        throw;
    }
}
} // namespace zipper_verify
