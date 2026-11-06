#include "engine.hpp"
#include <openssl/sha.h>
#include <iomanip>
#include <sstream>
namespace zipper {
namespace {
Json source_hashes = Json::object();
std::string catalog_hash;
std::string canonical(const Json &j, const std::string &key = "") {
    if (j.is_object()) {
        std::string s = "{";
        bool first = true;
        for (auto it = j.begin(); it != j.end(); ++it) {
            if (!first)
                s += ",";
            first = false;
            s += Json(it.key()).dump() + ":" + canonical(it.value(), it.key());
        }
        return s + "}";
    }
    if (j.is_array()) {
        std::string s = "[";
        bool first = true;
        for (const auto &item : j) {
            if (!first)
                s += ",";
            first = false;
            if (key == "coefficient" && item.is_array() && item.size() == 2 &&
                item[1].is_string()) {
                s += "[" + item[0].dump() + "," + item[1].get<std::string>() + "]";
            } else
                s += canonical(item);
        }
        return s + "]";
    }
    return j.dump(-1, ' ', true);
}
} // namespace
void set_reference_sources(const Json &sources, const std::string &sha) {
    if (!sources.is_object())
        throw EngineError("reference source map required");
    source_hashes = sources;
    catalog_hash = sha;
}
std::string reference_source(const std::string &name) {
    if (!source_hashes.contains(name))
        throw EngineError("missing reference source identity: " + name);
    return source_hashes.at(name).get<std::string>();
}
std::string reference_catalog_sha() {
    return catalog_hash;
}
std::string canonical_json(const Json &j) {
    return canonical(j);
}
std::string sha256_json(const Json &j) {
    std::string input = canonical_json(j);
    unsigned char bytes[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char *>(input.data()), input.size(), bytes);
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (unsigned char b : bytes)
        out << std::setw(2) << static_cast<unsigned int>(b);
    return out.str();
}
Json bound_rule_json(const Rule &rule, const Env &env) {
    Json result = rule.source;
    if (result.contains("semilinear")) {
        Int i = env.at("i");
        result["semilinear"]["parameters"] = env;
        result["semilinear"]["common_slider_unit"] =
            Json{{"integer", (i - 1) % 2 ? -1 : 1}, {"q_power", -checked_mul(i, i - 1) / 2}};
        result["semilinear"]["normalized_relator_unit"] =
            Json{{"integer", -1}, {"q_power", -(i - 1)}};
    }
    return result;
}
Json run_state_json(Int n, const std::vector<Branch> &branches) {
    Json bs = Json::array();
    for (const auto &b : branches)
        bs.push_back(branch_json(b));
    return Json{{"n", n}, {"branches", bs}};
}
Json contextual_transition(const Branch &before, const std::vector<Branch> &after, const Rule &rule,
                           const Env &env, Int offset, Int n, const Json &evidence,
                           const std::vector<Control> *controls, const NamedWords &words) {
    Json cs = nullptr;
    if (controls) {
        cs = Json::array();
        for (const auto &c : *controls)
            cs.push_back(control_json(c));
    }
    Json request{
        {"branch", 0},    {"offset", offset},     {"rule_id", rule.id},      {"parameters", env},
        {"words", words}, {"evidence", evidence}, {"successor_controls", cs}};
    Json successors = Json::array();
    for (const auto &b : after)
        successors.push_back(branch_json(b));
    return Json{{"request", request},
                {"rule_sha256", sha256_json(bound_rule_json(rule, env))},
                {"predecessor_state_sha256", sha256_json(run_state_json(n, {before}))},
                {"successor_state_sha256", sha256_json(run_state_json(n, after))},
                {"predecessor", branch_json(before)},
                {"successors", successors},
                {"equality_reference",
                 Json{{"catalog_rule", rule.id}, {"kind", "accepted_catalog_identity"}}}};
}
Json metric_delta(const Branch &old, const std::vector<Branch> &after) {
    Int word_sum = 0, support_sum = 0, word_max = 0, bits = 0, qmin = 0, qmax = 0;
    bool first = true;
    for (const auto &b : after) {
        word_sum = checked_add(word_sum, static_cast<Int>(b.word.size()));
        support_sum = checked_add(support_sum, static_cast<Int>(b.coefficient.size()));
        word_max = std::max(word_max, static_cast<Int>(b.word.size()));
        for (const auto &[q, c] : b.coefficient) {
            if (first) {
                qmin = qmax = q;
                first = false;
            } else {
                qmin = std::min(qmin, q);
                qmax = std::max(qmax, q);
            }
            bits = std::max(bits, static_cast<Int>(integer_bits(c)));
        }
    }
    return Json{{"word_tokens", word_sum - static_cast<Int>(old.word.size())},
                {"coefficient_support", support_sum - static_cast<Int>(old.coefficient.size())},
                {"max_word_tokens", word_max},
                {"q_exponent_min", qmin},
                {"q_exponent_max", qmax},
                {"max_integer_bits", bits}};
}
} // namespace zipper
