#include "engine.hpp"

#include <algorithm>
#include <limits>
#include <tuple>

namespace zipper {
namespace {
// All persistent braid words are compact vectors; these helpers concatenate
// native storage without a JSON conversion or a Python call.
void append_braid(Word &dst, const Braid &src) {
    for (auto x : src)
        dst.push_back(Token(x));
}
Braid braid_slice(const Word &word, std::size_t a, std::size_t b) {
    if (a > b || b > word.size())
        throw ControlError("invalid word boundaries");
    Braid out;
    out.reserve(b - a);
    for (std::size_t k = a; k < b; ++k) {
        if (word[k].is_slider())
            throw ControlError("exterior words must contain signed braid integers");
        out.push_back(word[k].generator);
    }
    return out;
}
bool starts(const Word &word, const Word &prefix) {
    return word.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), word.begin());
}
void validate_braid(const Braid &word, std::int64_t n) {
    for (auto x : word)
        if (x == 0 || x == std::numeric_limits<std::int64_t>::min() || std::abs(x) >= n)
            throw ControlError("word must contain signed braid generators in B_n");
}
void flatten_nodes(const Json &node, std::vector<Json> &out) {
    if (node.at("type") == "concat")
        for (const auto &child : node.at("items"))
            flatten_nodes(child, out);
    else
        out.push_back(node);
}
Word literal_nodes(const std::vector<Json> &nodes, std::size_t a, std::size_t b, const Env &env,
                   std::int64_t budget) {
    Json items = Json::array();
    for (std::size_t k = a; k < b; ++k)
        items.push_back(nodes[k]);
    return expand_word(Json{{"type", "concat"}, {"items", items}}, env, budget);
}
Control upward_control_impl(std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t left,
                            std::int64_t processed, std::int64_t tail) {
    Control c;
    c.phase = left ? "zip_up" : "upper_pass_complete";
    c.i = i;
    c.j = j;
    c.k = k;
    c.left_length = left;
    c.processed_length = processed;
    c.tail_length = tail;
    return c;
}
Control framed_control_impl(std::int64_t i, std::int64_t j, std::int64_t k, bool mirrored,
                            std::int64_t a, std::int64_t b, std::int64_t left,
                            std::int64_t processed, std::int64_t tail) {
    Control c = upward_control_impl(i, j, k, left, processed, tail);
    c.phase = left ? "framed_zip_up" : "framed_upper_pass_complete";
    c.mirrored = mirrored;
    c.output_power = a;
    c.input_power = b;
    return c;
}
Control successor_control(const Json &term, const Env &env, std::int64_t remaining,
                          std::int64_t processed, std::int64_t tail, std::int64_t budget,
                          const std::string &rule_id) {
    std::vector<Json> nodes;
    flatten_nodes(term.at("word"), nodes);
    auto at = std::find_if(nodes.begin(), nodes.end(),
                           [](const Json &v) { return v.at("type") == "slider"; });
    if (at == nodes.end())
        throw ControlError("catalog branch has no slider");
    const auto index = std::size_t(at - nodes.begin());
    auto i = integer(at->at("size"), env), j = integer(at->at("start"), env);
    if (at->value("orientation", std::string("normal")) == "mirrored" && i > 2) {
        Control c;
        c.phase = "unsupported_frame";
        c.rule_id = rule_id;
        c.message = "higher mirrored slider requires transport of retained half twists";
        return c;
    }
    std::int64_t k = 0;
    auto before = index;
    if (index && nodes[index - 1].at("type") == "block") {
        const auto &block = nodes[index - 1];
        if (integer(block.at("length"), env) == i - 1 && integer(block.at("start"), env) == j + 1 &&
            integer(block.value("sign", Json(1)), env) == 1 && block.at("order") == "descending") {
            k = integer(block.at("count"), env);
            before = index - 1;
        }
    }
    auto l = literal_nodes(nodes, 0, before, env, budget).size();
    auto p = literal_nodes(nodes, index + 1, nodes.size(), env, budget).size();
    return upward_control_impl(i, j, k, remaining + l, processed + p, tail);
}
} // namespace

Control upward_control(Int i, Int j, Int k, Int left, Int processed, Int tail) {
    return upward_control_impl(i, j, k, left, processed, tail);
}
Control framed_control(Int i, Int j, Int k, bool mirrored, Int a, Int b, Int left, Int processed,
                       Int tail) {
    return framed_control_impl(i, j, k, mirrored, a, b, left, processed, tail);
}

Braid friends(std::int64_t i, std::int64_t j, std::int64_t k) {
    Braid result;
    if (i < 1 || j < 1 || k < 0)
        throw ControlError("invalid friend interval");
    for (auto a = k; a > 0; --a) {
        const auto begin = checked_add(j, a);
        const auto end = checked_add(begin, i - 1);
        for (auto x = begin; x < end; ++x)
            result.push_back(x);
    }
    return result;
}
std::pair<std::int64_t, std::int64_t> anchors(std::int64_t i, std::int64_t j, std::int64_t k,
                                              bool mirrored) {
    (void)i;
    return mirrored ? std::make_pair(checked_sub(j, k), checked_add(j, 1))
                    : std::make_pair(checked_add(checked_add(j, k), 1), j);
}
Word core(std::int64_t n, std::int64_t i, std::int64_t j, std::int64_t k, bool mirrored,
          const std::string &policy) {
    auto f = friends(i, mirrored ? n - i - j + 2 : j, k);
    Word out;
    out.reserve(f.size() + 1);
    for (auto x : f)
        out.push_back(Token(mirrored ? (policy == "bar_mirror" ? -(n - x) : n - x) : x));
    out.push_back(Token::slider(i, j, mirrored));
    return out;
}
Braid friend_word(std::int64_t n, const Control &c) {
    auto word = core(n, c.i, c.j, c.k, c.mirrored, c.mirror_policy);
    return braid_slice(word, 0, word.size() - 1);
}
Braid power(std::int64_t size, std::int64_t start, std::int64_t exponent, std::int64_t budget) {
    if (size < 1 || start < 1 || budget < 0)
        throw LimitExceeded("expanded frame word limit exceeded");
    auto count = __int128(size) * (__int128(size) - 1) / 2;
    auto magnitude = exponent < 0 ? -__int128(exponent) : __int128(exponent);
    if (count == 0 || magnitude == 0)
        return {};
    if (count > budget || magnitude > __int128(budget) / count)
        throw LimitExceeded("expanded frame word limit exceeded");
    if (__int128(start) + size - 2 > std::numeric_limits<Int>::max())
        throw LimitExceeded("expanded frame coordinate limit exceeded");
    Braid half;
    for (auto length = size - 1; length > 0; --length)
        for (Int offset = 0; offset < length; ++offset)
            half.push_back(checked_add(start, offset));
    if (exponent < 0) {
        std::reverse(half.begin(), half.end());
        for (auto &x : half)
            x = -x;
    }
    Braid out;
    out.reserve(std::size_t(count * magnitude));
    for (__int128 k = 0; k < magnitude; ++k)
        append(out, half);
    return out;
}
Word expression(std::int64_t n, const Control &c, const Braid &u, const Braid &p, const Braid &tail,
                std::int64_t budget) {
    const auto i = c.i;
    auto magnitude = [](std::int64_t x) { return x < 0 ? -__int128(x) : __int128(x); };
    if (i < 2 || c.k < 0 || budget < 0)
        throw ControlError("invalid framed state boundaries");
    auto total = __int128(u.size()) + c.output_payload.size() + p.size() + tail.size() +
                 (__int128(i) - 1) * c.k + 1;
    auto frame_count = (__int128(i) - 1) * (i - 2) / 2;
    auto frame_power = magnitude(c.output_power) + magnitude(c.input_power);
    if (total > budget || (frame_count && frame_power > (__int128(budget) - total) / frame_count))
        throw LimitExceeded("expanded framed state word limit exceeded");
    total += frame_count * frame_power;
    auto [out, inp] = anchors(i, c.j, c.k, c.mirrored);
    Word result;
    result.reserve(std::size_t(total));
    append_braid(result, u);
    append_braid(result, c.output_payload);
    append_braid(result, power(i - 1, out, c.output_power, budget));
    append(result, core(n, i, c.j, c.k, c.mirrored, c.mirror_policy));
    append_braid(result, power(i - 1, inp, c.input_power, budget));
    append_braid(result, p);
    append_braid(result, tail);
    return result;
}
Branch initialize_upward(std::int64_t n, const Braid &left, std::int64_t size, std::int64_t start,
                         std::int64_t friend_count, const Braid &processed, const Braid &tail) {
    if (size < 2 || size > n || start < 1 || friend_count < 0 ||
        __int128(size) + start + friend_count - 1 > n)
        throw ControlError("valid upward slider and friend intervals required");
    validate_braid(left, n);
    validate_braid(processed, n);
    validate_braid(tail, n);
    Branch b;
    b.control =
        upward_control_impl(size, start, friend_count, left.size(), processed.size(), tail.size());
    append_braid(b.word, left);
    append_braid(b.word, friends(size, start, friend_count));
    b.word.push_back(Token::slider(size, start));
    append_braid(b.word, processed);
    append_braid(b.word, tail);
    b.coefficient = {{0, 1}};
    return b;
}
Branch initialize_framed(std::int64_t n, const Braid &left, std::int64_t size, std::int64_t start,
                         std::int64_t friend_count, bool mirrored, std::int64_t a, std::int64_t b,
                         const Braid &processed, const Braid &tail, const Braid &payload,
                         const std::string &policy, const Limits &limits) {
    if (size < 2 || size > n || start < 1 || friend_count < 0 || __int128(start) + size - 1 > n ||
        (mirrored ? start - friend_count < 1 : __int128(start) + size + friend_count - 1 > n) ||
        (policy != "M" && policy != "bar_mirror"))
        throw ControlError("valid framed slider, friends and integer powers required");
    if (size == 2) {
        if (mirrored || !payload.empty())
            throw ControlError("use the normal X2 convention for the two-strand delegate");
        return initialize_upward(n, left, size, start, friend_count, processed, tail);
    }
    validate_braid(left, n);
    validate_braid(processed, n);
    validate_braid(tail, n);
    Control c = framed_control_impl(size, start, friend_count, mirrored, a, b, left.size(),
                                    processed.size(), tail.size());
    c.mirror_policy = policy;
    c.output_payload = payload;
    auto [out, inp] = anchors(size, start, friend_count, mirrored);
    (void)inp;
    for (auto x : payload)
        if (x == std::numeric_limits<std::int64_t>::min() || std::abs(x) < out ||
            std::abs(x) >= checked_add(out, size - 2))
            throw ControlError("protected output payload must be supported on the marked tape");
    c.initial_slider = {start,    friend_count, a, b, std::int64_t(processed.size()),
                        mirrored, policy};
    Branch result{expression(n, c, left, processed, tail, limits.max_word_size), {{0, 1}}, c};
    check_branch(result, limits);
    return result;
}
Parts upward_parts(const Branch &branch, std::int64_t n) {
    const auto &c = branch.control;
    if ((c.phase != "zip_up" && c.phase != "upper_pass_complete") || c.i < 2 || c.i > n ||
        c.j < 1 || c.k < 0 || __int128(c.i) + c.j + c.k - 1 > n ||
        std::min({c.left_length, c.processed_length, c.tail_length}) < 0)
        throw ControlError("invalid upward state boundaries");
    auto a = std::size_t(c.left_length), b = std::size_t(c.processed_length),
         d = std::size_t(c.tail_length);
    const auto core_size = (__int128(c.i) - 1) * c.k + 1;
    if (__int128(a) + core_size + b + d != branch.word.size())
        throw ControlError("upward boundaries do not describe the full expression");
    auto middle = core(n, c.i, c.j, c.k, false, "M");
    if (__int128(a) + middle.size() + b + d != branch.word.size() ||
        !std::equal(middle.begin(), middle.end(), branch.word.begin() + a) ||
        c.phase != (a ? "zip_up" : "upper_pass_complete"))
        throw ControlError("upward boundaries do not describe the full expression");
    return {c, braid_slice(branch.word, 0, a),
            braid_slice(branch.word, a + middle.size(), a + middle.size() + b),
            braid_slice(branch.word, a + middle.size() + b, branch.word.size())};
}
Parts framed_parts(const Branch &branch, std::int64_t n, std::int64_t budget) {
    const auto &c = branch.control;
    if (c.i < 3 || c.i > n || c.j < 1 || c.k < 0 || __int128(c.j) + c.i - 1 > n ||
        (c.mirrored ? c.j - c.k < 1 : __int128(c.j) + c.i + c.k - 1 > n) ||
        std::min({c.left_length, c.processed_length, c.tail_length}) < 0 ||
        (c.mirror_policy != "M" && c.mirror_policy != "bar_mirror"))
        throw ControlError("invalid framed state boundaries");
    auto [out, inp] = anchors(c.i, c.j, c.k, c.mirrored);
    (void)inp;
    for (auto x : c.output_payload)
        if (x == std::numeric_limits<std::int64_t>::min() || std::abs(x) < out ||
            std::abs(x) >= checked_add(out, c.i - 2))
            throw ControlError("invalid protected output-tape payload");
    if (c.initial_slider.mirror_policy != c.mirror_policy)
        throw ControlError("invalid framed input-band record");
    for (const auto &event : c.processed_events)
        if (event.mirror_policy != c.mirror_policy || event.rhs_index < 0)
            throw ControlError("invalid processed-band event");
    if (c.phase != (c.left_length ? "framed_zip_up" : "framed_upper_pass_complete"))
        throw ControlError("invalid framed phase");
    if (__int128(c.processed_length) + c.tail_length > branch.word.size())
        throw ControlError("invalid framed state boundaries");
    auto p_start = branch.word.size() - c.processed_length - c.tail_length;
    Parts parts{c, braid_slice(branch.word, 0, c.left_length),
                braid_slice(branch.word, p_start, p_start + c.processed_length),
                braid_slice(branch.word, p_start + c.processed_length, branch.word.size())};
    if (branch.word != expression(n, c, parts.u, parts.p, parts.tail, budget))
        throw ControlError("framed boundaries do not describe the complete expression");
    return parts;
}
std::tuple<Int, Int, bool, Int, Int> branch_shape(const Json &term, const Env &env, bool mirrored) {
    std::vector<Json> nodes;
    flatten_nodes(term.at("word"), nodes);
    auto at = std::find_if(nodes.begin(), nodes.end(),
                           [](const Json &v) { return v.at("type") == "slider"; });
    if (at == nodes.end())
        throw ControlError("catalog branch has no slider");
    auto idx = std::size_t(at - nodes.begin());
    auto j = integer(at->at("start"), env);
    bool base_mirror = at->value("orientation", std::string("normal")) == "mirrored";
    std::int64_t k = 0;
    if (idx && nodes[idx - 1].at("type") == "block")
        k = integer(nodes[idx - 1].at("count"), env);
    auto exponent = [](const Json &node) -> std::int64_t {
        if (node.at("type") == "half_twist")
            return 1;
        if (node.at("type") == "inverse" && node.at("word").at("type") == "half_twist")
            return -1;
        return 0;
    };
    return {mirrored ? env.at("n") - env.at("i") - j + 2 : j, k, base_mirror != mirrored,
            idx ? exponent(nodes.front()) : 0,
            idx + 1 < nodes.size() ? exponent(nodes[idx + 1]) : 0};
}

BranchStep upward_step(const Branch &predecessor, Int n, const Catalog &catalog,
                       const Limits &limits, bool trace) {
    auto parts = upward_parts(predecessor, n);
    const auto &c = parts.control;
    const auto &u = parts.u;
    const auto &p = parts.p;
    const auto &t = parts.tail;
    if (u.empty())
        throw ControlError("upward pass is already complete");
    Env env{{"n", n},
            {"i", c.i},
            {"j", c.j},
            {"k", c.k},
            {"l", std::abs(u.back())},
            {"eps", u.back() > 0 ? 1 : -1}};
    auto prefix_end = checked_add(checked_add(Int(u.size()), checked_mul(c.i - 1, c.k)), 1);
    std::vector<std::pair<const Rule *, Int>> candidates;
    for (const auto *rule : catalog.select(env, "zip.", ".I", c.i == 2)) {
        auto lhs = expand_polynomial(rule->lhs, env, limits.expansion_budget);
        if (lhs.size() != 1)
            throw ControlError("contextual replay requires a monomial left side");
        const auto &occurrence = lhs.front().word;
        auto offset = prefix_end - Int(occurrence.size());
        if (offset >= 0 && offset < Int(u.size()) && prefix_end <= Int(predecessor.word.size()) &&
            std::equal(occurrence.begin(), occurrence.end(), predecessor.word.begin() + offset))
            candidates.emplace_back(rule, offset);
    }
    if (candidates.size() != 1)
        throw ControlError(candidates.empty() ? "missing full-tape rule"
                                              : "competing upward rules");
    const auto &rule = *candidates.front().first;
    auto offset = candidates.front().second;
    std::vector<std::pair<Word, Json>> terms;
    for (const auto &term : rule.rhs) {
        auto word = expand_word(term.at("word"), env, limits.expansion_budget);
        auto found = std::find_if(terms.begin(), terms.end(),
                                  [&](const auto &v) { return v.first == word; });
        if (found == terms.end())
            terms.emplace_back(std::move(word), term);
        else
            found->second = term;
    }
    std::vector<Control> controls;
    for (const auto &item : terms)
        controls.push_back(successor_control(item.second, env, offset, p.size(), t.size(),
                                             limits.expansion_budget, rule.id));
    auto following = apply_rule(predecessor, rule, env, offset, limits, &controls);
    Json equality;
    if (trace)
        equality = contextual_transition(predecessor, following, rule, env, offset, n,
                                         Json{{"controller", "normal_upward"},
                                              {"full_tape_literal_match", true},
                                              {"i2_internal_precedes_flype", c.i == 2}},
                                         &controls);
    Json normalizations = Json::array();
    for (std::size_t index = 0; index < following.size(); ++index) {
        auto &branch = following[index];
        bool changed = false;
        for (auto &token : branch.word)
            if (token.is_slider() && token.size == 2 && token.mirrored) {
                token.mirrored = false;
                changed = true;
            }
        if (changed && trace)
            normalizations.push_back({{"branch", index}, {"identity", "mirrored_X2_equals_X2"}});
        if (branch.control.phase != "unsupported_frame")
            upward_parts(branch, n);
    }
    BranchStep result{std::move(following), rule.id};
    result.evidence = {{"rule_id", rule.id},
                       {"metric_delta", metric_delta(predecessor, result.branches)}};
    if (trace) {
        result.evidence["catalog_transition"] = std::move(equality);
        result.evidence["normalizations"] = normalizations;
        result.evidence["successor_state"] = run_state_json(n, result.branches);
    }
    return result;
}

BranchStep framed_step(const Branch &old, Int n, const Catalog &catalog, const Limits &limits,
                       bool trace) {
    check_branch(old, limits);
    auto parts = framed_parts(old, n, limits.max_word_size);
    const auto &c = parts.control;
    const auto &u = parts.u;
    const auto &p = parts.p;
    const auto &tail = parts.tail;
    if (u.empty())
        throw ControlError("framed upward pass is already complete");
    auto i = c.i, j = c.j, k = c.k;
    bool mirror = c.mirrored;
    auto [out, inp] = anchors(i, j, k, mirror);
    const auto &policy = c.mirror_policy;
    bool bar_mirror = mirror && policy == "bar_mirror";
    if (out <= std::abs(u.back()) && std::abs(u.back()) < checked_add(out, i - 2)) {
        Branch next = old;
        auto &cc = next.control;
        --cc.left_length;
        cc.phase = cc.left_length ? "framed_zip_up" : "framed_upper_pass_complete";
        cc.output_payload.insert(cc.output_payload.begin(), u.back());
        framed_parts(next, n, limits.max_word_size);
        BranchStep result{{std::move(next)}, "control.absorb_internal_tape"};
        result.evidence = {{"rule_id", result.rule_id},
                           {"metric_delta", metric_delta(old, result.branches)}};
        if (trace) {
            result.evidence["equality_reference"] = "identical_literal_word";
            result.evidence["successor_state"] = run_state_json(n, result.branches);
        }
        return result;
    }
    Env env{{"n", n},
            {"i", i},
            {"j", mirror ? n - i - j + 2 : j},
            {"k", k},
            {"l", mirror ? n - std::abs(u.back()) : std::abs(u.back())},
            {"eps", (u.back() > 0 ? 1 : -1) * (bar_mirror ? -1 : 1)}};
    auto middle_core = core(n, i, j, k, mirror, policy);
    struct Candidate {
        const Rule *rule;
        Braid w;
        FrameTransport frame;
    };
    std::vector<Candidate> candidates;
    for (const auto *selected : catalog.select(env, "zip.", mirror ? ".M" : ".I")) {
        const Rule *rule = selected;
        if (bar_mirror)
            rule = &catalog.get(selected->id.substr(0, selected->id.size() - 2) + ".BM");
        auto lhs = expand_polynomial(rule->lhs, env, limits.expansion_budget);
        if (lhs.size() != 1)
            throw ControlError("contextual replay requires a monomial left side");
        const auto &occurrence = lhs.front().word;
        if (occurrence.size() < middle_core.size() ||
            !std::equal(middle_core.rbegin(), middle_core.rend(), occurrence.rbegin()))
            throw ControlError("catalog left side does not retain the framed core");
        auto w = braid_slice(occurrence, 0, occurrence.size() - middle_core.size());
        if (w.empty() || w.size() > u.size())
            continue;
        try {
            auto witness = transport_frame(slice(u, u.size() - w.size(), u.size()), n, i - 1, out,
                                           c.output_power, "left");
            if (witness.transformed == w)
                candidates.push_back({rule, std::move(w), std::move(witness)});
        } catch (const FrameError &) {
        }
    }
    if (candidates.size() != 1)
        throw ControlError(candidates.empty() ? "missing framed full-tape rule"
                                              : "competing framed rules");
    const auto &chosen = candidates.front();
    const auto &rule = *chosen.rule;
    const auto &before = chosen.frame;
    const auto &w = chosen.w;
    auto urest = slice(u, 0, u.size() - w.size());
    Braid payload;
    payload.reserve(c.output_payload.size());
    for (auto x : c.output_payload)
        payload.push_back((x > 0 ? 1 : -1) * checked_add(std::abs(x), before.end - out));
    auto leading = power(i - 1, before.end, c.output_power, limits.max_word_size);
    Word middle_word;
    append_braid(middle_word, urest);
    append_braid(middle_word, payload);
    append_braid(middle_word, leading);
    append_braid(middle_word, w);
    append(middle_word, middle_core);
    append_braid(middle_word, power(i - 1, inp, c.input_power, limits.max_word_size));
    append_braid(middle_word, p);
    append_braid(middle_word, tail);
    Branch middle{std::move(middle_word), old.coefficient, old.control};
    auto offset = Int(urest.size() + payload.size() + leading.size());
    auto raw = apply_rule(middle, rule, env, offset, limits);
    Json equality;
    if (trace)
        equality = contextual_transition(
            middle, raw, rule, env, offset, n,
            Json{{"controller", "framed_upward"}, {"marked_tape_transport", true}}, nullptr);
    const auto &base = catalog.get(rule.base_id + ".I");
    if (rule.rhs.size() != base.rhs.size() || raw.size() != rule.rhs.size())
        throw ControlError("catalog branch count disagrees with frame controls");
    std::vector<Branch> following;
    following.reserve(raw.size());
    Json witnesses = Json::array();
    for (std::size_t idx = 0; idx < rule.rhs.size(); ++idx) {
        const auto &term = rule.rhs[idx];
        auto [jj, kk, mm, aa, bb] = branch_shape(base.rhs[idx], env, mirror);
        if (bar_mirror) {
            aa = -aa;
            bb = -bb;
        }
        auto [out2, inp2] = anchors(i, jj, kk, mm);
        if (before.end != out2)
            throw ControlError("catalog output tape disagrees with frame transport");
        auto local = expand_word(term.at("word"), env, limits.expansion_budget);
        Word prefix = word_from_braid(power(i - 1, out2, aa, limits.max_word_size));
        append(prefix, core(n, i, jj, kk, mm, policy));
        append_braid(prefix, power(i - 1, inp2, bb, limits.max_word_size));
        if (!starts(local, prefix))
            throw ControlError("catalog branch is missing its retained flype frames");
        auto residual = braid_slice(local, prefix.size(), local.size());
        auto after = transport_frame(residual, n, i - 1, inp, c.input_power, "left");
        if (after.end != inp2 || after.transformed != residual)
            throw ControlError("catalog tail does not preserve the marked input tape");
        auto cc = framed_control(i, jj, kk, mm, checked_add(c.output_power, aa),
                                 checked_add(c.input_power, bb), urest.size(),
                                 residual.size() + p.size(), tail.size());
        cc.initial_slider = c.initial_slider;
        cc.mirror_policy = policy;
        cc.processed_events = c.processed_events;
        cc.processed_events.insert(cc.processed_events.begin(),
                                   ProcessedEvent{rule.id, env, Int(idx), policy});
        cc.output_payload = payload;
        auto word = expression(n, cc, urest, concat(residual, p), tail, limits.max_word_size);
        Branch branch{std::move(word), std::move(raw[idx].coefficient), std::move(cc)};
        check_branch(branch, limits);
        framed_parts(branch, n, limits.max_word_size);
        following.push_back(std::move(branch));
        if (trace)
            witnesses.push_back({{"tail_transport", frame_json(after)},
                                 {"emitted_output_power", aa},
                                 {"emitted_input_power", bb}});
    }
    BranchStep result{std::move(following), rule.id};
    result.evidence = {{"rule_id", rule.id}, {"metric_delta", metric_delta(old, result.branches)}};
    if (trace) {
        result.evidence["incoming_frame_transport"] = frame_json(before);
        result.evidence["catalog_transition"] = std::move(equality);
        result.evidence["payload_transport"] = {{"start", out},
                                                {"end", before.end},
                                                {"word", c.output_payload},
                                                {"transformed", payload}};
        result.evidence["branch_frames"] = witnesses;
        result.evidence["successor_state"] = run_state_json(n, result.branches);
        if (bar_mirror)
            result.evidence["semilinear_derivation"] = bound_rule_json(rule, env).at("semilinear");
    }
    return result;
}

} // namespace zipper
