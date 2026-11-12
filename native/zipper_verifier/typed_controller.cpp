#include "typed_controller.hpp"
#include "typed_transport.hpp"
#include <algorithm>
#include <set>
namespace zipper_verify::typed {
namespace {
std::vector<J> flatten_catalog(const J &node) {
    if (node.at("type") != "concat")
        return {node};
    std::vector<J> out;
    for (const auto &child : node.at("items")) {
        auto next = flatten_catalog(child);
        out.insert(out.end(), next.begin(), next.end());
    }
    return out;
}
bool fields(const J &v, std::initializer_list<const char *> names) {
    if (!v.is_object() || v.size() != names.size())
        return false;
    for (auto key : names)
        if (!v.contains(key))
            return false;
    return true;
}
bool pure(const Word &w) {
    return std::all_of(w.begin(), w.end(), [](const Token &t) { return is_int(t); });
}
const char *phase_name(ControlPhase p) {
    switch (p) {
    case ControlPhase::Up:
        return "zip_up";
    case ControlPhase::UpComplete:
        return "upper_pass_complete";
    case ControlPhase::FramedUp:
        return "framed_zip_up";
    case ControlPhase::FramedComplete:
        return "framed_upper_pass_complete";
    case ControlPhase::UnsupportedFrame:
        return "unsupported_frame";
    default:
        return "turnaround_complete";
    }
}
I extent(I i, I j, I k) {
    return checked_sub(checked_add(checked_add(i, j), k), 1);
}
J parameters_json(const std::map<std::string, I> &p) {
    J out = J::object();
    for (auto &[k, v] : p)
        out[k] = v;
    return out;
}
J changes(const Branch &old, const std::vector<Branch> &next, bool framed) {
    I words = 0, support = 0, max_word = 0, bits = 0, qmin = std::numeric_limits<I>::max(),
      qmax = std::numeric_limits<I>::min();
    if (next.empty())
        throw ControlError(framed ? "framed rule has no successors"
                                  : "upward rule has no surviving branches");
    for (const auto &b : next) {
        words = checked_add(words, length(b.word));
        support = checked_add(support, static_cast<I>(b.coefficient.size()));
        max_word = std::max(max_word, length(b.word));
        for (const auto &[q, c] : b.coefficient) {
            qmin = std::min(qmin, q);
            qmax = std::max(qmax, q);
            bits = std::max(bits, integer_bits(c));
        }
    }
    return {{"word_tokens", checked_sub(words, length(old.word))},
            {"coefficient_support", checked_sub(support, old.coefficient.size())},
            {"max_word_tokens", max_word},
            {"q_exponent_min", qmin},
            {"q_exponent_max", qmax},
            {"max_integer_bits", bits}};
}
struct Applied {
    std::vector<Branch> branches;
    J transition;
};
Applied apply(I n, const Branch &old, I offset, const J &rule, const J &env, const J &evidence,
              const std::vector<Control> *controls, Catalog &, const Limits &limits, bool trace) {
    check_state(n, {old}, limits);
    if (offset < 0)
        throw ReplayError("branch and offset must be nonnegative integers");
    if (!rule.contains("execution") || !rule["execution"].contains("algebra_verified") ||
        rule["execution"]["algebra_verified"] != J(true))
        throw ReplayError("rule lacks accepted algebraic verification");
    if (!condition(rule.at("condition"), env))
        throw ReplayError("rule guard is false");
    Polynomial lhs, rhs;
    try {
        lhs = expand_polynomial(rule.at("lhs"), env, J::object(), false, limits.expansion_budget);
        rhs = expand_polynomial(rule.at("rhs"), env, J::object(), false, limits.expansion_budget);
    } catch (const CatalogError &e) {
        std::string s = e.what();
        if (s.find("budget exceeded") != std::string::npos)
            throw LimitExceeded(s);
        throw ReplayError("invalid catalog instantiation: " + s);
    }
    if (lhs.size() != 1)
        throw ReplayError("contextual replay requires a monomial left side");
    if (lhs[0].q != 0 || lhs[0].coefficient != 1 || lhs[0].word.empty())
        throw ReplayError("left side must be a nonempty word with coefficient one");
    Word lhsword = decode_word(lhs[0].word);
    I end = checked_add(offset, length(lhsword));
    if (slice(old.word, offset, end) != lhsword)
        throw ReplayError("literal left side does not match selected occurrence");
    std::map<std::string, std::size_t> positions;
    std::vector<std::pair<Word, Laurent>> groups;
    for (const auto &term : rhs) {
        auto key = canonical_json(term.word);
        auto it = positions.find(key);
        if (it == positions.end()) {
            positions[key] = groups.size();
            groups.push_back({decode_word(term.word), {}});
            it = positions.find(key);
        }
        groups[it->second].second[term.q] += term.coefficient;
    }
    for (auto &g : groups)
        for (auto it = g.second.begin(); it != g.second.end();)
            if (it->second == 0)
                it = g.second.erase(it);
            else
                ++it;
    I count = 0;
    for (const auto &g : groups)
        if (!g.second.empty())
            ++count;
    if (controls && static_cast<I>(controls->size()) != count)
        throw ReplayError("one control update is required per distinct RHS word");
    Applied result;
    std::size_t k = 0;
    for (auto &g : groups)
        if (!g.second.empty()) {
            result.branches.push_back(
                {concat({slice(old.word, 0, offset), g.first, slice(old.word, end)}),
                 multiply_scalar(old.coefficient, g.second),
                 controls ? controls->at(k) : old.control});
            ++k;
        }
    check_state(n, result.branches, limits);
    if (trace) {
        J control_json = nullptr;
        if (controls) {
            control_json = J::array();
            for (auto &c : *controls)
                control_json.push_back(encode_control(c));
        }
        J req = request(0, offset, rule.at("id"), env, J::object(), evidence, control_json);
        J successors = J::array();
        for (auto &b : result.branches)
            successors.push_back(encode_branch(b));
        result.transition = {
            {"request", req},
            {"rule_sha256", digest(rule)},
            {"predecessor_state_sha256", digest(encode_state(n, {old}))},
            {"successor_state_sha256", digest(encode_state(n, result.branches))},
            {"predecessor", encode_branch(old)},
            {"successors", successors},
            {"equality_reference",
             {{"catalog_rule", rule.at("id")}, {"kind", "accepted_catalog_identity"}}}};
    }
    return result;
}
} // namespace
J encode_control(const Control &c) {
    if (c.phase == ControlPhase::UnsupportedFrame)
        return {{"phase", phase_name(c.phase)}, {"rule_id", c.rule_id}, {"message", c.message}};
    if (c.phase == ControlPhase::TurnaroundComplete)
        return {{"phase", "turnaround_complete"},
                {"i", c.i},
                {"j", c.j},
                {"mirrored", c.mirrored},
                {"mirror_policy", c.mirror_policy},
                {"output_power", c.output_power},
                {"input_power", c.input_power},
                {"exterior_length", c.exterior_length},
                {"payload_length", c.payload_length},
                {"routing_length", c.routing_length},
                {"transported_length", c.transported_length},
                {"processed_length", c.processed_length},
                {"tail_length", c.tail_length}};
    J out = {{"phase", phase_name(c.phase)},
             {"i", c.i},
             {"j", c.j},
             {"k", c.k},
             {"left_length", c.left_length},
             {"processed_length", c.processed_length},
             {"tail_length", c.tail_length}};
    if (c.phase == ControlPhase::Up || c.phase == ControlPhase::UpComplete)
        return out;
    out["mirrored"] = c.mirrored;
    out["output_power"] = c.output_power;
    out["input_power"] = c.input_power;
    out["output_payload"] = encode_word(c.output_payload);
    out["mirror_policy"] = c.mirror_policy;
    if (c.initial_slider) {
        auto &s = *c.initial_slider;
        out["initial_slider"] = {{"j", s.j},
                                 {"k", s.k},
                                 {"mirrored", s.mirrored},
                                 {"output_power", s.output_power},
                                 {"input_power", s.input_power},
                                 {"processed_length", s.processed_length},
                                 {"mirror_policy", s.mirror_policy}};
    }
    out["processed_events"] = J::array();
    for (auto &e : c.processed_events)
        out["processed_events"].push_back({{"rule_id", e.rule_id},
                                           {"parameters", parameters_json(e.parameters)},
                                           {"rhs_index", e.rhs_index},
                                           {"mirror_policy", e.mirror_policy}});
    return out;
}
Control decode_control(const J &v) {
    if (!v.is_object())
        throw ReplayError("branch control must be a JSON object");
    std::string phase = v.at("phase").get<std::string>();
    Control c;
    if (phase == "unsupported_frame") {
        c.phase = ControlPhase::UnsupportedFrame;
        c.rule_id = v.at("rule_id");
        c.message = v.at("message");
        return c;
    }
    if (phase == "turnaround_complete") {
        if (!fields(v, {"phase", "i", "j", "mirrored", "mirror_policy", "output_power",
                        "input_power", "exterior_length", "payload_length", "routing_length",
                        "transported_length", "processed_length", "tail_length"}) ||
            !v.at("mirrored").is_boolean() ||
            (v.at("mirror_policy") != "M" && v.at("mirror_policy") != "bar_mirror"))
            throw TurnaroundError("invalid turnaround state boundaries");
        for (auto k : {"i", "j", "output_power", "input_power", "exterior_length", "payload_length",
                       "routing_length", "transported_length", "processed_length", "tail_length"})
            if (!zipper_verify::is_int(v.at(k)))
                throw TurnaroundError("invalid turnaround state boundaries");
        c.phase = ControlPhase::TurnaroundComplete;
        c.i = zipper_verify::as_i(v.at("i"));
        c.j = zipper_verify::as_i(v.at("j"));
        c.mirrored = v.at("mirrored");
        c.mirror_policy = v.at("mirror_policy");
        c.output_power = zipper_verify::as_i(v.at("output_power"));
        c.input_power = zipper_verify::as_i(v.at("input_power"));
        c.exterior_length = zipper_verify::as_i(v.at("exterior_length"));
        c.payload_length = zipper_verify::as_i(v.at("payload_length"));
        c.routing_length = zipper_verify::as_i(v.at("routing_length"));
        c.transported_length = zipper_verify::as_i(v.at("transported_length"));
        c.processed_length = zipper_verify::as_i(v.at("processed_length"));
        c.tail_length = zipper_verify::as_i(v.at("tail_length"));
        return c;
    }
    const bool framed = phase == "framed_zip_up" || phase == "framed_upper_pass_complete";
    if (framed) {
        if (!fields(v, {"phase", "i", "j", "k", "mirrored", "output_power", "input_power",
                        "left_length", "processed_length", "tail_length", "initial_slider",
                        "processed_events", "output_payload", "mirror_policy"}) ||
            !v.at("mirrored").is_boolean() ||
            (v.at("mirror_policy") != "M" && v.at("mirror_policy") != "bar_mirror"))
            throw ControlError("invalid framed state boundaries");
        for (auto k : {"i", "j", "k", "output_power", "input_power", "left_length",
                       "processed_length", "tail_length"})
            if (!zipper_verify::is_int(v.at(k)))
                throw ControlError("invalid framed state boundaries");
        c.phase = phase == "framed_zip_up" ? ControlPhase::FramedUp : ControlPhase::FramedComplete;
        c.mirrored = v.at("mirrored");
        c.output_power = zipper_verify::as_i(v.at("output_power"));
        c.input_power = zipper_verify::as_i(v.at("input_power"));
        c.mirror_policy = v.at("mirror_policy");
        if (!v.at("output_payload").is_array())
            throw ControlError("invalid protected output-tape payload");
        try {
            c.output_payload = decode_word(v.at("output_payload"));
        } catch (const Error &) {
            throw ControlError("invalid protected output-tape payload");
        }
        const J &s = v.at("initial_slider");
        if (!fields(s, {"j", "k", "mirrored", "output_power", "input_power", "processed_length",
                        "mirror_policy"}) ||
            !s.at("mirrored").is_boolean() || s.at("mirror_policy") != v.at("mirror_policy") ||
            !v.at("processed_events").is_array())
            throw ControlError("invalid framed input-band record");
        for (auto k : {"j", "k", "output_power", "input_power", "processed_length"})
            if (!zipper_verify::is_int(s.at(k)))
                throw ControlError("invalid framed input-band record");
        c.initial_slider = InitialSlider{zipper_verify::as_i(s.at("j")),
                                         zipper_verify::as_i(s.at("k")),
                                         zipper_verify::as_i(s.at("output_power")),
                                         zipper_verify::as_i(s.at("input_power")),
                                         zipper_verify::as_i(s.at("processed_length")),
                                         s.at("mirrored").get<bool>(),
                                         s.at("mirror_policy").get<std::string>()};
        for (auto &e : v.at("processed_events")) {
            if (!fields(e, {"rule_id", "parameters", "rhs_index", "mirror_policy"}) ||
                !e.at("rule_id").is_string() || !e.at("parameters").is_object() ||
                e.at("mirror_policy") != v.at("mirror_policy") ||
                !zipper_verify::is_int(e.at("rhs_index")) ||
                zipper_verify::as_i(e.at("rhs_index")) < 0)
                throw ControlError("invalid processed-band event");
            ProcessedEvent event;
            event.rule_id = e.at("rule_id");
            event.rhs_index = zipper_verify::as_i(e.at("rhs_index"));
            event.mirror_policy = e.at("mirror_policy");
            for (auto it = e.at("parameters").begin(); it != e.at("parameters").end(); ++it) {
                if (!zipper_verify::is_int(it.value()))
                    throw ControlError("invalid processed-band event");
                event.parameters[it.key()] = zipper_verify::as_i(it.value());
            }
            c.processed_events.push_back(std::move(event));
        }
    } else {
        if (!fields(v, {"phase", "i", "j", "k", "left_length", "processed_length", "tail_length"}))
            throw ControlError("upward state boundaries are missing or unexpected");
        if (phase != "zip_up" && phase != "upper_pass_complete")
            throw ControlError("invalid upward state boundaries");
        for (auto k : {"i", "j", "k", "left_length", "processed_length", "tail_length"})
            if (!zipper_verify::is_int(v.at(k)))
                throw ControlError("invalid upward state boundaries");
        c.phase = phase == "zip_up" ? ControlPhase::Up : ControlPhase::UpComplete;
    }
    c.i = zipper_verify::as_i(v.at("i"));
    c.j = zipper_verify::as_i(v.at("j"));
    c.k = zipper_verify::as_i(v.at("k"));
    c.left_length = zipper_verify::as_i(v.at("left_length"));
    c.processed_length = zipper_verify::as_i(v.at("processed_length"));
    c.tail_length = zipper_verify::as_i(v.at("tail_length"));
    return c;
}
J encode_branch(const Branch &b) {
    return {{"word", encode_word(b.word)},
            {"coefficient", encoded_coefficient(b.coefficient)},
            {"control", encode_control(b.control)}};
}
Branch decode_branch(const J &b) {
    return {decode_word(b.at("word")), coefficient(b.at("coefficient"), false),
            decode_control(b.at("control"))};
}
J encode_state(I n, const std::vector<Branch> &b) {
    J branches = J::array();
    for (auto &x : b)
        branches.push_back(encode_branch(x));
    return {{"n", n}, {"branches", branches}};
}
Laurent multiply_scalar(const Laurent &a, const Laurent &b) {
    Laurent out;
    for (auto &[q, c] : a)
        for (auto &[r, d] : b)
            out[checked_add(q, r)] += c * d;
    for (auto it = out.begin(); it != out.end();)
        if (it->second == 0)
            it = out.erase(it);
        else
            ++it;
    return out;
}
void check_branch(const Branch &b, I rank, const Limits &limits) {
    if (b.coefficient.empty())
        throw ReplayError("a nonzero Laurent coefficient is required");
    validate_word(b.word, rank);
    if (length(b.word) > limits.max_word_size)
        throw LimitExceeded("word size limit exceeded");
    if (static_cast<I>(b.coefficient.size()) > limits.max_coefficient_support)
        throw LimitExceeded("coefficient support limit exceeded");
    for (auto &[q, c] : b.coefficient) {
        (void)q;
        if (integer_bits(c) > limits.max_integer_bits)
            throw LimitExceeded("coefficient integer bit limit exceeded");
    }
}
void check_state(I rank, const std::vector<Branch> &b, const Limits &limits) {
    if (rank < 1)
        throw ReplayError("rank must be a positive integer");
    for (auto &x : b) {
        if (x.coefficient.empty())
            throw ReplayError("a nonzero Laurent coefficient is required");
        validate_word(x.word, rank);
    }
    if (static_cast<I>(b.size()) > limits.max_terms)
        throw LimitExceeded("term limit exceeded");
    for (auto &x : b)
        check_branch(x, rank, limits);
}
std::pair<I, I> anchors(I i, I j, I k, bool m) {
    (void)i;
    return m ? std::make_pair(checked_sub(j, k), checked_add(j, 1))
             : std::make_pair(checked_add(checked_add(j, k), 1), j);
}
Word upward_friends(I size, I start, I friends) {
    Word out;
    if (size < 1 || friends < 0)
        throw ControlError("invalid friend interval");
    for (I a = friends; a > 0; --a) {
        I first = checked_add(start, a), stop = checked_add(first, checked_sub(size, 1));
        for (I x = first; x < stop; ++x)
            out.push_back(x);
    }
    return out;
}
Control upward_control(I i, I j, I k, I left, I processed, I tail) {
    Control c;
    c.phase = left ? ControlPhase::Up : ControlPhase::UpComplete;
    c.i = i;
    c.j = j;
    c.k = k;
    c.left_length = left;
    c.processed_length = processed;
    c.tail_length = tail;
    return c;
}
Branch upward_initialize(I n, const Word &left, I size, I start, I friends, const Word &processed,
                         const Word &tail) {
    if (size < 2 || size > n || start < 1 || friends < 0 || extent(size, start, friends) > n)
        throw ControlError("valid upward slider and friend intervals required");
    if (!pure(left) || !pure(processed) || !pure(tail))
        throw ControlError("exterior words must contain signed braid integers");
    Branch b{concat({left, upward_friends(size, start, friends), Word{slider(size, start)},
                     processed, tail}),
             {{0, Big(1)}},
             upward_control(size, start, friends, length(left), length(processed), length(tail))};
    validate_word(b.word, n);
    return b;
}
std::tuple<Control, Word, Word, Word> upward_parts(const Branch &b, I rank) {
    auto &c = b.control;
    I size = c.i, start = c.j, friends = c.k, a = c.left_length, p = c.processed_length,
      d = c.tail_length;
    if (c.phase != ControlPhase::Up && c.phase != ControlPhase::UpComplete)
        throw ControlError("invalid upward state boundaries");
    if (size < 2 || size > rank || start < 1 || friends < 0 ||
        extent(size, start, friends) > rank || std::min({a, p, d}) < 0)
        throw ControlError("invalid upward state boundaries");
    Word core = concat({upward_friends(size, start, friends), Word{slider(size, start)}});
    I end = checked_add(a, length(core));
    if (length(b.word) != checked_add(checked_add(end, p), d) || slice(b.word, a, end) != core ||
        c.phase != (a ? ControlPhase::Up : ControlPhase::UpComplete))
        throw ControlError("upward boundaries do not describe the full expression");
    Word u = slice(b.word, 0, a), pr = slice(b.word, end, checked_add(end, p)),
         t = slice(b.word, checked_add(end, p));
    if (!pure(u) || !pure(pr) || !pure(t))
        throw ControlError("upward exterior words must be pure braid words");
    return {c, u, pr, t};
}
Word upward_literal(const std::vector<J> &nodes, const J &p, I budget) {
    auto expanded =
        expand_word(J{{"type", "concat"}, {"items", nodes}}, p, J::object(), false, budget);
    if (expanded.empty())
        throw CatalogError("literal word expansion is empty");
    return decode_word(expanded.front().word);
}
Control upward_successor_control(const J &term, const J &env, I remaining, I processed, I tail,
                                 I budget, const std::string &id) {
    auto nodes = flatten_catalog(term.at("word"));
    auto it = std::find_if(nodes.begin(), nodes.end(),
                           [](const J &n) { return n.at("type") == "slider"; });
    if (it == nodes.end())
        throw ControlError("catalog branch has no slider");
    std::size_t at = std::distance(nodes.begin(), it);
    auto &token = nodes[at];
    I size = integer(token.at("size"), env), start = integer(token.at("start"), env);
    if (token.value("orientation", std::string("normal")) == "mirrored" && size > 2) {
        Control c;
        c.phase = ControlPhase::UnsupportedFrame;
        c.rule_id = id;
        c.message = "higher mirrored slider requires transport of retained half twists";
        return c;
    }
    I friends = 0;
    std::size_t before = at;
    if (at && nodes[at - 1].at("type") == "block") {
        auto &b = nodes[at - 1];
        if (integer(b.at("length"), env) == checked_sub(size, 1) &&
            integer(b.at("start"), env) == checked_add(start, 1) &&
            integer(b.value("sign", J(1)), env) == 1 && b.at("order") == "descending") {
            friends = integer(b.at("count"), env);
            before = at - 1;
        }
    }
    I left =
        length(upward_literal(std::vector<J>(nodes.begin(), nodes.begin() + before), env, budget));
    I emitted =
        length(upward_literal(std::vector<J>(nodes.begin() + at + 1, nodes.end()), env, budget));
    return upward_control(size, start, friends, checked_add(remaining, left),
                          checked_add(processed, emitted), tail);
}
std::pair<std::vector<Branch>, J> upward_step(I n, const Branch &old, Catalog &catalog,
                                              const Limits &limits, bool trace) {
    try {
        auto [c, u, p, t] = upward_parts(old, n);
        if (u.empty())
            throw ControlError("upward pass is already complete");
        I size = c.i;
        J env = {{"n", n},
                 {"i", size},
                 {"j", c.j},
                 {"k", c.k},
                 {"l", abs_i(as_i(u.back()))},
                 {"eps", as_i(u.back()) > 0 ? 1 : -1}};
        I end = checked_add(checked_add(length(u), checked_mul(checked_sub(size, 1), c.k)), 1);
        std::vector<std::pair<const J *, I>> candidates;
        std::vector<std::string> excludes =
            size == 2 ? std::vector<std::string>{"zip.flype."} : std::vector<std::string>{};
        for (auto &id : catalog.select_rule_ids(env, "zip.", ".I", excludes)) {
            auto &rule = catalog.get(id);
            auto lhs =
                expand_polynomial(rule.at("lhs"), env, J::object(), false, limits.expansion_budget);
            if (lhs.empty())
                throw CatalogError("upward rule has empty left side");
            Word occurrence = decode_word(lhs[0].word);
            I offset = checked_sub(end, length(occurrence));
            if (offset >= 0 && offset < length(u) && slice(old.word, offset, end) == occurrence)
                candidates.emplace_back(&rule, offset);
        }
        if (candidates.size() != 1) {
            if (candidates.empty())
                throw ControlError("missing full-tape rule");
            std::string message = "competing upward rules: ";
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                if (i)
                    message += ", ";
                message += candidates[i].first->at("id").get<std::string>();
            }
            throw ControlError(message);
        }
        auto [rule, offset] = candidates.front();
        std::map<std::string, std::size_t> positions;
        std::vector<J> terms;
        for (auto &term : rule->at("rhs")) {
            std::string key = canonical_json(
                encode_word(upward_literal({term.at("word")}, env, limits.expansion_budget)));
            auto it = positions.find(key);
            if (it == positions.end()) {
                positions.emplace(key, terms.size());
                terms.push_back(term);
            } else
                terms[it->second] = term;
        }
        std::vector<Control> controls;
        std::string id = rule->at("id");
        for (auto &term : terms)
            controls.push_back(upward_successor_control(term, env, offset, length(p), length(t),
                                                        limits.expansion_budget, id));
        auto applied = apply(n, old, offset, *rule, env,
                             {{"controller", "normal_upward"},
                              {"full_tape_literal_match", true},
                              {"i2_internal_precedes_flype", size == 2}},
                             &controls, catalog, limits, trace);
        J normalizations = J::array();
        for (std::size_t idx = 0; idx < applied.branches.size(); ++idx) {
            bool changed = false;
            for (auto &token : applied.branches[idx].word)
                if (auto *s = std::get_if<Slider>(&token); s && s->size == 2 && s->mirrored) {
                    *s = slider(2, s->start);
                    changed = true;
                }
            if (changed && trace)
                normalizations.push_back({{"branch", idx}, {"identity", "mirrored_X2_equals_X2"}});
        }
        for (auto &b : applied.branches)
            if (b.control.phase != ControlPhase::UnsupportedFrame)
                upward_parts(b, n);
        J event = {{"rule_id", id}, {"metric_delta", changes(old, applied.branches, false)}};
        if (trace) {
            event["catalog_transition"] = std::move(applied.transition);
            event["normalizations"] = std::move(normalizations);
            event["successor_state"] = encode_state(n, applied.branches);
        }
        return {std::move(applied.branches), std::move(event)};
    } catch (const CatalogError &e) {
        if (std::string(e.what()).find("budget") != std::string::npos)
            throw LimitExceeded(e.what());
        throw;
    }
}
Word frame_power(I size, I start, I exponent, I budget) {
    Big expanded = Big(size) * Big(checked_sub(size, 1));
    expanded /= 2;
    expanded *= Big(abs_i(exponent));
    if (expanded > Big(budget))
        throw LimitExceeded("expanded frame word limit exceeded");
    if (expanded == 0)
        return {};
    Word half;
    for (I width = checked_sub(size, 1); width > 0; --width)
        for (I off = 0; off < width; ++off)
            half.push_back(checked_add(start, off));
    if (exponent < 0)
        half = inverse(half);
    Word result;
    for (I copy = 0; copy < abs_i(exponent); ++copy)
        append(result, half);
    return result;
}
Word framed_core(I n, I i, I j, I k, bool mirrored, const std::string &policy) {
    I base = mirrored ? checked_add(checked_sub(checked_sub(n, i), j), 2) : j;
    Word word = upward_friends(i, base, k);
    if (mirrored)
        for (auto &x : word) {
            I v = checked_sub(n, as_i(x));
            x = policy == "bar_mirror" ? neg(v) : v;
        }
    word.push_back(slider(i, j, mirrored));
    return word;
}
Word friend_word(I n, const Control &c) {
    Word full = framed_core(n, c.i, c.j, c.k, c.mirrored, c.mirror_policy);
    return slice(full, 0, checked_sub(length(full), 1));
}
Control framed_control(I i, I j, I k, bool m, I a, I b, I left, I processed, I tail) {
    Control c = upward_control(i, j, k, left, processed, tail);
    c.phase = left ? ControlPhase::FramedUp : ControlPhase::FramedComplete;
    c.mirrored = m;
    c.output_power = a;
    c.input_power = b;
    return c;
}
Word framed_expression(I n, const Control &c, const Word &u, const Word &p, const Word &tail,
                       I budget) {
    I i = c.i, j = c.j, k = c.k, a = c.output_power, b = c.input_power;
    Big total = Big(length(u)) + Big(length(c.output_payload)) + Big(length(p)) + Big(length(tail));
    total += Big(checked_sub(i, 1)) * Big(k) + 1;
    Big half = Big(checked_sub(i, 1)) * Big(checked_sub(i, 2));
    half /= 2;
    total += half * (Big(abs_i(a)) + Big(abs_i(b)));
    if (total > Big(budget))
        throw LimitExceeded("expanded framed state word limit exceeded");
    auto [out, inp] = anchors(i, j, k, c.mirrored);
    return concat({u, c.output_payload, frame_power(i - 1, out, a, budget),
                   framed_core(n, i, j, k, c.mirrored, c.mirror_policy),
                   frame_power(i - 1, inp, b, budget), p, tail});
}
Branch framed_initialize(I n, const Word &left, I size, I start, I friends, bool mirrored, I a, I b,
                         const Word &processed, const Word &tail, const Word &payload,
                         const std::string &policy, const Limits &limits) {
    if ((policy != "M" && policy != "bar_mirror") || size < 2 || size > n || start < 1 ||
        friends < 0 || checked_sub(checked_add(start, size), 1) > n ||
        (mirrored ? checked_sub(start, friends) < 1
                  : checked_sub(checked_add(checked_add(start, size), friends), 1) > n))
        throw ControlError("valid framed slider, friends and integer powers required");
    if (size == 2) {
        if (mirrored || !payload.empty())
            throw ControlError("use the normal X2 convention for the two-strand delegate");
        return upward_initialize(n, left, size, start, friends, processed, tail);
    }
    if (!pure(left) || !pure(processed) || !pure(tail))
        throw ControlError("framed exterior must contain signed braid integers");
    auto [out, inp] = anchors(size, start, friends, mirrored);
    (void)inp;
    for (auto &x : payload)
        if (!is_int(x) || abs_i(as_i(x)) < out ||
            abs_i(as_i(x)) >= checked_sub(checked_add(out, size), 2))
            throw ControlError("protected output payload must be supported on the marked tape");
    Control c = framed_control(size, start, friends, mirrored, a, b, length(left),
                               length(processed), length(tail));
    c.mirror_policy = policy;
    c.output_payload = payload;
    c.initial_slider = InitialSlider{start, friends, a, b, length(processed), mirrored, policy};
    Branch result{
        framed_expression(n, c, left, processed, tail, limits.max_word_size), {{0, Big(1)}}, c};
    check_state(n, {result}, limits);
    return result;
}
std::tuple<Control, Word, Word, Word> framed_parts(const Branch &b, I n, I budget) {
    const Control &c = b.control;
    I i = c.i, j = c.j, k = c.k, a = c.left_length, p = c.processed_length, d = c.tail_length;
    if (c.phase != ControlPhase::FramedUp && c.phase != ControlPhase::FramedComplete)
        throw ControlError("invalid framed state boundaries");
    if ((c.mirror_policy != "M" && c.mirror_policy != "bar_mirror") || i < 3 || i > n || j < 1 ||
        k < 0 || checked_sub(checked_add(j, i), 1) > n ||
        (c.mirrored ? checked_sub(j, k) < 1
                    : checked_sub(checked_add(checked_add(j, i), k), 1) > n) ||
        std::min({a, p, d}) < 0)
        throw ControlError("invalid framed state boundaries");
    auto [out, inp] = anchors(i, j, k, c.mirrored);
    (void)inp;
    for (auto &x : c.output_payload)
        if (!is_int(x) || abs_i(as_i(x)) < out ||
            abs_i(as_i(x)) >= checked_sub(checked_add(out, i), 2))
            throw ControlError("invalid protected output-tape payload");
    if (!c.initial_slider || c.initial_slider->mirror_policy != c.mirror_policy)
        throw ControlError("invalid framed input-band record");
    for (auto &e : c.processed_events)
        if (e.rhs_index < 0 || e.mirror_policy != c.mirror_policy)
            throw ControlError("invalid processed-band event");
    if (c.phase != (a ? ControlPhase::FramedUp : ControlPhase::FramedComplete))
        throw ControlError("invalid framed phase");
    I begin = checked_sub(checked_sub(length(b.word), p), d);
    Word u = slice(b.word, 0, a), pr = slice(b.word, begin, checked_add(begin, p)),
         t = slice(b.word, checked_add(begin, p));
    if (!pure(u) || !pure(pr) || !pure(t))
        throw ControlError("framed exterior must be a pure braid word");
    if (b.word != framed_expression(n, c, u, pr, t, budget))
        throw ControlError("framed boundaries do not describe the complete expression");
    return {c, u, pr, t};
}
namespace {
std::tuple<I, I, bool, I, I> shape(const J &base, const J &env, bool mirrored) {
    auto nodes = flatten_catalog(base.at("word"));
    std::size_t at = 0;
    while (at < nodes.size() && nodes[at].at("type") != "slider")
        ++at;
    if (at == nodes.size())
        throw ControlError("catalog branch has no slider");
    auto &sl = nodes[at];
    I j = integer(sl.at("start"), env);
    bool mm = sl.value("orientation", std::string("normal")) == "mirrored";
    I k = 0;
    if (at && nodes[at - 1].at("type") == "block")
        k = integer(nodes[at - 1].at("count"), env);
    auto exponent = [](const J &node) -> I {
        if (node.at("type") == "half_twist")
            return 1;
        if (node.at("type") == "inverse" && node.at("word").at("type") == "half_twist")
            return -1;
        return 0;
    };
    I a = at ? exponent(nodes[0]) : 0, b = at + 1 < nodes.size() ? exponent(nodes[at + 1]) : 0;
    return {mirrored ? checked_add(checked_sub(checked_sub(zipper_verify::as_i(env.at("n")),
                                                           zipper_verify::as_i(env.at("i"))),
                                               j),
                                   2)
                     : j,
            k, mm != mirrored, a, b};
}
} // namespace
std::pair<std::vector<Branch>, J> framed_step(I n, const Branch &old, Catalog &catalog,
                                              const Limits &limits, bool trace) {
    try {
        check_state(n, {old}, limits);
        auto [c, u, p, tail] = framed_parts(old, n, limits.max_word_size);
        if (u.empty())
            throw ControlError("framed upward pass is already complete");
        I i = c.i, j = c.j, k = c.k;
        bool mirror = c.mirrored;
        auto [out, inp] = anchors(i, j, k, mirror);
        std::string policy = c.mirror_policy;
        bool bm = mirror && policy == "bar_mirror";
        I last = as_i(u.back());
        if (out <= abs_i(last) && abs_i(last) < checked_sub(checked_add(out, i), 2)) {
            Control cc = c;
            cc.left_length = checked_sub(length(u), 1);
            cc.phase = cc.left_length ? ControlPhase::FramedUp : ControlPhase::FramedComplete;
            cc.output_payload = concat({Word{u.back()}, c.output_payload});
            Branch b{old.word, old.coefficient, cc};
            framed_parts(b, n, limits.max_word_size);
            std::vector<Branch> next{b};
            J event = {{"rule_id", "control.absorb_internal_tape"},
                       {"metric_delta", changes(old, next, true)}};
            if (trace) {
                event["equality_reference"] = "identical_literal_word";
                event["successor_state"] = encode_state(n, next);
            }
            return {std::move(next), std::move(event)};
        }
        J env = {{"n", n},
                 {"i", i},
                 {"j", mirror ? checked_add(checked_sub(checked_sub(n, i), j), 2) : j},
                 {"k", k},
                 {"l", mirror ? checked_sub(n, abs_i(last)) : abs_i(last)},
                 {"eps", (last > 0 ? 1 : -1) * (bm ? -1 : 1)}};
        Word core = framed_core(n, i, j, k, mirror, policy);
        struct Candidate {
            std::shared_ptr<const J> derived;
            const J *ordinary;
            Word w;
            FrameWitness witness;
            const J &rule() const {
                return ordinary ? *ordinary : *derived;
            }
        };
        std::vector<Candidate> candidates;
        for (auto &id : catalog.select_rule_ids(env, "zip.", mirror ? ".M" : ".I")) {
            std::shared_ptr<const J> derived;
            const J &rule =
                bm ? *(derived = catalog.semilinear_rule_handle(id, env)) : catalog.get(id);
            auto lhs =
                expand_polynomial(rule.at("lhs"), env, J::object(), false, limits.expansion_budget);
            if (lhs.empty())
                throw ControlError("catalog left side does not retain the framed core");
            Word occurrence = decode_word(lhs[0].word);
            I boundary = checked_sub(length(occurrence), length(core));
            if (slice(occurrence, std::max<I>(0, boundary)) != core)
                throw ControlError("catalog left side does not retain the framed core");
            Word w = slice(occurrence, 0, boundary);
            if (w.empty() || length(w) > length(u))
                continue;
            try {
                auto witness = transport_frame(slice(u, checked_sub(length(u), length(w))), n,
                                               i - 1, out, c.output_power, "left");
                check_frame(witness);
                if (witness.transformed == w)
                    candidates.push_back(
                        {std::move(derived), bm ? nullptr : &rule, w, std::move(witness)});
            } catch (const FrameError &) {
            }
        }
        if (candidates.size() != 1)
            throw ControlError(candidates.empty() ? "missing framed full-tape rule"
                                                  : "competing framed rules");
        auto &chosen = candidates[0];
        auto &rule = chosen.rule();
        auto &w = chosen.w;
        auto &before = chosen.witness;
        Word rest = slice(u, 0, checked_sub(length(u), length(w)));
        I end = before.end, shift = checked_sub(end, out);
        Word payload;
        for (auto &x : c.output_payload) {
            I value = as_i(x), idx = checked_add(abs_i(value), shift);
            payload.push_back(value > 0 ? idx : neg(idx));
        }
        Word leading = frame_power(i - 1, end, c.output_power, limits.max_word_size);
        Branch middle{
            concat({rest, payload, leading, w, core,
                    frame_power(i - 1, inp, c.input_power, limits.max_word_size), p, tail}),
            old.coefficient, old.control};
        auto application = apply(
            n, middle, checked_add(checked_add(length(rest), length(payload)), length(leading)),
            rule, env, {{"controller", "framed_upward"}, {"marked_tape_transport", true}}, nullptr,
            catalog, limits, trace);
        auto &base = catalog.get(rule.at("base_id").get<std::string>() + ".I");
        if (rule.at("rhs").size() != base.at("rhs").size())
            throw ControlError("catalog branch counts differ");
        std::vector<Branch> next;
        J witnesses = J::array();
        for (std::size_t idx = 0; idx < rule.at("rhs").size(); ++idx) {
            auto &term = rule.at("rhs").at(idx), &base_term = base.at("rhs").at(idx);
            auto [jj, kk, mm, aa, bb] = shape(base_term, env, mirror);
            if (bm) {
                aa = neg(aa);
                bb = neg(bb);
            }
            auto [out2, inp2] = anchors(i, jj, kk, mm);
            if (end != out2)
                throw ControlError("catalog output tape disagrees with frame transport");
            Word local = upward_literal({term.at("word")}, env, limits.expansion_budget);
            Word prefix = concat({frame_power(i - 1, out2, aa, limits.max_word_size),
                                  framed_core(n, i, jj, kk, mm, policy)}),
                 right = frame_power(i - 1, inp2, bb, limits.max_word_size);
            I skip = checked_add(length(prefix), length(right));
            if (slice(local, 0, skip) != concat({prefix, right}))
                throw ControlError("catalog branch is missing its retained flype frames");
            Word residual = slice(local, skip);
            auto after = transport_frame(residual, n, i - 1, inp, c.input_power, "left");
            check_frame(after);
            if (after.end != inp2 || after.transformed != residual)
                throw ControlError("catalog tail does not preserve the marked input tape");
            Control cc = framed_control(i, jj, kk, mm, checked_add(c.output_power, aa),
                                        checked_add(c.input_power, bb), length(rest),
                                        checked_add(length(residual), length(p)), length(tail));
            cc.initial_slider = c.initial_slider;
            cc.mirror_policy = policy;
            cc.processed_events = c.processed_events;
            ProcessedEvent event;
            event.rule_id = rule.at("id");
            event.rhs_index = static_cast<I>(idx);
            event.mirror_policy = policy;
            for (auto it = env.begin(); it != env.end(); ++it)
                event.parameters[it.key()] = zipper_verify::as_i(it.value());
            cc.processed_events.insert(cc.processed_events.begin(), std::move(event));
            cc.output_payload = payload;
            next.push_back(
                {framed_expression(n, cc, rest, concat({residual, p}), tail, limits.max_word_size),
                 application.branches.at(idx).coefficient, cc});
            if (trace)
                witnesses.push_back({{"tail_transport", after.json()},
                                     {"emitted_output_power", aa},
                                     {"emitted_input_power", bb}});
        }
        check_state(n, next, limits);
        for (auto &b : next)
            framed_parts(b, n, limits.max_word_size);
        J event = {{"rule_id", rule.at("id")}, {"metric_delta", changes(old, next, true)}};
        if (trace) {
            event["incoming_frame_transport"] = before.json();
            event["catalog_transition"] = std::move(application.transition);
            event["payload_transport"] = {{"start", out},
                                          {"end", end},
                                          {"word", encode_word(c.output_payload)},
                                          {"transformed", encode_word(payload)}};
            event["branch_frames"] = witnesses;
            event["successor_state"] = encode_state(n, next);
            if (bm)
                event["semilinear_derivation"] = rule.at("semilinear");
        }
        return {std::move(next), std::move(event)};
    } catch (const CatalogError &e) {
        if (std::string(e.what()).find("budget") != std::string::npos)
            throw LimitExceeded(e.what());
        throw;
    }
}
Laurent add_scalar(const Laurent &a, const Laurent &b) {
    Laurent out = a;
    for (auto &[q, c] : b)
        out[q] += c;
    for (auto it = out.begin(); it != out.end();)
        if (it->second == 0)
            it = out.erase(it);
        else
            ++it;
    return out;
}
} // namespace zipper_verify::typed
