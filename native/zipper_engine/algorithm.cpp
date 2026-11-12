#include "algorithm.hpp"
#include <numeric>
#include <set>
#include <cstring>

namespace zipper {
namespace {
Work common(const Work &w) {
    Work x;
    x.rank = w.rank;
    x.depth = w.depth;
    x.coefficient = w.coefficient;
    x.contexts = w.contexts;
    x.world_phi = w.world_phi;
    x.boundary_left = w.boundary_left;
    x.boundary_right = w.boundary_right;
    x.boundary_left_factors = w.boundary_left_factors;
    x.boundary_right_factors = w.boundary_right_factors;
    x.bottom_labels = w.bottom_labels;
    x.top_labels = w.top_labels;
    return x;
}
Word phi(Word w, Int n) {
    std::reverse(w.begin(), w.end());
    for (auto &t : w)
        if (t.is_slider())
            t.start = n - t.size - t.start + 2;
        else
            t.generator = (t.generator > 0 ? 1 : -1) * (n - std::abs(t.generator));
    return w;
}
Braid phi(Braid w, Int n) {
    std::reverse(w.begin(), w.end());
    for (auto &x : w)
        x = (x > 0 ? 1 : -1) * (n - std::abs(x));
    return w;
}
Word shift(Word w, Int amount) {
    for (auto &t : w)
        if (t.is_slider())
            t.start = checked_add(t.start, amount);
        else
            t.generator = (t.generator > 0 ? 1 : -1) * checked_add(std::abs(t.generator), amount);
    return w;
}
Braid shift(Braid w, Int amount) {
    for (auto &x : w)
        x = (x > 0 ? 1 : -1) * checked_add(std::abs(x), amount);
    return w;
}
Word lift(Word word, const std::vector<Context> &contexts) {
    for (auto it = contexts.rbegin(); it != contexts.rend(); ++it) {
        Word next = word_from_braid(it->top_route);
        append(next, it->component_word);
        append(next, shift(word, it->size));
        append(next, word_from_braid(it->bottom_route));
        word = std::move(next);
    }
    return word;
}
Work single(const Work &c, const Braid &word) {
    auto w = common(c);
    w.phase = "single";
    w.left = word;
    w.position = w.rank;
    return w;
}
Work pass_work(const Work &c, const std::string &phase, Branch b) {
    auto w = common(c);
    w.phase = phase;
    w.coefficient = b.coefficient;
    w.branch = std::move(b);
    return w;
}
Work new_zip(const Work &c, const Braid &left, const Braid &right, Int size, Int start,
             const Limits &limits, bool mirrored = false, Int a = 0, Int b = 0,
             const Braid &payload = {}) {
    Branch branch;
    if (size == 2) {
        if (a || b)
            throw ControlError("two-strand framing must be the identity");
        branch = initialize_upward(c.rank, concat(left, payload), 2, start, 0, {}, right);
    } else
        branch = initialize_framed(c.rank, left, size, start, 0, mirrored, a, b, {}, right, payload,
                                   "bar_mirror", limits);
    branch.coefficient = c.coefficient;
    return pass_work(c, "zip_up", std::move(branch));
}
struct SliderParts {
    Braid left;
    Token slider;
    Braid right;
};
SliderParts slider_parts(const Word &word) {
    size_t at = 0, count = 0;
    for (size_t i = 0; i < word.size(); ++i)
        if (word[i].is_slider()) {
            at = i;
            ++count;
        }
    if (count != 1)
        throw ControlError("the active local diagram must have exactly one slider");
    return {braid_from_word(slice(word, 0, at)), word[at],
            braid_from_word(slice(word, at + 1, word.size()))};
}
Word joined(const Braid &left, const Token &slider, const Braid &right) {
    auto w = word_from_braid(left);
    w.push_back(slider);
    append(w, word_from_braid(right));
    return w;
}
Transition single_step(const Work &w, Catalog &catalog, const Limits &limits) {
    Transition r;
    if (w.left.empty()) {
        auto child = common(w);
        child.phase = "extract";
        child.word = joined({}, Token::slider(1, w.position), w.processed);
        r.children.push_back(child);
        r.event = {{"rule_id", "control.single_scan_complete"}};
        return r;
    }
    Int letter = w.left.back(), index = std::abs(letter);
    bool adjacent = w.position == index || w.position == index + 1;
    bool over = adjacent && ((w.position == index) == (letter > 0));
    bool offensive = adjacent && !over;
    std::string family = offensive ? "foe" : adjacent ? "close" : "far";
    std::string side = w.position == index + 1 ? "left" : "right";
    std::string rid = "single." + family + "." + (adjacent ? side + "." : "") +
                      (letter > 0 ? "positive" : "negative") + ".I";
    Env env{{"n", w.rank}, {"i", 1}, {"j", w.position}, {"k", 0}, {"l", index}};
    const Rule &rule = catalog.get(rid);
    if (!condition(rule.guard, env))
        throw CatalogError("rule guard is false");
    auto terms = expand_polynomial(rule.rhs, env, limits.expansion_budget);
    r.event = {{"rule_id", rid},         {"parameters", env},
               {"active_label", w.rank}, {"active_position", w.position},
               {"offensive", offensive}, {"active_over", adjacent ? Json(over) : Json(nullptr)}};
    if (!offensive) {
        auto child = w;
        child.left.pop_back();
        child.processed.insert(child.processed.begin(), letter);
        if (adjacent)
            child.position = w.position == index ? index + 1 : index;
        r.children.push_back(std::move(child));
        return r;
    }
    terms = serialize_order(std::move(terms));
    std::vector<std::pair<Word, Laurent>> grouped;
    for (const auto &t : terms) {
        auto it = std::find_if(grouped.begin(), grouped.end(),
                               [&](const auto &x) { return x.first == t.word; });
        if (it == grouped.end()) {
            grouped.push_back({t.word, {}});
            it = std::prev(grouped.end());
        }
        it->second[t.q] += t.coefficient;
    }
    for (auto &item : grouped) {
        auto c = common(w);
        c.coefficient = multiply(c.coefficient, item.second);
        auto p = slider_parts(item.first);
        auto left = slice(w.left, 0, w.left.size() - 1);
        append(left, p.left);
        auto right = concat(p.right, w.processed);
        if (p.slider.size == 2)
            r.children.push_back(new_zip(c, left, right, 2, p.slider.start, limits));
        else {
            auto child = single(c, left);
            child.processed = right;
            child.position = w.position;
            r.children.push_back(std::move(child));
        }
    }
    return r;
}
Transition zip_step(const Work &w, Catalog &catalog, const Limits &limits, bool trace,
                    bool want_evidence) {
    const auto &c = w.branch.control;
    Transition r;
    if (c.phase != "upper_pass_complete" && c.phase != "framed_upper_pass_complete") {
        auto step = c.i == 2 ? upward_step(w.branch, w.rank, catalog, limits, trace, want_evidence)
                             : framed_step(w.branch, w.rank, catalog, limits, trace, want_evidence);
        for (auto &b : step.branches)
            r.children.push_back(pass_work(w, w.phase, std::move(b)));
        r.event = step.evidence;
        if (!r.event.contains("rule_id"))
            r.event["rule_id"] = step.rule_id;
        return r;
    }
    auto turned = turnaround(w.branch, w.rank, catalog, limits, trace, trace);
    auto p = turnaround_parts(turned.first, w.rank, limits.max_word_size);
    const auto &cc = p.control;
    if (!p.p.empty())
        throw ControlError("a zipper traversal must begin with an empty processed word");
    r.event = {{"rule_id", "control.turnaround"},
               {"direction", w.phase},
               {"rotation_identity", "phi fixes slider orientation and sends each half twist to "
                                     "the canonical reflected half twist by Garside symmetry"},
               {"old_child_mirrored", bool(cc.mirrored ^ (cc.output_power % 2 != 0))}};
    if (trace)
        r.event["turnaround"] = turned.second;
    auto friend_braid = c.i == 2 ? friends(c.i, c.j, c.k) : friend_word(w.rank, c);
    auto prefix = concat(c.output_payload, friend_braid);
    if (p.u.size() < prefix.size() || !std::equal(prefix.begin(), prefix.end(), p.u.begin()))
        throw ControlError("turnaround exterior does not preserve its endpoint friend fan");
    auto body = shift(c.output_payload, c.mirrored ? c.k : -c.k);
    append(body, slice(p.u, prefix.size(), p.u.size()));
    auto next = common(w);
    if (w.phase == "zip_up") {
        append(next.boundary_left, friend_braid);
        if (!friend_braid.empty())
            next.boundary_left_factors.push_back({cc.i + 1, friend_braid});
        Int jj = w.rank - cc.i - cc.j + 2;
        auto initial =
            initialize_framed(w.rank, phi(p.tail, w.rank), cc.i, jj, 0, cc.mirrored, cc.input_power,
                              cc.output_power, {}, phi(body, w.rank), {},
                              c.mirror_policy == "" ? "bar_mirror" : c.mirror_policy, limits);
        initial.coefficient = w.coefficient;
        r.children.push_back(pass_work(next, "zip_down", std::move(initial)));
        return r;
    }
    next.boundary_right = concat(phi(friend_braid, w.rank), next.boundary_right);
    if (!friend_braid.empty())
        next.boundary_right_factors.insert(next.boundary_right_factors.begin(),
                                           {cc.i + 1, phi(friend_braid, w.rank)});
    auto local = word_from_braid(body);
    append(local, slice(turned.first.word, p.u.size(), turned.first.word.size()));
    next.phase = "layer";
    next.word = phi(local, w.rank);
    r.children.push_back(std::move(next));
    return r;
}
bool outside_over(const Braid &word, Int n, Int size, Int start) {
    std::vector<int> member(n);
    for (Int i = 0; i < n; ++i)
        member[i] = start <= i + 1 && i + 1 < start + size;
    for (auto it = word.rbegin(); it != word.rend(); ++it) {
        size_t at = std::abs(*it) - 1;
        auto a = member[at], b = member[at + 1];
        if (a != b && !(*it > 0 ? a : b))
            return true;
        std::swap(member[at], member[at + 1]);
    }
    return false;
}
std::vector<BoundaryRoute> phi_factors(const std::vector<BoundaryRoute> &f, Int n) {
    std::vector<BoundaryRoute> result;
    for (auto it = f.rbegin(); it != f.rend(); ++it)
        result.push_back({it->growth_size, phi(it->word, n)});
    return result;
}
Transition layer_step(const Work &w, const Limits &limits) {
    Transition r;
    auto p = slider_parts(w.word);
    auto s = p.slider;
    if (s.reversed || s.barred)
        throw ControlError("layering requires an explicit normal or mirrored coupon without "
                           "reversed or barred flags");
    size_t cut = p.left.size() - w.payload_length;
    auto unprocessed = slice(p.left, 0, cut), payload = slice(p.left, cut, p.left.size());
    if (!outside_over(unprocessed, w.rank, s.size, s.start)) {
        auto next = common(w);
        if (w.completed_sides == 0) {
            next.phase = "layer";
            next.word = phi(w.word, w.rank);
            next.completed_sides = 1;
            next.world_phi = !w.world_phi;
            next.boundary_left = phi(w.boundary_right, w.rank);
            next.boundary_right = phi(w.boundary_left, w.rank);
            next.boundary_left_factors = phi_factors(w.boundary_right_factors, w.rank);
            next.boundary_right_factors = phi_factors(w.boundary_left_factors, w.rank);
            r.event = {{"rule_id", "control.layer_change_side"}};
        } else {
            next.phase = "extract";
            next.word = w.word;
            r.event = {{"rule_id", "control.layer_complete"}};
        }
        r.children.push_back(std::move(next));
        return r;
    }
    Int letter = unprocessed.back(), index = std::abs(letter), sign = letter > 0 ? 1 : -1;
    if (s.start <= index && index < s.start + s.size - 1) {
        auto next = w;
        ++next.payload_length;
        r.children.push_back(next);
        r.event = {{"rule_id", "control.layer_absorb_internal"}};
        return r;
    }
    if (index < s.start - 1 || index >= s.start + s.size) {
        unprocessed.pop_back();
        append(unprocessed, payload);
        auto right = p.right;
        right.insert(right.begin(), letter);
        auto next = common(w);
        next.phase = "layer";
        next.word = joined(unprocessed, s, right);
        next.payload_length = w.payload_length;
        next.completed_sides = w.completed_sides;
        r.children.push_back(next);
        r.event = {{"rule_id", "control.layer_far"}};
        return r;
    }
    Int move = index == s.start + s.size - 1 ? 1 : -1;
    Braid crossing;
    if (move == 1)
        for (Int i = s.start; i < s.start + s.size; ++i)
            crossing.push_back(sign * i);
    else
        for (Int i = s.start + s.size - 2; i >= s.start - 1; --i)
            crossing.push_back(sign * i);
    if (unprocessed.size() < crossing.size() ||
        !std::equal(crossing.rbegin(), crossing.rend(), unprocessed.rbegin()))
        throw ControlError("layering requires the complete active-band crossing");
    unprocessed.resize(unprocessed.size() - crossing.size());
    auto shifted = shift(payload, move);
    bool offensive = (move == 1 && sign < 0) || (move == -1 && sign > 0);
    auto current = common(w);
    current.phase = "layer";
    current.payload_length = shifted.size();
    current.completed_sides = w.completed_sides;
    auto good = crossing;
    if (offensive)
        for (auto &x : good)
            x = -x;
    auto moved = s;
    moved.start += move;
    current.word = joined(concat(unprocessed, shifted), moved, concat(good, p.right));
    r.event = {{"rule_id", "control.layer_whole_band"},
               {"crossing", crossing},
               {"input_anchor", s.start},
               {"output_anchor", s.start + move},
               {"offensive", offensive},
               {"orientation", s.mirrored ? "mirrored" : "normal"}};
    if (!offensive) {
        r.children.push_back(current);
        return r;
    }
    Int start, a, b;
    bool mirror;
    Laurent scalar, growth;
    if (move == 1) {
        scalar = {{-s.size, 1}};
        growth = scalar;
        start = s.start;
        mirror = false;
        a = s.mirrored ? 1 : 0;
        b = s.mirrored ? -1 : 0;
    } else {
        scalar = {{s.size, 1}};
        growth = {{0, -1}};
        start = s.start - 1;
        mirror = true;
        a = s.mirrored ? 0 : -1;
        b = s.mirrored ? 0 : 1;
    }
    current.coefficient = multiply(w.coefficient, scalar);
    auto c = common(w);
    c.coefficient = multiply(w.coefficient, growth);
    r.children.push_back(
        new_zip(c, unprocessed, p.right, s.size + 1, start, limits, mirror, a, b, shifted));
    r.children.push_back(current);
    r.event["rule_id"] = move == 1 ? "layer.route.right" : "layer.route.left";
    r.event["growth_frames"] = {{"output_power", a}, {"input_power", b}};
    r.event["growth_child_type"] = {{"old_child_mirrored", s.mirrored},
                                    {"output_payload_mirrored", bool(mirror ^ (a % 2 != 0))},
                                    {"input_payload_mirrored", bool(mirror ^ (b % 2 != 0))}};
    r.event["equality_reference"] = {
        {"catalog_rule", move == 1 ? "layer.positive.I" : "layer.opposite_orientation.M"},
        {"full_band_garside_conjugation", move == 1 ? s.mirrored : !s.mirrored},
        {"scope", "braid-group-algebra slider recursion and whole-band naturality"}};
    return r;
}
Braid bottom_ports(const Word &word, Int n) {
    auto p = slider_parts(word);
    std::vector<int> member(n);
    for (Int i = 0; i < n; ++i)
        member[i] = p.slider.start <= i + 1 && i + 1 < p.slider.start + p.slider.size;
    for (auto x : p.right)
        std::swap(member[std::abs(x) - 1], member[std::abs(x)]);
    Braid result;
    for (Int i = 0; i < n; ++i)
        if (member[i])
            result.push_back(i + 1);
    return result;
}
std::vector<BoundaryRoute> project_factors(const std::vector<BoundaryRoute> &factors,
                                           std::vector<int> &member) {
    std::vector<BoundaryRoute> result;
    for (auto f = factors.rbegin(); f != factors.rend(); ++f) {
        Braid chronological;
        for (auto x = f->word.rbegin(); x != f->word.rend(); ++x) {
            size_t at = std::abs(*x) - 1;
            bool a = member[at], b = member[at + 1];
            if (a && b)
                chronological.push_back(
                    (*x > 0 ? 1 : -1) *
                    (std::accumulate(member.begin(), member.begin() + at, Int(0)) + 1));
            std::swap(member[at], member[at + 1]);
        }
        std::reverse(chronological.begin(), chronological.end());
        if (!chronological.empty())
            result.push_back({f->growth_size, chronological});
    }
    std::reverse(result.begin(), result.end());
    return result;
}
Int choose2(Int n) {
    return checked_mul(n, n - 1) / 2;
}
Int choose3(Int n) {
    return n < 3 ? 0 : checked_mul(checked_mul(n, n - 1), n - 2) / 6;
}
void check_fans(const std::vector<BoundaryRoute> &factors, Int size) {
    std::set<Int> tags;
    size_t total = 0;
    for (auto &f : factors) {
        Int growth = f.growth_size;
        if (growth < 3 || growth > size || tags.count(growth) ||
            f.word.size() > size_t(checked_mul(growth - 2, size - growth + 1)))
            throw ControlError("projected endpoint fan exceeds its recorded growth-size bound");
        tags.insert(growth);
        total += f.word.size();
    }
    if (total > size_t(choose3(size)))
        throw ControlError("projected endpoint collar exceeds its size bound");
}
Braid factor_word(const std::vector<BoundaryRoute> &f) {
    Braid result;
    for (auto &x : f)
        append(result, x.word);
    return result;
}
Braid select_labels(const Braid &labels, const Braid &ports, bool selected) {
    Braid result;
    for (size_t i = 0; i < labels.size(); ++i)
        if ((std::find(ports.begin(), ports.end(), Int(i + 1)) != ports.end()) == selected)
            result.push_back(labels[i]);
    return result;
}
Transition extract_step(const Work &w, const Limits &limits, bool trace) {
    Transition r;
    Word body = w.word, local = word_from_braid(w.boundary_left);
    append(local, body);
    append(local, word_from_braid(w.boundary_right));
    auto top = w.boundary_left_factors, bottom = w.boundary_right_factors;
    if (w.world_phi) {
        local = phi(local, w.rank);
        body = phi(body, w.rank);
        auto oldtop = top;
        top = phi_factors(bottom, w.rank);
        bottom = phi_factors(oldtop, w.rank);
    }
    auto source = slider_parts(local).slider;
    if (source.reversed || source.barred)
        throw ControlError(
            "component extraction requires a coupon without reversed or barred flags");
    auto selected = bottom_ports(local, w.rank);
    auto split = split_layered(local, w.rank, selected);
    check_split(split);
    auto component = split.component_word;
    Int size = selected.size();
    Json record = {
        {"bottom_ports", split.bottom_ports},
        {"top_ports", split.top_ports},
        {"original_bottom_ports", select_labels(w.bottom_labels, split.bottom_ports, true)},
        {"original_top_ports", select_labels(w.top_labels, split.top_ports, true)}};
    std::optional<CanonicalSlider> canonical;
    if (size <= 2) {
        Braid twists;
        if (size == 1)
            component.clear();
        else {
            Int power = 0;
            for (auto &t : component)
                if (!t.is_slider())
                    power = checked_add(power, t.generator);
            component = word_from_braid(Braid(std::abs(power), power >= 0 ? 1 : -1));
            component.push_back(Token::slider(2, 1));
            twists = {power};
        }
        record["twist_powers"] = twists;
        record["orientation"] = "normal";
        record["top_route"] = Json::array();
        record["bottom_route"] = Json::array();
        record["canonical_top_route"] = Json::array();
        record["canonical_bottom_route"] = Json::array();
        record["canonical_parities"] = {{"upper", Json::array()}, {"lower_inverse", Json::array()}};
    } else {
        auto inner = split_layered(body, w.rank, bottom_ports(body, w.rank));
        auto whole = slider_parts(component), inside = slider_parts(inner.component_word);
        if (whole.slider != inside.slider || whole.left.size() < inside.left.size() ||
            whole.right.size() < inside.right.size() ||
            !std::equal(inside.left.rbegin(), inside.left.rend(), whole.left.rbegin()) ||
            !std::equal(inside.right.begin(), inside.right.end(), whole.right.begin()))
            throw ControlError("inherited endpoint routes do not delimit the component body");
        auto route_left = slice(whole.left, 0, whole.left.size() - inside.left.size()),
             route_right = slice(whole.right, inside.right.size(), whole.right.size());
        std::vector<int> member(w.rank);
        for (auto x : selected)
            member[x - 1] = 1;
        auto inherited_bottom = project_factors(bottom, member);
        for (auto it = body.rbegin(); it != body.rend(); ++it)
            if (!it->is_slider())
                std::swap(member[std::abs(it->generator) - 1], member[std::abs(it->generator)]);
        auto inherited_top = project_factors(top, member);
        if (factor_word(inherited_top) != route_left ||
            factor_word(inherited_bottom) != route_right)
            throw ControlError("projected endpoint fan factors disagree with the full component");
        check_fans(inherited_top, size);
        check_fans(inherited_bottom, size);
        canonical = normalize_canonical_slider(component, inherited_top, inherited_bottom,
                                               limits.expansion_budget,
                                               checked_mul(100, limits.expansion_budget));
        const auto &witness = *canonical;
        Int bound = checked_add(choose3(size), choose2(size)) - 1;
        if (std::max(witness.top_route.size(), witness.bottom_route.size()) > size_t(bound))
            throw ControlError("terminal finite route exceeds its size bound");
        component = witness.normalized_word;
        record["twist_powers"] = witness.twist_powers;
        record["inherited_top_route"] = route_left;
        record["inherited_bottom_route"] = route_right;
        record["orientation"] = witness.slider.mirrored ? "mirrored" : "normal";
        record["top_route"] = witness.top_route;
        record["bottom_route"] = witness.bottom_route;
        record["canonical_top_route"] =
            slice(witness.top_route, route_left.size(), witness.top_route.size());
        record["canonical_bottom_route"] =
            slice(witness.bottom_route, 0, witness.bottom_route.size() - route_right.size());
        record["canonical_parities"] = {
            {"upper", witness.upper.boundary_parities()},
            {"lower_inverse", witness.lower_inverse.boundary_parities()}};
        record["output_hierarchy"] = witness.upper.labels;
        record["input_hierarchy"] = witness.lower_inverse.labels;
    }
    Context ctx{w.rank, size, component, split.top_route, split.bottom_route, record};
    auto contexts = w.contexts;
    contexts.push_back(ctx);
    Int remaining = w.rank - size;
    r.event = {{"rule_id", "control.extract"},
               {"component", record},
               {"rank_decrease", {w.rank, remaining}}};
    if (trace) {
        r.event["boundary_split"] = split_json(split);
        if (canonical)
            r.event["nested_slider"] = canonical_json(*canonical);
    }
    if (remaining) {
        Work next;
        next.rank = remaining;
        next.coefficient = w.coefficient;
        next.contexts = contexts;
        next.depth = w.depth;
        next.bottom_labels = select_labels(w.bottom_labels, split.bottom_ports, false);
        next.top_labels = select_labels(w.top_labels, split.top_ports, false);
        r.children.push_back(single(next, split.complement_word));
    } else
        r.output.push_back({w.coefficient, lift({}, contexts), contexts});
    return r;
}
Json context_json(const Context &c) {
    return {{"rank", c.rank},
            {"size", c.size},
            {"component_word", word_json(c.component_word)},
            {"top_route", c.top_route},
            {"bottom_route", c.bottom_route},
            {"record", c.record}};
}
Json contexts_json(const std::vector<Context> &contexts) {
    Json out = Json::array();
    for (auto &c : contexts)
        out.push_back(context_json(c));
    return out;
}
Json factors_json(const std::vector<BoundaryRoute> &factors) {
    Json out = Json::array();
    for (auto &f : factors)
        out.push_back({{"growth_size", f.growth_size}, {"word", f.word}});
    return out;
}
// Length-prefixed native structural keys. No JSON encoding of active states.
void key_int(std::string &key, Int v) {
    auto raw = static_cast<uint64_t>(v);
    for (int i = 0; i < 8; ++i)
        key.push_back(char((raw >> (i * 8)) & 255));
}
void key_text(std::string &key, const std::string &text) {
    key_int(key, text.size());
    key.append(text);
}
void key_braid(std::string &key, const Braid &word) {
    key_int(key, word.size());
    for (auto x : word)
        key_int(key, x);
}
void key_word(std::string &key, const Word &word) {
    key_int(key, word.size());
    for (auto &t : word) {
        key_int(key, t.generator);
        key_int(key, t.size);
        key_int(key, t.start);
        key.push_back(t.mirrored);
        key.push_back(t.reversed);
        key.push_back(t.barred);
    }
}
void key_json(std::string &key, const Json &j) {
    // JSON parsers distinguish positive unsigned values from native signed
    // values; mathematical integer metadata has one structural identity.
    key.push_back(char(j.is_number_integer() ? Json::value_t::number_integer : j.type()));
    if (j.is_object()) {
        key_int(key, j.size());
        for (auto it = j.begin(); it != j.end(); ++it) {
            key_text(key, it.key());
            key_json(key, it.value());
        }
    } else if (j.is_array()) {
        key_int(key, j.size());
        for (auto &x : j)
            key_json(key, x);
    } else if (j.is_string())
        key_text(key, j.get<std::string>());
    else if (j.is_number_integer())
        key_int(key, j.get<Int>());
    else if (j.is_boolean())
        key.push_back(j.get<bool>());
    else if (!j.is_null())
        throw EngineError("unsupported state metadata value");
}
void key_control(std::string &key, const Control &c) {
    key_text(key, c.phase);
    for (auto x : {c.i, c.j, c.k, c.left_length, c.processed_length, c.tail_length})
        key_int(key, x);
    if (c.phase == "zip_up" || c.phase == "upper_pass_complete")
        return;
    key.push_back(c.mirrored);
    key_int(key, c.output_power);
    key_int(key, c.input_power);
    key_braid(key, c.output_payload);
    key_text(key, c.mirror_policy);
    auto &s = c.initial_slider;
    for (auto x : {s.j, s.k, s.output_power, s.input_power, s.processed_length})
        key_int(key, x);
    key.push_back(s.mirrored);
    key_text(key, s.mirror_policy);
    key_int(key, c.processed_events.size());
    for (auto &e : c.processed_events) {
        key_text(key, e.rule_id);
        key_int(key, e.parameters.size());
        for (auto &p : e.parameters) {
            key_text(key, p.first);
            key_int(key, p.second);
        }
        key_int(key, e.rhs_index);
        key_text(key, e.mirror_policy);
    }
    for (auto x : {c.exterior_length, c.payload_length, c.routing_length, c.transported_length})
        key_int(key, x);
    key_text(key, c.rule_id);
    key_text(key, c.message);
}
void key_contexts(std::string &key, const std::vector<Context> &contexts) {
    key_int(key, contexts.size());
    for (auto &c : contexts) {
        key_int(key, c.rank);
        key_int(key, c.size);
        key_word(key, c.component_word);
        key_braid(key, c.top_route);
        key_braid(key, c.bottom_route);
        key_json(key, c.record);
    }
}
void key_factors(std::string &key, const std::vector<BoundaryRoute> &factors) {
    key_int(key, factors.size());
    for (auto &f : factors) {
        key_int(key, f.growth_size);
        key_braid(key, f.word);
    }
}
} // namespace
Work initial_work(Int n, const Braid &word) {
    Work w;
    w.rank = n;
    w.coefficient = {{0, 1}};
    for (Int i = 1; i <= n; ++i) {
        w.bottom_labels.push_back(i);
        w.top_labels.push_back(i);
    }
    return single(w, word);
}
Word work_expression(const Work &w) {
    Word local;
    if (w.phase == "single")
        local = joined(w.left, Token::slider(1, w.position), w.processed);
    else if (w.phase == "zip_up" || w.phase == "zip_down") {
        local = w.branch.word;
        if (w.phase == "zip_down")
            local = phi(local, w.rank);
    } else
        local = w.word;
    local = concat(word_from_braid(w.boundary_left), local);
    append(local, word_from_braid(w.boundary_right));
    if (w.world_phi)
        local = phi(local, w.rank);
    return lift(local, w.contexts);
}
size_t work_size(const Work &w) {
    size_t size = w.boundary_left.size() + w.boundary_right.size();
    if (w.phase == "single")
        size += w.left.size() + w.processed.size() + 1;
    else if (w.phase == "zip_up" || w.phase == "zip_down")
        size += w.branch.word.size();
    else
        size += w.word.size();
    for (auto &c : w.contexts)
        size += c.top_route.size() + c.bottom_route.size() + c.component_word.size();
    return size;
}
Transition algorithm_transition(const Work &work, Catalog &catalog, const Limits &limits,
                                bool trace, bool want_evidence) {
    auto unit = work;
    unit.coefficient = {{0, 1}};
    unit.branch.coefficient = unit.coefficient;
    auto local = limits;
    local.max_terms = std::max<Int>(4, local.max_terms);
    local.max_coefficient_support = std::max<Int>(2, local.max_coefficient_support);
    Transition r;
    try {
        if (unit.phase == "single")
            r = single_step(unit, catalog, local);
        else if (unit.phase == "zip_up" || unit.phase == "zip_down")
            r = zip_step(unit, catalog, local, trace, want_evidence);
        else if (unit.phase == "layer")
            r = layer_step(unit, local);
        else
            r = extract_step(unit, local, trace);
    } catch (const LimitExceeded &e) {
        r.error = "incomplete";
        r.message = e.what();
        return r;
    } catch (const RibbonLimit &e) {
        r.error = "incomplete";
        r.message = e.what();
        return r;
    } catch (const RibbonError &e) {
        r.error = "control_obligation";
        r.message = e.what();
        return r;
    } catch (const EngineError &e) {
        r.error = dynamic_cast<const CatalogError *>(&e) &&
                          std::string(e.what()).find("budget") != std::string::npos
                      ? "incomplete"
                      : "control_obligation";
        r.message = e.what();
        return r;
    }
    try {
        // A completed extraction still consumes one transition, even when it
        // emits only an output term and therefore has no child depth to store.
        const Int next_depth = checked_add(work.depth, 1);
        for (auto &child : r.children) {
            child.coefficient = multiply(work.coefficient, child.coefficient);
            if (child.phase == "zip_up" || child.phase == "zip_down")
                child.branch.coefficient = child.coefficient;
            child.depth = next_depth;
        }
        for (auto &output : r.output)
            output.coefficient = multiply(work.coefficient, output.coefficient);
    } catch (const LimitExceeded &error) {
        // Scaling belongs to the atomic transition, including overflow guards.
        Transition refused;
        refused.error = "incomplete";
        refused.message = error.what();
        return refused;
    }
    return r;
}
Json work_json(const Work &w) {
    Json j = {{"rank", w.rank},
              {"coefficient", coefficient_json(w.coefficient)},
              {"contexts", contexts_json(w.contexts)},
              {"world_phi", w.world_phi},
              {"boundary_left", w.boundary_left},
              {"boundary_right", w.boundary_right},
              {"boundary_left_factors", factors_json(w.boundary_left_factors)},
              {"boundary_right_factors", factors_json(w.boundary_right_factors)},
              {"depth", w.depth},
              {"bottom_labels", w.bottom_labels},
              {"top_labels", w.top_labels},
              {"phase", w.phase}};
    if (w.phase == "single") {
        j["left"] = w.left;
        j["processed"] = w.processed;
        j["position"] = w.position;
    } else if (w.phase == "zip_up" || w.phase == "zip_down")
        j["branch"] = branch_json(w.branch);
    else {
        j["word"] = word_json(w.word);
        if (w.phase == "layer") {
            j["payload_length"] = w.payload_length;
            j["completed_sides"] = w.completed_sides;
        }
    }
    return j;
}
Json output_json(const Output &o) {
    return {{"coefficient", coefficient_json(o.coefficient)},
            {"word", word_json(o.word)},
            {"components", contexts_json(o.components)}};
}
Json transition_json(const Transition &r) {
    if (!r.error.empty())
        return {{"error", r.error}, {"message", r.message}};
    Json children = Json::array(), output = Json::array();
    for (auto &c : r.children)
        children.push_back(work_json(c));
    for (auto &o : r.output)
        output.push_back(output_json(o));
    return {{"children", children}, {"output", output}, {"unit_evidence", r.event}};
}
std::string work_key(const Work &w) {
    std::string key;
    key_int(key, w.rank);
    key_contexts(key, w.contexts);
    key.push_back(w.world_phi);
    key_braid(key, w.boundary_left);
    key_braid(key, w.boundary_right);
    key_factors(key, w.boundary_left_factors);
    key_factors(key, w.boundary_right_factors);
    key_braid(key, w.bottom_labels);
    key_braid(key, w.top_labels);
    key_text(key, w.phase);
    if (w.phase == "single") {
        key_braid(key, w.left);
        key_braid(key, w.processed);
        key_int(key, w.position);
    } else if (w.phase == "zip_up" || w.phase == "zip_down") {
        key_word(key, w.branch.word);
        key_control(key, w.branch.control);
    } else {
        key_word(key, w.word);
        if (w.phase == "layer") {
            key_int(key, w.payload_length);
            key_int(key, w.completed_sides);
        }
    }
    return key;
}
std::string output_key(const Output &o) {
    std::string key;
    key_word(key, o.word);
    key_contexts(key, o.components);
    return key;
}
namespace {
void validate_saved_numbers(const Json &value) {
    if (value.is_number_float() ||
        (value.is_number_unsigned() && value.get<uint64_t>() > uint64_t(INT64_MAX)))
        throw EngineError("saved state metadata requires signed 64-bit integers");
    if (value.is_array() || value.is_object())
        for (auto it = value.begin(); it != value.end(); ++it) {
            // Laurent coefficients use their own arbitrary-precision parser.
            if (value.is_object() && it.key() == "coefficient")
                continue;
            validate_saved_numbers(it.value());
        }
}
void validate_saved_work(const Work &w) {
    if (w.rank < 1 || w.depth < 0 || w.coefficient.empty())
        throw EngineError("saved state requires positive rank and nonnegative depth");
    for (const Braid *labels : {&w.bottom_labels, &w.top_labels}) {
        if (labels->size() != size_t(w.rank))
            throw EngineError("saved endpoint label count differs from rank");
        std::set<Int> distinct(labels->begin(), labels->end());
        if (distinct.size() != labels->size() || *distinct.begin() < 1)
            throw EngineError("saved endpoint labels must be distinct and positive");
    }
    auto braid = [&](const Braid &b, Int rank) { validate_word(word_from_braid(b), rank); };
    braid(w.boundary_left, w.rank);
    braid(w.boundary_right, w.rank);
    for (const auto *factors : {&w.boundary_left_factors, &w.boundary_right_factors})
        for (const auto &f : *factors) {
            if (f.growth_size < 1)
                throw EngineError("invalid saved boundary factor");
            braid(f.word, w.rank);
        }
    for (const auto &c : w.contexts) {
        if (c.rank < 1 || c.size < 1 || c.size > c.rank)
            throw EngineError("invalid saved extraction context");
        validate_word(c.component_word, c.size);
        braid(c.top_route, c.rank);
        braid(c.bottom_route, c.rank);
        if (!c.record.is_object())
            throw EngineError("saved component record must be an object");
    }
    if (w.phase == "single") {
        braid(w.left, w.rank);
        braid(w.processed, w.rank);
        if (w.position < 1 || w.position > w.rank)
            throw EngineError("invalid saved active position");
    } else if (w.phase == "zip_up" || w.phase == "zip_down") {
        validate_word(w.branch.word, w.rank);
        const auto &c = w.branch.control;
        if (c.i < 2 || c.i > w.rank || c.j < 1 || c.j > w.rank || c.k < 0 || c.k > w.rank ||
            c.left_length < 0 || c.processed_length < 0 || c.tail_length < 0 ||
            c.left_length > Int(w.branch.word.size()) ||
            c.processed_length > Int(w.branch.word.size()) ||
            c.tail_length > Int(w.branch.word.size()))
            throw EngineError("invalid saved branch control boundaries");
        if (w.coefficient != w.branch.coefficient)
            throw EngineError("saved branch coefficient differs from term");
        if (c.i == 2)
            upward_parts(w.branch, w.rank);
        else {
            const Int tokens = static_cast<Int>(w.branch.word.size());
            if (c.output_power < -tokens || c.output_power > tokens || c.input_power < -tokens ||
                c.input_power > tokens)
                throw EngineError("saved frame powers exceed the literal branch length");
            framed_parts(w.branch, w.rank, tokens);
        }
    } else if (w.phase == "layer" || w.phase == "extract") {
        validate_word(w.word, w.rank);
        auto parts = slider_parts(w.word);
        if (w.phase == "layer" &&
            (w.payload_length < 0 || size_t(w.payload_length) > parts.left.size() ||
             (w.completed_sides != 0 && w.completed_sides != 1)))
            throw EngineError("invalid saved layering control");
    } else
        throw EngineError("unknown saved state phase");
}
} // namespace
Work work_from_json(const Json &j) {
    validate_saved_numbers(j);
    Work w;
    w.rank = j.at("rank");
    w.coefficient = coefficient_from_json(j.at("coefficient"));
    w.world_phi = j.at("world_phi");
    w.depth = j.at("depth");
    w.boundary_left = j.at("boundary_left").get<Braid>();
    w.boundary_right = j.at("boundary_right").get<Braid>();
    w.bottom_labels = j.at("bottom_labels").get<Braid>();
    w.top_labels = j.at("top_labels").get<Braid>();
    for (auto &c : j.at("contexts"))
        w.contexts.push_back({c.at("rank"), c.at("size"), word_from_json(c.at("component_word")),
                              c.at("top_route").get<Braid>(), c.at("bottom_route").get<Braid>(),
                              c.at("record")});
    for (auto &f : j.at("boundary_left_factors"))
        w.boundary_left_factors.push_back({f.at("growth_size"), f.at("word").get<Braid>()});
    for (auto &f : j.at("boundary_right_factors"))
        w.boundary_right_factors.push_back({f.at("growth_size"), f.at("word").get<Braid>()});
    w.phase = j.at("phase");
    if (w.phase == "single") {
        w.left = j.at("left").get<Braid>();
        w.processed = j.at("processed").get<Braid>();
        w.position = j.at("position");
    } else if (w.phase == "zip_up" || w.phase == "zip_down")
        w.branch = branch_from_json(j.at("branch"));
    else {
        w.word = word_from_json(j.at("word"));
        if (w.phase == "layer") {
            w.payload_length = j.at("payload_length");
            w.completed_sides = j.at("completed_sides");
        }
    }
    validate_saved_work(w);
    return w;
}
} // namespace zipper
