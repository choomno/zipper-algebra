// Independent owned controller. Catalog expansion and optional evidence use wire JSON.
#include "typed_state.hpp"
#include "typed_transport.hpp"
#include "catalog.hpp"
#include <algorithm>
#include <map>
namespace zipper_verify::typed {
namespace {
Transition result(std::vector<Work> children, J evidence) {
    Transition t;
    t.children = std::move(children);
    t.evidence = std::move(evidence);
    return t;
}
Transition unary_result(Work child, J evidence) {
    Transition t;
    t.children.push_back(std::move(child));
    t.evidence = std::move(evidence);
    return t;
}
Work common(const Work &w) {
    Work out;
    out.rank = w.rank;
    out.coefficient = w.coefficient;
    out.contexts = w.contexts;
    out.world_phi = w.world_phi;
    out.boundary_left = w.boundary_left;
    out.boundary_right = w.boundary_right;
    out.boundary_left_factors = w.boundary_left_factors;
    out.boundary_right_factors = w.boundary_right_factors;
    out.depth = w.depth;
    out.bottom_labels = w.bottom_labels;
    out.top_labels = w.top_labels;
    out.phase = Extract{};
    return out;
}
Work pass(Work w, bool down, Branch b) {
    w.coefficient = b.coefficient;
    if (down)
        w.phase = ZipDown{std::move(b)};
    else
        w.phase = ZipUp{std::move(b)};
    return w;
}
Work new_zip(const Work &base, const Word &left, const Word &right, I size, I start,
             const Limits &limits, bool mirrored = false, I a = 0, I b = 0,
             const Word &payload = {}) {
    Branch br;
    if (size == 2) {
        if (a || b)
            throw ControlObligation("two-strand framing must be the identity");
        br = upward_initialize(base.rank, concat({left, payload}), 2, start, 0, {}, right);
    } else
        br = framed_initialize(base.rank, left, size, start, 0, mirrored, a, b, {}, right, payload,
                               "bar_mirror", limits);
    br.coefficient = base.coefficient;
    return pass(base, false, std::move(br));
}
Transition single_step(const Work &w, Catalog &catalog, const Limits &limits) {
    const auto &s = std::get<Single>(w.phase);
    I n = w.rank, p = s.position;
    const Word &left = s.left, &processed = s.processed;
    if (left.empty()) {
        Work out = common(w);
        out.phase = Extract{concat({Word{slider(1, p)}, processed})};
        return unary_result(std::move(out), {{"rule_id", "control.single_scan_complete"}});
    }
    I letter = as_i(left.back()), index = abs_i(letter);
    bool adjacent = p == index || p == checked_add(index, 1),
         over = adjacent && ((p == index) == (letter > 0)), offensive = adjacent && !over;
    std::string sign = letter > 0 ? "positive" : "negative",
                side = p == checked_add(index, 1) ? "left" : "right",
                family = offensive  ? "foe"
                         : adjacent ? "close"
                                    : "far";
    std::string id = "single." + family + "." + (adjacent ? side + "." : "") + sign + ".I";
    J env = {{"n", n}, {"i", 1}, {"j", p}, {"k", 0}, {"l", index}};
    J instance = catalog.instantiate(id, env, limits.expansion_budget);
    J event = {{"rule_id", id},          {"parameters", env},
               {"active_label", n},      {"active_position", p},
               {"offensive", offensive}, {"active_over", adjacent ? J(over) : J(nullptr)}};
    if (!offensive) {
        Work out = w;
        auto &next = std::get<Single>(out.phase);
        next.left = slice(left, 0, length(left) - 1);
        next.processed = concat({Word{letter}, processed});
        if (adjacent)
            next.position = p == index ? checked_add(index, 1) : index;
        return unary_result(std::move(out), std::move(event));
    }
    struct Group {
        Word word;
        Laurent scalar;
    };
    std::vector<Group> groups;
    std::map<std::string, std::size_t> lookup;
    for (const J &term : instance.at("rhs")) {
        std::string key = canonical_json(term.at("word"));
        auto it = lookup.find(key);
        std::size_t at;
        if (it == lookup.end()) {
            at = groups.size();
            lookup.emplace(key, at);
            groups.push_back({decode_word(term.at("word")), {}});
        } else
            at = it->second;
        groups[at].scalar[zipper_verify::as_i(term.at("q_power"))] +=
            integer_value(term.at("integer"));
    }
    std::vector<Work> children;
    for (const auto &g : groups) {
        Work out = common(w);
        out.coefficient = multiply_scalar(w.coefficient, g.scalar);
        I at = 0;
        while (at < length(g.word) && is_int(g.word.at(at)))
            ++at;
        if (at == length(g.word))
            throw ControlObligation("single rewrite is missing its slider");
        const auto &coupon = std::get<Slider>(g.word.at(at));
        Word before = concat({slice(left, 0, length(left) - 1), slice(g.word, 0, at)}),
             after = concat({slice(g.word, at + 1), processed});
        if (coupon.size == 2)
            children.push_back(new_zip(out, before, after, 2, coupon.start, limits));
        else {
            out.phase = Single{std::move(before), std::move(after), p};
            children.push_back(std::move(out));
        }
    }
    return result(std::move(children), std::move(event));
}
bool outside_over(const Word &word, I n, I size, I start) {
    std::vector<bool> member;
    for (I p = 1; p <= n; ++p)
        member.push_back(start <= p && p < checked_add(start, size));
    for (auto it = word.rbegin(); it != word.rend(); ++it) {
        I letter = as_i(*it), at = checked_sub(abs_i(letter), 1);
        bool a = member.at(at), b = member.at(at + 1);
        if (a != b && !(letter > 0 ? a : b))
            return true;
        member.at(at) = b;
        member.at(at + 1) = a;
    }
    return false;
}
std::vector<Factor> reflected_factors(const std::vector<Factor> &f, I rank) {
    std::vector<Factor> out;
    for (auto it = f.rbegin(); it != f.rend(); ++it)
        out.push_back({it->growth_size, phi(it->word, rank)});
    return out;
}
Transition layer_step(const Work &w, const Limits &limits) {
    const auto &l = std::get<Layer>(w.phase);
    I n = w.rank;
    const Word &word = l.word;
    auto [left, coupon, right] = slider_parts(word);
    I size = coupon.size, start = coupon.start;
    bool mirrored = coupon.mirrored;
    if (coupon.reversed || coupon.barred)
        throw ControlObligation("layering requires an explicit normal or mirrored coupon without "
                                "reversed or barred flags");
    I plen = l.payload_length;
    Word unread = plen ? slice(left, 0, length(left) - plen) : left,
         payload = plen ? slice(left, length(left) - plen) : Word{};
    if (!outside_over(unread, n, size, start)) {
        Work out = common(w);
        if (l.completed_sides == 0) {
            out.phase = Layer{phi(word, n), 0, 1};
            out.world_phi = !w.world_phi;
            out.boundary_left = phi(w.boundary_right, n);
            out.boundary_right = phi(w.boundary_left, n);
            out.boundary_left_factors = reflected_factors(w.boundary_right_factors, n);
            out.boundary_right_factors = reflected_factors(w.boundary_left_factors, n);
            return unary_result(std::move(out), {{"rule_id", "control.layer_change_side"}});
        }
        out.phase = Extract{word};
        return unary_result(std::move(out), {{"rule_id", "control.layer_complete"}});
    }
    I letter = as_i(unread.back()), index = abs_i(letter), sign = letter > 0 ? 1 : -1;
    if (start <= index && index < checked_sub(checked_add(start, size), 1)) {
        Work out = w;
        std::get<Layer>(out.phase).payload_length = checked_add(plen, 1);
        return unary_result(std::move(out), {{"rule_id", "control.layer_absorb_internal"}});
    }
    if (index < checked_sub(start, 1) || index >= checked_add(start, size)) {
        Work out = common(w);
        out.phase = Layer{
            concat({slice(unread, 0, length(unread) - 1), payload, Word{coupon, letter}, right}),
            plen, l.completed_sides};
        return unary_result(std::move(out), {{"rule_id", "control.layer_far"}});
    }
    I move = index == checked_sub(checked_add(start, size), 1) ? 1 : -1;
    Word crossing;
    for (I k = 0; k < size; ++k) {
        I x = move == 1 ? checked_add(start, k)
                        : checked_sub(checked_sub(checked_add(start, size), 2), k);
        crossing.push_back(checked_mul(sign, x));
    }
    if (length(unread) < size || slice(unread, length(unread) - size) != crossing)
        throw ControlObligation("layering requires the complete active-band crossing");
    Word remaining = slice(unread, 0, length(unread) - size),
         shifted_payload = shifted(payload, move);
    bool offensive = (move == 1 && sign < 0) || (move == -1 && sign > 0);
    Work current = common(w);
    current.phase = Layer{
        concat({remaining, shifted_payload, Word{slider(size, checked_add(start, move), mirrored)},
                offensive ? bar(crossing) : crossing, right}),
        length(shifted_payload), l.completed_sides};
    J event = {{"rule_id", "control.layer_whole_band"},
               {"crossing", encode_word(crossing)},
               {"input_anchor", start},
               {"output_anchor", checked_add(start, move)},
               {"offensive", offensive},
               {"orientation", mirrored ? "mirrored" : "normal"}};
    if (!offensive)
        return unary_result(std::move(current), std::move(event));
    Laurent scalar, growth_scalar;
    I new_start, a, b;
    bool new_mirror;
    if (move == 1) {
        scalar = {{neg(size), 1}};
        growth_scalar = scalar;
        new_start = start;
        new_mirror = false;
        a = mirrored ? 1 : 0;
        b = mirrored ? -1 : 0;
    } else {
        scalar = {{size, 1}};
        growth_scalar = {{0, -1}};
        new_start = checked_sub(start, 1);
        new_mirror = true;
        a = mirrored ? 0 : -1;
        b = mirrored ? 0 : 1;
    }
    current.coefficient = multiply_scalar(w.coefficient, scalar);
    Work base = common(w);
    base.coefficient = multiply_scalar(w.coefficient, growth_scalar);
    Work growth = new_zip(base, remaining, right, checked_add(size, 1), new_start, limits,
                          new_mirror, a, b, shifted_payload);
    event["rule_id"] = move == 1 ? "layer.route.right" : "layer.route.left";
    event["growth_frames"] = {{"output_power", a}, {"input_power", b}};
    event["growth_child_type"] = {{"old_child_mirrored", mirrored},
                                  {"output_payload_mirrored", new_mirror != (a % 2 != 0)},
                                  {"input_payload_mirrored", new_mirror != (b % 2 != 0)}};
    event["equality_reference"] = {
        {"catalog_rule", move == 1 ? "layer.positive.I" : "layer.opposite_orientation.M"},
        {"full_band_garside_conjugation", move == 1 ? mirrored : !mirrored},
        {"scope", "braid-group-algebra slider recursion and whole-band naturality"}};
    std::vector<Work> children;
    children.push_back(std::move(growth));
    children.push_back(std::move(current));
    return result(std::move(children), std::move(event));
}
Transition zip_step(const Work &w, Catalog &catalog, const Limits &limits, bool trace) {
    I n = w.rank;
    bool down = std::holds_alternative<ZipDown>(w.phase);
    const Branch &br = down ? std::get<ZipDown>(w.phase).branch : std::get<ZipUp>(w.phase).branch;
    const Control &c = br.control;
    if (c.phase != ControlPhase::UpComplete && c.phase != ControlPhase::FramedComplete) {
        auto next = c.i == 2 ? upward_step(n, br, catalog, limits, trace)
                             : framed_step(n, br, catalog, limits, trace);
        std::vector<Work> children;
        for (auto &b : next.first)
            children.push_back(pass(common(w), down, std::move(b)));
        return result(std::move(children), std::move(next.second));
    }
    auto turned = turnaround(br, n, catalog, limits, trace, trace);
    auto [cc, exterior, protected_word, tail] =
        turnaround_parts(turned.first, n, limits.max_word_size);
    if (!protected_word.empty())
        throw ControlObligation("a zipper traversal must begin with an empty processed word");
    J event = {{"rule_id", "control.turnaround"},
               {"direction", down ? "zip_down" : "zip_up"},
               {"rotation_identity", "phi fixes slider orientation and sends each half twist to "
                                     "the canonical reflected half twist by Garside symmetry"},
               {"old_child_mirrored", cc.mirrored != (cc.output_power % 2 != 0)}};
    if (trace)
        event["turnaround"] = std::move(turned.second);
    Word payload = c.output_payload,
         friends = c.i == 2 ? upward_friends(c.i, c.j, c.k) : friend_word(n, c);
    I prefix = checked_add(length(payload), length(friends));
    if (slice(exterior, 0, prefix) != concat({payload, friends}))
        throw ControlObligation("turnaround exterior does not preserve its endpoint friend fan");
    I amount = c.mirrored ? c.k : neg(c.k);
    Word body = concat({shifted(payload, amount), slice(exterior, prefix)});
    Work out = common(w);
    if (!down) {
        append(out.boundary_left, friends);
        if (!friends.empty())
            out.boundary_left_factors.push_back({checked_add(cc.i, 1), friends});
        I start = checked_add(checked_sub(checked_sub(n, cc.i), cc.j), 2);
        Branch next =
            framed_initialize(n, phi(tail, n), cc.i, start, 0, cc.mirrored, cc.input_power,
                              cc.output_power, {}, phi(body, n), {}, c.mirror_policy, limits);
        next.coefficient = w.coefficient;
        return unary_result(pass(std::move(out), true, std::move(next)), std::move(event));
    }
    Word reflected = phi(friends, n);
    out.boundary_right = concat({reflected, out.boundary_right});
    if (!friends.empty())
        out.boundary_right_factors.insert(out.boundary_right_factors.begin(),
                                          {checked_add(cc.i, 1), reflected});
    out.phase = Layer{phi(concat({body, slice(turned.first.word, length(exterior))}), n), 0, 0};
    return unary_result(std::move(out), std::move(event));
}
} // namespace
std::tuple<Word, Slider, Word> slider_parts(const Word &word) {
    I at = -1;
    for (I i = 0; i < length(word); ++i)
        if (!is_int(word.at(i))) {
            if (at != -1)
                throw ControlObligation("the active local diagram must have exactly one slider");
            at = i;
        }
    if (at == -1)
        throw ControlObligation("the active local diagram must have exactly one slider");
    return {slice(word, 0, at), std::get<Slider>(word.at(at)), slice(word, at + 1)};
}
Word work_expression(const Work &w) {
    Word local;
    if (auto *s = std::get_if<Single>(&w.phase))
        local = concat({s->left, Word{slider(1, s->position)}, s->processed});
    else if (auto *s = std::get_if<ZipUp>(&w.phase))
        local = s->branch.word;
    else if (auto *s = std::get_if<ZipDown>(&w.phase))
        local = phi(s->branch.word, w.rank);
    else if (auto *s = std::get_if<Layer>(&w.phase))
        local = s->word;
    else
        local = std::get<Extract>(w.phase).word;
    local = concat({w.boundary_left, local, w.boundary_right});
    if (w.world_phi)
        local = phi(local, w.rank);
    return lift_contexts(local, w.contexts);
}
Work initial_work(I rank, const Word &word) {
    if (rank < 1)
        throw ReplayError("positive integer rank and literal word required");
    for (const auto &x : word)
        if (!is_int(x) || as_i(x) == 0 || abs_i(as_i(x)) >= rank)
            throw ReplayError("literal input must contain signed generators in the stated rank");
    Work w;
    w.rank = rank;
    w.coefficient = {{0, 1}};
    for (I i = 1; i <= rank; ++i) {
        w.bottom_labels.push_back(i);
        w.top_labels.push_back(i);
    }
    w.phase = Single{word, {}, rank};
    return w;
}
Transition algorithm_transition(const Work &work, Catalog &catalog, const Limits &limits,
                                bool trace) {
    static const Laurent one = {{0, 1}};
    const Branch *branch = nullptr;
    if (const auto *s = std::get_if<ZipUp>(&work.phase))
        branch = &s->branch;
    if (const auto *s = std::get_if<ZipDown>(&work.phase))
        branch = &s->branch;
    const bool unit_scalar = work.coefficient == one && (!branch || branch->coefficient == one);
    std::optional<Work> normalized;
    if (!unit_scalar) {
        normalized = work;
        normalized->coefficient = one;
        if (auto *s = std::get_if<ZipUp>(&normalized->phase))
            s->branch.coefficient = one;
        if (auto *s = std::get_if<ZipDown>(&normalized->phase))
            s->branch.coefficient = one;
    }
    const Work &unit = unit_scalar ? work : *normalized;
    Limits local = limits;
    local.max_terms = std::max<I>(4, limits.max_terms);
    local.max_coefficient_support = std::max<I>(2, limits.max_coefficient_support);
    Transition t;
    try {
        if (std::holds_alternative<Single>(unit.phase))
            t = single_step(unit, catalog, local);
        else if (std::holds_alternative<ZipUp>(unit.phase) ||
                 std::holds_alternative<ZipDown>(unit.phase))
            t = zip_step(unit, catalog, local, trace);
        else if (std::holds_alternative<Layer>(unit.phase))
            t = layer_step(unit, local);
        else
            t = extraction_transition(unit, local, trace);
    } catch (const LimitExceeded &e) {
        t.error = "incomplete";
        t.message = e.what();
        return t;
    } catch (const CatalogError &e) {
        t.message = e.what();
        t.error =
            t.message.find("budget") != std::string::npos ? "incomplete" : "control_obligation";
        return t;
    } catch (const Error &e) {
        t.error = "control_obligation";
        t.message = e.what();
        return t;
    }
    I depth = checked_add(work.depth, 1);
    for (auto &c : t.children) {
        c.coefficient = multiply_scalar(work.coefficient, c.coefficient);
        if (auto *s = std::get_if<ZipUp>(&c.phase))
            s->branch.coefficient = c.coefficient;
        if (auto *s = std::get_if<ZipDown>(&c.phase))
            s->branch.coefficient = c.coefficient;
        c.depth = depth;
    }
    for (auto &o : t.output)
        o.coefficient = multiply_scalar(work.coefficient, o.coefficient);
    return t;
}
} // namespace zipper_verify::typed
