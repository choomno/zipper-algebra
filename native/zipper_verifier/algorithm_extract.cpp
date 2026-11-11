// Independent transcription of Python algorithm._extract_step and its helpers.
// This file uses verifier JSON states and independent ribbon recognition only.
#include "algorithm.hpp"
#include "extraction.hpp"
#include <algorithm>
#include <set>
#include <vector>

namespace zipper_verify {
namespace {

using Membership = std::vector<unsigned char>;

std::size_t crossing_position(const J &letter, const Membership &member) {
    I at = checked_sub(abs_i(as_i(letter)), 1);
    if (at < 0 || static_cast<std::uint64_t>(at) + 1 >= member.size())
        throw ControlObligation("endpoint crossing lies outside its active rank");
    return static_cast<std::size_t>(at);
}

J bottom_ports(const Word &word, I rank) {
    auto [left, token, right] = slider_parts(word);
    (void)left;
    I start = as_i(token.at(2));
    I end = checked_add(start, as_i(token.at(1)));
    if (rank < 1 || start < 1 || end > checked_add(rank, 1))
        throw ControlObligation("component slider lies outside its active rank");
    Membership member(static_cast<std::size_t>(rank), 0);
    for (I p = 1; p <= rank; ++p)
        member[static_cast<std::size_t>(p - 1)] = (start <= p && p < end);
    for (const auto &letter : right) {
        auto at = crossing_position(letter, member);
        std::swap(member[at], member[at + 1]);
    }
    J result = J::array();
    for (std::size_t p = 0; p < member.size(); ++p)
        if (member[p]) result.push_back(static_cast<I>(p + 1));
    return result;
}

J projected_factors(const J &factors, Membership &member) {
    J reverse_factors = J::array();
    for (auto factor = factors.rbegin(); factor != factors.rend(); ++factor) {
        Word chronological = J::array();
        const Word &word = factor->at("word");
        for (auto letter = word.rbegin(); letter != word.rend(); ++letter) {
            auto at = crossing_position(*letter, member);
            bool a = member[at] != 0, b = member[at + 1] != 0;
            if (a && b) {
                I offset = 1;
                for (std::size_t p = 0; p < at; ++p)
                    if (member[p]) offset = checked_add(offset, 1);
                chronological.push_back(as_i(*letter) > 0 ? offset : neg(offset));
            }
            std::swap(member[at], member[at + 1]);
        }
        Word projected = reverse_word(chronological);
        if (!projected.empty())
            reverse_factors.push_back({{"growth_size", factor->at("growth_size")},
                                       {"word", projected}});
    }
    return reverse_word(reverse_factors);
}

Big choose(I n, unsigned r) {
    if (n < 0) throw ControlObligation("component size cannot be negative");
    if (n < static_cast<I>(r)) return Big(0);
    Big value(1);
    for (unsigned k = 1; k <= r; ++k) {
        value *= Big(std::to_string(n - static_cast<I>(k) + 1));
        value /= k;
    }
    return value;
}

void check_generated_fans(const J &factors, I size) {
    std::set<I> tags;
    Big total(0);
    for (const auto &factor : factors) {
        const J &tag = factor.at("growth_size");
        if (!is_int(tag))
            throw ControlObligation("projected endpoint fan exceeds its recorded growth-size bound");
        I growth = as_i(tag);
        if (growth < 3 || growth > size || tags.count(growth))
            throw ControlObligation("projected endpoint fan exceeds its recorded growth-size bound");
        Big bound = Big(std::to_string(growth - 2)) * Big(std::to_string(size - growth + 1));
        Big count(std::to_string(length(factor.at("word"))));
        if (count > bound)
            throw ControlObligation("projected endpoint fan exceeds its recorded growth-size bound");
        tags.insert(growth);
        total += count;
    }
    if (total > choose(size, 3))
        throw ControlObligation("projected endpoint collar exceeds its size bound");
}

J reflected_factors(const J &factors, I rank) {
    J result = J::array();
    for (auto it = factors.rbegin(); it != factors.rend(); ++it) {
        J factor = *it;
        factor["word"] = phi(factor.at("word"), rank);
        result.push_back(std::move(factor));
    }
    return result;
}

Word flatten_factors(const J &factors) {
    Word result = J::array();
    for (const auto &factor : factors) append(result, factor.at("word"));
    return result;
}

J parity_flags(const J &ribbon) {
    J result = J::array();
    const auto &counts = ribbon.at("crossing_counts");
    for (I p = 1; p < length(counts); ++p) {
        I value = as_i(counts.at(static_cast<std::size_t>(p)));
        result.push_back(value % 2 == 0 ? 0 : 1);
    }
    return result;
}

J selected_labels(const J &labels, const J &ports, bool complement) {
    J result = J::array();
    if (!complement) {
        for (const auto &port : ports) {
            I p = as_i(port);
            if (p < 1 || p > length(labels))
                throw ControlObligation("component endpoint label lies outside its active rank");
            result.push_back(labels.at(static_cast<std::size_t>(p - 1)));
        }
        return result;
    }
    std::set<I> selected;
    for (const auto &port : ports) selected.insert(as_i(port));
    for (I p = 1; p <= length(labels); ++p)
        if (!selected.count(p)) result.push_back(labels.at(static_cast<std::size_t>(p - 1)));
    return result;
}

void trivial_record(J &record) {
    record["orientation"] = "normal";
    record["top_route"] = J::array();
    record["bottom_route"] = J::array();
    record["canonical_top_route"] = J::array();
    record["canonical_bottom_route"] = J::array();
    record["canonical_parities"] = {{"upper", J::array()}, {"lower_inverse", J::array()}};
}

} // namespace

Word lift_contexts(const Word &input, const J &contexts) {
    Word word = input;
    for (auto context = contexts.rbegin(); context != contexts.rend(); ++context)
        word = concat({context->at("top_route"), context->at("component_word"),
                       shifted(word, as_i(context->at("size"))), context->at("bottom_route")});
    return word;
}

J extraction_transition(const J &work, const Limits &limits, bool trace) {
    I n = as_i(work.at("rank"));
    Word body = work.at("word");
    Word local = concat({work.at("boundary_left"), body, work.at("boundary_right")});
    J top_factors = work.at("boundary_left_factors");
    J bottom_factors = work.at("boundary_right_factors");
    if (work.at("world_phi").get<bool>()) {
        local = phi(local, n);
        body = phi(body, n);
        J previous_top = std::move(top_factors);
        top_factors = reflected_factors(bottom_factors, n);
        bottom_factors = reflected_factors(previous_top, n);
    }
    auto [source_left, source_slider, source_right] = slider_parts(local);
    (void)source_left;
    (void)source_right;
    if (!source_slider.at(4).is_boolean() || source_slider.at(4).get<bool>() ||
        !source_slider.at(5).is_boolean() || source_slider.at(5).get<bool>())
        throw ControlObligation("component extraction requires a coupon without reversed or barred flags");
    J selected = bottom_ports(local, n);
    J split = split_layered(local, n, selected);
    check_split(split);
    Word component = split.at("component_word");
    I size = length(selected);
    J record = {{"bottom_ports", split.at("bottom_ports")},
                {"top_ports", split.at("top_ports")},
                {"original_bottom_ports", selected_labels(work.at("bottom_labels"), split.at("bottom_ports"), false)},
                {"original_top_ports", selected_labels(work.at("top_labels"), split.at("top_ports"), false)}};
    J canonical = nullptr;
    if (size == 1) {
        component = J::array();
        record["twist_powers"] = J::array();
        trivial_record(record);
    } else if (size == 2) {
        I power = 0;
        for (const auto &token : component)
            if (is_int(token)) power = checked_add(power, as_i(token));
        component = J::array();
        I count = abs_i(power);
        for (I p = 0; p < count; ++p) component.push_back(power >= 0 ? 1 : -1);
        component.push_back(slider(2, 1));
        record["twist_powers"] = J::array({power});
        trivial_record(record);
    } else {
        J inner_split = split_layered(body, n, bottom_ports(body, n));
        auto [whole_left, whole_slider, whole_right] = slider_parts(component);
        auto [inner_left, inner_slider, inner_right] = slider_parts(inner_split.at("component_word"));
        I boundary = checked_sub(length(whole_left), length(inner_left));
        if (whole_slider != inner_slider || slice(whole_left, boundary) != inner_left ||
            slice(whole_right, 0, length(inner_right)) != inner_right)
            throw ControlObligation("inherited endpoint routes do not delimit the component body");
        Word route_left = slice(whole_left, 0, boundary);
        Word route_right = slice(whole_right, length(inner_right));
        Membership member(static_cast<std::size_t>(n), 0);
        for (const auto &p : selected) member.at(static_cast<std::size_t>(as_i(p) - 1)) = 1;
        J inherited_bottom = projected_factors(bottom_factors, member);
        for (auto token = body.rbegin(); token != body.rend(); ++token) {
            if (is_int(*token)) {
                auto at = crossing_position(*token, member);
                std::swap(member[at], member[at + 1]);
            }
        }
        J inherited_top = projected_factors(top_factors, member);
        if (flatten_factors(inherited_top) != route_left || flatten_factors(inherited_bottom) != route_right)
            throw ControlObligation("projected endpoint fan factors disagree with the full component");
        check_generated_fans(inherited_top, size);
        check_generated_fans(inherited_bottom, size);
        J witness = normalize_canonical_slider(component, inherited_top, inherited_bottom,
                                               limits.expansion_budget, checked_mul(100, limits.expansion_budget));
        Big route_bound = choose(size, 3) + choose(size, 2) - 1;
        I route_length = std::max(length(witness.at("top_route")), length(witness.at("bottom_route")));
        if (Big(std::to_string(route_length)) > route_bound)
            throw ControlObligation("terminal finite route exceeds its size bound");
        component = witness.at("normalized_word");
        record["twist_powers"] = witness.at("twist_powers");
        record["inherited_top_route"] = route_left;
        record["inherited_bottom_route"] = route_right;
        record["orientation"] = witness.at("slider").at(3).get<bool>() ? "mirrored" : "normal";
        record["top_route"] = witness.at("top_route");
        record["bottom_route"] = witness.at("bottom_route");
        record["canonical_top_route"] = slice(witness.at("top_route"), length(route_left));
        record["canonical_bottom_route"] = slice(witness.at("bottom_route"), 0,
                                                 checked_sub(length(witness.at("bottom_route")), length(route_right)));
        record["canonical_parities"] = {{"upper", parity_flags(witness.at("upper"))},
                                        {"lower_inverse", parity_flags(witness.at("lower_inverse"))}};
        record["output_hierarchy"] = witness.at("upper").at("labels");
        record["input_hierarchy"] = witness.at("lower_inverse").at("labels");
        if (trace) canonical = std::move(witness);
    }
    J context = J::object();
    context["rank"] = n;
    context["size"] = size;
    context["component_word"] = std::move(component);
    context["top_route"] = split.at("top_route");
    context["bottom_route"] = split.at("bottom_route");
    context["record"] = record;
    J contexts = work.at("contexts");
    contexts.push_back(std::move(context));
    I remaining = checked_sub(n, size);
    J event = J::object();
    event["rule_id"] = "control.extract";
    event["component"] = std::move(record);
    event["rank_decrease"] = J::array({n, remaining});
    if (trace) {
        event["boundary_split"] = split;
        if (!canonical.is_null()) event["nested_slider"] = std::move(canonical);
    }
    if (remaining != 0) {
        J following = {{"rank", remaining}, {"coefficient", work.at("coefficient")},
                       {"world_phi", false},
                       {"boundary_left", J::array()}, {"boundary_right", J::array()},
                       {"boundary_left_factors", J::array()}, {"boundary_right_factors", J::array()},
                       {"depth", work.at("depth")},
                       {"bottom_labels", selected_labels(work.at("bottom_labels"), split.at("bottom_ports"), true)},
                       {"top_labels", selected_labels(work.at("top_labels"), split.at("top_ports"), true)},
                       {"phase", "single"}, {"left", split.at("complement_word")},
                       {"processed", J::array()}, {"position", remaining}};
        following["contexts"] = std::move(contexts);
        J result = J::object();
        result["children"] = J::array();
        result["children"].push_back(std::move(following));
        result["output"] = J::array();
        result["unit_evidence"] = std::move(event);
        return result;
    }
    J output = J::object();
    output["coefficient"] = work.at("coefficient");
    output["word"] = lift_contexts(J::array(), contexts);
    output["components"] = std::move(contexts);
    J result = J::object();
    result["children"] = J::array();
    result["output"] = J::array();
    result["output"].push_back(std::move(output));
    result["unit_evidence"] = std::move(event);
    return result;
}

} // namespace zipper_verify
