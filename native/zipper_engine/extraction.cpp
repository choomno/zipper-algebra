#include "extraction.hpp"
#include <algorithm>
#include <numeric>
#include <set>
#include <map>
#include <sstream>
#include <iomanip>
#include <openssl/sha.h>

namespace zipper {
namespace {
int64_t floor_half(int64_t x) {
    return x / 2 - (x < 0 && x % 2 != 0);
}
int64_t parity(int64_t x) {
    return (x % 2 + 2) % 2;
}
int64_t sign(int64_t x) {
    return x > 0 ? 1 : -1;
}
int64_t magnitude(int64_t x) {
    if (x == INT64_MIN)
        throw RibbonError("signed braid generators inside the component required");
    return x < 0 ? -x : x;
}
void validate_word(const Braid &word, int64_t n) {
    for (auto x : word)
        if (x == 0 || magnitude(x) >= n)
            throw RibbonError("signed braid generators inside the component required");
}
size_t locate(const Braid &v, int64_t x) {
    auto it = std::find(v.begin(), v.end(), x);
    if (it == v.end())
        throw RibbonError("distinct integer tree labels and the same boundary labels required");
    return static_cast<size_t>(it - v.begin());
}
void validate_ports(const Braid &tree, const Braid &ports) {
    std::set<int64_t> a(tree.begin(), tree.end()), b(ports.begin(), ports.end());
    if (tree.empty() || a.size() != tree.size() || tree.size() != ports.size() || a != b)
        throw RibbonError("distinct integer tree labels and the same boundary labels required");
    size_t lo = locate(ports, tree[0]), hi = lo;
    for (size_t i = 1; i < tree.size(); ++i) {
        auto p = locate(ports, tree[i]);
        lo = std::min(lo, p);
        hi = std::max(hi, p);
        if (hi - lo + 1 != i + 1)
            throw RibbonError("every marked subtree must occupy a contiguous boundary interval");
    }
}
struct Children {
    size_t start, width;
    bool old_left;
};
Children children(const Braid &ports, const Braid &tree, size_t size) {
    size_t lo = ports.size(), hi = 0;
    for (size_t i = 0; i < size - 1; ++i) {
        auto p = locate(ports, tree[i]);
        lo = std::min(lo, p);
        hi = std::max(hi, p);
    }
    if (hi - lo + 1 != size - 1)
        throw RibbonError("the marked child ribbon is not contiguous");
    size_t singleton = locate(ports, tree[size - 1]);
    if (singleton == hi + 1)
        return {lo, size - 1, true};
    if (singleton + 1 == lo)
        return {singleton, size - 1, false};
    throw RibbonError("the marked children are not adjacent");
}
Braid canonical_block(Braid &ports, const Braid &tree, size_t size, int64_t s) {
    auto c = children(ports, tree, size);
    Braid result;
    result.reserve(c.width);
    for (size_t k = 0; k < c.width; ++k) {
        size_t x = c.old_left ? c.start + c.width - k : c.start + 1 + k;
        result.push_back(s * static_cast<int64_t>(x));
        std::swap(ports[x - 1], ports[x]);
    }
    return result;
}
Braid inverse(const Braid &w) {
    Braid r;
    r.reserve(w.size());
    for (auto i = w.rbegin(); i != w.rend(); ++i) {
        magnitude(*i);
        r.push_back(-*i);
    }
    return r;
}
void gather_limits(int64_t letters, int64_t ops) {
    if (letters < 1 || ops < 1)
        throw RibbonError("positive integer ribbon gathering limits are required");
}
std::vector<Braid> artin_action(const Braid &word, int64_t n, int64_t max_letters, int64_t max_ops,
                                int64_t &operations) {
    if (n > max_letters)
        throw RibbonLimit("ribbon gathering Artin image limit exceeded");
    std::vector<Braid> images(static_cast<size_t>(n));
    for (int64_t i = 0; i < n; ++i)
        images[i] = {i + 1};
    size_t total = static_cast<size_t>(n);
    auto op = [&]() {
        if (operations >= max_ops)
            throw RibbonLimit("ribbon gathering operation limit exceeded");
        ++operations;
    };
    for (auto letter : word) {
        op();
        size_t at = static_cast<size_t>(magnitude(letter) - 1);
        const auto a = images[at], b = images[at + 1];
        Braid stack;
        auto emit = [&](int64_t x) {
            op();
            if (!stack.empty() && stack.back() == -x)
                stack.pop_back();
            else
                stack.push_back(x);
            if (stack.size() > static_cast<size_t>(max_letters))
                throw RibbonLimit("ribbon gathering Artin image limit exceeded");
        };
        if (letter > 0) {
            for (auto x : a)
                emit(x);
            for (auto x : b)
                emit(x);
            for (auto it = a.rbegin(); it != a.rend(); ++it)
                emit(-*it);
        } else {
            for (auto it = b.rbegin(); it != b.rend(); ++it)
                emit(-*it);
            for (auto x : a)
                emit(x);
            for (auto x : b)
                emit(x);
        }
        total = total - (letter > 0 ? b.size() : a.size()) + stack.size();
        if (total > static_cast<size_t>(max_letters))
            throw RibbonLimit("ribbon gathering Artin image limit exceeded");
        if (letter > 0) {
            images[at] = std::move(stack);
            images[at + 1] = a;
        } else {
            images[at] = b;
            images[at + 1] = std::move(stack);
        }
    }
    return images;
}
RibbonGathering gathering(const Braid &word, const Braid &canonical, int64_t n, int64_t letters,
                          int64_t ops) {
    gather_limits(letters, ops);
    int64_t operations = 0;
    auto original = artin_action(word, n, letters, ops, operations);
    auto target = artin_action(canonical, n, letters, ops, operations);
    if (original != target)
        throw RibbonError("pair-count candidate is not equal to the original braid");
    std::ostringstream s;
    s << '[';
    for (size_t i = 0; i < original.size(); ++i) {
        if (i)
            s << ',';
        s << '[';
        for (size_t j = 0; j < original[i].size(); ++j) {
            if (j)
                s << ',';
            s << original[i][j];
        }
        s << ']';
    }
    s << ']';
    auto bytes = s.str();
    unsigned char digest[SHA256_DIGEST_LENGTH];
    if (!SHA256(reinterpret_cast<const unsigned char *>(bytes.data()), bytes.size(), digest))
        throw EngineError("SHA256 computation failed");
    std::ostringstream hex;
    hex << std::hex << std::setfill('0');
    for (auto c : digest)
        hex << std::setw(2) << unsigned(c);
    return {letters, ops, hex.str()};
}
} // namespace
Braid NestedRibbon::twist_powers() const {
    Braid r = crossing_counts;
    for (size_t i = 1; i < r.size(); ++i)
        r[i] = floor_half(r[i]);
    return r;
}
Braid NestedRibbon::boundary_parities() const {
    Braid r;
    for (size_t i = 1; i < crossing_counts.size(); ++i)
        r.push_back(parity(crossing_counts[i]));
    return r;
}
Braid reconstruct_nested(const Braid &tree, const Braid &counts, const Braid &boundary) {
    validate_ports(tree, boundary);
    if (counts.size() + 1 != tree.size())
        throw RibbonError("one integer crossing count per non-leaf node required");
    Braid ports = boundary, chronological;
    for (size_t i = 0; i < counts.size(); ++i)
        for (int64_t j = 0; j < magnitude(counts[i]); ++j) {
            auto block = canonical_block(ports, tree, i + 2, sign(counts[i]));
            chronological.insert(chronological.end(), block.begin(), block.end());
        }
    std::reverse(chronological.begin(), chronological.end());
    return chronological;
}
NestedRibbon recognize_nested_strict(const Braid &word, const Braid &tree, const Braid &boundary) {
    validate_ports(tree, boundary);
    validate_word(word, tree.size());
    Braid chronological(word.rbegin(), word.rend()), ports = boundary, counts(tree.size() - 1, 0);
    std::map<int64_t, size_t> order;
    for (size_t i = 0; i < tree.size(); ++i)
        order[tree[i]] = i;
    std::vector<RibbonBlock> blocks;
    size_t offset = 0;
    while (offset < chronological.size()) {
        auto at = static_cast<size_t>(magnitude(chronological[offset]) - 1);
        auto size = std::max(order.at(ports[at]), order.at(ports[at + 1])) + 1;
        children(ports, tree, size);
        std::set<int64_t> old(tree.begin(), tree.begin() + size - 1), crossed;
        auto singleton = tree[size - 1], s = sign(chronological[offset]);
        if (chronological.size() - offset < size - 1)
            throw RibbonError("incomplete crossing of the marked child ribbon");
        Braid block(chronological.begin() + offset, chronological.begin() + offset + size - 1);
        for (auto letter : block) {
            auto left = static_cast<size_t>(magnitude(letter) - 1);
            auto a = ports[left], b = ports[left + 1], other = a == singleton ? b : a;
            if (sign(letter) != s || (a != singleton && b != singleton) || !old.count(other) ||
                crossed.count(other))
                throw RibbonError(
                    "a whole ribbon crossing must meet every child leaf once with one sign");
            crossed.insert(other);
            std::swap(ports[left], ports[left + 1]);
        }
        counts[size - 2] += s;
        blocks.push_back({static_cast<int64_t>(size), s, std::move(block)});
        offset += size - 1;
    }
    return {tree,
            boundary,
            word,
            std::move(blocks),
            counts,
            ports,
            reconstruct_nested(tree, counts, boundary),
            std::nullopt};
}
NestedRibbon recognize_nested(const Braid &word, const Braid &tree, const Braid &boundary,
                              int64_t letters, int64_t ops) {
    gather_limits(letters, ops);
    validate_ports(tree, boundary);
    validate_word(word, tree.size());
    try {
        return recognize_nested_strict(word, tree, boundary);
    } catch (const RibbonError &) {
    }
    std::map<int64_t, size_t> order;
    for (size_t i = 0; i < tree.size(); ++i)
        order[tree[i]] = i;
    Braid ports = boundary, counts;
    std::map<std::pair<size_t, size_t>, int64_t> pairs;
    for (auto it = word.rbegin(); it != word.rend(); ++it) {
        size_t at = magnitude(*it) - 1, a = order.at(ports[at]), b = order.at(ports[at + 1]);
        pairs[{std::min(a, b), std::max(a, b)}] += sign(*it);
        std::swap(ports[at], ports[at + 1]);
    }
    for (size_t node = 1; node < tree.size(); ++node) {
        auto count = pairs[{0, node}];
        for (size_t leaf = 1; leaf < node; ++leaf)
            if (pairs[{leaf, node}] != count)
                throw RibbonError(
                    "leaf-pair crossing counts do not define a marked ribbon candidate");
        counts.push_back(count);
    }
    auto canonical = reconstruct_nested(tree, counts, boundary);
    auto action = gathering(word, canonical, tree.size(), letters, ops);
    auto target = recognize_nested_strict(canonical, tree, boundary);
    return {tree, boundary, word, std::move(target.blocks), counts, ports, canonical, action};
}
namespace {
Braid over_route(Int n, const Braid &selected) {
    Braid ports(static_cast<size_t>(n)), chronological;
    std::iota(ports.begin(), ports.end(), 1);
    for (size_t destination = 0; destination < selected.size(); ++destination) {
        size_t current = locate(ports, selected[destination]);
        while (current > destination) {
            chronological.push_back(-static_cast<Int>(current));
            std::swap(ports[current - 1], ports[current]);
            --current;
        }
    }
    std::reverse(chronological.begin(), chronological.end());
    return chronological;
}
void boundary_routes(const std::vector<BoundaryRoute> &factors, Int n) {
    std::set<Int> sizes;
    for (const auto &factor : factors) {
        if (factor.growth_size < 2 || factor.growth_size > n ||
            !sizes.insert(factor.growth_size).second)
            throw RibbonError(
                "distinct growth sizes and explicit immutable boundary routes required");
        validate_word(factor.word, n);
        Braid ports(static_cast<size_t>(n));
        std::iota(ports.begin(), ports.end(), 0);
        std::set<std::pair<Int, Int>> pairs;
        for (auto it = factor.word.rbegin(); it != factor.word.rend(); ++it) {
            size_t at = magnitude(*it) - 1;
            auto pair = std::minmax(ports[at], ports[at + 1]);
            if (!pairs.insert(pair).second)
                throw RibbonError(
                    "an inherited boundary route must cross each strand pair at most once");
            std::swap(ports[at], ports[at + 1]);
        }
    }
}
} // namespace
LayeredSplit split_layered(const Word &diagram, Int n, const Braid &selected) {
    if (n < 1 || selected.empty() || !std::is_sorted(selected.begin(), selected.end()) ||
        std::adjacent_find(selected.begin(), selected.end()) != selected.end() ||
        selected.front() < 1 || selected.back() > n)
        throw RibbonError("an increasing nonempty subset of the ambient bottom ports is required");
    Int size = static_cast<Int>(selected.size());
    std::vector<bool> member(static_cast<size_t>(n), false);
    for (auto x : selected)
        member[x - 1] = true;
    Word component;
    Braid complement;
    size_t coupons = 0;
    for (auto it = diagram.rbegin(); it != diagram.rend(); ++it) {
        const auto &token = *it;
        if (!token.is_slider()) {
            Int letter = token.generator;
            if (letter == 0 || magnitude(letter) >= n)
                throw RibbonError("ambient braid generator is out of range");
            size_t at = magnitude(letter) - 1;
            bool left = member[at], right = member[at + 1];
            Int preceding = std::count(member.begin(), member.begin() + at, true);
            if (left && right)
                component.emplace_back(sign(letter) * (preceding + 1));
            else if (!left && !right)
                complement.push_back(sign(letter) * (static_cast<Int>(at) - preceding + 1));
            else if (!(letter > 0 ? left : right))
                throw RibbonError("the complement crosses over the selected component");
            member[at] = right;
            member[at + 1] = left;
        } else {
            if (token.size != size || token.start < 1 || token.start > n - size + 1)
                throw RibbonError(
                    "one explicit slider on the complete selected component is required");
            if (!std::all_of(member.begin() + token.start - 1,
                             member.begin() + token.start + size - 1, [](bool x) { return x; }))
                throw RibbonError("the slider does not join exactly the selected component");
            ++coupons;
            component.push_back(
                Token::slider(size, 1, token.mirrored, token.reversed, token.barred));
        }
    }
    if (coupons > 1 || (size > 1 && coupons != 1))
        throw RibbonError(
            "the selected component requires exactly one slider, except for a singleton");
    Braid top;
    for (Int i = 0; i < n; ++i)
        if (member[i])
            top.push_back(i + 1);
    std::reverse(component.begin(), component.end());
    std::reverse(complement.begin(), complement.end());
    auto bottom = over_route(n, selected), upper = inverse(over_route(n, top));
    Word rebuilt = word_from_braid(upper);
    append(rebuilt, component);
    for (auto x : complement)
        rebuilt.emplace_back(sign(x) * checked_add(magnitude(x), size));
    append(rebuilt, word_from_braid(bottom));
    return {n, diagram, selected, top, component, complement, bottom, upper, rebuilt};
}
CanonicalSlider normalize_canonical_slider(const Word &diagram,
                                           const std::vector<BoundaryRoute> &inherited_top,
                                           const std::vector<BoundaryRoute> &inherited_bottom,
                                           Int letters, Int ops) {
    size_t at = 0, sliders = 0;
    for (size_t i = 0; i < diagram.size(); ++i)
        if (diagram[i].is_slider()) {
            ++sliders;
            at = i;
        }
    if (sliders != 1)
        throw RibbonError("one full-rank standard slider is required for twist collection");
    auto token = diagram[at];
    if (token.size < 1 || token.start != 1 || token.reversed || token.barred)
        throw RibbonError("normal or mirrored full-rank slider at the first port required");
    Int n = token.size;
    boundary_routes(inherited_top, n);
    boundary_routes(inherited_bottom, n);
    Braid top_boundary, bottom_boundary;
    for (const auto &f : inherited_top)
        append(top_boundary, f.word);
    for (const auto &f : inherited_bottom)
        append(bottom_boundary, f.word);
    if (bottom_boundary.size() > diagram.size())
        throw RibbonError("inherited boundary routes must be the exact exterior prefix and suffix");
    size_t end = diagram.size() - bottom_boundary.size();
    if (top_boundary.size() > at || end <= at ||
        slice(diagram, 0, top_boundary.size()) != word_from_braid(top_boundary) ||
        slice(diagram, end, diagram.size()) != word_from_braid(bottom_boundary))
        throw RibbonError("inherited boundary routes must be the exact exterior prefix and suffix");
    Braid ports(static_cast<size_t>(n));
    std::iota(ports.begin(), ports.end(), 1);
    Braid output_tree(ports.rbegin(), ports.rend());
    Braid upper_word = braid_from_word(slice(diagram, top_boundary.size(), at));
    Braid lower_word = inverse(braid_from_word(slice(diagram, at + 1, end)));
    validate_word(upper_word, n);
    validate_word(lower_word, n);
    if (token.mirrored) {
        Braid half;
        for (Int width = n - 1; width > 0; --width)
            for (Int x = 1; x <= width; ++x)
                half.push_back(x);
        append(upper_word, half);
        append(lower_word, half);
    }
    auto upper = recognize_nested(upper_word, output_tree, ports, letters, ops);
    auto lower = recognize_nested(lower_word, ports, ports, letters, ops);
    Braid top_flags, bottom_flags;
    if (n > 1) {
        top_flags.push_back(0);
        append(top_flags, upper.boundary_parities());
        bottom_flags.push_back(0);
        append(bottom_flags, lower.boundary_parities());
    }
    Braid top = concat(top_boundary, reconstruct_nested(output_tree, top_flags, ports));
    Braid bottom = concat(inverse(reconstruct_nested(ports, bottom_flags, ports)), bottom_boundary);
    Braid twists = upper.twist_powers(), lower_twists = lower.twist_powers();
    for (size_t i = 0; i < twists.size(); ++i)
        twists[i] = checked_sub(twists[i], lower_twists[i]);
    Braid combined_counts = twists;
    for (size_t i = 1; i < combined_counts.size(); ++i)
        combined_counts[i] = checked_mul(2, combined_counts[i]);
    auto combined = reconstruct_nested(output_tree, combined_counts, ports);
    auto normal = Token::slider(n, 1);
    Word normalized = word_from_braid(concat(top, combined));
    normalized.push_back(normal);
    append(normalized, word_from_braid(bottom));
    return {diagram,       normal,           upper,   lower, twists, top, bottom, normalized,
            inherited_top, inherited_bottom, letters, ops};
}
SliderExtraction extract_slider(const Word &diagram, const Braid &upper_labels,
                                const Braid &lower_labels, const Braid &upper_ports,
                                const Braid &lower_ports) {
    size_t at = 0, sliders = 0;
    for (size_t i = 0; i < diagram.size(); ++i)
        if (diagram[i].is_slider()) {
            ++sliders;
            at = i;
        }
    if (sliders != 1)
        throw RibbonError("an isolated component with exactly one explicit slider is required");
    auto slider = diagram[at];
    if (slider.size < 1 || slider.start != 1)
        throw RibbonError("a full-rank slider anchored at the first component port is required");
    auto upper =
        recognize_nested(braid_from_word(slice(diagram, 0, at)), upper_labels, upper_ports);
    auto lower = recognize_nested(braid_from_word(slice(diagram, at + 1, diagram.size())),
                                  lower_labels, lower_ports);
    if (upper.labels.size() != static_cast<size_t>(slider.size) ||
        lower.labels.size() != static_cast<size_t>(slider.size))
        throw RibbonError("both marked trees must have the slider's rank");
    Word rebuilt = word_from_braid(upper.reconstructed_word);
    rebuilt.push_back(slider);
    append(rebuilt, word_from_braid(lower.reconstructed_word));
    return {diagram, slider, upper, lower, rebuilt};
}
Json extraction_json(const SliderExtraction &w) {
    return {{"word", word_json(w.word)},
            {"slider", word_json(Word{w.slider})[0]},
            {"upper", nested_json(w.upper)},
            {"lower", nested_json(w.lower)},
            {"reconstructed_word", word_json(w.reconstructed_word)}};
}
Json nested_json(const NestedRibbon &w) {
    Json blocks = Json::array();
    for (const auto &b : w.blocks)
        blocks.push_back({{"node_size", b.node_size},
                          {"sign", b.sign},
                          {"chronological_word", b.chronological_word}});
    Json action = nullptr;
    if (w.gathering)
        action = {{"max_image_letters", w.gathering->max_image_letters},
                  {"max_operations", w.gathering->max_operations},
                  {"action_sha256", w.gathering->action_sha256}};
    return {{"labels", w.labels},
            {"initial_ports", w.initial_ports},
            {"word", w.word},
            {"blocks", blocks},
            {"crossing_counts", w.crossing_counts},
            {"final_ports", w.final_ports},
            {"reconstructed_word", w.reconstructed_word},
            {"gathering", action}};
}
Json split_json(const LayeredSplit &w) {
    return {{"n", w.n},
            {"word", word_json(w.word)},
            {"bottom_ports", w.bottom_ports},
            {"top_ports", w.top_ports},
            {"component_word", word_json(w.component_word)},
            {"complement_word", w.complement_word},
            {"bottom_route", w.bottom_route},
            {"top_route", w.top_route},
            {"reconstructed_word", word_json(w.reconstructed_word)}};
}
Json canonical_json(const CanonicalSlider &w) {
    Json top = Json::array(), bottom = Json::array();
    for (const auto &f : w.inherited_top)
        top.push_back({{"growth_size", f.growth_size}, {"word", f.word}});
    for (const auto &f : w.inherited_bottom)
        bottom.push_back({{"growth_size", f.growth_size}, {"word", f.word}});
    return {{"word", word_json(w.word)},
            {"slider", word_json(Word{w.slider})[0]},
            {"upper", nested_json(w.upper)},
            {"lower_inverse", nested_json(w.lower_inverse)},
            {"twist_powers", w.twist_powers},
            {"top_route", w.top_route},
            {"bottom_route", w.bottom_route},
            {"normalized_word", word_json(w.normalized_word)},
            {"inherited_top", top},
            {"inherited_bottom", bottom},
            {"max_image_letters", w.max_image_letters},
            {"max_operations", w.max_operations}};
}
namespace {
bool equal_nested(const NestedRibbon &a, const NestedRibbon &b) {
    if (a.labels != b.labels || a.initial_ports != b.initial_ports || a.word != b.word ||
        a.crossing_counts != b.crossing_counts || a.final_ports != b.final_ports ||
        a.reconstructed_word != b.reconstructed_word || a.blocks.size() != b.blocks.size() ||
        bool(a.gathering) != bool(b.gathering))
        return false;
    for (size_t i = 0; i < a.blocks.size(); ++i)
        if (a.blocks[i].node_size != b.blocks[i].node_size ||
            a.blocks[i].sign != b.blocks[i].sign ||
            a.blocks[i].chronological_word != b.blocks[i].chronological_word)
            return false;
    return !a.gathering || (a.gathering->max_image_letters == b.gathering->max_image_letters &&
                            a.gathering->max_operations == b.gathering->max_operations &&
                            a.gathering->action_sha256 == b.gathering->action_sha256);
}
bool equal_routes(const std::vector<BoundaryRoute> &a, const std::vector<BoundaryRoute> &b) {
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i].growth_size != b[i].growth_size || a[i].word != b[i].word)
            return false;
    return true;
}
} // namespace
void check_nested(const NestedRibbon &w) {
    auto expected = recognize_nested(w.word, w.labels, w.initial_ports,
                                     w.gathering ? w.gathering->max_image_letters : 100000,
                                     w.gathering ? w.gathering->max_operations : 10000000);
    if (!equal_nested(expected, w))
        throw RibbonError("nested ribbon witness does not match its complete braid word");
}
void check_extraction(const SliderExtraction &w) {
    check_nested(w.upper);
    check_nested(w.lower);
    auto e = extract_slider(w.word, w.upper.labels, w.lower.labels, w.upper.initial_ports,
                            w.lower.initial_ports);
    if (e.word != w.word || e.slider != w.slider || !equal_nested(e.upper, w.upper) ||
        !equal_nested(e.lower, w.lower) || e.reconstructed_word != w.reconstructed_word)
        throw RibbonError("extraction witness does not reconstruct the complete component");
}
void check_split(const LayeredSplit &w) {
    auto e = split_layered(w.word, w.n, w.bottom_ports);
    if (e.n != w.n || e.word != w.word || e.bottom_ports != w.bottom_ports ||
        e.top_ports != w.top_ports || e.component_word != w.component_word ||
        e.complement_word != w.complement_word || e.bottom_route != w.bottom_route ||
        e.top_route != w.top_route || e.reconstructed_word != w.reconstructed_word)
        throw RibbonError("layered component split witness does not match its source diagram");
}
void check_canonical_slider(const CanonicalSlider &w) {
    check_nested(w.upper);
    check_nested(w.lower_inverse);
    auto e = normalize_canonical_slider(w.word, w.inherited_top, w.inherited_bottom,
                                        w.max_image_letters, w.max_operations);
    if (e.word != w.word || e.slider != w.slider || !equal_nested(e.upper, w.upper) ||
        !equal_nested(e.lower_inverse, w.lower_inverse) || e.twist_powers != w.twist_powers ||
        e.top_route != w.top_route || e.bottom_route != w.bottom_route ||
        e.normalized_word != w.normalized_word || !equal_routes(e.inherited_top, w.inherited_top) ||
        !equal_routes(e.inherited_bottom, w.inherited_bottom) ||
        e.max_image_letters != w.max_image_letters || e.max_operations != w.max_operations)
        throw RibbonError("canonical slider witness does not reconstruct its source diagram");
}
} // namespace zipper
