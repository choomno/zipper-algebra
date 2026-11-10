// Independent typed component splitting and canonical slider normalization.
#include "typed_extraction.hpp"
#include <algorithm>
#include <set>
namespace zipper_verify::typed {
namespace {
I magnitude(I a) {
    if (a == std::numeric_limits<I>::min())
        throw RibbonError("signed braid generators inside the component required");
    return a < 0 ? -a : a;
}
Labels braid(const Word &a, I n) {
    Labels r;
    for (const auto &x : a) {
        if (!is_int(x) || magnitude(as_i(x)) < 1 || magnitude(as_i(x)) >= n)
            throw RibbonError("signed braid generators inside the component required");
        r.push_back(as_i(x));
    }
    return r;
}
Word word(const Labels &a) {
    Word r;
    for (I x : a)
        r.push_back(x);
    return r;
}
Labels sequence(I n) {
    Labels a;
    for (I i = 1; i <= n; i++)
        a.push_back(i);
    return a;
}
Word route(I n, const Labels &selected) {
    Labels ports = sequence(n);
    Word c;
    std::size_t dest = 0;
    for (I label : selected) {
        auto it = std::find(ports.begin(), ports.end(), label);
        std::size_t k = it - ports.begin();
        while (k > dest) {
            c.push_back(-static_cast<I>(k));
            std::swap(ports[k - 1], ports[k]);
            --k;
        }
        ++dest;
    }
    return reverse_word(c);
}
std::vector<Factor> boundary(const std::vector<Factor> &factors, I n) {
    std::set<I> sizes;
    for (const auto &f : factors) {
        I s = f.growth_size;
        if (s < 2 || s > n || !sizes.insert(s).second)
            throw RibbonError(
                "distinct growth sizes and explicit immutable boundary routes required");
        Labels ports = sequence(n);
        std::set<std::pair<I, I>> pairs;
        auto w = braid(f.word, n);
        for (auto it = w.rbegin(); it != w.rend(); ++it) {
            auto a = magnitude(*it) - 1;
            I x = ports[a], y = ports[a + 1];
            if (!pairs.insert(std::minmax(x, y)).second)
                throw RibbonError(
                    "an inherited boundary route must cross each strand pair at most once");
            std::swap(ports[a], ports[a + 1]);
        }
    }
    return factors;
}
Word flatten(const std::vector<Factor> &f) {
    Word r;
    for (const auto &a : f)
        append(r, a.word);
    return r;
}
I floorhalf(I x) {
    return x / 2 - (x < 0 && x % 2 != 0);
}
} // namespace
LayeredSplit split_layered(const Word &input, I n, const Labels &selected) {
    if (n < 1 || selected.empty())
        throw RibbonError("an increasing nonempty subset of the ambient bottom ports is required");
    I prev = 0;
    for (I x : selected) {
        if (x <= prev || x > n)
            throw RibbonError(
                "an increasing nonempty subset of the ambient bottom ports is required");
        prev = x;
    }
    I size = static_cast<I>(selected.size()), coupons = 0;
    std::vector<bool> m(static_cast<std::size_t>(n));
    for (I x : selected)
        m[x - 1] = true;
    Word component, complement;
    for (auto it = input.rbegin(); it != input.rend(); ++it) {
        const Token &t = *it;
        if (is_int(t)) {
            I g = as_i(t), at = magnitude(g) - 1, sign = g > 0 ? 1 : -1;
            if (at < 0 || at >= n - 1)
                throw RibbonError("ambient braid generator is out of range");
            bool l = m[at], r = m[at + 1];
            I before = std::count(m.begin(), m.begin() + at, true);
            if (l && r)
                component.push_back(sign * (before + 1));
            else if (!l && !r)
                complement.push_back(sign * (at - before + 1));
            else if (!(g > 0 ? l : r))
                throw RibbonError("the complement crosses over the selected component");
            m[at] = r;
            m[at + 1] = l;
        } else {
            const Slider &s = std::get<Slider>(t);
            if (s.size != size || s.start < 1 || s.start > n - size + 1)
                throw RibbonError(
                    "one explicit slider on the complete selected component is required");
            for (I k = s.start - 1; k < s.start + size - 1; k++)
                if (!m[k])
                    throw RibbonError("the slider does not join exactly the selected component");
            ++coupons;
            Slider nt = s;
            nt.start = 1;
            component.push_back(nt);
        }
    }
    if (coupons > 1 || (size > 1 && coupons != 1))
        throw RibbonError(
            "the selected component requires exactly one slider, except for a singleton");
    Labels top;
    for (I k = 0; k < n; k++)
        if (m[k])
            top.push_back(k + 1);
    component = reverse_word(component);
    complement = reverse_word(complement);
    Word br = route(n, selected), tr = inverse(route(n, top)), shift;
    for (const auto &x : complement) {
        I g = as_i(x);
        shift.push_back((g > 0 ? 1 : -1) * checked_add(magnitude(g), size));
    }
    Word rebuilt = concat({tr, component, shift, br});
    return {n,
            input,
            std::move(component),
            std::move(complement),
            std::move(br),
            std::move(tr),
            std::move(rebuilt),
            selected,
            std::move(top)};
}
CanonicalSlider normalize_canonical_slider(const Word &w, const std::vector<Factor> &top_routes,
                                           const std::vector<Factor> &bottom_routes, I images,
                                           I operations) {
    std::size_t at = 0, count = 0;
    for (std::size_t i = 0; i < w.size(); i++)
        if (!is_int(w[i])) {
            at = i;
            count++;
        }
    if (count != 1)
        throw RibbonError("one full-rank standard slider is required for twist collection");
    Slider t = std::get<Slider>(w[at]);
    if (t.size < 1 || t.start != 1 || t.reversed || t.barred)
        throw RibbonError("normal or mirrored full-rank slider at the first port required");
    I n = t.size;
    auto it = boundary(top_routes, n), ib = boundary(bottom_routes, n);
    Word tb = flatten(it), bb = flatten(ib);
    if (bb.size() > w.size())
        throw RibbonError("inherited boundary routes must be the exact exterior prefix and suffix");
    std::size_t end = w.size() - bb.size();
    if (tb.size() > at || end <= at || slice(w, 0, tb.size()) != tb || slice(w, end) != bb)
        throw RibbonError("inherited boundary routes must be the exact exterior prefix and suffix");
    Labels ports = sequence(n), output = ports;
    std::reverse(output.begin(), output.end());
    Word uw = word(braid(slice(w, tb.size(), at), n)),
         lw = inverse(word(braid(slice(w, at + 1, end), n)));
    if (t.mirrored) {
        Word half;
        for (I width = n - 1; width > 0; width--)
            for (I x = 1; x <= width; x++)
                half.push_back(x);
        append(uw, half);
        append(lw, half);
    }
    RibbonWitness u = recognize_nested(braid(uw, n), output, ports, images, operations),
                  l = recognize_nested(braid(lw, n), ports, ports, images, operations);
    Labels tf, bf, tw, cc;
    for (I k = 0; k < n - 1; k++) {
        I a = u.crossing_counts[k], b = l.crossing_counts[k];
        tf.push_back(k ? ((a % 2 + 2) % 2) : 0);
        bf.push_back(k ? ((b % 2 + 2) % 2) : 0);
        I power = checked_sub(k ? floorhalf(a) : a, k ? floorhalf(b) : b);
        tw.push_back(power);
        cc.push_back(k ? checked_mul(2, power) : power);
    }
    Word top = concat({tb, word(reconstruct_nested(output, tf, ports))}),
         bottom = concat({inverse(word(reconstruct_nested(ports, bf, ports))), bb});
    Slider normal{n, 1, false, false, false};
    Word body = word(reconstruct_nested(output, cc, ports));
    body.push_back(normal);
    Word normalized = concat({top, body, bottom});
    return {w,
            normal,
            std::move(u),
            std::move(l),
            std::move(tw),
            std::move(top),
            std::move(bottom),
            std::move(normalized),
            std::move(it),
            std::move(ib),
            images,
            operations};
}
J encode_split(const LayeredSplit &s) {
    return {{"n", s.n},
            {"word", encode_word(s.word)},
            {"bottom_ports", s.bottom_ports},
            {"top_ports", s.top_ports},
            {"component_word", encode_word(s.component_word)},
            {"complement_word", encode_word(s.complement_word)},
            {"bottom_route", encode_word(s.bottom_route)},
            {"top_route", encode_word(s.top_route)},
            {"reconstructed_word", encode_word(s.reconstructed_word)}};
}
J encode_canonical(const CanonicalSlider &s) {
    return {{"word", encode_word(s.word)},
            {"slider", encode_token(s.slider)},
            {"upper", encode_ribbon(s.upper)},
            {"lower_inverse", encode_ribbon(s.lower_inverse)},
            {"twist_powers", s.twist_powers},
            {"top_route", encode_word(s.top_route)},
            {"bottom_route", encode_word(s.bottom_route)},
            {"normalized_word", encode_word(s.normalized_word)},
            {"inherited_top", encode_factors(s.inherited_top)},
            {"inherited_bottom", encode_factors(s.inherited_bottom)},
            {"max_image_letters", s.max_image_letters},
            {"max_operations", s.max_operations}};
}
} // namespace zipper_verify::typed
