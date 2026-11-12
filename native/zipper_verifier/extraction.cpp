#include "extraction.hpp"
#include <algorithm>
#include <set>
namespace zipper_verify {
namespace {
J part(const J &a, size_t b, size_t e) {
  J r = J::array();
  for (size_t i = b; i < e; i++)
    r.push_back(a.at(i));
  return r;
}
J join(const J &a, const J &b) {
  J r = a;
  for (auto &x : b)
    r.push_back(x);
  return r;
}
J rev(J a) {
  std::reverse(a.begin(), a.end());
  return a;
}
I magnitude(I a) {
  if (a == std::numeric_limits<I>::min())
    throw RibbonError("signed braid generators inside the component required");
  return a < 0 ? -a : a;
}
J ribbon_inverse(const J &a) {
  J r = J::array();
  for (auto i = a.rbegin(); i != a.rend(); ++i)
    r.push_back(checked_sub(0, as_i(*i)));
  return r;
}
J braid(const J &a, I n) {
  for (auto &x : a)
    if (!is_int(x) || magnitude(as_i(x)) < 1 || magnitude(as_i(x)) >= n)
      throw RibbonError(
          "signed braid generators inside the component required");
  return a;
}
J sequence(I n) {
  J a = J::array();
  for (I i = 1; i <= n; i++)
    a.push_back(i);
  return a;
}
bool slider(const J &t) {
  return t.is_array() && t.size() == 6 && t[0] == "slider" && is_int(t[1]) &&
         is_int(t[2]) && t[3].is_boolean() && t[4].is_boolean() &&
         t[5].is_boolean();
}
J route(I n, const J &selected) {
  J ports = sequence(n), c = J::array();
  size_t dest = 0;
  for (auto &label : selected) {
    auto it = std::find(ports.begin(), ports.end(), label);
    size_t k = it - ports.begin();
    while (k > dest) {
      c.push_back(-static_cast<I>(k));
      std::swap(ports[k - 1], ports[k]);
      --k;
    }
    ++dest;
  }
  return rev(c);
}
J boundary(const J &factors, I n) {
  std::set<I> sizes;
  for (auto &f : factors) {
    I s = as_i(f.at("growth_size"));
    if (s < 2 || s > n || !sizes.insert(s).second || !f.at("word").is_array())
      throw RibbonError("distinct growth sizes and explicit immutable boundary "
                        "routes required");
    J ports = sequence(n);
    std::set<std::pair<I, I>> pairs;
    auto w = braid(f.at("word"), n);
    for (auto it = w.rbegin(); it != w.rend(); ++it) {
      auto a = magnitude(as_i(*it)) - 1;
      I x = as_i(ports[a]), y = as_i(ports[a + 1]);
      if (!pairs.insert(std::minmax(x, y)).second)
        throw RibbonError("an inherited boundary route must cross each strand "
                          "pair at most once");
      std::swap(ports[a], ports[a + 1]);
    }
  }
  return factors;
}
J flatten(const J &f) {
  J r = J::array();
  for (auto &a : f)
    r = join(r, a.at("word"));
  return r;
}
I floorhalf(I x) { return x / 2 - (x < 0 && x % 2 != 0); }
} // namespace
J split_layered(const J &word, I n, const J &selected) {
  if (n < 1 || !selected.is_array() || selected.empty())
    throw RibbonError("an increasing nonempty subset of the ambient bottom "
                      "ports is required");
  I prev = 0;
  for (auto &x : selected) {
    if (!is_int(x) || as_i(x) <= prev || as_i(x) > n)
      throw RibbonError("an increasing nonempty subset of the ambient bottom "
                        "ports is required");
    prev = as_i(x);
  }
  I size = selected.size(), coupons = 0;
  std::vector<bool> m(static_cast<size_t>(n));
  for (auto &x : selected)
    m[as_i(x) - 1] = true;
  J component = J::array(), complement = J::array();
  for (auto it = word.rbegin(); it != word.rend(); ++it) {
    const J &t = *it;
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
      if (!slider(t) || as_i(t[1]) != size || as_i(t[2]) < 1 ||
          as_i(t[2]) > n - size + 1)
        throw RibbonError("one explicit slider on the complete selected "
                          "component is required");
      for (I k = as_i(t[2]) - 1; k < as_i(t[2]) + size - 1; k++)
        if (!m[k])
          throw RibbonError(
              "the slider does not join exactly the selected component");
      ++coupons;
      J nt = t;
      nt[2] = 1;
      component.push_back(nt);
    }
  }
  if (coupons > 1 || (size > 1 && coupons != 1))
    throw RibbonError("the selected component requires exactly one slider, "
                      "except for a singleton");
  J top = J::array();
  for (I k = 0; k < n; k++)
    if (m[k])
      top.push_back(k + 1);
  component = rev(component);
  complement = rev(complement);
  J br = route(n, selected), tr = ribbon_inverse(route(n, top)),
    shift = J::array();
  for (auto &x : complement) {
    I g = as_i(x);
    shift.push_back((g > 0 ? 1 : -1) * checked_add(magnitude(g), size));
  }
  J rebuilt = join(join(join(tr, component), shift), br);
  return J{{"n", n},
           {"word", word},
           {"bottom_ports", selected},
           {"top_ports", top},
           {"component_word", component},
           {"complement_word", complement},
           {"bottom_route", br},
           {"top_route", tr},
           {"reconstructed_word", rebuilt}};
}
void check_split(const J &w) {
  if (split_layered(w.at("word"), as_i(w.at("n")), w.at("bottom_ports")) != w)
    throw RibbonError(
        "layered component split witness does not match its source diagram");
}
J extract_slider(const J &w, const J &ul, const J &ll, const J &up,
                 const J &lp) {
  size_t at = 0, count = 0;
  for (size_t i = 0; i < w.size(); i++)
    if (!is_int(w[i])) {
      at = i;
      count++;
    }
  if (count != 1)
    throw RibbonError(
        "an isolated component with exactly one explicit slider is required");
  J t = w[at];
  if (!slider(t) || as_i(t[1]) < 1 || as_i(t[2]) != 1)
    throw RibbonError(
        "a full-rank slider anchored at the first component port is required");
  J u = recognize_nested(part(w, 0, at), ul, up),
    l = recognize_nested(part(w, at + 1, w.size()), ll, lp);
  if (u["labels"].size() != static_cast<size_t>(as_i(t[1])) ||
      l["labels"].size() != static_cast<size_t>(as_i(t[1])))
    throw RibbonError("both marked trees must have the slider's rank");
  J mid = J::array();
  mid.push_back(t);
  return J{{"word", w},
           {"slider", t},
           {"upper", u},
           {"lower", l},
           {"reconstructed_word",
            join(join(u["reconstructed_word"], mid), l["reconstructed_word"])}};
}
void check_extraction(const J &w) {
  check_nested(w.at("upper"));
  check_nested(w.at("lower"));
  auto u = w.at("upper"), l = w.at("lower");
  if (extract_slider(w.at("word"), u.at("labels"), l.at("labels"),
                     u.at("initial_ports"), l.at("initial_ports")) != w)
    throw RibbonError(
        "extraction witness does not reconstruct the complete component");
}
J normalize_canonical_slider(const J &w, const J &top_routes,
                             const J &bottom_routes, I images, I operations) {
  size_t at = 0, count = 0;
  for (size_t i = 0; i < w.size(); i++)
    if (!is_int(w[i])) {
      at = i;
      count++;
    }
  if (count != 1)
    throw RibbonError(
        "one full-rank standard slider is required for twist collection");
  J t = w[at];
  if (!slider(t) || as_i(t[1]) < 1 || as_i(t[2]) != 1 || t[4].get<bool>() ||
      t[5].get<bool>())
    throw RibbonError(
        "normal or mirrored full-rank slider at the first port required");
  I n = as_i(t[1]);
  J it = boundary(top_routes, n), ib = boundary(bottom_routes, n),
    tb = flatten(it), bb = flatten(ib);
  if (bb.size() > w.size())
    throw RibbonError("inherited boundary routes must be the exact exterior "
                      "prefix and suffix");
  size_t end = w.size() - bb.size();
  if (tb.size() > at || end <= at || part(w, 0, tb.size()) != tb ||
      part(w, end, w.size()) != bb)
    throw RibbonError("inherited boundary routes must be the exact exterior "
                      "prefix and suffix");
  J ports = sequence(n), output = rev(ports),
    uw = braid(part(w, tb.size(), at), n),
    lw = ribbon_inverse(braid(part(w, at + 1, end), n));
  if (t[3].get<bool>()) {
    J half = J::array();
    for (I width = n - 1; width > 0; width--)
      for (I x = 1; x <= width; x++)
        half.push_back(x);
    uw = join(uw, half);
    lw = join(lw, half);
  }
  J u = recognize_nested(uw, output, ports, images, operations),
    l = recognize_nested(lw, ports, ports, images, operations), tf = J::array(),
    bf = J::array(), tw = J::array(), cc = J::array();
  for (I k = 0; k < n - 1; k++) {
    I a = as_i(u["crossing_counts"][k]), b = as_i(l["crossing_counts"][k]);
    tf.push_back(k ? ((a % 2 + 2) % 2) : 0);
    bf.push_back(k ? ((b % 2 + 2) % 2) : 0);
    I power = checked_sub(k ? floorhalf(a) : a, k ? floorhalf(b) : b);
    tw.push_back(power);
    cc.push_back(k ? checked_mul(2, power) : power);
  }
  J top = join(tb, reconstruct_nested(output, tf, ports)),
    bottom = join(ribbon_inverse(reconstruct_nested(ports, bf, ports)), bb),
    normal = J::array({"slider", n, 1, false, false, false}),
    body = reconstruct_nested(output, cc, ports);
  body.push_back(normal);
  return J{{"word", w},
           {"slider", normal},
           {"upper", u},
           {"lower_inverse", l},
           {"twist_powers", tw},
           {"top_route", top},
           {"bottom_route", bottom},
           {"normalized_word", join(join(top, body), bottom)},
           {"inherited_top", it},
           {"inherited_bottom", ib},
           {"max_image_letters", images},
           {"max_operations", operations}};
}
void check_canonical_slider(const J &w) {
  check_nested(w.at("upper"));
  check_nested(w.at("lower_inverse"));
  if (normalize_canonical_slider(
          w.at("word"), w.at("inherited_top"), w.at("inherited_bottom"),
          as_i(w.at("max_image_letters")), as_i(w.at("max_operations"))) != w)
    throw RibbonError(
        "canonical slider witness does not reconstruct its source diagram");
}
} // namespace zipper_verify
