#include "ribbons.hpp"
#include <algorithm>
#include <set>
namespace zipper_verify {
namespace {
using V = std::vector<I>;
V ints(const J &a) {
  V v;
  if (!a.is_array())
    throw RibbonError("integer sequence required");
  for (auto &x : a) {
    if (!is_int(x))
      throw RibbonError("integer sequence required");
    v.push_back(as_i(x));
  }
  return v;
}
std::pair<V, V> ports(const J &labels, const J &initial) {
  V t = ints(labels), p = initial.is_null() ? t : ints(initial);
  std::set<I> s(t.begin(), t.end()), u(p.begin(), p.end());
  if (t.empty() || s.size() != t.size() || p.size() != t.size() || s != u)
    throw RibbonError(
        "distinct integer tree labels and the same boundary labels required");
  for (size_t n = 2; n <= t.size(); ++n) {
    V pos;
    for (size_t k = 0; k < n; ++k)
      pos.push_back(std::find(p.begin(), p.end(), t[k]) - p.begin());
    auto mm = std::minmax_element(pos.begin(), pos.end());
    if (*mm.second - *mm.first + 1 != I(n))
      throw RibbonError(
          "every marked subtree must occupy a contiguous boundary interval");
  }
  return {t, p};
}
V braid(const J &w, I n) {
  V a = ints(w);
  for (I x : a)
    if (abs_i(x) < 1 || abs_i(x) >= n)
      throw RibbonError(
          "signed braid generators inside the component required");
  return a;
}
std::pair<I, bool> children(const V &p, const V &t, I n) {
  V pos;
  for (I k = 0; k < n - 1; ++k)
    pos.push_back(std::find(p.begin(), p.end(), t[k]) - p.begin());
  auto mm = std::minmax_element(pos.begin(), pos.end());
  if (*mm.second - *mm.first + 1 != n - 1)
    throw RibbonError("the marked child ribbon is not contiguous");
  I a = std::find(p.begin(), p.end(), t[n - 1]) - p.begin();
  if (a == *mm.second + 1)
    return {*mm.first, true};
  if (a == *mm.first - 1)
    return {a, false};
  throw RibbonError("the marked children are not adjacent");
}
void limits(I m, I o) {
  if (m < 1 || o < 1)
    throw RibbonError("positive integer ribbon gathering limits are required");
}
J action(const V &w, I n, I m, I maxop, I &op) {
  if (n > m)
    throw RibbonLimit("ribbon gathering Artin image limit exceeded");
  std::vector<V> images;
  for (I i = 1; i <= n; ++i)
    images.push_back({i});
  I total = n;
  auto tick = [&]() {
    if (op == maxop)
      throw RibbonLimit("ribbon gathering operation limit exceeded");
    ++op;
  };
  for (I g : w) {
    tick();
    I at = abs_i(g) - 1;
    const V a = images[at], b = images[at + 1];
    V stack;
    auto put = [&](I x) {
      tick();
      if (!stack.empty() && stack.back() == neg(x))
        stack.pop_back();
      else
        stack.push_back(x);
      if (stack.size() > size_t(m))
        throw RibbonLimit("ribbon gathering Artin image limit exceeded");
    };
    if (g > 0) {
      for (I x : a)
        put(x);
      for (I x : b)
        put(x);
      for (auto it = a.rbegin(); it != a.rend(); ++it)
        put(neg(*it));
    } else {
      for (auto it = b.rbegin(); it != b.rend(); ++it)
        put(neg(*it));
      for (I x : a)
        put(x);
      for (I x : b)
        put(x);
    }
    total = checked_add(
        total, checked_sub(I(stack.size()), I(g > 0 ? b.size() : a.size())));
    if (total > m)
      throw RibbonLimit("ribbon gathering Artin image limit exceeded");
    if (g > 0) {
      images[at] = std::move(stack);
      images[at + 1] = a;
    } else {
      images[at] = b;
      images[at + 1] = std::move(stack);
    }
  }
  return images;
}
J gathering(const V &w, const V &c, I n, I m, I maxop) {
  limits(m, maxop);
  I op = 0;
  J a = action(w, n, m, maxop, op), b = action(c, n, m, maxop, op);
  if (a != b)
    throw RibbonError(
        "pair-count candidate is not equal to the original braid");
  return J{{"max_image_letters", m},
           {"max_operations", maxop},
           {"action_sha256", sha256(a.dump())}};
}
J witness(const V &t, const V &boundary, const V &w, const J &blocks,
          const V &counts, const V &p, const J &canon, const J &g = nullptr) {
  return J{{"labels", t},
           {"initial_ports", boundary},
           {"word", w},
           {"blocks", blocks},
           {"crossing_counts", counts},
           {"final_ports", p},
           {"reconstructed_word", canon},
           {"gathering", g}};
}
} // namespace
J reconstruct_nested(const J &labels, const J &counts, const J &initial) {
  auto [t, p] = ports(labels, initial);
  V h = ints(counts), c;
  if (h.size() + 1 != t.size())
    throw RibbonError("one integer crossing count per non-leaf node required");
  for (I n = 2; n <= I(t.size()); ++n) {
    I sign = h[n - 2] > 0 ? 1 : -1;
    for (I k = 0; k < abs_i(h[n - 2]); ++k) {
      auto [start, left] = children(p, t, n);
      for (I j = 0; j < n - 1; ++j) {
        I g = left ? start + n - 1 - j : start + 1 + j;
        c.push_back(sign * g);
        std::swap(p[g - 1], p[g]);
      }
    }
  }
  std::reverse(c.begin(), c.end());
  return c;
}
J recognize_nested_strict(const J &word, const J &labels, const J &initial) {
  auto [t, b] = ports(labels, initial);
  V w = braid(word, t.size()), p = b, h(t.size() - 1);
  std::map<I, I> order;
  for (size_t i = 0; i < t.size(); ++i)
    order[t[i]] = i;
  J blocks = J::array();
  V c = w;
  std::reverse(c.begin(), c.end());
  size_t off = 0;
  while (off < c.size()) {
    I at = abs_i(c[off]) - 1, n = std::max(order[p[at]], order[p[at + 1]]) + 1;
    children(p, t, n);
    I sign = c[off] > 0 ? 1 : -1;
    if (c.size() - off < size_t(n - 1))
      throw RibbonError("incomplete crossing of the marked child ribbon");
    std::set<I> remaining(t.begin(), t.begin() + n - 1);
    V block;
    for (I k = 0; k < n - 1; ++k) {
      I g = c[off + k], a = abs_i(g) - 1, x = p[a], y = p[a + 1],
        other = x == t[n - 1] ? y : x;
      if ((g > 0 ? 1 : -1) != sign || (x != t[n - 1] && y != t[n - 1]) ||
          !remaining.erase(other))
        throw RibbonError("a whole ribbon crossing must meet every child leaf "
                          "once with one sign");
      std::swap(p[a], p[a + 1]);
      block.push_back(g);
    }
    h[n - 2] = checked_add(h[n - 2], sign);
    blocks.push_back(
        J{{"node_size", n}, {"sign", sign}, {"chronological_word", block}});
    off += n - 1;
  }
  return witness(t, b, w, blocks, h, p, reconstruct_nested(t, h, b));
}
J recognize_nested(const J &word, const J &labels, const J &initial, I m, I o) {
  limits(m, o);
  auto [t, b] = ports(labels, initial);
  V w = braid(word, t.size());
  try {
    return recognize_nested_strict(w, t, b);
  } catch (const RibbonError &) {
  }
  std::map<I, I> order;
  for (size_t k = 0; k < t.size(); ++k)
    order[t[k]] = k;
  V p = b, h;
  std::map<std::pair<I, I>, I> pairs;
  for (auto it = w.rbegin(); it != w.rend(); ++it) {
    I at = abs_i(*it) - 1, a = order[p[at]], z = order[p[at + 1]];
    auto key = std::minmax(a, z);
    pairs[key] = checked_add(pairs[key], *it > 0 ? 1 : -1);
    std::swap(p[at], p[at + 1]);
  }
  for (I n = 1; n < I(t.size()); ++n) {
    I v = pairs[{0, n}];
    for (I k = 1; k < n; ++k)
      if (pairs[{k, n}] != v)
        throw RibbonError("leaf-pair crossing counts do not define a marked "
                          "ribbon candidate");
    h.push_back(v);
  }
  J canon = reconstruct_nested(t, h, b);
  J g = gathering(w, ints(canon), t.size(), m, o);
  J target = recognize_nested_strict(canon, t, b);
  return witness(t, b, w, target["blocks"], h, p, canon, g);
}
void check_nested(const J &w) {
  if (!w.is_object())
    throw RibbonError("a nested ribbon witness is required");
  try {
    J expected;
    if (w.at("gathering").is_null())
      expected = recognize_nested_strict(w.at("word"), w.at("labels"),
                                         w.at("initial_ports"));
    else {
      auto &g = w.at("gathering");
      expected = recognize_nested_strict(w.at("reconstructed_word"),
                                         w.at("labels"), w.at("initial_ports"));
      expected["word"] = w.at("word");
      expected["gathering"] =
          gathering(braid(w.at("word"), I(w.at("labels").size())),
                    braid(w.at("reconstructed_word"), I(w.at("labels").size())),
                    I(w.at("labels").size()), as_i(g.at("max_image_letters")),
                    as_i(g.at("max_operations")));
    }
    if (!strict_equal(expected, w))
      throw RibbonError(
          "nested ribbon witness does not match its complete braid word");
  } catch (const nlohmann::json::exception &) {
    throw RibbonError("a nested ribbon witness is required");
  }
}
} // namespace zipper_verify
