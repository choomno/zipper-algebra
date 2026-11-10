// Independent typed ribbon recognition and reconstruction.
#include "typed_ribbons.hpp"
#include <algorithm>
#include <set>
namespace zipper_verify::typed {
namespace {
using V = std::vector<I>;
std::pair<V, V> ports(const V &labels, const V &initial) {
  V t = labels, p = initial;
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
V braid(const V &w, I n) {
  V a = w;
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
std::vector<V> action(const V &w, I n, I m, I maxop, I &op) {
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
Gathering gathering(const V &w, const V &c, I n, I m, I maxop) {
  limits(m, maxop);
  I op = 0;
  auto a = action(w, n, m, maxop, op), b = action(c, n, m, maxop, op);
  if (a != b)
    throw RibbonError(
        "pair-count candidate is not equal to the original braid");
  return {m,maxop,sha256(J(a).dump())};
}
RibbonWitness witness(const V &t,const V &boundary,const V &w,const std::vector<RibbonBlock> &blocks,
                      const V &counts,const V &p,const V &canon,const std::optional<Gathering> &g=std::nullopt) {
 return {t,boundary,w,counts,p,canon,blocks,g};
}
} // namespace
Labels reconstruct_nested(const Labels &labels,const Labels &counts,const Labels &initial) {
  auto [t, p] = ports(labels, initial);
  V h = counts, c;
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
RibbonWitness recognize_nested_strict(const Labels &word,const Labels &labels,const Labels &initial) {
  auto [t, b] = ports(labels, initial);
  V w = braid(word, t.size()), p = b, h(t.size() - 1);
  std::map<I, I> order;
  for (size_t i = 0; i < t.size(); ++i)
    order[t[i]] = i;
  std::vector<RibbonBlock> blocks;
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
    blocks.push_back({n,sign,block});
    off += n - 1;
  }
  return witness(t, b, w, blocks, h, p, reconstruct_nested(t, h, b));
}
RibbonWitness recognize_nested(const Labels &word,const Labels &labels,const Labels &initial,I m,I o) {
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
  V canon = reconstruct_nested(t, h, b);
  Gathering g = gathering(w, canon, t.size(), m, o);
  RibbonWitness target = recognize_nested_strict(canon, t, b);
  return witness(t, b, w, target.blocks, h, p, canon, g);
}
J encode_ribbon(const RibbonWitness &w) {
 J blocks=J::array();
 for(const auto &b:w.blocks) blocks.push_back({{"node_size",b.node_size},{"sign",b.sign},{"chronological_word",b.chronological_word}});
 J gathering=nullptr;
 if(w.gathering) gathering={{"max_image_letters",w.gathering->max_image_letters},{"max_operations",w.gathering->max_operations},{"action_sha256",w.gathering->action_sha256}};
 return {{"labels",w.labels},{"initial_ports",w.initial_ports},{"word",w.word},{"blocks",std::move(blocks)},
         {"crossing_counts",w.crossing_counts},{"final_ports",w.final_ports},{"reconstructed_word",w.reconstructed_word},{"gathering",std::move(gathering)}};
}
} // namespace zipper_verify::typed
