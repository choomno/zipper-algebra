#include "engine.hpp"
#include <limits>
#include <sstream>
namespace zipper {
Int checked_add(Int a, Int b) {
    Int c;
    if (__builtin_add_overflow(a, b, &c))
        throw LimitExceeded("signed integer coordinate overflow");
    return c;
}
Int checked_sub(Int a, Int b) {
    Int c;
    if (__builtin_sub_overflow(a, b, &c))
        throw LimitExceeded("signed integer coordinate overflow");
    return c;
}
Int checked_mul(Int a, Int b) {
    Int c;
    if (__builtin_mul_overflow(a, b, &c))
        throw LimitExceeded("signed integer coordinate overflow");
    return c;
}
namespace {
Int value_int(const Json &j) {
    if (!j.is_number_integer())
        throw EngineError("integer required");
    if (j.is_number_unsigned() &&
        j.get<std::uint64_t>() > static_cast<std::uint64_t>(std::numeric_limits<Int>::max()))
        throw LimitExceeded("signed integer coordinate overflow");
    return j.get<Int>();
}
BigInt big_from_json(const Json &j) {
    if (j.is_string()) {
        std::string s = j.get<std::string>();
        if (s.empty() || s == "-" || s.find_first_not_of("-0123456789") != std::string::npos ||
            s.find('-', 1) != std::string::npos)
            throw EngineError("invalid integer coefficient");
        BigInt x = 0;
        bool neg = s[0] == '-';
        for (std::size_t i = neg ? 1 : 0; i < s.size(); ++i)
            x = x * 10 + (s[i] - '0');
        return neg ? -x : x;
    }
    if (j.is_number_unsigned())
        return BigInt(j.get<std::uint64_t>());
    if (j.is_number_integer())
        return BigInt(j.get<Int>());
    throw EngineError("integer coefficient required");
}
Int magnitude(Int a) {
    if (a == std::numeric_limits<Int>::min())
        throw LimitExceeded("signed integer coordinate overflow");
    return a < 0 ? -a : a;
}
void require_fields(const Json &j, const std::vector<std::string> &required,
                    const std::vector<std::string> &optional = {}) {
    if (!j.is_object())
        throw CatalogError("AST object required");
    for (const auto &k : required)
        if (!j.contains(k))
            throw CatalogError("Missing AST field: " + k);
    for (auto it = j.begin(); it != j.end(); ++it)
        if (std::find(required.begin(), required.end(), it.key()) == required.end() &&
            std::find(optional.begin(), optional.end(), it.key()) == optional.end())
            throw CatalogError("Unexpected AST field: " + it.key());
}
Json negative(const Json &x) {
    if (x.is_number_integer())
        return checked_sub(0, value_int(x));
    return Json{{"op", "sub"}, {"args", Json::array({0, x})}};
}
Json semilinear_word(Json node) {
    std::string k = node.at("type");
    if (k == "identity" || k == "slider")
        return node;
    if (k == "sigma" || k == "fat" || k == "block") {
        node["sign"] = negative(node.value("sign", Json(1)));
        return node;
    }
    if (k == "half_twist")
        return Json{{"type", "inverse"}, {"word", node}};
    if (k == "concat") {
        for (auto &c : node["items"])
            c = semilinear_word(c);
        return node;
    }
    if (k == "inverse" || k == "bar" || k == "reverse" || k == "reflect") {
        node["word"] = semilinear_word(node["word"]);
        return node;
    }
    if (k == "word_ref")
        return Json{{"type", "bar"}, {"word", node}};
    throw CatalogError("unknown semilinear word");
}
} // namespace
Int integer(const Json &e, const Env &env) {
    if (e.is_number_integer())
        return value_int(e);
    if (e.is_string()) {
        auto it = env.find(e.get<std::string>());
        if (it == env.end())
            throw CatalogError("missing integer parameter");
        return it->second;
    }
    require_fields(e, {"op", "args"});
    std::string op = e.at("op");
    const auto &a = e.at("args");
    if (!a.is_array() || a.size() < 2)
        throw CatalogError("invalid integer operation");
    if (op == "sub") {
        if (a.size() != 2)
            throw CatalogError("sub has arity two");
        return checked_sub(integer(a[0], env), integer(a[1], env));
    }
    if (op != "add" && op != "mul")
        throw CatalogError("invalid integer operation");
    Int v = op == "add" ? 0 : 1;
    for (const auto &x : a)
        v = op == "add" ? checked_add(v, integer(x, env)) : checked_mul(v, integer(x, env));
    return v;
}
bool condition(const Json &e, const Env &env) {
    if (e.is_boolean())
        return e.get<bool>();
    if (!e.is_object())
        throw CatalogError("invalid condition");
    std::string op = e.at("op");
    if (op == "and" || op == "or") {
        require_fields(e, {"op", "args"});
        if (!e.at("args").is_array())
            throw CatalogError("condition arguments must be a list");
        bool out = op == "and";
        for (const auto &x : e.at("args")) {
            bool v = condition(x, env);
            out = op == "and" ? out && v : out || v;
        }
        return out;
    }
    if (op == "not") {
        require_fields(e, {"op", "arg"});
        return !condition(e.at("arg"), env);
    }
    require_fields(e, {"op", "left", "right"});
    Int a = integer(e.at("left"), env), b = integer(e.at("right"), env);
    if (op == "eq")
        return a == b;
    if (op == "ne")
        return a != b;
    if (op == "lt")
        return a < b;
    if (op == "le")
        return a <= b;
    if (op == "gt")
        return a > b;
    if (op == "ge")
        return a >= b;
    throw CatalogError("invalid condition operation");
}
Word expand_word(const Json &node, const Env &env, Int budget, const NamedWords &words) {
    if (budget < 1)
        throw CatalogError("Expansion budget must be positive");
    if (!node.is_object())
        throw CatalogError("Word must be an AST object");
    auto rank = env.find("n");
    if (rank == env.end() || rank->second < 1)
        throw CatalogError("n must be a positive integer");
    Int n = rank->second;
    std::string kind = node.at("type");
    Word result;
    auto iv = [&](const char *k) { return integer(node.at(k), env); };
    auto emit = [&](Int x) {
        if (x == 0 || magnitude(x) >= n)
            throw CatalogError("Braid generator outside B_n");
        if (result.size() >= static_cast<std::size_t>(budget))
            throw LimitExceeded("Expansion word budget exceeded");
        result.emplace_back(x);
    };
    if (kind == "identity") {
        require_fields(node, {"type"});
        return {};
    }
    if (kind == "sigma") {
        require_fields(node, {"type", "index"}, {"sign"});
        Int index = iv("index"), sign = integer(node.value("sign", Json(1)), env);
        if (index < 1 || (sign != 1 && sign != -1))
            throw CatalogError("Invalid signed generator");
        emit(checked_mul(index, sign));
    } else if (kind == "fat" || kind == "block") {
        if (kind == "fat")
            require_fields(node, {"type", "length", "start"}, {"sign"});
        else
            require_fields(node, {"type", "length", "start", "count", "order"}, {"sign"});
        Int length = iv("length"), start = iv("start"), count = kind == "block" ? iv("count") : 1,
            sign = integer(node.value("sign", Json(1)), env);
        std::string order = node.value("order", std::string("ascending"));
        if (length < 0 || count < 0 || start < 1 || (sign != 1 && sign != -1))
            throw CatalogError("Invalid block parameters");
        if (order != "ascending" && order != "descending")
            throw CatalogError("Invalid block order");
        if (length == 0 || count == 0)
            return {};
        if (checked_mul(length, count) > budget)
            throw LimitExceeded("Expansion word budget exceeded");
        for (Int t = 0; t < count; ++t) {
            Int offset = order == "ascending" ? t : count - t - 1;
            for (Int x = 0; x < length; ++x)
                emit(checked_mul(sign, checked_add(checked_add(start, offset), x)));
        }
    } else if (kind == "half_twist") {
        require_fields(node, {"type", "size", "start"});
        Int size = iv("size"), start = iv("start");
        if (size < 1 || start < 1 || checked_add(start, size - 1) > n)
            throw CatalogError("Invalid half-twist interval");
        if (checked_mul(size, size - 1) / 2 > budget)
            throw LimitExceeded("Expansion word budget exceeded");
        for (Int length = size - 1; length > 0; --length)
            for (Int x = 0; x < length; ++x)
                emit(checked_add(start, x));
    } else if (kind == "slider") {
        require_fields(node, {"type", "size", "start"}, {"orientation"});
        Int size = iv("size"), start = iv("start");
        std::string orientation = node.value("orientation", std::string("normal"));
        if (size < 1 || start < 1 || checked_add(start, size - 1) > n ||
            (orientation != "normal" && orientation != "mirrored"))
            throw CatalogError("Invalid slider");
        result.push_back(Token::slider(size, start, orientation == "mirrored"));
    } else if (kind == "concat") {
        require_fields(node, {"type", "items"});
        if (!node.at("items").is_array())
            throw CatalogError("concat.items must be a list");
        for (const auto &c : node.at("items")) {
            Word w = expand_word(c, env, budget, words);
            if (result.size() + w.size() > static_cast<std::size_t>(budget))
                throw LimitExceeded("Expansion word budget exceeded");
            append(result, w);
        }
    } else if (kind == "word_ref") {
        require_fields(node, {"type", "name"});
        auto it = words.find(node.at("name").get<std::string>());
        if (it == words.end())
            throw CatalogError("Unbound word reference");
        for (Int x : it->second)
            emit(x);
    } else if (kind == "inverse" || kind == "bar" || kind == "reverse" || kind == "reflect") {
        require_fields(node, {"type", "word"},
                       kind == "reflect" ? std::vector<std::string>{"n"}
                                         : std::vector<std::string>{});
        result = expand_word(node.at("word"), env, budget, words);
        if (kind == "reflect" && integer(node.value("n", Json(n)), env) != n)
            throw CatalogError("reflect.n must be ambient n");
        if (kind == "inverse" || kind == "reverse")
            std::reverse(result.begin(), result.end());
        for (auto &t : result) {
            if (!t.is_slider()) {
                if (kind == "inverse" || kind == "bar")
                    t.generator = checked_sub(0, t.generator);
                else if (kind == "reflect")
                    t.generator = checked_mul(t.generator > 0 ? 1 : -1, n - magnitude(t.generator));
            } else {
                if (kind == "inverse")
                    throw CatalogError("inverse is only defined for braid words");
                if (kind == "bar")
                    t.barred = !t.barred;
                else if (kind == "reverse")
                    t.reversed = !t.reversed;
                else {
                    t.start = checked_add(checked_sub(checked_sub(n, t.start), t.size), 2);
                    t.mirrored = !t.mirrored;
                }
            }
        }
    } else
        throw CatalogError("Unknown word node: " + kind);
    return result;
}
Polynomial expand_polynomial(const Json &terms, const Env &env, Int budget,
                             const NamedWords &words) {
    if (!terms.is_array())
        throw CatalogError("Polynomial must be a list");
    Polynomial out;
    std::map<std::pair<Int, Word>, std::size_t> indices;
    for (const auto &term : terms) {
        require_fields(term, {"coefficient", "word"});
        Word word = expand_word(term.at("word"), env, budget, words);
        if (!term.at("coefficient").is_array())
            throw CatalogError("Coefficient must be Laurent term list");
        for (const auto &coeff : term.at("coefficient")) {
            require_fields(coeff, {"integer", "q_power"});
            Int q = integer(coeff.at("q_power"), env);
            BigInt c = big_from_json(coeff.at("integer"));
            auto key = std::make_pair(q, word);
            auto it = indices.find(key);
            if (it == indices.end()) {
                indices.emplace(std::move(key), out.size());
                out.push_back({q, word, c});
            } else
                out[it->second].coefficient += c;
        }
        if (out.size() > static_cast<std::size_t>(budget))
            throw LimitExceeded("Expansion term budget exceeded");
    }
    out.erase(std::remove_if(out.begin(), out.end(),
                             [](const PolynomialTerm &t) { return t.coefficient == 0; }),
              out.end());
    return out;
}
namespace {
std::string python_word_repr(const Word &word) {
    std::string result = "(";
    for (std::size_t i = 0; i < word.size(); ++i) {
        if (i)
            result += ", ";
        const auto &t = word[i];
        if (!t.is_slider())
            result += std::to_string(t.generator);
        else
            result += "('slider', " + std::to_string(t.size) + ", " + std::to_string(t.start) +
                      ", " + (t.mirrored ? "True" : "False") + ", " +
                      (t.reversed ? "True" : "False") + ", " + (t.barred ? "True" : "False") + ")";
    }
    if (word.size() == 1)
        result += ",";
    return result + ")";
}
} // namespace
Polynomial serialize_order(Polynomial polynomial) {
    std::vector<std::pair<std::string, PolynomialTerm>> keyed;
    for (auto &t : polynomial) {
        std::string key = "(" + std::to_string(t.q) + ", " + python_word_repr(t.word) + ")";
        keyed.emplace_back(std::move(key), std::move(t));
    }
    std::stable_sort(keyed.begin(), keyed.end(),
                     [](const auto &a, const auto &b) { return a.first < b.first; });
    polynomial.clear();
    for (auto &item : keyed)
        polynomial.push_back(std::move(item.second));
    return polynomial;
}
std::vector<Json> flatten(const Json &node) {
    if (node.at("type") != "concat")
        return {node};
    std::vector<Json> result;
    for (const auto &c : node.at("items")) {
        auto a = flatten(c);
        append(result, a);
    }
    return result;
}
Word literal(const std::vector<Json> &nodes, const Env &env, Int budget) {
    return expand_word(Json{{"type", "concat"}, {"items", nodes}}, env, budget);
}
Catalog::Catalog(const Json &source) {
    if (!source.is_object() || !source.contains("rules") || !source.at("rules").is_array())
        throw CatalogError("catalog rules required");
    auto add = [&](const Json &j) {
        std::string id = j.at("id");
        if (rules_.count(id))
            throw CatalogError("duplicate catalog id");
        if (!j.contains("execution") || !j.at("execution").value("algebra_verified", false))
            throw CatalogError("rule lacks algebraic verification");
        rules_.emplace(id, Rule{id, j.value("base_id", std::string()), j.at("lhs"), j.at("rhs"),
                                j.at("condition"), j});
        ordered_ids_.push_back(id);
    };
    for (const auto &j : source.at("rules"))
        add(j);
    std::vector<Json> derived;
    for (const auto &id : ordered_ids_) {
        if (id.rfind("zip.", 0) != 0 || id.size() < 2 || id.substr(id.size() - 2) != ".M")
            continue;
        Json j = rules_.at(id).source;
        std::string source_hash = sha256_json(j);
        j["id"] = id.substr(0, id.size() - 2) + ".BM";
        j["source_record_id"] = id;
        j["source_rule_sha256"] = source_hash;
        j["semilinear"] =
            Json{{"proof_version", "zipper-semilinear-bar-v1"},
                 {"source_record_id", id},
                 {"source_rule_sha256", source_hash},
                 {"coefficient_map", "q -> q^-1"},
                 {"crossing_map", "sigma_l -> sigma_l^-1"},
                 {"preserves_product_order", true},
                 {"half_twist_identity", "bar(Delta) = inverse(Delta)"},
                 {"author_review",
                  "inherits partial source review; semilinear derivation recorded separately"}};
        j["symmetry"].update(Json{{"name", "BM"},
                                  {"invert_q", true},
                                  {"invert_crossing_sign", true},
                                  {"semilinear", true}});
        j["execution"]["derived_semilinear"] = true;
        for (const char *key :
             {"original_relator_cofactor", "cofactor_original_notation", "cofactor_transport"})
            if (j.contains(key)) {
                j[std::string("semilinear_source_") + key] = j.at(key);
                j.erase(key);
            }
        for (const char *side : {"lhs", "rhs"})
            for (auto &term : j[side]) {
                term["word"] = semilinear_word(term["word"]);
                for (auto &coeff : term["coefficient"])
                    coeff["q_power"] = negative(coeff["q_power"]);
            }
        derived.push_back(std::move(j));
    }
    for (const auto &j : derived)
        add(j);
}
const Rule &Catalog::get(const std::string &id) const {
    auto it = rules_.find(id);
    if (it == rules_.end())
        throw CatalogError("unknown rule: " + id);
    return it->second;
}
std::vector<const Rule *> Catalog::select(const Env &env, const std::string &prefix,
                                          const std::string &suffix, bool exclude_flype) const {
    std::vector<const Rule *> result;
    for (const auto &id : ordered_ids_) {
        if (id.rfind(prefix, 0) != 0 || id.size() < suffix.size() ||
            id.compare(id.size() - suffix.size(), suffix.size(), suffix) != 0 ||
            (exclude_flype && id.rfind("zip.flype.", 0) == 0))
            continue;
        const auto &r = rules_.at(id);
        if (condition(r.guard, env))
            result.push_back(&r);
    }
    return result;
}
Laurent multiply(const Laurent &a, const Laurent &b) {
    Laurent r;
    for (const auto &[q, c] : a)
        for (const auto &[p, d] : b)
            r[checked_add(q, p)] += c * d;
    for (auto it = r.begin(); it != r.end();)
        if (it->second == 0)
            it = r.erase(it);
        else
            ++it;
    return r;
}
void add_to(Laurent &a, const Laurent &b) {
    for (const auto &[q, c] : b) {
        auto it = a.find(q);
        if (it == a.end()) {
            if (c != 0)
                a.emplace(q, c);
        } else {
            it->second += c;
            if (it->second == 0)
                a.erase(it);
        }
    }
}
std::size_t integer_bits(const BigInt &c) {
    if (c == 0)
        return 0;
    BigInt positive = c < 0 ? -c : c;
    return boost::multiprecision::msb(positive) + 1;
}
void check_branch(const Branch &b, const Limits &l) {
    if (b.word.size() > static_cast<std::size_t>(l.max_word_size))
        throw LimitExceeded("word size limit exceeded");
    if (b.coefficient.size() > static_cast<std::size_t>(l.max_coefficient_support))
        throw LimitExceeded("coefficient support limit exceeded");
    for (const auto &[q, c] : b.coefficient) {
        (void)q;
        if (integer_bits(c) > static_cast<std::size_t>(l.max_integer_bits))
            throw LimitExceeded("coefficient integer bit limit exceeded");
    }
}
void validate_word(const Word &w, Int n) {
    if (n < 1)
        throw EngineError("positive rank required");
    for (const auto &t : w) {
        if (t.is_slider()) {
            if (t.generator || t.size < 1 || t.size > n || t.start < 1 || t.start > n - t.size + 1)
                throw EngineError("invalid slider token");
        } else if (!t.generator || magnitude(t.generator) >= n)
            throw EngineError("invalid generator token");
    }
}
Word word_from_braid(const Braid &b) {
    return Word(b.begin(), b.end());
}
Braid braid_from_word(const Word &w) {
    Braid out;
    out.reserve(w.size());
    for (const auto &t : w) {
        if (t.is_slider())
            throw ControlError("pure braid word required");
        out.push_back(t.generator);
    }
    return out;
}
Word slice(const Word &w, std::size_t a, std::size_t b) {
    if (a > b || b > w.size())
        throw ControlError("word boundary outside expression");
    return Word(w.begin() + a, w.begin() + b);
}
Braid slice(const Braid &w, std::size_t a, std::size_t b) {
    if (a > b || b > w.size())
        throw ControlError("word boundary outside expression");
    return Braid(w.begin() + a, w.begin() + b);
}
std::vector<Branch> apply_rule(const Branch &old, const Rule &rule, const Env &env, Int offset,
                               const Limits &limits, const std::vector<Control> *controls,
                               const NamedWords &words) {
    check_branch(old, limits);
    if (offset < 0)
        throw ControlError("negative occurrence offset");
    if (!condition(rule.guard, env))
        throw CatalogError("rule guard is false");
    Polynomial lhs = expand_polynomial(rule.lhs, env, limits.expansion_budget, words),
               rhs = expand_polynomial(rule.rhs, env, limits.expansion_budget, words);
    if (lhs.size() != 1 || lhs[0].q != 0 || lhs[0].coefficient != 1 || lhs[0].word.empty())
        throw CatalogError("nonempty monomial left side required");
    const auto &occurrence = lhs[0].word;
    std::size_t start = static_cast<std::size_t>(offset);
    if (start > old.word.size() || occurrence.size() > old.word.size() - start ||
        !std::equal(occurrence.begin(), occurrence.end(), old.word.begin() + start))
        throw ControlError("literal left side does not match selected occurrence");
    std::map<Word, std::size_t> index;
    std::vector<std::pair<Word, Laurent>> grouped;
    for (const auto &t : rhs) {
        auto it = index.find(t.word);
        if (it == index.end()) {
            index.emplace(t.word, grouped.size());
            grouped.emplace_back(t.word, Laurent{});
            it = index.find(t.word);
        }
        grouped[it->second].second[t.q] += t.coefficient;
    }
    for (auto &[word, coefficient] : grouped) {
        (void)word;
        for (auto it = coefficient.begin(); it != coefficient.end();)
            if (it->second == 0)
                it = coefficient.erase(it);
            else
                ++it;
    }
    grouped.erase(std::remove_if(grouped.begin(), grouped.end(),
                                 [](const auto &t) { return t.second.empty(); }),
                  grouped.end());
    if (controls && controls->size() != grouped.size())
        throw ControlError("one control update required per RHS word");
    if (grouped.size() > static_cast<std::size_t>(limits.max_terms))
        throw LimitExceeded("term limit exceeded");
    std::vector<Branch> result;
    for (std::size_t i = 0; i < grouped.size(); ++i) {
        Word word = slice(old.word, 0, start);
        append(word, grouped[i].first);
        append(word, slice(old.word, start + occurrence.size(), old.word.size()));
        Branch b{std::move(word), multiply(old.coefficient, grouped[i].second),
                 controls ? (*controls)[i] : old.control};
        check_branch(b, limits);
        result.push_back(std::move(b));
    }
    return result;
}
Json word_json(const Word &word) {
    Json out = Json::array();
    for (const auto &t : word) {
        if (t.is_slider())
            out.push_back(
                Json::array({"slider", t.size, t.start, t.mirrored, t.reversed, t.barred}));
        else
            out.push_back(t.generator);
    }
    return out;
}
Word word_from_json(const Json &j) {
    if (!j.is_array())
        throw EngineError("word array required");
    Word word;
    for (const auto &t : j) {
        if (t.is_array()) {
            if (t.size() != 6 || t[0] != "slider" || !t[3].is_boolean() || !t[4].is_boolean() ||
                !t[5].is_boolean())
                throw EngineError("invalid slider token");
            word.push_back(Token::slider(value_int(t[1]), value_int(t[2]), t[3].get<bool>(),
                                         t[4].get<bool>(), t[5].get<bool>()));
        } else
            word.emplace_back(value_int(t));
    }
    return word;
}
Json coefficient_json(const Laurent &coefficient) {
    Json j = Json::array();
    for (const auto &[q, c] : coefficient)
        j.push_back(Json::array({q, c.convert_to<std::string>()}));
    return j;
}
Laurent coefficient_from_json(const Json &j) {
    if (!j.is_array())
        throw EngineError("coefficient array required");
    Laurent out;
    Int previous = 0;
    bool first = true;
    for (const auto &t : j) {
        if (!t.is_array() || t.size() != 2)
            throw EngineError("coefficient pair required");
        Int q = value_int(t[0]);
        BigInt c = big_from_json(t[1]);
        if (c == 0 || (!first && q <= previous))
            throw EngineError("coefficient exponents must increase and values be nonzero");
        out.emplace(q, c);
        previous = q;
        first = false;
    }
    return out;
}
Json event_json(const ProcessedEvent &e) {
    return Json{{"rule_id", e.rule_id},
                {"parameters", e.parameters},
                {"rhs_index", e.rhs_index},
                {"mirror_policy", e.mirror_policy}};
}
Json control_json(const Control &c) {
    if (c.phase == "unsupported_frame")
        return Json{{"phase", c.phase}, {"rule_id", c.rule_id}, {"message", c.message}};
    Json j{{"phase", c.phase}, {"i", c.i}, {"j", c.j}};
    if (c.phase == "turnaround_complete") {
        j.update(Json{{"mirrored", c.mirrored},
                      {"mirror_policy", c.mirror_policy},
                      {"output_power", c.output_power},
                      {"input_power", c.input_power},
                      {"exterior_length", c.exterior_length},
                      {"payload_length", c.payload_length},
                      {"routing_length", c.routing_length},
                      {"transported_length", c.transported_length},
                      {"processed_length", c.processed_length},
                      {"tail_length", c.tail_length}});
        return j;
    }
    j.update(Json{{"k", c.k},
                  {"left_length", c.left_length},
                  {"processed_length", c.processed_length},
                  {"tail_length", c.tail_length}});
    if (c.phase == "framed_zip_up" || c.phase == "framed_upper_pass_complete") {
        const auto &a = c.initial_slider;
        Json events = Json::array();
        for (const auto &e : c.processed_events)
            events.push_back(event_json(e));
        j.update(Json{{"mirrored", c.mirrored},
                      {"output_power", c.output_power},
                      {"input_power", c.input_power},
                      {"output_payload", c.output_payload},
                      {"mirror_policy", c.mirror_policy},
                      {"initial_slider", Json{{"j", a.j},
                                              {"k", a.k},
                                              {"mirrored", a.mirrored},
                                              {"output_power", a.output_power},
                                              {"input_power", a.input_power},
                                              {"processed_length", a.processed_length},
                                              {"mirror_policy", a.mirror_policy}}},
                      {"processed_events", events}});
    }
    return j;
}
Control control_from_json(const Json &j) {
    if (!j.is_object() || !j.contains("phase") || !j.at("phase").is_string())
        throw EngineError("control object with phase required");
    const std::string phase = j.at("phase");
    if (phase == "unsupported_frame")
        require_fields(j, {"phase", "rule_id", "message"});
    else if (phase == "turnaround_complete")
        require_fields(j, {"phase", "i", "j", "mirrored", "mirror_policy", "output_power",
                           "input_power", "exterior_length", "payload_length", "routing_length",
                           "transported_length", "processed_length", "tail_length"});
    else if (phase == "framed_zip_up" || phase == "framed_upper_pass_complete")
        require_fields(j, {"phase", "i", "j", "k", "mirrored", "output_power", "input_power",
                           "left_length", "processed_length", "tail_length", "initial_slider",
                           "processed_events", "output_payload", "mirror_policy"});
    else if (phase == "zip_up" || phase == "upper_pass_complete")
        require_fields(j,
                       {"phase", "i", "j", "k", "left_length", "processed_length", "tail_length"});
    else
        throw EngineError("unknown control phase");
    Control c;
    c.phase = j.at("phase").get<std::string>();
    if (c.phase == "unsupported_frame") {
        c.rule_id = j.at("rule_id");
        c.message = j.at("message");
        return c;
    }
    c.i = value_int(j.at("i"));
    c.j = value_int(j.at("j"));
    auto get = [&](const char *key, Int def) {
        return j.contains(key) ? value_int(j.at(key)) : def;
    };
    c.k = get("k", 0);
    c.left_length = get("left_length", 0);
    c.processed_length = get("processed_length", 0);
    c.tail_length = get("tail_length", 0);
    c.mirrored = j.value("mirrored", false);
    c.output_power = get("output_power", 0);
    c.input_power = get("input_power", 0);
    if (j.contains("output_payload")) {
        if (!j.at("output_payload").is_array())
            throw EngineError("output payload array required");
        for (const auto &x : j.at("output_payload"))
            c.output_payload.push_back(value_int(x));
    }
    c.mirror_policy = j.value("mirror_policy", std::string("M"));
    c.exterior_length = get("exterior_length", 0);
    c.payload_length = get("payload_length", 0);
    c.routing_length = get("routing_length", 0);
    c.transported_length = get("transported_length", 0);
    if (j.contains("initial_slider")) {
        const auto &a = j.at("initial_slider");
        require_fields(a, {"j", "k", "mirrored", "output_power", "input_power", "processed_length",
                           "mirror_policy"});
        c.initial_slider = {value_int(a.at("j")),
                            value_int(a.at("k")),
                            value_int(a.at("output_power")),
                            value_int(a.at("input_power")),
                            value_int(a.at("processed_length")),
                            a.at("mirrored").get<bool>(),
                            a.at("mirror_policy").get<std::string>()};
    }
    if (j.contains("processed_events")) {
        if (!j.at("processed_events").is_array())
            throw EngineError("processed events array required");
        for (const auto &e : j.at("processed_events")) {
            require_fields(e, {"rule_id", "parameters", "rhs_index", "mirror_policy"});
            if (!e.at("parameters").is_object())
                throw EngineError("event parameters object required");
            ProcessedEvent ev;
            ev.rule_id = e.at("rule_id");
            for (auto it = e.at("parameters").begin(); it != e.at("parameters").end(); ++it)
                ev.parameters.emplace(it.key(), value_int(it.value()));
            ev.rhs_index = value_int(e.at("rhs_index"));
            ev.mirror_policy = e.at("mirror_policy");
            c.processed_events.push_back(std::move(ev));
        }
    }
    return c;
}
Json branch_json(const Branch &b) {
    return Json{{"word", word_json(b.word)},
                {"coefficient", coefficient_json(b.coefficient)},
                {"control", control_json(b.control)}};
}
Branch branch_from_json(const Json &j) {
    require_fields(j, {"word", "coefficient", "control"});
    if (j.at("coefficient").empty())
        throw EngineError("nonzero branch coefficient required");
    return Branch{word_from_json(j.at("word")), coefficient_from_json(j.at("coefficient")),
                  control_from_json(j.at("control"))};
}
Json limits_json(const Limits &l) {
    return Json{{"max_steps", l.max_steps},
                {"max_terms", l.max_terms},
                {"max_word_size", l.max_word_size},
                {"max_coefficient_support", l.max_coefficient_support},
                {"max_integer_bits", l.max_integer_bits},
                {"expansion_budget", l.expansion_budget}};
}
Limits limits_from_json(const Json &j) {
    if (!j.is_object())
        throw EngineError("limits object required");
    Limits l;
    for (auto it = j.begin(); it != j.end(); ++it) {
        Int v = value_int(it.value());
        if (v < 1)
            throw EngineError("positive limits required");
        if (it.key() == "max_steps")
            l.max_steps = v;
        else if (it.key() == "max_terms")
            l.max_terms = v;
        else if (it.key() == "max_word_size")
            l.max_word_size = v;
        else if (it.key() == "max_coefficient_support")
            l.max_coefficient_support = v;
        else if (it.key() == "max_integer_bits")
            l.max_integer_bits = v;
        else if (it.key() == "expansion_budget")
            l.expansion_budget = v;
        else
            throw EngineError("unknown execution limit");
    }
    return l;
}
} // namespace zipper
