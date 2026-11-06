#include "engine.hpp"
#include <algorithm>
#include <limits>
#include <sstream>

namespace zipper {
namespace {
int64_t magnitude(int64_t value) {
    if (value == std::numeric_limits<int64_t>::min())
        throw FrameError("generator index overflow");
    return value < 0 ? -value : value;
}
void validate_braid(const Braid &word, int64_t n, int64_t size, int64_t start) {
    if (n < 1 || size < 1 || size > n || start < 1 || start > n - size + 1)
        throw FrameError("valid rank and tape interval required");
    for (auto letter : word)
        if (letter == 0 || magnitude(letter) >= n)
            throw FrameError("word must contain signed braid generators in B_n");
}
bool segment(const Braid &word, int64_t offset, const Braid &block) {
    return offset >= 0 && static_cast<uint64_t>(offset) <= word.size() &&
           block.size() <= word.size() - static_cast<size_t>(offset) &&
           std::equal(block.begin(), block.end(), word.begin() + offset);
}
FrameStep local_frame(int64_t letter, int64_t size, int64_t start, int64_t power, bool left) {
    const Int index = magnitude(letter), sign = letter > 0 ? 1 : -1;
    FrameStep result{"", 0, start, start, {letter}, {letter}};
    if (index < start - 1 || (index >= start && index - start >= size)) {
        result.kind = "far";
        return result;
    }
    if (start <= index && index < start + (size - 1)) {
        result.kind = "internal";
        // Written without 2*start to avoid overflowing even for valid int64 ranks.
        const auto reflected = power % 2 ? start + (start + (size - 2) - index) : index;
        result.transformed = {sign * reflected};
        return result;
    }
    if (size == 1) {
        result.kind = "singleton";
        result.start_after += index == start ? 1 : -1;
        return result;
    }
    result.block.clear();
    if (index == (left ? start + (size - 1) : start - 1)) {
        result.start_after += left ? 1 : -1;
        const auto low = left ? start : start - 1;
        for (int64_t k = 0; k < size; ++k)
            result.block.push_back(sign * (low + k));
    } else {
        result.start_after += left ? -1 : 1;
        const auto low = left ? start - 1 : start;
        for (int64_t k = size; k > 0; --k)
            result.block.push_back(sign * (low + (k - 1)));
    }
    result.kind = result.start_after > start ? "right" : "left";
    result.transformed = result.block;
    return result;
}
Braid joined_pieces(const std::vector<Braid> &pieces, bool reverse) {
    Braid result;
    if (reverse) {
        for (auto it = pieces.rbegin(); it != pieces.rend(); ++it)
            result.insert(result.end(), it->begin(), it->end());
    } else {
        for (const auto &piece : pieces)
            result.insert(result.end(), piece.begin(), piece.end());
    }
    return result;
}
std::string transport_message(const std::string &reason, int64_t offset, int64_t start,
                              const Braid &word) {
    std::ostringstream result;
    result << reason << " at offset " << offset << ", band starts at " << start << "; remaining=(";
    for (size_t k = static_cast<size_t>(offset); k < word.size(); ++k) {
        if (k > static_cast<size_t>(offset))
            result << ", ";
        result << word[k];
    }
    if (static_cast<size_t>(offset) + 1 == word.size())
        result << ",";
    result << ")";
    return result.str();
}
} // namespace
FrameTransport transport_frame(const Braid &word, int64_t n, int64_t size, int64_t start,
                               int64_t power, const std::string &side) {
    if (side != "left" && side != "right")
        throw FrameError("valid tape transport side required");
    const bool left = side == "left";
    validate_braid(word, n, size, start);
    FrameTransport result{n, size, start, start, power, side, word, {}, {}};
    int64_t cursor = left ? static_cast<int64_t>(word.size()) : 0;
    std::vector<Braid> pieces;
    while (left ? cursor > 0 : static_cast<size_t>(cursor) < word.size()) {
        auto step = local_frame(word[static_cast<size_t>(left ? cursor - 1 : cursor)], size,
                                result.end, power, left);
        step.offset = left ? cursor - static_cast<int64_t>(step.block.size()) : cursor;
        if (!segment(word, step.offset, step.block) || step.start_after < 1 ||
            step.start_after > n - size + 1) {
            throw FrameError("whole marked tape required at offset " + std::to_string(step.offset) +
                             "; tape starts at " + std::to_string(result.end) +
                             ", side=" + (left ? "left" : "right"));
        }
        pieces.push_back(step.transformed);
        cursor = left ? step.offset : cursor + static_cast<int64_t>(step.block.size());
        result.end = step.start_after;
        result.steps.push_back(std::move(step));
    }
    result.transformed = joined_pieces(pieces, left);
    return result;
}
void check_frame(const FrameTransport &witness) {
    if (witness.side != "left" && witness.side != "right")
        throw FrameError("valid tape transport side required");
    validate_braid(witness.word, witness.n, witness.size, witness.start);
    if (witness.end < 1 || witness.end > witness.n - witness.size + 1)
        throw FrameError("invalid final tape interval");
    int64_t position = witness.start,
            cursor = (witness.side == "left") ? static_cast<int64_t>(witness.word.size()) : 0;
    std::vector<Braid> pieces;
    for (const auto &step : witness.steps) {
        if (step.block.empty() || step.start_before != position)
            throw FrameError("invalid frame step");
        for (auto letter : step.block)
            if (!letter || magnitude(letter) >= witness.n)
                throw FrameError("invalid frame step");
        for (auto letter : step.transformed)
            if (!letter || magnitude(letter) >= witness.n)
                throw FrameError("invalid frame step");
        const auto expected =
            local_frame((witness.side == "left") ? step.block.back() : step.block.front(),
                        witness.size, position, witness.power, (witness.side == "left"));
        const auto offset =
            (witness.side == "left") ? cursor - static_cast<int64_t>(step.block.size()) : cursor;
        if (step.kind != expected.kind || step.start_after != expected.start_after ||
            step.block != expected.block || step.transformed != expected.transformed ||
            step.offset != offset || !segment(witness.word, offset, step.block) ||
            step.start_after < 1 || step.start_after > witness.n - witness.size + 1)
            throw FrameError("invalid local frame identity");
        pieces.push_back(step.transformed);
        cursor =
            (witness.side == "left") ? offset : cursor + static_cast<int64_t>(step.block.size());
        position = step.start_after;
    }
    if (cursor != ((witness.side == "left") ? 0 : static_cast<int64_t>(witness.word.size())) ||
        position != witness.end ||
        joined_pieces(pieces, (witness.side == "left")) != witness.transformed)
        throw FrameError("frame witness does not cover the entire identity");
}
BandTransport transport_band(const Braid &word, int64_t n, int64_t size, int64_t start,
                             bool mirrored) {
    validate_braid(word, n, size, start);
    BandTransport result{n, size, start, start, mirrored, word, {}};
    size_t offset = 0;
    while (offset < word.size()) {
        const Int letter = word[offset], index = magnitude(letter), sign = letter > 0 ? 1 : -1;
        TransportStep step{"", static_cast<int64_t>(offset), result.end, result.end, {letter}};
        if (index < result.end - 1 || (index >= result.end && index - result.end >= size))
            step.kind = "far";
        else if (size == 1) {
            step.kind = "singleton";
            step.start_after += index == result.end ? 1 : -1;
        } else if (size == 2 && index == result.end)
            step.kind = "internal_i2";
        else if (index == result.end + (size - 1)) {
            step.kind = "right";
            ++step.start_after;
            step.block.clear();
            for (int64_t k = size; k > 0; --k)
                step.block.push_back(sign * (result.end + (k - 1)));
        } else if (index == result.end - 1) {
            step.kind = "left";
            --step.start_after;
            step.block.clear();
            for (int64_t k = 0; k < size; ++k)
                step.block.push_back(sign * (result.end - 1 + k));
        } else
            throw TransportError(
                transport_message("internal crossing requires a flype or another explicit witness",
                                  offset, result.end, word));
        if (!segment(word, offset, step.block) || step.start_after < 1 ||
            step.start_after > n - size + 1)
            throw TransportError(
                transport_message("whole-band crossing required", offset, result.end, word));
        offset += step.block.size();
        result.end = step.start_after;
        result.steps.push_back(std::move(step));
    }
    return result;
}
void check_transport(const BandTransport &witness) {
    validate_braid(witness.word, witness.n, witness.size, witness.start);
    if (witness.end < 1 || witness.end > witness.n - witness.size + 1)
        throw TransportError("invalid final band interval");
    int64_t position = witness.start, offset = 0;
    for (const auto &step : witness.steps) {
        if (step.offset != offset || step.start_before != position || step.block.empty())
            throw TransportError("invalid transport step location or block");
        for (auto letter : step.block)
            if (!letter || magnitude(letter) >= witness.n)
                throw TransportError("invalid transport step location or block");
        const Int index = magnitude(step.block.front()), sign = step.block.front() > 0 ? 1 : -1;
        Braid expected;
        bool valid = false;
        if (step.kind == "far")
            valid =
                step.block.size() == 1 && step.start_after == position &&
                (index < position - 1 || (index >= position && index - position >= witness.size));
        else if (step.kind == "internal_i2")
            valid = witness.size == 2 && step.block.size() == 1 && index == position &&
                    step.start_after == position;
        else if (step.kind == "singleton")
            valid = witness.size == 1 && step.block.size() == 1 &&
                    (index == position - 1 || index == position) &&
                    step.start_after == position + (index == position ? 1 : -1);
        else if (step.kind == "right") {
            for (int64_t k = witness.size; k > 0; --k)
                expected.push_back(sign * (position + (k - 1)));
            valid =
                step.block == expected && position < witness.n && step.start_after == position + 1;
        } else if (step.kind == "left") {
            for (int64_t k = 0; k < witness.size; ++k)
                expected.push_back(sign * (position - 1 + k));
            valid = step.block == expected && step.start_after == position - 1;
        }
        if (!valid || !segment(witness.word, offset, step.block) || step.start_after < 1 ||
            step.start_after > witness.n - witness.size + 1)
            throw TransportError("invalid local band transport identity");
        offset += static_cast<int64_t>(step.block.size());
        position = step.start_after;
    }
    if (static_cast<size_t>(offset) != witness.word.size() || position != witness.end)
        throw TransportError("transport witness does not cover the entire word");
}
Json frame_json(const FrameTransport &witness) {
    Json steps = Json::array();
    for (const auto &step : witness.steps)
        steps.push_back({{"kind", step.kind},
                         {"offset", step.offset},
                         {"start_before", step.start_before},
                         {"start_after", step.start_after},
                         {"block", step.block},
                         {"emitted", step.transformed}});
    return {{"n", witness.n},         {"size", witness.size},
            {"start", witness.start}, {"end", witness.end},
            {"power", witness.power}, {"side", witness.side},
            {"word", witness.word},   {"transformed", witness.transformed},
            {"steps", steps}};
}
Json band_json(const BandTransport &witness) {
    Json steps = Json::array();
    for (const auto &step : witness.steps)
        steps.push_back({{"kind", step.kind},
                         {"offset", step.offset},
                         {"start_before", step.start_before},
                         {"start_after", step.start_after},
                         {"block", step.block}});
    return {{"n", witness.n},     {"size", witness.size},         {"start", witness.start},
            {"end", witness.end}, {"mirrored", witness.mirrored}, {"word", witness.word},
            {"steps", steps}};
}
namespace {
bool suffix(const std::string &word, const std::string &ending) {
    return word.size() >= ending.size() &&
           word.compare(word.size() - ending.size(), ending.size(), ending) == 0;
}
Word monomial_word(const Json &terms, const Env &env, Int budget) {
    auto polynomial = expand_polynomial(terms, env, budget);
    if (polynomial.size() != 1)
        throw TurnaroundError("a monomial packet word is required");
    if (polynomial[0].q != 0 || polynomial[0].coefficient != 1)
        throw TurnaroundError("packet word must have coefficient one");
    return std::move(polynomial[0].word);
}
bool prefix_word(const Word &word, const Word &prefix) {
    return prefix.size() <= word.size() && std::equal(prefix.begin(), prefix.end(), word.begin());
}
struct TypedTail {
    Env env;
    Int pre_j, post_j, output_delta, input_delta;
    bool pre_mirrored, post_mirrored;
    Braid tail;
    std::string mirror_policy;
};
TypedTail typed_event(const ProcessedEvent &event, Int n, Int i, const Catalog &catalog,
                      const Limits &limits) {
    const auto &env = event.parameters;
    const auto &name = event.rule_id;
    if (env.size() != 6 || !env.count("n") || !env.count("i") || !env.count("j") ||
        !env.count("k") || !env.count("l") || !env.count("eps") || event.rhs_index < 0 ||
        (event.mirror_policy != "M" && event.mirror_policy != "bar_mirror"))
        throw TurnaroundError("invalid processed-tail event");
    if (env.at("n") != n || env.at("i") != i || name.compare(0, 4, "zip.") != 0 ||
        !(suffix(name, ".I") || suffix(name, ".M") || suffix(name, ".BM")) ||
        (suffix(name, ".BM") && event.mirror_policy != "bar_mirror") ||
        (suffix(name, ".M") && event.mirror_policy != "M"))
        throw TurnaroundError("tail event has incompatible rank or orientation");
    const auto &rule = catalog.get(name);
    if (!condition(rule.guard, env) || static_cast<size_t>(event.rhs_index) >= rule.rhs.size())
        throw TurnaroundError("tail event lacks a valid accepted rule branch");
    const bool mirror = suffix(name, ".M") || suffix(name, ".BM"), barred = suffix(name, ".BM");
    const auto old_j =
        mirror ? checked_add(checked_sub(checked_sub(n, i), env.at("j")), 2) : env.at("j");
    const auto &base = catalog.get(rule.base_id + ".I");
    if (static_cast<size_t>(event.rhs_index) >= base.rhs.size())
        throw TurnaroundError("base tail branch is absent");
    auto [jj, kk, mm, aa, bb] =
        branch_shape(base.rhs[static_cast<size_t>(event.rhs_index)], env, mirror);
    if (barred) {
        aa = checked_sub(0, aa);
        bb = checked_sub(0, bb);
    }
    const auto [output, inp] = anchors(i, jj, kk, mm);
    auto local = expand_word(rule.rhs[static_cast<size_t>(event.rhs_index)].at("word"), env,
                             limits.expansion_budget);
    auto prefix = word_from_braid(power(i - 1, output, aa, limits.max_word_size));
    append(prefix, core(n, i, jj, kk, mm, event.mirror_policy));
    append(prefix, word_from_braid(power(i - 1, inp, bb, limits.max_word_size)));
    if (!prefix_word(local, prefix))
        throw TurnaroundError("tail event does not have the catalog frame shape");
    return {env,
            old_j,
            jj,
            aa,
            bb,
            mirror,
            mm,
            braid_from_word(slice(local, prefix.size(), local.size())),
            event.mirror_policy};
}
Word turnaround_core(const Control &c, Int budget) {
    const auto [output, inp] = anchors(c.i, c.j, 0, c.mirrored);
    auto result = word_from_braid(power(c.i - 1, output, c.output_power, budget));
    result.push_back(Token::slider(c.i, c.j, c.mirrored));
    append(result, word_from_braid(power(c.i - 1, inp, c.input_power, budget)));
    return result;
}
} // namespace
Parts turnaround_parts(const Branch &branch, Int n, Int budget) {
    validate_word(branch.word, n);
    const auto &c = branch.control;
    if (c.phase != "turnaround_complete" ||
        (c.mirror_policy != "M" && c.mirror_policy != "bar_mirror") || c.i < 2 || c.i > n ||
        c.j < 1 || c.j > n - c.i + 1 || c.exterior_length < 0 || c.payload_length < 0 ||
        c.routing_length < 0 || c.transported_length < 0 || c.processed_length < 0 ||
        c.tail_length < 0 ||
        c.exterior_length !=
            checked_add(checked_add(c.payload_length, c.routing_length), c.transported_length))
        throw TurnaroundError("invalid turnaround output boundaries");
    const auto right_length = checked_add(c.tail_length, c.processed_length);
    if (static_cast<uint64_t>(right_length) > branch.word.size() ||
        static_cast<uint64_t>(c.exterior_length) > branch.word.size())
        throw TurnaroundError("turnaround output does not match its full expression");
    const auto end = branch.word.size() - static_cast<size_t>(right_length);
    auto exterior = braid_from_word(slice(branch.word, 0, c.exterior_length));
    auto processed =
        braid_from_word(slice(branch.word, end, end + static_cast<size_t>(c.processed_length)));
    auto tail = braid_from_word(
        slice(branch.word, end + static_cast<size_t>(c.processed_length), branch.word.size()));
    auto expected = word_from_braid(exterior);
    append(expected, turnaround_core(c, budget));
    append(expected, word_from_braid(processed));
    append(expected, word_from_braid(tail));
    if (expected != branch.word)
        throw TurnaroundError("turnaround output does not match its full expression");
    return {c, std::move(exterior), std::move(processed), std::move(tail)};
}
std::pair<Branch, Json> turnaround(const Branch &branch, Int n, const Catalog &catalog,
                                   const Limits &limits, bool trace, bool want_certificate) {
    check_branch(branch, limits);
    validate_word(branch.word, n);
    Json evidence = Json::array();
    Branch following;
    following.coefficient = branch.coefficient;
    Int event_count = 0;
    if (branch.control.i == 2) {
        const auto parts = upward_parts(branch, n);
        const auto &c = parts.control;
        if (!parts.u.empty())
            throw TurnaroundError("completed upward pass required");
        const auto band = transport_band(parts.p, n, 2, c.j);
        check_transport(band);
        auto exterior = friends(2, c.j, c.k);
        append(exterior, parts.p);
        auto &cc = following.control;
        cc.phase = "turnaround_complete";
        cc.i = 2;
        cc.j = band.end;
        cc.mirrored = false;
        cc.mirror_policy = "M";
        cc.exterior_length = static_cast<Int>(exterior.size());
        cc.payload_length = 0;
        cc.routing_length = c.k;
        cc.transported_length = static_cast<Int>(parts.p.size());
        cc.processed_length = 0;
        cc.tail_length = static_cast<Int>(parts.tail.size());
        following.word = word_from_braid(exterior);
        following.word.push_back(Token::slider(2, band.end));
        append(following.word, word_from_braid(parts.tail));
        if (trace)
            evidence.push_back({{"band_transport", band_json(band)}});
        event_count = static_cast<Int>(band.steps.size());
    } else {
        auto parts = framed_parts(branch, n, limits.max_word_size);
        const auto &c = parts.control;
        if (!parts.u.empty() || c.phase != "framed_upper_pass_complete")
            throw TurnaroundError("completed framed upward pass required");
        const auto i = c.i;
        auto j = c.j, a = c.output_power, b = c.input_power;
        auto mirror = c.mirrored;
        const auto &policy = c.mirror_policy;
        auto exterior = c.output_payload;
        auto core_word = core(n, i, j, c.k, mirror, policy);
        core_word.pop_back();
        append(exterior, braid_from_word(core_word));
        const auto payload_length = static_cast<Int>(c.output_payload.size());
        const auto routing_length = checked_mul(i - 1, c.k);
        if (c.processed_events.size() > static_cast<uint64_t>(limits.max_steps))
            throw LimitExceeded("turnaround event limit exceeded");
        auto p = parts.p;
        for (const auto &event : c.processed_events) {
            auto typed = typed_event(event, n, i, catalog, limits);
            if (typed.post_j != j || typed.post_mirrored != mirror ||
                typed.mirror_policy != policy || !segment(p, 0, typed.tail))
                throw TurnaroundError("tail event disagrees with current slider or word");
            auto T = std::move(typed.tail);
            p.erase(p.begin(), p.begin() + static_cast<std::ptrdiff_t>(T.size()));
            const auto aa = typed.output_delta, bb = typed.input_delta;
            Json packet_id = nullptr;
            Braid W, R;
            if (aa || bb) {
                if ((aa != -1 && aa != 1) || bb != -aa || typed.pre_mirrored == mirror)
                    throw TurnaroundError("invalid emitted flype frame increments");
                const bool barred = typed.pre_mirrored && policy == "bar_mirror";
                const std::string sign = aa == (barred ? 1 : -1) ? "positive" : "negative";
                const std::string variant = barred ? ".BM" : typed.pre_mirrored ? ".M" : ".I";
                const std::string packet_name = "zip.flype." + sign + variant;
                packet_id = packet_name;
                auto env = typed.env;
                env["j"] =
                    typed.pre_mirrored ? checked_add(checked_sub(checked_sub(n, i), j), 2) : j;
                env["k"] = 0;
                const auto &packet = catalog.get(packet_name);
                auto incoming = monomial_word(packet.lhs, env, limits.expansion_budget);
                if (incoming.empty() || incoming.back() != Token::slider(i, j, typed.pre_mirrored))
                    throw TurnaroundError("flype packet has an invalid input slider");
                incoming.pop_back();
                W = braid_from_word(incoming);
                const auto packet_word = monomial_word(packet.rhs, env, limits.expansion_budget);
                const auto [output, inp] = anchors(i, j, 0, mirror);
                auto packet_core = word_from_braid(power(i - 1, output, aa, limits.max_word_size));
                packet_core.push_back(Token::slider(i, j, mirror));
                append(packet_core, word_from_braid(power(i - 1, inp, bb, limits.max_word_size)));
                if (!prefix_word(packet_word, packet_core))
                    throw TurnaroundError("flype packet has an invalid retained frame");
                const auto C =
                    braid_from_word(slice(packet_word, packet_core.size(), packet_word.size()));
                if (C.size() != static_cast<uint64_t>(i - 1) || !segment(T, 0, C))
                    throw TurnaroundError("tail is missing its complete outgoing flype tape");
                R = slice(T, C.size(), T.size());
                a = checked_sub(a, aa);
                b = checked_sub(b, bb);
                mirror = typed.pre_mirrored;
            } else {
                if (typed.pre_mirrored != mirror)
                    throw TurnaroundError("orientation change requires an explicit flype packet");
                R = std::move(T);
            }
            const auto band = transport_band(R, n, i, j, mirror);
            check_transport(band);
            if (band.end != typed.pre_j)
                throw TurnaroundError("residual tail does not return to the previous band anchor");
            append(exterior, W);
            append(exterior, R);
            j = band.end;
            if (trace)
                evidence.push_back({{"event", event_json(event)},
                                    {"packet_rule", packet_id},
                                    {"band_transport", band_json(band)},
                                    {"output_power", a},
                                    {"input_power", b},
                                    {"j", j},
                                    {"mirrored", mirror}});
        }
        const auto &initial = c.initial_slider;
        if (j != initial.j || mirror != initial.mirrored || a != initial.output_power ||
            b != initial.input_power ||
            p.size() != static_cast<uint64_t>(initial.processed_length) ||
            policy != initial.mirror_policy)
            throw TurnaroundError("typed tails do not restore the initial slider data");
        auto &cc = following.control;
        cc.phase = "turnaround_complete";
        cc.i = i;
        cc.j = j;
        cc.mirrored = mirror;
        cc.mirror_policy = policy;
        cc.output_power = a;
        cc.input_power = b;
        cc.exterior_length = static_cast<Int>(exterior.size());
        cc.payload_length = payload_length;
        cc.routing_length = routing_length;
        cc.transported_length =
            checked_sub(checked_sub(cc.exterior_length, payload_length), routing_length);
        cc.processed_length = static_cast<Int>(p.size());
        cc.tail_length = static_cast<Int>(parts.tail.size());
        following.word = word_from_braid(exterior);
        append(following.word, turnaround_core(cc, limits.max_word_size));
        append(following.word, word_from_braid(p));
        append(following.word, word_from_braid(parts.tail));
        check_branch(following, limits);
        turnaround_parts(following, n, limits.max_word_size);
        event_count = static_cast<Int>(c.processed_events.size());
    }
    if (!want_certificate)
        return {std::move(following), Json::object()};
    Json dependencies = Json::object();
    for (const std::string name :
         {"framed.py", "frames.py", "transport.py", "catalog.py", "upward.py", "semilinear.py"})
        dependencies[name] = reference_source(name);
    Json certificate = {{"schema", "zipper-turnaround-v1"},
                        {"n", n},
                        {"trace_enabled", trace},
                        {"source_sha256", reference_source("turnaround.py")},
                        {"dependency_sha256", dependencies},
                        {"catalog_sha256", reference_catalog_sha()},
                        {"initial", branch_json(branch)},
                        {"final", branch_json(following)},
                        {"events", evidence},
                        {"summary",
                         {{"events", event_count},
                          {"word_tokens_before", branch.word.size()},
                          {"word_tokens_after", following.word.size()}}},
                        {"limits", limits_json(limits)}};
    return {std::move(following), std::move(certificate)};
}
} // namespace zipper
