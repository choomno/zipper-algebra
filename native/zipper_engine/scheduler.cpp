#include "scheduler.hpp"
#include "ordered.hpp"
#include <algorithm>
#include <chrono>
#include <optional>
#include <unordered_set>

namespace zipper {
namespace {
struct SchedulerLimit : LimitExceeded {
    using LimitExceeded::LimitExceeded;
};
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
struct PhaseTimer {
    EngineSession *session;
    const char *phase;
    Clock::time_point start = Clock::now();
    bool active = true;
    void stop() {
        if (active && session)
            session->record_timing(phase, elapsed(start));
        active = false;
    }
    ~PhaseTimer() {
        stop();
    }
};
using Removed = std::unordered_set<std::string>;

template <class Item> using Updates = OrderedMap<std::optional<Item>>;

void set_coefficient(Work &item, Laurent coefficient) {
    item.coefficient = std::move(coefficient);
    if (item.phase == "zip_up" || item.phase == "zip_down")
        item.branch.coefficient = item.coefficient;
}
void set_coefficient(Output &item, Laurent coefficient) {
    item.coefficient = std::move(coefficient);
}
void merge_depth(Work &item, const Work &previous) {
    item.depth = std::max(item.depth, previous.depth);
}
void merge_depth(Output &, const Output &) {}
std::string key_of(const Work &item) {
    return work_key(item);
}
std::string key_of(const Output &item) {
    return output_key(item);
}
size_t size_of(const Work &item) {
    return work_size(item);
}
size_t size_of(const Output &item) {
    return item.word.size();
}
Json json_of(const Work &item) {
    return work_json(item);
}
Json json_of(const Output &item) {
    return output_json(item);
}

int64_t scalar_bits(const BigInt &value) {
    if (value == 0)
        return 0;
    BigInt absolute = value < 0 ? -value : value;
    return static_cast<int64_t>(boost::multiprecision::msb(absolute)) + 1;
}
Laurent add(const Laurent &a, const Laurent &b) {
    Laurent result = a;
    for (const auto &[q, c] : b)
        result[q] += c;
    for (auto it = result.begin(); it != result.end();) {
        if (it->second == 0)
            it = result.erase(it);
        else
            ++it;
    }
    return result;
}

template <class Item>
void stage_one(const OrderedMap<Item> &current, Updates<Item> &updates, const Item &item,
               const Removed &removed, int64_t &merges) {
    auto key = key_of(item);
    const Item *previous = nullptr;
    auto updated = updates.get(key);
    if (updated) {
        if (*updated)
            previous = &(**updated);
    } else if (!removed.count(key))
        previous = current.get(key);
    Item combined = item;
    if (previous) {
        ++merges;
        set_coefficient(combined, add(previous->coefficient, item.coefficient));
        merge_depth(combined, *previous);
    }
    if (combined.coefficient.empty())
        updates.set(std::move(key), std::nullopt);
    else
        updates.set(std::move(key), std::move(combined));
}

template <class Item> int64_t zero_count(const Updates<Item> &updates) {
    return std::count_if(updates.begin(), updates.end(),
                         [](const auto &entry) { return !entry.second; });
}
template <class Item>
int64_t count_after(const OrderedMap<Item> &current, const Updates<Item> &updates,
                    const Removed &removed) {
    int64_t result = static_cast<int64_t>(current.size()) - static_cast<int64_t>(removed.size());
    for (const auto &[key, value] : updates)
        result += int(bool(value)) - int(current.contains(key) && !removed.count(key));
    return result;
}
template <class Item>
void commit(OrderedMap<Item> &current, Updates<Item> &updates, const Removed &removed) {
    for (const auto &key : removed)
        current.erase(key);
    for (auto &[key, value] : updates) {
        if (value)
            current.set(key, std::move(*value));
        else
            current.erase(key);
    }
}
template <class Item> void check_item(const Item &item, const Limits &limits) {
    if (size_of(item) > static_cast<size_t>(limits.max_word_size))
        throw SchedulerLimit("word size limit exceeded");
    if (item.coefficient.size() > static_cast<size_t>(limits.max_coefficient_support))
        throw SchedulerLimit("coefficient support limit exceeded");
    for (const auto &[q, c] : item.coefficient) {
        (void)q;
        if (scalar_bits(c) > limits.max_integer_bits)
            throw SchedulerLimit("coefficient integer bit limit exceeded");
    }
}

struct Metrics {
    int64_t steps = 0, attempted_steps = 0, batches = 0, max_terms = 1, max_local_word_tokens = 1;
    int64_t max_coefficient_support = 1, max_integer_bits = 1, q_exponent_min = 0,
            q_exponent_max = 0;
    int64_t maximum_depth = 0, max_branch_count = 0, merged_terms = 0, zero_terms = 0;
    std::map<std::string, int64_t> counts;
    template <class Item> void observe(const Item &item) {
        max_local_word_tokens =
            std::max(max_local_word_tokens, static_cast<int64_t>(size_of(item)));
        max_coefficient_support =
            std::max(max_coefficient_support, static_cast<int64_t>(item.coefficient.size()));
        for (const auto &[q, c] : item.coefficient) {
            q_exponent_min = std::min(q_exponent_min, q);
            q_exponent_max = std::max(q_exponent_max, q);
            max_integer_bits = std::max(max_integer_bits, scalar_bits(c));
        }
    }
    Json json(size_t pending, size_t output) const {
        return {{"steps", steps},
                {"attempted_steps", attempted_steps},
                {"batches", batches},
                {"max_terms", max_terms},
                {"max_local_word_tokens", max_local_word_tokens},
                {"max_coefficient_support", max_coefficient_support},
                {"max_integer_bits", max_integer_bits},
                {"q_exponent_min", q_exponent_min},
                {"q_exponent_max", q_exponent_max},
                {"maximum_depth", maximum_depth},
                {"max_branch_count", max_branch_count},
                {"merged_terms", merged_terms},
                {"zero_terms", zero_terms},
                {"rule_applications", counts},
                {"output_terms", output},
                {"pending_terms", pending}};
    }
};
std::string event_rule(const Json &event) {
    if (event.contains("rule_id"))
        return event.at("rule_id").get<std::string>();
    if (event.contains("catalog_transition") &&
        event.at("catalog_transition").contains("request")) {
        const auto &request = event.at("catalog_transition").at("request");
        if (request.contains("rule_id"))
            return request.at("rule_id").get<std::string>();
    }
    return "None";
}
// Python's trace dictionary uses canonical textual continuation keys. The
// native hot path uses binary keys; textual keys are materialized only here.
std::string trace_key(Json item, bool pending) {
    item.erase("coefficient");
    if (pending) {
        item.erase("depth");
        if (item.contains("branch"))
            item["branch"].erase("coefficient");
    }
    return item.dump();
}
template <class Item>
Json trace_updates(const Updates<Item> &updates, const OrderedMap<Item> &current, bool pending,
                   const std::unordered_map<std::string, Json> &representatives) {
    Json result = Json::object();
    for (const auto &[key, value] : updates) {
        Json record;
        if (value)
            record = json_of(*value);
        else {
            auto found = representatives.find(key);
            if (found != representatives.end())
                record = found->second;
            else if (auto item = current.get(key))
                record = json_of(*item);
            else
                throw std::logic_error("missing canceled trace representative");
        }
        result[trace_key(record, pending)] = value ? record : Json(nullptr);
    }
    return result;
}
} // namespace

Json reduce_collected(const Work &initial, Catalog &catalog, const Limits &limits, bool trace,
                      size_t batch_size, EngineSession *session) {
    if (batch_size == 0)
        throw std::invalid_argument("batch size must be positive");
    OrderedMap<Work> frontier;
    frontier.set(work_key(initial), initial);
    OrderedMap<Output> output;
    Metrics metrics;
    metrics.max_local_word_tokens = static_cast<int64_t>(initial.left.size()) + 1;
    Json transitions = Json::array(), diagnostic = nullptr;
    std::string status = "complete";
    const Removed none;
    while (!frontier.empty()) {
        std::vector<std::pair<std::string, Work>> selected;
        selected.reserve(std::min(batch_size, frontier.size()));
        for (const auto &item : frontier) {
            selected.push_back(item);
            if (selected.size() == batch_size)
                break;
        }
        const Work *before = &selected.front().second;
        Updates<Work> pending_updates;
        Updates<Output> output_updates;
        Removed removed;
        std::vector<Transition> results;
        int64_t pm = 0, om = 0, term_count = 0;
        Int next_steps = 0, next_batches = 0, next_merged = 0, next_zeroes = 0;
        Int next_maximum_depth = metrics.maximum_depth;
        Json batch_trace;
        try {
            if (static_cast<uint64_t>(limits.max_steps - metrics.steps) < selected.size())
                throw SchedulerLimit("step limit would be exceeded by the next atomic batch");
            for (const auto &[key, item] : selected) {
                (void)key;
                check_item(item, limits);
            }
            metrics.attempted_steps += static_cast<int64_t>(selected.size());
            results.reserve(selected.size());
            auto transition_start = Clock::now();
            // Every transition in an admitted batch is evaluated even if an
            // earlier one refuses, matching the reference's atomic batch.
            for (const auto &[key, item] : selected) {
                (void)key;
                results.push_back(session
                                      ? session->transition(item, limits, trace)
                                      : algorithm_transition(item, catalog, limits, trace, trace));
            }
            if (session)
                session->record_timing("transition_seconds", elapsed(transition_start));
            bool failed = false;
            for (size_t index = 0; index < results.size(); ++index) {
                const auto &result = results[index];
                if (!result.error.empty()) {
                    before = &selected[index].second;
                    status = result.error;
                    diagnostic = {{"message", result.message},
                                  {"phase", before->phase},
                                  {"state", work_json(*before)},
                                  {"batch_rolled_back", true}};
                    failed = true;
                    break;
                }
            }
            if (failed)
                break;
            PhaseTimer collection_timer{session, "collection_seconds"};
            for (const auto &[key, item] : selected) {
                (void)item;
                removed.insert(key);
            }
            std::unordered_map<std::string, Json> pending_representatives, output_representatives;
            for (const auto &result : results) {
                for (const auto &item : result.children) {
                    stage_one(frontier, pending_updates, item, removed, pm);
                    if (trace)
                        pending_representatives[key_of(item)] = work_json(item);
                }
                for (const auto &item : result.output) {
                    stage_one(output, output_updates, item, none, om);
                    if (trace)
                        output_representatives[key_of(item)] = output_json(item);
                }
            }
            term_count = checked_add(count_after(frontier, pending_updates, removed),
                                     count_after(output, output_updates, none));
            if (term_count > limits.max_terms)
                throw SchedulerLimit("collected term limit exceeded");
            for (const auto &[key, item] : pending_updates) {
                (void)key;
                if (item)
                    check_item(*item, limits);
            }
            for (const auto &[key, item] : output_updates) {
                (void)key;
                if (item)
                    check_item(*item, limits);
            }
            // All signed diagnostic arithmetic is staged before committing
            // the batch, including output-only transitions at the depth limit.
            next_steps = checked_add(metrics.steps, static_cast<Int>(selected.size()));
            next_batches = checked_add(metrics.batches, 1);
            next_merged = checked_add(metrics.merged_terms, checked_add(pm, om));
            next_zeroes = checked_add(metrics.zero_terms, checked_add(zero_count(pending_updates),
                                                                      zero_count(output_updates)));
            for (const auto &selected_item : selected)
                next_maximum_depth =
                    std::max(next_maximum_depth, checked_add(selected_item.second.depth, 1));
            collection_timer.stop();
            if (trace) {
                auto trace_start = Clock::now();
                batch_trace = {{"batch", metrics.batches},
                               {"predecessors", Json::array()},
                               {"rewrites", Json::array()},
                               {"pending_updates", trace_updates(pending_updates, frontier, true,
                                                                 pending_representatives)},
                               {"output_updates", trace_updates(output_updates, output, false,
                                                                output_representatives)}};
                for (const auto &[key, item] : selected) {
                    (void)key;
                    batch_trace["predecessors"].push_back(work_json(item));
                }
                for (const auto &result : results)
                    batch_trace["rewrites"].push_back(transition_json(result));
                if (session)
                    session->record_timing("trace_materialization_seconds", elapsed(trace_start));
            }
        } catch (const LimitExceeded &exc) {
            status = "incomplete";
            diagnostic = {{"message", exc.what()},
                          {"phase", "batch"},
                          {"batch_representative", work_json(*before)},
                          {"batch_rolled_back", true}};
            break;
        }
        metrics.steps = next_steps;
        metrics.batches = next_batches;
        metrics.merged_terms = next_merged;
        metrics.zero_terms = next_zeroes;
        metrics.maximum_depth = next_maximum_depth;
        metrics.max_terms = std::max(metrics.max_terms, term_count);
        for (size_t i = 0; i < results.size(); ++i) {
            const auto &result = results[i];
            ++metrics.counts[event_rule(result.event)];
            metrics.max_branch_count =
                std::max(metrics.max_branch_count,
                         static_cast<int64_t>(result.children.size() + result.output.size()));
        }
        for (const auto &[key, item] : pending_updates) {
            (void)key;
            if (item)
                metrics.observe(*item);
        }
        for (const auto &[key, item] : output_updates) {
            (void)key;
            if (item)
                metrics.observe(*item);
        }
        commit(frontier, pending_updates, removed);
        commit(output, output_updates, none);
        if (trace)
            transitions.push_back(std::move(batch_trace));
    }
    auto boundary_start = Clock::now();
    Json pending = Json::array(), completed = Json::array();
    for (const auto &[key, item] : frontier) {
        (void)key;
        pending.push_back(work_json(item));
    }
    for (const auto &[key, item] : output) {
        (void)key;
        completed.push_back(output_json(item));
    }
    if (session)
        session->record_timing("boundary_materialization_seconds", elapsed(boundary_start));
    return {{"status", status},
            {"diagnostic", diagnostic},
            {"pending", pending},
            {"output", completed},
            {"transitions", transitions},
            {"metrics", metrics.json(frontier.size(), output.size())},
            {"execution",
             {{"batch_size", batch_size},
              {"schedule", "stable-prefix-atomic-v1"},
              {"collection", "full-continuation-except-scalar-and-depth-v1"},
              {"local_term_floor", 4},
              {"local_support_floor", 2}}}};
}
namespace {
void append_number(std::string &key, uint64_t value) {
    for (unsigned i = 0; i < 8; ++i)
        key.push_back(static_cast<char>((value >> (8 * i)) & 255));
}
std::string cache_key(const Work &item, const Limits &limits, bool trace) {
    std::string key = work_key(item);
    for (Int value :
         {limits.max_steps, limits.max_terms, limits.max_word_size, limits.max_coefficient_support,
          limits.max_integer_bits, limits.expansion_budget})
        append_number(key, static_cast<uint64_t>(value));
    key.push_back(trace ? 1 : 0);
    return key;
}
size_t scalar_charge(const Laurent &coefficient) {
    size_t total = 0;
    for (const auto &[q, c] : coefficient) {
        (void)q;
        total += 8 + (scalar_bits(c) + 7) / 8 + 1;
    }
    return total;
}
size_t json_charge(const Json &value) {
    size_t total = sizeof(Json);
    if (value.is_object())
        for (auto it = value.begin(); it != value.end(); ++it)
            total += it.key().size() + json_charge(it.value());
    else if (value.is_array())
        for (const auto &item : value)
            total += json_charge(item);
    else if (value.is_string())
        total += value.get_ref<const std::string &>().size();
    return total;
}
size_t transition_charge(const Transition &result) {
    size_t total = sizeof(Transition) + json_charge(result.event);
    for (const auto &child : result.children)
        total += sizeof(Work) + work_key(child).size() + scalar_charge(child.coefficient);
    for (const auto &item : result.output)
        total += sizeof(Output) + output_key(item).size() + scalar_charge(item.coefficient);
    return total;
}
// The unit transition is local to its active component. Extracted outer
// components and original endpoint names are restored only on return. Every
// geometric boundary, phase, frame, and processed-event field remains keyed.
Work intrinsic_work(const Work &work) {
    Work unit;
    unit.rank = work.rank;
    unit.coefficient = {{0, 1}};
    unit.world_phi = work.world_phi;
    unit.boundary_left = work.boundary_left;
    unit.boundary_right = work.boundary_right;
    unit.boundary_left_factors = work.boundary_left_factors;
    unit.boundary_right_factors = work.boundary_right_factors;
    unit.bottom_labels.resize(work.rank);
    for (Int i = 0; i < work.rank; ++i)
        unit.bottom_labels[i] = i + 1;
    unit.top_labels = unit.bottom_labels;
    unit.phase = work.phase;
    if (work.phase == "single") {
        unit.left = work.left;
        unit.processed = work.processed;
        unit.position = work.position;
    } else if (work.phase == "zip_up" || work.phase == "zip_down") {
        unit.branch = work.branch;
        unit.branch.coefficient = unit.coefficient;
    } else {
        unit.word = work.word;
        unit.payload_length = work.payload_length;
        unit.completed_sides = work.completed_sides;
    }
    return unit;
}
Braid relabel_ports(const Braid &ports, const Braid &labels) {
    Braid restored;
    restored.reserve(ports.size());
    for (Int port : ports) {
        if (port < 1 || static_cast<size_t>(port) > labels.size())
            throw EngineError("cached endpoint port lies outside the caller labels");
        restored.push_back(labels[port - 1]);
    }
    return restored;
}
void relabel_record(Json &record, const Work &work) {
    record["original_bottom_ports"] =
        relabel_ports(record.at("original_bottom_ports").get<Braid>(), work.bottom_labels);
    record["original_top_ports"] =
        relabel_ports(record.at("original_top_ports").get<Braid>(), work.top_labels);
}
void restore_contexts(std::vector<Context> &contexts, const Work &work) {
    for (auto &context : contexts)
        relabel_record(context.record, work);
    contexts.insert(contexts.begin(), work.contexts.begin(), work.contexts.end());
}
Word lift_contexts(Word word, const std::vector<Context> &contexts) {
    for (auto context = contexts.rbegin(); context != contexts.rend(); ++context) {
        Word next;
        next.reserve(context->top_route.size() + context->component_word.size() + word.size() +
                     context->bottom_route.size());
        for (Int letter : context->top_route)
            next.emplace_back(letter);
        append(next, context->component_word);
        for (auto token : word) {
            if (token.is_slider())
                token.start = checked_add(token.start, context->size);
            else
                token.generator = (token.generator > 0 ? 1 : -1) *
                                  checked_add(std::abs(token.generator), context->size);
            next.push_back(token);
        }
        for (Int letter : context->bottom_route)
            next.emplace_back(letter);
        word = std::move(next);
    }
    return word;
}
void rebase(Transition &result, const Work &work) {
    const Int next_depth = checked_add(work.depth, 1);
    for (auto &child : result.children) {
        restore_contexts(child.contexts, work);
        child.bottom_labels = relabel_ports(child.bottom_labels, work.bottom_labels);
        child.top_labels = relabel_ports(child.top_labels, work.top_labels);
        set_coefficient(child, multiply(work.coefficient, child.coefficient));
        child.depth = next_depth;
    }
    for (auto &item : result.output) {
        restore_contexts(item.components, work);
        item.word = lift_contexts(std::move(item.word), work.contexts);
        set_coefficient(item, multiply(work.coefficient, item.coefficient));
    }
    if (result.event.contains("component"))
        relabel_record(result.event["component"], work);
}
} // namespace
struct EngineSession::Impl {
    struct Entry {
        Transition result;
        size_t bytes;
        uint64_t input;
    };
    Catalog &catalog;
    size_t max_entries, max_bytes, bytes = 0, peak_bytes = 0, peak_entries = 0;
    uint64_t hits = 0, misses = 0, evictions = 0, oversized = 0, evaluations = 0;
    uint64_t input = 0, cross_input_hits = 0, complement_hits = 0;
    std::map<std::string, double> timings;
    OrderedMap<Entry> cache;
    Impl(Catalog &c, size_t entries, size_t memory)
        : catalog(c), max_entries(entries), max_bytes(memory) {}
};
EngineSession::EngineSession(Catalog &catalog, size_t cache_entries, size_t cache_bytes)
    : impl_(std::make_unique<Impl>(catalog, cache_entries, cache_bytes)) {}
EngineSession::~EngineSession() = default;
void EngineSession::clear_cache() {
    auto &p = *impl_;
    p.cache.clear();
    p.bytes = p.peak_bytes = p.peak_entries = 0;
    p.hits = p.misses = p.evictions = p.oversized = p.evaluations = 0;
    p.cross_input_hits = p.complement_hits = 0;
}
Json EngineSession::cache_info() const {
    const auto &p = *impl_;
    return {
        {"hits", p.hits},
        {"misses", p.misses},
        {"cross_input_hits", p.cross_input_hits},
        {"complement_hits", p.complement_hits},
        {"evictions", p.evictions},
        {"oversized", p.oversized},
        {"evaluations", p.evaluations},
        {"entries", p.cache.size()},
        {"bytes", p.bytes},
        {"peak_bytes", p.peak_bytes},
        {"peak_entries", p.peak_entries},
        {"max_entries", p.max_entries},
        {"max_bytes", p.max_bytes},
        {"byte_measure",
         "estimated-typed-state-and-structural-key-payload; excludes allocator overhead"},
        {"key_scope", "intrinsic-unit-transition-except-outer-contexts-labels-scalar-and-depth"}};
}
Transition EngineSession::transition(const Work &work, const Limits &limits, bool trace) {
    auto &p = *impl_;
    if (!p.max_entries || !p.max_bytes) {
        ++p.misses;
        ++p.evaluations;
        return algorithm_transition(work, p.catalog, limits, trace, trace);
    }
    auto key_start = Clock::now();
    auto unit = intrinsic_work(work);
    auto key = cache_key(unit, limits, trace);
    record_timing("cache_key_seconds", elapsed(key_start));
    Transition result;
    if (auto entry = p.cache.get(key)) {
        ++p.hits;
        p.cross_input_hits += entry->input != p.input;
        p.complement_hits += !work.contexts.empty();
        // Move the existing entry to the LRU tail. Copy the return template
        // first so no pointer survives the erase.
        result = entry->result;
        auto retained = std::move(*entry);
        p.cache.erase(key);
        p.cache.set(key, std::move(retained));
    } else {
        ++p.misses;
        ++p.evaluations;
        result = algorithm_transition(unit, p.catalog, limits, trace, trace);
        if (!result.error.empty())
            return result;
        const size_t charge = key.size() + transition_charge(result);
        if (charge > p.max_bytes)
            ++p.oversized;
        else {
            while (!p.cache.empty() &&
                   (p.cache.size() >= p.max_entries || charge > p.max_bytes - p.bytes)) {
                const auto &oldest = *p.cache.begin();
                auto old_key = oldest.first;
                p.bytes -= oldest.second.bytes;
                p.cache.erase(old_key);
                ++p.evictions;
            }
            p.cache.set(key, Impl::Entry{result, charge, p.input});
            p.bytes += charge;
            p.peak_bytes = std::max(p.peak_bytes, p.bytes);
            p.peak_entries = std::max(p.peak_entries, p.cache.size());
        }
    }
    try {
        rebase(result, work);
    } catch (const LimitExceeded &error) {
        Transition refused;
        refused.error = "incomplete";
        refused.message = error.what();
        return refused;
    }
    return result;
}
Json EngineSession::last_telemetry() const {
    return {{"seconds", impl_->timings},
            {"timing_scope", "native reduce; transition includes cache lookup and keys; collection "
                             "includes staging and checks; excludes process IO"}};
}
void EngineSession::record_timing(const char *phase, double seconds) {
    impl_->timings[phase] += seconds;
}
Json EngineSession::reduce(const Work &initial, const Limits &limits, bool trace,
                           size_t batch_size) {
    impl_->timings.clear();
    ++impl_->input;
    auto start = Clock::now();
    auto result = reduce_collected(initial, impl_->catalog, limits, trace, batch_size, this);
    record_timing("total_seconds", elapsed(start));
    return result;
}
} // namespace zipper
