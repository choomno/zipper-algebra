#include "scheduler.hpp"
#include "algorithm.hpp"
#include "typed_work_size.hpp"
#include "typed_canonical.hpp"
#include "catalog.hpp"
#include <algorithm>
#include <memory>
#include <set>
#include <unordered_set>

namespace zipper_verify {
namespace {
using Work = typed::Work;
using Output = typed::Output;
using Transition = typed::Transition;
using Context = typed::Context;
using Labels = typed::Labels;
J encode(const Work &w) {
    return typed::encode_work(w);
}
J encode(const Output &o) {
    return typed::encode_output(o);
}
std::string state_key(const Work &w) {
    return typed::continuation_key(w);
}
std::string state_key(const Output &o) {
    return typed::output_collection_key(o);
}
void set_scalar(Work &w, Laurent scalar) {
    w.coefficient = std::move(scalar);
    if (auto *s = std::get_if<typed::ZipUp>(&w.phase))
        s->branch.coefficient = w.coefficient;
    if (auto *s = std::get_if<typed::ZipDown>(&w.phase))
        s->branch.coefficient = w.coefficient;
}
void set_scalar(Output &o, Laurent scalar) {
    o.coefficient = std::move(scalar);
}
template <class Value> J encode_values(const std::vector<Value> &values) {
    J result = J::array();
    for (const auto &value : values)
        result.push_back(encode(value));
    return result;
}
// A sequence determines the schedule; a separate index gives exact lookup.
// Replacing an existing element leaves its position unchanged, as in dict.
template <class Value> class OrderedStates {
    using List = std::list<std::pair<std::string, Value>>;
    List sequence_;
    std::unordered_map<std::string, typename List::iterator> lookup_;

  public:
    const Value *get(const std::string &key) const {
        auto it = lookup_.find(key);
        return it == lookup_.end() ? nullptr : &it->second->second;
    }
    bool contains(const std::string &key) const {
        return lookup_.count(key) != 0;
    }
    I size() const {
        return static_cast<I>(sequence_.size());
    }
    bool empty() const {
        return sequence_.empty();
    }
    void set(const std::string &key, Value value) {
        auto found = lookup_.find(key);
        if (found != lookup_.end())
            found->second->second = std::move(value);
        else {
            sequence_.emplace_back(key, std::move(value));
            auto last = sequence_.end();
            --last;
            lookup_.emplace(key, last);
        }
    }
    void erase(const std::string &key) {
        auto it = lookup_.find(key);
        if (it != lookup_.end()) {
            sequence_.erase(it->second);
            lookup_.erase(it);
        }
    }
    std::vector<const std::pair<std::string, Value> *> prefix(I count) const {
        std::vector<const std::pair<std::string, Value> *> result;
        for (const auto &item : sequence_) {
            if (static_cast<I>(result.size()) == count)
                break;
            result.push_back(&item);
        }
        return result;
    }
    std::vector<Value> values() const {
        std::vector<Value> result;
        for (const auto &item : sequence_)
            result.push_back(item.second);
        return result;
    }
    J object() const {
        J result = J::object();
        for (const auto &item : sequence_)
            result[item.first] = encode(item.second);
        return result;
    }
    const List &items() const {
        return sequence_;
    }
    List &items() {
        return sequence_;
    }
};
template <class Value> struct Staged {
    // A disengaged value records a cancellation, preserving update order.
    using Update = std::pair<std::string, std::optional<Value>>;
    std::vector<Update> updates;
    I merged = 0, zeroes = 0;
    J object() const {
        J out = J::object();
        for (const auto &[key, value] : updates)
            out[key] = value ? encode(*value) : J(nullptr);
        return out;
    }
};
template <class Value>
Staged<Value> stage(const OrderedStates<Value> &current, const std::vector<const Value *> &incoming,
                    const std::unordered_set<std::string> &removed) {
    Staged<Value> result;
    std::unordered_map<std::string, std::size_t> index;
    for (const Value *pointer : incoming) {
        const Value &item = *pointer;
        auto key = state_key(item);
        auto found = index.find(key);
        const Value *previous = nullptr;
        if (found != index.end()) {
            const auto &old = result.updates[found->second].second;
            if (old)
                previous = &*old;
        } else if (!removed.count(key))
            previous = current.get(key);
        Value combined = item;
        Laurent scalar = item.coefficient;
        if (previous) {
            result.merged = checked_add(result.merged, 1);
            scalar = typed::add_scalar(previous->coefficient, scalar);
            if constexpr (std::is_same_v<Value, Work>)
                combined.depth = std::max(previous->depth, item.depth);
        }
        set_scalar(combined, std::move(scalar));
        std::optional<Value> value;
        if (!combined.coefficient.empty())
            value = std::move(combined);
        if (found == index.end()) {
            index.emplace(key, result.updates.size());
            result.updates.emplace_back(std::move(key), std::move(value));
        } else
            result.updates[found->second].second = std::move(value);
    }
    for (const auto &[key, value] : result.updates)
        if (!value)
            ++result.zeroes;
    return result;
}
template <class Value>
I count_after(const OrderedStates<Value> &current, const Staged<Value> &updates,
              const std::unordered_set<std::string> &removed) {
    I value = checked_sub(current.size(), static_cast<I>(removed.size()));
    for (const auto &[key, item] : updates.updates) {
        value = checked_add(value, item ? 1 : 0);
        if (current.contains(key) && !removed.count(key))
            value = checked_sub(value, 1);
    }
    return value;
}
template <class Value>
void commit(OrderedStates<Value> &current, Staged<Value> &updates,
            const std::unordered_set<std::string> &removed = {}) {
    for (const auto &key : removed)
        current.erase(key);
    for (auto &[key, item] : updates.updates) {
        if (item)
            current.set(key, std::move(*item));
        else
            current.erase(key);
    }
}
I item_size(const Output &o) {
    return typed::length(o.word);
}
I item_size(const Work &w) {
    I count = std::visit(
        [](const auto &phase) -> I {
            using P = std::decay_t<decltype(phase)>;
            if constexpr (std::is_same_v<P, typed::Single>)
                return checked_add(checked_add(typed::length(phase.left), 1),
                                   typed::length(phase.processed));
            else if constexpr (std::is_same_v<P, typed::ZipUp> || std::is_same_v<P, typed::ZipDown>)
                return typed::length(phase.branch.word);
            else
                return typed::length(phase.word);
        },
        w.phase);
    count = checked_add(count, typed::length(w.boundary_left));
    count = checked_add(count, typed::length(w.boundary_right));
    for (const auto &c : w.contexts)
        for (const auto *word : {&c.top_route, &c.component_word, &c.bottom_route})
            count = checked_add(count, typed::length(*word));
    return count;
}
void check_scalar(const Laurent &scalar, const Limits &limits) {
    if (static_cast<I>(scalar.size()) > limits.max_coefficient_support)
        throw LimitExceeded("coefficient support limit exceeded");
    for (const auto &[q, value] : scalar)
        if (integer_bits(value) > limits.max_integer_bits)
            throw LimitExceeded("coefficient integer bit limit exceeded");
}
void check_item(const Work &w, const Limits &limits) {
    if (typed::work_expression_size(w) > limits.max_word_size)
        throw LimitExceeded("word size limit exceeded");
    check_scalar(w.coefficient, limits);
}
void check_item(const Output &o, const Limits &limits) {
    if (typed::length(o.word) > limits.max_word_size)
        throw LimitExceeded("word size limit exceeded");
    check_scalar(o.coefficient, limits);
}
void update_max(J &metrics, const char *key, I value) {
    metrics[key] = std::max(as_i(metrics.at(key)), value);
}
template <class Value> void update_item_metrics(J &metrics, const Value &item) {
    update_max(metrics, "max_local_word_tokens", item_size(item));
    update_max(metrics, "max_coefficient_support", static_cast<I>(item.coefficient.size()));
    for (const auto &[q, c] : item.coefficient) {
        metrics["q_exponent_min"] = std::min(as_i(metrics.at("q_exponent_min")), q);
        update_max(metrics, "q_exponent_max", q);
        update_max(metrics, "max_integer_bits", integer_bits(c));
    }
}
std::string event_name(const J &event) {
    J name = nullptr;
    if (event.contains("rule_id"))
        name = event.at("rule_id");
    else if (event.contains("catalog_transition")) {
        const auto &transition = event.at("catalog_transition");
        if (transition.contains("request") && transition.at("request").contains("rule_id"))
            name = transition.at("request").at("rule_id");
    }
    if (name.is_null())
        return "None";
    if (!name.is_string())
        throw ReplayError("transition rule identifier must be a string");
    return name.get<std::string>();
}

// Retain the whole control state. Only the scalar, depth, extracted outer
// components, and the names of the local endpoint ports are external data.
Work intrinsic_state(Work state) {
    set_scalar(state, {{0, 1}});
    state.depth = 0;
    state.contexts.clear();
    state.bottom_labels.clear();
    for (I k = 1; k <= state.rank; ++k)
        state.bottom_labels.push_back(k);
    state.top_labels = state.bottom_labels;
    return state;
}
Labels endpoint_names(const Labels &ports, const Labels &names) {
    Labels result;
    for (I p : ports) {
        if (p < 1 || p > static_cast<I>(names.size()))
            throw ReplayError("cached endpoint port is outside the caller labels");
        result.push_back(names.at(static_cast<std::size_t>(p - 1)));
    }
    return result;
}
struct Restoration {
    I rank;
    Laurent coefficient;
    std::vector<Context> contexts;
    Labels bottom_labels, top_labels;
    explicit Restoration(Work &w)
        : rank(w.rank), coefficient(std::move(w.coefficient)), contexts(std::move(w.contexts)),
          bottom_labels(std::move(w.bottom_labels)), top_labels(std::move(w.top_labels)) {}
};
struct RestorationView {
    I rank;
    const Laurent &coefficient;
    const std::vector<Context> &contexts;
    const Labels &bottom_labels, &top_labels;
    RestorationView(const Work &w)
        : rank(w.rank), coefficient(w.coefficient), contexts(w.contexts),
          bottom_labels(w.bottom_labels), top_labels(w.top_labels) {}
    RestorationView(const Restoration &w)
        : rank(w.rank), coefficient(w.coefficient), contexts(w.contexts),
          bottom_labels(w.bottom_labels), top_labels(w.top_labels) {}
};
void name_record(typed::ComponentRecord &r, const RestorationView &caller) {
    r.original_bottom_ports = endpoint_names(r.original_bottom_ports, caller.bottom_labels);
    r.original_top_ports = endpoint_names(r.original_top_ports, caller.top_labels);
}
void restore_components(std::vector<Context> &components, const RestorationView &caller) {
    for (auto &component : components)
        name_record(component.record, caller);
    std::vector<Context> all = caller.contexts;
    for (auto &component : components)
        all.push_back(std::move(component));
    components = std::move(all);
}
Output restore_output(Output output, const RestorationView &caller) {
    restore_components(output.components, caller);
    output.word = typed::lift_contexts(output.word, caller.contexts);
    output.coefficient = typed::multiply_scalar(caller.coefficient, output.coefficient);
    return output;
}
void collect_output(OrderedStates<Output> &output, Output term, const Limits &limits) {
    check_item(term, limits);
    const std::string key = state_key(term);
    if (const auto *old = output.get(key))
        term.coefficient = typed::add_scalar(old->coefficient, term.coefficient);
    if (term.coefficient.empty())
        output.erase(key);
    else {
        check_item(term, limits);
        output.set(key, std::move(term));
    }
    if (output.size() > limits.max_terms)
        throw LimitExceeded("completed local output term limit exceeded");
}
} // namespace

std::string continuation_key(const J &item, bool pending) {
    return continuation_json(item, pending);
}
ReplaySession::ReplaySession(Catalog &catalog, std::size_t entries, std::size_t bytes,
                             bool unary_fusion)
    : catalog_(catalog), unary_fusion_(unary_fusion), max_entries_(entries), max_bytes_(bytes) {
    catalog_.configure_selection_cache(std::min<std::size_t>(entries, 4096),
                                       std::min<std::size_t>(bytes, 4 * 1024 * 1024));
    catalog_.configure_bound_cache(std::min<std::size_t>(entries, 4096),
                                   std::min<std::size_t>(bytes, 4 * 1024 * 1024));
}
void ReplaySession::clear_cache() {
    catalog_.clear_selection_cache();
    catalog_.clear_bound_cache();
    index_.clear();
    cache_.clear();
    bytes_ = 0;
    input_ = 0;
    hits_ = misses_ = evictions_ = oversized_ = cross_input_hits_ = complement_hits_ =
        peak_entries_ = peak_bytes_ = 0;
    for (auto &entry : complete_cache_)
        entry.second.result->cache_references = 0;
    complete_index_.clear();
    complete_cache_.clear();
    complete_bytes_ = 0;
    complete_key_bytes_ = complete_payload_bytes_ = complete_unique_outputs_ = 0;
    complete_hits_ = complete_misses_ = complete_evictions_ = complete_oversized_ = 0;
    complete_peak_bytes_ = complete_peak_entries_ = complete_cross_input_hits_ = 0;
}
J ReplaySession::cache_info() const {
    return {{"entries", cache_.size()},
            {"bytes", bytes_},
            {"max_entries", max_entries_},
            {"max_bytes", max_bytes_},
            {"hits", hits_},
            {"misses", misses_},
            {"evaluations", misses_},
            {"evictions", evictions_},
            {"oversized", oversized_},
            {"cross_input_hits", cross_input_hits_},
            {"complement_hits", complement_hits_},
            {"peak_entries", peak_entries_},
            {"peak_bytes", peak_bytes_},
            {"byte_measure", "canonical-key-and-payload-estimate"},
            {"mode", "independent-cpp-intrinsic-unit-replay-v2"},
            {"completed",
             {{"entries", complete_cache_.size()},
              {"bytes", complete_bytes_},
              {"key_bytes", complete_key_bytes_},
              {"payload_bytes", complete_payload_bytes_},
              {"unique_outputs", complete_unique_outputs_},
              {"max_entries", max_entries_},
              {"max_bytes", max_bytes_},
              {"hits", complete_hits_},
              {"misses", complete_misses_},
              {"evictions", complete_evictions_},
              {"oversized", complete_oversized_},
              {"cross_input_hits", complete_cross_input_hits_},
              {"peak_entries", complete_peak_entries_},
              {"peak_bytes", complete_peak_bytes_},
              {"byte_measure", "canonical keys per entry plus each shared immutable payload once; "
                               "excludes container overhead"}}},
            {"selection", catalog_.selection_cache_info()},
            {"bound_rules", catalog_.bound_cache_info()},
            {"memory_scope", "separate estimated payload budget per cache; excludes container and "
                             "active-frame overhead"}};
}
Transition ReplaySession::evaluate(const Work &work, const Limits &limits, bool trace) {
    Work unit = trace ? work : intrinsic_state(work);
    set_scalar(unit, {{0, 1}});
    unit.depth = 0;
    std::string key = "[" + canonical_json(limits.json()) + "," + (trace ? "true" : "false") + "," +
                      typed::canonical_work(unit) + "]";
    auto found = index_.find(key);
    Transition result;
    if (found != index_.end()) {
        auto item = found->second;
        ++hits_;
        if (item->second.input != input_)
            ++cross_input_hits_;
        if (std::holds_alternative<typed::Extract>(work.phase))
            ++complement_hits_;
        result = item->second.result;
        cache_.splice(cache_.end(), cache_, item);
    } else {
        ++misses_;
        result = typed::algorithm_transition(unit, catalog_, limits, trace);
        if (max_bytes_ && max_entries_) {
            std::size_t amount =
                key.size() + canonical_json(typed::encode_transition(result)).size();
            if (amount > max_bytes_)
                ++oversized_;
            else {
                while (!cache_.empty() &&
                       (cache_.size() >= max_entries_ || bytes_ > max_bytes_ - amount)) {
                    auto &old = cache_.front();
                    bytes_ -= old.second.bytes;
                    index_.erase(old.first);
                    cache_.pop_front();
                    ++evictions_;
                }
                cache_.emplace_back(key, Entry{result, amount, input_});
                auto item = cache_.end();
                --item;
                try {
                    index_.emplace(key, item);
                } catch (...) {
                    cache_.pop_back();
                    throw;
                }
                bytes_ += amount;
                peak_bytes_ = std::max(peak_bytes_, static_cast<I>(bytes_));
                peak_entries_ = std::max(peak_entries_, static_cast<I>(cache_.size()));
            }
        }
    }
    if (result.error.empty()) {
        const I next_depth = checked_add(work.depth, 1);
        for (auto &item : result.children) {
            if (!trace) {
                restore_components(item.contexts, work);
                item.bottom_labels = endpoint_names(item.bottom_labels, work.bottom_labels);
                item.top_labels = endpoint_names(item.top_labels, work.top_labels);
            }
            set_scalar(item, typed::multiply_scalar(work.coefficient, item.coefficient));
            item.depth = next_depth;
        }
        for (auto &item : result.output) {
            if (!trace)
                item = restore_output(std::move(item), work);
            else
                item.coefficient = typed::multiply_scalar(work.coefficient, item.coefficient);
        }
        if (!trace && result.evidence.contains("component")) {
            auto record = typed::decode_record(result.evidence.at("component"));
            name_record(record, work);
            result.evidence["component"] = typed::encode_record(record);
        }
    }
    return result;
}
void ReplaySession::retain_complete(const std::string &key,
                                    const std::shared_ptr<const CompleteOutput> &output) {
    if (!max_bytes_ || !max_entries_)
        return;
    const std::size_t key_bytes = key.size();
    if (output->serialized_bytes > max_bytes_ ||
        key_bytes > max_bytes_ - output->serialized_bytes) {
        ++complete_oversized_;
        return;
    }
    // Recompute after every eviction: removing the last cached alias of this
    // very payload changes the cost of retaining it from key-only to key+sum.
    std::size_t payload_bytes = output->cache_references ? 0 : output->serialized_bytes;
    while (true) {
        if (payload_bytes > max_bytes_ || key_bytes > max_bytes_ - payload_bytes) {
            ++complete_oversized_;
            return;
        }
        if (complete_cache_.size() < max_entries_ &&
            complete_bytes_ <= max_bytes_ - (key_bytes + payload_bytes))
            break;
        const auto &old = complete_cache_.front();
        complete_bytes_ -= old.second.bytes;
        complete_key_bytes_ -= old.second.bytes;
        if (--old.second.result->cache_references == 0) {
            complete_bytes_ -= old.second.result->serialized_bytes;
            complete_payload_bytes_ -= old.second.result->serialized_bytes;
            --complete_unique_outputs_;
        }
        complete_index_.erase(old.first);
        complete_cache_.pop_front();
        ++complete_evictions_;
        payload_bytes = output->cache_references ? 0 : output->serialized_bytes;
    }
    complete_cache_.emplace_back(key, CompleteEntry{output, key_bytes, input_});
    auto added = complete_cache_.end();
    --added;
    try {
        complete_index_.emplace(key, added);
    } catch (...) {
        complete_cache_.pop_back();
        throw;
    }
    if (output->cache_references++ == 0) {
        complete_payload_bytes_ += output->serialized_bytes;
        ++complete_unique_outputs_;
    }
    complete_key_bytes_ += key_bytes;
    complete_bytes_ += key_bytes + payload_bytes;
    complete_peak_bytes_ = std::max(complete_peak_bytes_, static_cast<I>(complete_bytes_));
    complete_peak_entries_ =
        std::max(complete_peak_entries_, static_cast<I>(complete_cache_.size()));
}

J ReplaySession::replay_reused(I n, const Word &word, const Limits &limits) {
    if (n < 1 || !word.is_array())
        throw ReplayError("positive integer rank and literal word required");
    for (const J &token : word)
        if (!is_int(token) || as_i(token) == 0 || abs_i(as_i(token)) >= n)
            throw ReplayError("literal input must contain signed generators in the stated rank");
    Work initial = typed::initial_work(n, typed::decode_word(word));
    input_ = checked_add(input_, 1);
    // A fused link is a logical frame: its intrinsic key stays active, and
    // completion still restores and checks this exact link before cache admission.
    struct UnaryLink {
        Restoration caller;
        std::string key;
    };
    struct Frame {
        std::vector<UnaryLink> chain;
        Restoration caller;
        Work unit;
        std::string key;
        std::vector<Work> children;
        std::size_t next = 0;
        OrderedStates<Output> output;
        std::shared_ptr<const CompleteOutput> sole_result;
        bool expanded = false;
        explicit Frame(Restoration value) : caller(std::move(value)) {}
    };
    std::vector<std::unique_ptr<Frame>> stack;
    std::unordered_set<std::string> active;
    I expanded = 0, peak_stack = 0, references = 0, cache_hits = 0, live_terms = 0,
      peak_live_terms = 0, logical_stack = 0;
    J counts = J::object(), answer = J::array(), diagnostic = nullptr;
    std::string status = "complete";
    const std::string limits_key = canonical_json(limits.json());
    auto push = [&](Work caller, bool fuse = false) {
        check_item(caller, limits);
        Restoration restoration(caller);
        Work unit = intrinsic_state(std::move(caller));
        std::string key = limits_key + "\n" + state_key(unit);
        if (active.count(key))
            throw ControlObligation("intrinsic completion encountered a control cycle");
        active.insert(key);
        if (fuse) {
            auto &frame = *stack.back();
            frame.chain.push_back(UnaryLink{std::move(frame.caller), std::move(frame.key)});
            // Run the next unary microstep in this same scheduler frame.
            frame.caller = std::move(restoration);
            frame.unit = std::move(unit);
            frame.key = std::move(key);
            frame.children.clear();
            frame.next = 0;
            frame.expanded = false;
        } else {
            auto frame = std::make_unique<Frame>(std::move(restoration));
            frame->unit = std::move(unit);
            frame->key = std::move(key);
            stack.push_back(std::move(frame));
        }
        logical_stack = checked_add(logical_stack, 1);
        peak_stack = std::max(peak_stack, logical_stack);
        references = checked_add(references, 1);
    };
    auto account = [&](I change) {
        live_terms = checked_add(live_terms, change);
        peak_live_terms = std::max(peak_live_terms, live_terms);
        if (live_terms > limits.max_terms)
            throw LimitExceeded("active completion terms limit exceeded");
    };
    auto identity_restoration = [](const Restoration &caller) {
        bool identity = caller.contexts.empty() && caller.coefficient == Laurent{{0, 1}} &&
                        caller.top_labels == caller.bottom_labels &&
                        static_cast<I>(caller.bottom_labels.size()) == caller.rank;
        I port = 1;
        for (I name : caller.bottom_labels)
            identity = identity && name == port++;
        return identity;
    };
    auto make_output = [&](const OrderedStates<Output> &output) {
        auto terms = output.values();
        const std::size_t bytes =
            (max_entries_ && max_bytes_) ? canonical_json(encode_values(terms)).size() : 0;
        return std::make_shared<const CompleteOutput>(CompleteOutput{std::move(terms), bytes});
    };
    auto finish = [&](std::shared_ptr<const CompleteOutput> local) {
        Restoration caller = std::move(stack.back()->caller);
        account(-stack.back()->output.size());
        active.erase(stack.back()->key);
        --logical_stack;
        // Each compact parent completes in exactly the old postorder. In
        // particular, scalar/context maps are applied separately: cancellation
        // cannot hide an intermediate overflow or resource refusal.
        auto chain = std::move(stack.back()->chain);
        stack.pop_back();
        while (!chain.empty()) {
            auto link = std::move(chain.back());
            chain.pop_back();
            I collected = 0;
            OrderedStates<Output> restored;
            if (identity_restoration(caller)) {
                for (const auto &term : local->terms)
                    check_item(term, limits);
            } else {
                for (const auto &term : local->terms) {
                    const I before = restored.size();
                    collect_output(restored, restore_output(term, caller), limits);
                    account(checked_sub(restored.size(), before));
                }
                collected = restored.size();
            }
            account(-1);
            if (!identity_restoration(caller))
                local = make_output(restored);
            retain_complete(link.key, local);
            account(-collected);
            active.erase(link.key);
            --logical_stack;
            caller = std::move(link.caller);
        }
        if (stack.empty()) {
            OrderedStates<Output> root;
            for (const auto &term : local->terms)
                collect_output(root, restore_output(term, caller), limits);
            answer = encode_values(root.values());
        } else {
            auto &parent = *stack.back();
            if (identity_restoration(caller) && parent.children.size() == 1 &&
                parent.output.empty()) {
                for (const auto &term : local->terms)
                    check_item(term, limits);
                parent.sole_result = local;
            } else {
                for (const auto &term : local->terms) {
                    const I before = parent.output.size();
                    collect_output(parent.output, restore_output(term, caller), limits);
                    account(checked_sub(parent.output.size(), before));
                }
            }
            account(-1);
        }
    };
    try {
        push(std::move(initial));
        while (!stack.empty()) {
            Frame &frame = *stack.back();
            if (!frame.expanded) {
                auto found = complete_index_.find(frame.key);
                if (found != complete_index_.end()) {
                    auto entry = found->second;
                    ++complete_hits_;
                    ++cache_hits;
                    if (entry->second.input != input_)
                        ++complete_cross_input_hits_;
                    auto local = entry->second.result;
                    complete_cache_.splice(complete_cache_.end(), complete_cache_, entry);
                    finish(local);
                    continue;
                }
                ++complete_misses_;
                if (expanded >= limits.max_steps)
                    throw LimitExceeded("completion transition evaluation limit exceeded");
                check_item(frame.unit, limits);
                ++expanded;
                Transition result =
                    typed::algorithm_transition(frame.unit, catalog_, limits, false);
                if (!result.error.empty()) {
                    status = result.error;
                    diagnostic = {{"message", result.message},
                                  {"phase", typed::phase_name(frame.unit)}};
                    break;
                }
                const std::string name = event_name(result.evidence);
                counts[name] = checked_add(counts.contains(name) ? as_i(counts[name]) : 0, 1);
                // Collect sibling continuations before recursively evaluating
                // them, including exact coefficient cancellations.
                OrderedStates<Work> pending;
                std::vector<const Work *> children;
                for (const auto &child : result.children)
                    children.push_back(&child);
                auto staged = stage(pending, children, {});
                commit(pending, staged);
                frame.children = pending.values();
                account(static_cast<I>(frame.children.size()));
                for (const auto &child : frame.children)
                    check_item(child, limits);
                for (auto &term : result.output) {
                    const I before = frame.output.size();
                    collect_output(frame.output, std::move(term), limits);
                    account(checked_sub(frame.output.size(), before));
                }
                if (checked_add(static_cast<I>(frame.children.size()), frame.output.size()) >
                    limits.max_terms)
                    throw LimitExceeded("local transition term limit exceeded");
                // The complete unit state has been consumed; descendants need
                // its key, restoration map, and staged children only.
                frame.unit = Work{};
                frame.expanded = true;
            }
            if (frame.next < frame.children.size()) {
                Work child = std::move(frame.children.at(frame.next++));
                push(std::move(child),
                     unary_fusion_ && frame.children.size() == 1 && frame.output.empty());
            } else {
                std::shared_ptr<const CompleteOutput> local = frame.sole_result;
                if (!local) {
                    local = make_output(frame.output);
                }
                retain_complete(frame.key, local);
                finish(local);
            }
        }
    } catch (const LimitExceeded &error) {
        status = "incomplete";
        diagnostic = {{"message", error.what()}};
    } catch (const ControlObligation &error) {
        status = "control_obligation";
        diagnostic = {{"message", error.what()}};
    }
    if (status != "complete")
        answer = J::array();
    return {{"status", status},
            {"diagnostic", diagnostic},
            {"output", answer},
            {"metrics",
             {{"expanded_nodes", expanded},
              {"node_references", references},
              {"completed_cache_hits", cache_hits},
              {"peak_stack", peak_stack},
              {"peak_live_terms", peak_live_terms},
              {"output_terms", length(answer)},
              {"rule_applications", counts}}},
            {"cache", cache_info()}};
}

J ReplaySession::replay(I n, const Word &word, const Limits &limits, I batch_size, bool trace) {
    if (n < 1 || batch_size < 1 || !word.is_array())
        throw ReplayError("positive integer rank and batch size required");
    for (const auto &token : word)
        if (!is_int(token) || as_i(token) == 0 || abs_i(as_i(token)) >= n)
            throw ReplayError("a literal signed-generator list in the stated rank is required");
    input_ = checked_add(input_, 1);
    auto initial = typed::initial_work(n, typed::decode_word(word));
    OrderedStates<Work> frontier;
    OrderedStates<Output> output;
    frontier.set(state_key(initial), initial);
    J transitions = J::array(), counts = J::object();
    J metrics = {{"steps", 0},
                 {"attempted_steps", 0},
                 {"batches", 0},
                 {"max_terms", 1},
                 {"max_local_word_tokens", checked_add(length(word), 1)},
                 {"max_coefficient_support", 1},
                 {"max_integer_bits", 1},
                 {"q_exponent_min", 0},
                 {"q_exponent_max", 0},
                 {"maximum_depth", 0},
                 {"max_branch_count", 0},
                 {"merged_terms", 0},
                 {"zero_terms", 0}};
    std::string status = "complete";
    J diagnostic = nullptr;
    while (!frontier.empty()) {
        auto selected = frontier.prefix(batch_size);
        const Work *before = &selected.front()->second;
        std::vector<Transition> results;
        Staged<Work> pending;
        Staged<Output> completed;
        I terms = 0;
        std::unordered_set<std::string> removed;
        try {
            if (checked_add(as_i(metrics["steps"]), static_cast<I>(selected.size())) >
                limits.max_steps)
                throw LimitExceeded("step limit would be exceeded by the next atomic batch");
            for (const auto &item : selected)
                check_item(item->second, limits);
            metrics["attempted_steps"] =
                checked_add(as_i(metrics["attempted_steps"]), static_cast<I>(selected.size()));
            for (const auto &item : selected)
                results.push_back(evaluate(item->second, limits, trace));
            bool failed = false;
            for (std::size_t i = 0; i < selected.size(); ++i)
                if (!results[i].error.empty()) {
                    before = &selected[i]->second;
                    status = results[i].error;
                    diagnostic = {{"message", results[i].message},
                                  {"phase", typed::phase_name(*before)},
                                  {"state", encode(*before)},
                                  {"batch_rolled_back", true}};
                    failed = true;
                    break;
                }
            if (failed)
                break;
            std::vector<const Work *> children;
            std::vector<const Output *> added;
            for (const auto &item : selected)
                removed.insert(item->first);
            for (const auto &result : results) {
                for (const auto &child : result.children)
                    children.push_back(&child);
                for (const auto &term : result.output)
                    added.push_back(&term);
            }
            pending = stage(frontier, children, removed);
            completed = stage(output, added, {});
            terms = checked_add(count_after(frontier, pending, removed),
                                count_after(output, completed, {}));
            if (terms > limits.max_terms)
                throw LimitExceeded("collected term limit exceeded");
            for (const auto &[key, item] : pending.updates)
                if (item)
                    check_item(*item, limits);
            for (const auto &[key, item] : completed.updates)
                if (item)
                    check_item(*item, limits);
        } catch (const LimitExceeded &error) {
            status = "incomplete";
            diagnostic = {{"message", error.what()},
                          {"phase", "batch"},
                          {"batch_representative", encode(*before)},
                          {"batch_rolled_back", true}};
            break;
        }
        metrics["steps"] = checked_add(as_i(metrics["steps"]), static_cast<I>(selected.size()));
        metrics["batches"] = checked_add(as_i(metrics["batches"]), 1);
        metrics["merged_terms"] = checked_add(as_i(metrics["merged_terms"]),
                                              checked_add(pending.merged, completed.merged));
        metrics["zero_terms"] =
            checked_add(as_i(metrics["zero_terms"]), checked_add(pending.zeroes, completed.zeroes));
        update_max(metrics, "max_terms", terms);
        for (std::size_t i = 0; i < selected.size(); ++i) {
            auto name = event_name(results[i].evidence);
            counts[name] = checked_add(counts.contains(name) ? as_i(counts[name]) : 0, 1);
            update_max(metrics, "maximum_depth", checked_add(selected[i]->second.depth, 1));
            update_max(metrics, "max_branch_count",
                       checked_add(static_cast<I>(results[i].children.size()),
                                   static_cast<I>(results[i].output.size())));
        }
        for (const auto &[key, item] : pending.updates)
            if (item)
                update_item_metrics(metrics, *item);
        for (const auto &[key, item] : completed.updates)
            if (item)
                update_item_metrics(metrics, *item);
        if (trace) {
            J predecessors = J::array();
            for (const auto &item : selected)
                predecessors.push_back(encode(item->second));
            J rewrites = J::array();
            for (const auto &result : results)
                rewrites.push_back(typed::encode_transition(result));
            transitions.push_back({{"batch", checked_sub(as_i(metrics["batches"]), 1)},
                                   {"predecessors", predecessors},
                                   {"rewrites", rewrites},
                                   {"pending_updates", pending.object()},
                                   {"output_updates", completed.object()}});
        }
        commit(frontier, pending, removed);
        commit(output, completed);
    }
    metrics["rule_applications"] = counts;
    metrics["output_terms"] = output.size();
    metrics["pending_terms"] = frontier.size();
    J certificate = {{"schema", "zipper-collected-run-v1"},
                     {"scope", "literal_word_collected_spanning_procedure"},
                     {"n", n},
                     {"input_word", word},
                     {"coefficient_ring", "Z[q,q^-1]"},
                     {"source_sha256", reference_sources()},
                     {"catalog_sha256", catalog_.sha256},
                     {"trace_enabled", trace},
                     {"limits", limits.json()},
                     {"execution",
                      {{"batch_size", batch_size},
                       {"schedule", "stable-prefix-atomic-v1"},
                       {"collection", "full-continuation-except-scalar-and-depth-v1"},
                       {"local_term_floor", 4},
                       {"local_support_floor", 2}}},
                     {"status", status}};
    certificate["diagnostic"] = std::move(diagnostic);
    certificate["pending"] = encode_values(frontier.values());
    certificate["output"] = encode_values(output.values());
    certificate["transitions"] = std::move(transitions);
    certificate["metrics"] = std::move(metrics);
    return certificate;
}
} // namespace zipper_verify
