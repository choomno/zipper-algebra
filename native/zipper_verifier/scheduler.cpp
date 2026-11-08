#include "scheduler.hpp"
#include "algorithm.hpp"
#include "catalog.hpp"
#include <algorithm>
#include <memory>
#include <set>
#include <unordered_set>

namespace zipper_verify {
namespace {
// A sequence determines the schedule; a separate index gives exact lookup.
// Replacing an existing element leaves its position unchanged, as in dict.
class OrderedStates {
    using List = std::list<std::pair<std::string, J>>;
    List sequence_;
    std::unordered_map<std::string, List::iterator> lookup_;
public:
    const J *get(const std::string &key) const {
        auto it = lookup_.find(key);
        return it == lookup_.end() ? nullptr : &it->second->second;
    }
    bool contains(const std::string &key) const { return lookup_.count(key) != 0; }
    I size() const { return static_cast<I>(sequence_.size()); }
    bool empty() const { return sequence_.empty(); }
    void set(const std::string &key, J value) {
        auto found = lookup_.find(key);
        if (found != lookup_.end()) found->second->second = std::move(value);
        else {
            sequence_.emplace_back(key, std::move(value));
            auto last = sequence_.end(); --last;
            lookup_.emplace(key, last);
        }
    }
    void erase(const std::string &key) {
        auto it = lookup_.find(key);
        if (it != lookup_.end()) { sequence_.erase(it->second); lookup_.erase(it); }
    }
    std::vector<std::pair<std::string,J>> prefix(I count) const {
        std::vector<std::pair<std::string,J>> result;
        for (const auto &item : sequence_) {
            if (static_cast<I>(result.size()) == count) break;
            result.push_back(item);
        }
        return result;
    }
    J values() const {
        J result = J::array();
        for (const auto &item : sequence_) result.push_back(item.second);
        return result;
    }
    J object() const {
        J result = J::object();
        for (const auto &item : sequence_) result[item.first] = item.second;
        return result;
    }
    const List &items() const { return sequence_; }
    List &items() { return sequence_; }
};
struct Staged {
    OrderedStates updates;
    I merged = 0, zeroes = 0;
};
Staged stage(const OrderedStates &current, const std::vector<const J *> &incoming, bool pending,
             const std::unordered_set<std::string> &removed) {
    Staged result;
    for (const J *pointer : incoming) {
        const J &item=*pointer;
        auto key = continuation_key(item, pending);
        const J *previous = result.updates.get(key);
        if (previous == nullptr && !removed.count(key)) previous = current.get(key);
        J combined = item;
        J coef = item.at("coefficient");
        if (previous != nullptr && !previous->is_null()) {
            result.merged = checked_add(result.merged, 1);
            coef = add_coefficients(previous->at("coefficient"), coef);
            if (pending) combined["depth"] = std::max(as_i(previous->at("depth")), as_i(item.at("depth")));
        }
        set_coefficient(combined, coef);
        result.updates.set(key, coef.empty() ? J(nullptr) : std::move(combined));
    }
    for (const auto &item : result.updates.items()) if (item.second.is_null()) ++result.zeroes;
    return result;
}
I count_after(const OrderedStates &current, const OrderedStates &updates,
              const std::unordered_set<std::string> &removed) {
    I value = checked_sub(current.size(), static_cast<I>(removed.size()));
    for (const auto &item : updates.items()) {
        value = checked_add(value, item.second.is_null() ? 0 : 1);
        if (current.contains(item.first) && !removed.count(item.first)) value = checked_sub(value, 1);
    }
    return value;
}
void commit(OrderedStates &current, OrderedStates &updates,
            const std::unordered_set<std::string> &removed={}) {
    for (const auto &key : removed) current.erase(key);
    for (auto &item : updates.items()) {
        if (item.second.is_null()) current.erase(item.first);
        else current.set(item.first, std::move(item.second));
    }
}
I item_size(const J &item, bool pending) {
    if (!pending) return length(item.at("word"));
    // Internal states have validated coordinates. Mirror and shift maps preserve
    // token count; count the concatenations without constructing their words.
    const auto &phase=item.at("phase");
    I count;
    if (phase=="single") count=checked_add(checked_add(length(item.at("left")),1),length(item.at("processed")));
    else if (phase=="zip_up" || phase=="zip_down") count=length(item.at("branch").at("word"));
    else count=length(item.at("word"));
    count=checked_add(count,length(item.at("boundary_left")));
    count=checked_add(count,length(item.at("boundary_right")));
    for (const auto &context:item.at("contexts"))
        for (const char *field:{"top_route","component_word","bottom_route"})
            count=checked_add(count,length(context.at(field)));
    return count;
}
void check_item(const J &item, const Limits &limits, bool pending) {
    if ((pending ? length(work_expression(item)) : length(item.at("word"))) > limits.max_word_size) throw LimitExceeded("word size limit exceeded");
    auto coef = coefficient(item.at("coefficient"));
    if (static_cast<I>(coef.size()) > limits.max_coefficient_support) throw LimitExceeded("coefficient support limit exceeded");
    for (const auto &term : coef) if (integer_bits(term.second) > limits.max_integer_bits)
        throw LimitExceeded("coefficient integer bit limit exceeded");
}
void update_max(J &metrics, const char *key, I value) {
    metrics[key] = std::max(as_i(metrics.at(key)), value);
}
void update_item_metrics(J &metrics, const J &item, bool pending) {
    update_max(metrics, "max_local_word_tokens", item_size(item, pending));
    auto coef = coefficient(item.at("coefficient"));
    update_max(metrics, "max_coefficient_support", static_cast<I>(coef.size()));
    for (const auto &[q,c] : coef) {
        metrics["q_exponent_min"] = std::min(as_i(metrics.at("q_exponent_min")), q);
        update_max(metrics, "q_exponent_max", q);
        update_max(metrics, "max_integer_bits", integer_bits(c));
    }
}
std::string event_name(const J &event) {
    J name = nullptr;
    if (event.contains("rule_id")) name = event.at("rule_id");
    else if (event.contains("catalog_transition")) {
        const auto &transition = event.at("catalog_transition");
        if (transition.contains("request") && transition.at("request").contains("rule_id"))
            name = transition.at("request").at("rule_id");
    }
    if (name.is_null()) return "None";
    if (!name.is_string()) throw ReplayError("transition rule identifier must be a string");
    return name.get<std::string>();
}

// Retain the whole control state. Only the scalar, depth, extracted outer
// components, and the names of the local endpoint ports are external data.
J intrinsic_state(const J &work) {
    J state=work;
    set_coefficient(state,one());
    state["depth"]=0;
    state["contexts"]=J::array();
    J ports=J::array();
    for (I k=1;k<=as_i(work.at("rank"));++k) ports.push_back(k);
    state["bottom_labels"]=ports;
    state["top_labels"]=std::move(ports);
    return state;
}
J endpoint_names(const J &ports,const J &names) {
    J result=J::array();
    for (const J &port:ports) {
        const I p=as_i(port);
        if (p<1 || p>length(names)) throw ReplayError("cached endpoint port is outside the caller labels");
        result.push_back(names.at(static_cast<std::size_t>(p-1)));
    }
    return result;
}
void name_record(J &record,const J &caller) {
    record["original_bottom_ports"]=endpoint_names(record.at("original_bottom_ports"),caller.at("bottom_labels"));
    record["original_top_ports"]=endpoint_names(record.at("original_top_ports"),caller.at("top_labels"));
}
void restore_components(J &components,const J &caller) {
    for (J &component:components) name_record(component.at("record"),caller);
    J all=caller.at("contexts");
    for (J &component:components) all.push_back(std::move(component));
    components=std::move(all);
}
J restore_output(J output,const J &caller) {
    restore_components(output["components"],caller);
    output["word"]=lift_contexts(output.at("word"),caller.at("contexts"));
    set_coefficient(output,multiply(caller.at("coefficient"),output.at("coefficient")));
    return output;
}
void collect_output(OrderedStates &output,J term,const Limits &limits) {
    check_item(term,limits,false);
    const std::string key=continuation_key(term,false);
    if (const J *old=output.get(key)) set_coefficient(term,add_coefficients(old->at("coefficient"),term.at("coefficient")));
    if (term.at("coefficient").empty()) output.erase(key);
    else {
        check_item(term,limits,false);
        output.set(key,std::move(term));
    }
    if (output.size()>limits.max_terms) throw LimitExceeded("completed local output term limit exceeded");
}
} // namespace

std::string continuation_key(const J &item, bool pending) {
    J identity = item;
    identity.erase("coefficient");
    if (pending) {
        identity.erase("depth");
        if (identity.contains("branch")) identity["branch"].erase("coefficient");
    }
    return canonical_json(identity);
}
ReplaySession::ReplaySession(Catalog &catalog, std::size_t entries, std::size_t bytes)
    : catalog_(catalog), max_entries_(entries), max_bytes_(bytes) {
    catalog_.configure_selection_cache(std::min<std::size_t>(entries,4096),
                                       std::min<std::size_t>(bytes,4*1024*1024));
    catalog_.configure_bound_cache(std::min<std::size_t>(entries,4096),
                                   std::min<std::size_t>(bytes,4*1024*1024));
}
void ReplaySession::clear_cache() {
    catalog_.clear_selection_cache();
    catalog_.clear_bound_cache();
    index_.clear(); cache_.clear(); bytes_=0; input_=0;
    hits_=misses_=evictions_=oversized_=cross_input_hits_=complement_hits_=peak_entries_=peak_bytes_=0;
    for (auto &entry:complete_cache_) entry.second.result->cache_references=0;
    complete_index_.clear(); complete_cache_.clear(); complete_bytes_=0;
    complete_key_bytes_=complete_payload_bytes_=complete_unique_outputs_=0;
    complete_hits_=complete_misses_=complete_evictions_=complete_oversized_=0;
    complete_peak_bytes_=complete_peak_entries_=complete_cross_input_hits_=0;
}
J ReplaySession::cache_info() const {
    return {{"entries",cache_.size()},{"bytes",bytes_},{"max_entries",max_entries_},{"max_bytes",max_bytes_},
            {"hits",hits_},{"misses",misses_},{"evaluations",misses_},{"evictions",evictions_},
            {"oversized",oversized_},{"cross_input_hits",cross_input_hits_},{"complement_hits",complement_hits_},
            {"peak_entries",peak_entries_},{"peak_bytes",peak_bytes_},
            {"byte_measure","canonical-key-and-payload-estimate"},
            {"mode","independent-cpp-intrinsic-unit-replay-v2"},
            {"completed",{{"entries",complete_cache_.size()},{"bytes",complete_bytes_},
                {"key_bytes",complete_key_bytes_},{"payload_bytes",complete_payload_bytes_},
                {"unique_outputs",complete_unique_outputs_},
                {"max_entries",max_entries_},{"max_bytes",max_bytes_},
                {"hits",complete_hits_},{"misses",complete_misses_},{"evictions",complete_evictions_},
                {"oversized",complete_oversized_},{"cross_input_hits",complete_cross_input_hits_},
                {"peak_entries",complete_peak_entries_},{"peak_bytes",complete_peak_bytes_},
                {"byte_measure","canonical keys per entry plus each shared immutable payload once; excludes container overhead"}}},
            {"selection",catalog_.selection_cache_info()},
            {"bound_rules",catalog_.bound_cache_info()},
            {"memory_scope","separate estimated payload budget per cache; excludes container and active-frame overhead"}};
}
J ReplaySession::evaluate(const J &work, const Limits &limits, bool trace) {
    J unit=trace?work:intrinsic_state(work);
    set_coefficient(unit, one()); unit["depth"]=0;
    std::string key=canonical_json(J::array({limits.json(),trace,unit}));
    auto found=index_.find(key);
    J result;
    if (found!=index_.end()) {
        auto item=found->second;
        ++hits_;
        if (item->second.input!=input_) ++cross_input_hits_;
        if (work.at("phase")=="extract") ++complement_hits_;
        result=item->second.result;
        cache_.splice(cache_.end(),cache_,item);
    } else {
        ++misses_;
        result=algorithm_transition(unit,catalog_,limits,trace);
        if (max_bytes_ && max_entries_) {
            std::size_t amount=key.size()+canonical_json(result).size();
            if (amount>max_bytes_) ++oversized_;
            else {
                while (!cache_.empty() && (cache_.size()>=max_entries_ || bytes_>max_bytes_-amount)) {
                    auto &old=cache_.front(); bytes_-=old.second.bytes;
                    index_.erase(old.first); cache_.pop_front(); ++evictions_;
                }
                cache_.emplace_back(key,Entry{result,amount,input_});
                auto item=cache_.end(); --item;
                try { index_.emplace(key,item); }
                catch (...) { cache_.pop_back(); throw; }
                bytes_+=amount;
                peak_bytes_=std::max(peak_bytes_,static_cast<I>(bytes_));
                peak_entries_=std::max(peak_entries_,static_cast<I>(cache_.size()));
            }
        }
    }
    if (!result.contains("error")) {
        const I next_depth=checked_add(as_i(work.at("depth")),1);
        for (auto &item : result.at("children")) {
            if (!trace) {
                restore_components(item["contexts"],work);
                item["bottom_labels"]=endpoint_names(item.at("bottom_labels"),work.at("bottom_labels"));
                item["top_labels"]=endpoint_names(item.at("top_labels"),work.at("top_labels"));
            }
            set_coefficient(item,multiply(work.at("coefficient"),item.at("coefficient")));
            item["depth"]=next_depth;
        }
        for (auto &item : result.at("output")) {
            if (!trace) item=restore_output(std::move(item),work);
            else set_coefficient(item,multiply(work.at("coefficient"),item.at("coefficient")));
        }
        if (!trace && result.at("unit_evidence").contains("component"))
            name_record(result["unit_evidence"]["component"],work);
    }
    return result;
}
void ReplaySession::retain_complete(const std::string &key,const std::shared_ptr<const CompleteOutput> &output) {
    if (!max_bytes_ || !max_entries_) return;
    const std::size_t key_bytes=key.size();
    if (output->serialized_bytes>max_bytes_ || key_bytes>max_bytes_-output->serialized_bytes) {
        ++complete_oversized_; return;
    }
    // Recompute after every eviction: removing the last cached alias of this
    // very payload changes the cost of retaining it from key-only to key+sum.
    std::size_t payload_bytes=output->cache_references?0:output->serialized_bytes;
    while (true) {
        if (payload_bytes>max_bytes_ || key_bytes>max_bytes_-payload_bytes) {
            ++complete_oversized_; return;
        }
        if (complete_cache_.size()<max_entries_ &&
            complete_bytes_<=max_bytes_-(key_bytes+payload_bytes)) break;
        const auto &old=complete_cache_.front();
        complete_bytes_-=old.second.bytes;
        complete_key_bytes_-=old.second.bytes;
        if (--old.second.result->cache_references==0) {
            complete_bytes_-=old.second.result->serialized_bytes;
            complete_payload_bytes_-=old.second.result->serialized_bytes;
            --complete_unique_outputs_;
        }
        complete_index_.erase(old.first);
        complete_cache_.pop_front();
        ++complete_evictions_;
        payload_bytes=output->cache_references?0:output->serialized_bytes;
    }
    complete_cache_.emplace_back(key,CompleteEntry{output,key_bytes,input_});
    auto added=complete_cache_.end(); --added;
    try { complete_index_.emplace(key,added); }
    catch (...) { complete_cache_.pop_back(); throw; }
    if (output->cache_references++==0) {
        complete_payload_bytes_+=output->serialized_bytes;
        ++complete_unique_outputs_;
    }
    complete_key_bytes_+=key_bytes;
    complete_bytes_+=key_bytes+payload_bytes;
    complete_peak_bytes_=std::max(complete_peak_bytes_,static_cast<I>(complete_bytes_));
    complete_peak_entries_=std::max(complete_peak_entries_,static_cast<I>(complete_cache_.size()));
}

J ReplaySession::replay_reused(I n,const Word &word,const Limits &limits) {
    J initial=initial_work(n,word);
    input_=checked_add(input_,1);
    struct Frame {
        J caller;
        J unit;
        std::string key;
        J children=J::array();
        std::size_t next=0;
        OrderedStates output;
        std::shared_ptr<const CompleteOutput> sole_result;
        bool expanded=false;
    };
    std::vector<std::unique_ptr<Frame>> stack;
    std::unordered_set<std::string> active;
    I expanded=0, peak_stack=0, references=0, cache_hits=0, live_terms=0, peak_live_terms=0;
    J counts=J::object(), answer=J::array(), diagnostic=nullptr;
    std::string status="complete";
    const std::string limits_key=canonical_json(limits.json());
    auto push=[&](const J &caller) {
        check_item(caller,limits,true);
        J unit=intrinsic_state(caller);
        std::string key=limits_key+"\n"+continuation_key(unit);
        if (active.count(key)) throw ControlObligation("intrinsic completion encountered a control cycle");
        auto frame=std::make_unique<Frame>();
        frame->caller=caller;
        frame->unit=std::move(unit);
        frame->key=std::move(key);
        active.insert(frame->key);
        stack.push_back(std::move(frame));
        peak_stack=std::max(peak_stack,static_cast<I>(stack.size()));
        references=checked_add(references,1);
    };
    auto account=[&](I change) {
        live_terms=checked_add(live_terms,change);
        peak_live_terms=std::max(peak_live_terms,live_terms);
        if (live_terms>limits.max_terms) throw LimitExceeded("active completion terms limit exceeded");
    };
    auto finish=[&](const std::shared_ptr<const CompleteOutput> &local) {
        // Detach the caller before popping. Every local result is rechecked
        // after exact scalar multiplication and endpoint/context restoration.
        J caller=std::move(stack.back()->caller);
        account(-stack.back()->output.size());
        active.erase(stack.back()->key);
        stack.pop_back();
        if (stack.empty()) {
            OrderedStates root;
            for (const J &term:local->terms) collect_output(root,restore_output(term,caller),limits);
            answer=root.values();
        } else {
            auto &parent=*stack.back();
            bool identity=caller.at("contexts").empty() && caller.at("coefficient")==one() &&
                caller.at("top_labels")==caller.at("bottom_labels") &&
                length(caller.at("bottom_labels"))==as_i(caller.at("rank"));
            I port=1;
            for (const J &name:caller.at("bottom_labels")) identity=identity && as_i(name)==port++;
            if (identity && parent.children.size()==1 && parent.output.empty()) {
                for (const J &term:local->terms) check_item(term,limits,false);
                parent.sole_result=local;
            } else {
                for (const J &term:local->terms) {
                    const I before=parent.output.size();
                    collect_output(parent.output,restore_output(term,caller),limits);
                    account(checked_sub(parent.output.size(),before));
                }
            }
            account(-1);
        }
    };
    try {
        push(initial);
        while (!stack.empty()) {
            Frame &frame=*stack.back();
            if (!frame.expanded) {
                auto found=complete_index_.find(frame.key);
                if (found!=complete_index_.end()) {
                    auto entry=found->second;
                    ++complete_hits_; ++cache_hits;
                    if (entry->second.input!=input_) ++complete_cross_input_hits_;
                    auto local=entry->second.result;
                    complete_cache_.splice(complete_cache_.end(),complete_cache_,entry);
                    finish(local);
                    continue;
                }
                ++complete_misses_;
                if (expanded>=limits.max_steps) throw LimitExceeded("completion transition evaluation limit exceeded");
                check_item(frame.unit,limits,true);
                ++expanded;
                J result=algorithm_transition(frame.unit,catalog_,limits,false);
                if (result.contains("error")) {
                    status=result.at("error").get<std::string>();
                    diagnostic={{"message",result.at("message")},{"phase",frame.unit.at("phase")}};
                    break;
                }
                const std::string name=event_name(result.at("unit_evidence"));
                counts[name]=checked_add(counts.contains(name)?as_i(counts[name]):0,1);
                // Collect sibling continuations before recursively evaluating
                // them, including exact coefficient cancellations.
                OrderedStates pending;
                std::vector<const J *> children;
                for (const J &child:result.at("children")) children.push_back(&child);
                auto staged=stage(pending,children,true,{});
                commit(pending,staged.updates);
                frame.children=pending.values();
                account(length(frame.children));
                for (const J &child:frame.children) check_item(child,limits,true);
                for (J &term:result.at("output")) {
                    const I before=frame.output.size();
                    collect_output(frame.output,std::move(term),limits);
                    account(checked_sub(frame.output.size(),before));
                }
                if (checked_add(length(frame.children),frame.output.size())>limits.max_terms)
                    throw LimitExceeded("local transition term limit exceeded");
                frame.expanded=true;
            }
            if (frame.next<frame.children.size()) {
                J child=std::move(frame.children.at(frame.next++));
                push(child);
            } else {
                std::shared_ptr<const CompleteOutput> local=frame.sole_result;
                if (!local) {
                    J terms=frame.output.values();
                    const std::size_t bytes=(max_entries_ && max_bytes_)?canonical_json(terms).size():0;
                    local=std::make_shared<const CompleteOutput>(CompleteOutput{std::move(terms),bytes});
                }
                retain_complete(frame.key,local);
                finish(local);
            }
        }
    } catch (const LimitExceeded &error) {
        status="incomplete";
        diagnostic={{"message",error.what()}};
    } catch (const ControlObligation &error) {
        status="control_obligation";
        diagnostic={{"message",error.what()}};
    }
    if (status!="complete") answer=J::array();
    return {{"status",status},{"diagnostic",diagnostic},{"output",answer},
            {"metrics",{{"expanded_nodes",expanded},{"node_references",references},
                {"completed_cache_hits",cache_hits},{"peak_stack",peak_stack},{"peak_live_terms",peak_live_terms},
                {"output_terms",length(answer)},{"rule_applications",counts}}},
            {"cache",cache_info()}};
}

J ReplaySession::replay(I n, const Word &word, const Limits &limits, I batch_size, bool trace) {
    if (n<1 || batch_size<1 || !word.is_array()) throw ReplayError("positive integer rank and batch size required");
    for (const auto &token : word) if (!is_int(token) || as_i(token)==0 || abs_i(as_i(token))>=n)
        throw ReplayError("a literal signed-generator list in the stated rank is required");
    input_=checked_add(input_,1);
    auto initial=initial_work(n,word);
    OrderedStates frontier, output;
    frontier.set(continuation_key(initial),initial);
    J transitions=J::array(), counts=J::object();
    J metrics={{"steps",0},{"attempted_steps",0},{"batches",0},{"max_terms",1},
               {"max_local_word_tokens",checked_add(length(word),1)},{"max_coefficient_support",1},
               {"max_integer_bits",1},{"q_exponent_min",0},{"q_exponent_max",0},
               {"maximum_depth",0},{"max_branch_count",0},{"merged_terms",0},{"zero_terms",0}};
    std::string status="complete"; J diagnostic=nullptr;
    while (!frontier.empty()) {
        auto selected=frontier.prefix(batch_size);
        J before=selected.front().second, results=J::array();
        Staged pending, completed;
        I terms=0;
        std::unordered_set<std::string> removed;
        try {
            if (checked_add(as_i(metrics["steps"]),static_cast<I>(selected.size()))>limits.max_steps)
                throw LimitExceeded("step limit would be exceeded by the next atomic batch");
            for (const auto &item : selected) check_item(item.second,limits,true);
            metrics["attempted_steps"]=checked_add(as_i(metrics["attempted_steps"]),static_cast<I>(selected.size()));
            for (const auto &item : selected) results.push_back(evaluate(item.second,limits,trace));
            bool failed=false;
            for (std::size_t i=0;i<selected.size();++i) if (results[i].contains("error")) {
                before=selected[i].second;
                status=results[i].at("error").get<std::string>();
                diagnostic={{"message",results[i].at("message")},{"phase",before.at("phase")},
                            {"state",before},{"batch_rolled_back",true}};
                failed=true; break;
            }
            if (failed) break;
            std::vector<const J *> children, added;
            for (const auto &item : selected) removed.insert(item.first);
            for (const auto &result : results) {
                for (const auto &child:result.at("children")) children.push_back(&child);
                for (const auto &term:result.at("output")) added.push_back(&term);
            }
            pending=stage(frontier,children,true,removed);
            completed=stage(output,added,false,{});
            terms=checked_add(count_after(frontier,pending.updates,removed),count_after(output,completed.updates,{}));
            if (terms>limits.max_terms) throw LimitExceeded("collected term limit exceeded");
            for (const auto &item : pending.updates.items()) if (!item.second.is_null()) check_item(item.second,limits,true);
            for (const auto &item : completed.updates.items()) if (!item.second.is_null()) check_item(item.second,limits,false);
        } catch (const LimitExceeded &error) {
            status="incomplete";
            diagnostic={{"message",error.what()},{"phase","batch"},{"batch_representative",before},{"batch_rolled_back",true}};
            break;
        }
        metrics["steps"]=checked_add(as_i(metrics["steps"]),static_cast<I>(selected.size()));
        metrics["batches"]=checked_add(as_i(metrics["batches"]),1);
        metrics["merged_terms"]=checked_add(as_i(metrics["merged_terms"]),checked_add(pending.merged,completed.merged));
        metrics["zero_terms"]=checked_add(as_i(metrics["zero_terms"]),checked_add(pending.zeroes,completed.zeroes));
        update_max(metrics,"max_terms",terms);
        for (std::size_t i=0;i<selected.size();++i) {
            auto name=event_name(results[i].at("unit_evidence"));
            counts[name]=checked_add(counts.contains(name)?as_i(counts[name]):0,1);
            update_max(metrics,"maximum_depth",checked_add(as_i(selected[i].second.at("depth")),1));
            update_max(metrics,"max_branch_count",checked_add(length(results[i].at("children")),length(results[i].at("output"))));
        }
        for (const auto &item : pending.updates.items()) if (!item.second.is_null()) update_item_metrics(metrics,item.second,true);
        for (const auto &item : completed.updates.items()) if (!item.second.is_null()) update_item_metrics(metrics,item.second,false);
        if (trace) {
            J predecessors=J::array(); for (const auto &item : selected) predecessors.push_back(item.second);
            transitions.push_back({{"batch",checked_sub(as_i(metrics["batches"]),1)},{"predecessors",predecessors},
                {"rewrites",results},{"pending_updates",pending.updates.object()},{"output_updates",completed.updates.object()}});
        }
        commit(frontier,pending.updates,removed); commit(output,completed.updates);
    }
    metrics["rule_applications"]=counts; metrics["output_terms"]=output.size(); metrics["pending_terms"]=frontier.size();
    J certificate={{"schema","zipper-collected-run-v1"},{"scope","literal_word_collected_spanning_procedure"},
            {"n",n},{"input_word",word},{"coefficient_ring","Z[q,q^-1]"},
            {"source_sha256",reference_sources()},{"catalog_sha256",catalog_.sha256},
            {"trace_enabled",trace},{"limits",limits.json()},
            {"execution",{{"batch_size",batch_size},{"schedule","stable-prefix-atomic-v1"},
                          {"collection","full-continuation-except-scalar-and-depth-v1"},
                          {"local_term_floor",4},{"local_support_floor",2}}},
            {"status",status}};
    certificate["diagnostic"]=std::move(diagnostic);
    certificate["pending"]=frontier.values();
    certificate["output"]=output.values();
    certificate["transitions"]=std::move(transitions);
    certificate["metrics"]=std::move(metrics);
    return certificate;
}
} // namespace zipper_verify
