// Independent JSON-state turnaround replay, derived from the Python specification.
#include "transport.hpp"
#include "catalog.hpp"
#include "controllers.hpp"
#include <algorithm>
#include <set>

namespace zipper_verify {
namespace {
bool fields(const J &value, const std::set<std::string> &names) {
    if (!value.is_object() || value.size() != names.size()) return false;
    for (const auto &name : names) if (!value.contains(name)) return false;
    return true;
}
bool pure(const Word &word) {
    return word.is_array() && std::all_of(word.begin(), word.end(), [](const J &x){ return is_int(x); });
}
bool suffix(const std::string &text, const std::string &ending) {
    return text.size() >= ending.size() && text.compare(text.size()-ending.size(), ending.size(), ending)==0;
}
Word packet_word(const J &terms, const J &env, I budget) {
    const auto polynomial = expand_polynomial(terms, env, J::object(), false, budget);
    if (polynomial.size()!=1) throw TurnaroundError("a monomial packet word is required");
    const auto &term=polynomial.front();
    if (term.q!=0 || term.coefficient!=1) throw TurnaroundError("packet word must have coefficient one");
    return term.word;
}
J typed_event(const J &event, I n, I i, Catalog &catalog, const Limits &limits) {
    if (!fields(event,{"rule_id","parameters","rhs_index","mirror_policy"}) ||
        !event["rule_id"].is_string() ||
        !fields(event["parameters"],{"n","i","j","k","l","eps"}) ||
        !is_int(event["rhs_index"]) ||
        (event["mirror_policy"]!="M" && event["mirror_policy"]!="bar_mirror"))
        throw TurnaroundError("invalid processed-tail event");
    const J &env=event["parameters"];
    for (const auto &value:env) if (!is_int(value)) throw TurnaroundError("invalid processed-tail event");
    const auto name=event["rule_id"].get<std::string>();
    const bool mirrored=suffix(name,".M") || suffix(name,".BM");
    const bool barred=suffix(name,".BM");
    if (as_i(env["n"])!=n || as_i(env["i"])!=i || name.rfind("zip.",0)!=0 ||
        !(suffix(name,".I") || mirrored) ||
        (barred && event["mirror_policy"]!="bar_mirror") ||
        (suffix(name,".M") && event["mirror_policy"]!="M"))
        throw TurnaroundError("tail event has incompatible rank or orientation");
    const J rule=rule_for_event(catalog,event);
    const I index=as_i(event["rhs_index"]);
    if (!rule.contains("execution") || !rule["execution"].is_object() ||
        !rule["execution"].contains("algebra_verified") ||
        !rule["execution"]["algebra_verified"].is_boolean() ||
        !rule["execution"]["algebra_verified"].get<bool>() ||
        !condition(rule["condition"],env) || index<0 || index>=length(rule["rhs"]))
        throw TurnaroundError("tail event lacks a valid accepted rule branch");
    const I old_j=mirrored ? checked_add(checked_sub(checked_sub(n,i),as_i(env["j"])),2) : as_i(env["j"]);
    const J base=catalog.get(rule["base_id"].get<std::string>()+".I");
    auto [jj,kk,mm,aa,bb]=branch_shape(base["rhs"][index],env,mirrored);
    if (barred) {aa=neg(aa);bb=neg(bb);}
    const auto [output,inp]=anchors(i,jj,kk,mm);
    const J &term=rule["rhs"][index];
    const Word local=packet_word(J::array({J{{"coefficient",J::array({J{{"integer",1},{"q_power",0}}})},{"word",term["word"]}}}),env,limits.expansion_budget);
    const Word prefix=concat({frame_power(checked_sub(i,1),output,aa,limits.max_word_size),
        framed_core(n,i,jj,kk,mm,event["mirror_policy"].get<std::string>()),
        frame_power(checked_sub(i,1),inp,bb,limits.max_word_size)});
    if (slice(local,0,length(prefix))!=prefix) throw TurnaroundError("tail event does not have the catalog frame shape");
    return J{{"env",env},{"pre_j",old_j},{"pre_mirrored",mirrored},
        {"post_j",jj},{"post_mirrored",mm},{"output_delta",aa},{"input_delta",bb},
        {"tail",slice(local,length(prefix))},{"mirror_policy",event["mirror_policy"]}};
}
J turnaround_certificate(const J &before,const J &after,I n,Catalog &catalog,
                        const Limits &limits,bool trace,const J &events,I event_count) {
    return J{{"schema","zipper-turnaround-v1"},{"n",n},{"trace_enabled",trace},
        {"source_sha256",reference_sources().at("turnaround.py")},
        {"dependency_sha256",source_subset({"framed.py","frames.py","transport.py","catalog.py","upward.py","semilinear.py"})},
        {"catalog_sha256",catalog.sha256},{"initial",before},{"final",after},{"events",events},
        {"summary",J{{"events",event_count},{"word_tokens_before",length(before["word"])},{"word_tokens_after",length(after["word"])}}},
        {"limits",limits.json()}};
}
} // namespace

TurnParts turnaround_parts(const J &value,I n,I budget) {
    validate_branch(value,n);
    const J &c=value.at("control");
    const std::set<std::string> names={"phase","i","j","mirrored","mirror_policy",
        "output_power","input_power","exterior_length","payload_length","routing_length",
        "transported_length","processed_length","tail_length"};
    if (!fields(c,names) || c["phase"]!="turnaround_complete" || !c["mirrored"].is_boolean() ||
        (c["mirror_policy"]!="M" && c["mirror_policy"]!="bar_mirror"))
        throw TurnaroundError("invalid turnaround output boundaries");
    for (const auto &key:names)
        if (key!="phase" && key!="mirrored" && key!="mirror_policy" && !is_int(c[key]))
            throw TurnaroundError("invalid turnaround output boundaries");
    const I i=as_i(c["i"]),j=as_i(c["j"]);
    if (i<2 || i>n || j<1 || j>checked_add(checked_sub(n,i),1))
        throw TurnaroundError("invalid turnaround output boundaries");
    for (const auto &key:{"exterior_length","payload_length","routing_length","transported_length","processed_length","tail_length"})
        if (as_i(c[key])<0) throw TurnaroundError("invalid turnaround output boundaries");
    if (as_i(c["exterior_length"])!=checked_add(checked_add(as_i(c["payload_length"]),as_i(c["routing_length"])),as_i(c["transported_length"])))
        throw TurnaroundError("invalid turnaround output boundaries");
    const I left=as_i(c["exterior_length"]),plen=as_i(c["processed_length"]),tlen=as_i(c["tail_length"]);
    const Word &word=value["word"];
    const I end=checked_sub(checked_sub(length(word),tlen),plen);
    const Word exterior=slice(word,0,left),processed=slice(word,end,checked_add(end,plen)),tail=slice(word,checked_add(end,plen));
    const auto [output,inp]=anchors(i,j,0,c["mirrored"].get<bool>());
    const Word core=concat({frame_power(checked_sub(i,1),output,as_i(c["output_power"]),budget),
        J::array({slider(i,j,c["mirrored"].get<bool>())}),
        frame_power(checked_sub(i,1),inp,as_i(c["input_power"]),budget)});
    if (word!=concat({exterior,core,processed,tail}) || !pure(concat({exterior,processed,tail})))
        throw TurnaroundError("turnaround output does not match its full expression");
    return {c,exterior,processed,tail};
}

std::pair<J,J> turnaround(const J &initial,I n,Catalog &catalog,const Limits &limits,
                         bool trace,bool want_certificate) {
    limits.check(run_state(n,J::array({initial})));
    const J &start_control=initial["control"];
    if (start_control.contains("i") && start_control["i"]==2) {
        const auto [c,u,p,tail]=upward_parts(initial,n);
        if (!u.empty()) throw TurnaroundError("completed upward pass required");
        const J band=transport_band(p,n,2,as_i(c["j"]));
        check_transport(band);
        const Word exterior=concat({upward_friends(2,as_i(c["j"]),as_i(c["k"])),p});
        const I end=as_i(band["end"]);
        const J cc={{"phase","turnaround_complete"},{"i",2},{"j",end},{"mirrored",false},
            {"mirror_policy","M"},{"output_power",0},{"input_power",0},
            {"exterior_length",length(exterior)},{"payload_length",0},{"routing_length",c["k"]},
            {"transported_length",length(p)},{"processed_length",0},{"tail_length",length(tail)}};
        const J following=branch(concat({exterior,J::array({slider(2,end)}),tail}),initial["coefficient"],cc);
        if (!(trace || want_certificate)) return {following,J::object()};
        const J events=trace ? J::array({J{{"band_transport",band}}}) : J::array();
        return {following,turnaround_certificate(initial,following,n,catalog,limits,trace,events,length(band["steps"]))};
    }
    auto [c,u,p,tail]=framed_parts(initial,n,limits.max_word_size);
    if (!u.empty() || c["phase"]!="framed_upper_pass_complete") throw TurnaroundError("completed framed upward pass required");
    const I i=as_i(c["i"]);
    I j=as_i(c["j"]),a=as_i(c["output_power"]),b=as_i(c["input_power"]);
    bool mirror=c["mirrored"].get<bool>();
    const std::string policy=c["mirror_policy"].get<std::string>();
    const Word payload=c.value("output_payload",J::array());
    const Word friends=framed_core(n,i,j,as_i(c["k"]),mirror,policy);
    Word exterior=concat({payload,slice(friends,0,checked_sub(length(friends),1))});
    const I payload_length=length(payload),routing_length=checked_mul(checked_sub(i,1),as_i(c["k"]));
    const J &processed_events=c["processed_events"];
    if (length(processed_events)>limits.max_steps) throw LimitExceeded("turnaround event limit exceeded");
    J evidence=J::array();
    for (const auto &event:processed_events) {
        const J typed=typed_event(event,n,i,catalog,limits);
        const Word T=typed["tail"];
        if (as_i(typed["post_j"])!=j || typed["post_mirrored"]!=mirror || typed["mirror_policy"]!=policy ||
            slice(p,0,length(T))!=T)
            throw TurnaroundError("tail event disagrees with current slider or word");
        p=slice(p,length(T));
        const I aa=as_i(typed["output_delta"]),bb=as_i(typed["input_delta"]);
        const bool pre_mirrored=typed["pre_mirrored"].get<bool>();
        J packet_id=nullptr;
        Word W=J::array(),R;
        if (aa!=0 || bb!=0) {
            if ((aa!=-1 && aa!=1) || bb!=neg(aa) || pre_mirrored==mirror)
                throw TurnaroundError("invalid emitted flype frame increments");
            const bool barred=pre_mirrored && policy=="bar_mirror";
            const std::string sign=aa==(barred?1:-1) ? "positive" : "negative";
            const std::string variant=barred?".BM":pre_mirrored?".M":".I";
            packet_id="zip.flype."+sign+variant;
            J env=typed["env"];
            env["j"]=pre_mirrored?checked_add(checked_sub(checked_sub(n,i),j),2):j;
            env["k"]=0;
            const J packet=rule_for_event(catalog,J{{"rule_id",packet_id},{"parameters",env},{"rhs_index",0},{"mirror_policy",policy}});
            const Word incoming=packet_word(packet["lhs"],env,limits.expansion_budget);
            if (slice(incoming,checked_sub(length(incoming),1))!=J::array({slider(i,j,pre_mirrored)}))
                throw TurnaroundError("flype packet has an invalid input slider");
            W=slice(incoming,0,checked_sub(length(incoming),1));
            const Word outgoing=packet_word(packet["rhs"],env,limits.expansion_budget);
            const auto [output,inp]=anchors(i,j,0,mirror);
            const Word packet_core=concat({frame_power(checked_sub(i,1),output,aa,limits.max_word_size),
                J::array({slider(i,j,mirror)}),frame_power(checked_sub(i,1),inp,bb,limits.max_word_size)});
            if (slice(outgoing,0,length(packet_core))!=packet_core)
                throw TurnaroundError("flype packet has an invalid retained frame");
            const Word C=slice(outgoing,length(packet_core));
            if (length(C)!=checked_sub(i,1) || slice(T,0,length(C))!=C)
                throw TurnaroundError("tail is missing its complete outgoing flype tape");
            R=slice(T,length(C));
            a=checked_sub(a,aa);b=checked_sub(b,bb);mirror=pre_mirrored;
        } else {
            if (pre_mirrored!=mirror) throw TurnaroundError("orientation change requires an explicit flype packet");
            R=T;
        }
        const J band=transport_band(R,n,i,j,mirror);
        check_transport(band);
        if (band["end"]!=typed["pre_j"])
            throw TurnaroundError("residual tail does not return to the previous band anchor");
        append(exterior,W);append(exterior,R);j=as_i(band["end"]);
        if (trace) evidence.push_back(J{{"event",event},{"packet_rule",packet_id},{"band_transport",band},
            {"output_power",a},{"input_power",b},{"j",j},{"mirrored",mirror}});
    }
    const J &original=c["initial_slider"];
    if (j!=as_i(original["j"]) || mirror!=original["mirrored"] || a!=as_i(original["output_power"]) ||
        b!=as_i(original["input_power"]) || length(p)!=as_i(original["processed_length"]) || policy!=original["mirror_policy"])
        throw TurnaroundError("typed tails do not restore the initial slider data");
    const auto [output,inp]=anchors(i,j,0,mirror);
    const Word word=concat({exterior,frame_power(checked_sub(i,1),output,a,limits.max_word_size),
        J::array({slider(i,j,mirror)}),frame_power(checked_sub(i,1),inp,b,limits.max_word_size),p,tail});
    const J cc={{"phase","turnaround_complete"},{"i",i},{"j",j},{"mirrored",mirror},{"mirror_policy",policy},
        {"output_power",a},{"input_power",b},{"exterior_length",length(exterior)},
        {"payload_length",payload_length},{"routing_length",routing_length},
        {"transported_length",checked_sub(checked_sub(length(exterior),payload_length),routing_length)},
        {"processed_length",length(p)},{"tail_length",length(tail)}};
    const J following=branch(word,initial["coefficient"],cc);
    limits.check(run_state(n,J::array({following})));
    turnaround_parts(following,n,limits.max_word_size);
    if (!(trace || want_certificate)) return {following,J::object()};
    return {following,turnaround_certificate(initial,following,n,catalog,limits,trace,evidence,length(processed_events))};
}
} // namespace zipper_verify
