// Independent JSON-state framed controller, derived from the Python specification.
#include "controllers.hpp"
#include "transport.hpp"
#include <algorithm>
#include <set>

namespace zipper_verify {
namespace {
bool fields(const J &value,std::initializer_list<const char*> names) {
    if(!value.is_object() || value.size()!=names.size()) return false;
    for(auto name:names) if(!value.contains(name)) return false;
    return true;
}
bool pure(const J &word) {
    if(!word.is_array()) return false;
    return std::all_of(word.begin(),word.end(),[](const J &x){return is_int(x);});
}
I field(const J &object,const char *name) {return as_i(object.at(name));}
J replace_branches(const J &state,I index,const J &replacement,I remove=1) {
    J result=state;
    result["branches"]=concat({slice(state.at("branches"),0,index),replacement,
                               slice(state.at("branches"),checked_add(index,remove))});
    return result;
}
J changes(const J &predecessor,const J &successors) {
    I total_words=0,total_support=0,max_word=0,bits=0;
    I qmin=std::numeric_limits<I>::max(),qmax=std::numeric_limits<I>::min();
    for(const auto &b:successors) {
        total_words=checked_add(total_words,length(b.at("word")));
        total_support=checked_add(total_support,length(b.at("coefficient")));
        max_word=std::max(max_word,length(b.at("word")));
        for(const auto &[q,c]:coefficient(b.at("coefficient"))) {
            qmin=std::min(qmin,q);qmax=std::max(qmax,q);bits=std::max(bits,integer_bits(c));
        }
    }
    if(successors.empty()) throw ControlError("framed rule has no successors");
    return {{"word_tokens",checked_sub(total_words,length(predecessor.at("word")))},
            {"coefficient_support",checked_sub(total_support,length(predecessor.at("coefficient")))},
            {"max_word_tokens",max_word},{"q_exponent_min",qmin},{"q_exponent_max",qmax},
            {"max_integer_bits",bits}};
}
}

std::pair<I,I> anchors(I i,I j,I k,bool mirrored) {
    (void)i;
    return mirrored?std::make_pair(checked_sub(j,k),checked_add(j,1)):
                    std::make_pair(checked_add(checked_add(j,k),1),j);
}
Word frame_power(I size,I start,I exponent,I budget) {
    Big expanded=Big(size)*Big(checked_sub(size,1));expanded/=2;expanded*=Big(abs_i(exponent));
    if(expanded>Big(budget)) throw LimitExceeded("expanded frame word limit exceeded");
    if(expanded==0) return J::array();
    Word half=J::array();
    for(I width=checked_sub(size,1);width>0;--width)
        for(I off=0;off<width;++off) half.push_back(checked_add(start,off));
    if(exponent<0) half=inverse(half);
    Word result=J::array();
    for(I copy=0;copy<abs_i(exponent);++copy) append(result,half);
    return result;
}
Word framed_core(I n,I i,I j,I k,bool mirrored,const std::string &policy) {
    I base=mirrored?checked_add(checked_sub(checked_sub(n,i),j),2):j;
    Word word=upward_friends(i,base,k);
    if(mirrored) for(auto &x:word) {
        I value=checked_sub(n,as_i(x));
        x=policy=="bar_mirror"?neg(value):value;
    }
    word.push_back(slider(i,j,mirrored,false,false));
    return word;
}
Word friend_word(I n,const J &control) {
    Word full=framed_core(n,field(control,"i"),field(control,"j"),field(control,"k"),
                          control.at("mirrored").get<bool>(),control.at("mirror_policy").get<std::string>());
    return slice(full,0,checked_sub(length(full),1));
}
EventRule held_rule_for_event(Catalog &catalog,const J &event) {
    std::string rid=event.at("rule_id").get<std::string>();
    if(rid.size()>=3 && rid.compare(rid.size()-3,3,".BM")==0)
        return {catalog.semilinear_rule_handle(rid.substr(0,rid.size()-3)+".M",event.at("parameters")),nullptr};
    return {nullptr,&catalog.get(rid)};
}
J rule_for_event(Catalog &catalog,const J &event) {
    return held_rule_for_event(catalog,event).get();
}
J framed_control(I i,I j,I k,bool mirrored,I a,I b,I left,I processed,I tail) {
    return {{"phase",left?"framed_zip_up":"framed_upper_pass_complete"},
            {"i",i},{"j",j},{"k",k},{"mirrored",mirrored},{"output_power",a},{"input_power",b},
            {"left_length",left},{"processed_length",processed},{"tail_length",tail},{"output_payload",J::array()}};
}
Word framed_expression(I n,const J &c,const Word &u,const Word &p,const Word &tail,I budget) {
    I i=field(c,"i"),j=field(c,"j"),k=field(c,"k");
    I a=field(c,"output_power"),b=field(c,"input_power");
    const J &payload=c.at("output_payload");
    Big total=Big(length(u))+Big(length(payload))+Big(length(p))+Big(length(tail));
    total+=Big(checked_sub(i,1))*Big(k)+1;
    Big half=Big(checked_sub(i,1))*Big(checked_sub(i,2));half/=2;
    total+=half*(Big(abs_i(a))+Big(abs_i(b)));
    if(total>Big(budget)) throw LimitExceeded("expanded framed state word limit exceeded");
    auto [out,inp]=anchors(i,j,k,c.at("mirrored").get<bool>());
    return concat({u,payload,frame_power(i-1,out,a,budget),
                   framed_core(n,i,j,k,c.at("mirrored").get<bool>(),c.at("mirror_policy").get<std::string>()),
                   frame_power(i-1,inp,b,budget),p,tail});
}
J framed_initialize(I n,const Word &left,I size,I start,I friends,bool mirrored,I output_power,
                    I input_power,const Word &processed,const Word &tail,const Word &payload,
                    const std::string &policy,const Limits &limits) {
    if((policy!="M" && policy!="bar_mirror") || size<2 || size>n || start<1 || friends<0 ||
       checked_sub(checked_add(start,size),1)>n ||
       (mirrored?checked_sub(start,friends)<1:checked_sub(checked_add(checked_add(start,size),friends),1)>n))
        throw ControlError("valid framed slider, friends and integer powers required");
    if(size==2) {
        if(mirrored || !payload.empty()) throw ControlError("use the normal X2 convention for the two-strand delegate");
        return upward_initialize(n,left,size,start,friends,processed,tail);
    }
    if(!pure(left)||!pure(processed)||!pure(tail))
        throw ControlError("framed exterior must contain signed braid integers");
    auto [out,unused]=anchors(size,start,friends,mirrored);(void)unused;
    if(!payload.is_array()) throw ControlError("protected output payload must be supported on the marked tape");
    for(const auto &x:payload)
        if(!is_int(x)||abs_i(as_i(x))<out || abs_i(as_i(x))>=checked_sub(checked_add(out,size),2))
            throw ControlError("protected output payload must be supported on the marked tape");
    J c=framed_control(size,start,friends,mirrored,output_power,input_power,length(left),length(processed),length(tail));
    c["mirror_policy"]=policy;c["output_payload"]=payload;
    c["initial_slider"]={{"j",start},{"k",friends},{"mirrored",mirrored},{"output_power",output_power},
                         {"input_power",input_power},{"processed_length",length(processed)},{"mirror_policy",policy}};
    c["processed_events"]=J::array();
    J result=run_state(n,J::array({branch(framed_expression(n,c,left,processed,tail,limits.max_word_size),one(),c)}));
    limits.check(result);return result;
}
std::tuple<J,Word,Word,Word> framed_parts(const J &b,I n,I budget) {
    const J &c=b.at("control");
    if(!fields(c,{"phase","i","j","k","mirrored","output_power","input_power","left_length",
                  "processed_length","tail_length","initial_slider","processed_events","output_payload","mirror_policy"})
       ||!c.at("mirrored").is_boolean() || (c.at("mirror_policy")!="M" && c.at("mirror_policy")!="bar_mirror"))
        throw ControlError("invalid framed state boundaries");
    for(auto key:{"i","j","k","output_power","input_power","left_length","processed_length","tail_length"})
        if(!is_int(c.at(key))) throw ControlError("invalid framed state boundaries");
    I i=field(c,"i"),j=field(c,"j"),k=field(c,"k");bool mirror=c.at("mirrored").get<bool>();
    I a=field(c,"left_length"),bsize=field(c,"processed_length"),d=field(c,"tail_length");
    if(i<3||i>n||j<1||k<0||checked_sub(checked_add(j,i),1)>n ||
       (mirror?checked_sub(j,k)<1:checked_sub(checked_add(checked_add(j,i),k),1)>n) ||
       std::min({a,bsize,d})<0) throw ControlError("invalid framed state boundaries");
    auto [out,unused]=anchors(i,j,k,mirror);(void)unused;
    const auto &payload=c.at("output_payload");
    if(!payload.is_array()) throw ControlError("invalid protected output-tape payload");
    for(const auto &x:payload)
        if(!is_int(x)||abs_i(as_i(x))<out||abs_i(as_i(x))>=checked_sub(checked_add(out,i),2))
            throw ControlError("invalid protected output-tape payload");
    const auto &initial=c.at("initial_slider");
    if(!fields(initial,{"j","k","mirrored","output_power","input_power","processed_length","mirror_policy"}) ||
       !initial.at("mirrored").is_boolean() || initial.at("mirror_policy")!=c.at("mirror_policy") ||
       !c.at("processed_events").is_array()) throw ControlError("invalid framed input-band record");
    for(auto key:{"j","k","output_power","input_power","processed_length"})
        if(!is_int(initial.at(key))) throw ControlError("invalid framed input-band record");
    for(const auto &event:c.at("processed_events")) {
        if(!fields(event,{"rule_id","parameters","rhs_index","mirror_policy"}) ||
           !event.at("rule_id").is_string() || !event.at("parameters").is_object() ||
           event.at("mirror_policy")!=c.at("mirror_policy") || !is_int(event.at("rhs_index")) || field(event,"rhs_index")<0)
            throw ControlError("invalid processed-band event");
        for(const auto &x:event.at("parameters")) if(!is_int(x)) throw ControlError("invalid processed-band event");
    }
    if(c.at("phase")!=(a?"framed_zip_up":"framed_upper_pass_complete")) throw ControlError("invalid framed phase");
    const auto &word=b.at("word");
    I p_start=checked_sub(checked_sub(length(word),bsize),d);
    Word u=slice(word,0,a),p=slice(word,p_start,checked_add(p_start,bsize)),t=slice(word,checked_add(p_start,bsize));
    if(!pure(u)||!pure(p)||!pure(t)) throw ControlError("framed exterior must be a pure braid word");
    if(word!=framed_expression(n,c,u,p,t,budget)) throw ControlError("framed boundaries do not describe the complete expression");
    return {c,u,p,t};
}
std::tuple<I,I,bool,I,I> branch_shape(const J &base_term,const J &env,bool mirrored) {
    auto nodes=upward_flatten(base_term.at("word"));
    std::size_t at=0;while(at<nodes.size() && nodes[at].at("type")!="slider") ++at;
    if(at==nodes.size()) throw ControlError("catalog branch has no slider");
    const J &sl=nodes[at];I j=integer(sl.at("start"),env);
    bool base_mirror=sl.value("orientation",std::string("normal"))=="mirrored";
    I k=0;if(at && nodes[at-1].at("type")=="block") k=integer(nodes[at-1].at("count"),env);
    auto exponent=[](const J &node)->I {
        if(node.at("type")=="half_twist") return 1;
        if(node.at("type")=="inverse" && node.at("word").at("type")=="half_twist") return -1;
        return 0;
    };
    I aa=at?exponent(nodes[0]):0,bb=at+1<nodes.size()?exponent(nodes[at+1]):0;
    return {mirrored?checked_add(checked_sub(checked_sub(field(env,"n"),field(env,"i")),j),2):j,
            k,base_mirror!=mirrored,aa,bb};
}

std::pair<J,J> framed_step(const J &state,I selected,Catalog &catalog,const Limits &limits,bool trace) {
    try {
        limits.check(state);
        const auto &branches=state.at("branches");I n=field(state,"n");
        if(selected<0||selected>=length(branches)) throw ControlError("selected framed branch does not exist");
        const J &old=branches.at(static_cast<std::size_t>(selected));
        auto [c,u,p,tail]=framed_parts(old,n,limits.max_word_size);
        if(u.empty()) throw ControlError("framed upward pass is already complete");
        I i=field(c,"i"),j=field(c,"j"),k=field(c,"k");bool mirror=c.at("mirrored").get<bool>();
        auto [out,inp]=anchors(i,j,k,mirror);std::string policy=c.at("mirror_policy").get<std::string>();
        bool bm=mirror&&policy=="bar_mirror";I last=as_i(u.back());
        if(out<=abs_i(last) && abs_i(last)<checked_sub(checked_add(out,i),2)) {
            J cc=c;cc["left_length"]=checked_sub(length(u),1);
            cc["phase"]=as_i(cc["left_length"])?"framed_zip_up":"framed_upper_pass_complete";
            cc["output_payload"]=concat({J::array({u.back()}),c.at("output_payload")});
            J new_branch=branch(old.at("word"),old.at("coefficient"),cc);
            J following=replace_branches(state,selected,J::array({new_branch}));
            framed_parts(new_branch,n,limits.max_word_size);
            J event={{"rule_id","control.absorb_internal_tape"},{"metric_delta",changes(old,J::array({new_branch}))}};
            if(trace) {event["equality_reference"]="identical_literal_word";event["successor_state"]=following;}
            return {following,event};
        }
        J env={{"n",n},{"i",i},{"j",mirror?checked_add(checked_sub(checked_sub(n,i),j),2):j},
               {"k",k},{"l",mirror?checked_sub(n,abs_i(last)):abs_i(last)},
               {"eps",(last>0?1:-1)*(bm?-1:1)}};
        Word core=framed_core(n,i,j,k,mirror,policy);
        struct Candidate {
            std::shared_ptr<const J> derived;
            const J *ordinary;
            Word w;
            FrameWitness witness;
            const J &rule() const {return ordinary?*ordinary:*derived;}
        };
        std::vector<Candidate> candidates;
        for(const auto &id:catalog.select_rule_ids(env,"zip.",mirror?".M":".I")) {
            std::shared_ptr<const J> derived;
            const J &rule=bm?*(derived=catalog.semilinear_rule_handle(id,env)):catalog.get(id);
            auto lhs=expand_polynomial(rule.at("lhs"),env,J::object(),false,limits.expansion_budget);
            if(lhs.empty()) throw ControlError("catalog left side does not retain the framed core");
            const Word &occurrence=lhs.front().word;I boundary=checked_sub(length(occurrence),length(core));
            if(slice(occurrence,std::max<I>(0,boundary))!=core)
                throw ControlError("catalog left side does not retain the framed core");
            Word w=slice(occurrence,0,boundary);
            if(w.empty()||length(w)>length(u)) continue;
            try {
                auto witness=typed_transport_frame(slice(u,checked_sub(length(u),length(w))),n,i-1,out,field(c,"output_power"),"left");
                check_frame(witness);
                if(witness.transformed==w)
                    candidates.push_back({std::move(derived),bm?nullptr:&rule,w,std::move(witness)});
            } catch(const FrameError &) {}
        }
        if(candidates.size()!=1) throw ControlError(candidates.empty()?"missing framed full-tape rule":"competing framed rules");
        const auto &chosen=candidates.front();const J &rule=chosen.rule();const Word &w=chosen.w;const auto &before=chosen.witness;
        Word urest=slice(u,0,checked_sub(length(u),length(w)));
        I end=before.end,shift=checked_sub(end,out);
        Word payload=J::array();for(const auto &x:c.at("output_payload")) {
            I value=as_i(x);I idx=checked_add(abs_i(value),shift);payload.push_back(value>0?idx:neg(idx));
        }
        Word leading=frame_power(i-1,end,field(c,"output_power"),limits.max_word_size);
        Word middle_word=concat({urest,payload,leading,w,core,frame_power(i-1,inp,field(c,"input_power"),limits.max_word_size),p,tail});
        J middle=replace_branches(state,selected,J::array({branch(middle_word,old.at("coefficient"),old.at("control"))}));
        J req=request(selected,checked_add(checked_add(length(urest),length(payload)),length(leading)),
                      rule.at("id").get<std::string>(),env,J::object(),
                      {{"controller","framed_upward"},{"marked_tape_transport",true}});
        J application=apply_rule(middle,req,catalog,limits,trace);J raw=application.at("state");
        const J &base=catalog.get(rule.at("base_id").get<std::string>()+".I");
        if(rule.at("rhs").size()!=base.at("rhs").size()) throw ControlError("catalog branch counts differ");
        J new_branches=J::array(),witnesses=J::array();
        for(std::size_t idx=0;idx<rule.at("rhs").size();++idx) {
            const auto &term=rule.at("rhs").at(idx),&base_term=base.at("rhs").at(idx);
            auto [jj,kk,mm,aa,bb]=branch_shape(base_term,env,mirror);if(bm){aa=neg(aa);bb=neg(bb);}
            auto [out2,inp2]=anchors(i,jj,kk,mm);
            if(end!=out2) throw ControlError("catalog output tape disagrees with frame transport");
            Word local=upward_literal({term.at("word")},env,limits.expansion_budget);
            Word prefix=concat({frame_power(i-1,out2,aa,limits.max_word_size),framed_core(n,i,jj,kk,mm,policy)});
            Word right=frame_power(i-1,inp2,bb,limits.max_word_size);I skip=checked_add(length(prefix),length(right));
            if(slice(local,0,skip)!=concat({prefix,right})) throw ControlError("catalog branch is missing its retained flype frames");
            Word residual=slice(local,skip);
            auto after=typed_transport_frame(residual,n,i-1,inp,field(c,"input_power"),"left");check_frame(after);
            if(after.end!=inp2||after.transformed!=residual)
                throw ControlError("catalog tail does not preserve the marked input tape");
            J cc=framed_control(i,jj,kk,mm,checked_add(field(c,"output_power"),aa),checked_add(field(c,"input_power"),bb),
                                length(urest),checked_add(length(residual),length(p)),length(tail));
            cc["initial_slider"]=c.at("initial_slider");cc["mirror_policy"]=policy;
            cc["processed_events"]=concat({J::array({J{{"rule_id",rule.at("id")},{"parameters",env},{"rhs_index",idx},{"mirror_policy",policy}}}),c.at("processed_events")});
            cc["output_payload"]=payload;
            const J &emitted=raw.at("branches").at(static_cast<std::size_t>(selected)+idx);
            Word word=framed_expression(n,cc,urest,concat({residual,p}),tail,limits.max_word_size);
            new_branches.push_back(branch(word,emitted.at("coefficient"),cc));
            if(trace) witnesses.push_back({{"tail_transport",after.json()},{"emitted_output_power",aa},{"emitted_input_power",bb}});
        }
        J following=replace_branches(raw,selected,new_branches,length(new_branches));limits.check(following);
        for(const auto &b:new_branches) framed_parts(b,n,limits.max_word_size);
        J event={{"rule_id",rule.at("id")},{"metric_delta",changes(old,new_branches)}};
        if(trace) {
            event["incoming_frame_transport"]=before.json();event["catalog_transition"]=application.at("transition");
            event["payload_transport"]={{"start",out},{"end",end},{"word",c.at("output_payload")},{"transformed",payload}};
            event["branch_frames"]=witnesses;event["successor_state"]=following;
            if(bm) event["semilinear_derivation"]=rule.at("semilinear");
        }
        return {following,event};
    } catch(const CatalogError &e) {
        if(std::string(e.what()).find("budget")!=std::string::npos) throw LimitExceeded(e.what());
        throw;
    }
}
} // namespace zipper_verify
