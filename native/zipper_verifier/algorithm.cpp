#include "algorithm.hpp"
#include "catalog.hpp"
#include "controllers.hpp"
#include "transport.hpp"
#include <algorithm>
#include <map>
#include <vector>

namespace zipper_verify {
namespace {
J single_work(const J &common,const Word &word) {
    J out=common;
    out["phase"]="single"; out["left"]=word; out["processed"]=J::array();
    out["position"]=common.at("rank");
    return out;
}
J pass_work(J out,const std::string &phase,J branch_value) {
    out["phase"]=phase;
    out["coefficient"]=branch_value.at("coefficient");
    out["branch"]=std::move(branch_value);
    return out;
}
J first_branch_scaled(J &&state,const J &coef) {
    J out=std::move(state.at("branches").at(0));
    out["coefficient"]=coef;
    return out;
}
J phase_result(J children,J evidence) {
    J result=J::object();
    result["children"]=std::move(children);
    result["output"]=J::array();
    result["unit_evidence"]=std::move(evidence);
    return result;
}
J only_child(J child) {
    J children=J::array();
    children.push_back(std::move(child));
    return children;
}
J new_zip(const J &common,const Word &left,const Word &right,I size,I start,
          const Limits &limits,bool mirrored=false,I output_power=0,I input_power=0,
          const Word &output_payload=J::array()) {
    J state;
    I n=as_i(common.at("rank"));
    if(size==2) {
        if(output_power || input_power)
            throw ControlObligation("two-strand framing must be the identity");
        state=upward_initialize(n,concat({left,output_payload}),2,start,0,J::array(),right);
    } else {
        state=framed_initialize(n,left,size,start,0,mirrored,output_power,input_power,
                               J::array(),right,output_payload,"bar_mirror",limits);
    }
    return pass_work(common,"zip_up",first_branch_scaled(std::move(state),common.at("coefficient")));
}
J single_step(const J &work,Catalog &catalog,const Limits &limits) {
    I n=as_i(work.at("rank"));
    const Word &left=work.at("left"), &processed=work.at("processed");
    I position=as_i(work.at("position"));
    if(left.empty()) {
        J out=common_work(work);
        out["phase"]="extract";
        out["word"]=concat({J::array({slider(1,position)}),processed});
        return phase_result(only_child(std::move(out)),{{"rule_id","control.single_scan_complete"}});
    }
    I letter=as_i(left.back()),index=abs_i(letter);
    bool adjacent=position==index || position==checked_add(index,1);
    bool active_over=adjacent && ((position==index)==(letter>0));
    bool offensive=adjacent && !active_over;
    std::string sign=letter>0?"positive":"negative";
    std::string side=position==checked_add(index,1)?"left":"right";
    std::string family=offensive?"foe":adjacent?"close":"far";
    std::string id="single."+family+"."+(adjacent?side+".":"")+sign+".I";
    J env={{"n",n},{"i",1},{"j",position},{"k",0},{"l",index}};
    J instance=catalog.instantiate(id,env,limits.expansion_budget);
    J event={{"rule_id",id},{"parameters",env},{"active_label",n},
             {"active_position",position},{"offensive",offensive},
             {"active_over",adjacent?J(active_over):J(nullptr)}};
    if(!offensive) {
        J out=work;
        out["left"]=slice(left,0,length(left)-1);
        out["processed"]=concat({J::array({letter}),processed});
        if(adjacent) out["position"]=position==index?checked_add(index,1):index;
        return phase_result(only_child(std::move(out)),std::move(event));
    }
    struct Group { Word word; Laurent scalar; };
    std::vector<Group> groups;
    std::map<std::string,std::size_t> index_by_word;
    for(const J &term:instance.at("rhs")) {
        const Word &word=term.at("word");
        std::string key=canonical_json(word);
        auto found=index_by_word.find(key);
        std::size_t at;
        if(found==index_by_word.end()) {
            at=groups.size();index_by_word.emplace(key,at);groups.push_back({word,{}});
        } else at=found->second;
        groups.at(at).scalar[as_i(term.at("q_power"))]+=integer_value(term.at("integer"));
    }
    J children=J::array();
    for(const Group &group:groups) {
        J common=common_work(work);
        common["coefficient"]=multiply(common.at("coefficient"),encoded_coefficient(group.scalar));
        I at=0;
        while(at<length(group.word) && is_int(group.word.at(at))) ++at;
        if(at==length(group.word)) throw ControlObligation("single rewrite is missing its slider");
        const J &coupon=group.word.at(at);
        Word before=concat({slice(left,0,length(left)-1),slice(group.word,0,at)});
        Word after=concat({slice(group.word,at+1),processed});
        if(as_i(coupon.at(1))==2) {
            children.push_back(new_zip(common,before,after,2,as_i(coupon.at(2)),limits));
        } else {
            common["phase"]="single";common["left"]=std::move(before);
            common["processed"]=std::move(after);common["position"]=position;
            children.push_back(std::move(common));
        }
    }
    return phase_result(std::move(children),std::move(event));
}
J zip_step(const J &work,Catalog &catalog,const Limits &limits,bool trace) {
    I n=as_i(work.at("rank"));
    const J &br=work.at("branch"), &c=br.at("control");
    std::string phase=c.at("phase").get<std::string>();
    std::string direction=work.at("phase").get<std::string>();
    if(phase!="upper_pass_complete" && phase!="framed_upper_pass_complete") {
        J state=run_state(n,J::array({br}));
        auto result=as_i(c.at("i"))==2 ? upward_step(state,0,catalog,limits,trace)
                                               : framed_step(state,0,catalog,limits,trace);
        J children=J::array();
        for(J &b:result.first.at("branches"))
            children.push_back(pass_work(common_work(work),direction,std::move(b)));
        return phase_result(std::move(children),std::move(result.second));
    }
    auto turned=turnaround(br,n,catalog,limits,trace,trace);
    auto [cc,exterior,protected_word,tail]=turnaround_parts(turned.first,n,limits.max_word_size);
    if(!protected_word.empty())
        throw ControlObligation("a zipper traversal must begin with an empty processed word");
    J event={{"rule_id","control.turnaround"},{"direction",direction},
             {"rotation_identity","phi fixes slider orientation and sends each half twist to the canonical reflected half twist by Garside symmetry"},
             {"old_child_mirrored",cc.at("mirrored").get<bool>() != (as_i(cc.at("output_power"))%2!=0)}};
    if(trace) event["turnaround"]=std::move(turned.second);
    Word payload=c.value("output_payload",J::array());
    Word friends=as_i(c.at("i"))==2 ? upward_friends(as_i(c.at("i")),as_i(c.at("j")),as_i(c.at("k")))
                                           : friend_word(n,c);
    I prefix_length=checked_add(length(payload),length(friends));
    if(slice(exterior,0,prefix_length)!=concat({payload,friends}))
        throw ControlObligation("turnaround exterior does not preserve its endpoint friend fan");
    I amount=c.value("mirrored",false)?as_i(c.at("k")):neg(as_i(c.at("k")));
    Word body=concat({shifted(payload,amount),slice(exterior,prefix_length)});
    J common=common_work(work);
    if(direction=="zip_up") {
        append(common["boundary_left"],friends);
        if(!friends.empty()) common["boundary_left_factors"].push_back(
            {{"growth_size",checked_add(as_i(cc.at("i")),1)},{"word",friends}});
        I start=checked_add(checked_sub(checked_sub(n,as_i(cc.at("i"))),as_i(cc.at("j"))),2);
        J state=framed_initialize(n,phi(tail,n),as_i(cc.at("i")),start,0,
                                 cc.at("mirrored").get<bool>(),as_i(cc.at("input_power")),
                                 as_i(cc.at("output_power")),J::array(),phi(body,n),J::array(),
                                 c.value("mirror_policy",std::string("bar_mirror")),limits);
        return phase_result(only_child(pass_work(std::move(common),"zip_down",first_branch_scaled(std::move(state),work.at("coefficient")))),std::move(event));
    }
    Word reflected_friends=phi(friends,n);
    common["boundary_right"]=concat({reflected_friends,common.at("boundary_right")});
    if(!friends.empty()) common["boundary_right_factors"]=concat({
        J::array({J{{"growth_size",checked_add(as_i(cc.at("i")),1)},{"word",reflected_friends}}}),
        common.at("boundary_right_factors")});
    Word local=phi(concat({body,slice(turned.first.at("word"),length(exterior))}),n);
    common["phase"]="layer";common["word"]=std::move(local);
    common["payload_length"]=0;common["completed_sides"]=0;
    return phase_result(only_child(std::move(common)),std::move(event));
}
bool outside_over(const Word &word,I n,I size,I start) {
    std::vector<bool> member;
    for(I p=1;p<=n;++p) member.push_back(start<=p && p<checked_add(start,size));
    for(auto it=word.rbegin();it!=word.rend();++it) {
        I letter=as_i(*it),at=checked_sub(abs_i(letter),1);
        bool a=member.at(at),b=member.at(at+1);
        if(a!=b && !(letter>0?a:b)) return true;
        member.at(at)=b;member.at(at+1)=a;
    }
    return false;
}
J reflected_factors(const J &factors,I n) {
    J result=J::array();
    for(auto it=factors.rbegin();it!=factors.rend();++it) {
        J f=*it;f["word"]=phi(f.at("word"),n);result.push_back(std::move(f));
    }
    return result;
}
J layer_step(const J &work,const Limits &limits) {
    I n=as_i(work.at("rank"));
    const Word &word=work.at("word");
    auto [left,coupon,right]=slider_parts(word);
    I size=as_i(coupon.at(1)),start=as_i(coupon.at(2));
    bool mirrored=coupon.at(3).get<bool>();
    if(coupon.at(4)!=J(false) || coupon.at(5)!=J(false))
        throw ControlObligation("layering requires an explicit normal or mirrored coupon without reversed or barred flags");
    I payload_length=as_i(work.at("payload_length"));
    Word unprocessed=payload_length?slice(left,0,length(left)-payload_length):left;
    Word payload=payload_length?slice(left,length(left)-payload_length):J::array();
    if(!outside_over(unprocessed,n,size,start)) {
        J common=common_work(work);
        if(as_i(work.at("completed_sides"))==0) {
            common["phase"]="layer";common["word"]=phi(word,n);
            common["payload_length"]=0;common["completed_sides"]=1;
            common["world_phi"]=!work.at("world_phi").get<bool>();
            common["boundary_left"]=phi(work.at("boundary_right"),n);
            common["boundary_right"]=phi(work.at("boundary_left"),n);
            common["boundary_left_factors"]=reflected_factors(work.at("boundary_right_factors"),n);
            common["boundary_right_factors"]=reflected_factors(work.at("boundary_left_factors"),n);
            return phase_result(only_child(std::move(common)),{{"rule_id","control.layer_change_side"}});
        }
        common["phase"]="extract";common["word"]=word;
        return phase_result(only_child(std::move(common)),{{"rule_id","control.layer_complete"}});
    }
    I letter=as_i(unprocessed.back()),index=abs_i(letter),sign=letter>0?1:-1;
    if(start<=index && index<checked_sub(checked_add(start,size),1)) {
        J out=work;out["payload_length"]=checked_add(payload_length,1);
        return phase_result(only_child(std::move(out)),{{"rule_id","control.layer_absorb_internal"}});
    }
    if(index<checked_sub(start,1) || index>=checked_add(start,size)) {
        Word remaining=slice(unprocessed,0,length(unprocessed)-1);
        J out=common_work(work);out["phase"]="layer";
        out["word"]=concat({remaining,payload,J::array({coupon,letter}),right});
        out["payload_length"]=payload_length;out["completed_sides"]=work.at("completed_sides");
        return phase_result(only_child(std::move(out)),{{"rule_id","control.layer_far"}});
    }
    I move=index==checked_sub(checked_add(start,size),1)?1:-1;
    Word crossing=J::array();
    for(I k=0;k<size;++k) {
        I x=move==1?checked_add(start,k):checked_sub(checked_sub(checked_add(start,size),2),k);
        crossing.push_back(checked_mul(sign,x));
    }
    if(length(unprocessed)<size || slice(unprocessed,length(unprocessed)-size)!=crossing)
        throw ControlObligation("layering requires the complete active-band crossing");
    Word remaining=slice(unprocessed,0,length(unprocessed)-size);
    Word shifted_payload=shifted(payload,move);
    bool offensive=(move==1 && sign<0) || (move==-1 && sign>0);
    J current=common_work(work);current["phase"]="layer";
    current["payload_length"]=length(shifted_payload);current["completed_sides"]=work.at("completed_sides");
    Word good=offensive?bar(crossing):crossing;
    J moved=slider(size,checked_add(start,move),mirrored);
    current["word"]=concat({remaining,shifted_payload,J::array({moved}),good,right});
    J event={{"rule_id","control.layer_whole_band"},{"crossing",crossing},
             {"input_anchor",start},{"output_anchor",checked_add(start,move)},
             {"offensive",offensive},{"orientation",mirrored?"mirrored":"normal"}};
    if(!offensive) return phase_result(only_child(std::move(current)),std::move(event));
    J scalar,growth_scalar;I new_start,a,b;bool new_mirror;
    if(move==1) {
        scalar=J::array({J::array({neg(size),"1"})});growth_scalar=scalar;
        new_start=start;new_mirror=false;a=mirrored?1:0;b=mirrored?-1:0;
    } else {
        scalar=J::array({J::array({size,"1"})});growth_scalar=J::array({J::array({0,"-1"})});
        new_start=checked_sub(start,1);new_mirror=true;a=mirrored?0:-1;b=mirrored?0:1;
    }
    current["coefficient"]=multiply(work.at("coefficient"),scalar);
    J common=common_work(work);common["coefficient"]=multiply(work.at("coefficient"),growth_scalar);
    J growth=new_zip(common,remaining,right,checked_add(size,1),new_start,limits,new_mirror,a,b,shifted_payload);
    event["rule_id"]=move==1?"layer.route.right":"layer.route.left";
    event["growth_frames"]={{"output_power",a},{"input_power",b}};
    event["growth_child_type"]={{"old_child_mirrored",mirrored},
        {"output_payload_mirrored",new_mirror!=(a%2!=0)},
        {"input_payload_mirrored",new_mirror!=(b%2!=0)}};
    event["equality_reference"]={{"catalog_rule",move==1?"layer.positive.I":"layer.opposite_orientation.M"},
        {"full_band_garside_conjugation",move==1?mirrored:!mirrored},
        {"scope","braid-group-algebra slider recursion and whole-band naturality"}};
    J children=J::array();
    children.push_back(std::move(growth));
    children.push_back(std::move(current));
    return phase_result(std::move(children),std::move(event));
}
} // namespace

J common_work(const J &work) {
    J result=J::object();
    for(const char *key:{"rank","coefficient","contexts","world_phi","boundary_left","boundary_right",
                         "boundary_left_factors","boundary_right_factors","depth","bottom_labels","top_labels"})
        result[key]=work.at(key);
    return result;
}
std::tuple<Word,J,Word> slider_parts(const Word &word) {
    I at=-1;
    for(I i=0;i<length(word);++i) if(!is_int(word.at(i))) {
        if(at!=-1) throw ControlObligation("the active local diagram must have exactly one slider");
        at=i;
    }
    if(at==-1) throw ControlObligation("the active local diagram must have exactly one slider");
    return {slice(word,0,at),word.at(at),slice(word,at+1)};
}
Word work_expression(const J &work) {
    I n=as_i(work.at("rank"));
    std::string phase=work.at("phase").get<std::string>();
    Word local;
    if(phase=="single") local=concat({work.at("left"),J::array({slider(1,as_i(work.at("position")))}),work.at("processed")});
    else if(phase=="zip_up" || phase=="zip_down") {
        local=work.at("branch").at("word");
        if(phase=="zip_down") local=phi(local,n);
    } else local=work.at("word");
    local=concat({work.at("boundary_left"),local,work.at("boundary_right")});
    if(work.at("world_phi").get<bool>()) local=phi(local,n);
    return lift_contexts(local,work.at("contexts"));
}
J initial_work(I rank,const Word &word) {
    if(rank<1 || !word.is_array()) throw ReplayError("positive integer rank and literal word required");
    for(const J &x:word)
        if(!is_int(x) || as_i(x)==0 || abs_i(as_i(x))>=rank)
            throw ReplayError("literal input must contain signed generators in the stated rank");
    J labels=J::array();for(I i=1;i<=rank;++i) labels.push_back(i);
    J common={{"rank",rank},{"coefficient",one()},{"contexts",J::array()},{"world_phi",false},
              {"boundary_left",J::array()},{"boundary_right",J::array()},
              {"boundary_left_factors",J::array()},{"boundary_right_factors",J::array()},
              {"depth",0},{"bottom_labels",labels},{"top_labels",labels}};
    return single_work(common,word);
}
J algorithm_transition(const J &work,Catalog &catalog,const Limits &limits,bool trace) {
    // Replay cache callers already normalize the scalar. Borrow that full state
    // when both scalar copies are unit, retaining normalization for diagnostics.
    static const J unit_coefficient=one();
    const bool already_unit=work.at("coefficient")==unit_coefficient &&
        (!work.contains("branch") || work.at("branch").at("coefficient")==unit_coefficient);
    J normalized;
    if(!already_unit) { normalized=work;set_coefficient(normalized,unit_coefficient); }
    const J &unit=already_unit?work:normalized;
    Limits local=limits;
    local.max_terms=std::max<I>(4,limits.max_terms);
    local.max_coefficient_support=std::max<I>(2,limits.max_coefficient_support);
    J result;
    try {
        std::string phase=unit.at("phase").get<std::string>();
        if(phase=="single") result=single_step(unit,catalog,local);
        else if(phase=="zip_up" || phase=="zip_down") result=zip_step(unit,catalog,local,trace);
        else if(phase=="layer") result=layer_step(unit,local);
        else result=extraction_transition(unit,local,trace);
    } catch(const LimitExceeded &e) {
        return {{"error","incomplete"},{"message",e.what()}};
    } catch(const CatalogError &e) {
        std::string message=e.what();
        return {{"error",message.find("budget")!=std::string::npos?"incomplete":"control_obligation"},{"message",message}};
    } catch(const Error &e) {
        return {{"error","control_obligation"},{"message",e.what()}};
    }
    const I next_depth=checked_add(as_i(work.at("depth")),1);
    for(J &child:result.at("children")) {
        set_coefficient(child,multiply(work.at("coefficient"),child.at("coefficient")));
        child["depth"]=next_depth;
    }
    for(J &out:result.at("output")) set_coefficient(out,multiply(work.at("coefficient"),out.at("coefficient")));
    return result;
}
} // namespace zipper_verify
