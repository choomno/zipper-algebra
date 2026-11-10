#include "typed_transport.hpp"
#include "catalog.hpp"
#include "controllers.hpp"
#include <algorithm>
#include <sstream>
namespace zipper_verify::typed {
namespace {
bool pure(const Word&w) {return std::all_of(w.begin(),w.end(),[](const Token&x){return is_int(x);});}
bool valid_braid(const Word&w,I n) {for(const auto& x:w) if(!is_int(x)||as_i(x)==0||abs_i(as_i(x))>=n)return false;return true;}
void validate_braid(const Word&w,I n) {if(!valid_braid(w,n)) throw Error("word must contain signed braid generators in B_n");}
I last_start(I n,I size){return checked_add(checked_sub(n,size),1);}
void validate_interval(I n,I size,I start,const std::string&m){if(size<1||size>n||start<1||start>last_start(n,size))throw Error(m);}
Word consecutive(I low,I count,I sign,bool descending){Word r;for(I k=0;k<count;++k){I offset=descending?checked_sub(checked_sub(count,1),k):k;r.push_back(checked_mul(sign,checked_add(low,offset)));}return r;}
std::string tuple_text(const Word&w){std::ostringstream s;s<<'(';for(std::size_t k=0;k<w.size();++k){if(k)s<<", ";s<<as_i(w[k]);}if(w.size()==1)s<<',';s<<')';return s.str();}
[[noreturn]] void band_failure(const std::string&r,I off,I start,const Word&word){throw TransportError(r+" at offset "+std::to_string(off)+", band starts at "+std::to_string(start)+"; remaining="+tuple_text(word));}
}
namespace {
struct LocalFrame {std::string kind; I after; Word block,emitted;};
I frame_block_size(I letter,I size,I start) {
    const I index=abs_i(letter), right=checked_add(start,size);
    return index<checked_sub(start,1)||index>=right||
           (start<=index&&index<checked_sub(right,1))||size==1 ? 1 : size;
}

LocalFrame frame_relation(I letter, I size, I start, I power, const std::string &side) {
    const I index=abs_i(letter), sign=letter>0 ? 1 : -1;
    const I left=checked_sub(start,1), right=checked_add(start,size);
    const Word singleton{letter};
    if (index < left || index >= right) return {"far",start,singleton,singleton};
    if (start <= index && index < checked_sub(right,1)) {
        const I reflected = power % 2 != 0
            ? checked_sub(checked_sub(checked_add(checked_mul(2,start),size),2),index) : index;
        return {"internal",start,singleton,Word{checked_mul(sign,reflected)}};
    }
    if (size == 1) return {"singleton",checked_add(start,index==start ? 1 : -1),singleton,singleton};
    const bool forward = index == (side=="left" ? checked_sub(right,1) : left);
    I after,low; Word block;
    if (forward) {
        after=checked_add(start,side=="left" ? 1 : -1);
        low=side=="left" ? start : left;
        block=consecutive(low,size,sign,false);
    } else {
        after=checked_add(start,side=="left" ? -1 : 1);
        low=side=="left" ? left : start;
        block=consecutive(low,size,sign,true);
    }
    return {after>start ? "right" : "left",after,block,block};
}

Word flatten(std::vector<Word> pieces,bool reversed){if(reversed)std::reverse(pieces.begin(),pieces.end());Word r;for(const auto&p:pieces)append(r,p);return r;}
}
FrameWitness transport_frame(const Word &word,I n,I size,I start,I power,const std::string &side) {
    const std::string invalid="valid rank, tape interval, integer power and side required";
    validate_interval(n,size,start,invalid);
    if(side!="left"&&side!="right") throw Error(invalid);
    validate_braid(word,n);
    I position=start,cursor=side=="left" ? length(word) : 0;
    const I total=length(word), last=last_start(n,size);
    std::vector<Word> pieces;
    std::vector<FrameStep> steps;
    while(side=="left" ? cursor>0 : cursor<total) {
        const I letter=as_i(word.at(static_cast<std::size_t>(side=="left" ? cursor-1 : cursor)));
        const I expected_size=frame_block_size(letter,size,position);
        if(expected_size>(side=="left" ? cursor : checked_sub(total,cursor))) {
            const I offset=side=="left" ? checked_sub(cursor,expected_size) : cursor;
            throw FrameError("whole marked tape required at offset "+std::to_string(offset)+
                             "; tape starts at "+std::to_string(position)+", side="+side);
        }
        const LocalFrame relation=frame_relation(letter,size,position,power,side);
        const I count=length(relation.block);
        const I offset=side=="left" ? checked_sub(cursor,count) : cursor;
        if(offset<0 || slice(word,offset,checked_add(offset,count))!=relation.block ||
           relation.after<1 || relation.after>last)
            throw FrameError("whole marked tape required at offset "+std::to_string(offset)+
                             "; tape starts at "+std::to_string(position)+", side="+side);
        steps.push_back({relation.kind,offset,position,relation.after,relation.block,relation.emitted});
        pieces.push_back(relation.emitted);
        cursor=side=="left" ? offset : checked_add(cursor,count);
        position=relation.after;
    }
    return {n,size,start,position,power,side,word,flatten(pieces,side=="left"),std::move(steps)};
}
void check_frame(const FrameWitness &w) {
    validate_interval(w.n,w.size,w.start,"valid rank, tape interval, integer power and side required");
    if(w.side!="left"&&w.side!="right") throw Error("valid rank, tape interval, integer power and side required");
    validate_braid(w.word,w.n);
    if(!valid_braid(w.transformed,w.n)) throw Error("frame witness does not cover the entire identity");
    const I last=last_start(w.n,w.size);
    if(w.end<1||w.end>last) throw Error("invalid final tape interval");
    I position=w.start,cursor=w.side=="left"?length(w.word):0;
    std::vector<Word> pieces;
    for(const auto &s:w.steps) {
        if(!valid_braid(s.block,w.n)||s.block.empty()||!valid_braid(s.emitted,w.n)||s.start_before!=position)
            throw Error("invalid frame step");
        const I letter=as_i(w.side=="left"?s.block.back():s.block.front());
        if(length(s.block)!=frame_block_size(letter,w.size,position)) throw Error("invalid local frame identity");
        const auto expected=frame_relation(letter,w.size,position,w.power,w.side);
        const I offset=w.side=="left"?checked_sub(cursor,length(s.block)):cursor;
        if(s.kind!=expected.kind||s.start_after!=expected.after||s.block!=expected.block||s.emitted!=expected.emitted||
           s.offset!=offset||offset<0||slice(w.word,offset,checked_add(offset,length(s.block)))!=s.block||
           s.start_after<1||s.start_after>last) throw Error("invalid local frame identity");
        pieces.push_back(s.emitted);
        cursor=w.side=="left"?offset:checked_add(cursor,length(s.block));
        position=s.start_after;
    }
    if(cursor!=(w.side=="left"?0:length(w.word))||position!=w.end||flatten(pieces,w.side=="left")!=w.transformed)
        throw Error("frame witness does not cover the entire identity");
}
J FrameWitness::json() const {J encoded=J::array();for(const auto&s:steps)encoded.push_back({{"kind",s.kind},{"offset",s.offset},{"start_before",s.start_before},{"start_after",s.start_after},{"block",encode_word(s.block)},{"emitted",encode_word(s.emitted)}});return {{"n",n},{"size",size},{"start",start},{"end",end},{"power",power},{"side",side},{"word",encode_word(word)},{"transformed",encode_word(transformed)},{"steps",std::move(encoded)}};}

BandWitness transport_band(const Word &word,I n,I size,I start,bool mirrored) {
    validate_interval(n,size,start,"valid rank, band interval and boolean orientation required");
    validate_braid(word,n);
    I position=start,offset=0; std::vector<BandStep> steps;
    while(offset<length(word)) {
        const I letter=as_i(word.at(static_cast<std::size_t>(offset))), index=abs_i(letter),sign=letter>0 ? 1 : -1;
        I after=position; Word block; std::string kind;
        const I left=checked_sub(position,1),right=checked_add(position,size);
        if(index<left||index>=right) {kind="far";block=Word{letter};}
        else if(size==1) {kind="singleton";block=Word{letter};after=checked_add(position,index==position ? 1 : -1);}
        else if(size==2&&index==position) {kind="internal_i2";block=Word{letter};}
        else if(index==checked_sub(right,1)) {
            if(size>checked_sub(length(word),offset)) band_failure("whole-band crossing required",offset,position,slice(word,offset));
            kind="right";block=consecutive(position,size,sign,true);after=checked_add(position,1);}
        else if(index==left) {
            if(size>checked_sub(length(word),offset)) band_failure("whole-band crossing required",offset,position,slice(word,offset));
            kind="left";block=consecutive(left,size,sign,false);after=checked_sub(position,1);}
        else band_failure("internal crossing requires a flype or another explicit witness",offset,position,slice(word,offset));
        if(slice(word,offset,checked_add(offset,length(block)))!=block)
            band_failure("whole-band crossing required",offset,position,slice(word,offset));
        steps.push_back({kind,offset,position,after,std::move(block)});
        position=after;offset=checked_add(offset,length(steps.back().block));
    }
    return {n,size,start,position,mirrored,word,std::move(steps)};
}

void check_transport(const BandWitness &w) {
    validate_interval(w.n,w.size,w.start,"valid rank, band interval and boolean orientation required");
    validate_braid(w.word,w.n);
    const I last=last_start(w.n,w.size);
    if(w.end<1||w.end>last) throw Error("invalid final band interval");
    I position=w.start,offset=0;
    for(const auto &s:w.steps) {
        if(!valid_braid(s.block,w.n)||s.block.empty()||s.offset!=offset||s.start_before!=position)
            throw Error("invalid transport step location or block");
        const I letter=as_i(s.block.front()),sign=letter>0?1:-1,index=abs_i(letter),after=s.start_after;
        bool valid=false;
        if(s.kind=="far") valid=s.block.size()==1&&after==position&&(index<checked_sub(position,1)||index>=checked_add(position,w.size));
        else if(s.kind=="internal_i2") valid=w.size==2&&s.block.size()==1&&index==position&&after==position;
        else if(s.kind=="singleton") valid=w.size==1&&s.block.size()==1&&(index==checked_sub(position,1)||index==position)&&after==checked_add(position,index==position?1:-1);
        else if(s.kind=="right") valid=length(s.block)==w.size&&s.block==consecutive(position,w.size,sign,true)&&after==checked_add(position,1);
        else if(s.kind=="left") valid=length(s.block)==w.size&&s.block==consecutive(checked_sub(position,1),w.size,sign,false)&&after==checked_sub(position,1);
        if(!valid||slice(w.word,offset,checked_add(offset,length(s.block)))!=s.block||after<1||after>last)
            throw Error("invalid local band transport identity");
        position=after;offset=checked_add(offset,length(s.block));
    }
    if(offset!=length(w.word)||position!=w.end) throw Error("transport witness does not cover the entire word");
}
J BandWitness::json() const {J encoded=J::array();for(const auto&s:steps)encoded.push_back({{"kind",s.kind},{"offset",s.offset},{"start_before",s.start_before},{"start_after",s.start_after},{"block",encode_word(s.block)}});return {{"n",n},{"size",size},{"start",start},{"end",end},{"mirrored",mirrored},{"word",encode_word(word)},{"steps",std::move(encoded)}};}


namespace {
J encode_event(const ProcessedEvent&e){return {{"rule_id",e.rule_id},{"parameters",e.parameters},{"rhs_index",e.rhs_index},{"mirror_policy",e.mirror_policy}};}
bool suffix(const std::string &text, const std::string &ending) {
    return text.size() >= ending.size() && text.compare(text.size()-ending.size(), ending.size(), ending)==0;
}
Word packet_word(const J &terms, const J &env, I budget) {
    const auto polynomial = expand_polynomial(terms, env, J::object(), false, budget);
    if (polynomial.size()!=1) throw TurnaroundError("a monomial packet word is required");
    const auto &term=polynomial.front();
    if (term.q!=0 || term.coefficient!=1) throw TurnaroundError("packet word must have coefficient one");
    return decode_word(term.word);
}
struct ReplayEvent {
    std::map<std::string,I> parameters;
    I pre_j, post_j, output_delta, input_delta;
    bool pre_mirrored, post_mirrored;
    Word tail;
    std::string mirror_policy;
};
ReplayEvent typed_event(const ProcessedEvent &record, I n, I i, Catalog &catalog, const Limits &limits) {
    if(record.parameters.size()!=6 || !record.parameters.count("n") || !record.parameters.count("i") || !record.parameters.count("j") || !record.parameters.count("k") || !record.parameters.count("l") || !record.parameters.count("eps") || (record.mirror_policy!="M" && record.mirror_policy!="bar_mirror")) throw TurnaroundError("invalid processed-tail event");
    const J env=record.parameters; // catalog expression boundary
    const auto &name=record.rule_id;
    const bool mirrored=suffix(name,".M") || suffix(name,".BM");
    const bool barred=suffix(name,".BM");
    if (record.parameters.at("n")!=n || record.parameters.at("i")!=i || name.rfind("zip.",0)!=0 ||
        !(suffix(name,".I") || mirrored) ||
        (barred && record.mirror_policy!="bar_mirror") ||
        (suffix(name,".M") && record.mirror_policy!="M"))
        throw TurnaroundError("tail event has incompatible rank or orientation");
    const auto held_rule=held_rule_for_event(catalog,J{{"rule_id",name},{"parameters",env}});
    const J &rule=held_rule.get();
    const I index=record.rhs_index;
    if (!rule.contains("execution") || !rule["execution"].is_object() ||
        !rule["execution"].contains("algebra_verified") ||
        !rule["execution"]["algebra_verified"].is_boolean() ||
        !rule["execution"]["algebra_verified"].get<bool>() ||
        !condition(rule["condition"],env) || index<0 || index>=zipper_verify::length(rule["rhs"]))
        throw TurnaroundError("tail event lacks a valid accepted rule branch");
    const I old_j=mirrored ? checked_add(checked_sub(checked_sub(n,i),record.parameters.at("j")),2) : record.parameters.at("j");
    const J &base=catalog.get(rule["base_id"].get<std::string>()+".I");
    auto [jj,kk,mm,aa,bb]=branch_shape(base["rhs"][index],env,mirrored);
    if (barred) {aa=neg(aa);bb=neg(bb);}
    const auto [output,inp]=anchors(i,jj,kk,mm);
    const J &term=rule["rhs"][index];
    const Word local=packet_word(J::array({J{{"coefficient",J::array({J{{"integer",1},{"q_power",0}}})},{"word",term["word"]}}}),env,limits.expansion_budget);
    const Word prefix=concat({frame_power(checked_sub(i,1),output,aa,limits.max_word_size),
        framed_core(n,i,jj,kk,mm,record.mirror_policy),
        frame_power(checked_sub(i,1),inp,bb,limits.max_word_size)});
    if (slice(local,0,length(prefix))!=prefix) throw TurnaroundError("tail event does not have the catalog frame shape");
    return {record.parameters,old_j,jj,aa,bb,mirrored,mm,slice(local,length(prefix)),record.mirror_policy};
}

}
namespace {
J turnaround_certificate(const Branch&before,const Branch&after,I n,Catalog&catalog,const Limits&limits,bool trace,const J&events,I count){return {{"schema","zipper-turnaround-v1"},{"n",n},{"trace_enabled",trace},{"source_sha256",reference_sources().at("turnaround.py")},{"dependency_sha256",source_subset({"framed.py","frames.py","transport.py","catalog.py","upward.py","semilinear.py"})},{"catalog_sha256",catalog.sha256},{"initial",encode_branch(before)},{"final",encode_branch(after)},{"events",events},{"summary",{{"events",count},{"word_tokens_before",length(before.word)},{"word_tokens_after",length(after.word)}}},{"limits",limits.json()}};}
Control completed(I i,I j,bool mirror,const std::string&policy,I a,I b,I exterior,I payload,I routing,I transported,I processed,I tail){Control c;c.phase=ControlPhase::TurnaroundComplete;c.i=i;c.j=j;c.mirrored=mirror;c.mirror_policy=policy;c.output_power=a;c.input_power=b;c.exterior_length=exterior;c.payload_length=payload;c.routing_length=routing;c.transported_length=transported;c.processed_length=processed;c.tail_length=tail;return c;}
}
TurnParts turnaround_parts(const Branch&value,I n,I budget){
    if(value.coefficient.empty()) throw ReplayError("a nonzero Laurent coefficient is required");
    for(const auto &[q,c]:value.coefficient) { (void)q; if(c==0)throw ReplayError("a nonzero Laurent coefficient is required"); }
    validate_word(value.word,n);
    const Control&c=value.control;
    if(c.phase!=ControlPhase::TurnaroundComplete || (c.mirror_policy!="M"&&c.mirror_policy!="bar_mirror"))throw TurnaroundError("invalid turnaround output boundaries");
    I i=c.i,j=c.j;
    if(i<2||i>n||j<1||j>checked_add(checked_sub(n,i),1))throw TurnaroundError("invalid turnaround output boundaries");
    if(std::min({c.exterior_length,c.payload_length,c.routing_length,c.transported_length,c.processed_length,c.tail_length})<0 || c.exterior_length!=checked_add(checked_add(c.payload_length,c.routing_length),c.transported_length))throw TurnaroundError("invalid turnaround output boundaries");
    I end=checked_sub(checked_sub(length(value.word),c.tail_length),c.processed_length);
    Word exterior=slice(value.word,0,c.exterior_length),processed=slice(value.word,end,checked_add(end,c.processed_length)),tail=slice(value.word,checked_add(end,c.processed_length));
    auto [output,inp]=anchors(i,j,0,c.mirrored);
    Word core=concat({frame_power(checked_sub(i,1),output,c.output_power,budget),Word{slider(i,j,c.mirrored)},frame_power(checked_sub(i,1),inp,c.input_power,budget)});
    if(value.word!=concat({exterior,core,processed,tail})||!pure(concat({exterior,processed,tail})))throw TurnaroundError("turnaround output does not match its full expression");
    return {c,exterior,processed,tail};
}
std::pair<Branch,J> turnaround(const Branch&initial,I n,Catalog&catalog,const Limits&limits,bool trace,bool want_certificate){
    check_state(n,{initial},limits);
    if(initial.control.i==2){
        auto [c,u,p,tail]=upward_parts(initial,n);
        if(!u.empty())throw TurnaroundError("completed upward pass required");
        auto band=transport_band(p,n,2,c.j);check_transport(band);
        Word exterior=concat({upward_friends(2,c.j,c.k),p});
        Control cc=completed(2,band.end,false,"M",0,0,length(exterior),0,c.k,length(p),0,length(tail));
        Branch following{concat({exterior,Word{slider(2,band.end)},tail}),initial.coefficient,cc};
        if(!(trace||want_certificate))return {following,J::object()};
        J events=trace?J::array({J{{"band_transport",band.json()}}}):J::array();
        return {following,turnaround_certificate(initial,following,n,catalog,limits,trace,events,static_cast<I>(band.steps.size()))};
    }
    auto [c,u,p,tail]=framed_parts(initial,n,limits.max_word_size);
    if(!u.empty()||c.phase!=ControlPhase::FramedComplete)throw TurnaroundError("completed framed upward pass required");
    I i=c.i,j=c.j,a=c.output_power,b=c.input_power;bool mirror=c.mirrored;
    const std::string policy=c.mirror_policy;const Word payload=c.output_payload;
    Word friends=framed_core(n,i,j,c.k,mirror,policy);
    Word exterior=concat({payload,slice(friends,0,checked_sub(length(friends),1))});
    I payload_length=length(payload),routing_length=checked_mul(checked_sub(i,1),c.k);
    if(static_cast<I>(c.processed_events.size())>limits.max_steps)throw LimitExceeded("turnaround event limit exceeded");
    J evidence=J::array();
    for(const auto&event:c.processed_events){
        auto typed=typed_event(event,n,i,catalog,limits);const Word&T=typed.tail;
        if(typed.post_j!=j||typed.post_mirrored!=mirror||typed.mirror_policy!=policy||slice(p,0,length(T))!=T)throw TurnaroundError("tail event disagrees with current slider or word");
        p=slice(p,length(T));I aa=typed.output_delta,bb=typed.input_delta;bool pre_mirrored=typed.pre_mirrored;
        J packet_id=nullptr;Word W,R;
        if(aa!=0||bb!=0){
            if((aa!=-1&&aa!=1)||bb!=neg(aa)||pre_mirrored==mirror)throw TurnaroundError("invalid emitted flype frame increments");
            bool barred=pre_mirrored&&policy=="bar_mirror";
            std::string sign=aa==(barred?1:-1)?"positive":"negative",variant=barred?".BM":pre_mirrored?".M":".I";
            packet_id="zip.flype."+sign+variant;
            J env=typed.parameters; // catalog expression boundary
            env["j"]=pre_mirrored?checked_add(checked_sub(checked_sub(n,i),j),2):j;env["k"]=0;
            auto held_packet=held_rule_for_event(catalog,J{{"rule_id",packet_id},{"parameters",env},{"rhs_index",0},{"mirror_policy",policy}});const J&packet=held_packet.get();
            Word incoming=packet_word(packet["lhs"],env,limits.expansion_budget);
            if(slice(incoming,checked_sub(length(incoming),1))!=Word{slider(i,j,pre_mirrored)})throw TurnaroundError("flype packet has an invalid input slider");
            W=slice(incoming,0,checked_sub(length(incoming),1));Word outgoing=packet_word(packet["rhs"],env,limits.expansion_budget);
            auto [output,inp]=anchors(i,j,0,mirror);
            Word packet_core=concat({frame_power(checked_sub(i,1),output,aa,limits.max_word_size),Word{slider(i,j,mirror)},frame_power(checked_sub(i,1),inp,bb,limits.max_word_size)});
            if(slice(outgoing,0,length(packet_core))!=packet_core)throw TurnaroundError("flype packet has an invalid retained frame");
            Word C=slice(outgoing,length(packet_core));
            if(length(C)!=checked_sub(i,1)||slice(T,0,length(C))!=C)throw TurnaroundError("tail is missing its complete outgoing flype tape");
            R=slice(T,length(C));a=checked_sub(a,aa);b=checked_sub(b,bb);mirror=pre_mirrored;
        }else{if(pre_mirrored!=mirror)throw TurnaroundError("orientation change requires an explicit flype packet");R=T;}
        auto band=transport_band(R,n,i,j,mirror);check_transport(band);
        if(band.end!=typed.pre_j)throw TurnaroundError("residual tail does not return to the previous band anchor");
        append(exterior,W);append(exterior,R);j=band.end;
        if(trace)evidence.push_back({{"event",encode_event(event)},{"packet_rule",packet_id},{"band_transport",band.json()},{"output_power",a},{"input_power",b},{"j",j},{"mirrored",mirror}});
    }
    if(!c.initial_slider)throw TurnaroundError("typed tails do not restore the initial slider data");
    const auto&original=*c.initial_slider;
    if(j!=original.j||mirror!=original.mirrored||a!=original.output_power||b!=original.input_power||length(p)!=original.processed_length||policy!=original.mirror_policy)throw TurnaroundError("typed tails do not restore the initial slider data");
    auto [output,inp]=anchors(i,j,0,mirror);
    Word word=concat({exterior,frame_power(checked_sub(i,1),output,a,limits.max_word_size),Word{slider(i,j,mirror)},frame_power(checked_sub(i,1),inp,b,limits.max_word_size),p,tail});
    Control cc=completed(i,j,mirror,policy,a,b,length(exterior),payload_length,routing_length,checked_sub(checked_sub(length(exterior),payload_length),routing_length),length(p),length(tail));
    Branch following{std::move(word),initial.coefficient,std::move(cc)};check_state(n,{following},limits);turnaround_parts(following,n,limits.max_word_size);
    if(!(trace||want_certificate))return {following,J::object()};
    return {following,turnaround_certificate(initial,following,n,catalog,limits,trace,evidence,static_cast<I>(c.processed_events.size()))};
}
} // namespace zipper_verify::typed
