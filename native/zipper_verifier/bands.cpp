// Independent local frame and whole-band relations, from the Python specification.
#include "transport.hpp"
#include <algorithm>
#include <sstream>

namespace zipper_verify {
namespace {
struct LocalFrame { std::string kind; I after; J block, emitted; };

void validate_braid(const J &word, I rank) {
    if (!word.is_array()) throw Error("word must contain signed braid generators in B_n");
    for (const auto &token : word)
        if (!is_int(token) || as_i(token) == 0 || abs_i(as_i(token)) >= rank)
            throw Error("word must contain signed braid generators in B_n");
}
I last_start(I n, I size) { return checked_add(checked_sub(n, size), 1); }
void validate_interval(I n, I size, I start, const std::string &message) {
    if (size < 1 || size > n || start < 1 || start > last_start(n, size)) throw Error(message);
}
I frame_block_size(I letter,I size,I start) {
    const I index=abs_i(letter), right=checked_add(start,size);
    return index<checked_sub(start,1)||index>=right||
           (start<=index&&index<checked_sub(right,1))||size==1 ? 1 : size;
}
J consecutive(I low, I count, I sign, bool descending) {
    J result = J::array();
    for (I k = 0; k < count; ++k) {
        const I offset = descending ? checked_sub(checked_sub(count, 1), k) : k;
        result.push_back(checked_mul(sign, checked_add(low, offset)));
    }
    return result;
}
LocalFrame frame_relation(I letter, I size, I start, I power, const std::string &side) {
    const I index=abs_i(letter), sign=letter>0 ? 1 : -1;
    const I left=checked_sub(start,1), right=checked_add(start,size);
    const J singleton=J::array({letter});
    if (index < left || index >= right) return {"far",start,singleton,singleton};
    if (start <= index && index < checked_sub(right,1)) {
        const I reflected = power % 2 != 0
            ? checked_sub(checked_sub(checked_add(checked_mul(2,start),size),2),index) : index;
        return {"internal",start,singleton,J::array({checked_mul(sign,reflected)})};
    }
    if (size == 1) return {"singleton",checked_add(start,index==start ? 1 : -1),singleton,singleton};
    const bool forward = index == (side=="left" ? checked_sub(right,1) : left);
    I after,low; J block;
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
std::string tuple_text(const J &word) {
    std::ostringstream text;
    text << '(';
    for (std::size_t k=0;k<word.size();++k) { if(k) text << ", "; text << as_i(word[k]); }
    if(word.size()==1) text << ',';
    text << ')';
    return text.str();
}
[[noreturn]] void band_failure(const std::string &reason,I offset,I start,const J &remaining) {
    throw TransportError(reason+" at offset "+std::to_string(offset)+", band starts at "+
                         std::to_string(start)+"; remaining="+tuple_text(remaining));
}
bool ints(const J &object,std::initializer_list<const char*> keys) {
    if(!object.is_object()) return false;
    for(const auto *key:keys) if(!object.contains(key)||!is_int(object.at(key))) return false;
    return true;
}
bool braid_array(const J &value,I rank) {
    if(!value.is_array()) return false;
    for(const auto &x:value)
        if(!is_int(x)||as_i(x)==0||abs_i(as_i(x))>=rank) return false;
    return true;
}
J flatten(std::vector<J> pieces,bool reversed) {
    if(reversed) std::reverse(pieces.begin(),pieces.end());
    J result=J::array();
    for(const auto &piece:pieces) for(const auto &token:piece) result.push_back(token);
    return result;
}
} // namespace

J transport_frame(const J &word,I n,I size,I start,I power,const std::string &side) {
    const std::string invalid="valid rank, tape interval, integer power and side required";
    validate_interval(n,size,start,invalid);
    if(side!="left"&&side!="right") throw Error(invalid);
    validate_braid(word,n);
    I position=start,cursor=side=="left" ? length(word) : 0;
    const I total=length(word), last=last_start(n,size);
    std::vector<J> pieces;
    J steps=J::array();
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
        steps.push_back({{"kind",relation.kind},{"offset",offset},{"start_before",position},
                         {"start_after",relation.after},{"block",relation.block},{"emitted",relation.emitted}});
        pieces.push_back(relation.emitted);
        cursor=side=="left" ? offset : checked_add(cursor,count);
        position=relation.after;
    }
    return {{"n",n},{"size",size},{"start",start},{"end",position},{"power",power},{"side",side},
            {"word",word},{"transformed",flatten(pieces,side=="left")},{"steps",steps}};
}

void check_frame(const J &w) {
    if(!w.is_object() || !w.contains("word") || !w["word"].is_array() ||
       !w.contains("transformed") || !w["transformed"].is_array() ||
       !w.contains("steps") || !w["steps"].is_array()) throw Error("immutable frame witness required");
    if(!ints(w,{"n","size","start","power"}) || !w.contains("side") || !w["side"].is_string())
        throw Error("valid rank, tape interval, integer power and side required");
    I n=as_i(w["n"]),size=as_i(w["size"]),start=as_i(w["start"]),power=as_i(w["power"]);
    const std::string side=w["side"].get<std::string>();
    validate_interval(n,size,start,"valid rank, tape interval, integer power and side required");
    if(side!="left"&&side!="right") throw Error("valid rank, tape interval, integer power and side required");
    validate_braid(w["word"],n);
    if(!braid_array(w["transformed"],n)) throw Error("frame witness does not cover the entire identity");
    const I last=last_start(n,size);
    if(!ints(w,{"end"})||as_i(w["end"])<1||as_i(w["end"])>last) throw Error("invalid final tape interval");
    const J &word=w["word"];
    I position=start,cursor=side=="left" ? length(word) : 0;
    std::vector<J> pieces;
    for(const auto &step:w["steps"]) {
        if(!ints(step,{"offset","start_before","start_after"}) ||
           !step.contains("kind") || !step["kind"].is_string() ||
           !step.contains("block") || !braid_array(step["block"],n) || step["block"].empty() ||
           !step.contains("emitted") || !braid_array(step["emitted"],n) ||
           as_i(step["start_before"])!=position) throw Error("invalid frame step");
        const J &block=step["block"];
        const I letter=as_i(side=="left" ? block.back() : block.front());
        if(length(block)!=frame_block_size(letter,size,position)) throw Error("invalid local frame identity");
        const auto expected=frame_relation(letter,size,position,power,side);
        const I offset=side=="left" ? checked_sub(cursor,length(block)) : cursor;
        const I after=as_i(step["start_after"]);
        if(step["kind"]!=expected.kind || after!=expected.after || block!=expected.block ||
           step["emitted"]!=expected.emitted || as_i(step["offset"])!=offset || offset<0 ||
           slice(word,offset,checked_add(offset,length(block)))!=block || after<1 || after>last)
            throw Error("invalid local frame identity");
        pieces.push_back(step["emitted"]);
        cursor=side=="left" ? offset : checked_add(cursor,length(block));
        position=after;
    }
    if(cursor!=(side=="left" ? 0 : length(word)) || position!=as_i(w["end"]) ||
       flatten(pieces,side=="left")!=w["transformed"])
        throw Error("frame witness does not cover the entire identity");
}

J transport_band(const J &word,I n,I size,I start,bool mirrored) {
    validate_interval(n,size,start,"valid rank, band interval and boolean orientation required");
    validate_braid(word,n);
    I position=start,offset=0; J steps=J::array();
    while(offset<length(word)) {
        const I letter=as_i(word.at(static_cast<std::size_t>(offset))), index=abs_i(letter),sign=letter>0 ? 1 : -1;
        I after=position; J block; std::string kind;
        const I left=checked_sub(position,1),right=checked_add(position,size);
        if(index<left||index>=right) {kind="far";block=J::array({letter});}
        else if(size==1) {kind="singleton";block=J::array({letter});after=checked_add(position,index==position ? 1 : -1);}
        else if(size==2&&index==position) {kind="internal_i2";block=J::array({letter});}
        else if(index==checked_sub(right,1)) {
            if(size>checked_sub(length(word),offset)) band_failure("whole-band crossing required",offset,position,slice(word,offset));
            kind="right";block=consecutive(position,size,sign,true);after=checked_add(position,1);}
        else if(index==left) {
            if(size>checked_sub(length(word),offset)) band_failure("whole-band crossing required",offset,position,slice(word,offset));
            kind="left";block=consecutive(left,size,sign,false);after=checked_sub(position,1);}
        else band_failure("internal crossing requires a flype or another explicit witness",offset,position,slice(word,offset));
        if(slice(word,offset,checked_add(offset,length(block)))!=block)
            band_failure("whole-band crossing required",offset,position,slice(word,offset));
        steps.push_back({{"kind",kind},{"offset",offset},{"start_before",position},{"start_after",after},{"block",block}});
        position=after;offset=checked_add(offset,length(block));
    }
    return {{"n",n},{"size",size},{"start",start},{"end",position},{"mirrored",mirrored},{"word",word},{"steps",steps}};
}

void check_transport(const J &w) {
    if(!w.is_object() || !w.contains("word") || !w["word"].is_array() ||
       !w.contains("steps") || !w["steps"].is_array()) throw Error("immutable band transport witness required");
    if(!ints(w,{"n","size","start"}) || !w.contains("mirrored") || !w["mirrored"].is_boolean())
        throw Error("valid rank, band interval and boolean orientation required");
    const I n=as_i(w["n"]),size=as_i(w["size"]),start=as_i(w["start"]);
    validate_interval(n,size,start,"valid rank, band interval and boolean orientation required");
    validate_braid(w["word"],n);
    const I last=last_start(n,size);
    if(!ints(w,{"end"})||as_i(w["end"])<1||as_i(w["end"])>last) throw Error("invalid final band interval");
    I position=start,offset=0;
    for(const auto &step:w["steps"]) {
        if(!step.is_object()) throw Error("transport step record required");
        if(!ints(step,{"offset","start_before","start_after"}) ||
           !step.contains("block") || !braid_array(step["block"],n) || step["block"].empty() ||
           as_i(step["offset"])!=offset || as_i(step["start_before"])!=position)
            throw Error("invalid transport step location or block");
        const J &block=step["block"];
        I letter=as_i(block.front()),sign=letter>0 ? 1 : -1,index=abs_i(letter),after=as_i(step["start_after"]);
        const std::string kind=step.contains("kind")&&step["kind"].is_string() ? step["kind"].get<std::string>() : "";
        bool valid=false;
        if(kind=="far") valid=block.size()==1&&after==position&&(index<checked_sub(position,1)||index>=checked_add(position,size));
        else if(kind=="internal_i2") valid=size==2&&block.size()==1&&index==position&&after==position;
        else if(kind=="singleton") valid=size==1&&block.size()==1&&(index==checked_sub(position,1)||index==position)&&after==checked_add(position,index==position ? 1 : -1);
        else if(kind=="right") valid=length(block)==size&&block==consecutive(position,size,sign,true)&&after==checked_add(position,1);
        else if(kind=="left") valid=length(block)==size&&block==consecutive(checked_sub(position,1),size,sign,false)&&after==checked_sub(position,1);
        if(!valid||slice(w["word"],offset,checked_add(offset,length(block)))!=block||after<1||after>last)
            throw Error("invalid local band transport identity");
        position=after;offset=checked_add(offset,length(block));
    }
    if(offset!=length(w["word"])||position!=as_i(w["end"])) throw Error("transport witness does not cover the entire word");
}
FrameWitness typed_transport_frame(const J &word,I n,I size,I start,I power,const std::string &side) {
    const std::string invalid="valid rank, tape interval, integer power and side required";
    validate_interval(n,size,start,invalid);
    if(side!="left"&&side!="right") throw Error(invalid);
    validate_braid(word,n);
    I position=start,cursor=side=="left" ? length(word) : 0;
    const I total=length(word), last=last_start(n,size);
    std::vector<J> pieces;
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

BandWitness typed_transport_band(const J &word,I n,I size,I start,bool mirrored) {
    validate_interval(n,size,start,"valid rank, band interval and boolean orientation required");
    validate_braid(word,n);
    I position=start,offset=0; std::vector<BandStep> steps;
    while(offset<length(word)) {
        const I letter=as_i(word.at(static_cast<std::size_t>(offset))), index=abs_i(letter),sign=letter>0 ? 1 : -1;
        I after=position; J block; std::string kind;
        const I left=checked_sub(position,1),right=checked_add(position,size);
        if(index<left||index>=right) {kind="far";block=J::array({letter});}
        else if(size==1) {kind="singleton";block=J::array({letter});after=checked_add(position,index==position ? 1 : -1);}
        else if(size==2&&index==position) {kind="internal_i2";block=J::array({letter});}
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

J FrameWitness::json() const {
    J encoded=J::array();
    for(const auto &s:steps)
        encoded.push_back({{"kind",s.kind},{"offset",s.offset},{"start_before",s.start_before},
                           {"start_after",s.start_after},{"block",s.block},{"emitted",s.emitted}});
    return {{"n",n},{"size",size},{"start",start},{"end",end},{"power",power},{"side",side},
            {"word",word},{"transformed",transformed},{"steps",std::move(encoded)}};
}
J BandWitness::json() const {
    J encoded=J::array();
    for(const auto &s:steps)
        encoded.push_back({{"kind",s.kind},{"offset",s.offset},{"start_before",s.start_before},
                           {"start_after",s.start_after},{"block",s.block}});
    return {{"n",n},{"size",size},{"start",start},{"end",end},{"mirrored",mirrored},
            {"word",word},{"steps",std::move(encoded)}};
}
void check_frame(const FrameWitness &w) {
    validate_interval(w.n,w.size,w.start,"valid rank, tape interval, integer power and side required");
    if(w.side!="left"&&w.side!="right") throw Error("valid rank, tape interval, integer power and side required");
    validate_braid(w.word,w.n);
    if(!braid_array(w.transformed,w.n)) throw Error("frame witness does not cover the entire identity");
    const I last=last_start(w.n,w.size);
    if(w.end<1||w.end>last) throw Error("invalid final tape interval");
    I position=w.start,cursor=w.side=="left"?length(w.word):0;
    std::vector<J> pieces;
    for(const auto &s:w.steps) {
        if(!braid_array(s.block,w.n)||s.block.empty()||!braid_array(s.emitted,w.n)||s.start_before!=position)
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
void check_transport(const BandWitness &w) {
    validate_interval(w.n,w.size,w.start,"valid rank, band interval and boolean orientation required");
    validate_braid(w.word,w.n);
    const I last=last_start(w.n,w.size);
    if(w.end<1||w.end>last) throw Error("invalid final band interval");
    I position=w.start,offset=0;
    for(const auto &s:w.steps) {
        if(!braid_array(s.block,w.n)||s.block.empty()||s.offset!=offset||s.start_before!=position)
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
} // namespace zipper_verify
