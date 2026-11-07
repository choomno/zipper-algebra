#include "core.hpp"
#include "catalog.hpp"
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <openssl/evp.h>

namespace zipper_verify {
namespace {
J sources = J::object();
bool decimal(const std::string &s) {
    if (s.empty()) return false;
    std::size_t i = s[0]=='-' ? 1 : 0;
    if (i==s.size() || (s[i]=='0' && i+1<s.size()) || s=="-0") return false;
    for (;i<s.size();++i) if(s[i]<'0'||s[i]>'9') return false;
    return true;
}
void fields(const J &v, std::initializer_list<const char *> names) {
    if(!v.is_object() || v.size()!=names.size()) throw ReplayError("unexpected or missing record fields");
    for(const char *name:names) if(!v.contains(name)) throw ReplayError("unexpected or missing record fields");
}
bool coefficient_shape(const J &v) {
    if(!v.is_array()) return false;
    for(const auto &pair:v)
        if(!pair.is_array()||pair.size()!=2||!is_int(pair[0])||
           !(pair[1].is_string()||is_int(pair[1]))) return false;
    return true;
}
void emit_json(std::string &out,const J &v,bool coef=false) {
    if(coef && coefficient_shape(v)) {
        out.push_back('['); bool first=true;
        for(const auto &pair:v) {
            if(!first) out.push_back(',');
            first=false;
            out.push_back('['); out+=std::to_string(as_i(pair[0])); out.push_back(',');
            out+=integer_value(pair[1]).get_str(); out.push_back(']');
        }
        out.push_back(']'); return;
    }
    if(v.is_object()) {
        out.push_back('{'); bool first=true;
        for(auto it=v.begin();it!=v.end();++it) {
            if(!first) out.push_back(',');
            first=false;
            out+=J(it.key()).dump(-1,' ',true); out.push_back(':');
            emit_json(out,it.value(),it.key()=="coefficient");
        }
        out.push_back('}');
    } else if(v.is_array()) {
        out.push_back('['); bool first=true;
        for(const auto &item:v) { if(!first) out.push_back(',');
            first=false; emit_json(out,item); }
        out.push_back(']');
    } else {
        if(v.is_binary()||v.is_discarded()) throw ReplayError("unsupported JSON value");
        out+=v.dump(-1,' ',true);
    }
}
}
bool is_int(const J &v) { return v.is_number_integer() || v.is_number_unsigned(); }
I as_i(const J &v) {
    if(!is_int(v)) throw ReplayError("integer required");
    if(v.is_number_unsigned()) {
        auto x=v.get<std::uint64_t>();
        if(x>static_cast<std::uint64_t>(std::numeric_limits<I>::max())) throw LimitExceeded("signed integer coordinate overflow");
        return static_cast<I>(x);
    }
    return v.get<I>();
}
I checked_add(I a,I b) { I r; if(__builtin_add_overflow(a,b,&r)) throw LimitExceeded("signed integer coordinate overflow"); return r; }
I checked_sub(I a,I b) { I r; if(__builtin_sub_overflow(a,b,&r)) throw LimitExceeded("signed integer coordinate overflow"); return r; }
I checked_mul(I a,I b) { I r; if(__builtin_mul_overflow(a,b,&r)) throw LimitExceeded("signed integer coordinate overflow"); return r; }
I neg(I a) { return checked_sub(0,a); }
I abs_i(I a) { return a<0?neg(a):a; }
I floor_div(I a,I b) {
    if(b==0) throw CatalogError("division by zero");
    if(a==std::numeric_limits<I>::min()&&b==-1) throw LimitExceeded("signed integer coordinate overflow");
    I q=a/b,r=a%b; return r && ((r<0)!=(b<0))?checked_sub(q,1):q;
}
I length(const J &v) {
    if(v.size()>static_cast<std::size_t>(std::numeric_limits<I>::max())) throw LimitExceeded("signed integer coordinate overflow");
    return static_cast<I>(v.size());
}
Big integer_value(const J &v) {
    std::string text;
    if(v.is_string()) text=v.get<std::string>();
    else if(is_int(v)) text=v.dump();
    else throw ReplayError("exact integer coefficient required");
    if(!decimal(text)) throw ReplayError("canonical decimal integer coefficient required");
    Big result;
    if(mpz_set_str(result.get_mpz_t(),text.c_str(),10)!=0) throw ReplayError("invalid integer coefficient");
    return result;
}
I integer_bits(const Big &value) {
    if(value==0) return 0;
    const auto bits=mpz_sizeinbase(value.get_mpz_t(),2);
    if(bits>static_cast<std::size_t>(std::numeric_limits<I>::max())) throw LimitExceeded("coefficient integer bit limit exceeded");
    return static_cast<I>(bits);
}
Laurent coefficient(const J &value,bool allow_zero) {
    if(!value.is_array()||(!allow_zero&&value.empty())) throw ReplayError("a nonzero Laurent coefficient is required");
    Laurent result; I last=0;bool first=true;
    for(const auto &pair:value) {
        if(!pair.is_array()||pair.size()!=2||!is_int(pair[0])) throw ReplayError("a nonzero Laurent coefficient is required");
        I q=as_i(pair[0]); Big c=integer_value(pair[1]);
        if(c==0) throw ReplayError("a nonzero Laurent coefficient is required");
        if(!first&&q<=last) throw ReplayError("Laurent exponents must be distinct and increasing");
        result.emplace(q,std::move(c)); first=false;last=q;
    }
    return result;
}
J encoded_coefficient(const Laurent &value) {
    J result=J::array();
    for(const auto &[q,c]:value) if(c!=0) result.push_back(J::array({q,c.get_str()}));
    return result;
}
J one() { return J::array({J::array({0,"1"})}); }
J add_coefficients(const J &a,const J &b) {
    Laurent out=coefficient(a);
    for(const auto &[q,c]:coefficient(b)) out[q]+=c;
    return encoded_coefficient(out);
}
J multiply(const J &a,const J &b) {
    Laurent out;const auto aa=coefficient(a),bb=coefficient(b);
    for(const auto &[q,c]:aa) for(const auto &[p,d]:bb) out[checked_add(q,p)]+=c*d;
    return encoded_coefficient(out);
}
void set_coefficient(J &item,const J &value) {
    item["coefficient"]=value;
    if(item.contains("branch")) item["branch"]["coefficient"]=value;
}
J normalize_coefficients(const J &value) {
    if(value.is_array()) {
        J out=J::array();for(const auto &item:value)out.push_back(normalize_coefficients(item));return out;
    }
    if(value.is_object()) {
        J out=J::object();for(auto it=value.begin();it!=value.end();++it)
            out[it.key()]=it.key()=="coefficient" && coefficient_shape(it.value())?
                encoded_coefficient(coefficient(it.value())):normalize_coefficients(it.value());
        return out;
    }
    return value;
}
Word slice(const Word &word,I begin,I end) {
    if(!word.is_array()) throw ReplayError("word must be an array");
    const I size=length(word);
    auto position=[size](I x) { if(x<0) return x < -size ? I(0) : size+x; return std::min(x,size); };
    begin=position(begin);end=position(end); Word out=J::array();
    for(I k=begin;k<end;++k)out.push_back(word[static_cast<std::size_t>(k)]);
    return out;
}
void append(Word &target,const Word &word) {
    if(!target.is_array()||!word.is_array()) throw ReplayError("word must be an array");
    for(const auto &token:word) target.push_back(token);
}
Word concat(std::initializer_list<Word> words) { Word out=J::array();for(const auto &word:words)append(out,word);return out; }
Word inverse(const Word &word) {
    Word out=J::array();for(auto it=word.rbegin();it!=word.rend();++it) {
        if(!is_int(*it)) throw ReplayError("inverse requires a pure braid word");
        out.push_back(neg(as_i(*it)));
    } return out;
}
Word bar(const Word &word) {
    Word out=J::array();for(auto token:word) {
        if(is_int(token))token=neg(as_i(token));
        else if(is_slider(token)) token[5]=!token[5].get<bool>();
        else throw ReplayError("word contains an invalid braid or slider token");
        out.push_back(std::move(token));
    }return out;
}
Word reverse_word(const Word &word) { Word out=J::array();for(auto it=word.rbegin();it!=word.rend();++it)out.push_back(*it);return out; }
Word phi(const Word &word,I rank) {
    Word out=J::array();for(auto it=word.rbegin();it!=word.rend();++it) {
        J token=*it;
        if(is_int(token)) {I g=as_i(token);token=checked_mul(g>0?1:-1,checked_sub(rank,abs_i(g)));}
        else if(is_slider(token)) token[2]=checked_add(checked_sub(checked_sub(rank,as_i(token[1])),as_i(token[2])),2);
        else throw ReplayError("word contains an invalid braid or slider token");
        out.push_back(std::move(token));
    }return out;
}
Word shifted(const Word &word,I shift) {
    Word out=J::array();for(auto token:word) {
        if(is_int(token)) {I g=as_i(token);token=checked_mul(g>0?1:-1,checked_add(abs_i(g),shift));}
        else if(is_slider(token)) token[2]=checked_add(as_i(token[2]),shift);
        else throw ReplayError("word contains an invalid braid or slider token");
        out.push_back(std::move(token));
    }return out;
}
J slider(I size,I start,bool mirrored,bool reversed,bool barred) {return J::array({"slider",size,start,mirrored,reversed,barred});}
bool is_slider(const J &t) {return t.is_array()&&t.size()==6&&t[0]=="slider"&&is_int(t[1])&&is_int(t[2])&&t[3].is_boolean()&&t[4].is_boolean()&&t[5].is_boolean();}
void validate_word(const Word &word,I rank) {
    if(!word.is_array())throw ReplayError("word contains an invalid braid or slider token");
    for(const auto &token:word) {
        bool valid=false;
        if(is_int(token)) {I g=as_i(token);valid=g!=0&&abs_i(g)<rank;}
        else if(is_slider(token)) {I size=as_i(token[1]),start=as_i(token[2]);valid=size>=1&&size<=rank&&start>=1&&start<=checked_add(checked_sub(rank,size),1);}
        if(!valid)throw ReplayError("word contains an invalid braid or slider token");
    }
}
void validate_branch(const J &b,I rank) {
    fields(b,{"word","coefficient","control"});
    if(!b["control"].is_object())throw ReplayError("branch control must be a JSON object");
    coefficient(b["coefficient"],false);validate_word(b["word"],rank);
}
J branch(const Word &word,const J &coef,const J &control) {
    if(!control.is_object())throw ReplayError("branch control must be a JSON object");
    return J{{"word",word},{"coefficient",encoded_coefficient(coefficient(coef,false))},{"control",control}};
}
J run_state(I rank,const J &branches) {
    if(rank<1)throw ReplayError("rank must be a positive integer");
    if(!branches.is_array())throw ReplayError("state entries must be Branch records");
    for(const auto &b:branches)validate_branch(b,rank);
    return J{{"n",rank},{"branches",branches}};
}
Limits::Limits(const J &value) {
    fields(value,{"max_steps","max_terms","max_word_size","max_coefficient_support","max_integer_bits","expansion_budget"});
    max_steps=as_i(value["max_steps"]);max_terms=as_i(value["max_terms"]);max_word_size=as_i(value["max_word_size"]);
    max_coefficient_support=as_i(value["max_coefficient_support"]);max_integer_bits=as_i(value["max_integer_bits"]);expansion_budget=as_i(value["expansion_budget"]);
    for(I x:{max_steps,max_terms,max_word_size,max_coefficient_support,max_integer_bits,expansion_budget})
        if(x<1)throw ReplayError("execution limits must be positive integers");
}
J Limits::json() const {return J{{"max_steps",max_steps},{"max_terms",max_terms},{"max_word_size",max_word_size},{"max_coefficient_support",max_coefficient_support},{"max_integer_bits",max_integer_bits},{"expansion_budget",expansion_budget}};}
void Limits::check_branch(const J &b) const {
    if(length(b.at("word"))>max_word_size)throw LimitExceeded("word size limit exceeded");
    if(length(b.at("coefficient"))>max_coefficient_support)throw LimitExceeded("coefficient support limit exceeded");
    for(const auto &[q,c]:coefficient(b.at("coefficient"))) { (void)q;if(integer_bits(c)>max_integer_bits)throw LimitExceeded("coefficient integer bit limit exceeded"); }
}
void Limits::check(const J &state) const {
    if(length(state.at("branches"))>max_terms)throw LimitExceeded("term limit exceeded");
    for(const auto &b:state.at("branches"))check_branch(b);
}
J metrics(const J &state) {
    std::map<std::string,Laurent> totals;I tokens=0,max_tokens=0,support=0,max_bits=0,min_q=0,max_q=0;bool seen=false;
    for(const auto &b:state.at("branches")) {
        const auto c=coefficient(b.at("coefficient")); auto &sum=totals[canonical_json(b.at("word"))];
        for(const auto &[q,scalar]:c) {sum[q]+=scalar; max_bits=std::max(max_bits,integer_bits(scalar)); if(!seen){min_q=max_q=q;seen=true;}else{min_q=std::min(min_q,q);max_q=std::max(max_q,q);}}
        tokens=checked_add(tokens,length(b.at("word")));max_tokens=std::max(max_tokens,length(b.at("word")));support=checked_add(support,static_cast<I>(c.size()));
    }
    I count=0;for(const auto &[key,c]:totals) { (void)key;if(std::any_of(c.begin(),c.end(),[](const auto &p){return p.second!=0;}))++count; }
    return J{{"uncollected_terms",length(state.at("branches"))},{"collected_terms",count},{"word_tokens",tokens},{"max_word_tokens",max_tokens},{"coefficient_support",support},{"q_exponent_min",min_q},{"q_exponent_max",max_q},{"max_integer_bits",max_bits}};
}
std::string canonical_json(const J &value) {std::string out;emit_json(out,value);return out;}
std::string sha256(const std::string &text) {
    unsigned char value[EVP_MAX_MD_SIZE];unsigned int size=0;
    if(EVP_Digest(text.data(),text.size(),value,&size,EVP_sha256(),nullptr)!=1)throw Error("SHA256 digest failed");
    std::ostringstream result;result<<std::hex<<std::setfill('0');for(unsigned int i=0;i<size;++i)result<<std::setw(2)<<static_cast<unsigned>(value[i]);return result.str();
}
std::string digest(const J &value) {return sha256(canonical_json(value));}
bool strict_equal(const J &a,const J &b) {
    if(a.is_object()&&b.is_object()) {if(a.size()!=b.size())return false;for(auto it=a.begin();it!=a.end();++it)if(!b.contains(it.key())||!strict_equal(it.value(),b.at(it.key())))return false;return true;}
    if(a.is_array()&&b.is_array()) {if(a.size()!=b.size())return false;for(std::size_t i=0;i<a.size();++i)if(!strict_equal(a[i],b[i]))return false;return true;}
    if(is_int(a)&&is_int(b))return integer_value(a)==integer_value(b);
    return a.type()==b.type()&&a==b;
}
void set_reference_sources(const J &value) {if(!value.is_object())throw ReplayError("source hashes must be an object");for(auto it=value.begin();it!=value.end();++it)if(!it.value().is_string())throw ReplayError("source hash must be a string");sources=value;}
const J &reference_sources() {return sources;}
J source_subset(std::initializer_list<std::string> names) {J result=J::object();for(const auto &name:names){if(!sources.contains(name))throw ReplayError("missing reference source hash: "+name);result[name]=sources.at(name);}return result;}
std::string source_digest(std::initializer_list<std::string> names) {return digest(source_subset(names));}
J request(I branch_index,I offset,const std::string &rule_id,const J &parameters,const J &words,const J &evidence,const J &controls) {
    return J{{"branch",branch_index},{"offset",offset},{"rule_id",rule_id},{"parameters",parameters},{"words",words},{"evidence",evidence},{"successor_controls",controls}};
}
J apply_rule(const J &state,const J &req,Catalog &catalog,const Limits &limits,bool trace) {
    fields(state,{"n","branches"}); const I n=as_i(state.at("n"));run_state(n,state.at("branches"));limits.check(state);
    fields(req,{"branch","offset","rule_id","parameters","words","evidence","successor_controls"});
    I index=as_i(req["branch"]),offset=as_i(req["offset"]);
    if(index<0||offset<0)throw ReplayError("branch and offset must be nonnegative integers");
    if(!req["rule_id"].is_string()||req["rule_id"].get<std::string>().empty())throw ReplayError("rule id must be nonempty");
    if(!req["parameters"].is_object())throw ReplayError("parameters must map names to integers");
    for(const auto &v:req["parameters"])if(!is_int(v))throw ReplayError("parameters must map names to integers");
    if(!req["words"].is_object())throw ReplayError("named words must contain signed braid integers");
    for(const auto &v:req["words"]) {if(!v.is_array())throw ReplayError("named words must contain signed braid integers");for(const auto &x:v)if(!is_int(x))throw ReplayError("named words must contain signed braid integers");}
    if(!req["evidence"].is_object())throw ReplayError("evidence must be an object");
    if(index>=length(state["branches"]))throw ReplayError("selected branch does not exist");
    const J &env=req["parameters"];
    if(!env.contains("n")||as_i(env["n"])!=n)throw ReplayError("rule ambient rank differs from state rank");
    std::string rid=req["rule_id"];
    J rule=rid.size()>=3 && rid.compare(rid.size()-3,3,".BM")==0 ? catalog.semilinear_rule(rid.substr(0,rid.size()-3)+".M",env) : catalog.get(rid);
    if(!rule.contains("execution")||!rule["execution"].contains("algebra_verified")||rule["execution"]["algebra_verified"]!=J(true))throw ReplayError("rule lacks accepted algebraic verification");
    if(!condition(rule.at("condition"),env))throw ReplayError("rule guard is false");
    Polynomial lhs,rhs;
    try {lhs=expand_polynomial(rule.at("lhs"),env,req["words"],false,limits.expansion_budget);rhs=expand_polynomial(rule.at("rhs"),env,req["words"],false,limits.expansion_budget);}
    catch(const CatalogError &e) {std::string s=e.what();if(s.find("budget exceeded")!=std::string::npos)throw LimitExceeded(s);throw ReplayError("invalid catalog instantiation: "+s);}
    if(lhs.size()!=1)throw ReplayError("contextual replay requires a monomial left side");
    if(lhs[0].q!=0||lhs[0].coefficient!=1||lhs[0].word.empty())throw ReplayError("left side must be a nonempty word with coefficient one");
    const J &before=state["branches"][static_cast<std::size_t>(index)];I end=checked_add(offset,length(lhs[0].word));
    if(!strict_equal(slice(before["word"],offset,end),lhs[0].word))throw ReplayError("literal left side does not match selected occurrence");
    std::map<std::string,std::size_t> positions;std::vector<std::pair<Word,Laurent>> groups;
    for(const auto &term:rhs) {std::string key=canonical_json(term.word);auto found=positions.find(key);if(found==positions.end()){positions[key]=groups.size();groups.emplace_back(term.word,Laurent{});found=positions.find(key);}groups[found->second].second[term.q]+=term.coefficient;}
    std::vector<std::pair<Word,J>> terms;
    for(const auto &[word,c]:groups){J encoded=encoded_coefficient(c);if(!encoded.empty())terms.emplace_back(word,std::move(encoded));}
    const J &controls=req["successor_controls"];
    if(!controls.is_null()&&(!controls.is_array()||controls.size()!=terms.size()))throw ReplayError("one control update is required per distinct RHS word");
    J successors=J::array();
    for(std::size_t k=0;k<terms.size();++k) successors.push_back(branch(concat({slice(before["word"],0,offset),terms[k].first,slice(before["word"],end)}),multiply(before["coefficient"],terms[k].second),controls.is_null()?before["control"]:controls[k]));
    J branches=concat({slice(state["branches"],0,index),successors,slice(state["branches"],checked_add(index,1))});J result=run_state(n,branches);limits.check(result);
    J transition=J::object();if(trace)transition=J{{"request",req},{"rule_sha256",digest(rule)},{"predecessor_state_sha256",digest(state)},{"successor_state_sha256",digest(result)},{"predecessor",before},{"successors",successors},{"equality_reference",J{{"catalog_rule",rid},{"kind","accepted_catalog_identity"}}}};
    return J{{"state",result},{"transition",transition}};
}
} // namespace zipper_verify
