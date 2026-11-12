// Independent catalog evaluator, implemented from the declarative AST and Python reference.
#include "catalog.hpp"
#include <algorithm>
#include <set>
#include <sstream>

namespace zipper_verify {
namespace {
bool starts(const std::string &s, const std::string &prefix) {
    return s.compare(0, prefix.size(), prefix) == 0;
}
bool ends(const std::string &s, const std::string &suffix) {
    return s.size() >= suffix.size() && s.compare(s.size()-suffix.size(), suffix.size(), suffix)==0;
}
void keys(const J &node, std::initializer_list<const char *> required,
          std::initializer_list<const char *> optional={}) {
    if (!node.is_object()) throw CatalogError("Unexpected/missing fields in word AST");
    std::set<std::string> allowed;
    for (const auto *key: required) {
        allowed.insert(key);
        if (!node.contains(key)) throw CatalogError("Unexpected/missing fields in word AST");
    }
    for (const auto *key: optional) allowed.insert(key);
    for (auto it=node.begin(); it!=node.end(); ++it)
        if (!allowed.count(it.key())) throw CatalogError("Unexpected/missing fields in word AST");
}
std::string poly_key(I q, const Word &word) {
    return std::to_string(q) + ":" + word.dump();
}
class Accumulator {
    std::map<std::string, size_t> index_;
    Polynomial terms_;
public:
    void add(I q, const Word &word, const Big &coefficient) {
        auto key=poly_key(q,word);
        auto inserted=index_.emplace(key,terms_.size());
        if (inserted.second) terms_.push_back({q,word,coefficient});
        else terms_[inserted.first->second].coefficient += coefficient;
    }
    size_t size() const { return terms_.size(); }
    Polynomial clean() {
        Polynomial result;
        for (auto &term:terms_) if (term.coefficient != 0) result.push_back(std::move(term));
        return result;
    }
};
Polynomial unit() { return {{0,J::array(),Big(1)}}; }
Polynomial product(const Polynomial &a,const Polynomial &b,I budget) {
    if (a.size() && b.size() > static_cast<size_t>(budget)/a.size())
        throw CatalogError("Expansion term budget exceeded");
    size_t left=0,right=0;
    for (const auto &t:a) left=std::max(left,t.word.size());
    for (const auto &t:b) right=std::max(right,t.word.size());
    if (!a.empty() && !b.empty() && (left>static_cast<size_t>(budget) || right>static_cast<size_t>(budget)-left))
        throw CatalogError("Expansion word budget exceeded");
    Accumulator out;
    for (const auto &x:a) for (const auto &y:b)
        out.add(checked_add(x.q,y.q),concat({x.word,y.word}),x.coefficient*y.coefficient);
    return out.clean();
}
Polynomial expand_slider(I size,I start,I budget) {
    if (size==1) return unit();
    Polynomial result{{1,J::array({neg(start)}),Big(1)}, {0,J::array(),Big(1)},
                      {1,J::array(),Big(-1)}, {0,J::array({start}),Big(-1)}};
    for (I s=3; s<=size; ++s) {
        Word positive=J::array(),negative=J::array();
        for (I x=0;x<s-1;++x) {
            I g=checked_add(start,x);positive.push_back(g);negative.push_back(neg(g));
        }
        result=product({{s-1,negative,Big(1)},{0,positive,Big(-1)}},result,budget);
    }
    return result;
}
J negative_expression(const J &value) {
    return is_int(value) ? J(neg(as_i(value))) : J{{"op","sub"},{"args",J::array({0,value})}};
}
std::vector<J> slider_ranks(const J &word,bool inverted=false) {
    std::string kind=word.value("type","");
    if (kind=="slider") {
        if (inverted) throw CatalogError("semilinear zipper transformation cannot invert a slider");
        return {word.at("size")};
    }
    std::vector<J> result;
    if (word.contains("word")) {
        auto child=slider_ranks(word.at("word"),inverted||kind=="inverse");
        result.insert(result.end(),child.begin(),child.end());
    }
    if (word.contains("items")) for (const auto &item:word.at("items")) {
        auto child=slider_ranks(item,inverted);result.insert(result.end(),child.begin(),child.end());
    }
    return result;
}
J semilinear_template(const J &source) {
    auto rid=source.at("id").get<std::string>();
    auto source_hash=digest(source);
    J result=source;
    result["id"]=rid.substr(0,rid.size()-2)+".BM";
    result["source_record_id"]=rid;
    result["source_rule_sha256"]=source_hash;
    for (auto side:{"lhs","rhs"}) for (auto &term:result[side]) {
        term["word"]=semilinear_word(term.at("word"));
        for (auto &coefficient:term["coefficient"])
            coefficient["q_power"]=negative_expression(coefficient.at("q_power"));
    }
    result["semilinear"]={{"proof_version","zipper-semilinear-bar-v1"},
        {"source_record_id",rid},{"source_rule_sha256",source_hash},
        {"coefficient_map","q -> q^-1"},{"crossing_map","sigma_l -> sigma_l^-1"},
        {"preserves_product_order",true},{"half_twist_identity","bar(Delta) = inverse(Delta)"},
        {"author_review","inherits partial source review; semilinear derivation recorded separately"}};
    result["symmetry"]["name"]="BM";result["symmetry"]["invert_q"]=true;
    result["symmetry"]["invert_crossing_sign"]=true;result["symmetry"]["semilinear"]=true;
    result["execution"]["derived_semilinear"]=true;
    for (auto key:{"original_relator_cofactor","cofactor_original_notation","cofactor_transport"}) {
        if (result.contains(key)) {
            result[std::string("semilinear_source_")+key]=result.at(key);result.erase(key);
        }
    }
    return result;
}
J bind_semilinear(J result,const J &env) {
    I i=as_i(env.at("i"));
    // Divide before multiplying to retain every representable triangular number.
    I a=i,b=i-1;if (a%2==0) a/=2;else b/=2;
    result["semilinear"]["parameters"]=env;
    result["semilinear"]["common_slider_unit"]={{"integer",(i-1)%2 ? -1:1},
        {"q_power",neg(checked_mul(a,b))}};
    result["semilinear"]["normalized_relator_unit"]={{"integer",-1},{"q_power",neg(i-1)}};
    return result;
}
std::string token_repr(const J &token) {
    if (is_int(token)) return std::to_string(as_i(token));
    return "('slider', "+std::to_string(as_i(token[1]))+", "+std::to_string(as_i(token[2]))+", "+
        (token[3].get<bool>()?"True":"False")+", "+(token[4].get<bool>()?"True":"False")+", "+
        (token[5].get<bool>()?"True":"False")+")";
}
std::string python_key_repr(const Monomial &term) {
    std::string word="(";
    for (size_t i=0;i<term.word.size();++i) {
        if (i) word+=", ";
        word+=token_repr(term.word[i]);
    }
    if (term.word.size()==1) word+=",";
    word+=")";
    return "("+std::to_string(term.q)+", "+word+")";
}
} // namespace

I integer(const J &expression,const J &env) {
    if (is_int(expression)) return as_i(expression);
    if (expression.is_string()) {
        static const std::set<std::string> parameters{"n","i","j","k","l","eps"};
        auto name=expression.get<std::string>();
        if (parameters.count(name)) {
            if (!env.contains(name)||!is_int(env.at(name))) throw CatalogError("Missing integer parameter: "+name);
            return as_i(env.at(name));
        }
    }
    if (expression.is_object() && expression.size()==2 && expression.contains("op") && expression.contains("args")) {
        auto op=expression.at("op");const auto &args=expression.at("args");
        if (!op.is_string() || (op!="add" && op!="sub" && op!="mul") || !args.is_array() || args.size()<2)
            throw CatalogError("Invalid integer operation");
        std::vector<I> values;for (const auto &arg:args) values.push_back(integer(arg,env));
        if (op=="sub") {
            if (values.size()!=2) throw CatalogError("sub has arity two");
            return checked_sub(values[0],values[1]);
        }
        I result=op=="add"?0:1;
        for (I value:values) result=op=="add"?checked_add(result,value):checked_mul(result,value);
        return result;
    }
    throw CatalogError("Invalid integer expression: "+expression.dump());
}
bool condition(const J &expression,const J &env) {
    if (expression.is_boolean()) return expression.get<bool>();
    if (!expression.is_object()) throw CatalogError("Invalid condition");
    auto op=expression.value("op",J());
    if ((op=="and"||op=="or") && expression.size()==2 && expression.contains("args") && expression.at("args").is_array()) {
        bool result=op=="and";
        // Eager evaluation is intentional: invalid operands cannot hide behind short-circuiting.
        for (const auto &arg:expression.at("args")) {
            bool value=condition(arg,env);result=op=="and"?(result&&value):(result||value);
        }
        return result;
    }
    if (op=="not" && expression.size()==2 && expression.contains("arg")) return !condition(expression.at("arg"),env);
    if ((op=="eq"||op=="ne"||op=="lt"||op=="le"||op=="gt"||op=="ge") && expression.size()==3 && expression.contains("left") && expression.contains("right")) {
        I a=integer(expression.at("left"),env),b=integer(expression.at("right"),env);
        if (op=="eq") return a==b;
        if (op=="ne") return a!=b;
        if (op=="lt") return a<b;
        if (op=="le") return a<=b;
        if (op=="gt") return a>b;
        return a>=b;
    }
    throw CatalogError("Invalid condition operation");
}

Polynomial expand_word(const J &node,const J &env,const J &words,bool expand_sliders,I max_terms) {
    if (max_terms<1) throw CatalogError("Expansion budget must be positive");
    if (!node.is_object()) throw CatalogError("Word must be an AST object");
    std::string kind=node.value("type","");
    if (!env.contains("n")||!is_int(env.at("n"))||as_i(env.at("n"))<1) throw CatalogError("n must be a positive integer");
    I n=as_i(env.at("n"));
    auto iv=[&](const char *key){return integer(node.at(key),env);};
    auto expand=[&](const J &child){return expand_word(child,env,words,expand_sliders,max_terms);};
    auto braid=[&](const Word &word)->Polynomial {
        if (!word.is_array()) throw CatalogError("Braid generator outside B_n");
        if (word.size()>static_cast<size_t>(max_terms)) throw CatalogError("Expansion word budget exceeded");
        for (const auto &token:word) if (!is_int(token)||as_i(token)==0||abs_i(as_i(token))>=n)
            throw CatalogError("Braid generator outside B_n");
        return {{0,word,Big(1)}};
    };
    if (kind=="identity") {keys(node,{"type"});return unit();}
    if (kind=="sigma") {
        keys(node,{"type","index"},{"sign"});I index=iv("index"),sign=integer(node.value("sign",J(1)),env);
        if (index<1||(sign!=-1&&sign!=1)) throw CatalogError("Invalid signed generator");
        return braid(J::array({checked_mul(sign,index)}));
    }
    if (kind=="fat"||kind=="block") {
        if (kind=="block") keys(node,{"type","length","start","count","order"},{"sign"});
        else keys(node,{"type","length","start"},{"sign"});
        I length_=iv("length"),start=iv("start"),count=kind=="block"?iv("count"):1;
        I sign=integer(node.value("sign",J(1)),env);
        if (length_<0||count<0||start<1||(sign!=-1&&sign!=1)) throw CatalogError("Invalid block parameters");
        std::string order=node.value("order","ascending");
        if (order!="ascending"&&order!="descending") throw CatalogError("Invalid block order");
        if (!length_||!count) return unit();
        if (length_>max_terms/count) throw CatalogError("Expansion word budget exceeded");
        Word word=J::array();
        for (I t=0;t<count;++t) {
            I offset=order=="ascending"?t:count-1-t;
            for (I x=0;x<length_;++x) word.push_back(checked_mul(sign,checked_add(checked_add(start,offset),x)));
        }
        return braid(word);
    }
    if (kind=="half_twist") {
        keys(node,{"type","size","start"});I size=iv("size"),start=iv("start");
        if (size<1||start<1||checked_add(start,size-1)>n) throw CatalogError("Invalid half-twist interval");
        I a=size,b=size-1;if (a%2==0) a/=2;else b/=2;
        if (a && b>max_terms/a) throw CatalogError("Expansion word budget exceeded");
        Word word=J::array();for (I width=size-1;width>0;--width)
            for (I x=0;x<width;++x) word.push_back(checked_add(start,x));
        return braid(word);
    }
    if (kind=="slider") {
        keys(node,{"type","size","start"},{"orientation"});I size=iv("size"),start=iv("start");
        std::string orientation=node.value("orientation","normal");
        if (size<1||start<1||checked_add(start,size-1)>n||(orientation!="normal"&&orientation!="mirrored"))
            throw CatalogError("Invalid slider");
        bool mirrored=orientation=="mirrored";
        if (expand_sliders&&size>max_terms) throw CatalogError("Expansion word budget exceeded");
        if (!expand_sliders) return {{0,J::array({slider(size,start,mirrored,false,false)}),Big(1)}};
        auto result=expand_slider(size,start,max_terms);
        if (mirrored) for (auto &term:result) for (auto &token:term.word) {
            I g=as_i(token);I reflection=checked_sub(checked_add(checked_mul(2,start),size-2),abs_i(g));
            token=g>0?reflection:neg(reflection);
        }
        return result;
    }
    if (kind=="concat") {
        keys(node,{"type","items"});if (!node.at("items").is_array()) throw CatalogError("concat.items must be a list");
        auto out=unit();for (const auto &child:node.at("items")) out=product(out,expand(child),max_terms);return out;
    }
    if (kind=="word_ref") {
        keys(node,{"type","name"});if (!node.at("name").is_string()||!words.contains(node.at("name").get<std::string>()))
            throw CatalogError("Unbound word reference");
        return braid(words.at(node.at("name").get<std::string>()));
    }
    if (kind=="inverse"||kind=="bar"||kind=="reverse"||kind=="reflect") {
        if (kind=="reflect") keys(node,{"type","word"},{"n"});else keys(node,{"type","word"});
        auto original=expand(node.at("word"));
        if (kind=="inverse" && (original.size()!=1||original[0].coefficient!=1||original[0].q!=0||
            std::any_of(original[0].word.begin(),original[0].word.end(),[](const J &t){return !is_int(t);})))
            throw CatalogError("inverse is only defined for braid words, not slider polynomials");
        I width=integer(node.value("n",J(n)),env);
        if (kind=="reflect"&&width!=n) throw CatalogError("reflect.n must be the ambient n");
        Accumulator out;
        for (auto &term:original) {
            Word input=term.word;if (kind=="inverse"||kind=="reverse") std::reverse(input.begin(),input.end());
            Word transformed=J::array();
            for (const auto &token:input) {
                if (is_int(token)) {
                    I value=as_i(token);
                    if (kind=="inverse"||kind=="bar") value=neg(value);
                    else if (kind=="reflect") value=value>0?checked_sub(n,abs_i(value)):neg(checked_sub(n,abs_i(value)));
                    transformed.push_back(value);
                } else {
                    J token_=token;
                    if (kind=="bar") token_[5]=!token_[5].get<bool>();
                    else if (kind=="reverse") token_[4]=!token_[4].get<bool>();
                    else {
                        token_[2]=checked_add(checked_sub(checked_sub(n,as_i(token[2])),as_i(token[1])),2);
                        token_[3]=!token[3].get<bool>();
                    }
                    transformed.push_back(token_);
                }
            }
            out.add(term.q,transformed,term.coefficient);
        }
        return out.clean();
    }
    throw CatalogError("Unknown word node: "+kind);
}
Polynomial expand_polynomial(const J &terms,const J &env,const J &words,bool expand_sliders,I max_terms) {
    if (!terms.is_array()) throw CatalogError("Polynomial must be a list");
    Accumulator out;
    for (const auto &term:terms) {
        if (!term.is_object()) throw CatalogError("Invalid polynomial term");
        keys(term,{"coefficient","word"});
        if (!term.at("coefficient").is_array()) throw CatalogError("Coefficient must be Laurent term list");
        auto word=expand_word(term.at("word"),env,words,expand_sliders,max_terms);
        for (const auto &coeff:term.at("coefficient")) {
            keys(coeff,{"integer","q_power"});
            if (!is_int(coeff.at("integer"))) throw CatalogError("Coefficient numerator must be integer");
            Big scalar=integer_value(coeff.at("integer"));I power=integer(coeff.at("q_power"),env);
            for (const auto &piece:word) out.add(checked_add(piece.q,power),piece.word,scalar*piece.coefficient);
        }
        if (out.size()>static_cast<size_t>(max_terms)) throw CatalogError("Expansion term budget exceeded");
    }
    return out.clean();
}
J serialize_polynomial(const Polynomial &poly) {
    std::vector<const Monomial *> order;for (const auto &term:poly) order.push_back(&term);
    std::stable_sort(order.begin(),order.end(),[](const Monomial *a,const Monomial *b){return python_key_repr(*a)<python_key_repr(*b);});
    J result=J::array();for (auto *term:order) result.push_back({{"integer",term->coefficient.get_str()},{"q_power",term->q},{"word",term->word}});
    return result;
}
J semilinear_word(const J &word) {
    J node=word;std::string kind=node.value("type","");
    if (kind=="identity"||kind=="slider") return node;
    if (kind=="sigma"||kind=="fat"||kind=="block") {node["sign"]=negative_expression(node.value("sign",J(1)));return node;}
    if (kind=="half_twist") return {{"type","inverse"},{"word",node}};
    if (kind=="concat") {for (auto &child:node["items"]) child=semilinear_word(child);return node;}
    if (kind=="inverse"||kind=="bar"||kind=="reverse"||kind=="reflect") {
        if (kind=="inverse") slider_ranks(node);
        node["word"]=semilinear_word(node.at("word"));return node;
    }
    if (kind=="word_ref") return {{"type","bar"},{"word",node}};
    throw CatalogError("unknown word in semilinear zipper transformation");
}
void validate_semilinear_bindings(const J &rule,const J &env) {
    const auto rid=rule.value("id",J());
    if (!rid.is_string()||!starts(rid.get<std::string>(),"zip.")||!ends(rid.get<std::string>(),".M")||
        !rule.contains("execution")||rule.at("execution").value("algebra_verified",J())!=J(true))
        throw CatalogError("an accepted frozen mirrored zipper rule is required");
    if (!env.contains("i")||!is_int(env.at("i"))||as_i(env.at("i"))<2||!condition(rule.at("condition"),env))
        throw CatalogError("valid original bindings and a true zipper guard are required");
    for (const auto side:{"lhs","rhs"}) {
        if (rule.at(side).empty()) throw CatalogError("nonempty zipper polynomial required");
        for (const auto &term:rule.at(side)) {
            auto ranks=slider_ranks(term.at("word"));
            if (ranks.size()!=1||integer(ranks[0],env)!=as_i(env.at("i")))
                throw CatalogError("every zipper term must contain one slider of the same rank");
        }
    }
}
J semilinear_zip_rule(const J &rule,const J &env) {
    validate_semilinear_bindings(rule,env);return bind_semilinear(semilinear_template(rule),env);
}
J instantiate_rule(const J &rule,const J &env,const J &words,bool expand_sliders,I max_terms,bool require_verified) {
    bool verified=rule.value("execution",J::object()).value("algebra_verified",J())==J(true);
    if (require_verified&&!verified) throw CatalogError("Rule lacks reviewed algebraic acceptance");
    if (!condition(rule.at("condition"),env)) throw CatalogError("Rule condition is false");
    J result{{"id",rule.at("id")},{"parameters",env},{"algebra_verified",verified},
        {"scope","Formal instantiation only; control preconditions are not checked."}};
    for (const auto side:{"lhs","rhs"}) result[side]=serialize_polynomial(expand_polynomial(rule.at(side),env,words,expand_sliders,max_terms));
    return result;
}
Catalog::Catalog(const J &catalog,std::string source_sha256):catalog_(catalog),sha256(std::move(source_sha256)) {
    if (!catalog_.is_object()||!catalog_.contains("rules")||!catalog_.at("rules").is_array())
        throw CatalogError("catalog: expected rule array");
    for (const auto &rule:catalog_.at("rules")) {
        if (!rule.is_object()||!rule.contains("id")||!rule.at("id").is_string()||
            !rule.contains("lhs")||!rule.contains("rhs")||!rule.contains("condition"))
            throw CatalogError("catalog: invalid rule record");
        auto id=rule.at("id").get<std::string>();
        if (!rules_.emplace(id,rule).second) throw CatalogError("Duplicate rule id");
        order_.push_back(id);
    }
}
J Catalog::get(const std::string &rule_id) const {
    auto found=rules_.find(rule_id);if (found==rules_.end()) throw CatalogError("Unknown rule id: "+rule_id);
    return found->second;
}
std::vector<J> Catalog::rules() const {
    std::vector<J> result;for (const auto &id:order_) result.push_back(rules_.at(id));return result;
}
std::vector<std::string> Catalog::select_rule_ids(const J &env,const std::string &prefix,const std::string &suffix,
        const std::vector<std::string> &exclude_prefixes) const {
    std::vector<std::string> result;
    for (const auto &id:order_) {
        if (!starts(id,prefix)||!ends(id,suffix)) continue;
        if (std::any_of(exclude_prefixes.begin(),exclude_prefixes.end(),[&](const std::string &x){return starts(id,x);})) continue;
        if (condition(rules_.at(id).at("condition"),env)) result.push_back(id);
    }
    return result;
}
std::vector<J> Catalog::select_rules(const J &env,const std::string &prefix,const std::string &suffix,
        const std::vector<std::string> &exclude_prefixes) const {
    std::vector<J> result;for (const auto &id:select_rule_ids(env,prefix,suffix,exclude_prefixes)) result.push_back(get(id));return result;
}
J Catalog::semilinear_rule(const std::string &source_id,const J &env) {
    J rule=get(source_id);validate_semilinear_bindings(rule,env);
    auto found=semilinear_templates_.find(source_id);
    if (found==semilinear_templates_.end()) found=semilinear_templates_.emplace(source_id,semilinear_template(rule)).first;
    return bind_semilinear(found->second,env);
}
J Catalog::instantiate(const std::string &rule_id,const J &env,const J &words,bool expand_sliders,I max_terms,bool require_verified) const {
    return instantiate_rule(get(rule_id),env,words,expand_sliders,max_terms,require_verified);
}
} // namespace zipper_verify
