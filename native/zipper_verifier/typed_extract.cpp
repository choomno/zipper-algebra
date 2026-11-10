// Independent typed extraction controller with owned typed ribbon witnesses.
#include "typed_state.hpp"
#include "typed_extraction.hpp"
#include <algorithm>
#include <set>
namespace zipper_verify::typed {
namespace {
using Membership=std::vector<unsigned char>;
std::size_t crossing_position(const Token &letter,const Membership &member) {
 I at=checked_sub(abs_i(as_i(letter)),1);
 if(at<0 || static_cast<std::uint64_t>(at)+1>=member.size())throw ControlObligation("endpoint crossing lies outside its active rank");
 return static_cast<std::size_t>(at);
}
Labels bottom_ports(const Word &word,I rank) {
 auto [left,token,right]=slider_parts(word);(void)left;
 I start=token.start,end=checked_add(start,token.size);
 if(rank<1 || start<1 || end>checked_add(rank,1))throw ControlObligation("component slider lies outside its active rank");
 Membership member(static_cast<std::size_t>(rank),0);
 for(I p=1;p<=rank;++p)member[static_cast<std::size_t>(p-1)]=(start<=p && p<end);
 for(const auto &letter:right){auto at=crossing_position(letter,member);std::swap(member[at],member[at+1]);}
 Labels result;for(std::size_t p=0;p<member.size();++p)if(member[p])result.push_back(static_cast<I>(p+1));return result;
}
std::vector<Factor> projected_factors(const std::vector<Factor> &factors,Membership &member) {
 std::vector<Factor> reverse_factors;
 for(auto factor=factors.rbegin();factor!=factors.rend();++factor){
  Word chronological;
  for(auto letter=factor->word.rbegin();letter!=factor->word.rend();++letter){
   auto at=crossing_position(*letter,member);bool a=member[at]!=0,b=member[at+1]!=0;
   if(a && b){I offset=1;for(std::size_t p=0;p<at;++p)if(member[p])offset=checked_add(offset,1);chronological.push_back(as_i(*letter)>0?offset:neg(offset));}
   std::swap(member[at],member[at+1]);
  }
  Word projected=reverse_word(chronological);if(!projected.empty())reverse_factors.push_back({factor->growth_size,std::move(projected)});
 }
 std::reverse(reverse_factors.begin(),reverse_factors.end());return reverse_factors;
}
Big choose(I n,unsigned r){if(n<0)throw ControlObligation("component size cannot be negative");if(n<static_cast<I>(r))return Big(0);Big value(1);for(unsigned k=1;k<=r;++k){value*=Big(std::to_string(n-static_cast<I>(k)+1));value/=k;}return value;}
void check_generated_fans(const std::vector<Factor> &factors,I size){std::set<I> tags;Big total(0);for(const auto &factor:factors){I growth=factor.growth_size;if(growth<3 || growth>size || tags.count(growth))throw ControlObligation("projected endpoint fan exceeds its recorded growth-size bound");Big bound=Big(std::to_string(growth-2))*Big(std::to_string(size-growth+1));Big count(std::to_string(length(factor.word)));if(count>bound)throw ControlObligation("projected endpoint fan exceeds its recorded growth-size bound");tags.insert(growth);total+=count;}if(total>choose(size,3))throw ControlObligation("projected endpoint collar exceeds its size bound");}
std::vector<Factor> reflected_factors(const std::vector<Factor> &factors,I rank){std::vector<Factor> result;for(auto it=factors.rbegin();it!=factors.rend();++it)result.push_back({it->growth_size,phi(it->word,rank)});return result;}
Word flatten_factors(const std::vector<Factor> &factors){Word result;for(const auto &f:factors)append(result,f.word);return result;}
Labels parity_flags(const RibbonWitness &ribbon){Labels result;for(std::size_t p=1;p<ribbon.crossing_counts.size();++p){I value=ribbon.crossing_counts[p];result.push_back(value%2==0?0:1);}return result;}
Labels selected_labels(const Labels &labels,const Labels &ports,bool complement){Labels result;if(!complement){for(I p:ports){if(p<1 || p>static_cast<I>(labels.size()))throw ControlObligation("component endpoint label lies outside its active rank");result.push_back(labels.at(static_cast<std::size_t>(p-1)));}return result;}std::set<I> selected(ports.begin(),ports.end());for(I p=1;p<=static_cast<I>(labels.size());++p)if(!selected.count(p))result.push_back(labels.at(static_cast<std::size_t>(p-1)));return result;}
}
Word lift_contexts(const Word &input,const std::vector<Context> &contexts){Word word=input;for(auto c=contexts.rbegin();c!=contexts.rend();++c)word=concat({c->top_route,c->component_word,shifted(word,c->size),c->bottom_route});return word;}
Transition extraction_transition(const Work &work,const Limits &limits,bool trace){
 const I n=work.rank;Word body=std::get<Extract>(work.phase).word;
 Word local=concat({work.boundary_left,body,work.boundary_right});
 auto top_factors=work.boundary_left_factors,bottom_factors=work.boundary_right_factors;
 if(work.world_phi){local=phi(local,n);body=phi(body,n);auto previous_top=std::move(top_factors);top_factors=reflected_factors(bottom_factors,n);bottom_factors=reflected_factors(previous_top,n);}
 auto [source_left,source_slider,source_right]=slider_parts(local);(void)source_left;(void)source_right;
 if(source_slider.reversed || source_slider.barred)throw ControlObligation("component extraction requires a coupon without reversed or barred flags");
 Labels selected=bottom_ports(local,n);
 LayeredSplit split=split_layered(local,n,selected);
 Word component=split.component_word;I size=static_cast<I>(selected.size());
 Labels bp=split.bottom_ports,tp=split.top_ports;
 ComponentRecord record;record.bottom_ports=bp;record.top_ports=tp;
 record.original_bottom_ports=selected_labels(work.bottom_labels,bp,false);record.original_top_ports=selected_labels(work.top_labels,tp,false);
 J canonical=nullptr;
 if(size==1)component.clear();
 else if(size==2){I power=0;for(const auto &token:component)if(is_int(token))power=checked_add(power,as_i(token));component.clear();I count=abs_i(power);for(I p=0;p<count;++p)component.push_back(power>=0?1:-1);component.push_back(slider(2,1));record.twist_powers={power};}
 else {
  LayeredSplit inner_split=split_layered(body,n,bottom_ports(body,n));
  auto [whole_left,whole_slider,whole_right]=slider_parts(component);
  auto [inner_left,inner_slider,inner_right]=slider_parts(inner_split.component_word);
  I boundary=checked_sub(length(whole_left),length(inner_left));
  if(!(whole_slider==inner_slider) || slice(whole_left,boundary)!=inner_left || slice(whole_right,0,length(inner_right))!=inner_right)throw ControlObligation("inherited endpoint routes do not delimit the component body");
  Word route_left=slice(whole_left,0,boundary),route_right=slice(whole_right,length(inner_right));
  Membership member(static_cast<std::size_t>(n),0);for(I p:selected)member.at(static_cast<std::size_t>(p-1))=1;
  auto inherited_bottom=projected_factors(bottom_factors,member);
  for(auto token=body.rbegin();token!=body.rend();++token)if(is_int(*token)){auto at=crossing_position(*token,member);std::swap(member[at],member[at+1]);}
  auto inherited_top=projected_factors(top_factors,member);
  if(flatten_factors(inherited_top)!=route_left || flatten_factors(inherited_bottom)!=route_right)throw ControlObligation("projected endpoint fan factors disagree with the full component");
  check_generated_fans(inherited_top,size);check_generated_fans(inherited_bottom,size);
  CanonicalSlider witness=normalize_canonical_slider(component,inherited_top,inherited_bottom,limits.expansion_budget,checked_mul(100,limits.expansion_budget));
  Big route_bound=choose(size,3)+choose(size,2)-1;
  I route_length=std::max(length(witness.top_route),length(witness.bottom_route));
  if(Big(std::to_string(route_length))>route_bound)throw ControlObligation("terminal finite route exceeds its size bound");
  component=witness.normalized_word;record.twist_powers=witness.twist_powers;record.inherited_top_route=route_left;record.inherited_bottom_route=route_right;
  record.mirrored=witness.slider.mirrored;record.top_route=witness.top_route;record.bottom_route=witness.bottom_route;
  record.canonical_top_route=slice(record.top_route,length(route_left));record.canonical_bottom_route=slice(record.bottom_route,0,checked_sub(length(record.bottom_route),length(route_right)));
  record.upper_parities=parity_flags(witness.upper);record.lower_inverse_parities=parity_flags(witness.lower_inverse);record.output_hierarchy=witness.upper.labels;record.input_hierarchy=witness.lower_inverse.labels;
  if(trace)canonical=encode_canonical(witness);
 }
 Context context{n,size,std::move(component),split.top_route,split.bottom_route,record};
 auto contexts=work.contexts;contexts.push_back(std::move(context));I remaining=checked_sub(n,size);
 J event={{"rule_id","control.extract"},{"component",encode_record(record)},{"rank_decrease",J::array({n,remaining})}};
 if(trace){event["boundary_split"]=encode_split(split);if(!canonical.is_null())event["nested_slider"]=std::move(canonical);}
 Transition result;result.evidence=std::move(event);
 if(remaining!=0){Work following;following.rank=remaining;following.coefficient=work.coefficient;following.depth=work.depth;following.bottom_labels=selected_labels(work.bottom_labels,bp,true);following.top_labels=selected_labels(work.top_labels,tp,true);following.phase=Single{split.complement_word,{},remaining};following.contexts=std::move(contexts);result.children.push_back(std::move(following));}
 else {Output out;out.coefficient=work.coefficient;out.word=lift_contexts({},contexts);out.components=std::move(contexts);result.output.push_back(std::move(out));}
 return result;
}
}
