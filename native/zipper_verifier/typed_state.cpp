#include "typed_state.hpp"
namespace zipper_verify::typed {
Labels decode_labels(const J &j) { Labels r; for(const auto &v:j) r.push_back(zipper_verify::as_i(v)); return r; }
J encode_labels(const Labels &v) {return J(v);}
ComponentRecord decode_record(const J &j) {
 ComponentRecord r;
 r.bottom_ports=decode_labels(j.at("bottom_ports"));r.top_ports=decode_labels(j.at("top_ports"));
 r.original_bottom_ports=decode_labels(j.at("original_bottom_ports"));r.original_top_ports=decode_labels(j.at("original_top_ports"));
 r.twist_powers=decode_labels(j.at("twist_powers"));
 const auto orientation=j.at("orientation").get<std::string>();
 if(orientation!="normal" && orientation!="mirrored") throw ReplayError("invalid component orientation");
 r.mirrored=orientation=="mirrored";
 r.top_route=decode_word(j.at("top_route"));r.bottom_route=decode_word(j.at("bottom_route"));
 r.canonical_top_route=decode_word(j.at("canonical_top_route"));r.canonical_bottom_route=decode_word(j.at("canonical_bottom_route"));
 r.upper_parities=decode_labels(j.at("canonical_parities").at("upper"));r.lower_inverse_parities=decode_labels(j.at("canonical_parities").at("lower_inverse"));
 if(j.contains("inherited_top_route")) r.inherited_top_route=decode_word(j.at("inherited_top_route"));
 if(j.contains("inherited_bottom_route")) r.inherited_bottom_route=decode_word(j.at("inherited_bottom_route"));
 if(j.contains("output_hierarchy")) r.output_hierarchy=decode_labels(j.at("output_hierarchy"));
 if(j.contains("input_hierarchy")) r.input_hierarchy=decode_labels(j.at("input_hierarchy"));
 return r;
}
J encode_record(const ComponentRecord &r) {
 J j={{"bottom_ports",r.bottom_ports},{"top_ports",r.top_ports},{"original_bottom_ports",r.original_bottom_ports},{"original_top_ports",r.original_top_ports},
 {"twist_powers",r.twist_powers},{"orientation",r.mirrored?"mirrored":"normal"},{"top_route",encode_word(r.top_route)},
 {"bottom_route",encode_word(r.bottom_route)},{"canonical_top_route",encode_word(r.canonical_top_route)},{"canonical_bottom_route",encode_word(r.canonical_bottom_route)},
 {"canonical_parities",{{"upper",r.upper_parities},{"lower_inverse",r.lower_inverse_parities}}}};
 if(r.inherited_top_route)j["inherited_top_route"]=encode_word(*r.inherited_top_route);
 if(r.inherited_bottom_route)j["inherited_bottom_route"]=encode_word(*r.inherited_bottom_route);
 if(r.output_hierarchy)j["output_hierarchy"]=*r.output_hierarchy;
 if(r.input_hierarchy)j["input_hierarchy"]=*r.input_hierarchy;
 return j;
}
Context decode_context(const J &j) {return {zipper_verify::as_i(j.at("rank")),zipper_verify::as_i(j.at("size")),decode_word(j.at("component_word")),decode_word(j.at("top_route")),decode_word(j.at("bottom_route")),decode_record(j.at("record"))};}
J encode_context(const Context &c) {return {{"rank",c.rank},{"size",c.size},{"component_word",encode_word(c.component_word)},{"top_route",encode_word(c.top_route)},{"bottom_route",encode_word(c.bottom_route)},{"record",encode_record(c.record)}};}
std::vector<Factor> decode_factors(const J &j) {std::vector<Factor> r; for(const auto &f:j)r.push_back({zipper_verify::as_i(f.at("growth_size")),decode_word(f.at("word"))});return r;}
J encode_factors(const std::vector<Factor> &fs) {J r=J::array();for(const auto &f:fs)r.push_back({{"growth_size",f.growth_size},{"word",encode_word(f.word)}});return r;}
std::string phase_name(const Work &w) {if(std::holds_alternative<Single>(w.phase))return "single";if(std::holds_alternative<ZipUp>(w.phase))return "zip_up";if(std::holds_alternative<ZipDown>(w.phase))return "zip_down";if(std::holds_alternative<Layer>(w.phase))return "layer";return "extract";}
Work decode_work(const J &j) {
 Work w; w.rank=zipper_verify::as_i(j.at("rank"));w.coefficient=zipper_verify::coefficient(j.at("coefficient"));
 for(const auto &c:j.at("contexts"))w.contexts.push_back(decode_context(c));
 w.world_phi=j.at("world_phi").get<bool>();w.boundary_left=decode_word(j.at("boundary_left"));w.boundary_right=decode_word(j.at("boundary_right"));
 w.boundary_left_factors=decode_factors(j.at("boundary_left_factors"));w.boundary_right_factors=decode_factors(j.at("boundary_right_factors"));
 w.depth=zipper_verify::as_i(j.at("depth"));w.bottom_labels=decode_labels(j.at("bottom_labels"));w.top_labels=decode_labels(j.at("top_labels"));
 const auto p=j.at("phase").get<std::string>();
 if(p=="single")w.phase=Single{decode_word(j.at("left")),decode_word(j.at("processed")),zipper_verify::as_i(j.at("position"))};
 else if(p=="zip_up")w.phase=ZipUp{decode_branch(j.at("branch"))};
 else if(p=="zip_down")w.phase=ZipDown{decode_branch(j.at("branch"))};
 else if(p=="layer")w.phase=Layer{decode_word(j.at("word")),zipper_verify::as_i(j.at("payload_length")),zipper_verify::as_i(j.at("completed_sides"))};
 else w.phase=Extract{decode_word(j.at("word"))};
 return w;
}
J encode_work(const Work &w) {
 J cs=J::array();for(const auto &c:w.contexts)cs.push_back(encode_context(c));
 J j={{"rank",w.rank},{"coefficient",encoded_coefficient(w.coefficient)},{"contexts",std::move(cs)},{"world_phi",w.world_phi},
 {"boundary_left",encode_word(w.boundary_left)},{"boundary_right",encode_word(w.boundary_right)},
 {"boundary_left_factors",encode_factors(w.boundary_left_factors)},{"boundary_right_factors",encode_factors(w.boundary_right_factors)},
 {"depth",w.depth},{"bottom_labels",w.bottom_labels},{"top_labels",w.top_labels},{"phase",phase_name(w)}};
 if(const auto *p=std::get_if<Single>(&w.phase)){j["left"]=encode_word(p->left);j["processed"]=encode_word(p->processed);j["position"]=p->position;}
 else if(const auto *p=std::get_if<ZipUp>(&w.phase))j["branch"]=encode_branch(p->branch);
 else if(const auto *p=std::get_if<ZipDown>(&w.phase))j["branch"]=encode_branch(p->branch);
 else if(const auto *p=std::get_if<Layer>(&w.phase)){j["word"]=encode_word(p->word);j["payload_length"]=p->payload_length;j["completed_sides"]=p->completed_sides;}
 else j["word"]=encode_word(std::get<Extract>(w.phase).word);
 return j;
}
Output decode_output(const J &j) {Output o;o.coefficient=zipper_verify::coefficient(j.at("coefficient"));o.word=decode_word(j.at("word"));for(const auto &c:j.at("components"))o.components.push_back(decode_context(c));return o;}
J encode_output(const Output &o) {J cs=J::array();for(const auto &c:o.components)cs.push_back(encode_context(c));return {{"coefficient",encoded_coefficient(o.coefficient)},{"word",encode_word(o.word)},{"components",std::move(cs)}};}
J encode_transition(const Transition &t) {if(!t.error.empty())return {{"error",t.error},{"message",t.message}};J children=J::array(),output=J::array();for(const auto &w:t.children)children.push_back(encode_work(w));for(const auto &o:t.output)output.push_back(encode_output(o));return {{"children",std::move(children)},{"output",std::move(output)},{"unit_evidence",t.evidence}};}
}
