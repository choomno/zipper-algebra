#pragma once
#include "typed_controller.hpp"
#include <variant>
namespace zipper_verify::typed {
using Labels=std::vector<I>;
struct Factor { I growth_size=0; Word word; };
struct ComponentRecord {
    Labels bottom_ports,top_ports,original_bottom_ports,original_top_ports,twist_powers;
    bool mirrored=false;
    Word top_route,bottom_route,canonical_top_route,canonical_bottom_route;
    Labels upper_parities,lower_inverse_parities;
    std::optional<Word> inherited_top_route,inherited_bottom_route;
    std::optional<Labels> output_hierarchy,input_hierarchy;
};
struct Context { I rank=0,size=0; Word component_word,top_route,bottom_route; ComponentRecord record; };
struct Single { Word left,processed; I position=0; };
struct ZipUp { Branch branch; };
struct ZipDown { Branch branch; };
struct Layer { Word word; I payload_length=0,completed_sides=0; };
struct Extract { Word word; };
using Phase=std::variant<Single,ZipUp,ZipDown,Layer,Extract>;
struct Work {
    I rank=0; Laurent coefficient; std::vector<Context> contexts; bool world_phi=false;
    Word boundary_left,boundary_right; std::vector<Factor> boundary_left_factors,boundary_right_factors;
    I depth=0; Labels bottom_labels,top_labels; Phase phase;
};
struct Output { Laurent coefficient; Word word; std::vector<Context> components; };
struct Transition { std::vector<Work> children; std::vector<Output> output; J evidence; std::string error,message; };
Labels decode_labels(const J &); J encode_labels(const Labels &);
ComponentRecord decode_record(const J &); J encode_record(const ComponentRecord &);
Context decode_context(const J &); J encode_context(const Context &);
std::vector<Factor> decode_factors(const J &); J encode_factors(const std::vector<Factor> &);
Work decode_work(const J &); J encode_work(const Work &);
Output decode_output(const J &); J encode_output(const Output &);
J encode_transition(const Transition &);
std::string phase_name(const Work &);
Word lift_contexts(const Word &,const std::vector<Context> &);
Transition extraction_transition(const Work &,const Limits &,bool trace);
std::tuple<Word,Slider,Word> slider_parts(const Word &);
Word work_expression(const Work &);
Work initial_work(I,const Word &);
Transition algorithm_transition(const Work &,Catalog &,const Limits &,bool trace);
}
