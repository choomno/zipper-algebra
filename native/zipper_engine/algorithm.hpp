#pragma once
#include "core.hpp"
#include "extraction.hpp"

namespace zipper {
struct Context {
    int64_t rank = 0, size = 0;
    Word component_word;
    Braid top_route, bottom_route;
    Json record;
};
struct Work {
    int64_t rank = 0, depth = 0;
    Laurent coefficient;
    std::vector<Context> contexts;
    bool world_phi = false;
    Braid boundary_left, boundary_right;
    std::vector<BoundaryRoute> boundary_left_factors, boundary_right_factors;
    Braid bottom_labels, top_labels;
    std::string phase;
    Braid left, processed;
    int64_t position = 0;
    Branch branch;
    Word word;
    int64_t payload_length = 0, completed_sides = 0;
};
struct Output {
    Laurent coefficient;
    Word word;
    std::vector<Context> components;
};
struct Transition {
    std::vector<Work> children;
    std::vector<Output> output;
    Json event;
    std::string error, message;
};
Work initial_work(int64_t n, const Braid &word);
Transition algorithm_transition(const Work &work, Catalog &catalog, const Limits &limits,
                                bool trace, bool want_evidence = true);
Word work_expression(const Work &work);
size_t work_size(const Work &work);
Json work_json(const Work &work);
Json output_json(const Output &output);
Json transition_json(const Transition &result);
Work work_from_json(const Json &value);
std::string work_key(const Work &work);
std::string output_key(const Output &output);
} // namespace zipper
