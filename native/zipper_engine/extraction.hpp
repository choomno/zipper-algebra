#pragma once
#include "core.hpp"
#include <optional>
#include <vector>

namespace zipper {
struct RibbonError : EngineError {
    using EngineError::EngineError;
};
struct RibbonLimit : RibbonError {
    using RibbonError::RibbonError;
};
struct RibbonGathering {
    int64_t max_image_letters, max_operations;
    std::string action_sha256;
};
struct RibbonBlock {
    int64_t node_size, sign;
    Braid chronological_word;
};
struct NestedRibbon {
    Braid labels, initial_ports, word;
    std::vector<RibbonBlock> blocks;
    Braid crossing_counts, final_ports, reconstructed_word;
    std::optional<RibbonGathering> gathering;
    Braid twist_powers() const;
    Braid boundary_parities() const;
};
Braid reconstruct_nested(const Braid &labels, const Braid &counts, const Braid &initial_ports);
NestedRibbon recognize_nested_strict(const Braid &word, const Braid &labels,
                                     const Braid &initial_ports);
NestedRibbon recognize_nested(const Braid &word, const Braid &labels, const Braid &initial_ports,
                              int64_t max_image_letters = 100000,
                              int64_t max_operations = 10000000);
struct SliderExtraction {
    Word word;
    Token slider;
    NestedRibbon upper, lower;
    Word reconstructed_word;
};
SliderExtraction extract_slider(const Word &, const Braid &upper_labels, const Braid &lower_labels,
                                const Braid &upper_ports, const Braid &lower_ports);
Json extraction_json(const SliderExtraction &);
void check_extraction(const SliderExtraction &);
struct BoundaryRoute {
    int64_t growth_size;
    Braid word;
};
struct LayeredSplit {
    int64_t n;
    Word word;
    Braid bottom_ports, top_ports;
    Word component_word;
    Braid complement_word, bottom_route, top_route;
    Word reconstructed_word;
};
void check_nested(const NestedRibbon &);
Json nested_json(const NestedRibbon &);
Json split_json(const LayeredSplit &);
void check_split(const LayeredSplit &);
LayeredSplit split_layered(const Word &word, int64_t n, const Braid &bottom_ports);
struct CanonicalSlider {
    Word word;
    Token slider;
    NestedRibbon upper, lower_inverse;
    Braid twist_powers, top_route, bottom_route;
    Word normalized_word;
    std::vector<BoundaryRoute> inherited_top, inherited_bottom;
    int64_t max_image_letters, max_operations;
};
Json canonical_json(const CanonicalSlider &);
void check_canonical_slider(const CanonicalSlider &);
CanonicalSlider normalize_canonical_slider(const Word &word,
                                           const std::vector<BoundaryRoute> &top_routes = {},
                                           const std::vector<BoundaryRoute> &bottom_routes = {},
                                           int64_t max_image_letters = 100000,
                                           int64_t max_operations = 10000000);
} // namespace zipper
