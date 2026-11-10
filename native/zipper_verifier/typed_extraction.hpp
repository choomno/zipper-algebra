#pragma once
#include "typed_state.hpp"
#include "typed_ribbons.hpp"
namespace zipper_verify::typed {
struct LayeredSplit {
    I n = 0;
    Word word, component_word, complement_word, bottom_route, top_route, reconstructed_word;
    Labels bottom_ports, top_ports;
};
struct CanonicalSlider {
    Word word;
    Slider slider;
    RibbonWitness upper, lower_inverse;
    Labels twist_powers;
    Word top_route, bottom_route, normalized_word;
    std::vector<Factor> inherited_top, inherited_bottom;
    I max_image_letters = 0, max_operations = 0;
};
LayeredSplit split_layered(const Word &, I, const Labels &);
CanonicalSlider normalize_canonical_slider(const Word &, const std::vector<Factor> &,
                                           const std::vector<Factor> &,
                                           I max_image_letters = 100000,
                                           I max_operations = 10000000);
J encode_split(const LayeredSplit &);
J encode_canonical(const CanonicalSlider &);
} // namespace zipper_verify::typed
