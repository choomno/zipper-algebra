#pragma once
#include "core.hpp"
#include "ribbons.hpp"
namespace zipper_verify {
J split_layered(const J &word, I n, const J &bottom_ports);
void check_split(const J &witness);
J normalize_canonical_slider(const J &word, const J &top_routes = J::array(),
                             const J &bottom_routes = J::array(),
                             I max_image_letters = 100000,
                             I max_operations = 10000000);
void check_canonical_slider(const J &witness);
J extract_slider(const J &word, const J &upper_labels, const J &lower_labels,
                 const J &upper_ports = J(), const J &lower_ports = J());
void check_extraction(const J &witness);
} // namespace zipper_verify
