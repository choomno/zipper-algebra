#pragma once
#include "core.hpp"
namespace zipper_verify {
J reconstruct_nested(const J &labels, const J &counts,
                     const J &initial_ports = nullptr);
J recognize_nested_strict(const J &word, const J &labels,
                          const J &initial_ports = nullptr);
J recognize_nested(const J &word, const J &labels,
                   const J &initial_ports = nullptr,
                   I max_image_letters = 100000, I max_operations = 10000000);
void check_nested(const J &witness);
} // namespace zipper_verify
