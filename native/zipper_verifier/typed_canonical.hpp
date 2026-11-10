#pragma once
#include "typed_state.hpp"
namespace zipper_verify::typed {
// Exact sorted, ensure_ascii JSON identity; exact coefficient scalars are numbers.
std::string canonical_work(const Work &);
std::string continuation_key(const Work &);        // omit scalar, depth and branch scalar
std::string work_collection_key(const Work &);     // omit top-level scalar
std::string output_collection_key(const Output &); // omit scalar
} // namespace zipper_verify::typed
