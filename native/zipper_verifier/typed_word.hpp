#pragma once
// Independent verifier words. JSON is decoded at the protocol/catalog boundary.
#include "core.hpp"
#include <variant>

namespace zipper_verify::typed {
struct Slider {
    I size, start;
    bool mirrored = false, reversed = false, barred = false;
    bool operator==(const Slider &other) const {
        return size == other.size && start == other.start && mirrored == other.mirrored &&
               reversed == other.reversed && barred == other.barred;
    }
    bool operator!=(const Slider &other) const {
        return !(*this == other);
    }
};
using Token = std::variant<I, Slider>;
using Word = std::vector<Token>;
Token decode_token(const J &token);
Word decode_word(const J &word);
J encode_token(const Token &token);
J encode_word(const Word &word);
bool is_int(const Token &token);
I as_i(const Token &token);
bool is_slider(const Token &token);
I length(const Word &word);
Slider slider(I size, I start, bool mirrored = false, bool reversed = false, bool barred = false);
struct WordView {
    const Word &value;
    WordView(const Word &word) : value(word) {}
};
Word slice(const Word &word, I begin, I end = std::numeric_limits<I>::max());
void append(Word &target, const Word &word);
Word concat(std::initializer_list<WordView> words);
Word inverse(const Word &word);
Word bar(const Word &word);
Word reverse_word(const Word &word);
Word phi(const Word &word, I rank);
Word shifted(const Word &word, I shift);
void validate_word(const Word &word, I rank);
} // namespace zipper_verify::typed
