#include "typed_word.hpp"
#include <algorithm>

namespace zipper_verify::typed {
Token decode_token(const J &token) {
    if (zipper_verify::is_int(token))
        return zipper_verify::as_i(token);
    if (!zipper_verify::is_slider(token))
        throw ReplayError("word contains an invalid braid or slider token");
    return Slider{zipper_verify::as_i(token[1]), zipper_verify::as_i(token[2]),
                  token[3].get<bool>(), token[4].get<bool>(), token[5].get<bool>()};
}
Word decode_word(const J &word) {
    if (!word.is_array())
        throw ReplayError("word must be an array");
    Word result;
    result.reserve(word.size());
    for (const J &token : word)
        result.push_back(decode_token(token));
    return result;
}
J encode_token(const Token &token) {
    if (const auto *integer = std::get_if<I>(&token))
        return *integer;
    const Slider &coupon = std::get<Slider>(token);
    return J::array(
        {"slider", coupon.size, coupon.start, coupon.mirrored, coupon.reversed, coupon.barred});
}
J encode_word(const Word &word) {
    J result = J::array();
    for (const Token &token : word)
        result.push_back(encode_token(token));
    return result;
}
bool is_int(const Token &token) {
    return std::holds_alternative<I>(token);
}
I as_i(const Token &token) {
    if (const auto *integer = std::get_if<I>(&token))
        return *integer;
    throw ReplayError("integer required");
}
bool is_slider(const Token &token) {
    return std::holds_alternative<Slider>(token);
}
I length(const Word &word) {
    if (word.size() > static_cast<std::size_t>(std::numeric_limits<I>::max()))
        throw LimitExceeded("signed integer coordinate overflow");
    return static_cast<I>(word.size());
}
Slider slider(I size, I start, bool mirrored, bool reversed, bool barred) {
    return {size, start, mirrored, reversed, barred};
}
Word slice(const Word &word, I begin, I end) {
    const I size = length(word);
    auto position = [size](I x) {
        if (x < 0)
            return x < -size ? I(0) : size + x;
        return std::min(x, size);
    };
    begin = position(begin);
    end = position(end);
    if (begin >= end)
        return {};
    return Word(word.begin() + begin, word.begin() + end);
}
void append(Word &target, const Word &word) {
    // Explicitly support aliasing: insertion through invalidated iterators is undefined.
    if (&target == &word) {
        Word copy = word;
        target.insert(target.end(), copy.begin(), copy.end());
    } else
        target.insert(target.end(), word.begin(), word.end());
}
Word concat(std::initializer_list<WordView> words) {
    Word result;
    for (const WordView &word : words)
        append(result, word.value);
    return result;
}
Word inverse(const Word &word) {
    Word result;
    result.reserve(word.size());
    for (auto it = word.rbegin(); it != word.rend(); ++it) {
        if (!is_int(*it))
            throw ReplayError("inverse requires a pure braid word");
        result.push_back(neg(as_i(*it)));
    }
    return result;
}
Word bar(const Word &word) {
    Word result;
    result.reserve(word.size());
    for (Token token : word) {
        if (auto *integer = std::get_if<I>(&token))
            *integer = neg(*integer);
        else {
            auto &coupon = std::get<Slider>(token);
            coupon.barred = !coupon.barred;
        }
        result.push_back(std::move(token));
    }
    return result;
}
Word reverse_word(const Word &word) {
    return Word(word.rbegin(), word.rend());
}
Word phi(const Word &word, I rank) {
    Word result;
    result.reserve(word.size());
    for (auto it = word.rbegin(); it != word.rend(); ++it) {
        Token token = *it;
        if (auto *integer = std::get_if<I>(&token)) {
            const I g = *integer;
            *integer = checked_mul(g > 0 ? 1 : -1, checked_sub(rank, abs_i(g)));
        } else {
            auto &coupon = std::get<Slider>(token);
            coupon.start =
                checked_add(checked_sub(checked_sub(rank, coupon.size), coupon.start), 2);
        }
        result.push_back(std::move(token));
    }
    return result;
}
Word shifted(const Word &word, I shift) {
    Word result;
    result.reserve(word.size());
    for (Token token : word) {
        if (auto *integer = std::get_if<I>(&token)) {
            const I g = *integer;
            *integer = checked_mul(g > 0 ? 1 : -1, checked_add(abs_i(g), shift));
        } else {
            auto &coupon = std::get<Slider>(token);
            coupon.start = checked_add(coupon.start, shift);
        }
        result.push_back(std::move(token));
    }
    return result;
}
void validate_word(const Word &word, I rank) {
    for (const Token &token : word) {
        bool valid;
        if (const auto *integer = std::get_if<I>(&token))
            valid = *integer != 0 && abs_i(*integer) < rank;
        else {
            const auto &coupon = std::get<Slider>(token);
            valid = coupon.size >= 1 && coupon.size <= rank && coupon.start >= 1 &&
                    coupon.start <= checked_add(checked_sub(rank, coupon.size), 1);
        }
        if (!valid)
            throw ReplayError("word contains an invalid braid or slider token");
    }
}
} // namespace zipper_verify::typed
