// Exact single-thread CPU execution of the packed rewrite instruction format.
// The frozen Python/OpenCL implementation is the independent transport peer;
// mathematical coefficients remain arbitrary-width unsigned-limb bit vectors.
#include <cstdint>
#include <limits>
using i64 = std::int64_t;
using u64 = std::uint64_t;
using u32 = std::uint32_t;
static_assert(sizeof(i64) == 8 && sizeof(u64) == 8 && sizeof(u32) == 4);
static inline void atomic_or(u32 *target, u32 value) { *target |= value; }
extern "C" u64 zipper_rewrite_abi() { return 0x5A495050455201ULL; }


#define LMAX 9223372036854775807L
#define LMIN (-9223372036854775807L - 1L)
i64 checked_add(i64 a, i64 b, u32 *error) {
    if ((b > 0 && a > LMAX - b) || (b < 0 && a < LMIN - b)) {
        *error = 1; return 0;
    }
    return a + b;
}
i64 checked_sub(i64 a, i64 b, u32 *error) {
    if ((b < 0 && a > LMAX + b) || (b > 0 && a < LMIN + b)) {
        *error = 1; return 0;
    }
    return a - b;
}
u64 magnitude(i64 value) {
    return value < 0 ? (u64)(-(value + 1)) + 1UL : (u64)value;
}
i64 checked_mul(i64 a, i64 b, u32 *error) {
    const u32 negative = (a < 0) != (b < 0);
    const u64 aa = magnitude(a), bb = magnitude(b);
    const u64 bound = negative ? (1UL << 63) : (u64)LMAX;
    if (bb && aa > bound / bb) { *error = 1; return 0; }
    const u64 value = aa * bb;
    if (!negative) return (i64)value;
    return value == (1UL << 63) ? LMIN : -(i64)value;
}
extern "C" void zipper_evaluate_registers(
    const i64 *programs, const i64 *expressions,
    const i64 *cases, const i64 *parameters,
    i64 *values, u32 *status, const i64 case_count) {
    for (i64 c = 0; c < case_count; ++c) {
    const i64 program = cases[c * 10], begin = programs[program * 4];
    const i64 count = programs[program * 4 + 1], out = cases[c * 10 + 1];
    u32 error = 0;
    for (i64 j = 0; j < count; ++j) {
        const i64 at = (begin + j) * 4, op = expressions[at];
        const i64 left = expressions[at + 1], right = expressions[at + 2];
        i64 result = 0;
        if (op == 0) result = expressions[at + 3];
        else if (op == 1) result = parameters[c * 6 + expressions[at + 3]];
        else {
            const i64 a = values[out + left], b = values[out + right];
            if (op == 2) result = checked_add(a, b, &error);
            else if (op == 3) result = checked_sub(a, b, &error);
            else result = checked_mul(a, b, &error);
        }
        values[out + j] = result;
    }
    if (error) atomic_or(status + c, 1U);

    }
}
void emit_token(
    i64 *output, i64 base, i64 capacity, i64 *cursor,
    i64 kind, i64 first, i64 second, i64 mirror, i64 reverse, i64 bar,
    u32 *error) {
    if (*cursor >= capacity) { *error = 1; return; }
    const i64 target = (base + *cursor) * 6;
    output[target] = kind; output[target + 1] = first;
    output[target + 2] = second; output[target + 3] = mirror;
    output[target + 4] = reverse; output[target + 5] = bar;
    *cursor += 1;
}
extern "C" void zipper_expand_words(
    const i64 *cases, const i64 *values,
    const i64 *terms, const i64 *primitives,
    const i64 *words, const i64 *refs,
    const i64 *named_words, const i64 *jobs,
    i64 *output, u32 *status, const i64 job_count) {
    for (i64 job = 0; job < job_count; ++job) {
    const i64 c = jobs[job * 4], term = jobs[job * 4 + 1];
    const i64 base = jobs[job * 4 + 2], capacity = jobs[job * 4 + 3];
    const i64 cv = c * 10, vo = cases[cv + 1], old = cases[cv + 2];
    const i64 old_len = cases[cv + 3], prefix = cases[cv + 4];
    const i64 end = cases[cv + 5], n = cases[cv + 9];
    i64 cursor = 0;
    u32 error = 0;
    for (i64 p = 0; p < prefix; ++p) {
        const i64 src = (old + p) * 6;
        emit_token(output, base, capacity, &cursor, words[src], words[src + 1],
                   words[src + 2], words[src + 3], words[src + 4], words[src + 5], &error);
    }
    const i64 begin = terms[term * 4], count = terms[term * 4 + 1];
    for (i64 p = 0; p < count; ++p) {
        const i64 at = (begin + p) * 6, kind = primitives[at];
        const i64 ia = primitives[at + 1], ib = primitives[at + 2];
        const i64 ic = primitives[at + 3], id = primitives[at + 4];
        const i64 flags = primitives[at + 5];
        const i64 a = kind == 5 ? ia : values[vo + ia];
        const i64 b = values[vo + ib], cc = values[vo + ic], d = values[vo + id];
        if (kind == 4) {
            i64 start = b, mirror = !!(flags & 8);
            if (flags & 4) { start = n - start - a + 2; mirror = !mirror; }
            emit_token(output, base, capacity, &cursor, 1, a, start, mirror,
                       !!(flags & 16), !!(flags & 32), &error);
            continue;
        }
        i64 total = 1;
        if (kind == 2) total = checked_mul(a, cc, &error);
        else if (kind == 3) total = checked_mul(a, a - 1, &error) / 2;
        else if (kind == 5) total = refs[(cases[cv + 6] + a) * 2 + 1];
        i64 half_length = flags & 1 ? 1 : a - 1, half_offset = 0;
        for (i64 pos = 0; pos < total; ++pos) {
            const i64 index = flags & 1 ? total - pos - 1 : pos;
            i64 generator = 0;
            if (kind == 1) generator = checked_mul(a, b, &error);
            else if (kind == 2) {
                i64 block = index / a;
                if (flags & 64) block = cc - block - 1;
                generator = checked_mul(d, b + block + index % a, &error);
            } else if (kind == 3) {
                generator = b + half_offset;
                if (flags & 1) {
                    --half_offset;
                    if (half_offset < 0) { ++half_length; half_offset = half_length - 1; }
                } else {
                    ++half_offset;
                    if (half_offset == half_length) { --half_length; half_offset = 0; }
                }
            } else {
                generator = named_words[refs[(cases[cv + 6] + a) * 2] + index];
            }
            if (flags & 2) generator = -generator;
            if (flags & 4) {
                generator = generator > 0 ? n - generator : -(n + generator);
            }
            emit_token(output, base, capacity, &cursor, 0, generator, 0, 0, 0, 0, &error);
        }
    }
    for (i64 p = end; p < old_len; ++p) {
        const i64 src = (old + p) * 6;
        emit_token(output, base, capacity, &cursor, words[src], words[src + 1],
                   words[src + 2], words[src + 3], words[src + 4], words[src + 5], &error);
    }
    if (cursor != capacity || error) atomic_or(status + c, 2U);

    }
}
extern "C" void zipper_propagate_coefficients(
    const i64 *cases, const i64 *values,
    const i64 *factors, const i64 *old_q,
    const u64 *input, const i64 *jobs,
    i64 *new_q, u64 *output,
    u32 *status, const i64 limbs, const i64 job_count) {
    for (i64 job = 0; job < job_count; ++job) {
    const i64 c = jobs[job * 4], factor = jobs[job * 4 + 1];
    const i64 source = jobs[job * 4 + 2], target = jobs[job * 4 + 3];
    const i64 power = values[cases[c * 10 + 1] + factors[factor * 2]];
    const i64 sign = factors[factor * 2 + 1];
    u32 error = 0;
    new_q[target] = sign == 0 ? power : checked_add(old_q[source], power, &error);
    u64 carry = sign < 0 ? 1UL : 0UL;
    for (i64 j = 0; j < limbs; ++j) {
        const u64 value = input[source * limbs + j];
        if (sign == 0) output[target * limbs + j] = 0;
        else if (sign > 0) output[target * limbs + j] = value;
        else {
            const u64 inverted = ~value, result = inverted + carry;
            carry = result < inverted;
            output[target * limbs + j] = result;
        }
    }
    if (error) atomic_or(status + c, 4U);

    }
}
