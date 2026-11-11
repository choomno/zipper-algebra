"""Packed OpenCL contextual catalog expansion with exact signed coefficients.

The CPU admits guarded literal occurrences and computes allocation lengths.
Reusable programs retain symbolic integer expressions and repeat descriptors;
the device evaluates their registers and emits successor words itself. No
CPU-expanded RHS word is an input buffer. Coefficients in the accepted bundled
catalog are +/-1 monomials, so arbitrary-width input integer propagation needs
exact copy/negation plus checked exponent shifts. Final duplicate collection
and geometric scheduler control remain on the CPU.
"""

from __future__ import annotations

from array import array
from collections import Counter
from dataclasses import dataclass
import hashlib
import importlib
import json
import sys
from time import perf_counter
from typing import Any, Iterable

from .gpu_collection import GPUInitializationError, GPUUnavailableError
from .gpu_rewrite_program import CompiledRule, compile_rule
from .gpu_rewrite_host import AdmittedCase, admit, expand_cpu
from .interpreter import Branch, RunState, ReplayError


OPENCL_SOURCE = r"""
typedef unsigned long ulong64;
#define LMAX 9223372036854775807L
#define LMIN (-9223372036854775807L - 1L)
long checked_add(long a, long b, __private uint *error) {
    if ((b > 0 && a > LMAX - b) || (b < 0 && a < LMIN - b)) {
        *error = 1; return 0;
    }
    return a + b;
}
long checked_sub(long a, long b, __private uint *error) {
    if ((b < 0 && a > LMAX + b) || (b > 0 && a < LMIN + b)) {
        *error = 1; return 0;
    }
    return a - b;
}
ulong magnitude(long value) {
    return value < 0 ? (ulong)(-(value + 1)) + 1UL : (ulong)value;
}
long checked_mul(long a, long b, __private uint *error) {
    const uint negative = (a < 0) != (b < 0);
    const ulong aa = magnitude(a), bb = magnitude(b);
    const ulong bound = negative ? (1UL << 63) : (ulong)LMAX;
    if (bb && aa > bound / bb) { *error = 1; return 0; }
    const ulong value = aa * bb;
    if (!negative) return (long)value;
    return value == (1UL << 63) ? LMIN : -(long)value;
}
__kernel void evaluate_registers(
    __global const long *programs, __global const long *expressions,
    __global const long *cases, __global const long *parameters,
    __global long *values, __global uint *status, const long case_count) {
    const long c = get_global_id(0);
    if (c >= case_count) return;
    const long program = cases[c * 10], begin = programs[program * 4];
    const long count = programs[program * 4 + 1], out = cases[c * 10 + 1];
    uint error = 0;
    for (long j = 0; j < count; ++j) {
        const long at = (begin + j) * 4, op = expressions[at];
        const long left = expressions[at + 1], right = expressions[at + 2];
        long result = 0;
        if (op == 0) result = expressions[at + 3];
        else if (op == 1) result = parameters[c * 6 + expressions[at + 3]];
        else {
            const long a = values[out + left], b = values[out + right];
            if (op == 2) result = checked_add(a, b, &error);
            else if (op == 3) result = checked_sub(a, b, &error);
            else result = checked_mul(a, b, &error);
        }
        values[out + j] = result;
    }
    if (error) atomic_or(status + c, 1U);
}
void emit_token(
    __global long *output, long base, long capacity, __private long *cursor,
    long kind, long first, long second, long mirror, long reverse, long bar,
    __private uint *error) {
    if (*cursor >= capacity) { *error = 1; return; }
    const long target = (base + *cursor) * 6;
    output[target] = kind; output[target + 1] = first;
    output[target + 2] = second; output[target + 3] = mirror;
    output[target + 4] = reverse; output[target + 5] = bar;
    *cursor += 1;
}
__kernel void expand_words(
    __global const long *cases, __global const long *values,
    __global const long *terms, __global const long *primitives,
    __global const long *words, __global const long *refs,
    __global const long *named_words, __global const long *jobs,
    __global long *output, __global uint *status, const long job_count) {
    const long job = get_global_id(0);
    if (job >= job_count) return;
    const long c = jobs[job * 4], term = jobs[job * 4 + 1];
    const long base = jobs[job * 4 + 2], capacity = jobs[job * 4 + 3];
    const long cv = c * 10, vo = cases[cv + 1], old = cases[cv + 2];
    const long old_len = cases[cv + 3], prefix = cases[cv + 4];
    const long end = cases[cv + 5], n = cases[cv + 9];
    long cursor = 0;
    uint error = 0;
    for (long p = 0; p < prefix; ++p) {
        const long src = (old + p) * 6;
        emit_token(output, base, capacity, &cursor, words[src], words[src + 1],
                   words[src + 2], words[src + 3], words[src + 4], words[src + 5], &error);
    }
    const long begin = terms[term * 4], count = terms[term * 4 + 1];
    for (long p = 0; p < count; ++p) {
        const long at = (begin + p) * 6, kind = primitives[at];
        const long ia = primitives[at + 1], ib = primitives[at + 2];
        const long ic = primitives[at + 3], id = primitives[at + 4];
        const long flags = primitives[at + 5];
        const long a = kind == 5 ? ia : values[vo + ia];
        const long b = values[vo + ib], cc = values[vo + ic], d = values[vo + id];
        if (kind == 4) {
            long start = b, mirror = !!(flags & 8);
            if (flags & 4) { start = n - start - a + 2; mirror = !mirror; }
            emit_token(output, base, capacity, &cursor, 1, a, start, mirror,
                       !!(flags & 16), !!(flags & 32), &error);
            continue;
        }
        long total = 1;
        if (kind == 2) total = checked_mul(a, cc, &error);
        else if (kind == 3) total = checked_mul(a, a - 1, &error) / 2;
        else if (kind == 5) total = refs[(cases[cv + 6] + a) * 2 + 1];
        long half_length = flags & 1 ? 1 : a - 1, half_offset = 0;
        for (long pos = 0; pos < total; ++pos) {
            const long index = flags & 1 ? total - pos - 1 : pos;
            long generator = 0;
            if (kind == 1) generator = checked_mul(a, b, &error);
            else if (kind == 2) {
                long block = index / a;
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
    for (long p = end; p < old_len; ++p) {
        const long src = (old + p) * 6;
        emit_token(output, base, capacity, &cursor, words[src], words[src + 1],
                   words[src + 2], words[src + 3], words[src + 4], words[src + 5], &error);
    }
    if (cursor != capacity || error) atomic_or(status + c, 2U);
}
__kernel void propagate_coefficients(
    __global const long *cases, __global const long *values,
    __global const long *factors, __global const long *old_q,
    __global const ulong *input, __global const long *jobs,
    __global long *new_q, __global ulong *output,
    __global uint *status, const long limbs, const long job_count) {
    const long job = get_global_id(0);
    if (job >= job_count) return;
    const long c = jobs[job * 4], factor = jobs[job * 4 + 1];
    const long source = jobs[job * 4 + 2], target = jobs[job * 4 + 3];
    const long power = values[cases[c * 10 + 1] + factors[factor * 2]];
    const long sign = factors[factor * 2 + 1];
    uint error = 0;
    new_q[target] = sign == 0 ? power : checked_add(old_q[source], power, &error);
    ulong carry = sign < 0 ? 1UL : 0UL;
    for (long j = 0; j < limbs; ++j) {
        const ulong value = input[source * limbs + j];
        if (sign == 0) output[target * limbs + j] = 0;
        else if (sign > 0) output[target * limbs + j] = value;
        else {
            const ulong inverted = ~value, result = inverted + carry;
            carry = result < inverted;
            output[target * limbs + j] = result;
        }
    }
    if (error) atomic_or(status + c, 4U);
}
"""


@dataclass(frozen=True)
class ExpansionResult:
    output: tuple[dict[str, Any], ...]
    timings: dict[str, float]
    metadata: dict[str, Any]


@dataclass
class _Packed:
    admitted: list[AdmittedCase]
    arrays: dict[str, bytes]
    word_jobs: list[tuple[int, int, int, int]]
    coefficient_jobs: list[tuple[int, int, int, int]]
    case_word_jobs: list[list[int]]
    term_factor_indices: list[list[int]]
    word_tokens: int
    coefficient_records: int
    expression_values: int
    limbs: int
    source_bits: int


def _i64(values: Iterable[int]) -> bytes:
    data = array("q", values)
    if sys.byteorder != "little":
        data.byteswap()
    return data.tobytes()


def _decode_i64(data: bytes) -> list[int]:
    result = array("q")
    result.frombytes(data)
    if sys.byteorder != "little":
        result.byteswap()
    return list(result)


def _token(token: Any) -> tuple[int, int, int, int, int, int]:
    if type(token) is int:
        return (0, token, 0, 0, 0, 0)
    return (1, token[1], token[2], int(token[3]), int(token[4]), int(token[5]))


class OpenClRewriteExpander:
    """Retain compiled symbolic templates and execute genuine packed expansion."""

    def __init__(self, device: int = 0):
        start = perf_counter()
        selected = False
        if type(device) is not int or device < 0:
            raise ValueError("device must be a nonnegative integer")
        try:
            self.cl: Any = importlib.import_module("pyopencl")
            self.np: Any = importlib.import_module("numpy")
            devices = [item for platform in self.cl.get_platforms()
                       for item in platform.get_devices(device_type=self.cl.device_type.GPU)]
            if device >= len(devices):
                raise GPUUnavailableError(f"OpenCL GPU index {device} is unavailable")
            self.device = devices[device]
            if not self.device.endian_little:
                raise GPUUnavailableError("Packed prototype requires a little-endian GPU")
            selected = True
            self.context = self.cl.Context([self.device])
            self.queue = self.cl.CommandQueue(self.context)
            self.program = self.cl.Program(self.context, OPENCL_SOURCE).build()
            self.evaluate_kernel = self.cl.Kernel(self.program, "evaluate_registers")
            self.word_kernel = self.cl.Kernel(self.program, "expand_words")
            self.coefficient_kernel = self.cl.Kernel(self.program, "propagate_coefficients")
            self.queue.finish()
            self.metadata = {
                "backend": "opencl", "device": device,
                "name": self.device.name.strip(), "vendor": self.device.vendor.strip(),
                "driver_version": self.device.driver_version,
                "opencl_version": self.device.version,
                "pyopencl_version": self.cl.VERSION_TEXT,
                "kernel_source_sha256": hashlib.sha256(OPENCL_SOURCE.encode()).hexdigest(),
                "arithmetic": "exact signed multi-limb integers and checked i64 coordinates",
                "successor_word_generation": "GPU symbolic program execution",
                "host_work": "guard/LHS admission, allocation sizes, final duplicate collection and controls",
            }
        except GPUUnavailableError:
            raise
        except Exception as error:
            if selected:
                raise GPUInitializationError(f"OpenCL rewrite kernel initialization failed: {error}") from error
            raise GPUUnavailableError(f"Cannot initialize OpenCL GPU: {error}") from error
        self.initialization_seconds = perf_counter() - start
        self._templates: list[CompiledRule] = []
        self._template_keys: dict[str, int] = {}
        self._template_arrays: dict[str, bytes] = {}
        self._template_buffers: dict[str, Any] = {}
        self._program_term_indices: list[list[int]] = []
        self._term_factor_indices: list[list[int]] = []
        self._templates_dirty = False

    def _template(self, rule: dict[str, Any]) -> tuple[int, CompiledRule]:
        key = json.dumps(rule["rhs"], sort_keys=True, separators=(",", ":"), allow_nan=False)
        if key not in self._template_keys:
            compiled = compile_rule(rule)
            self._template_keys[key] = len(self._templates)
            self._templates.append(compiled)
            self._templates_dirty = True
        index = self._template_keys[key]
        return index, self._templates[index]

    def _rebuild_templates(self) -> None:
        programs: list[int] = []
        expressions: list[int] = []
        primitives: list[int] = []
        terms: list[int] = []
        factors: list[int] = []
        self._program_term_indices = []
        self._term_factor_indices = []
        for program in self._templates:
            programs.extend((len(expressions) // 4, len(program.expressions),
                             len(terms) // 4, len(program.terms)))
            expressions.extend(value for row in program.expressions for value in row)
            term_indices = []
            for term in program.terms:
                term_indices.append(len(terms) // 4)
                start = len(factors) // 2
                factor_indices = list(range(start, start + len(term.coefficients)))
                self._term_factor_indices.append(factor_indices)
                terms.extend((len(primitives) // 6, len(term.primitives),
                              start, len(term.coefficients)))
                primitives.extend(value for row in term.primitives for value in row)
                factors.extend(value for row in term.coefficients for value in row)
            self._program_term_indices.append(term_indices)
        self._template_arrays = {
            "programs": _i64(programs), "expressions": _i64(expressions),
            "primitives": _i64(primitives), "terms": _i64(terms),
            "factors": _i64(factors),
        }

    def _pack(self, raw_cases: Iterable[dict[str, Any]]) -> _Packed:
        raw_cases = list(raw_cases)
        seen = set()
        admitted: list[AdmittedCase] = []
        programs: list[int] = []
        for raw_case in raw_cases:
            if (not isinstance(raw_case, dict)
                    or any(type(raw_case.get(name)) is not int for name in ("root", "id"))):
                raise ReplayError("packed cases require integer root and case identities")
            identity = (raw_case["root"], raw_case["id"])
            if identity in seen:
                raise ReplayError("duplicate root/case identity in packed rewrite batch")
            seen.add(identity)
            program_index, program = self._template(raw_case["rule"])
            admitted.append(admit(raw_case, program))
            programs.append(program_index)
        if self._templates_dirty:
            self._rebuild_templates()
        bits = max((abs(c).bit_length() for case in admitted for _, c in case.branch.coefficient), default=0)
        limbs = max(1, (bits + 2 + 63) // 64)
        width = limbs * 8
        cases: list[int] = []
        parameters: list[int] = []
        words: list[int] = []
        refs: list[int] = []
        named_words: list[int] = []
        old_q: list[int] = []
        coefficients = bytearray()
        word_jobs: list[tuple[int, int, int, int]] = []
        coefficient_jobs: list[tuple[int, int, int, int]] = []
        case_word_jobs: list[list[int]] = []
        expression_values = word_tokens = 0
        for index, (case, program_index) in enumerate(zip(admitted, programs, strict=True)):
            program = case.program
            word_offset = len(words) // 6
            words.extend(value for token in case.branch.word for value in _token(token))
            ref_offset = len(refs) // 2
            for name in program.refs:
                reference = tuple((case.request.words or {})[name])
                refs.extend((len(named_words), len(reference)))
                named_words.extend(reference)
            coefficient_start = len(old_q)
            for power, value in case.branch.coefficient:
                old_q.append(power)
                coefficients.extend(value.to_bytes(width, "little", signed=True))
            cases.extend((program_index, expression_values, word_offset,
                          len(case.branch.word), case.request.offset, case.end,
                          ref_offset, coefficient_start, len(case.branch.coefficient),
                          case.raw["n"]))
            expression_values += len(program.expressions)
            used_slots = {node[3] for node in program.expressions if node[0] == 1}
            parameters.extend(case.request.parameters.get(name, 0) if slot in used_slots else 0
                              for slot, name in enumerate(("n", "i", "j", "k", "l", "eps")))
            current_word_jobs = []
            for term_index, capacity in enumerate(case.lengths):
                global_term = self._program_term_indices[program_index][term_index]
                word_job = len(word_jobs)
                current_word_jobs.append(word_job)
                word_jobs.append((index, global_term, word_tokens, capacity))
                word_tokens += capacity
                for factor in self._term_factor_indices[global_term]:
                    for source in range(coefficient_start, len(old_q)):
                        coefficient_jobs.append((index, factor, source, len(coefficient_jobs)))
            case_word_jobs.append(current_word_jobs)
        arrays = {
            "cases": _i64(cases), "parameters": _i64(parameters), "words": _i64(words),
            "refs": _i64(refs), "named_words": _i64(named_words), "old_q": _i64(old_q),
            "coefficients": bytes(coefficients),
            "word_jobs": _i64(value for job in word_jobs for value in job),
            "coefficient_jobs": _i64(value for job in coefficient_jobs for value in job),
        }
        return _Packed(admitted, arrays, word_jobs, coefficient_jobs, case_word_jobs,
                       self._term_factor_indices, word_tokens, len(coefficient_jobs),
                       expression_values, limbs, bits)

    def _execute(self, packed: _Packed) -> tuple[dict[str, bytes], dict[str, float], int]:
        cl, np = self.cl, self.np
        start = perf_counter()
        upload_bytes = 0
        if self._templates_dirty:
            for name, payload in self._template_arrays.items():
                data = payload or bytes(8)
                self._template_buffers[name] = cl.Buffer(
                    self.context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                    hostbuf=data)
                upload_bytes += len(data)
            self._templates_dirty = False
        buffers = {}
        for name, payload in packed.arrays.items():
            data = payload or bytes(8)
            buffers[name] = cl.Buffer(self.context, cl.mem_flags.READ_ONLY, size=len(data))
            cl.enqueue_copy(self.queue, buffers[name], data, is_blocking=True)
            upload_bytes += len(data)
        output_sizes = {
            "word_output": packed.word_tokens * 6 * 8,
            "coefficient_output": packed.coefficient_records * packed.limbs * 8,
            "q_output": packed.coefficient_records * 8,
            "status": len(packed.admitted) * 4,
        }
        output = {name: cl.Buffer(self.context, cl.mem_flags.READ_WRITE,
                                 size=max(size, 8))
                  for name, size in output_sizes.items()}
        cl.enqueue_copy(self.queue, output["status"],
                        bytes(max(output_sizes["status"], 8)), is_blocking=True)
        upload_bytes += max(output_sizes["status"], 8)
        values = cl.Buffer(self.context, cl.mem_flags.READ_WRITE,
                           size=max(packed.expression_values * 8, 8))
        self.queue.finish()
        transferred = perf_counter()
        templates = self._template_buffers
        self.evaluate_kernel(self.queue, (len(packed.admitted),), None,
                             templates["programs"], templates["expressions"],
                             buffers["cases"], buffers["parameters"], values,
                             output["status"], np.int64(len(packed.admitted)))
        self.queue.finish()
        evaluated = perf_counter()
        if packed.word_jobs:
            self.word_kernel(self.queue, (len(packed.word_jobs),), None,
                             buffers["cases"], values, templates["terms"],
                             templates["primitives"], buffers["words"], buffers["refs"],
                             buffers["named_words"], buffers["word_jobs"],
                             output["word_output"], output["status"], np.int64(len(packed.word_jobs)))
        self.queue.finish()
        expanded = perf_counter()
        if packed.coefficient_jobs:
            self.coefficient_kernel(
                self.queue, (len(packed.coefficient_jobs),), None, buffers["cases"],
                values, templates["factors"], buffers["old_q"], buffers["coefficients"],
                buffers["coefficient_jobs"], output["q_output"], output["coefficient_output"],
                output["status"], np.int64(packed.limbs), np.int64(len(packed.coefficient_jobs)))
        self.queue.finish()
        propagated = perf_counter()
        result = {}
        for name, size in output_sizes.items():
            if not size:
                result[name] = b""
                continue
            destination = np.empty(size, dtype=np.uint8)
            cl.enqueue_copy(self.queue, destination, output[name], is_blocking=True)
            result[name] = destination.tobytes()
        self.queue.finish()
        returned = perf_counter()
        statuses = np.frombuffer(result["status"], dtype="<u4")
        if any(statuses):
            raise ReplayError("GPU rewrite checked arithmetic/allocation failure: " + repr(statuses.tolist()))
        timings = {
            "h2d_and_allocation": transferred - start, "evaluate": evaluated - transferred,
            "word_expansion": expanded - evaluated,
            "coefficient_propagation": propagated - expanded, "d2h": returned - propagated,
        }
        return result, timings, upload_bytes

    def _reconstruct(self, packed: _Packed, raw: dict[str, bytes]) -> tuple[dict[str, Any], ...]:
        tokens = _decode_i64(raw["word_output"])
        q_values = _decode_i64(raw["q_output"])
        width = packed.limbs * 8
        coefficient_values = [
            int.from_bytes(raw["coefficient_output"][j * width:(j + 1) * width],
                           "little", signed=True)
            for j in range(packed.coefficient_records)]
        result = []
        coefficient_job = 0
        for case_index, case in enumerate(packed.admitted):
            # Keep raw (q,word) insertion order until zero factors are removed;
            # the interpreter groups distinct words only after this cleaning.
            factors: dict[tuple[int, tuple[Any, ...]], int] = {}
            products: dict[tuple[int, tuple[Any, ...]], Counter[int]] = {}
            for local_term, word_job_index in enumerate(packed.case_word_jobs[case_index]):
                _, _, offset, length = packed.word_jobs[word_job_index]
                word_tokens: list[Any] = []
                for pos in range(offset, offset + length):
                    row = tokens[pos * 6:(pos + 1) * 6]
                    if row[0] == 0:
                        word_tokens.append(row[1])
                    else:
                        word_tokens.append(("slider", row[1], row[2], bool(row[3]), bool(row[4]), bool(row[5])))
                word_tuple = tuple(word_tokens)
                for _, factor in case.program.terms[local_term].coefficients:
                    q_rhs = (q_values[coefficient_job] if factor == 0
                             else q_values[coefficient_job] - case.branch.coefficient[0][0])
                    key = q_rhs, word_tuple
                    factors[key] = factors.get(key, 0) + factor
                    product = products.setdefault(key, Counter())
                    for _ in case.branch.coefficient:
                        product[q_values[coefficient_job]] += coefficient_values[coefficient_job]
                        coefficient_job += 1
            grouped: dict[tuple[Any, ...], Counter[int]] = {}
            for key, factor in factors.items():
                if factor:
                    grouped.setdefault(key[1], Counter()).update(products[key])
            controls = case.request._controls_json
            if controls is not None and len(controls) != len(grouped):
                raise ReplayError("one control update is required per distinct RHS word")
            branches = []
            for index, (following_word, coefficients) in enumerate(grouped.items()):
                coefficient = tuple((q, c) for q, c in sorted(coefficients.items()) if c)
                control = case.branch.control_json if controls is None else controls[index]
                branches.append(Branch(following_word, coefficient, control))
            state = RunState(case.raw["n"], tuple(branches))
            case.limits.check(state)
            result.append({"root": case.raw["root"], "id": case.raw["id"],
                           "branches": [branch.to_json() for branch in branches]})
        if coefficient_job != packed.coefficient_records:
            raise ReplayError("GPU coefficient reconstruction did not consume every emitted monomial")
        return tuple(result)

    def expand(self, cases: Iterable[dict[str, Any]], *, verify: bool = True) -> ExpansionResult:
        if type(verify) is not bool:
            raise ValueError("verify must be a boolean")
        started = perf_counter()
        raw_cases = list(cases)
        templates_before = len(self._templates)
        packed = self._pack(raw_cases)
        packed_at = perf_counter()
        template_upload = bool(packed.admitted and self._templates_dirty)
        template_upload_bytes = (sum(max(len(data), 8) for data in self._template_arrays.values())
                                 if template_upload else 0)
        if packed.admitted:
            raw, timings, h2d_bytes = self._execute(packed)
            output = self._reconstruct(packed, raw)
            d2h_bytes = sum(map(len, raw.values()))
        else:
            timings = {key: 0.0 for key in (
                "h2d_and_allocation", "evaluate", "word_expansion", "coefficient_propagation", "d2h")}
            output, h2d_bytes, d2h_bytes = (), 0, 0
        reconstructed = perf_counter()
        device_total = sum(timings.values())
        timings["packing_and_admission"] = packed_at - started
        timings["reconstruction"] = max(0.0, reconstructed - packed_at - device_total)
        if verify and output != expand_cpu(raw_cases):
            raise ReplayError("GPU rewrite differs from independent contextual CPU replay")
        checked = perf_counter()
        timings["verification"] = checked - reconstructed
        timings["total"] = checked - started
        metadata = dict(self.metadata, input_cases=len(packed.admitted),
                        template_count=len(self._templates),
                        new_templates=len(self._templates) - templates_before,
                        emitted_word_tokens=packed.word_tokens,
                        emitted_coefficient_records=packed.coefficient_records,
                        output_word_tokens=packed.word_tokens,
                        output_coefficient_monomials=packed.coefficient_records,
                        emission_count_scope="raw emitted buffers before host collection; duplicates and zero-factor ordering records included",
                        template_upload_happened=template_upload,
                        template_upload_bytes=template_upload_bytes,
                        coefficient_limbs=packed.limbs, max_input_coefficient_bits=packed.source_bits,
                        h2d_bytes=h2d_bytes, d2h_bytes=d2h_bytes, verified=verify,
                        h2d_byte_measure="logical payload including template initialization, padded empty inputs and zeroed status buffer",
                        initialization_seconds=self.initialization_seconds,
                        unsupported_factor_policy="reject factors outside zero/plus-or-minus one")
        return ExpansionResult(output, timings, metadata)
