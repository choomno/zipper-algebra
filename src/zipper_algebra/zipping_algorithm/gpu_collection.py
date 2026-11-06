"""Optional exact GPU summation of captured zipper Laurent contributions.

This prototype accelerates coefficient addition only. Sorting/grouping keys,
packing, transfers, reconstruction, and independent verification are timed.
The active zipper scheduler is untouched. Keys include an input root identity,
so coefficients from different literal input words can never cancel.

Coefficients use signed two's-complement vectors of unsigned 64-bit limbs.
For ``m`` inputs of absolute bit length at most ``b``, the allocated width is
at least ``b + ceil(log2(m)) + 2`` bits. Every partial group sum has absolute
value below ``m * 2**b``, so its mathematical signed value fits that width.
Unsigned limb addition computes the residue modulo the width; signed decoding
therefore gives the exact integer, even with cancellation or >64-bit inputs.
No floating-point arithmetic or silent CPU arithmetic fallback is used.

CUDA needs CuPy; OpenCL needs PyOpenCL. Imports and device initialization occur
only when their collector is constructed. OpenCL selects GPU devices only.
Constructor/compilation time is reported separately as ``initialization_seconds``;
``timings['total']`` includes the complete collect call and its verification.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import importlib
import struct
from time import perf_counter
from typing import Any, Iterable


Contribution = tuple[int, int, int, int]
Key = tuple[int, int, int]


class GPUUnavailableError(RuntimeError):
    """The requested GPU backend or device is unavailable."""


class CollectionVerificationError(RuntimeError):
    """A GPU result disagrees with independent Python integer arithmetic."""


class GPUInitializationError(RuntimeError):
    """A selected GPU could not compile or initialize the exact kernel."""


@dataclass(frozen=True)
class CollectionResult:
    output: tuple[Contribution, ...]
    timings: dict[str, float]
    metadata: dict[str, Any]


@dataclass(frozen=True)
class PackedContributions:
    """Host-sorted keys and little-endian coefficient/offset byte buffers."""

    records: tuple[Contribution, ...]
    keys: tuple[Key, ...]
    offsets: tuple[int, ...]
    coefficients: bytes
    offset_bytes: bytes
    limbs: int
    max_input_bits: int


def _normalize(records: Iterable[Contribution]) -> tuple[Contribution, ...]:
    result = []
    for row in records:
        if not isinstance(row, (tuple, list)) or len(row) != 4:
            raise ValueError("Each contribution must contain root, state, exponent, coefficient")
        if any(type(value) is not int for value in row):
            raise ValueError("Contribution fields must be Python integers")
        result.append((row[0], row[1], row[2], row[3]))
    return tuple(result)


def collect_cpu(records: Iterable[Contribution]) -> tuple[Contribution, ...]:
    """Independent exact Python-integer collection, in canonical key order."""
    sums: dict[Key, int] = {}
    for root, state, exponent, coefficient in _normalize(records):
        key = (root, state, exponent)
        sums[key] = sums.get(key, 0) + coefficient
    return tuple((*key, value) for key, value in sorted(sums.items()) if value)


def pack_contributions(
    records: Iterable[Contribution],
) -> tuple[PackedContributions, dict[str, float]]:
    """Prepare exact GPU buffers without adding coefficients on the host."""
    start = perf_counter()
    normalized = _normalize(records)
    validated = perf_counter()
    ordered = sorted(normalized, key=lambda row: row[:3])
    keys: list[Key] = []
    offsets: list[int] = []
    for index, row in enumerate(ordered):
        key = row[:3]
        if not keys or key != keys[-1]:
            keys.append(key)
            offsets.append(index)
    offsets.append(len(ordered))
    sorted_at = perf_counter()
    bits = max((abs(row[3]).bit_length() for row in ordered), default=0)
    required_bits = bits + max(0, len(ordered) - 1).bit_length() + 2
    limbs = max(1, (required_bits + 63) // 64)
    width = limbs * 8
    coefficient_bytes = b"".join(row[3].to_bytes(width, "little", signed=True) for row in ordered)
    offset_bytes = b"".join(struct.pack("<Q", value) for value in offsets)
    packed = PackedContributions(normalized, tuple(keys), tuple(offsets),
                                 coefficient_bytes, offset_bytes, limbs, bits)
    return packed, {
        "validation": validated - start,
        "sort": sorted_at - validated,
        "packing": perf_counter() - sorted_at,
    }


def unpack_sums(packed: PackedContributions, result_bytes: bytes) -> tuple[Contribution, ...]:
    """Decode signed exact sums, removing zero groups."""
    width = packed.limbs * 8
    if len(result_bytes) != len(packed.keys) * width:
        raise ValueError("GPU output buffer has the wrong size")
    result = []
    for index, key in enumerate(packed.keys):
        value = int.from_bytes(result_bytes[index * width:(index + 1) * width],
                               "little", signed=True)
        if value:
            result.append((*key, value))
    return tuple(result)


_BODY = r"""
    if (group >= groups) return;
    const u64 base = group * limbs;
    for (u64 j = 0; j < limbs; ++j) output[base + j] = 0;
    for (u64 row = offsets[group]; row < offsets[group + 1]; ++row) {
        u64 carry = 0;
        for (u64 j = 0; j < limbs; ++j) {
            const u64 old = output[base + j];
            const u64 value = input[row * limbs + j];
            const u64 first = old + value;
            const u64 carry_first = (first < old);
            const u64 second = first + carry;
            const u64 carry_second = (second < first);
            output[base + j] = second;
            carry = carry_first | carry_second;
        }
    }
"""

CUDA_SOURCE = r"""
typedef unsigned long long u64;
extern "C" __global__ void sum_groups(
    const u64 *input, const u64 *offsets, u64 *output,
    const u64 groups, const u64 limbs) {
    const u64 group = (u64)blockIdx.x * blockDim.x + threadIdx.x;
""" + _BODY + "\n}\n"

OPENCL_SOURCE = r"""
typedef ulong u64;
__kernel void sum_groups(
    __global const u64 *input, __global const u64 *offsets,
    __global u64 *output, const u64 groups, const u64 limbs) {
    const u64 group = get_global_id(0);
""" + _BODY + "\n}\n"


class _Collector:
    metadata: dict[str, Any]
    initialization_seconds: float

    def _execute(self, packed: PackedContributions) -> tuple[bytes, dict[str, float]]:
        raise NotImplementedError

    def collect(self, records: Iterable[Contribution], *, verify: bool = True) -> CollectionResult:
        if type(verify) is not bool:
            raise ValueError("verify must be a boolean")
        start = perf_counter()
        packed, timings = pack_contributions(records)
        if packed.keys:
            raw, gpu_timings = self._execute(packed)
        else:
            raw, gpu_timings = b"", {name: 0.0 for name in ("h2d", "reduce", "d2h")}
        timings.update(gpu_timings)
        before_decode = perf_counter()
        output = unpack_sums(packed, raw)
        timings["reconstruction"] = perf_counter() - before_decode
        before_verify = perf_counter()
        if verify and output != collect_cpu(packed.records):
            raise CollectionVerificationError("GPU collection differs from independent exact CPU collection")
        timings["verification"] = perf_counter() - before_verify
        timings["total"] = perf_counter() - start
        metadata = dict(self.metadata)
        metadata.update({
            "input_contributions": len(packed.records), "groups": len(packed.keys),
            "output_nonzero_groups": len(output), "limbs": packed.limbs,
            "max_input_bits": packed.max_input_bits,
            "h2d_bytes": len(packed.coefficients) + len(packed.offset_bytes),
            "d2h_bytes": len(raw), "verified": verify,
            "sort_location": "host", "arithmetic": "exact signed multi-limb integers",
            "initialization_seconds": self.initialization_seconds,
        })
        return CollectionResult(output, timings, metadata)


class CudaCollector(_Collector):
    """Exact CUDA collection using lazy optional CuPy and a retained kernel."""

    def __init__(self, device: int = 0):
        start = perf_counter()
        selected = False
        if type(device) is not int or device < 0:
            raise ValueError("device must be a nonnegative integer")
        try:
            self.cp: Any = importlib.import_module("cupy")
            self.np: Any = importlib.import_module("numpy")
            if self.cp.cuda.runtime.getDeviceCount() <= device:
                raise GPUUnavailableError(f"CUDA GPU index {device} is unavailable")
            self.device = self.cp.cuda.Device(device)
            selected = True
            with self.device:
                self.kernel = self.cp.RawKernel(CUDA_SOURCE, "sum_groups")
                self.kernel.compile()
                self.device.synchronize()
                props = self.cp.cuda.runtime.getDeviceProperties(device)
            name = props["name"]
            if isinstance(name, bytes):
                name = name.decode(errors="replace")
            self.metadata = {"backend": "cuda", "device": device, "name": name,
                             "cupy_version": self.cp.__version__,
                             "driver_version": self.cp.cuda.runtime.driverGetVersion(),
                             "runtime_version": self.cp.cuda.runtime.runtimeGetVersion(),
                             "kernel_source_sha256": hashlib.sha256(CUDA_SOURCE.encode()).hexdigest()}
        except GPUUnavailableError:
            raise
        except Exception as error:
            if selected:
                raise GPUInitializationError(f"CUDA kernel initialization failed: {error}") from error
            raise GPUUnavailableError(f"Cannot initialize CUDA GPU: {error}") from error
        self.initialization_seconds = perf_counter() - start

    def _execute(self, packed: PackedContributions) -> tuple[bytes, dict[str, float]]:
        cp, np = self.cp, self.np
        with self.device:
            start = perf_counter()
            coefficients = cp.asarray(np.frombuffer(packed.coefficients, dtype="<u8"))
            offsets = cp.asarray(np.frombuffer(packed.offset_bytes, dtype="<u8"))
            output = cp.empty(len(packed.keys) * packed.limbs, dtype=cp.uint64)
            self.device.synchronize()
            transferred = perf_counter()
            self.kernel(((len(packed.keys) + 127) // 128,), (128,),
                        (coefficients, offsets, output, np.uint64(len(packed.keys)),
                         np.uint64(packed.limbs)))
            self.device.synchronize()
            reduced = perf_counter()
            raw = cp.asnumpy(output).astype("<u8", copy=False).tobytes()
            self.device.synchronize()
            end = perf_counter()
        return raw, {"h2d": transferred - start, "reduce": reduced - transferred,
                     "d2h": end - reduced}


class OpenClCollector(_Collector):
    """Exact OpenCL collection; device indices enumerate GPU devices only."""

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
                raise GPUUnavailableError("Prototype requires a little-endian OpenCL GPU")
            selected = True
            self.context = self.cl.Context([self.device])
            self.queue = self.cl.CommandQueue(self.context)
            self.program = self.cl.Program(self.context, OPENCL_SOURCE).build()
            self.kernel = self.cl.Kernel(self.program, "sum_groups")
            self.queue.finish()
            self.metadata = {"backend": "opencl", "device": device,
                             "name": self.device.name.strip(),
                             "vendor": self.device.vendor.strip(),
                             "driver_version": self.device.driver_version,
                             "opencl_version": self.device.version,
                             "pyopencl_version": self.cl.VERSION_TEXT,
                             "global_memory_bytes": self.device.global_mem_size,
                             "compute_units": self.device.max_compute_units,
                             "kernel_source_sha256": hashlib.sha256(OPENCL_SOURCE.encode()).hexdigest()}
        except GPUUnavailableError:
            raise
        except Exception as error:
            if selected:
                raise GPUInitializationError(f"OpenCL kernel initialization failed: {error}") from error
            raise GPUUnavailableError(f"Cannot initialize OpenCL GPU: {error}") from error
        self.initialization_seconds = perf_counter() - start

    def _execute(self, packed: PackedContributions) -> tuple[bytes, dict[str, float]]:
        cl, np = self.cl, self.np
        start = perf_counter()
        flags = cl.mem_flags
        coefficients = cl.Buffer(self.context, flags.READ_ONLY, size=len(packed.coefficients))
        offsets = cl.Buffer(self.context, flags.READ_ONLY, size=len(packed.offset_bytes))
        output = cl.Buffer(self.context, flags.READ_WRITE, size=len(packed.keys) * packed.limbs * 8)
        cl.enqueue_copy(self.queue, coefficients, packed.coefficients, is_blocking=True)
        cl.enqueue_copy(self.queue, offsets, packed.offset_bytes, is_blocking=True)
        self.queue.finish()
        transferred = perf_counter()
        self.kernel(self.queue, (len(packed.keys),), None, coefficients, offsets, output,
                    np.uint64(len(packed.keys)), np.uint64(packed.limbs))
        self.queue.finish()
        reduced = perf_counter()
        destination = np.empty(len(packed.keys) * packed.limbs, dtype="<u8")
        cl.enqueue_copy(self.queue, destination, output, is_blocking=True)
        self.queue.finish()
        raw = destination.tobytes()
        end = perf_counter()
        return raw, {"h2d": transferred - start, "reduce": reduced - transferred,
                     "d2h": end - reduced}
