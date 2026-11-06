"""Exact packing and independent CPU oracle; actual GPU checks when available."""

import random
import struct
from types import SimpleNamespace

import pytest

from zipper_algebra.zipping_algorithm import gpu_collection as gpu


def emulate_limbs(packed):
    """Unsigned-limb carry recurrence, explicitly not a hardware GPU test."""
    values = struct.unpack("<" + "Q" * (len(packed.coefficients) // 8), packed.coefficients)
    output = []
    mask = (1 << 64) - 1
    for start, stop in zip(packed.offsets, packed.offsets[1:]):
        accumulator = [0] * packed.limbs
        for row in range(start, stop):
            carry = 0
            for limb in range(packed.limbs):
                old = accumulator[limb]
                first = (old + values[row * packed.limbs + limb]) & mask
                carry_first = first < old
                second = (first + carry) & mask
                carry_second = second < first
                accumulator[limb] = second
                carry = carry_first | carry_second
        output.extend(accumulator)
    return b"".join(struct.pack("<Q", value) for value in output)


def test_root_state_and_exponent_all_belong_to_the_identity():
    records = [(0, 0, 0, 7), (1, 0, 0, -7), (0, 1, 0, -7), (0, 0, 1, -7),
               (0, 0, 0, -2), (0, 0, 1, 7)]
    assert gpu.collect_cpu(records) == ((0, 0, 0, 5), (0, 1, 0, -7), (1, 0, 0, -7))


@pytest.mark.parametrize("bits", [0, 1, 63, 64, 65, 127, 128, 1024, 32769])
def test_limb_arithmetic_preserves_huge_signed_cancellation(bits):
    large = 1 << bits
    records = [(0, 0, 0, large)] * 7 + [(0, 0, 0, -large)] * 6
    records += [(0, 1, -large, -large), (0, 1, -large, large - 3)]
    packed, timings = gpu.pack_contributions(records)
    assert packed.limbs * 64 >= bits + 1 + (len(records) - 1).bit_length() + 2
    assert gpu.unpack_sums(packed, emulate_limbs(packed)) == gpu.collect_cpu(records)
    assert all(value >= 0 for value in timings.values())


def test_random_packing_matches_independent_arbitrary_precision_collection():
    rng = random.Random(90210)
    records = [(rng.randrange(4), rng.randrange(8), rng.randrange(-5, 6),
                rng.getrandbits(rng.randrange(1, 2048)) * rng.choice([-1, 1]))
               for _ in range(1000)]
    packed, _ = gpu.pack_contributions(records)
    assert gpu.unpack_sums(packed, emulate_limbs(packed)) == gpu.collect_cpu(records)


def test_empty_and_zero_groups():
    for records in [[], [(0, 0, 0, 0)], [(0, 0, 0, 11), (0, 0, 0, -11)]]:
        packed, _ = gpu.pack_contributions(records)
        assert gpu.unpack_sums(packed, emulate_limbs(packed)) == ()
    packed, _ = gpu.pack_contributions([])
    with pytest.raises(ValueError, match="wrong size"):
        gpu.unpack_sums(packed, b"wrong")


@pytest.mark.parametrize("record", [(0, 0, 0), (0, 0, 0, 1, 2), (0, 0, 0, 1.0),
                                    (False, 0, 0, 1), (0, 0, "0", 1), "abcd"])
def test_malformed_records_are_rejected(record):
    with pytest.raises(ValueError):
        gpu.collect_cpu([record])
    with pytest.raises(ValueError):
        gpu.pack_contributions([record])


class EmulatedCollector(gpu._Collector):
    metadata = {"backend": "test-only-limb-emulation"}
    initialization_seconds = 0

    def _execute(self, packed):
        return emulate_limbs(packed), {"h2d": 0, "reduce": 0, "d2h": 0}


def test_result_timing_and_explicit_independent_verification(monkeypatch):
    collector = EmulatedCollector()
    records = [(0, 0, 0, 2**128), (0, 0, 0, -2**128 + 3)]
    result = collector.collect(records)
    assert result.output == ((0, 0, 0, 3),)
    assert result.metadata["verified"]
    assert result.metadata["sort_location"] == "host"
    assert result.metadata["h2d_bytes"] > 0
    assert set(result.timings) == {"validation", "packing", "sort", "h2d", "reduce",
                                  "d2h", "reconstruction", "verification", "total"}
    assert all(value >= 0 for value in result.timings.values())
    assert sum(value for key, value in result.timings.items() if key != "total") <= result.timings["total"]
    assert collector.collect([]).output == ()
    with pytest.raises(ValueError, match="boolean"):
        collector.collect(records, verify=1)
    monkeypatch.setattr(gpu, "collect_cpu", lambda _: ())
    with pytest.raises(gpu.CollectionVerificationError):
        collector.collect(records)
    assert collector.collect(records, verify=False).output == ((0, 0, 0, 3),)


@pytest.mark.parametrize("collector_type", [gpu.CudaCollector, gpu.OpenClCollector])
def test_invalid_device_selection(collector_type):
    for value in [-1, True, 0.5]:
        with pytest.raises(ValueError):
            collector_type(value)


@pytest.mark.parametrize("collector_type", [gpu.CudaCollector, gpu.OpenClCollector])
def test_missing_backend_does_not_silently_execute_on_cpu(monkeypatch, collector_type):
    def missing(_):
        raise ImportError("missing test backend")
    monkeypatch.setattr(gpu.importlib, "import_module", missing)
    with pytest.raises(gpu.GPUUnavailableError, match="missing test backend"):
        collector_type()


@pytest.mark.parametrize("collector_type", [gpu.CudaCollector, gpu.OpenClCollector])
def test_actual_gpu_matches_exact_oracle(collector_type):
    try:
        collector = collector_type()
    except gpu.GPUUnavailableError as error:
        pytest.skip(str(error))
    rng = random.Random(173)
    records = [(rng.randrange(3), rng.randrange(5), rng.randrange(-3, 4),
                rng.choice([-1, 1]) * rng.getrandbits(rng.randrange(0, 2048)))
               for _ in range(513)]
    result = collector.collect(records)
    assert result.output == gpu.collect_cpu(records)
    assert result.metadata["backend"] in {"cuda", "opencl"}
    assert result.metadata["verified"]


def fake_backend(monkeypatch, backend, *, device_present=True, compile_error=False,
                 little_endian=True, cuda_bytes_name=True):
    """Driver-contract doubles, never presented as GPU hardware evidence."""
    np = pytest.importorskip("numpy")

    def reduce_buffers(coefficients, offsets, groups, limbs):
        rows = np.frombuffer(coefficients, dtype="<u8").tobytes()
        offsets = np.frombuffer(offsets, dtype="<u8")
        width = int(limbs) * 8
        raw = []
        for group in range(int(groups)):
            total = sum(int.from_bytes(rows[row * width:(row + 1) * width],
                                       "little", signed=True)
                        for row in range(int(offsets[group]), int(offsets[group + 1])))
            raw.append(total.to_bytes(width, "little", signed=True))
        return b"".join(raw)

    if backend == "cuda":
        class Device:
            def __init__(self, index):
                self.index = index
            def __enter__(self):
                return self
            def __exit__(self, *args):
                return None
            def synchronize(self):
                return None

        class Kernel:
            def __init__(self, source, name):
                assert source == gpu.CUDA_SOURCE and name == "sum_groups"
            def compile(self):
                if compile_error:
                    raise RuntimeError("deliberate CUDA compilation failure")
            def __call__(self, grid, block, args):
                assert block == (128,) and grid[0] >= 1
                coefficients, offsets, output, groups, limbs = args
                assert coefficients.dtype == offsets.dtype == output.dtype == np.dtype("uint64")
                raw = reduce_buffers(coefficients.tobytes(), offsets.tobytes(), groups, limbs)
                output[:] = np.frombuffer(raw, dtype="<u8")

        runtime = SimpleNamespace(
            getDeviceCount=lambda: int(device_present),
            getDeviceProperties=lambda _: {"name": b"contract-device" if cuda_bytes_name else "contract-device"},
            driverGetVersion=lambda: 123, runtimeGetVersion=lambda: 456)
        package = SimpleNamespace(cuda=SimpleNamespace(runtime=runtime, Device=Device),
                                  RawKernel=Kernel, asarray=np.asarray, empty=np.empty,
                                  uint64=np.uint64, asnumpy=np.asarray, __version__="test")
        package_name = "cupy"
    else:
        class Buffer:
            def __init__(self, context, flags, size):
                self.flags, self.data = flags, bytearray(size)

        class Queue:
            def __init__(self, context):
                pass
            def finish(self):
                pass

        class Program:
            def __init__(self, context, source):
                assert source == gpu.OPENCL_SOURCE
            def build(self):
                if compile_error:
                    raise RuntimeError("deliberate OpenCL compilation failure")
                return self

        class Kernel:
            def __init__(self, program, name):
                assert name == "sum_groups"
            def __call__(self, queue, shape, local, coefficients, offsets, output, groups, limbs):
                # The real kernel reads/writes its accumulator. A WRITE_ONLY
                # allocation violates OpenCL's access contract even if a driver
                # happens to tolerate that undefined operation on one device.
                if output.flags != 3:
                    raise RuntimeError("kernel cannot read a WRITE_ONLY accumulator")
                assert shape == (int(groups),) and local is None
                output.data[:] = reduce_buffers(coefficients.data, offsets.data, groups, limbs)

        def copy(queue, destination, source, *, is_blocking):
            assert is_blocking
            if isinstance(destination, Buffer):
                destination.data[:] = source
            else:
                destination[:] = np.frombuffer(source.data, dtype="<u8")

        device = SimpleNamespace(endian_little=little_endian, name=" contract-device ",
                                 vendor=" test ", driver_version="test", version="test",
                                 global_mem_size=100000, max_compute_units=2)
        platform = SimpleNamespace(get_devices=lambda **_: [device] if device_present else [])
        package = SimpleNamespace(get_platforms=lambda: [platform], device_type=SimpleNamespace(GPU=4),
                                  Context=lambda _: object(), CommandQueue=Queue, Program=Program,
                                  Kernel=Kernel, mem_flags=SimpleNamespace(READ_ONLY=1, WRITE_ONLY=2,
                                                                          READ_WRITE=3),
                                  Buffer=Buffer, enqueue_copy=copy, VERSION_TEXT="test")
        package_name = "pyopencl"
    def load(name):
        return package if name == package_name else np
    monkeypatch.setattr(gpu.importlib, "import_module", load)
    return package


@pytest.mark.parametrize("backend", ["cuda", "opencl"])
def test_backend_copy_launch_and_buffer_contracts(monkeypatch, backend):
    fake_backend(monkeypatch, backend)
    collector_type = gpu.CudaCollector if backend == "cuda" else gpu.OpenClCollector
    collector = collector_type()
    rows = [(0, 0, 0, 2**129), (0, 0, 0, 5 - 2**129), (1, 1, -2, -7)]
    result = collector.collect(rows)
    assert result.output == gpu.collect_cpu(rows)
    assert result.metadata["name"] == "contract-device"
    assert result.metadata["kernel_source_sha256"]
    assert result.metadata["backend"] == backend


@pytest.mark.parametrize("backend", ["cuda", "opencl"])
def test_missing_device_and_kernel_compile_failures_are_distinct(monkeypatch, backend):
    collector_type = gpu.CudaCollector if backend == "cuda" else gpu.OpenClCollector
    fake_backend(monkeypatch, backend, device_present=False)
    with pytest.raises(gpu.GPUUnavailableError):
        collector_type()
    fake_backend(monkeypatch, backend, compile_error=True)
    with pytest.raises(gpu.GPUInitializationError, match="compilation failure"):
        collector_type()


def test_cuda_string_device_name_and_opencl_endianness_guard(monkeypatch):
    fake_backend(monkeypatch, "cuda", cuda_bytes_name=False)
    assert gpu.CudaCollector().metadata["name"] == "contract-device"
    fake_backend(monkeypatch, "opencl", little_endian=False)
    with pytest.raises(gpu.GPUUnavailableError, match="little-endian"):
        gpu.OpenClCollector()
    with pytest.raises(NotImplementedError):
        gpu._Collector()._execute(None)
