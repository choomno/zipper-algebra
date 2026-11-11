"""Adversarial packed-rewrite admission and exact contextual output contracts."""

from copy import deepcopy
from types import SimpleNamespace

import pytest

from zipper_algebra.zipping_algorithm.catalog import CatalogError
from zipper_algebra.zipping_algorithm.interpreter import Branch, Limits, ReplayError
from zipper_algebra.zipping_algorithm import gpu_rewrite_host as host
from zipper_algebra.zipping_algorithm import gpu_rewrite_program as program
from zipper_algebra.zipping_algorithm import gpu_rewrite as gpu


def term(word, coefficient=1, power=0):
    return {"word": word, "coefficient": [{"integer": coefficient, "q_power": power}]}


def sigma(index=1, sign=1):
    return {"type": "sigma", "index": index, "sign": sign}


def case(rhs=None, *, coefficient=((0, 7), (2, -3)), words=None, controls=None):
    """An explicit contextual primitive, separate from geometric admission."""
    return {
        "root": 3, "id": 9, "batch": 4, "n": 3,
        "rule": {"id": "test.contextual", "execution": {"algebra_verified": True},
                 "condition": True, "lhs": [term(sigma())],
                 "rhs": [term(sigma(2))] if rhs is None else rhs},
        "parameters": {"n": 3, "i": 1, "j": 1, "k": 0, "l": 1, "eps": 1},
        "words": {} if words is None else words,
        "offset": 1,
        "predecessor": Branch((-1, 1, -2), coefficient).to_json(),
        "successor_controls": controls,
        "limits": vars(Limits(expansion_budget=10000)),
    }


def admit(value):
    return host.admit(value, program.compile_rule(value["rule"]))


@pytest.mark.parametrize("word", [
    sigma(0), sigma(3), sigma(1, 0), sigma(1, 2),
    {"type": "fat", "length": -1, "start": 1},
    {"type": "fat", "length": 2, "start": 2},
    {"type": "block", "length": 1, "start": 1, "count": -1, "order": "ascending"},
    {"type": "half_twist", "size": 0, "start": 1},
    {"type": "half_twist", "size": 3, "start": 2},
    {"type": "slider", "size": 0, "start": 1},
    {"type": "slider", "size": 2, "start": 0},
    {"type": "slider", "size": 3, "start": 2},
])
def test_invalid_primitive_geometry_is_rejected_before_device_allocation(word):
    with pytest.raises((CatalogError, ReplayError)):
        admit(case([term(word)]))


@pytest.mark.parametrize("words", [{}, {"U": [0]}, {"U": [3]},
                                   {"U": [True]}, {"U": [("slider", 2, 1, False, False, False)]}])
def test_named_refs_require_bound_pure_braid_words(words):
    with pytest.raises((CatalogError, ReplayError)):
        admit(case([term({"type": "word_ref", "name": "U"})], words=words))


@pytest.mark.parametrize("old_power,shift", [
    (program.I64_MAX, 1), (program.I64_MIN, -1),
])
def test_q_shift_rejects_intermediate_signed_64_bit_overflow(old_power, shift):
    with pytest.raises((CatalogError, ReplayError)):
        admit(case([term(sigma(2), power=shift)], coefficient=((old_power, 1),)))


def test_zero_monomial_removal_precedes_word_order_and_control_assignment():
    value = case([term(sigma(2)), term(sigma(1)), term(sigma(2), -1),
                  term(sigma(2), power=1)], controls=[{"label": "B"}, {"label": "A"}])
    outputs = host.expand_cpu([value])
    assert outputs == ({"root": 3, "id": 9, "branches": [
        Branch((-1, 1, -2), ((0, 7), (2, -3)), '{"label":"B"}').to_json(),
        Branch((-1, 2, -2), ((1, 7), (3, -3)), '{"label":"A"}').to_json(),
    ]},)


def test_contextual_cpu_reference_preserves_inputs_and_independent_roots():
    first = case()
    second = deepcopy(first)
    second["root"] = 4
    before = deepcopy([first, second])
    result = host.expand_cpu([first, second])
    assert [item["root"] for item in result] == [3, 4]
    assert result[0]["branches"] == result[1]["branches"]
    assert [first, second] == before


def test_direct_cpu_reference_handles_coefficients_beyond_json_decimal_limit():
    huge = 1 << 20000
    value = case([term(sigma(2), -1, power=-3)], coefficient=((-2, huge), (4, -huge + 1)))
    result = host.expand_cpu([value])
    assert Branch.from_json(result[0]["branches"][0]).coefficient == ((-5, -huge), (1, huge - 1))


def test_controls_are_counted_after_actual_rhs_word_collection():
    value = case([term(sigma(2)), term(sigma(2), power=1)], controls=[{}, {}])
    with pytest.raises((CatalogError, ReplayError)):
        host.expand_cpu([value])


def test_failed_compilation_does_not_poison_reusable_template_identity():
    expander = object.__new__(gpu.OpenClRewriteExpander)
    expander._templates = []
    expander._template_keys = {}
    expander._templates_dirty = False
    invalid = case([term({"type": "unknown"})])["rule"]
    for _ in range(2):
        with pytest.raises(CatalogError, match="Unsupported"):
            expander._template(invalid)
        assert expander._templates == []
        assert expander._template_keys == {}
    valid = case()["rule"]
    index, compiled = expander._template(valid)
    assert index == 0 and compiled == program.compile_rule(valid)
    with pytest.raises(CatalogError, match="Unsupported"):
        expander._template(invalid)
    assert expander._template(valid) == (index, compiled)
    assert len(expander._templates) == len(expander._template_keys) == 1


class FakeBytes(bytearray):
    def tobytes(self):
        return bytes(self)


class FakeValues(list):
    def tolist(self):
        return list(self)


def fake_opencl(monkeypatch, *, present=True, little_endian=True,
                compilation_error=False, status_failure=False, corrupt=False):
    """Transport/launch contract double, not evidence of GPU execution."""
    import struct

    def ints(payload):
        return list(struct.unpack("<" + "q" * (len(payload) // 8), payload))

    def write_ints(destination, values):
        raw = struct.pack("<" + "q" * len(values), *values)
        destination.data[:len(raw)] = raw

    class Buffer:
        def __init__(self, context, flags, size=None, hostbuf=None):
            self.flags = flags
            if hostbuf is None:
                assert size >= 8
                self.data = bytearray(size)
            else:
                assert flags == 5
                self.data = bytearray(hostbuf)

    class Queue:
        def __init__(self, context):
            self.finished = 0
        def finish(self):
            self.finished += 1

    class Program:
        def __init__(self, context, source):
            assert source == gpu.OPENCL_SOURCE
        def build(self):
            if compilation_error:
                raise RuntimeError("deliberate driver compilation failure")
            return self

    class Kernel:
        def __init__(self, compiled, name):
            self.name = name

        def __call__(self, queue, shape, local, *arguments):
            assert local is None
            if self.name == "evaluate_registers":
                programs, expressions, cases, parameters, values, status, count = arguments
                assert shape == (count,) and values.flags == status.flags == 3
                programs, expressions, cases, parameters = (
                    ints(b.data) for b in (programs, expressions, cases, parameters))
                output = ints(values.data)
                for c in range(count):
                    begin, length = programs[cases[c * 10] * 4:][:2]
                    base = cases[c * 10 + 1]
                    for j in range(length):
                        op, left, right, immediate = expressions[(begin + j) * 4:][:4]
                        if op == program.CONST:
                            result = immediate
                        elif op == program.PARAM:
                            result = parameters[c * 6 + immediate]
                        else:
                            a, b = output[base + left], output[base + right]
                            result = {program.ADD: lambda: a + b,
                                      program.SUB: lambda: a - b,
                                      program.MUL: lambda: a * b}[op]()
                        output[base + j] = result
                write_ints(values, output)
                if status_failure:
                    status.data[:4] = struct.pack("<I", 1)
            elif self.name == "expand_words":
                cases, values, terms, primitives, words, refs, named_words, jobs, output, status, count = arguments
                assert shape == (count,) and output.flags == status.flags == 3
                cases, values, terms, primitives, words, refs, named_words, jobs = (
                    ints(b.data) for b in (cases, values, terms, primitives, words, refs, named_words, jobs))
                result = ints(output.data)
                for j in range(count):
                    c, t, base, capacity = jobs[j * 4:][:4]
                    _, vo, old, oldlen, prefix, end = cases[c * 10:][:6]
                    emitted = words[old * 6:(old + prefix) * 6]
                    begin, length = terms[t * 4:][:2]
                    for p in range(begin, begin + length):
                        kind, a, b, _, _, flags = primitives[p * 6:][:6]
                        if kind == program.SIGMA:
                            assert flags == 0
                            emitted += [0, values[vo + a] * values[vo + b], 0, 0, 0, 0]
                        elif kind == program.SLIDER:
                            assert flags == 0
                            emitted += [1, values[vo + a], values[vo + b], 0, 0, 0]
                        else:
                            assert kind == program.WORD_REF and flags == 0
                            reference_index = cases[c * 10 + 6] + a
                            start, count = refs[reference_index * 2:][:2]
                            for generator in named_words[start:start + count]:
                                emitted += [0, generator, 0, 0, 0, 0]
                    emitted += words[(old + end) * 6:(old + oldlen) * 6]
                    assert len(emitted) == capacity * 6
                    result[base * 6:(base + capacity) * 6] = emitted
                write_ints(output, result)
            else:
                assert self.name == "propagate_coefficients"
                cases, values, factors, old_q, source, jobs, q_output, output, status, limbs, count = arguments
                assert shape == (count,) and output.flags == status.flags == 3
                cases, values, factors, old_q, jobs = (
                    ints(b.data) for b in (cases, values, factors, old_q, jobs))
                powers = ints(q_output.data)
                width = limbs * 8
                for j in range(count):
                    c, f, old, target = jobs[j * 4:][:4]
                    power_reg, sign = factors[f * 2:][:2]
                    power = values[cases[c * 10 + 1] + power_reg]
                    powers[target] = old_q[old] + power if sign else power
                    value = int.from_bytes(source.data[old * width:(old + 1) * width],
                                           "little", signed=True) * sign
                    if corrupt:
                        value += 1
                    output.data[target * width:(target + 1) * width] = value.to_bytes(
                        width, "little", signed=True)
                write_ints(q_output, powers)

    def copy(queue, destination, source, *, is_blocking):
        assert is_blocking
        if isinstance(destination, Buffer):
            destination.data[:len(source)] = source
        else:
            destination[:] = source.data[:len(destination)]

    device = SimpleNamespace(endian_little=little_endian, name=" contract GPU ",
                             vendor=" test ", driver_version="test", version="test")
    platform = SimpleNamespace(get_devices=lambda **_: [device] if present else [])
    cl = SimpleNamespace(get_platforms=lambda: [platform], device_type=SimpleNamespace(GPU=4),
                         Context=lambda _: object(), CommandQueue=Queue, Program=Program,
                         Kernel=Kernel, mem_flags=SimpleNamespace(READ_ONLY=1, READ_WRITE=3,
                                                                 COPY_HOST_PTR=4),
                         Buffer=Buffer, enqueue_copy=copy, VERSION_TEXT="test")
    def frombuffer(payload, dtype):
        assert dtype == "<u4"
        return FakeValues(struct.unpack("<" + "I" * (len(payload) // 4), payload))
    np = SimpleNamespace(int64=int, uint8="uint8", empty=lambda size, **_: FakeBytes(size),
                         frombuffer=frombuffer)
    monkeypatch.setattr(gpu.importlib, "import_module", lambda name: cl if name == "pyopencl" else np)
    return cl


def test_mock_driver_transport_and_template_cache_are_exact(monkeypatch):
    fake_opencl(monkeypatch)
    expander = gpu.OpenClRewriteExpander()
    first = case([term(sigma(2), -1, -3)])
    expected = host.expand_cpu([first])
    result = expander.expand([first])
    assert result.output == expected
    assert result.metadata["new_templates"] == 1 and result.metadata["verified"]
    second = deepcopy(first)
    second["root"] += 1
    second["parameters"]["i"] = 2
    assert expander.expand([second]).output == host.expand_cpu([second])
    assert expander.expand([second]).metadata["new_templates"] == 0
    assert expander.expand([]).output == ()
    assert expander.expand([case([])]).output[0]["branches"] == []
    assert expander.expand([case([{"word": {"type": "identity"}, "coefficient": []}])]).output[0]["branches"] == []


def test_mock_device_orders_surviving_monomials_and_rejects_bad_controls(monkeypatch):
    fake_opencl(monkeypatch)
    expander = gpu.OpenClRewriteExpander()
    ordered = case([term(sigma(2)), term(sigma(1)), term(sigma(2), -1),
                    term(sigma(2), power=1)], controls=[{"label": "B"}, {"label": "A"}])
    assert expander.expand([ordered]).output == host.expand_cpu([ordered])
    bad = case([term(sigma(2)), term(sigma(2), power=1)], controls=[{}, {}])
    with pytest.raises(ReplayError, match="control update"):
        expander.expand([bad])


def test_zero_factor_resurrection_retains_its_original_word_position(monkeypatch):
    fake_opencl(monkeypatch)
    value = case([term(sigma(2), 0), term(sigma(1)), term(sigma(2))],
                 coefficient=((5, 9),), controls=[{"label": "A"}, {"label": "B"}])
    result = gpu.OpenClRewriteExpander().expand([value])
    assert result.output == host.expand_cpu([value])
    assert result.output[0]["branches"][0]["word"] == [-1, 2, -2]
    assert result.output[0]["branches"][0]["control"] == {"label": "A"}


@pytest.mark.parametrize("old_power,shift", [(program.I64_MAX, 1), (program.I64_MIN, -1)])
def test_zero_factor_never_evaluates_an_irrelevant_overflowing_product_shift(monkeypatch, old_power, shift):
    fake_opencl(monkeypatch)
    value = case([term(sigma(2), 0, shift), term(sigma(1))],
                 coefficient=((old_power, 7),))
    admit(value)
    assert gpu.OpenClRewriteExpander().expand([value]).output == host.expand_cpu([value])


def test_zero_only_word_still_checks_syntax_but_does_not_limit_surviving_words(monkeypatch):
    fake_opencl(monkeypatch)
    value = case([term({"type": "concat", "items": [sigma(), sigma()]}, 0),
                  term(sigma())])
    value["limits"]["max_word_size"] = 3
    assert admit(value).lengths == (4, 3)
    assert gpu.OpenClRewriteExpander().expand([value]).output == host.expand_cpu([value])


def test_cached_program_uses_current_bindings_and_rhs_value_semantics(monkeypatch):
    fake_opencl(monkeypatch)
    expander = gpu.OpenClRewriteExpander()
    original = case([term(sigma("l"))])
    assert expander.expand([original]).output == host.expand_cpu([original])


    rebound = deepcopy(original)
    rebound["parameters"]["l"] = 2
    result = expander.expand([rebound])
    assert result.output == host.expand_cpu([rebound])
    assert result.metadata["new_templates"] == 0
    changed = deepcopy(original)
    changed["rule"]["rhs"][0]["word"] = sigma(2)
    result = expander.expand([changed])
    assert result.output == host.expand_cpu([changed])
    assert result.metadata["new_templates"] == 1
    assert expander.expand([original]).output == host.expand_cpu([original])


def test_named_word_and_slider_transport_round_trip_without_gpu_dependencies(monkeypatch):
    fake_opencl(monkeypatch)
    value = case([term({"type": "concat", "items": [
        {"type": "word_ref", "name": "U"},
        {"type": "slider", "size": 2, "start": 1},
    ]})], words={"U": [1, -2]})
    old = Branch.from_json(value["predecessor"])
    value["predecessor"] = Branch((("slider", 2, 1, False, False, False), *old.word),
                                   old.coefficient).to_json()
    value["offset"] += 1
    assert gpu.OpenClRewriteExpander().expand([value]).output == host.expand_cpu([value])


def test_unused_large_parameters_never_overflow_packed_register_slots(monkeypatch):
    fake_opencl(monkeypatch)
    value = case()
    value["parameters"]["eps"] = 1 << 20000
    assert gpu.OpenClRewriteExpander().expand([value]).output == host.expand_cpu([value])


def test_unconsumed_device_monomial_records_are_rejected(monkeypatch):
    fake_opencl(monkeypatch)
    expander = gpu.OpenClRewriteExpander()
    packed = expander._pack([case()])
    raw, _, _ = expander._execute(packed)
    packed.coefficient_records += 1
    raw["coefficient_output"] += bytes(packed.limbs * 8)
    raw["q_output"] += bytes(8)
    with pytest.raises(ReplayError, match="consume every emitted monomial"):
        expander._reconstruct(packed, raw)


def test_boolean_verification_contract_and_big_endian_host_codec(monkeypatch):
    fake_opencl(monkeypatch)
    expander = gpu.OpenClRewriteExpander()
    with pytest.raises(ValueError, match="boolean"):
        expander.expand([], verify=1)

    class BigEndianArray:
        """Simulate native big-endian array storage, not a GPU device."""
        def __init__(self, code, values=()):
            assert code == "q"
            self.values = list(values)
        def byteswap(self):
            self.values = [int.from_bytes(v.to_bytes(8, "big", signed=True), "little", signed=True)
                           for v in self.values]
        def tobytes(self):
            return b"".join(v.to_bytes(8, "big", signed=True) for v in self.values)
        def frombytes(self, payload):
            self.values.extend(int.from_bytes(payload[i:i + 8], "big", signed=True)
                               for i in range(0, len(payload), 8))
        def __iter__(self):
            return iter(self.values)

    monkeypatch.setattr(gpu, "array", BigEndianArray)
    monkeypatch.setattr(gpu, "sys", SimpleNamespace(byteorder="big"))
    values = [-12345678, 0, 2**60 + 1]
    raw = gpu._i64(values)
    assert raw == b"".join(v.to_bytes(8, "little", signed=True) for v in values)
    assert gpu._decode_i64(raw) == values


def test_checked_device_status_does_not_return_a_partial_answer(monkeypatch):
    fake_opencl(monkeypatch, status_failure=True)
    with pytest.raises(ReplayError, match="checked arithmetic/allocation failure"):
        gpu.OpenClRewriteExpander().expand([case()])


def test_independent_verification_detects_corrupt_device_coefficients(monkeypatch):
    fake_opencl(monkeypatch, corrupt=True)
    expander = gpu.OpenClRewriteExpander()
    with pytest.raises(ReplayError, match="independent contextual"):
        expander.expand([case()])
    unverified = expander.expand([case()], verify=False)
    assert not unverified.metadata["verified"]
    assert unverified.output != host.expand_cpu([case()])


def test_gpu_backend_presence_endianness_and_compilation_contracts(monkeypatch):
    for index in [-1, True, 0.5]:
        with pytest.raises(ValueError):
            gpu.OpenClRewriteExpander(index)
    fake_opencl(monkeypatch, present=False)
    with pytest.raises(gpu.GPUUnavailableError):
        gpu.OpenClRewriteExpander()
    fake_opencl(monkeypatch, little_endian=False)
    with pytest.raises(gpu.GPUUnavailableError, match="little-endian"):
        gpu.OpenClRewriteExpander()
    fake_opencl(monkeypatch, compilation_error=True)
    with pytest.raises(gpu.GPUInitializationError, match="compilation failure"):
        gpu.OpenClRewriteExpander()
    monkeypatch.setattr(gpu.importlib, "import_module", lambda _: (_ for _ in ()).throw(ImportError("missing driver")))
    with pytest.raises(gpu.GPUUnavailableError, match="missing driver"):
        gpu.OpenClRewriteExpander()
