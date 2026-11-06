"""Reusable primitive programs for exact packed GPU catalog expansion.

Compilation preserves symbolic integer expressions and repeat descriptors.
It does not expand successor braid words on the CPU. Integer coordinates use
checked signed 64-bit registers; unsupported or overflowing programs fail
closed before dispatch. Coefficient integers in the bundled catalog are 1
and -1; arbitrary-width incoming coefficient propagation is handled by the
GPU executor.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import AbstractSet, Any, Mapping

from .catalog import CatalogError

PARAMETER_NAMES = ("n", "i", "j", "k", "l", "eps")
I64_MIN = -(1 << 63)
I64_MAX = (1 << 63) - 1

CONST, PARAM, ADD, SUB, MUL = range(5)
SIGMA, BLOCK, HALF_TWIST, SLIDER, WORD_REF = range(1, 6)
REVERSE = 1
NEGATE = 2
REFLECT = 4
SLIDER_MIRRORED = 8
SLIDER_REVERSED = 16
SLIDER_BARRED = 32
BLOCK_DESCENDING = 64

Expression = tuple[int, int, int, int]
Primitive = tuple[int, int, int, int, int, int]


@dataclass(frozen=True)
class CompiledTerm:
    primitives: tuple[Primitive, ...]
    coefficients: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class CompiledRule:
    expressions: tuple[Expression, ...]
    terms: tuple[CompiledTerm, ...]
    refs: tuple[str, ...]
    reflection_registers: tuple[int, ...] = ()


def checked_i64(value: Any) -> int:
    if type(value) is not int or not I64_MIN <= value <= I64_MAX:
        raise CatalogError("GPU coordinate/exponent is outside signed 64-bit range")
    return value


def _fields(node: Any, required: set[str], optional: AbstractSet[str] = frozenset()) -> None:
    if not isinstance(node, dict) or not required <= node.keys() or node.keys() - required - optional:
        raise CatalogError("Unexpected/missing packed-program AST fields")


def _contains_slider(node: Any) -> bool:
    if not isinstance(node, dict):
        return False
    return (node.get("type") == "slider"
            or ("word" in node and _contains_slider(node["word"]))
            or (isinstance(node.get("items"), list)
                and any(_contains_slider(child) for child in node["items"])))


class _Compiler:
    def __init__(self) -> None:
        self.expressions: list[Expression] = []
        self.registers: dict[Expression, int] = {}
        self.refs: list[str] = []
        self.reflections: list[int] = []
        self.zero = self.expression(0)
        self.one = self.expression(1)

    def append(self, instruction: Expression) -> int:
        if instruction not in self.registers:
            self.registers[instruction] = len(self.expressions)
            self.expressions.append(instruction)
        return self.registers[instruction]

    def expression(self, expr: Any) -> int:
        if type(expr) is int:
            return self.append((CONST, 0, 0, checked_i64(expr)))
        if isinstance(expr, str) and expr in PARAMETER_NAMES:
            return self.append((PARAM, 0, 0, PARAMETER_NAMES.index(expr)))
        _fields(expr, {"op", "args"})
        op, args = expr["op"], expr["args"]
        if (op not in ("add", "sub", "mul") or not isinstance(args, list)
                or len(args) < 2 or (op == "sub" and len(args) != 2)):
            raise CatalogError("Invalid packed integer operation")
        indices = [self.expression(arg) for arg in args]
        opcode = {"add": ADD, "sub": SUB, "mul": MUL}[op]
        result = indices[0]
        for index in indices[1:]:
            result = self.append((opcode, result, index, 0))
        return result

    def word(self, node: Any, flags: int = 0) -> list[Primitive]:
        if not isinstance(node, dict):
            raise CatalogError("Packed word must be an AST object")
        kind = node.get("type")
        z = self.zero
        if kind == "identity":
            _fields(node, {"type"})
            return []
        if kind == "concat":
            _fields(node, {"type", "items"})
            if not isinstance(node["items"], list):
                raise CatalogError("Packed concat.items must be a list")
            children = node["items"][::-1] if flags & REVERSE else node["items"]
            return [primitive for child in children for primitive in self.word(child, flags)]
        if kind in ("inverse", "bar", "reverse", "reflect"):
            _fields(node, {"type", "word"}, {"n"} if kind == "reflect" else frozenset())
            if kind == "inverse":
                if _contains_slider(node["word"]):
                    raise CatalogError("Packed inverse cannot contain a slider")
                flags ^= REVERSE | NEGATE
            elif kind == "bar":
                flags ^= NEGATE | SLIDER_BARRED
            elif kind == "reverse":
                flags ^= REVERSE | SLIDER_REVERSED
            else:
                self.reflections.append(self.expression(node.get("n", "n")))
                flags ^= REFLECT
            return self.word(node["word"], flags)
        if kind == "sigma":
            _fields(node, {"type", "index"}, {"sign"})
            return [(SIGMA, self.expression(node["index"]),
                     self.expression(node.get("sign", 1)), z, z, flags)]
        if kind in ("fat", "block"):
            required = {"type", "length", "start"}
            if kind == "block":
                required |= {"count", "order"}
            _fields(node, required, {"sign"})
            if kind == "block" and node["order"] not in ("ascending", "descending"):
                raise CatalogError("Invalid packed block order")
            if kind == "block" and node["order"] == "descending":
                flags |= BLOCK_DESCENDING
            return [(BLOCK, self.expression(node["length"]),
                     self.expression(node["start"]),
                     self.expression(node["count"]) if kind == "block" else self.one,
                     self.expression(node.get("sign", 1)), flags)]
        if kind == "half_twist":
            _fields(node, {"type", "size", "start"})
            return [(HALF_TWIST, self.expression(node["size"]),
                     self.expression(node["start"]), z, z, flags)]
        if kind == "slider":
            _fields(node, {"type", "size", "start"}, {"orientation"})
            orientation = node.get("orientation", "normal")
            if orientation not in ("normal", "mirrored"):
                raise CatalogError("Invalid packed slider orientation")
            if orientation == "mirrored":
                flags |= SLIDER_MIRRORED
            return [(SLIDER, self.expression(node["size"]),
                     self.expression(node["start"]), z, z, flags)]
        if kind == "word_ref":
            _fields(node, {"type", "name"})
            name = node["name"]
            if not isinstance(name, str):
                raise CatalogError("Packed word_ref name must be a string")
            if name not in self.refs:
                self.refs.append(name)
            return [(WORD_REF, self.refs.index(name), z, z, z, flags)]
        raise CatalogError("Unsupported packed word node: " + repr(kind))


def compile_rule(rule: Mapping[str, Any]) -> CompiledRule:
    """Compile a rule's RHS into symbolic emit programs, never emitted words."""
    rhs = rule.get("rhs")
    if not isinstance(rhs, list):
        raise CatalogError("Packed rule RHS must be a polynomial list")
    compiler = _Compiler()
    terms = []
    for term in rhs:
        _fields(term, {"coefficient", "word"})
        coefficients = term["coefficient"]
        if not isinstance(coefficients, list):
            raise CatalogError("Packed coefficients must be a Laurent list")
        primitive = tuple(compiler.word(term["word"]))
        compiled_coefficients = []
        for coefficient in coefficients:
            _fields(coefficient, {"integer", "q_power"})
            value = coefficient["integer"]
            if type(value) is not int or value not in (-1, 0, 1):
                raise CatalogError("Packed bundled-rule factors must be zero, one, or minus one")
            power = compiler.expression(coefficient["q_power"])
            # A zero entry can reserve an insertion position later resurrected
            # by another term, so preserve it until monomial collection.
            compiled_coefficients.append((power, value))
        terms.append(CompiledTerm(primitive, tuple(compiled_coefficients)))
    return CompiledRule(tuple(compiler.expressions), tuple(terms),
                        tuple(compiler.refs), tuple(compiler.reflections))


def evaluate(program: CompiledRule, env: Mapping[str, int]) -> tuple[int, ...]:
    """Evaluate allocation/admission coordinates exactly, rejecting overflow."""
    values: list[int] = []
    for opcode, left, right, immediate in program.expressions:
        if opcode == CONST:
            result = immediate
        elif opcode == PARAM:
            name = PARAMETER_NAMES[immediate]
            if name not in env:
                raise CatalogError("Missing packed integer parameter: " + name)
            result = env[name]
        else:
            if opcode not in (ADD, SUB, MUL):
                raise CatalogError("Unknown packed expression opcode")
            a, b = values[left], values[right]
            if opcode == ADD:
                result = a + b
            elif opcode == SUB:
                result = a - b
            else:
                result = a * b
        values.append(checked_i64(result))
    for index in program.reflection_registers:
        if values[index] != env.get("n"):
            raise CatalogError("Packed reflect.n must equal ambient n")
    return tuple(values)
