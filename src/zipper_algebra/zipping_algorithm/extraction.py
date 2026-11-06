"""Exact recognition of braid words carried by an explicit nested ribbon tree.

The introduction order ``labels`` specifies a comb: at level r the first r-1
labels form the old ribbon and label r is the new singleton. This finite
combinatorial witness is required. Connected-component incidence does not
supply it. Words use algebraic order; the parser reads bottom to top.

Whole crossings of the two children collapse to generators of B_2. Naturality
moves every child braid through such a crossing without changing it, so the
signed crossing count at each node determines the braid in this groupoid.
The canonical reconstruction executes inner nodes before outer nodes. This
module normalizes an already certified nested ribbon; it does not prove that
an arbitrary completed zipper phase has the required hierarchy.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from typing import Any, Iterable


class RibbonError(ValueError):
    """The supplied word fails the explicit nested ribbon witness."""


class RibbonLimit(RibbonError):
    """Exact ribbon gathering exceeded its explicit action or operation limit."""


@dataclass(frozen=True)
class RibbonGathering:
    max_image_letters: int
    max_operations: int
    action_sha256: str


@dataclass(frozen=True)
class RibbonBlock:
    node_size: int
    sign: int
    chronological_word: tuple[int, ...]


@dataclass(frozen=True)
class NestedRibbon:
    labels: tuple[int, ...]
    initial_ports: tuple[int, ...]
    word: tuple[int, ...]
    blocks: tuple[RibbonBlock, ...]
    crossing_counts: tuple[int, ...]
    final_ports: tuple[int, ...]
    reconstructed_word: tuple[int, ...]
    gathering: RibbonGathering | None = None

    @property
    def twist_powers(self) -> tuple[int, ...]:
        """Base sigma power and higher full-lasso powers, with parity separate."""
        return self.crossing_counts[:1] + tuple(h // 2 for h in self.crossing_counts[1:])

    @property
    def boundary_parities(self) -> tuple[int, ...]:
        return tuple(h % 2 for h in self.crossing_counts[1:])


def _ports(labels: Iterable[int], initial_ports: Iterable[int] | None) -> tuple[tuple[int, ...], tuple[int, ...]]:
    tree = tuple(labels)
    ports = tree if initial_ports is None else tuple(initial_ports)
    if (not tree or any(type(x) is not int for x in tree)
            or len(set(tree)) != len(tree) or len(ports) != len(tree)
            or any(type(x) is not int for x in ports) or set(ports) != set(tree)):
        raise RibbonError("distinct integer tree labels and the same boundary labels required")
    for size in range(2, len(tree) + 1):
        positions = sorted(ports.index(x) for x in tree[:size])
        if positions[-1] - positions[0] + 1 != size:
            raise RibbonError("every marked subtree must occupy a contiguous boundary interval")
    return tree, ports


def _word(word: Iterable[int], n: int) -> tuple[int, ...]:
    braid = tuple(word)
    if any(type(x) is not int or not 1 <= abs(x) < n for x in braid):
        raise RibbonError("signed braid generators inside the component required")
    return braid


def _children(ports: list[int], tree: tuple[int, ...], size: int) -> tuple[int, int, bool]:
    old = set(tree[:size - 1])
    positions = sorted(ports.index(x) for x in old)
    if positions[-1] - positions[0] + 1 != size - 1:
        raise RibbonError("the marked child ribbon is not contiguous")
    singleton = ports.index(tree[size - 1])
    if singleton == positions[-1] + 1:
        return positions[0], size - 1, True
    if singleton == positions[0] - 1:
        return singleton, size - 1, False
    raise RibbonError("the marked children are not adjacent")


def _canonical_block(ports: list[int], tree: tuple[int, ...], size: int, sign: int) -> tuple[int, ...]:
    start, width, old_left = _children(ports, tree, size)
    # Zero-based ports become the one-based indices of Artin generators.
    block = tuple(sign * x for x in (range(start + width, start, -1) if old_left
                                    else range(start + 1, start + width + 1)))
    for letter in block:
        at = abs(letter) - 1
        ports[at], ports[at + 1] = ports[at + 1], ports[at]
    return block


def reconstruct_nested(labels: Iterable[int], counts: Iterable[int], *,
                       initial_ports: Iterable[int] | None = None) -> tuple[int, ...]:
    """Reconstruct the normal word from node counts and finite boundary data."""
    tree, boundary = _ports(labels, initial_ports)
    powers = tuple(counts)
    if len(powers) != len(tree) - 1 or any(type(x) is not int for x in powers):
        raise RibbonError("one integer crossing count per non-leaf node required")
    ports = list(boundary)
    chronological: list[int] = []
    for size, count in enumerate(powers, 2):
        sign = 1 if count > 0 else -1
        for _ in range(abs(count)):
            chronological.extend(_canonical_block(ports, tree, size, sign))
    return tuple(reversed(chronological))


def recognize_nested_strict(word: Iterable[int], *, labels: Iterable[int],
                     initial_ports: Iterable[int] | None = None) -> NestedRibbon:
    """Require every inter-child interaction to be a full uniform ribbon crossing.

    A complete crossing may use any reduced shuffle of the child leaves. The
    parser checks every pair once and checks a single over/under choice; it does
    not silently commute incomplete blocks or infer a hierarchy from writhe.
    """
    tree, boundary = _ports(labels, initial_ports)
    braid = _word(word, len(tree))
    chronological = tuple(reversed(braid))
    order = {label: idx for idx, label in enumerate(tree)}
    ports = list(boundary)
    counts = [0] * (len(tree) - 1)
    blocks: list[RibbonBlock] = []
    offset = 0
    while offset < len(chronological):
        at = abs(chronological[offset]) - 1
        size = max(order[ports[at]], order[ports[at + 1]]) + 1
        _children(ports, tree, size)
        old = set(tree[:size - 1])
        singleton = tree[size - 1]
        sign = 1 if chronological[offset] > 0 else -1
        block = chronological[offset:offset + size - 1]
        if len(block) != size - 1:
            raise RibbonError("incomplete crossing of the marked child ribbon")
        crossed: set[int] = set()
        for letter in block:
            left = abs(letter) - 1
            a, b = ports[left], ports[left + 1]
            other = b if a == singleton else a
            if ((1 if letter > 0 else -1) != sign or singleton not in (a, b)
                    or other not in old or other in crossed):
                raise RibbonError("a whole ribbon crossing must meet every child leaf once with one sign")
            crossed.add(other)
            ports[left], ports[left + 1] = b, a
        counts[size - 2] += sign
        blocks.append(RibbonBlock(size, sign, block))
        offset += size - 1
    canonical = reconstruct_nested(tree, counts, initial_ports=boundary)
    return NestedRibbon(tree, boundary, braid, tuple(blocks), tuple(counts), tuple(ports), canonical)


def _gather_limits(max_image_letters: int, max_operations: int) -> None:
    if (type(max_image_letters) is not int or max_image_letters < 1
            or type(max_operations) is not int or max_operations < 1):
        raise RibbonError("positive integer ribbon gathering limits are required")


def _artin_action(word: tuple[int, ...], n: int, *, max_image_letters: int,
                  max_operations: int, operations: list[int]) -> tuple[tuple[int, ...], ...]:
    if n > max_image_letters:
        raise RibbonLimit("ribbon gathering Artin image limit exceeded")
    images: list[tuple[int, ...]] = [(idx,) for idx in range(1, n + 1)]
    total = n
    for letter in word:
        operations[0] += 1
        if operations[0] > max_operations:
            raise RibbonLimit("ribbon gathering operation limit exceeded")
        at = abs(letter) - 1
        a, b = images[at:at + 2]
        stack: list[int] = []
        # Consume inverse images lazily, without allocating an unbounded
        # intermediate substituted free word.
        pieces = (iter(a), iter(b), (-x for x in reversed(a))) if letter > 0 else (
            (-x for x in reversed(b)), iter(a), iter(b))
        for piece in pieces:
            for x in piece:
                operations[0] += 1
                if operations[0] > max_operations:
                    raise RibbonLimit("ribbon gathering operation limit exceeded")
                if stack and stack[-1] == -x:
                    stack.pop()
                else:
                    stack.append(x)
                if len(stack) > max_image_letters:
                    raise RibbonLimit("ribbon gathering Artin image limit exceeded")
        changed = tuple(stack)
        total += len(changed) - (len(b) if letter > 0 else len(a))
        if total > max_image_letters:
            raise RibbonLimit("ribbon gathering Artin image limit exceeded")
        if letter > 0:
            images[at:at + 2] = [changed, a]
        else:
            images[at:at + 2] = [b, changed]
    return tuple(images)


def _gathering(word: tuple[int, ...], canonical: tuple[int, ...], n: int,
               max_image_letters: int, max_operations: int) -> RibbonGathering:
    _gather_limits(max_image_letters, max_operations)
    operations = [0]
    original_action = _artin_action(word, n, max_image_letters=max_image_letters,
                                    max_operations=max_operations, operations=operations)
    canonical_action = _artin_action(canonical, n, max_image_letters=max_image_letters,
                                     max_operations=max_operations, operations=operations)
    if original_action != canonical_action:
        raise RibbonError("pair-count candidate is not equal to the original braid")
    digest = hashlib.sha256(json.dumps(original_action, separators=(",", ":")).encode()).hexdigest()
    return RibbonGathering(max_image_letters, max_operations, digest)


def recognize_nested(word: Iterable[int], *, labels: Iterable[int],
                     initial_ports: Iterable[int] | None = None,
                     max_image_letters: int = 100000,
                     max_operations: int = 10000000) -> NestedRibbon:
    """Recognize a marked ribbon element, allowing exact braid gathering.

    Literal whole blocks have a direct structural witness. Otherwise the
    signed crossings of actual leaf pairs give only a candidate node vector.
    Equality to its reconstructed braid is checked by the faithful Artin
    action with explicit allocation and operation limits. Pair counts alone
    never establish membership, and a limit hit leaves the word unresolved.
    """
    _gather_limits(max_image_letters, max_operations)
    tree, boundary = _ports(labels, initial_ports)
    braid = _word(word, len(tree))
    try:
        return recognize_nested_strict(braid, labels=tree, initial_ports=boundary)
    except RibbonError:
        pass
    order = {label: idx for idx, label in enumerate(tree)}
    ports = list(boundary)
    pairs: dict[tuple[int, int], int] = {}
    for letter in reversed(braid):
        at = abs(letter) - 1
        a, b = order[ports[at]], order[ports[at + 1]]
        pair = min(a, b), max(a, b)
        pairs[pair] = pairs.get(pair, 0) + (1 if letter > 0 else -1)
        ports[at], ports[at + 1] = ports[at + 1], ports[at]
    counts = []
    for node in range(1, len(tree)):
        count = pairs.get((0, node), 0)
        if any(pairs.get((leaf, node), 0) != count for leaf in range(1, node)):
            raise RibbonError("leaf-pair crossing counts do not define a marked ribbon candidate")
        counts.append(count)
    canonical = reconstruct_nested(tree, counts, initial_ports=boundary)
    action = _gathering(braid, canonical, len(tree), max_image_letters, max_operations)
    target = recognize_nested_strict(canonical, labels=tree, initial_ports=boundary)
    return NestedRibbon(tree, boundary, braid, target.blocks, tuple(counts), tuple(ports), canonical, action)


def check_nested(witness: NestedRibbon) -> None:
    """Check the saved full-block witness and its canonical reconstruction."""
    if not isinstance(witness, NestedRibbon):
        raise RibbonError("a nested ribbon witness is required")
    tree, boundary = _ports(witness.labels, witness.initial_ports)
    braid = _word(witness.word, len(tree))
    if witness.gathering is not None:
        if not isinstance(witness.gathering, RibbonGathering):
            raise RibbonError("an explicit exact braid gathering witness is required")
        action = witness.gathering
        check_nested(NestedRibbon(witness.labels, witness.initial_ports,
                                  witness.reconstructed_word, witness.blocks,
                                  witness.crossing_counts, witness.final_ports,
                                  witness.reconstructed_word))
        expected = _gathering(braid, witness.reconstructed_word, len(tree),
                             action.max_image_letters, action.max_operations)
        if action != expected:
            raise RibbonError("exact gathering action does not match its saved identity")
        return
    if not isinstance(witness.blocks, tuple):
        raise RibbonError("immutable crossing blocks are required")
    ports = list(boundary)
    counts = [0] * (len(tree) - 1)
    chronological: list[int] = []
    for block in witness.blocks:
        if (not isinstance(block, RibbonBlock) or type(block.node_size) is not int
                or not 2 <= block.node_size <= len(tree)
                or type(block.sign) is not int or block.sign not in (-1, 1)
                or not isinstance(block.chronological_word, tuple)
                or len(block.chronological_word) != block.node_size - 1):
            raise RibbonError("invalid marked crossing block")
        _children(ports, tree, block.node_size)
        singleton = tree[block.node_size - 1]
        remaining = set(tree[:block.node_size - 1])
        for letter in _word(block.chronological_word, len(tree)):
            at = abs(letter) - 1
            a, b = ports[at], ports[at + 1]
            other = b if a == singleton else a
            if (singleton not in (a, b) or other not in remaining
                    or (1 if letter > 0 else -1) != block.sign):
                raise RibbonError("invalid local whole-ribbon crossing")
            remaining.remove(other)
            ports[at], ports[at + 1] = b, a
        counts[block.node_size - 2] += block.sign
        chronological.extend(block.chronological_word)
    if (tuple(reversed(chronological)) != braid or tuple(counts) != witness.crossing_counts
            or tuple(ports) != witness.final_ports
            or reconstruct_nested(tree, counts, initial_ports=boundary) != witness.reconstructed_word):
        raise RibbonError("nested ribbon witness does not match its complete braid word")



@dataclass(frozen=True)
class SliderExtraction:
    """A component with its coupon and both independently marked boundary trees."""
    word: tuple[Any, ...]
    slider: tuple[Any, ...]
    upper: NestedRibbon
    lower: NestedRibbon
    reconstructed_word: tuple[Any, ...]


def extract_slider(word: Iterable[Any], *, upper_labels: Iterable[int],
                   lower_labels: Iterable[int], upper_ports: Iterable[int] | None = None,
                   lower_ports: Iterable[int] | None = None) -> SliderExtraction:
    """Normalize the two nested ribbon sides, retaining the complete slider token.

    Input is the isolated component, reindexed to its own rank. The caller
    supplies the two marked hierarchies obtained from its zipping certificate.
    The two label sets independently name their own boundary ports.
    No internal permutation of the slider is assigned.
    """
    diagram = tuple(word)
    sliders = [idx for idx, token in enumerate(diagram) if type(token) is not int]
    if len(sliders) != 1:
        raise RibbonError("an isolated component with exactly one explicit slider is required")
    at = sliders[0]
    slider = diagram[at]
    if (not isinstance(slider, tuple) or len(slider) != 6 or slider[0] != "slider"
            or type(slider[1]) is not int or slider[1] < 1
            or type(slider[2]) is not int or slider[2] != 1
            or any(type(x) is not bool for x in slider[3:])):
        raise RibbonError("a full-rank slider anchored at the first component port is required")
    upper = recognize_nested(diagram[:at], labels=upper_labels, initial_ports=upper_ports)
    lower = recognize_nested(diagram[at + 1:], labels=lower_labels, initial_ports=lower_ports)
    if len(upper.labels) != slider[1] or len(lower.labels) != slider[1]:
        raise RibbonError("both marked trees must have the slider's rank")
    reconstructed = upper.reconstructed_word + (slider,) + lower.reconstructed_word
    return SliderExtraction(diagram, slider, upper, lower, reconstructed)


def check_extraction(witness: SliderExtraction) -> None:
    """Check both sides and the retained coupon without expanding that coupon."""
    if not isinstance(witness, SliderExtraction):
        raise RibbonError("a slider extraction witness is required")
    check_nested(witness.upper)
    check_nested(witness.lower)
    expected = extract_slider(witness.word, upper_labels=witness.upper.labels,
                              lower_labels=witness.lower.labels,
                              upper_ports=witness.upper.initial_ports,
                              lower_ports=witness.lower.initial_ports)
    if witness != expected:
        raise RibbonError("extraction witness does not reconstruct the complete component")


@dataclass(frozen=True)
class LayeredSplit:
    """Deletion and explicit overcrossing reconstruction of one component."""
    n: int
    word: tuple[Any, ...]
    bottom_ports: tuple[int, ...]
    top_ports: tuple[int, ...]
    component_word: tuple[Any, ...]
    complement_word: tuple[int, ...]
    bottom_route: tuple[int, ...]
    top_route: tuple[int, ...]
    reconstructed_word: tuple[Any, ...]


def _inverse_braid(word: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-x for x in reversed(word))


def _over_route(n: int, selected: tuple[int, ...]) -> tuple[int, ...]:
    """Stable partition to the left, with the selected strands always over."""
    ports = list(range(1, n + 1))
    chronological = []
    for destination, label in enumerate(selected):
        current = ports.index(label)
        while current > destination:
            chronological.append(-current)
            ports[current - 1], ports[current] = ports[current], ports[current - 1]
            current -= 1
    return tuple(reversed(chronological))


def split_layered(word: Iterable[Any], *, n: int, bottom_ports: Iterable[int]) -> LayeredSplit:
    """Isolate an overcrossing component and retain canonical boundary routes.

    The selected component has one explicit slider joining all of its strands,
    or is a single strand with no slider. Every mixed crossing must have that
    component over the complement. Internal slider permutations are never used:
    only the boolean component membership of each boundary port is transported.

    The exact reconstruction is top_route * (component tensor complement) *
    bottom_route, with the component placed at the left. Straightening an arc
    that lies over the complete complement proves this braid-diagram equality.
    """
    selected = tuple(bottom_ports)
    if (type(n) is not int or n < 1 or not selected
            or any(type(x) is not int or not 1 <= x <= n for x in selected)
            or selected != tuple(sorted(set(selected)))):
        raise RibbonError("an increasing nonempty subset of the ambient bottom ports is required")
    diagram = tuple(word)
    size = len(selected)
    membership = [idx in selected for idx in range(1, n + 1)]
    component: list[Any] = []
    complement: list[int] = []
    coupons = 0
    for token in reversed(diagram):
        if type(token) is int:
            if not 1 <= abs(token) < n:
                raise RibbonError("ambient braid generator is out of range")
            at = abs(token) - 1
            left, right = membership[at:at + 2]
            sign = 1 if token > 0 else -1
            if left and right:
                component.append(sign * (sum(membership[:at]) + 1))
            elif not left and not right:
                complement.append(sign * (at - sum(membership[:at]) + 1))
            elif not (left if token > 0 else right):
                raise RibbonError("the complement crosses over the selected component")
            membership[at], membership[at + 1] = right, left
        else:
            if (not isinstance(token, tuple) or len(token) != 6 or token[0] != "slider"
                    or type(token[1]) is not int or token[1] != size
                    or type(token[2]) is not int or not 1 <= token[2] <= n - size + 1
                    or any(type(x) is not bool for x in token[3:])):
                raise RibbonError("one explicit slider on the complete selected component is required")
            if not all(membership[token[2] - 1:token[2] + size - 1]):
                raise RibbonError("the slider does not join exactly the selected component")
            coupons += 1
            component.append(("slider", size, 1, *token[3:]))
    if coupons > 1 or (size > 1 and coupons != 1):
        raise RibbonError("the selected component requires exactly one slider, except for a singleton")
    top = tuple(idx for idx, active in enumerate(membership, 1) if active)
    active_word = tuple(reversed(component))
    other_word = tuple(reversed(complement))
    bottom_route = _over_route(n, selected)
    top_route = _inverse_braid(_over_route(n, top))
    shifted_complement = tuple((1 if x > 0 else -1) * (abs(x) + size) for x in other_word)
    rebuilt = top_route + active_word + shifted_complement + bottom_route
    return LayeredSplit(n, diagram, selected, top, active_word, other_word,
                        bottom_route, top_route, rebuilt)


def check_split(witness: LayeredSplit) -> None:
    """Check a saved deletion, all overcrossing guards, and its full reconstruction."""
    if not isinstance(witness, LayeredSplit):
        raise RibbonError("a layered component split witness is required")
    expected = split_layered(witness.word, n=witness.n, bottom_ports=witness.bottom_ports)
    if witness != expected:
        raise RibbonError("layered component split witness does not match its source diagram")


@dataclass(frozen=True)
class BoundaryRoute:
    """An exterior signed permutation factor, kept at its literal boundary.

    The scheduler names a fan saved at active size s by growth_size=s+1.
    This public witness admits more general signed permutation factors; its
    validator does not establish the scheduler's sharper fan-origin bound.
    """
    growth_size: int
    word: tuple[int, ...]


def _boundary_routes(routes: Iterable[BoundaryRoute], n: int) -> tuple[BoundaryRoute, ...]:
    factors = tuple(routes)
    sizes: set[int] = set()
    for factor in factors:
        if (not isinstance(factor, BoundaryRoute) or type(factor.growth_size) is not int
                or not 2 <= factor.growth_size <= n or factor.growth_size in sizes
                or not isinstance(factor.word, tuple)):
            raise RibbonError("distinct growth sizes and explicit immutable boundary routes required")
        sizes.add(factor.growth_size)
        ports = list(range(n))
        pairs: set[tuple[int, int]] = set()
        for letter in reversed(_word(factor.word, n)):
            at = abs(letter) - 1
            pair = (min(ports[at], ports[at + 1]), max(ports[at], ports[at + 1]))
            if pair in pairs:
                raise RibbonError("an inherited boundary route must cross each strand pair at most once")
            pairs.add(pair)
            ports[at], ports[at + 1] = ports[at + 1], ports[at]
    return factors


@dataclass(frozen=True)
class CanonicalSlider:
    """Exact component value with one twist vector and a default normal coupon.

    The finite routes are component braid words, not ambient endpoint subsets
    or recoverable history flags. The source word and inherited factors are
    retained separately for witness checking.
    """
    word: tuple[Any, ...]
    slider: tuple[Any, ...]
    upper: NestedRibbon
    lower_inverse: NestedRibbon
    twist_powers: tuple[int, ...]
    top_route: tuple[int, ...]
    bottom_route: tuple[int, ...]
    normalized_word: tuple[Any, ...]
    inherited_top: tuple[BoundaryRoute, ...] = ()
    inherited_bottom: tuple[BoundaryRoute, ...] = ()
    max_image_letters: int = 100000
    max_operations: int = 10000000


def normalize_canonical_slider(word: Iterable[Any], *,
                               top_routes: Iterable[BoundaryRoute] = (),
                               bottom_routes: Iterable[BoundaryRoute] = (),
                               max_image_letters: int = 100000,
                               max_operations: int = 10000000) -> CanonicalSlider:
    """Collect canonical input twists through the coupon into its output chain.

    A mirrored coupon becomes Delta X_i Delta^-1, transporting both marked
    hierarchies before recognition. The final normal coupon has the prefix
    input chain and suffix output chain. The inverse lower word is read
    outward from the coupon, keeping its finite routing exterior. Inherited
    collars stay at their literal exterior cuts; they need not be comb
    morphisms. The value record does not recover separate base parities.
    """
    diagram = tuple(word)
    at_list = [idx for idx, token in enumerate(diagram) if type(token) is not int]
    if len(at_list) != 1:
        raise RibbonError("one full-rank standard slider is required for twist collection")
    at = at_list[0]
    token = diagram[at]
    if (not isinstance(token, tuple) or len(token) != 6 or token[0] != "slider"
            or type(token[1]) is not int or token[1] < 1
            or type(token[2]) is not int or token[2] != 1
            or any(type(x) is not bool for x in token[3:]) or token[4] or token[5]):
        raise RibbonError("normal or mirrored full-rank slider at the first port required")
    n = token[1]
    inherited_top = _boundary_routes(top_routes, n)
    inherited_bottom = _boundary_routes(bottom_routes, n)
    top_boundary = tuple(x for factor in inherited_top for x in factor.word)
    bottom_boundary = tuple(x for factor in inherited_bottom for x in factor.word)
    end = len(diagram) - len(bottom_boundary)
    if (len(top_boundary) > at or end <= at or diagram[:len(top_boundary)] != top_boundary
            or diagram[end:] != bottom_boundary):
        raise RibbonError("inherited boundary routes must be the exact exterior prefix and suffix")
    ports = tuple(range(1, n + 1))
    prefix, suffix = ports, tuple(reversed(ports))
    input_tree, output_tree = prefix, suffix
    upper_word = _word(diagram[len(top_boundary):at], n)
    lower_inverse = _inverse_braid(_word(diagram[at + 1:end], n))
    if token[3]:
        # U (Delta X Delta^-1) V: both outward body readings acquire Delta.
        # Its reversal transports every descendant label, not only the root.
        half_twist = tuple(x for width in range(n - 1, 0, -1)
                           for x in range(1, width + 1))
        upper_word += half_twist
        lower_inverse += half_twist
    normal_token = ("slider", n, 1, False, False, False)
    upper = recognize_nested(upper_word, labels=output_tree, initial_ports=ports,
                             max_image_letters=max_image_letters, max_operations=max_operations)
    lower = recognize_nested(lower_inverse, labels=input_tree, initial_ports=ports,
                             max_image_letters=max_image_letters, max_operations=max_operations)
    top_flags = (0,) + upper.boundary_parities if n > 1 else ()
    bottom_flags = (0,) + lower.boundary_parities if n > 1 else ()
    top = top_boundary + reconstruct_nested(output_tree, top_flags, initial_ports=ports)
    bottom = _inverse_braid(reconstruct_nested(input_tree, bottom_flags, initial_ports=ports)) + bottom_boundary
    twists = tuple(a - b for a, b in zip(upper.twist_powers, lower.twist_powers))
    combined_counts = twists[:1] + tuple(2 * power for power in twists[1:])
    combined = reconstruct_nested(output_tree, combined_counts, initial_ports=ports)
    normalized = top + combined + (normal_token,) + bottom
    return CanonicalSlider(diagram, normal_token, upper, lower, twists, top, bottom, normalized,
                           inherited_top, inherited_bottom, max_image_letters, max_operations)


def check_canonical_slider(witness: CanonicalSlider) -> None:
    """Check the two ribbon witnesses, finite routing and coupon twist collection."""
    if not isinstance(witness, CanonicalSlider):
        raise RibbonError("a canonical slider witness is required")
    check_nested(witness.upper)
    check_nested(witness.lower_inverse)
    if normalize_canonical_slider(witness.word, top_routes=witness.inherited_top,
                                  bottom_routes=witness.inherited_bottom,
                                  max_image_letters=witness.max_image_letters,
                                  max_operations=witness.max_operations) != witness:
        raise RibbonError("canonical slider witness does not reconstruct its source diagram")
