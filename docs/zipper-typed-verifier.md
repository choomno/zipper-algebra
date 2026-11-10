# Independent typed verifier

The independent verifier uses owned C++ values throughout collected replay and
completed-continuation replay. It shares no mathematical implementation with the
producer. `typed_word` represents a token as a signed generator or a slider with
all five slider fields. `typed_state` represents the five algorithm phases as a
variant and stores exact Laurent coefficients, endpoint labels, boundary factors,
extracted contexts and conditional component records as concrete fields.
`typed_controller` owns upward and framed branch controls, initial slider records,
processed rule events and frame geometry. `typed_transport`, `typed_ribbons` and
`typed_extraction` implement transport and extraction using typed witnesses.

`ReplaySession` stores typed work, transitions and output in its frontier and
caches. Its restoration records and unary continuation links also own typed
values. The direct canonical encoder preserves the complete continuation
identity, excluding precisely the scalar and depth fields specified by the
collection contract. Optional component fields retain their presence. The word
size observer performs every reflection and shift in the original order,
including operations that cancel in the final coordinates.

The unary scheduler stores a compact sequence of continuation links for a
single-child transition with no direct output. Each logical transition retains
its own key, active-cycle membership, restoration, resource checks, rule count
and cache admission. Unwinding applies each restoration in order. Standalone
atomic batches and optional traces retain their existing schedule.
`ZIPPER_DISABLE_UNARY_FUSION` selects the same typed implementation with fusion
disabled for controlled comparisons; the wire protocol does not select this flag.

## JSON boundaries

JSON remains in the following explicit places:

- The native request and response protocol, source manifests, diagnostics,
  metrics and optional trace/evidence encoders.
- Catalog rule syntax, parameter environments, conditions and expanded catalog
  packets. Catalog words are decoded into typed words before owned branch state
  is constructed. Catalog code remains independent of the producer.
- Exact Artin image digest serialization in ribbon recognition, preserving the
  established certificate digest bytes.
- String escaping in the direct canonical encoder. It does not construct a JSON
  work-state tree.

The historical `transition` operation retains the unrestricted JSON diagnostic
contract, including malformed or extended records and their original rejection
order. This is an explicit compatibility operation. Collected `replay` and
`replay_reused` do not call that dispatcher or fall back to it on errors.
`typed_transition` decodes the owned typed schema and exposes the production
controller for differential testing. It rejects malformed typed containers and
unknown phase names; its accepted input schema does not include arbitrary
historical extensions.

## Independent verification

The test fixture under `tests/native/frozen_verifier_9c477dd` authenticates a full
frozen source archive and builds a separate JSON verifier executable. Expected
transitions, traces, resource refusals and completed outputs come from that
executable. The producer supplies build identity metadata only. Local native
harnesses additionally compare typed word operations, canonical identities,
framed and band transport, turnaround, ribbon recognition and extraction against
legacy independent implementations.

Formatting of the new typed sources uses clang-format 23.1.5 with the existing
native engine formatting configuration. Performance measurements and their
source identities are recorded separately from these structural contracts.
