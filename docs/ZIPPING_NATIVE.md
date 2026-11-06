# Native zipper state engine

The optional C++ engine executes complete literal-word zipper calculations.
It retains native braid-token vectors, Laurent polynomials over arbitrary-precision
integers, control states, ordered frontiers, and transition caches throughout a
calculation. Python sends one input request and receives one final certificate.
Optional traces retain the committed batches and their local rewrite evidence.
Small event and terminal-record objects remain native JSON values; they do not
require a Python callback or a text round trip between transitions.

The existing Python collected executor is the independent replay reference.
Native arithmetic results are accepted by comparing the full output, pending
frontier, diagnostics, deterministic metrics, and enabled trace with that replay.
The native backend has its own source and executable provenance.

## Build and run

A GCC-compatible C++17 compiler, Boost headers, and OpenSSL development
headers/library are required. The Python process wrapper currently targets POSIX. The pinned JSON header and its license are included. Importing the
module does not compile code; constructing a session builds a content-addressed
executable when no matching build is present.

```python
from zipper_algebra.zipping_algorithm.native import NativeSession
from zipper_algebra.zipping_algorithm.interpreter import Limits

with NativeSession(cache_entries=262_144, cache_bytes=256 * 1024**2) as session:
    result = session.reduce(
        [1, -2] * 5, n=3,
        limits=Limits(max_steps=100_000, max_terms=100_000),
        trace=False,
    )
    session.check(result)
```

Input words include inverse pairs and are not identified up to braid isotopy.
All eligible slider sizes are available. Resource limits yield an exact stopped
frontier with atomic batch rollback. Coordinate and Laurent-exponent arithmetic
uses checked signed 64-bit integers; integer coefficients use arbitrary precision.

Before running pytest against a cold checkout, build the engine outside the
per-test timeout:

```sh
PYTHONPATH=src python -c 'from zipper_algebra.zipping_algorithm.native import build_engine; build_engine()'
```

A session retains the parsed catalog and native transition cache across inputs.
Cache keys retain rank, active geometry, control, boundary factors, limits, and
trace policy. Already extracted components, endpoint label values, coefficients,
and depth are restored from the caller when a unit transition is reused. This
includes complement steps; it does not identify different literal input words.
The cache byte counter estimates typed payload and keys and excludes allocator
overhead, frontiers, optional traces, and other process memory.

Each input has its own formal sum. Independent word processes can run separate
sessions; the first implementation executes each native batch in one thread.
Cache budgets and timing telemetry are execution data, not mathematical output.

## Benchmark protocol

`experiments/zipper_native_engine_benchmark.py` compares full reductions with
both Python collected and memoized sessions. Its SQLite database retains source
identities, compressed certificates, reduction timings, independent verification
timings, and their sums. Native reduction time includes request encoding, process
communication, native computation, final certificate encoding, and response
decoding. Compilation and initialization are recorded separately.

Measurements must identify the source epoch, limits, batch size, trace policy,
and cache policy. Complete outputs and bounded prefixes are distinguished.
Concurrent census load and repeated-input cache reuse affect timings; neither
should be interpreted as a universal speedup or complexity theorem.

## Implementation

- `native/zipper_engine/engine.hpp`: native tokens, coefficients and pass controls.
- `core.cpp`: exact catalog evaluation and contextual rewrites.
- `framed.cpp`, `transport.cpp`: traversal, frame transport and turnaround.
- `extraction.cpp`: component splitting and nested ribbon normalization.
- `algorithm.cpp`: complete algorithm phase transitions.
- `scheduler.cpp`: stable atomic batches, exact collection and transition reuse.
- `main.cpp`: persistent process interface.

The frozen Python engine and running census remain independently usable.
