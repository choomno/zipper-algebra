# Independent native zipper verifier

For census results, `CppVerifier.check_reused` independently evaluates the
completed continuation graph and compares exact formal outputs. It maintains
its own completed-result cache. `VerifierIncomplete` distinguishes a resource
refusal from a different completed replay. Standalone `check` continues to
verify full execution certificates. See [census reuse](zipper-census-reuse.md).

The C++ verifier replays literal-word calculations with a separate implementation
of the catalog evaluator, control phases, transport and ribbon extraction,
collection, and scheduling. It compares the entire result: terminal terms,
pending frontier, diagnostics, deterministic metrics, and enabled trace records.
The Python reference remains available for differential checks.

## Independence

The producer retains typed packed states and uses Boost arbitrary-precision
integers. The verifier owns typed token, controller, context and formal-output records and
uses GMP integers. Its complete typed state and remaining JSON boundaries are
described in [the typed verifier contract](zipper-typed-verifier.md). Its
source and object files are separate; it does not include or link the producer's
rewriting, controller, extraction, collection, or scheduler implementation.
Only the mathematical catalog and generic libraries, including the pinned JSON
parser, are shared. Each bundle includes its own authenticated copy of that
header and its MIT license.

A verifier cache contains only transitions computed by the verifier itself.
Producer results and optional witnesses are comparison inputs. They do not
supply transition evaluations or populate the verifier's cache. The setup
binds producer source/build identity as provenance; acceptance additionally
requires independent exact replay. The finite differential test suite checks
implementation agreement and does not establish a general complexity theorem.

## Build and check

A C++17 compiler, GMP C++ development headers/libraries, and OpenSSL development
headers/library are required. On Debian or Ubuntu, the relevant packages are
`build-essential`, `libgmp-dev`, and `libssl-dev`. The producer additionally
requires Boost headers. Each native bundle is authenticated in both the wheel
and the source distribution. Builds use a separate content-addressed cache.

```python
from zipper_algebra.zipping_algorithm.native import NativeSession
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier

with NativeSession() as producer:
    result = producer.reduce([1, -2, 1], n=3, trace=True, batch_size=4)
    with CppVerifier(producer_metadata=producer.metadata) as verifier:
        assert verifier.check(result)
```

Before running tests against a cold installation, compile outside the per-test
timeout:

```sh
python -c 'from zipper_algebra.zipping_algorithm.native_verifier import build_verifier; print(build_verifier()[0])'
```

Native verification is an explicit API. The running census configuration retains
its recorded checker and source identity until an explicit, versioned migration.
The census runner accepts `--verifier cpp` and optionally
`--verifier-allocator /path/to/libtcmalloc_minimal.so.4`. Changing an existing
campaign's checker requires draining its writer, then `--migrate` with a fresh
`--migration-backup` path. Historical certificates retain their original epochs;
the new epoch binds verifier code, build, and allocator identity. The Python
checker remains available with `--verifier python`.

`check()` accepts an exactly reproduced bounded prefix as well as a completed
calculation. Inspect the certificate's `status` when completion is required.

An installed allocator library can be selected explicitly with
`CppVerifier(..., allocator_library=Path("/path/to/libtcmalloc_minimal.so.4"))`.
The wrapper records its resolved path and SHA-256, verifies that the verifier
process loaded that file, and checks for changes before subsequent requests.
Only the verifier child receives the preload setting. Without this argument,
the system allocator remains the default.

## Benchmark interpretation

Compare both verifiers on the same saved census inputs and resource limits.
Include request encoding, process transfer, independent replay, strict comparison,
and response decoding in verification time. Record cold versus retained cache
policy, source and executable identities, compilation/setup separately, memory,
and complete versus bounded-prefix outcomes. Repeated measurements under an
active census may vary with contention. A speedup measured for one input family
is not a bound for arbitrary strand counts or word lengths.

The transition cache's byte budget measures canonical keys and canonical
payload sizes; resident typed records and container allocations require additional
memory. Compare measured
process memory when sizing workers.

The benchmark captures certificates from the census through a read-only
connection and saves its frozen cases, measurements, and provenance in a separate
SQLite database:

```sh
PYTHONPATH=src:. python experiments/zipper_cpp_verifier_benchmark.py \
  /tmp/verifier-benchmark.sqlite3 --source /path/to/words.sqlite3 --passes 2
```

Use `--order cpp-first` to reverse backend order. Each backend gets a fresh
process; subsequent passes retain that backend's transition cache. Reusing the
benchmark database preserves its original cases and requires the same source
database and recorded source identities.
Pass `--allocator-library /path/to/libtcmalloc_minimal.so.4` to select and
authenticate that allocator for the C++ verifier measurements.
