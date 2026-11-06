# Packed exact GPU rewrite expansion

This prototype executes contextual catalog expansion on an OpenCL GPU. It
generates successor braid words from reusable symbolic programs and propagates
exact Laurent coefficients. The CPU admits catalog occurrences and handles
geometric transport, recursive scheduling, and final collection.

The running CPU census uses its frozen CPU implementation. Importing the
prototype does not initialize a GPU or change the scheduler.

## What executes on the device

The compiler retains integer-expression bytecode and primitive word emitters
for generators, repeated blocks, half twists, named braid words, and slider
tokens. Concatenation and inverse, reverse, bar, and reflection operations are
compiled with their literal ordering and orientation data. Parameters remain
separate from the retained program, allowing the same program to run at
different slider widths and positions.

The GPU evaluates the symbolic registers, generates each successor's local
word, and copies the original exterior words around it. No CPU-expanded RHS
word is supplied as a device input. Sliders remain explicit tokens with size,
position, mirrored, reversed, and barred data; there is no fixed slider-size
cutoff.

All bundled catalog factors are plus or minus one. The device copies or
negates arbitrarily wide signed coefficients using unsigned 64-bit limbs and
applies checked Laurent-exponent shifts. Coordinate and exponent arithmetic
uses checked signed 64-bit values; an unsupported or overflowing input fails
before acceptance. It does not wrap, approximate, or silently fall back to CPU
arithmetic. Final duplicate collection uses exact Python integers and retains
the reference interpreter's monomial and successor ordering.

## Reproduce the workload comparison

Install the GPU vendor's OpenCL runtime and the optional Python dependencies:

```sh
python -m pip install -e '.[gpu-opencl]'
```

Capture actual local catalog calls from bounded zipper calculations into a new
SQLite store, then run the same stage on each backend:

```sh
python experiments/zipper_gpu_rewrite_benchmark.py capture /tmp/zipper-rewrites.sqlite3 --batch-size 1024 --max-steps 20000 --max-terms 10000
python experiments/zipper_gpu_rewrite_benchmark.py benchmark /tmp/zipper-rewrites.sqlite3 --backend cpu --repeats 5
python experiments/zipper_gpu_rewrite_benchmark.py benchmark /tmp/zipper-rewrites.sqlite3 --backend opencl --repeats 5
```

With `g++` available, the same packed programs can run on one native CPU thread:

```sh
python experiments/zipper_gpu_rewrite_benchmark.py benchmark /tmp/zipper-rewrites.sqlite3 --backend cpp --repeats 5
```

The native library builds into a source/compiler/flags-bound local cache.
Its manifest checks the source, binary, and ABI identities. Build and load
times are reported separately. This C++ comparator shares Python admission,
packing, reconstruction, collection, and verification with the GPU path; it
does not implement the complete scheduler in C++.

The default workloads are an alternating three-strand word of length ten, a
positive five-strand word of length six, and a positive eight-strand word of
length seven. Complete calculations and checked incomplete prefixes are
reported separately. `--workloads` accepts a JSON list of `name`, `n`, and
`word` objects. `--roots` selects saved workload IDs.

The default `--kind contextual` measures actual contextual interpreter calls.
The capture also records single-strand instantiations: `--kind used` includes
used single-strand expansions, and `--kind all` includes current instantiations
whose results the controller discards. These single-strand cases are expressed
as equivalent contextual expansions for the comparison. They are labeled
separately from the actual interpreter-call measurements.

Sequential mode preserves the batches available during one calculation.
`--mode pooled` is an offline throughput diagnostic that combines different
scheduler times.

## Exactness and timing

Every case retains input-root, case, scheduler-batch, source-template,
parameter, predecessor, exterior, control, and coefficient-weight identities.
The collector's primitive calls use unit coefficients; capture restores the
actual selected term's coefficient using coefficient linearity and checks the
scaled result against independent contextual replay. Only committed batches
are retained. Final shadow states reconcile with the full scheduler
certificate, which is JSON-reloaded and independently replayed.

Both timed paths use the same independent CPU interpreter for verification and
compare the complete selected output with the captured result. GPU totals
include RHS-program compilation/cache preparation on first use, host admission and
packing, template uploads, allocation, all transfers, register evaluation,
word expansion, coefficient propagation, reconstruction, collection, and
verification. Device initialization and kernel compilation are separately
reported so cold-start totals can include them. Capture, database loading,
integrity checking, and original full-calculation replay are separate costs.

The CPU comparator is the existing, unmemoized contextual interpreter. A stage
ratio includes the effects of reusable compilation as well as GPU execution.
It does not establish a full-census speedup: the production memoized executor
can avoid local expansions altogether, and geometric control still runs on CPU.
The prototype currently targets OpenCL; the older coefficient-collection
prototype's CUDA backend does not provide CUDA rewrite expansion.

Captured JSON artifacts retain the existing decimal-integer conversion limit.
Direct API wide-integer tests avoid decimal serialization. No global integer
conversion setting is changed.
