# Exact GPU coefficient-collection prototype

The optional prototype collects Laurent coefficients from actual zipper rewrite
batches on an OpenCL or CUDA GPU. The existing CPU census continues to use its
frozen CPU implementation. Full zipper rewriting and continuation-state
canonicalization remain on the CPU.

## Arithmetic and identities

Each contribution contains an input-root ID, complete continuation-state ID,
Laurent exponent, and signed integer coefficient. Input roots and scheduler
batch namespaces are distinct; different input words or different scheduler
times cannot cancel each other's coefficients. The capture includes incoming
terms and affected prior accumulators, and compares its expected answer with
the scheduler's actual collection operation.

The host sorts the keys and packs coefficients into variable-width unsigned
64-bit limbs using signed two's-complement encoding. The width is chosen from
the largest coefficient bit length and contribution count, with enough space
for every intermediate partial sum. One GPU work item accumulates each group
with explicit carry propagation. Signed decoding recovers the exact integer.
No floating-point coefficients, fixed-precision truncation, or silent CPU
arithmetic fallback is used. Missing devices and compilation failures raise
distinct errors.

## Run a comparison

Install the appropriate GPU vendor runtime first. The OpenCL backend enumerates
GPU devices only. Install the optional Python dependencies in a separate
environment:

```sh
python -m pip install -e '.[gpu-opencl]'
```

Capture real workloads into a new SQLite database, then benchmark the same
captured batches on both backends:

```sh
python experiments/zipper_gpu_benchmark.py capture /tmp/zipper-frontiers.sqlite3 --max-steps 20000
python experiments/zipper_gpu_benchmark.py benchmark /tmp/zipper-frontiers.sqlite3 --backend cpu --repeats 5
python experiments/zipper_gpu_benchmark.py benchmark /tmp/zipper-frontiers.sqlite3 --backend opencl --device 0 --repeats 5
```

The default workloads include branching sign patterns and several strand
counts. A bounded workload can produce a verified incomplete prefix; status
is retained explicitly. Supply `--workloads` with a JSON list of objects with
`name`, `n`, and `word` fields to select other literal inputs. `--roots`
selects saved workload IDs. Capturing refuses an existing nonempty database.

The CUDA backend uses the `gpu-cuda` extra and `--backend cuda`. Actual hardware
validation currently covers the local Intel Arc OpenCL GPU. The CUDA transport
contract has unit tests; no NVIDIA hardware benchmark has been performed.

## Timing and interpretation

The default `sequential` mode processes each captured scheduler batch in order.
Reported totals include input validation, host sorting and packing, allocation,
host-to-device copies, exact reduction, device-to-host copies, signed
reconstruction, the same independent Python-integer verifier on both paths,
and comparison with the captured scheduler output. Initialization and kernel
compilation, database loading and integrity checks, and original rewrite/capture
time are reported separately. GPU operations synchronize before their timing
interval ends. Constructor metadata includes device, runtime, and kernel hashes.

The timed stage starts with complete continuation identities already encoded.
Its result is a coefficient-collection comparison; it does not establish a
speedup of a full zipper calculation. `--mode pooled` combines different saved
scheduler times and is an offline throughput diagnostic. Such frontiers are not
all simultaneously available during one live calculation.

The Python API accepts a retained collector so paired CPU/GPU trials can
alternate order without recompiling the kernel. SQLite stores compressed
captures and certificates, source identities, output hashes, per-batch timings,
and benchmark records. Do not run the benchmark against the live census store.

## Validation

CPU packing and carry tests cover cancellation, distinct roots/states/exponents,
malformed input, and coefficients wider than a machine word. Actual OpenCL tests
exercise the same kernel on the GPU, including 32,769-bit cancellation. Mock
driver tests separately check transport, allocation, launch, and error contracts;
they are not evidence of NVIDIA hardware execution. Independent certificate
replay validates the original zipper captures.

Capture uses the existing JSON certificate representation. Extremely large
decimal integers can reach Python's JSON conversion limit and fail closed;
the direct limb-collection API has no decimal conversion step. The wide-integer
hardware probes exercise that direct API. No global integer-conversion setting
is changed.
