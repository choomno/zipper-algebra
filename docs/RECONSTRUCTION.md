# Certified matrix reconstruction with FiniteFlow

The reconstruction interface produces a complete right kernel over the rationals or the declared
rational-function field. It uses [FiniteFlow](https://github.com/peraro/finiteflow),
Tiziano Peraro's framework for finite-field computation and multivariate rational
reconstruction ([paper](https://doi.org/10.1007/JHEP07(2019)031)). The
[maintained fork](https://github.com/Choomno/finiteflow) supplies the strict
original-system checks and reconstruction-state controls required here.

The interface accepts a sparse `ExactMatrix`, reconstructs its kernel, and
returns a `KernelCertificate` only after independent exact verification. It accepts matrices arising from
action, trace-form, and contextual cofactor systems once those matrices are formed. The installed `discover` command forms a contextual cofactor matrix from a fresh
f4ncgb proof graph and uses this interface for its scalar reconstruction stage.
The [reproduction guide](X3_REPRODUCTION.md) describes that complete route.

## Installation

Install the symbolic extra in the Python environment used for the computation.
FiniteFlow is an external native dependency; installing zipper-algebra does not
build it. Its Python interface additionally needs CFFI, CMake 3.18 or newer,
C/C++ compilers, Python development headers, and GMP development headers.
For example, the Python headers are provided by `python3-devel` on Fedora or
`python3-dev` on Debian/Ubuntu. FLINT is recommended. See the
[upstream build instructions](https://github.com/peraro/finiteflow#installation)
for system dependencies.

The required native revision is `a2a5abb9bf90d15416116c0a36a8821e04cbdb7b` from the
fork's `build/python-module-integration` branch, developed in
[draft PR #5](https://github.com/choomno/finiteflow/pull/5). Check out this exact
revision so the installed APIs match the integration tests.

From the directory containing your zipper-algebra checkout, with the intended
Python environment active:

```sh
git clone https://github.com/choomno/finiteflow.git
git -C finiteflow checkout --detach a2a5abb9bf90d15416116c0a36a8821e04cbdb7b
python -m pip install cffi
python -m pip install -e './zipper-algebra[symbolic]'
cmake -S finiteflow -B finiteflow/build \
  -DCMAKE_BUILD_TYPE=Release -DFFLOW_MATHEMATICA=OFF -DFFLOW_PYTHON=ON \
  -DPython3_EXECUTABLE="$(command -v python)" \
  -DCMAKE_INSTALL_PREFIX="$PWD/finiteflow/install" \
  -DFFLOW_PYTHON_INSTALL="$(python -c 'import sysconfig; print(sysconfig.get_paths()["platlib"])')"
cmake --build finiteflow/build --parallel 4
cmake --install finiteflow/build
python -c 'from zipper_algebra.execution.dependencies import reconstruction_provenance; print(reconstruction_provenance()["implementation"])'
```

The producer rejects an `fflow` module without `LSolveValidateFullSystem`. Checking an already supplied matrix/kernel certificate does not
load FiniteFlow. Rational-field checking with a witness prime below `2**32`
uses the standard library; parameter fields and larger witness primes use SymPy.

## Python interface

```python
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.matrices import ExactMatrix
from zipper_algebra.reconstruction import reconstruct_kernel
from zipper_algebra.verification.linear import KernelCertificate, verify_kernel

field = RationalFunctionField(("q", "t"))
q, t = (field.parameter(name) for name in field.parameters)
matrix = ExactMatrix.from_rows(field, [[1, q, t], [0, 1, q]])
certificate = reconstruct_kernel(matrix, threads=4)
assert certificate.rank == 2
assert certificate.nullity == 1

# Reload the mathematical witness, then check against the original matrix.
reloaded = KernelCertificate.from_json(certificate.to_json())
verify_kernel(matrix, reloaded)
```

This small linear example has kernel vector `(q*q-t, -q, 1)`. It illustrates the
interface rather than a new algebra computation. `ExactMatrix` stores entries
as `(row, column, coefficient)` triples with zero-based coordinates, declared
dimensions, and a complete coefficient field. Repeated coordinates and mixed
fields are rejected; the matrix acts on column vectors. `from_rows` accepts
exact scalar entries and field elements.

The producer supports the rationals and rational-function fields over the
rationals. Parameter names and coefficients are passed as structured monomial
data, without expression parsing. Its settings include `threads`, `max_primes`,
`max_degree`, `start_prime`, and an optional explicit `shift`. `max_products`
bounds scalar products in independent verification. The modular inverse check
multiplies only nonzero entries of the inverse and specialized minor. The
identity block is checked from its stored coefficients. Search bounds do not impose
wall-time or resident-memory limits; census workers need separate process
resource limits. Exhausted reconstruction or witness searches raise an error
rather than return a partial result.

The verifier checks the complete generic product of the original matrix with
the proposed kernel, an identity block proving kernel-column independence, and
a denominator-safe invertible modular minor proving the opposite rank bound.
It checks every stored input coefficient at the minor's specialization point.
No library success status or collection of matching samples substitutes for
these checks. The conclusion is over the generic declared field; ranks may
drop and specialized kernels may be larger at particular parameter values.
`verify_kernel` returns the supplied certificate on success and raises an
exception for invalid data, failed identities, or exhausted operation budgets. The [compiled mathematical reference](mathematical-specification/mathematical-specification.pdf)
states the rank argument precisely.

`zipper-matrix-v1` and `zipper-kernel-certificate-v1` JSON data preserve the field,
shape, exact coefficients and witness. The verifier checks the certificate's SHA-256 against
the complete canonical data of the supplied matrix. These mathematical encodings do not write
files. CMake installs native files separately from pip; retain the configured
installation prefix so the Python extension can find its library. The fork's
CMake `uninstall` target removes its installed files. Saved computation workflows should also retain source and dependency
identities through the project's [provenance facilities](PROVENANCE.md).

For census concurrency, use independent worker processes created with the
`spawn` start method. The interface serializes its own simultaneous Python
callers because FiniteFlow uses a process-global session. `threads` controls
parallel evaluations within one native call. Raw `fflow` calls must not
interleave with an active reconstruction. Strict original-system validation
materializes dense original-system and dependent-column buffers and performs a
dense rank calculation at each modular evaluation, including for sparse input.
Each native evaluator clone retains its own checking buffers. Account for this
per-thread memory and computation when choosing worker and thread counts.

## Discovery from a proof graph

`discover` uses the chosen anchor only to obtain candidate word contexts. It
expands the proof graph to triples `(left, relation index, right)`, retaining
contextual occurrences before deduplication and ignoring their numerical
coefficients. For each distinct triple it expands the original symbolic
relation between those words. The resulting columns form `A`; the generic target
is the last column of the augmented matrix `[A|-h]`.

The complete reconstructed kernel is checked independently against that matrix.
A kernel column with nonzero last coordinate provides a cofactor representation
after division by that coordinate. A separate free-algebra checker then expands
the original relations with those coefficients and compares the result with the
target. Both checks are required before a verified bundle is published.

A sampled proof graph can omit contexts needed generically. Failure to find a
kernel vector with nonzero target coordinate is inconclusive for membership;
changing the anchor can supply a different support. The command does not claim a
complete method for finding supports or a completed Gröbner basis. Native prime,
degree and thread controls govern the scalar computation; `--max-contexts`
bounds proof expansion. `--max-products` separately bounds support construction
(context substitutions and matrix placements) and scalar products in exact
kernel verification. It does not bound native FiniteFlow work or the final
free-algebra check. These are separate from the f4ncgb process timeout and
iteration limits.

The published bundle contains the context/word lists, the augmented matrix,
the complete kernel witness, and the independently checked cofactor. It also
hashes the actual loaded FiniteFlow Python module, CFFI extension and shared
library. These byte identities describe the installed dependencies; they do not
assert that an arbitrary local build corresponds to a claimed source revision.
