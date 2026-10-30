# Certified matrix reconstruction with FiniteFlow

The reconstruction interface produces a complete right kernel over the rationals or the declared
rational-function field. It uses [FiniteFlow](https://github.com/peraro/finiteflow),
Tiziano Peraro's framework for finite-field computation and multivariate rational
reconstruction ([paper](https://doi.org/10.1007/JHEP07(2019)031)). The
[maintained fork](https://github.com/Choomno/finiteflow) supplies the strict
original-system checks and reconstruction-state controls required here.

The interface accepts a sparse `ExactMatrix`, reconstructs its kernel, and
returns a `KernelCertificate` only after independent exact verification. It accepts matrices arising from
action, trace-form, and contextual cofactor systems once those matrices are formed. Matrix formation from
a presentation is a separate operation. This first interface does not replace
the existing proof-sampling and coefficient-fitting discovery command.

## Installation

Install the symbolic extra in the Python environment used for the computation.
FiniteFlow is an external native dependency; installing zipper-algebra does not
build it. Its Python interface additionally needs CFFI, CMake 3.18 or newer,
C/C++ compilers, Python development headers, and GMP development headers.
For example, the Python headers are provided by `python3-devel` on Fedora or
`python3-dev` on Debian/Ubuntu. FLINT is recommended. See the
[upstream build instructions](https://github.com/peraro/finiteflow#installation)
for system dependencies.

The required native revision is `c0db1ed610b2f3a341ec0b247ec06067516ea575` from the
fork's `fix/reconstruction-reliability` branch. These native changes are
currently local; use the matching local checkout until that revision is
published. Publish that dependency before distributing this integration.

From a directory containing the two matching checkouts, with the intended Python
environment active:

```sh
python -m pip install cffi
python -m pip install -e './zipper-algebra[symbolic]'
cmake -S finiteflow -B finiteflow/build \
  -DCMAKE_BUILD_TYPE=Release -DFFLOW_MATHEMATICA=OFF -DFFLOW_PYTHON=ON \
  -DPython3_EXECUTABLE="$(command -v python)" \
  -DCMAKE_INSTALL_PREFIX="$PWD/finiteflow/install" \
  -DFFLOW_PYTHON_INSTALL="$(python -c 'import sysconfig; print(sysconfig.get_paths()["platlib"])')"
cmake --build finiteflow/build --parallel 4
cmake --install finiteflow/build
python -c 'import fflow; assert hasattr(fflow, "LSolveValidateFullSystem")'
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
bounds scalar products in independent verification. Search bounds do not impose
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
retains a scratch matrix in each native solver-data clone; account for that
per-thread storage when choosing worker and thread counts.
