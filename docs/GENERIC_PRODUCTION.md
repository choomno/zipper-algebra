# Generic finite-algebra production

The production interface starts with a literal `AlgebraPresentation` over
`Q(parameters)`, acquires a numerical f4ncgb proof, reconstructs generic shared
proofs with FiniteFlow, and independently certifies a finite word basis and
generator actions. Complete Gröbner reduction is an optional additional result. The parameter names and polynomial untwistings are
unrestricted by the historical `q,t` certificate format.

Install the package with its symbolic extra, the patched
[f4ncgb solver](../README.md#algebra-solver-and-citation), and the pinned
[FiniteFlow bindings](RECONSTRUCTION.md). f4ncgb remains the basis producer.
FiniteFlow solves the local coefficient systems; it does not perform
noncommutative completion.

## From a presentation to certificates

```python
from pathlib import Path
from zipper_algebra.algebra import AlgebraPresentation
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.execution.config import RunConfig
from zipper_algebra.execution.production import produce_algebra

K = RationalFunctionField(("q", "c0", "c1", "r"))
q, c0, c1, r = (K.parameter(name) for name in K.parameters)
P = AlgebraPresentation.zipper(
    3, field=K, q=q, untwistings={2: (c0, c1, 1), 3: (-r, 1)},
)
config = RunConfig(
    binary=Path("/absolute/path/to/f4ncgb"),
    output_root=Path("runs/quadratic-three-strand"),
    timeout_seconds=60,
    solver_io="sqlite",  # Requires the solver's optional SQLite build.
)
result = produce_algebra(P, {"q": 5, "c0": 7, "c1": 11, "r": 13}, config, threads=2)
```

Here the defining untwisting is `z^2+c1*z+c0` at level two, and `z-r` at level three. The result concerns the algebra over **Q(q,c0,c1,r)**.
The rational point proposes a proof shape; it does not change the domain of an
accepted result. A pole, unsupported proof shape, failed finite-basis check or
exhausted bound produces `status="inconclusive"`, with no dimension claim.

For a runnable linear or quadratic example:

```sh
python examples/generic_production.py --solver /absolute/path/to/f4ncgb \
  --output runs/generic-three-strand --degree 1 --threads 2
```

The linear example imposes only the level-two polynomial `z+c0`.
Use `--degree 2` for the quadratic/linear family above. Fresh independently
certified runs give dimensions 19 and 24 respectively. These are local research computations,
not CI workloads. A single numerical proof shape is not guaranteed to lift. Different native
pivot choices can produce different valid certificates and generation hashes.

## What is checked

1. **Numerical acquisition.** `run_algebra_point` substitutes every named
   parameter by an exact integer or `Fraction`, serializes rational coefficients,
   and independently checks every local numerical proof identity. It rechecks
   cache hits. Zero specialized input rows are omitted from native input and an
   explicit index map restores their positions in the original presentation.
2. **Generic rewrite proposals.** `reconstruct_basis` processes the shared
   acyclic proof graph in dependency order. Columns are contextual original
   relations or already reconstructed parent polynomials, followed by negative
   candidate monomials. `reconstruct_kernel_candidate` supplies native vector
   proposals for a polynomial and its derivation together. No local exact-kernel,
   rank or kernel-completeness claim is inferred from these vectors. Numerical tails are tried first; if necessary, contextual source
   words up to the proposed leading word enlarge the tails. Higher words remain
   cancellation equations. Numerical node coefficients select compatible generic
   kernel combinations at the anchor when needed. They guide the search; they
   do not establish generic identities from one evaluation. The returned proposal
   makes no completion or dimension assertion.
3. **Generator actions and exact certification.** `word_basis_from_rewrites`
   enumerates all words irreducible under the proposed leading words. Degree-lex
   orientation guarantees terminating reductions; confluence is not assumed.
   Shared identities express each rewritten word minus its selected reduction
   as an element of the original ideal. The independent finite-basis checker
   checks the complete generic proof graph once, every generator-column identity,
   every original relation on the full coordinate space, and cyclic independence.
   These checks prove both spanning and independence of the actual presented
   algebra. A bounded enumeration that does not finish gives no infinity verdict.
   The certified actions feed the existing exact operations and
   [minimal-polynomial computations](MINIMAL_POLYNOMIALS.md).

The default route omits local rank witnesses and kernel certification: only the
resulting local ideal identities are needed, and the final finite-basis checker
verifies them exactly. It also omits a preliminary generic-graph check and
Gröbner ambiguity tests. The separate `reconstruct_kernel` API still certifies
complete exact kernels for matrix invariants, minimal polynomials and cofactor
discovery. To request complete Gröbner reduction as well, pass
`verify_completion=True` to `produce_algebra`. This runs one additional
independent completion check. The separate `reconstruct_groebner` API likewise
reconstructs a proposal and explicitly checks both ideal inclusions and all
inclusion and overlap ambiguities.

A verified default report has `word_basis_verified=True` and
`groebner_complete=None`, meaning completion was not separately checked. An
explicitly completed report has `groebner_complete=True`. Finite-algebra
certification is exact in both cases.

The dimension and certificates are statements over the declared rational-function
field. They do not assert a uniform dimension at specializations or describe all
exceptional fibers. Shared derivations remain shared; no expansion into original
cofactor contexts or cubic multiplication table is required.

## Outputs and independent audit

One run folder contains `samples.sqlite3`, `result.json`, and a selected immutable
`algebras/<generation>/` bundle with:

- `presentation.json`: literal coefficient field, alphabet and defining relations;
- `rewrite-candidate.json`: generic rewrite proposal and shared ideal derivations;
- `word-basis-certificate.json`: finite words, generator actions and spanning proofs;
- `production.json`: solver generation, reconstruction controls/dependency hashes
  stage timings, verification counts, and optional completion results;
- `provenance.json`: producing commit, uncommitted changes and source/output hashes.

The folder receives an ignore policy automatically. Numerical samples share one
SQLite store. Certificate files can be checked without f4ncgb or FiniteFlow:

```python
import json
from zipper_algebra.verification.finite import WordBasisCertificate, verify_word_basis

folder = Path(result["artifacts"]["presentation.json"]).parent
P = AlgebraPresentation.from_json(json.loads((folder / "presentation.json").read_text()))
W = WordBasisCertificate.from_json(json.loads((folder / "word-basis-certificate.json").read_text()))
A = verify_word_basis(P, W)
assert A.dimension == result["dimension"]
```

The finite-basis certificate includes its entire ideal derivation graph, so the
independent audit above needs only the literal presentation and that certificate.
To audit the additional Gröbner claim separately, decode `rewrite-candidate.json`
with `GroebnerCertificate.from_json` and call `verify_groebner(P, candidate)`.
That separate operation checks completion even when the production run did not
request it.

The source path consists of `execution/general.py` and the common
`execution/runner.py`, `native_io.py`, `store.py`, `artifacts.py` for acquisition;
`production.py` for shared generic reconstruction; `reconstruction.py` and
`verification/linear.py` for the separate certified kernel interface; `action_production.py` for action
construction; and `execution/production.py` for orchestration/publication.
`algebra.py`, `fields.py`, `matrices.py` and `provenance.py` supply the underlying
data. `verification/f4ncgb.py`, `derivations.py`, `groebner.py`, `finite.py` and
`working_field.py` decode and independently check the mathematics.

## Limits

The solver has a subprocess timeout and explicit iteration/degree limits.
`max_operations` separately bounds numerical graph checking, generic matrix
formation, action construction and finite-basis verification. The final finite
check includes the generic proof graph; it is not an additional preliminary pass. `max_nodes` bounds retained proof nodes; `max_words` bounds word
basis enumeration. Native prime/degree/thread controls apply per solve. Final
generic identity and representation products are covered by `max_operations`;
local candidate reconstruction has no kernel-verification product budget. These operation/search bounds do not impose an operating-system
RAM cap or a wall-clock limit on reconstruction. A census supervisor must supply
those, along with bounded retries and resumable case records.

Empty alphabets and inputs whose specialized relations are all zero currently
require a direct certificate route rather than native acquisition.

`result.json` records acquisition, reconstruction, actions/certification and
publication timings. Its immutable `production.json` records the stages completed
before publication. Optional completion has its own timing and reduction budget.
