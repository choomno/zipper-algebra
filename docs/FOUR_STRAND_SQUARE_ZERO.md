# A four-strand square-zero quotient over Q(q,t)

This computation treats the four-strand zipper algebra over the generic field
`Q(q,t)`, with polynomial untwistings

\[
p_2(\lambda)=(\lambda-t)^2,\qquad
p_3(\lambda)=\lambda-q^2,\qquad
p_4(\lambda)=\lambda-q^3.
\]

All valid shifts are imposed. The literal input has 31 ordered defining
relations, including inverse relations and the canonical four-strand
untwisting relation `T4`. It is constructed by `Presentation(4, "square-zero")`
and converted to the exact field `Q(q,t)`. In particular, the higher relations
use the canonical free-algebra representatives `T3` and `T4`. Literal input
from `from_untwistings` uses polynomial expressions in the lasso full-twists
instead; the inverse relations identify their ideals by invertible left
factors, but their free polynomials and presentation hashes differ. The
[mathematical specification](mathematical-specification/mathematical-specification.pdf)
defines the sliders and the free-algebra representatives of these relations.

The quotient `A` has dimension **186** over `Q(q,t)`. Its three-strand slider
`X3` vanishes. For the first two-strand slider, put

\[
w=(x_1-t)X_{2,1}.
\]

Exact reduction proves `w != 0` and `w A w = 0`. Consequently the two-sided
ideal `J = A w A` is nonzero and satisfies `J^2 = 0`, so `A` is not semisimple.
The ideal assertion follows by checking `w b w = 0` for every one of the 186
standard words `b`, then extending by linearity.

## Reproduction and evidence

Use Linux, the installed Python package, Singular 4.4.1 with letterplace
support, and GAP with GBNP 1.1.0. The recorded computation used Singular
4.4.1p5, GAP 4.15.1 and GBNP 1.1.0. Candidate production and independent verification are local research
computations; they are outside CI.

From the repository root, produce the presentation and its derivations in a
fresh, empty output directory:

```sh
python examples/four_strand_square_zero.py \
  --singular Singular \
  --output runs/four-strand-square-zero \
  --timeout 180 --memory-mib 6144
```

Then run the independent checker:

```sh
gap -A -q -r --quitonbreak examples/verify_four_strand_square_zero.g
```

A successful check prints `PASS` records for the local identities, both ideal
inclusions, completion, dimension and ideal calculation, ending with
`PASS: exact four-strand square-zero verification over Q(q,t)`. It exits with
a failure status if a required check fails.

The checker reads `runs/four-strand-square-zero/proof.txt` by default. For a
different output directory, set its input explicitly:

```sh
gap -A -q -r --quitonbreak \
  -c 'SquareZeroProofFile:="runs/another-run/proof.txt";;' \
  examples/verify_four_strand_square_zero.g
```

Singular computes an auxiliary basis with `T4` omitted, then derives a compact
candidate basis. The final proof is bound to the full 31-relation input.
The GAP checker requires every input reference in every derivation to lie
among the first 30 relations, so `T4` is never used in proving the candidate
polynomials. It then checks that all 31 input relations, including `T4`,
reduce to zero. Together these checks prove that adjoining `T4` does not
change the quotient.

The candidate is accompanied by a directed acyclic graph of local polynomial
identities. The verification chain is:

| Data or check | Mathematical role |
| --- | --- |
| 31 original relations | Identify the algebra being studied. |
| 299 auxiliary derived polynomials | Supply explicitly input-derived proof nodes. |
| 148 further derived candidate polynomials | Supply the completed reduction basis. |
| All 447 local derivation identities, using only the first 30 inputs | Prove that the candidate basis belongs to the ideal generated without `T4`. |
| Zero normal forms of all 31 input relations | Prove the reverse ideal inclusion. |
| Resolution of all critical ambiguities | Prove that the candidate is a complete Gröbner basis. |
| 186 irreducible words | Give a basis and the exact dimension of the original quotient. |
| Normal forms of `X3`, `T4`, `w` and all `w b w` | Establish the slider vanishing, redundancy and nonzero square-zero ideal. |

The two ideal inclusions identify the computed ideal with the defining ideal.
GAP/GBNP independently checks completion, enumeration and the quotient
calculations using exact rational-function arithmetic. These checks concern
`Q(q,t)` throughout, rather than a collection of numerical parameter values.
The general proof criterion is explained in the
[completed-basis guide](GROEBNER_VERIFICATION.md).

## Files and outputs

The executable entry point is
[`examples/four_strand_square_zero.py`](../examples/four_strand_square_zero.py).
Its adjacent
[`four_strand_square_zero.sing`](../examples/four_strand_square_zero.sing)
implements compact completion and local derivation recording. The Python entry
point constructs the exact presentation using the maintained field and algebra
interfaces and emits the Singular calculation and local proof data. The
independent checker is
[`examples/verify_four_strand_square_zero.g`](../examples/verify_four_strand_square_zero.g):
it reconstructs the defining relations independently and checks all derivations,
completion, standard words and the stated quotient identities.

Generated files belong under the ignored `runs/` directory. A run retains:

| File | Purpose |
| --- | --- |
| `presentation.json` | The full domain-bound 31-relation input. |
| `certificate.json` | The candidate polynomials and 447-node derivation graph. |
| `proof.txt` | Framed exact data for the independent GAP check. |
| `w-square-identity.json` | An exact cofactor identity for `w^2`, verified by the producer. |
| `environment.sing`, `environment.log` | Runtime/library identification, version information and the `freegb.lib` source hash. |
| `candidate.log`, `derivation.log`, `bridge.log` | Solver output for candidate production, the 299 auxiliary derivations, and derivation of all 148 candidate polynomials from those nodes. |
| `report.json` | Recorded production results and counts. |
| `provenance.json` | Source revision, local changes and relevant source/input hashes. |
| Generated `.sing` scripts | The exact Singular programs for the run. |

The producer probes the selected Singular runtime and resolves `freegb.lib`
before the computation. Generated scripts load that resolved library path; the
report and provenance record its version and source hash.

The maintained Python completion verifier accepts the domain-bound certificate
format, but full verification of this 447-node certificate did not finish in a
60-second local probe. The documented completed verification route is
Singular production followed by the independent GAP check. A successful producer
exit or a parsed certificate alone supplies no completion verdict.

## References

1. Wolfram Decker, Gert-Martin Greuel, Gerhard Pfister and Hans Schönemann,
   *Singular*, version 4.4.1, computer algebra system.
   [Official project site](https://www.singular.uni-kl.de/);
   [citation information](https://www.singular.uni-kl.de/index.php/how-to-cite-singular).
2. A. M. Cohen and J. W. Knopper, *GBNP*, version 1.1.0, GAP package.
   [Official function manual](https://gap-packages.github.io/gbnp/doc/chap3.html).
3. George M. Bergman, *The diamond lemma for ring theory*, Advances in
   Mathematics **29** (1978), 178–218.
   [DOI](https://doi.org/10.1016/0001-8708(78)90010-5);
   [author’s corrections](https://math.berkeley.edu/~gbergman/papers/updates/diamond.html).
