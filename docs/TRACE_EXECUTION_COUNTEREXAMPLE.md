# Zipper execution does not descend to trace classes

**Status: exact computed counterexample, November 7, 2026.** Distant-generator
commutation does not preserve the existing procedure's formal terminal sum,
transition count, or bounded-resource completion status. This is not a
counterexample to equality in the zipper algebra or to eventual termination.

## A connected, four-letter counterexample

In B4, take

    u = sigma_1 sigma_3^-1 sigma_2 sigma_3^-1,
    v = sigma_3^-1 sigma_1 sigma_2 sigma_3^-1.

They differ by one permitted distant commutation and both have full connected
generator support. The signed integer inputs are `(1,-3,2,-3)` and
`(-3,1,2,-3)`.

| Input | Logical transitions | Saved terminal terms |
|---|---:|---:|
| u | 342 | 31 |
| v | 503 | 31 |

The formal outputs differ even after forgetting routing words, orientations,
and normalization provenance. Retain only the ordered component sizes,
original top/bottom endpoint sets, and integer twist vectors, and collect
coefficients in Z[q,q^-1]. For the signature consisting of one size-four
component, top and bottom endpoint sets `{1,2,3,4}`, and twist vector
`(0,0,-1)`, the coefficient of `q^-2` is **0 for u and +1 for v**.
There are twelve differing Laurent monomials in this coarser comparison.
Thus the discrepancy is not just a different spelling of ambient routes.

The terminal family has not been established as a linearly independent basis
or a unique algebraic normal form. Different formal sums here therefore do
not imply different elements of the zipper algebra. Both inputs are equal
already in the braid group.

## Other exact differences

The full-support pair `(-2,3,-1)` and `(-2,-1,3)` takes 38 and 39 transitions.
At `max_steps=38`, the first completes while the second remains incomplete:
its next atomic batch would exceed the budget. This is a bounded-status
difference, not evidence of nontermination.

The pair `(1,2,3,1,2,3)` and its lexicographic trace representative
`(1,2,1,3,2,3)` takes 35,308 and 7,728 transitions, with 218 and 147 saved
terminal terms. Lexicographic selection can therefore substantially change
cost. The four-letter example shows that the lexicographic representative
can also be the more expensive member.

An exhaustive B4 comparison of noncanonical trace representatives at lengths
two and three found no intrinsic-output mismatch (4 and 48 comparisons,
respectively). Length four is minimal for this output discrepancy in B4;
the search stopped after finding the counterexample and did not classify all
length-four words. Cost differences already occur at length two.

## Verification and reproduction

All six certificates for the primary length-two, length-four, and length-six
pairs passed both independent C++ replay and independent Python replay.
The connected length-three and budget-38 examples additionally passed C++
replay. Arithmetic is exact over Z[q,q^-1].

The calculations use the producer and verifier implementations at commit
`5b08c087b71e1b4ebfb811f69f9c5ca077d38d78`, unchanged by the enumeration work,
catalog SHA-256
`27c5dd71d09718a5561d18c802f714e07abac87890dc87169215369541462470`,
batch size 64, one worker, and traces disabled. Limits are one million steps,
100,000 terms, and expansion budget one million; other `Limits` defaults apply.
The archived length-four certificate-file SHA-256 values are:

| Input | Certificate SHA-256 |
|---|---|
| u | `72dd8120e693fdb014efb3de1aaa85dfe71ac20fe510c9dd353f3fb5d579d6be` |
| v | `4f668c5681f078170e2e2fdcc185dd7a89e3e47966fe0e8eccd344bdc175656c` |

The following reproduces the distinguishing coefficient:

```python
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeSession
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier

limits = Limits(max_steps=1_000_000, max_terms=100_000,
                expansion_budget=1_000_000)
producer = NativeSession()
verifier = CppVerifier(producer_metadata=producer.metadata)
try:
    for word in ((1, -3, 2, -3), (-3, 1, 2, -3)):
        result = producer.reduce(word, n=4, limits=limits, batch_size=64)
        assert result["status"] == "complete" and verifier.check(result)
        coefficient = 0
        for term in result["output"]:
            if len(term["components"]) != 1:
                continue
            component = term["components"][0]
            record = component["record"]
            if (component["size"] == 4
                    and record["original_top_ports"] == [1, 2, 3, 4]
                    and record["original_bottom_ports"] == [1, 2, 3, 4]
                    and record["twist_powers"] == [0, 0, -1]):
                coefficient += sum(c for q, c in term["coefficient"] if q == -2)
        print(word, result["metrics"]["steps"], coefficient)
finally:
    verifier.close()
    producer.close()
```

Expected output: 342 transitions and coefficient 0 for u; 503 transitions and
coefficient 1 for v. Future changes to rewrite selection or collection may
change these execution-specific figures; the archived certificates bind the
counterexample to the recorded implementation.

## Consequence for enumeration

The [lexicographic trace enumerator](essential-enumeration.md) correctly counts
trace classes. A checked derivation for its representative can be preceded by
the input's checked commutations to give an algebraic derivation for another
member. That is a different execution policy.

One representative does not certify every existing literal execution, its
resource status, its worst-case cost, or the distribution of those costs.
Such conclusions require a separate bound over all linear extensions of the
trace. This obstruction is separate from [disconnected-block
separation](zipper-block-separation.md): the counterexample already has
connected support.
