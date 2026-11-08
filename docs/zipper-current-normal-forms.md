# Current-algorithm input and output conventions

This interface preserves the current zipper procedure's full formal outputs
under proved identifications. It does not choose a normal form for braid-group
isotopy. A literal word or chosen diagram remains part of the provenance:
the current procedure can give different formal outputs for braid-equivalent
words. The [modified representative procedures](zipper-census-baselines.md)
are separate experiments and are not used here.

A small obstruction already occurs in B3: `(1,2,1)` and `(2,1,2)`
represent the same braid, but current execution gives 18 and 13 fully collected
terminal terms respectively, in 158 and 127 logical transitions. Both native
runs passed independent C++ replay and matched fresh reference executions.
Thus a braid-group normal form cannot also classify this procedure's exact
outputs without changing which representative-producing procedure is used.

## Geometric and algebraic frame

- A rank-n input is a list of nonzero signed integers with absolute value <n.
- The list is in algebraic multiplication order; the algorithm traverses it
  **right to left**, from bottom to top. Reverse a list recorded chronologically
  from bottom to top before importing it.
- Number bottom slots and top slots independently 1,...,n from left to right.
  Transport original bottom strand labels through the diagram when following
  strands. A top port number denotes a top slot, not its incoming strand label.
- At positive sigma_i the strand entering position i passes over the strand
  entering i+1; the negative generator reverses that over/under choice.
- Coefficients lie in Z[q,q^-1], with q formal. Output coefficient pairs are
  `(q exponent, integer coefficient)`; no specialization is implicit.

A diagram from another convention needs an explicit projection, endpoint
numbering, and sign conversion. These conventions do not infer isotopy from
an image, and do not silently cancel inverse pairs in a higher-rank connected
word. Preserve any external endpoint names as a separate map to these slots.

## Input normal form: requests plus reconstruction

`input_normal_form(word,n=n)` produces schema `zipper-current-input-v1`,
policy `literal-support-peeling-and-rank-two-v1`. Its checked tree retains the
exact original word. It performs only:

1. Literal support separation and index translation. Each block retains its
   literal subsequence; disconnected blocks are processed right to left.
2. Certified initial singleton peeling. The peeled original bottom label,
   its original top port, and the literal complement are retained.
3. Rank-two exponent normalization. An isolated two-strand word or compatible
   recursive two-strand restart is replaced by the pure-sign word with the
   same exponent sum. [Full-output invariance is proved](zipper-rank-two-identification.md).

`check_input_normal_form` independently checks the decomposition and deletion
witnesses rather than calling the constructor. It rejects unknown policies,
malformed port maps, and changed literal subsequences. `canonical_input_json`
serializes the complete checked witness deterministically.

`leaf_requests(form)` gives the distinct primitive local requests, each with
rank, literal word, and a key. `reconstruction_plan(form)` flattens the tree
into an extraction-ordered sequence of singleton and leaf components with
original top/bottom embeddings. `reconstruction_key(form)` hashes this plan,
excluding the original word and incidental tree shape.

**Equal plans imply equal reconstructed current outputs**, under the source
separation, peeling, and rank-two theorems. The converse is not asserted.
Equal primitive request keys alone do not imply equal ambient outputs: the
endpoint maps and surrounding components must be restored. Neither equality
identifies execution traces, transition counts, or resource-limited statuses.
Store canonical request/plan data alongside hashes rather than using hashes
as mathematical certificates.

For example, the B4 words `(1,-1,1,3,3)` and `(3,1,3)` have the same plan:
a rank-two exponent-2 request on original slots {3,4}, followed by a rank-two
exponent-1 request on {1,2}. This uses disconnected support and the proved
rank-two result; it is not higher-rank free reduction.

## Output normal form: an exact formal sum

`canonical_output(output,n=n,complete=True)` produces schema
`zipper-formal-output-normal-form-v1`. `from_certificate(run)` imports supported
current certificates, recording whether the emitted sum is complete or partial.
A certificate marked complete with a nonempty pending frontier is rejected.

Normalization collects coefficients of identical full terms, removes zero
terms, sorts the terms and coefficient powers, and canonicalizes dictionary
key order. It retains:

- the complete routed output word;
- the extraction order and all component words;
- original and local endpoint sets;
- twist powers, orientation, hierarchy, parity, and routing fields;
- additional JSON term/record fields, as part of term identity.

No braid relation, free cancellation, component permutation, or further slider
identity is used. In particular the routed identity word `(1,-1)` is retained
when that is the current output encoding. This is equality of exact formal
outputs, not a decision procedure for equality in the zipper algebra.

Coefficients use canonical decimal strings in the interchange format, avoiding
fixed-width or floating-point loss. `engine_output(normal)` converts those
coefficients back to arbitrary-precision integers while preserving the terms.
`canonical_bytes(normal)` re-normalizes a compatible envelope and gives the deterministic encoding;
`output_fingerprint(normal)` hashes it. Completion and ambient rank are included.

Structural checks cover ranks, endpoint partitions, hierarchy/parity shapes,
orientation consistency, deterministic ambient routes for the recorded ports,
and literal reconstruction of the full word from the contexts. They do not verify every internal algebraic identity or prove that
the output equals the original input. Independently replay the original
certificate to establish that; see the [independent verifier](ZIPPING_NATIVE_VERIFIER.md).
A partial-output fingerprint identifies only
the emitted sum, not its pending frontier or the complete computation; retain
the original checkpoint to resume.

## Execution binding and reconstruction

`execution_identity(run)` extracts the source/catalog/execution identity.
`execution_cache_key(request,identity)` binds a primitive request to that
identity. An output fingerprint is deliberately mathematical and independent
of source provenance; it must not be substituted for an execution cache key.

`reconstruct_output(form,leaf_runs,source_identity=identity)` restores the full
ambient canonical output. It requires completed leaf certificates bound to
the exact requested rank and word and the explicitly supplied execution
identity. It rebuilds local endpoint sets, stable-partition routes, ordered
contexts, and the full word; multiplies exact coefficients across independent
blocks; and collects the resulting identical full terms.

The result is a reconstructed formal sum, not a fabricated execution trace or
run certificate. Store the checked input witness and verified leaf certificates
as its proof dependencies. Unknown term/record extensions are preserved by the
output encoder, but rejected by reconstruction when no embedding rule for those
fields is defined. Incomplete or mixed-version leaf runs are not accepted.

```python
from zipper_algebra.zipping_algorithm.collected import CollectedSession
from zipper_algebra.zipping_algorithm.input_normal_form import (
    input_normal_form, leaf_requests, reconstruction_key,
)
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import (
    execution_identity, execution_cache_key, reconstruct_output,
)
from zipper_algebra.zipping_algorithm.output_normal_form import output_fingerprint

form = input_normal_form((1, -1, 1, 3, 3), n=4)
requests = leaf_requests(form)
session = CollectedSession()
runs = {r["key"]: session.reduce(r["word"], n=r["n"]) for r in requests}
identity = execution_identity(next(iter(runs.values())))
# Independently replay imported/stored certificates before trusting them.
cache_keys = [execution_cache_key(r, identity) for r in requests]
normal = reconstruct_output(form, runs, source_identity=identity)
print(reconstruction_key(form), output_fingerprint(normal))
```

For an identity-only plan there are no leaf runs; pass the execution identity
from the trusted result store and an empty mapping. Existing certificate and
checkpoint formats remain the execution/replay protocols.

## Current baseline and proof status

Through ranks 2–8 and lengths at most 10, the support/peeling baseline was
42,454,830,949 residual literal requests. Rank-two normalization replaces
2,036 of them by 19 exponent representatives: exponents 1,...,10 and -2,...,-10.
Exponent zero is identity and exponent -1 is handled by peeling. The new
sufficient current-procedure domain is **42,454,828,932 requests**, plus one
identity. This is a small reduction; savings from repeated rank-two complements
inside larger runs have not been measured.

All 2,047 rank-two words through length ten matched full outputs by exponent
and passed independent C++ replay; 24 also passed cold Python replay. These
checks support the general rank-two proof. General-rank inverse cancellation
and unique-index far-swap invariance remain open. No modified Artin or reduced-
trace preprocessing is enabled by this interface.

The [inverse-cancellation audit](zipper-inverse-cancellation.md) records the
remaining proof obligations, additional bounded evidence, and conditional
counts. Free reduction at arbitrary rank is not enabled.
