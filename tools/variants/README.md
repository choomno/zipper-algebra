# Experimental variant-quotient calculations

These scripts explore dimensions, trace forms, radical candidates, tower maps and
Markov-trace conditions for specified untwisting quotients. They are research
programs retained with their earlier inputs and outputs; they are outside the
maintained installed package.

| Script | Computation |
| --- | --- |
| `compute_dim_aho.py` | Standard words and dimension from a supplied basis. |
| `check_x3_nilpotent.py` | Specialized nilpotence search. |
| `semisimple_Q.py`, `semisimple_v3.py` | Trace-form calculations over the rationals or a finite field. |
| `trace_space.py`, `trace_space_fast.py` | Trace-space linear algebra in their stated coefficient domain. |
| `tower_check.py` | Tests of proposed maps at chosen parameters. |
| `markov_search.py` | Linear systems for candidate Markov traces. |
| `square_zero_*`, `explicit_rep_n3.py`, `generic_non_ss.py` | Specialized representation and ideal calculations. |

A dimension or normal-form conclusion requires a completed basis of the intended
ideal, with the coefficient field and specialization stated. Raising an engine
iteration limit does not prove completion. A calculation at several parameter
values does not establish the corresponding generic statement. In particular,
rank agreement modulo primes and generic non-semisimplicity require separate
justification.

The historical [claim ledger](../../LEDGER_claims.md) and
[specialization audit](../../AUDIT_spec_vs_generic.md) record earlier claims,
evidence and gaps. They have not all been reestablished by the current package
tests. Their manuscript references name research sources outside the software
distribution.

The older solver integration suite is `tools/tests.py`. See
[TESTING](../../docs/TESTING.md) for its status and for the maintained CI suite.
