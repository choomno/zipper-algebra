# Experimental transition DAG certificates

`experiments.zipper_derivation_dag` captures a bounded graph of complete
continuation states, records the transition from each state, and checks every
reachable node using the independent native verifier. It is an experiment,
not a replacement for production census certification.

The experimental wire schema is `zipper-experimental-transition-dag-v2`;
the earlier exploratory shape is not accepted as a v2 certificate.

The checker binds the expected literal input, source/catalog identities and
per-transition limits. Implementation identity retains authenticated source/build
fields and excludes session timings, cache-hit build flags, and executable paths.
Node identity includes every continuation field except diagnostic depth. Coefficients and outer contexts remain in the state, so this experiment
uses less sharing than the production intrinsic-continuation evaluator.
Children must reference present nodes, every node must be reachable from the
input root, and the graph must be acyclic. Every recorded transition, including
its evidence, must agree with the independent native transition implementation.
Exact terminal output is assembled separately using integer path
multiplicities and Laurent coefficients, then compared with the claimed output.
Node, accumulated terminal-key, and integer-bit guards bound the experiment;
resource refusal supplies no completed certificate. These are separate graph
assembly guards: terminal-key and coefficient-bit limits apply during assembly,
including contributions that may later cancel. The experiment does not certify
the standalone scheduler's step count, frontier limits, depth-overflow behavior,
or success under its whole-run budget. The certificate records this scope as
`transition_limits` and an explicit semantics record.

The root uses the versioned `literal-singleton-v1` construction in `initial()`;
it is not obtained through a native root-construction API. Source inspection
shows that the current producer and checker carry and increment depth rather
than using it to select geometric rewrites. Normalizing depth deliberately
omits diagnostic depth growth and overflow. The benchmark compares every graph
answer with a full independently replayed native run, checking root/output
agreement on those inputs. Those finite checks do not replace an arbitrary-rank
proof of the initialization and depth-normalization conventions.

The existing native transition diagnostic still discovers the applicable rule
and performs the whole transition. Consequently this prototype does **not**
implement a smaller local-rule proof checker. It measures the graph transport,
storage, sharing and assembly boundary available through the current API.
One source audit surrounds the entire frozen-source graph check; the checker
does not rehash every source file separately for every node.

Run a bounded experiment with:

```sh
python -m experiments.zipper_derivation_dag /path/to/results --max-nodes 20000 \
  --cases '[[3,[1,2,1]],[3,[-1,-2,-1,-2]],[4,[1,3,2]]]'
```

The output directory contains compressed exact graph certificates and timing
records for ordinary native production and independent replay, graph capture,
storage/loading, and graph verification. The benchmark separately compares the
graph's terminal sum with the independently replayed native result. Process
initialization and compilation are excluded; graph capture and checking include
per-transition process transfers. These small experiments are not population
throughput estimates.

Initial measurements on frozen revision `357b55e` showed that streamed JSON
transition graphs were slower than native whole-run replay. On a negative
three-strand length-ten input, the graph had 13,640 nodes and 84 terminal terms,
and required approximately 14.3 MB of JSON (1.7 MB compressed). Graph creation,
storage and checking exceeded the cost of ordinary production plus replay.
This route therefore remains experimental. A useful replacement would need
native batched checking and a certificate interface that validates an indicated
local rule without repeating geometric rule discovery. Those changes require
separate correctness and performance evidence before production adoption.
