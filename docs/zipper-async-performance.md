# Asynchronous census measurements

Separate computation and certification improve control of the census backlog;
they do not automatically reduce total mathematical work. On the measured
low-rank cohort, one producer with two certifiers performs about as well as two
synchronous producer/checker pairs, with lower peak sampled memory. One producer
and one certifier leaves the certification stage underprovisioned.

The maintained driver is `experiments.zipper_async_benchmark`. It exercises the
production queue APIs, exact native engine and independent native verifier. Each
arm uses fresh sessions and a new SQLite sidecar. A barrier separates native
initialization from measured work. Work wall time includes transfers, exact
candidate normalization, compression, committed storage, certification, queue
coordination and worker teardown. Worker work time and setup are also recorded.
Total CPU includes setup and must not be confused with work-only CPU time.

All workers inherited CPU affinity 10,12 on the development host; no CPU was
reserved exclusively. The existing census and quality checks continued on the
same host. Cache budgets were 256 MiB per native cache, and both routes used the
same authenticated tcmalloc configuration. Every comparison required identical
inputs, successful independent certification and equality of complete canonical
outputs. Output and certificate digests were checked, but digest equality alone
was not the output comparison. Source manifests and exact policies accompany
each SQLite measurement.

## Hard inputs

The cohort consists of B5 `[-4,-3,1,2,3,4]` and B8 `[4,1,6,7,2,5,6,3]`.
Arms ran in order sequential, split, paired2, paired2, split, sequential.

| Scheduling | Two elapsed times | Median | Maximum sampled tree RSS |
|---|---:|---:|---:|
| One synchronous producer/checker pair | 47.805, 47.436 s | 47.621 s | 3.77 GB |
| One producer and one certifier | 42.168, 38.694 s | 40.431 s | 2.44 GB |
| Two synchronous producer/checker pairs | 29.696, 29.867 s | 29.781 s | 4.07 GB |

Both inputs completed and were independently certified in every arm; all 15
pairwise comparisons of their full outputs agreed. The two compressed candidate
certificates occupied approximately 412 KB. First-arm compilation took about
104 seconds and is excluded from the work-wall column. Subsequent initialization
took 0.31–0.45 seconds. Peak RSS includes initialization and can double-count
shared pages; it is not unique resident memory or a cache-budget measurement.

## Longer low-rank inputs

This cohort contains 16 distinct primitive inputs: four each in B3 at lengths
8, 10 and 12, and B4 at length 6. For a layer with N primitive inputs, its four
ordinals are floor((2j+1)N/8), for j=0,1,2,3. This deterministic diagnostic sample
is not a random population sample or a census-wide runtime estimate.

Arms ran in order paired2, split, split12, split12, split, paired2.

| Scheduling | Two elapsed times | Median | Maximum sampled tree RSS |
|---|---:|---:|---:|
| Two synchronous producer/checker pairs | 68.625, 72.199 s | 70.412 s | 2.66 GB |
| One producer and one certifier | 108.904, 102.624 s | 105.764 s | 1.16 GB |
| One producer and two certifiers | 72.363, 72.714 s | 72.538 s | 2.07 GB |

All 16 inputs completed and were independently certified in every arm; all 15
pairwise full-output comparisons agreed. Compressed candidates occupied about
75 KB per arm. Initialization took 0.37–0.51 seconds. Total CPU including setup
was 135.6/141.5 seconds for paired2 and 136.0/137.9 seconds for split12. The
separation of roles did not reduce the measured CPU work.

The split12 arm has three role processes and three native sessions, while
paired2 has two role processes and four native sessions. Both inherit the same
two-CPU affinity. Thus these are comparisons under a common CPU affinity,
not equal process counts or equal aggregate cache capacity. At approximately
3% slower elapsed time and 22% lower peak sampled RSS than paired2, split12 is
a reasonable configuration for the durable staged queue on this cohort. Two
observations do not establish a population confidence interval or a universal
best allocation.

The hard run preceded final queue robustness fixes; the low-rank run included
those queue fixes. Neither changed the mathematical backend. A final-source
three-request smoke check used the shared production verification adapter and
confirmed exact certified equality between sequential and split12 execution.
This smoke check is functional validation, not a throughput measurement.

See [the queue contract](zipper-async-census.md) for provisional versus certified
coverage and [the census plan](zipper-census-plan.md) for the chosen workload.
