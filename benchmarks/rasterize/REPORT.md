# Expressing a pileup with DISJOIN: what it costs

Benchmark evidence for [#246](https://github.com/abdenlab/giql/issues/246) (add a `RASTERIZE` operator).

## Question

Coverage is already expressible in GIQL today. Over 1 bp cut sites, `DISJOIN` plus `GROUP BY` reproduces a pileup exactly; over wide intervals it reproduces true per-base depth as RLE runs. Both identities were verified against DuckDB before benchmarking. So the question is not whether `RASTERIZE` adds expressiveness, but whether the existing spelling is fast enough that a dedicated operator would be redundant.

To answer it on real data, three arms run cherimoya's whole ATAC pre-processing path, from BAM to the tensors `Cherimoya.fit` moves to the GPU:

| arm | pileup | locus extraction |
|---|---|---|
| `reference` | `bam2bw` to bigWig | `tangermeme.io.extract_loci` |
| `giql-disjoin` | GIQL `DISJOIN` + `GROUP BY`, in memory | GIQL `INTERSECTS` |
| `giql-groupby` | hand-written `GROUP BY chrom, pos`, in memory | GIQL `INTERSECTS` |

The third arm is the plan a `RASTERIZE` sweep-line fast path would be expected to produce, so it turns "`DISJOIN` is slow" into a measured gap.

## Answer

**`DISJOIN` costs a flat 20-fold penalty on the pileup stage and 2.5 times the memory.** It does not blow up asymptotically, because polars-bio's rule does fire on the breakpoint join, but the constant factor is large enough to invert the end-to-end result: the GIQL path is **2.4 times faster** than the standard pipeline when the pileup is written by hand, and **1.75 times slower** when it is expressed with `DISJOIN`.

That is the case for the fast path in #246. The operator is not needed for expressiveness. It is needed so that the expressible spelling is also the fast one.

## Setup

Run on `atlas`, AMD EPYC 7763, 128 cores, 503 GB RAM, Linux 6.8, no GPU, bare metal at load average 0.2 (not through SLURM, matching how the existing `~/giql-bench` harness was run). Python 3.12.14 with polars-bio 0.35.1, polars 1.44.2, datafusion 53.0.0, pyarrow 24.0.0, torch 2.14.0, tangermeme 1.4.1, giql 0.7+12e772d, cherimoya 0.2.1, bam2bw 0.5.1.

Data is ENCODE [ENCSR483RKN](https://www.encodeproject.org/experiments/ENCSR483RKN/), deeply sequenced K562 ATAC-seq: filtered alignments `ENCFF512VEZ` and `ENCFF987XOV` (99.4 M cut sites), IDR peaks `ENCFF925CYR`, exclusion list `ENCFF356LFX`, UCSC hg38. Tn5 shift `+4 / -4`, unstranded, matching the recipe in `docs/recipes/atacseq.rst`. Every arm ran as its own subprocess under a 3600 s cap; none timed out. Whole ladder: 15 min 57 s wall, 47.9 GB peak RSS.

Reproduce with:

```bash
uv run benchmarks/rasterize/pileup_bench.py run --data DIR --explain
uv run benchmarks/rasterize/pileup_bench.py gate --data DIR
```

## Correctness gate

Timings are only meaningful if the arms agree, so each rung is gated before it is reported. At every scale all three arms produced **bit-identical** sequence, signal and mask tensors, and identical sha256 digests over every batch at the GPU-ingestion boundary:

| scale | digest | peaks kept | batches |
|---|---|---|---|
| chr21 | `40e3e3b85459` | 1 400 | 27 |
| chr8 | `861d1aac1c0c` | 4 796 | 93 |
| chr1 | `e897938b1e36` | 14 134 | 276 |
| genome | `afa0f272f665` | 85 065 | 1 661 |

## Pileup stage

The stage where the arms actually differ. Seconds.

| scale | cut sites | positions | `DISJOIN` | `GROUP BY` | penalty |
|---|---|---|---|---|---|
| chr21 | 1.28 M | 0.87 M | 2.28 | 0.10 | 23× |
| chr8 | 4.12 M | 3.03 M | 7.99 | 0.31 | 26× |
| chr1 | 12.09 M | 7.65 M | 20.97 | 0.95 | 22× |
| genome | 99.39 M | 68.03 M | 264.87 | 13.73 | 19× |

The penalty is roughly constant across two orders of magnitude of input. Both spellings scale close to linearly from chr1 to genome: 8.2 times the data costs `DISJOIN` 12.6 times the time and `GROUP BY` 14.4 times.

## Whole path

Seconds, BAM to the GPU boundary. `cut sites` is the BAM scan and Tn5 transform, identical across the GIQL arms; the reference arm has no equivalent because `bam2bw` fuses scanning and aggregation.

| scale | arm | cut sites | pileup | loader | epoch | total | peak RSS |
|---|---|---|---|---|---|---|---|
| chr21 | `reference` | - | 67.3 | 0.4 | 0.1 | **67.9** | 0.9 GB |
| chr21 | `giql-disjoin` | 17.5 | 2.3 | 16.8 | 0.1 | **36.6** | 10.7 GB |
| chr21 | `giql-groupby` | 17.4 | 0.1 | 11.3 | 0.1 | **28.9** | 10.5 GB |
| chr8 | `reference` | - | 70.5 | 1.4 | 0.3 | **72.1** | 1.1 GB |
| chr8 | `giql-disjoin` | 17.5 | 8.0 | 14.3 | 0.3 | **40.1** | 11.3 GB |
| chr8 | `giql-groupby` | 17.5 | 0.3 | 10.8 | 0.3 | **28.9** | 10.8 GB |
| chr1 | `reference` | - | 77.7 | 3.7 | 0.8 | **82.2** | 1.8 GB |
| chr1 | `giql-disjoin` | 17.7 | 21.0 | 16.7 | 0.8 | **56.2** | 12.8 GB |
| chr1 | `giql-groupby` | 17.6 | 1.0 | 17.0 | 0.9 | **36.5** | 11.5 GB |
| genome | `reference` | - | 151.0 | 26.2 | 5.0 | **182.2** | 7.1 GB |
| genome | `giql-disjoin` | 19.0 | 264.9 | 29.9 | 4.8 | **318.5** | 46.9 GB |
| genome | `giql-groupby` | 19.1 | 13.7 | 38.3 | 5.4 | **76.5** | 18.8 GB |

The crossover is the point of the table. At every sub-genome rung `DISJOIN` still beats the reference pipeline, because `bam2bw` pays a near-fixed 67 to 78 s scanning both BAMs in pure Python whatever the rung. At genome scale the pileup stage dominates and the ordering flips.

## Why it is a constant factor and not a blowup

`DISJOIN` expands to a five-CTE chain whose third `UNION` branch joins targets against breakpoints on one equality key plus two strict range comparisons:

```sql
SELECT ... FROM __giql_dj_tgt AS t
JOIN __giql_dj_bp AS bp ON bp.chrom = t."chrom" AND bp.pos > t."start" AND bp.pos < t."end"
```

On 1 bp cut sites that branch returns nothing, but it is still planned and executed. `EXPLAIN` confirms polars-bio rewrites it rather than falling back:

```
IntervalJoinExec: mode=CollectLeft, join_type=Inner, on=[(chrom@0, chrom@0)],
  filter=pos@2 > start@0 AND pos@2 < end@1
HashJoinExec: mode=CollectLeft, join_type=Inner,
  on=[(chrom@0, kc@0), (start@1, ks@1), (end@2, ke@2)]
WindowAggExec
```

So the cost is not a quadratic fallback. It is the sum of a `UNION` that deduplicates twice the input, an interval join that finds nothing, a `LEAD` window over the union, and a three-key hash join back to the targets, in place of one hash aggregate. Memory follows: 46.9 GB against 18.8 GB at genome scale.

This matters for the #246 design. The plan's proposed **sweep-line fast path** (`+1` / `-1` deltas, pre-aggregated by position, cumulative `SUM` with `LEAD` for run ends) targets exactly the work this measurement shows is wasted, and the `giql-groupby` column is a lower bound on what it should achieve for the self-grid, invertible-aggregate case.

## Caveats

- **Thread counts differ and are not normalized.** polars-bio pins `datafusion.execution.target_partitions=1`, so both GIQL arms are effectively single-threaded per query, while `bam2bw` was given `-p -1` and forks one process per BAM. The GIQL arms therefore win the end-to-end comparison from behind.
- **Sub-genome rungs do not reduce IO.** Both pipelines scan both BAMs in full at every rung; only the aggregation shrinks. That is deliberate, so the ladder isolates aggregation scaling, but it means the sub-genome totals are dominated by a fixed scan.
- **The GIQL arms hold the pileup in memory,** which is why their RSS is an order of magnitude above the reference at small rungs. `--materialize` writes Parquet instead.
- **One defect found while building this.** polars-bio's interval-join rule raises `Invalid arithmetic operation: Int64 - Int32` when the strict comparison mixes widths, so the cut-site view must cast both bounds to Int32 *after* their arithmetic. Casting first and adding 1 afterwards silently promotes `end` back to Int64 and reintroduces the crash. Any `RASTERIZE` expander that emits strict comparisons will meet the same defect.
