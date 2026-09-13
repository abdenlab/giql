# Cherimoya ATAC pre-processing: standard pipeline versus GIQL on DataFusion

Wall-clock comparison of cherimoya's ATAC-seq pre-processing, from BAM to the tensors `Cherimoya.fit` moves to the GPU, run on real deeply sequenced data. A third arm measures what the same path costs when the pileup is expressed with `DISJOIN` instead of hand-written SQL, as evidence for [giql#246](https://github.com/abdenlab/giql/issues/246).

## Result

**Replacing the pileup and locus extraction with GIQL on DataFusion makes the whole pre-GPU path 2.4 times faster.** The margin is stable across two orders of magnitude of input, and the tensors are bit-identical.

| scale | cut sites | standard pipeline | GIQL | speedup |
|---|---|---|---|---|
| chr21 | 1.28 M | 67.9 s | 28.9 s | 2.35x |
| chr8 | 4.12 M | 72.1 s | 28.9 s | 2.50x |
| chr1 | 12.09 M | 82.2 s | 36.5 s | 2.25x |
| genome | 99.39 M | 182.2 s | 76.5 s | **2.38x** |

The GIQL path also writes nothing to disk. The standard path's pileup exists only as a bigWig on the way to being read back one locus at a time.

## What the two paths do

| stage | standard | GIQL |
|---|---|---|
| BAM to pileup | `bam2bw`, a single-threaded Python loop over every record into per-chromosome dicts, then a bigWig | a DataFusion view over the BAM table provider deriving the Tn5 cut site, then one hash aggregate, held in memory |
| loci to tensors | `tangermeme.io.extract_loci`, per-locus pybigtools reads plus pyfaidx | GIQL `INTERSECTS` joins for the exclusion list and the per-locus signal, sequences from the FASTA provider |
| batching | `PeakGenerator` into a DataLoader | unchanged, through an injected extractor |

Only the first two stages differ. Sampling, augmentation and batching are the same code in both arms, which is what makes the tensor comparison meaningful.

## Setup

`atlas`, AMD EPYC 7763, 128 cores, 503 GB RAM, Linux 6.8, no GPU, bare metal at load average 0.2. Python 3.12.14, polars-bio 0.35.1, polars 1.44.2, datafusion 53.0.0, pyarrow 24.0.0, torch 2.14.0, tangermeme 1.4.1, giql 0.7+12e772d, cherimoya 0.2.1, bam2bw 0.5.1.

Data is ENCODE [ENCSR483RKN](https://www.encodeproject.org/experiments/ENCSR483RKN/), K562 ATAC-seq: filtered alignments `ENCFF512VEZ` and `ENCFF987XOV` totalling 99.4 M cut sites, IDR peaks `ENCFF925CYR`, exclusion list `ENCFF356LFX`, UCSC hg38. Tn5 shift `+4 / -4`, unstranded, per `docs/recipes/atacseq.rst`. Each arm ran as its own subprocess under a 3600 s cap; none timed out. Whole ladder: 15 min 57 s.

```bash
uv run benchmarks/rasterize/pileup_bench.py run --data DIR --explain
uv run benchmarks/rasterize/pileup_bench.py gate --data DIR
```

## Correctness

Timings mean nothing if the arms disagree, so each rung is gated before it is reported. At every scale all three arms produced **bit-identical** sequence, signal and mask tensors, and identical sha256 digests over every batch at the GPU-ingestion boundary.

| scale | digest | peaks kept | batches |
|---|---|---|---|
| chr21 | `40e3e3b85459` | 1 400 | 27 |
| chr8 | `861d1aac1c0c` | 4 796 | 93 |
| chr1 | `e897938b1e36` | 14 134 | 276 |
| genome | `afa0f272f665` | 85 065 | 1 661 |

## Where the time goes

Genome scale, seconds. `cut sites` is the BAM scan and Tn5 transform; the standard arm has no equivalent line because `bam2bw` fuses scanning and aggregation.

| stage | standard | GIQL |
|---|---|---|
| cut sites | - | 19.1 |
| pileup | 151.0 | 13.7 |
| loader (both extractions) | 26.2 | 38.3 |
| epoch through the sink | 5.0 | 5.4 |
| **total** | **182.2** | **76.5** |
| peak RSS | 7.1 GB | 18.8 GB |

The win is concentrated entirely in one stage. BAM to pileup takes 151.0 s in the standard path against 32.8 s in GIQL, a 4.6-fold difference, and that single stage is more than the entire GIQL path.

Two honest counterpoints. GIQL's locus extraction is **slower**, 38.3 s against 26.2 s, so the INTERSECTS join and the FASTA provider give back about a third of what the pileup gains; that is where the next optimization belongs, not in the pileup. And GIQL uses 2.6 times the memory, because it holds the pileup and the extracted tensors in RAM where the standard path streams a bigWig from disk.

A shared GC-matched negatives step, identical in both arms, runs once per rung and is charged to neither total: 28.9 s at genome scale.

## The DISJOIN arm, and giql#246

Coverage is already expressible in GIQL, since `DISJOIN` plus `GROUP BY` reproduces a pileup over 1 bp cut sites exactly and true per-base depth over wide intervals. The question this benchmark was built to answer is whether that spelling is fast enough to make a dedicated `RASTERIZE` operator redundant. It is not.

| genome scale | pileup | total path |
|---|---|---|
| GIQL, hand-written aggregate | 13.7 s | 76.5 s |
| GIQL, `DISJOIN` spelling | 264.9 s | 318.5 s |
| standard pipeline | 151.0 s | 182.2 s |

Expressing the pileup with `DISJOIN` costs about 20-fold on that stage and 2.5 times the memory, which is enough to turn a 2.4x speedup into a 1.75x slowdown. `EXPLAIN` confirms polars-bio genuinely rewrites `DISJOIN`'s breakpoint join into `IntervalJoinExec`, so this is not a quadratic fallback but the cost of the surrounding CTE chain: a deduplicating `UNION` over twice the input, an interval join that finds nothing on point input, a `LEAD` window, and a three-key hash join back to the targets, standing in for one aggregate.

Plan comparison on identical input, 99.4 M intervals:

| input | hash aggregate | sweep-line | `DISJOIN` |
|---|---|---|---|
| 1 bp points | **13.2 s** | 155.5 s | 264.9 s |
| 92.5 bp reads | not a valid plan | **143.9 s** | 7x the point cost (at chr21) |

So a self-grid `RASTERIZE` needs more than one fast path. Point input wants the hash aggregate, wide input wants the sweep-line, which is flat in width where `DISJOIN` degrades sevenfold with it. Since GIQL cannot infer interval width at transpile time, the choice has to be declared. `pileup_sql` in the benchmark script carries the full argument.

## Caveats

- **Thread counts are not normalized.** polars-bio pins `datafusion.execution.target_partitions=1`, so the GIQL arms are effectively single-threaded per query, while `bam2bw` was given `-p -1` and forked one process per BAM. GIQL wins the comparison from behind.
- **Sub-genome rungs do not reduce IO.** Both pipelines scan both BAMs in full at every rung; only the aggregation shrinks. That is deliberate, so the ladder isolates aggregation scaling, but it means the sub-genome totals are dominated by a fixed scan.
- **The GIQL arms hold the pileup in memory,** which is most of the RSS difference. `--materialize` writes Parquet instead, which only helps reruns.
- **One defect found in passing.** polars-bio's interval-join rule raises `Invalid arithmetic operation: Int64 - Int32` when a strict comparison mixes widths, so the cut-site view casts both bounds to Int32 *after* their arithmetic. Casting first and adding 1 afterwards silently promotes `end` back to Int64 and reintroduces the crash.
