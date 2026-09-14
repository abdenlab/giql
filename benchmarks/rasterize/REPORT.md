# Cherimoya ATAC pre-processing: standard pipeline versus GIQL on DataFusion

Wall-clock comparison of cherimoya's ATAC-seq pre-processing, from BAM to the tensors `Cherimoya.fit` moves to the GPU, run on real deeply sequenced data. A third arm measures what the same path costs when the pileup is expressed with `DISJOIN` instead of hand-written SQL, as evidence for [giql#246](https://github.com/abdenlab/giql/issues/246).

## Result

**Replacing the pileup and locus extraction with GIQL on DataFusion makes the whole pre-GPU path 2.4 times faster.** The margin is stable across two orders of magnitude of input, and the tensors are bit-identical.

| scale | cut sites | standard pipeline | GIQL | speedup |
|---|---|---|---|---|
| chr21 | 1.28 M | 67.6 s | 30.8 s | 2.19x |
| chr8 | 4.12 M | 72.1 s | 32.7 s | 2.20x |
| chr1 | 12.09 M | 82.9 s | 37.1 s | 2.24x |
| genome | 99.39 M | 183.6 s | 76.5 s | **2.40x** |

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

Data is ENCODE [ENCSR483RKN](https://www.encodeproject.org/experiments/ENCSR483RKN/), K562 ATAC-seq: filtered alignments `ENCFF512VEZ` and `ENCFF987XOV` totalling 99.4 M cut sites, IDR peaks `ENCFF925CYR`, exclusion list `ENCFF356LFX`, UCSC hg38. Tn5 shift `+4 / -4`, unstranded, per `docs/recipes/atacseq.rst`. Each arm ran as its own subprocess under a 3600 s cap; none timed out.

**cherimoya is used unmodified.** The benchmark depends on upstream `jmschrei/cherimoya` and needs no patches to it. Swapping the loci extractor would normally require one, since `PeakGenerator` hard-codes tangermeme's, so `peak_generator` in the script reproduces its body against the same public pieces and takes the extractor as an argument. `PeakNegativeSampler` is reused unchanged, which is what keeps the two arms comparable, and the `wrapper` subcommand pins the equivalence: over the chr21 reference bigWig the wrapper and cherimoya's own `PeakGenerator` produce the same sha256 across all 27 batches.

```bash
uv run benchmarks/rasterize/pileup_bench.py run --data DIR --explain
uv run benchmarks/rasterize/pileup_bench.py gate --data DIR
```

## Correctness

Timings mean nothing if the arms disagree, so each rung is gated before it is reported. At every scale all three arms produced **bit-identical** sequence, signal and mask tensors, and identical sha256 digests over every batch at the GPU-ingestion boundary.

| scale | digest | peaks kept | batches |
|---|---|---|---|
| chr21 | `523b0d9ed24a` | 1 400 | 27 |
| chr8 | `3c8f906063ee` | 4 796 | 93 |
| chr1 | `184b1cb72d26` | 14 134 | 276 |
| genome | `65c5cb2098ec` | 85 065 | 1 661 |

## Where the time goes

Genome scale, seconds. `cut sites` is the BAM scan and Tn5 transform; the standard arm has no equivalent line because `bam2bw` fuses scanning and aggregation.

| stage | standard | GIQL |
|---|---|---|
| cut sites | - | 18.8 |
| pileup | 152.3 | 13.9 |
| loader (both extractions) | 26.3 | 39.0 |
| epoch through the sink | 5.0 | 4.9 |
| **total** | **183.6** | **76.5** |
| peak RSS | 6.8 GB | 18.3 GB |

The win is concentrated entirely in one stage. BAM to pileup takes 152.3 s in the standard path against 32.7 s in GIQL, a 4.7-fold difference, and that single stage alone is more than the entire GIQL path.

Two honest counterpoints. GIQL's locus extraction is **slower**, 39.0 s against 26.3 s, so the INTERSECTS join, the FASTA provider and per-window one-hot encoding give back about a third of what the pileup gains; that is where the next optimization belongs, not in the pileup. And GIQL uses 2.7 times the memory, because it holds the pileup and the extracted tensors in RAM where the standard path streams a bigWig from disk.

A shared GC-matched negatives step, identical in both arms, runs once per rung and is charged to neither total: 28.8 s at genome scale.

## The DISJOIN arm, and giql#246

Coverage is already expressible in GIQL, since `DISJOIN` plus `GROUP BY` reproduces a pileup over 1 bp cut sites exactly and true per-base depth over wide intervals. The question this benchmark was built to answer is whether that spelling is fast enough to make a dedicated `RASTERIZE` operator redundant. It is not.

| genome scale | pileup | total path |
|---|---|---|
| GIQL, hand-written aggregate | 13.9 s | 76.5 s |
| GIQL, `DISJOIN` spelling | 263.1 s | 318.0 s |
| standard pipeline | 152.3 s | 183.6 s |

Expressing the pileup with `DISJOIN` costs about 19-fold on that stage and 2.6 times the memory, which is enough to turn a 2.4x speedup into a 1.73x slowdown. `EXPLAIN` confirms polars-bio genuinely rewrites `DISJOIN`'s breakpoint join into `IntervalJoinExec`, so this is not a quadratic fallback but the cost of the surrounding CTE chain: a deduplicating `UNION` over twice the input, an interval join that finds nothing on point input, a `LEAD` window, and a three-key hash join back to the targets, standing in for one aggregate.

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
- **One defect found in passing, and it is polars-bio's, not GIQL's.** `DISJOIN` emits a strict range comparison for its breakpoint join, which is exactly the shape polars-bio rewrites into `IntervalJoinExec`. That rewrite converts strict to non-strict by adjusting by one, and builds the adjustment as an uncoerced Int32 literal, so it raises `Invalid arithmetic operation: Int64 - Int32` against Int64 coordinates. The benchmark works around it by narrowing coordinates to Int32 for that plan alone; hg38 fits comfortably. Nothing else needs it, since the hash aggregate never reaches the rule and the `INTERSECTS` extraction is safe because the datafusion-bio target already emits the closed non-strict form. Were the coercion fixed upstream, the cut-site query would need no casts at all.
