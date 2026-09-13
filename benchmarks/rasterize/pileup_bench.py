# /// script
# requires-python = ">=3.11,<3.15"
# dependencies = [
#   "polars-bio>=0.35.1",
#   "polars",
#   "numpy",
#   "pandas",
#   "torch",
#   "tangermeme",
#   "pysam",
#   "pyBigWig",
#   "pyfaidx",
#   "pybigtools",
#   "bam2bw",
#   "giql @ git+https://github.com/abdenlab/giql@cherimoya-mvp",
#   "cherimoya @ git+https://github.com/conradbzura/cherimoya@giql-datafusion-pipeline",
# ]
# ///

"""Benchmark a DISJOIN-expressed pileup against cherimoya's standard pre-GPU path.

Evidence for GIQL issue #246 (add a RASTERIZE operator). Coverage is expressible
today as ``DISJOIN`` plus ``GROUP BY`` -- verified identical to a plain pileup on
1 bp cut sites, and to true per-base depth on wide intervals -- but ``DISJOIN``
expands to a CTE chain whose breakpoint join is one equality key plus two strict
range comparisons. This measures what that costs on real, deeply sequenced ATAC
data against the pipeline it would replace.

Three arms run the whole path from BAM to the tensors ``Cherimoya.fit`` would move
to the GPU, all executing on polars-bio / DataFusion:

    reference     bam2bw -> bigWig, then tangermeme.io.extract_loci
    giql-disjoin  GIQL DISJOIN + GROUP BY pileup, then GIQL INTERSECTS extraction
    giql-groupby  hand-written GROUP BY pileup, then the same GIQL extraction

The third arm is the pileup a RASTERIZE fast path would be expected to match, so
it turns "DISJOIN is slow" into a quantified gap. See `pileup_sql` for what that
fast implementation is, why it is legal, and what it would require of a
RASTERIZE operator; see REPORT.md for the measured comparison.

Every (scale, arm) pair runs as its own subprocess so a blowup is capped and
recorded as ``DNF`` rather than hanging the ladder. The GIQL arms hold the pileup
in memory; ``--materialize`` writes it to Parquet instead, which only helps reruns.

Usage::

    uv run benchmarks/rasterize/pileup_bench.py run --data DIR
    uv run benchmarks/rasterize/pileup_bench.py run --data DIR --scales chr21 --explain
    uv run benchmarks/rasterize/pileup_bench.py report --data DIR
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

os.environ.setdefault("TQDM_DISABLE", "1")
os.environ.setdefault("TORCH_COMPILE_DISABLE", "1")
os.environ.setdefault("TORCHDYNAMO_DISABLE", "1")

IGNORE = list("QWERYUIOPSDFHJKLZXVBNM")
ARMS = ("reference", "giql-disjoin", "giql-groupby")
GIQL_ARMS = ("giql-disjoin", "giql-groupby")

#: Ladder rungs, smallest first. ``None`` means every contig in the FASTA index,
#: which is what the reference pipeline does by default.
SCALES = {
    "chr21": ["chr21"],
    "chr8": ["chr8"],
    "chr1": ["chr1"],
    "genome": None,
}

#: Tn5 shift from the cherimoya ATAC recipe (`bam2bw -ps 4 -ns -4`).
POS_SHIFT, NEG_SHIFT = 4, -4


###
# Paths
###


def add_common_args(parser):
    parser.add_argument(
        "--data",
        default=".",
        help="Directory holding the BAMs, peaks, FASTA and exclusion list.",
    )
    parser.add_argument(
        "--bams", nargs="+", default=["ENCFF512VEZ.bam", "ENCFF987XOV.bam"]
    )
    parser.add_argument("--peaks", default="ENCFF925CYR.bed.gz")
    parser.add_argument("--fasta", default="hg38.fa")
    parser.add_argument("--exclusion", default="ENCFF356LFX.bed.gz")
    parser.add_argument(
        "--results",
        default=None,
        help="Directory for per-run JSON (default <data>/results).",
    )
    parser.add_argument("--in-window", type=int, default=2114)
    parser.add_argument("--out-window", type=int, default=1000)
    parser.add_argument("--max-jitter", type=int, default=500)
    parser.add_argument("--negative-ratio", type=float, default=0.25)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument(
        "--materialize",
        action="store_true",
        help="Write the GIQL pileup to Parquet instead of holding it in memory.",
    )
    parser.add_argument(
        "--explain",
        action="store_true",
        help="Record the physical plan of the pileup query.",
    )
    parser.add_argument(
        "--bam2bw",
        default=None,
        help="Path to the bam2bw executable (default: next to this interpreter, then PATH).",
    )


def resolve(args):
    args.data = os.path.abspath(args.data)
    args.results = os.path.abspath(args.results or os.path.join(args.data, "results"))
    os.makedirs(args.results, exist_ok=True)

    def under(p):
        return p if os.path.isabs(p) else os.path.join(args.data, p)

    args.bams = [under(b) for b in args.bams]
    args.peaks = under(args.peaks)
    args.fasta = under(args.fasta)
    args.exclusion = under(args.exclusion)
    return args


def scale_chroms(scale, fasta):
    """Contigs a scale covers; the full FASTA index for the genome rung."""

    if SCALES[scale] is not None:
        return list(SCALES[scale])

    with open(fasta + ".fai") as f:
        return [line.split("\t")[0] for line in f]


def training_chroms_for(scale):
    """Which contigs the loader draws training loci from.

    Sub-genome rungs train on the rung itself; the genome rung uses cherimoya's
    default split, which is what the verified ENCSR483RKN run used.
    """

    if SCALES[scale] is not None:
        return list(SCALES[scale])

    from cherimoya_cli.defaults import training_chroms

    return list(training_chroms)


def scale_paths(args, scale):
    prefix = os.path.join(args.results, "_" + scale)
    return {
        "sizes": prefix + ".chrom.sizes",
        "peaks": prefix + ".peaks.bed",
        "negatives": prefix + ".negatives.bed",
        "bigwig_prefix": prefix + ".ref",
        "bigwig": prefix + ".ref.bw",
        "parquet": prefix + ".cuts.parquet",
    }


###
# Shared per-scale inputs
###


def prepare_scale(args, scale):
    """Build the rung's chrom-sizes file, its peak subset, and its negatives.

    All three are inputs both arms consume, so this work is timed once and
    charged to both totals rather than being attributed to either engine.
    """

    import pandas

    paths = scale_paths(args, scale)
    chroms = scale_chroms(scale, args.fasta)
    report = {"n_chroms": len(chroms)}

    sizes = {}
    with open(args.fasta + ".fai") as f:
        for line in f:
            name, length = line.split("\t")[:2]
            sizes[name] = int(length)

    keep = set(chroms)
    with open(paths["sizes"], "w") as f:
        for name in chroms:
            f.write("{}\t{}\n".format(name, sizes[name]))

    peaks = pandas.read_csv(
        args.peaks,
        sep="\t",
        header=None,
        usecols=(0, 1, 2),
        names=["chrom", "start", "end"],
    )
    peaks = peaks[peaks["chrom"].isin(keep)]
    peaks.to_csv(paths["peaks"], sep="\t", header=False, index=False)
    report["n_peaks"] = len(peaks)

    t = time.time()
    if not os.path.exists(paths["negatives"]):
        from tangermeme.match import extract_matching_loci

        negatives = extract_matching_loci(
            loci=paths["peaks"],
            fasta=args.fasta,
            gc_bin_width=0.02,
            max_n_perc=0.1,
            bigwig=None,
            signal_beta=0.5,
            in_window=args.in_window,
            out_window=args.out_window,
            chroms=None,
            verbose=False,
            n_jobs=1,
        )
        negatives.to_csv(paths["negatives"], header=False, sep="\t", index=False)
    report["seconds_negatives"] = round(time.time() - t, 2)
    report["n_negatives"] = sum(1 for _ in open(paths["negatives"]))
    return report


###
# The GIQL pileup, both spellings
###


def cut_site_sql(bam_tables, chroms, restrict):
    """Strand-aware Tn5 cut sites, matching `bam2bw -ps 4 -ns -4` exactly.

    Forward reads cut at ``start + pos_shift``, reverse reads at
    ``end - 1 + neg_shift``; unmapped records are dropped, as are contigs the
    rung does not cover. Coordinates are Int32: hg38 fits comfortably, and
    polars-bio's interval-join rule builds an uncoerced Int32 literal that
    raises on Int64 operands, which DISJOIN's strict breakpoint join would
    otherwise hit. Both bounds are cast *after* their arithmetic: casting
    first and adding 1 afterwards silently promotes `end` back to Int64 and
    reintroduces the crash.
    """

    where = "(CAST(b.flags AS BIGINT) & 4) = 0"
    if restrict:
        where += " AND b.chrom IN ({})".format(
            ", ".join("'{}'".format(c) for c in chroms)
        )

    branches = [
        "SELECT CAST(b.chrom AS VARCHAR) AS chrom, "
        "CAST(CASE WHEN (CAST(b.flags AS BIGINT) & 16) = 0 "
        "THEN CAST(b.start AS BIGINT) + {pos} "
        'ELSE CAST(b."end" AS BIGINT) + {neg} - 1 END AS INT) AS start, '
        "CAST(CASE WHEN (CAST(b.flags AS BIGINT) & 16) = 0 "
        "THEN CAST(b.start AS BIGINT) + {pos} + 1 "
        'ELSE CAST(b."end" AS BIGINT) + {neg} END AS INT) AS "end" '
        "FROM {table} AS b WHERE {where}".format(
            pos=POS_SHIFT, neg=NEG_SHIFT, table=table, where=where
        )
        for table in bam_tables
    ]
    return " UNION ALL ".join(branches)


def pileup_sql(arm):
    """The pileup query for an arm, over a registered `cuts` table.

    Both spellings return one row per occupied position, and the benchmark
    gates on them producing bit-identical tensors downstream. They differ by
    roughly 20x, which is the whole point of the exercise.

    THE FAST IMPLEMENTATION
    -----------------------
    `giql-groupby` is a single hash aggregate::

        SELECT chrom, start AS pos, COUNT(*) FROM cuts GROUP BY chrom, start

    That is a *correct pileup only because a Tn5 cut site is one base wide*.
    Coverage over width-1 intervals degenerates to counting occurrences per
    position: no interval logic is required, no sort, no window, one pass.
    Measured genome-wide over 99.4 M cut sites it runs in 13.2 s, against
    155.5 s for a sweep-line and 264.9 s for the DISJOIN spelling below.

    The same query over *wide* intervals is silently WRONG. It counts
    interval starts per position, not depth: on the same input read as
    92.5 bp alignments it returns 66.7 M rows where true coverage has
    119.5 M runs. Width is what makes this plan legal, and nothing in the
    SQL says so.

    WHAT THIS REQUIRES OF A RASTERIZE OPERATOR (giql#246)
    -----------------------------------------------------
    For `RASTERIZE(cuts)` to emit the fast plan automatically it must know
    the input is point-like, and that is exactly what GIQL cannot work out
    for itself. GIQL is a transpiler with no access to table statistics, and
    DataFusion's own statistics would not help either: they carry per-column
    min, max, null and distinct counts, whereas interval width is a
    two-column derived property (`max(end - start)`), and `max(end) -
    min(start)` is merely the chromosome span.

    So the requirement is a *declaration*, not an inference, and because
    declaring it wrongly yields wrong answers rather than slow ones, it
    belongs with the other load-bearing schema claims (`coordinate_system`,
    `interval_type`) rather than being offered as a performance knob. The
    corollary for the operator's expansion is that a self-grid RASTERIZE
    needs at least three plans, not the two originally proposed:

      * width-1 input, invertible aggregate  -> this hash aggregate
      * wide input, invertible aggregate     -> sweep-line, flat in width
                                                (143.9 s on 92.5 bp reads)
      * anything else                        -> the general cells + join plan

    A fourth, a dense per-contig array as `bam2bw` and `mosdepth` use, is
    plausible where chrom sizes are declared and has not been measured.
    """

    from giql import Table
    from giql import transpile

    # The operator-level spelling: correct on any width, and the identity a
    # RASTERIZE self grid is defined by. EXPLAIN shows polars-bio does rewrite
    # its breakpoint join into IntervalJoinExec, so the ~20x penalty is not a
    # quadratic fallback: it is a deduplicating UNION over twice the input, an
    # interval join that finds nothing on point input, a LEAD window and a
    # three-key hash join back to the targets, in place of one aggregate.
    if arm == "giql-disjoin":
        return transpile(
            "SELECT disjoin_chrom AS chrom, disjoin_start AS pos, "
            "COUNT(*) AS value FROM DISJOIN(cuts) "
            "GROUP BY disjoin_chrom, disjoin_start",
            tables=[
                Table(
                    "cuts",
                    chrom_col="chrom",
                    start_col="start",
                    end_col="end",
                    strand_col=None,
                )
            ],
            dialect="datafusion-bio",
        )

    # The fast plan. Legal only for width-1 intervals; see the docstring.
    return (
        "SELECT chrom, start AS pos, CAST(COUNT(*) AS DOUBLE) AS value "
        "FROM cuts GROUP BY chrom, start"
    )


def build_pileup(args, scale, arm, record):
    """BAM to pileup, timed as two stages.

    `cut_sites` is the BAM scan and the Tn5 transform, identical across the
    GIQL arms. `pileup` is the aggregation, and is the only stage where the
    two GIQL arms differ -- which is the comparison issue #246 turns on.
    """

    import polars
    import polars_bio as pb

    chroms = scale_chroms(scale, args.fasta)
    restrict = SCALES[scale] is not None

    tables = []
    for i, bam in enumerate(args.bams):
        name = "reads_{}".format(i)
        pb.register_bam(bam, name)
        tables.append(name)

    t = time.time()
    cuts = pb.sql(cut_site_sql(tables, chroms, restrict)).collect()
    record["seconds_cut_sites"] = round(time.time() - t, 2)
    record["n_cut_sites"] = cuts.height

    pb.from_polars("cuts", cuts)
    sql = pileup_sql(arm)
    record["pileup_sql"] = sql

    if args.explain:
        plan = pb.sql("EXPLAIN " + sql).collect()
        text = "\n".join(
            str(v)
            for col in plan.get_columns()
            for v in col.to_list()
            if isinstance(v, str)
        )
        record["plan_operators"] = sorted(
            {
                op
                for op in (
                    "IntervalJoinExec",
                    "HashJoinExec",
                    "NestedLoopJoinExec",
                    "SortMergeJoinExec",
                    "WindowAggExec",
                )
                if op in text
            }
        )
        record["plan"] = text

    t = time.time()
    pileup = pb.sql(sql).collect()
    record["seconds_pileup"] = round(time.time() - t, 2)
    record["n_positions"] = pileup.height
    record["total_counts"] = float(pileup["value"].sum())

    pileup = pileup.select(
        polars.col("chrom").cast(polars.Utf8),
        polars.col("pos").cast(polars.Int64),
        polars.col("value").cast(polars.Float64),
    )

    if args.materialize:
        path = scale_paths(args, scale)["parquet"]
        pileup.sort("chrom", "pos").write_parquet(path, statistics=True)
        record["materialized"] = path
        return path

    return pileup


def find_bam2bw(override=None):
    """Locate the bam2bw executable.

    It is installed beside the running interpreter in a venv, but `~/.local/bin`
    is not on PATH for a non-interactive ssh command, so look next to
    `sys.executable` before falling back to PATH and to `~/.local/bin`.
    """

    import shutil

    if override:
        return override

    candidates = [
        os.path.join(os.path.dirname(sys.executable), "bam2bw"),
        shutil.which("bam2bw"),
        os.path.expanduser("~/.local/bin/bam2bw"),
    ]
    for candidate in candidates:
        if candidate and os.path.exists(candidate):
            return candidate

    raise FileNotFoundError("bam2bw not found; pass --bam2bw")


def build_bigwig(args, scale, record):
    """The reference pileup: bam2bw, invoked exactly as `cherimoya pipeline` does.

    The rung is enforced through the chrom-sizes file, which is how bam2bw
    already drops contigs it was not given.
    """

    paths = scale_paths(args, scale)
    cmd = [
        find_bam2bw(getattr(args, "bam2bw", None)),
        "-s",
        paths["sizes"],
        "-n",
        paths["bigwig_prefix"],
        "-ps",
        str(POS_SHIFT),
        "-ns",
        str(NEG_SHIFT),
        "-sf",
        "1",
        "-p",
        "-1",
        "-u",
    ] + args.bams

    env = dict(os.environ)
    env.setdefault("JOBLIB_START_METHOD", "fork")
    record["bam2bw_command"] = " ".join(cmd)

    t = time.time()
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    record["seconds_pileup"] = round(time.time() - t, 2)
    if proc.returncode != 0:
        raise RuntimeError("bam2bw failed: {}".format(proc.stderr[-2000:]))

    return paths["bigwig"]


###
# The pre-GPU path
###


def sink_unpack(data, batch_size, device="cpu"):
    """What `Cherimoya.fit` does to a batch before the forward pass.

    Mirrors cherimoya/cherimoya.py: unpack, drop a ragged final batch, move and
    cast. Returns None for a batch the fit loop would skip.
    """

    X, y, labels = data[0], data[-2], data[-1]
    X_ctl = data[1].to(device) if len(data) == 4 else None
    if X.shape[0] != batch_size:
        return None

    return X.to(device).float(), X_ctl, y.to(device), labels


def run_arm(args, scale, arm):
    """Run one arm of one rung, from BAM to the GPU-ingestion boundary."""

    import torch
    from cherimoya.io import PeakGenerator

    paths = scale_paths(args, scale)
    record = {"scale": scale, "arm": arm, "pid": os.getpid()}

    if arm == "reference":
        signals, extractor = [build_bigwig(args, scale, record)], None
    else:
        from cherimoya import giql_io

        signals = [build_pileup(args, scale, arm, record)]
        extractor = giql_io.extract_loci

    t = time.time()
    loader = PeakGenerator(
        peaks=paths["peaks"],
        negatives=paths["negatives"],
        sequences=args.fasta,
        signals=signals,
        chroms=training_chroms_for(scale),
        in_window=args.in_window,
        out_window=args.out_window,
        max_jitter=args.max_jitter,
        negative_ratio=args.negative_ratio,
        reverse_complement=True,
        summits=False,
        exclusion_lists=[args.exclusion],
        random_state=args.random_state,
        batch_size=args.batch_size,
        num_workers=0,
        pin_memory=False,
        verbose=False,
        signal_groups=[1],
        extractor=extractor,
    )
    record["seconds_loader"] = round(time.time() - t, 2)
    record["n_peaks_kept"] = int(loader.dataset.peak_sequences.shape[0])
    record["n_negatives_kept"] = int(loader.dataset.negative_sequences.shape[0])

    digest = hashlib.sha256()
    n_batches = n_skipped = 0
    t = time.time()
    for data in loader:
        unpacked = sink_unpack(data, args.batch_size)
        if unpacked is None:
            n_skipped += 1
            continue

        X, _, y, labels = unpacked
        digest.update(X.numpy().tobytes())
        digest.update(y.numpy().tobytes())
        digest.update(labels.numpy().tobytes())
        n_batches += 1
    record["seconds_epoch"] = round(time.time() - t, 2)
    record["n_batches"] = n_batches
    record["n_skipped_ragged"] = n_skipped
    record["sha256"] = digest.hexdigest()

    stages = [
        record.get(k, 0.0)
        for k in (
            "seconds_cut_sites",
            "seconds_pileup",
            "seconds_loader",
            "seconds_epoch",
        )
    ]
    record["seconds_total"] = round(sum(stages), 2)
    record["torch_threads"] = torch.get_num_threads()
    return record


def peak_rss_mb():
    import resource

    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Linux reports kB, macOS bytes.
    return round(rss / (1024 if sys.platform.startswith("linux") else 1048576), 1)


###
# Commands
###


def cmd_once(args):
    """Worker: run a single (scale, arm) and write its JSON. Internal."""

    args = resolve(args)
    record = run_arm(args, args.scale, args.arm)
    record["peak_rss_mb"] = peak_rss_mb()
    record["status"] = "ok"
    with open(args.out, "w") as f:
        json.dump(record, f, indent=2, sort_keys=True)
    print(
        "[{} {}] {:.1f}s total".format(args.scale, args.arm, record["seconds_total"]),
        flush=True,
    )
    return 0


def cmd_run(args):
    """Orchestrator: prepare each rung, then run every arm in its own process."""

    args = resolve(args)
    for scale in args.scales:
        print("=== scale {}".format(scale), flush=True)
        shared = prepare_scale(args, scale)
        shared["scale"] = scale
        with open(os.path.join(args.results, "{}_shared.json".format(scale)), "w") as f:
            json.dump(shared, f, indent=2, sort_keys=True)
        print(
            "  {} peaks, {} negatives, negatives took {:.1f}s".format(
                shared["n_peaks"], shared["n_negatives"], shared["seconds_negatives"]
            ),
            flush=True,
        )

        for arm in args.arms:
            out = os.path.join(args.results, "{}_{}.json".format(scale, arm))
            cmd = (
                [
                    sys.executable,
                    os.path.abspath(__file__),
                    "once",
                    "--scale",
                    scale,
                    "--arm",
                    arm,
                    "--out",
                    out,
                    "--data",
                    args.data,
                    "--results",
                    args.results,
                    "--peaks",
                    args.peaks,
                    "--fasta",
                    args.fasta,
                    "--exclusion",
                    args.exclusion,
                    "--bams",
                ]
                + args.bams
                + [
                    "--in-window",
                    str(args.in_window),
                    "--out-window",
                    str(args.out_window),
                    "--max-jitter",
                    str(args.max_jitter),
                    "--batch-size",
                    str(args.batch_size),
                    "--random-state",
                    str(args.random_state),
                ]
            )
            if args.materialize:
                cmd.append("--materialize")
            if args.explain:
                cmd.append("--explain")
            if args.bam2bw:
                cmd += ["--bam2bw", args.bam2bw]

            started = time.time()
            try:
                proc = subprocess.run(
                    cmd, timeout=args.timeout, capture_output=True, text=True
                )
                status = "ok" if proc.returncode == 0 else "error"
                detail = proc.stderr[-4000:] if proc.returncode else ""
            except subprocess.TimeoutExpired:
                status, detail = "timeout", ""

            if status != "ok":
                with open(out, "w") as f:
                    json.dump(
                        {
                            "scale": scale,
                            "arm": arm,
                            "status": status,
                            "seconds_wall": round(time.time() - started, 1),
                            "timeout": args.timeout,
                            "stderr": detail,
                        },
                        f,
                        indent=2,
                        sort_keys=True,
                    )
                print(
                    "  {}: {} after {:.0f}s".format(
                        arm, status.upper(), time.time() - started
                    ),
                    flush=True,
                )
                if detail:
                    print(detail[-1500:], flush=True)
            else:
                print("  " + proc.stdout.strip(), flush=True)

    return cmd_report(args)


def load_results(args):
    records = []
    for scale in args.scales:
        for arm in args.arms:
            path = os.path.join(args.results, "{}_{}.json".format(scale, arm))
            if os.path.exists(path):
                with open(path) as f:
                    records.append(json.load(f))
    return records


def cmd_gate(args):
    """Assert every completed arm of a rung produced identical tensors."""

    args = resolve(args)
    records = [r for r in load_results(args) if r.get("status") == "ok"]
    by_scale = {}
    for r in records:
        by_scale.setdefault(r["scale"], []).append(r)

    ok = True
    for scale, group in sorted(by_scale.items()):
        digests = {r["arm"]: r.get("sha256") for r in group}
        agree = len(set(digests.values())) == 1
        ok &= agree
        print(
            "{}: {} ({})".format(
                scale,
                "identical" if agree else "DIFFERENT",
                ", ".join(
                    "{}={}".format(a, (d or "?")[:12])
                    for a, d in sorted(digests.items())
                ),
            )
        )
    return 0 if ok else 1


def _fmt(record, key):
    if record.get("status") != "ok":
        return record.get("status", "?").upper()
    value = record.get(key)
    return "-" if value is None else "{:.1f}".format(value)


def cmd_report(args):
    args = resolve(args)
    records = load_results(args)
    if not records:
        print("no results in {}".format(args.results))
        return 1

    lines = [
        "| scale | arm | cut sites | pileup | loader | epoch | total | status |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in records:
        lines.append(
            "| {} | `{}` | {} | {} | {} | {} | {} | {} |".format(
                r["scale"],
                r["arm"],
                _fmt(r, "seconds_cut_sites"),
                _fmt(r, "seconds_pileup"),
                _fmt(r, "seconds_loader"),
                _fmt(r, "seconds_epoch"),
                _fmt(r, "seconds_total"),
                r.get("status", "?"),
            )
        )

    table = "\n".join(lines)
    print(table)
    out = os.path.join(args.results, "TIMINGS.md")
    with open(out, "w") as f:
        f.write("# Pileup benchmark timings\n\nAll times in seconds.\n\n")
        f.write(table + "\n")
    print("\nwrote {}".format(out))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="command")

    run = sub.add_parser("run", help="Run the ladder and render timings.")
    add_common_args(run)
    run.add_argument(
        "--scales",
        nargs="+",
        default=["chr21", "chr8", "chr1", "genome"],
        choices=list(SCALES),
    )
    run.add_argument("--arms", nargs="+", default=list(ARMS), choices=list(ARMS))
    run.add_argument(
        "--timeout",
        type=float,
        default=1800.0,
        help="Per (scale, arm) wall-clock cap in seconds.",
    )
    run.set_defaults(func=cmd_run)

    once = sub.add_parser("once", help="Run one (scale, arm). Used internally by run.")
    add_common_args(once)
    once.add_argument("--scale", required=True, choices=list(SCALES))
    once.add_argument("--arm", required=True, choices=list(ARMS))
    once.add_argument("--out", required=True)
    once.set_defaults(func=cmd_once)

    for name, fn, helptext in (
        ("report", cmd_report, "Render timings from existing results."),
        ("gate", cmd_gate, "Check every arm of a rung agrees bit for bit."),
    ):
        p = sub.add_parser(name, help=helptext)
        add_common_args(p)
        p.add_argument("--scales", nargs="+", default=list(SCALES), choices=list(SCALES))
        p.add_argument("--arms", nargs="+", default=list(ARMS), choices=list(ARMS))
        p.set_defaults(func=fn)

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
