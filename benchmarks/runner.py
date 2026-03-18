"""Timing harness and result formatting for GIQL benchmarks."""

from __future__ import annotations

import math
import tempfile
from collections.abc import Callable
from pathlib import Path
from time import perf_counter
from typing import NamedTuple

from benchmarks.data import GenomicDataset, generate_dataset, load_from_bed, load_from_parquet
from benchmarks.engines import (
    make_datafusion_engine,
    run_bedtools_cluster,
    run_bedtools_cluster_unsorted,
    run_bedtools_intersect_filter,
    run_bedtools_intersect_join,
    run_bedtools_merge,
    run_bedtools_merge_unsorted,
    run_datafusion_cluster,
    run_datafusion_cluster_unsorted,
    run_datafusion_intersect_filter,
    run_datafusion_intersect_join,
    run_datafusion_merge,
    run_datafusion_merge_unsorted,
)
from benchmarks.operations import ALL_OPS, QUERY_REGION_BED


class BenchmarkResult(NamedTuple):
    operation: str
    engine: str    # "datafusion" | "bedtools"
    size: int
    mean_s: float
    std_s: float
    n_reps: int
    status: str    # "ok" | "na" | "error:<msg>"


def time_reps(fn: Callable[[], None], n_reps: int) -> tuple[float, float]:
    """Time n_reps calls to fn, returning (mean_s, std_s).

    No warmup rep is discarded — I/O-bound operations are affected by the OS
    page cache, and including the cold-start rep reflects real-world behavior.
    """
    times: list[float] = []
    for _ in range(n_reps):
        t0 = perf_counter()
        fn()
        times.append(perf_counter() - t0)
    mean = sum(times) / len(times)
    variance = sum((t - mean) ** 2 for t in times) / len(times)
    return mean, math.sqrt(variance)


def _make_query_bed(workdir: Path) -> Path:
    """Write a single-interval BED file for the intersect filter region."""
    chrom, start, end = QUERY_REGION_BED
    path = workdir / "query_region.bed"
    path.write_text(f"{chrom}\t{start}\t{end}\n")
    return path


def run_benchmark(
    sizes: list[int],
    ops: list[str],
    n_reps: int,
    input_path: Path | None,
    has_datafusion: bool,
    has_bedtools: bool,
    *,
    row_group_per_chrom: bool = False,
    unsorted: bool = False,
) -> list[BenchmarkResult]:
    """Run benchmarks and return results.

    Parameters
    ----------
    sizes :
        Dataset sizes to benchmark. Ignored when input_path is given.
    ops :
        Operation keys to run.
    n_reps :
        Number of timed repetitions per (engine, op, size).
    input_path :
        Path to a .bed or .parquet file to use as the primary dataset.
        When given, runs once at the file's actual size; sizes is ignored.
    has_datafusion :
        Whether the DataFusion backend is available.
    has_bedtools :
        Whether the bedtools backend is available.
    row_group_per_chrom :
        If True, write one Parquet row group per chromosome so DataFusion
        can prune row groups via min/max statistics. Default False (single
        row group — worst case for pruning).
    unsorted :
        If True, intervals within each chromosome are shuffled (chromosomes
        remain contiguous). DataFusion merge/cluster queries include an
        ORDER BY subquery; bedtools handles unsorted data natively.
    """
    results: list[BenchmarkResult] = []

    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)
        query_bed = _make_query_bed(workdir)

        if input_path is not None:
            if input_path.suffix == ".bed":
                primary = load_from_bed(input_path, workdir, row_group_per_chrom=row_group_per_chrom)
            else:
                primary = load_from_parquet(input_path, workdir, row_group_per_chrom=row_group_per_chrom)
            secondary = generate_dataset(
                primary.n_intervals, workdir, name="secondary", seed=43,
                row_group_per_chrom=row_group_per_chrom, unsorted=unsorted,
            )
            dataset_pairs: list[tuple[GenomicDataset, GenomicDataset]] = [(primary, secondary)]
        else:
            dataset_pairs = []
            for size in sizes:
                primary = generate_dataset(
                    size, workdir, name=f"primary_{size}", seed=42,
                    row_group_per_chrom=row_group_per_chrom, unsorted=unsorted,
                )
                secondary = generate_dataset(
                    size, workdir, name=f"secondary_{size}", seed=43,
                    row_group_per_chrom=row_group_per_chrom, unsorted=unsorted,
                )
                dataset_pairs.append((primary, secondary))

        for primary, secondary in dataset_pairs:
            size = primary.n_intervals

            # Create DataFusion engines once per table-set, reused across all ops and reps.
            df_engine_single = None
            df_engine_pair = None
            if has_datafusion:
                df_engine_single = make_datafusion_engine(primary)
                df_engine_pair = make_datafusion_engine(primary, secondary)

            for op_key in ops:
                if op_key not in ALL_OPS:
                    continue
                op = ALL_OPS[op_key]

                # --- DataFusion ---
                if has_datafusion and op.datafusion_supported:
                    engine = df_engine_pair if op.needs_secondary else df_engine_single
                    try:
                        if op_key == "intersect_filter":
                            fn: Callable[[], None] = lambda e=engine: run_datafusion_intersect_filter(e)
                        elif op_key == "intersect_join":
                            fn = lambda e=engine: run_datafusion_intersect_join(e)
                        elif op_key == "merge":
                            merge_fn = run_datafusion_merge_unsorted if unsorted else run_datafusion_merge
                            fn = lambda e=engine, f=merge_fn: f(e)
                        elif op_key == "cluster":
                            cluster_fn = run_datafusion_cluster_unsorted if unsorted else run_datafusion_cluster
                            fn = lambda e=engine, f=cluster_fn: f(e)
                        else:
                            fn = None  # type: ignore[assignment]

                        if fn is not None:
                            mean_s, std_s = time_reps(fn, n_reps)
                            results.append(BenchmarkResult(
                                operation=op_key,
                                engine="datafusion",
                                size=size,
                                mean_s=mean_s,
                                std_s=std_s,
                                n_reps=n_reps,
                                status="ok",
                            ))
                        else:
                            results.append(BenchmarkResult(
                                operation=op_key, engine="datafusion", size=size,
                                mean_s=0.0, std_s=0.0, n_reps=0, status="na",
                            ))
                    except Exception as exc:
                        results.append(BenchmarkResult(
                            operation=op_key, engine="datafusion", size=size,
                            mean_s=0.0, std_s=0.0, n_reps=0,
                            status=f"error:{exc!s}",
                        ))
                else:
                    results.append(BenchmarkResult(
                        operation=op_key, engine="datafusion", size=size,
                        mean_s=0.0, std_s=0.0, n_reps=0, status="na",
                    ))

                # --- bedtools ---
                if has_bedtools:
                    try:
                        import pybedtools  # noqa: F401

                        if op_key == "intersect_filter":
                            fn = lambda p=primary, q=query_bed: run_bedtools_intersect_filter(p.bed_path, q)
                        elif op_key == "intersect_join":
                            fn = lambda p=primary, s=secondary: run_bedtools_intersect_join(p.bed_path, s.bed_path)
                        elif op_key == "merge":
                            bt_merge = run_bedtools_merge_unsorted if unsorted else run_bedtools_merge
                            fn = lambda p=primary, f=bt_merge: f(p.bed_path)
                        elif op_key == "cluster":
                            bt_cluster = run_bedtools_cluster_unsorted if unsorted else run_bedtools_cluster
                            fn = lambda p=primary, f=bt_cluster: f(p.bed_path)
                        else:
                            fn = None  # type: ignore[assignment]

                        if fn is not None:
                            mean_s, std_s = time_reps(fn, n_reps)
                            results.append(BenchmarkResult(
                                operation=op_key,
                                engine="bedtools",
                                size=size,
                                mean_s=mean_s,
                                std_s=std_s,
                                n_reps=n_reps,
                                status="ok",
                            ))
                            pybedtools.cleanup(remove_all=False)
                        else:
                            results.append(BenchmarkResult(
                                operation=op_key, engine="bedtools", size=size,
                                mean_s=0.0, std_s=0.0, n_reps=0, status="na",
                            ))
                    except Exception as exc:
                        results.append(BenchmarkResult(
                            operation=op_key, engine="bedtools", size=size,
                            mean_s=0.0, std_s=0.0, n_reps=0,
                            status=f"error:{exc!s}",
                        ))
                else:
                    results.append(BenchmarkResult(
                        operation=op_key, engine="bedtools", size=size,
                        mean_s=0.0, std_s=0.0, n_reps=0, status="na",
                    ))

    return results


def print_results_table(results: list[BenchmarkResult]) -> None:
    """Print an aligned console table of benchmark results."""
    # Column widths
    col_op = 20
    col_size = 12
    col_df = 20
    col_bt = 20
    col_sp = 10

    header = (
        "Operation".ljust(col_op)
        + "Size".rjust(col_size)
        + "DataFusion (s)".rjust(col_df)
        + "bedtools (s)".rjust(col_bt)
        + "Speedup".rjust(col_sp)
    )
    sep = "-" * len(header)

    print(sep)
    print(header)
    print(sep)

    # Group by (operation, size) preserving insertion order
    seen: list[tuple[str, int]] = []
    grouped: dict[tuple[str, int], dict[str, BenchmarkResult]] = {}
    for r in results:
        key = (r.operation, r.size)
        if key not in grouped:
            seen.append(key)
            grouped[key] = {}
        grouped[key][r.engine] = r

    for key in seen:
        op_key, size = key
        engines = grouped[key]
        label = ALL_OPS[op_key].label if op_key in ALL_OPS else op_key
        df_r = engines.get("datafusion")
        bt_r = engines.get("bedtools")

        def fmt_time(r: BenchmarkResult | None) -> str:
            if r is None or r.status == "na":
                return "N/A"
            if r.status.startswith("error:"):
                return "ERROR"
            return f"{r.mean_s:.4f} ± {r.std_s:.4f}"

        df_str = fmt_time(df_r)
        bt_str = fmt_time(bt_r)

        speedup_str = ""
        if (
            df_r and bt_r
            and df_r.status == "ok" and bt_r.status == "ok"
            and df_r.mean_s > 0
        ):
            speedup_str = f"{bt_r.mean_s / df_r.mean_s:.2f}x"

        print(
            label.ljust(col_op)
            + f"{size:,}".rjust(col_size)
            + df_str.rjust(col_df)
            + bt_str.rjust(col_bt)
            + speedup_str.rjust(col_sp)
        )

    print(sep)
    print("Speedup > 1 means DataFusion is faster.")


def write_results_csv(results: list[BenchmarkResult], path: Path) -> None:
    """Write benchmark results to a CSV file."""
    with open(path, "w") as f:
        f.write("operation,engine,size,mean_s,std_s,n_reps,status\n")
        for r in results:
            f.write(
                f"{r.operation},{r.engine},{r.size},"
                f"{r.mean_s:.6f},{r.std_s:.6f},{r.n_reps},{r.status}\n"
            )
