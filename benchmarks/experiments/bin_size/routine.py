"""Wool routine for bin-size sweep experiment."""

from __future__ import annotations

import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

import wool

_REPO_ROOT = str(Path(__file__).resolve().parent.parent.parent.parent)

_TIMEOUT_S = 300


@dataclass(frozen=True)
class TrialResult:
    """Result of a single timed trial."""

    median_length: int
    sigma: float
    n_intervals: int
    bin_size: int
    rep: int
    time_s: float
    n_pairs: int


def build_binned_join_sql(bin_size: int) -> str:
    """Build INTERSECTS full-pair join SQL with parameterized bin size.

    Extracted from ``IntersectsJoinTransformer._build_binned_query``
    in ``src/giql/transformer.py``, simplified for the ``a``/``b``
    table convention used by the experiment.
    """
    b = bin_size
    return (
        "WITH __giql_left AS ("
        "SELECT *, "
        f'UNNEST(range(CAST("start" / {b} AS BIGINT), '
        f'CAST(("end" - 1) / {b} + 1 AS BIGINT))) '
        "AS __giql_bin "
        "FROM a), "
        "__giql_right AS ("
        "SELECT *, "
        f'UNNEST(range(CAST("start" / {b} AS BIGINT), '
        f'CAST(("end" - 1) / {b} + 1 AS BIGINT))) '
        "AS __giql_bin "
        "FROM b) "
        "SELECT DISTINCT "
        'l."chrom", l."start", l."end", '
        'r."chrom" AS chrom_r, '
        'r."start" AS start_r, '
        'r."end" AS end_r '
        "FROM __giql_left AS l "
        "JOIN __giql_right AS r "
        'ON l."chrom" = r."chrom" '
        "AND l.__giql_bin = r.__giql_bin "
        'WHERE l."start" < r."end" '
        'AND l."end" > r."start"'
    )


def _make_engine(primary_path: Path, secondary_path: Path):
    """Create a fresh DataFusion SessionContext with tables a/b."""
    import datafusion

    ctx = datafusion.SessionContext()
    ctx.register_parquet("a", str(primary_path))
    ctx.register_parquet("b", str(secondary_path))
    return ctx


@wool.routine
async def bin_size_sweep(
    median_length: int,
    sigma: float,
    n_intervals: int,
    bin_sizes: list[int],
    n_reps: int,
    seed: int = 42,
):
    """Sweep bin sizes for one dataset configuration.

    Generates one dataset pair per (median_length, sigma, n_intervals),
    then sweeps all bin_sizes x n_reps, yielding TrialResult
    incrementally.
    """
    if _REPO_ROOT not in sys.path:
        sys.path.insert(0, _REPO_ROOT)
    from benchmarks.data import generate_dataset

    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)

        primary = generate_dataset(
            n_intervals,
            workdir,
            name="primary",
            seed=seed,
            median_length=median_length,
            sigma=sigma,
        )
        secondary = generate_dataset(
            n_intervals,
            workdir,
            name="secondary",
            seed=seed + 1,
            median_length=median_length,
            sigma=sigma,
        )

        for bin_size in bin_sizes:
            if median_length / bin_size > 100:
                continue
            # Avoid OOM: skip small bins for large datasets.
            if n_intervals >= 1_000_000 and bin_size < 10000:
                continue
            if n_intervals >= 500_000 and bin_size < 5000:
                continue

            sql = build_binned_join_sql(bin_size)

            for rep in range(n_reps):
                ctx = _make_engine(
                    primary.parquet_path,
                    secondary.parquet_path,
                )
                t0 = perf_counter()
                try:
                    batches = ctx.sql(sql).collect()
                    elapsed = perf_counter() - t0
                    n_pairs = sum(
                        b.num_rows for b in batches
                    )
                except Exception:
                    elapsed = perf_counter() - t0
                    n_pairs = -1

                yield TrialResult(
                    median_length=median_length,
                    sigma=sigma,
                    n_intervals=n_intervals,
                    bin_size=bin_size,
                    rep=rep,
                    time_s=elapsed,
                    n_pairs=n_pairs,
                )

                if elapsed > _TIMEOUT_S:
                    break
