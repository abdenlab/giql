"""CLI entry point for the GIQL benchmarking suite.

Usage examples::

    pixi run python benchmarks/run.py --sizes 10000 --reps 2
    pixi run python benchmarks/run.py --input peaks.bed --reps 1
    pixi run python benchmarks/run.py --ops merge,cluster --output /tmp/results.csv
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

# Allow running as `python benchmarks/run.py` from the project root.
_repo_root = Path(__file__).parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from benchmarks.operations import ALL_OPS


def _detect_backends() -> tuple[bool, bool]:
    """Detect available backends. Returns (has_datafusion, has_bedtools)."""
    has_datafusion = False
    has_bedtools = False

    try:
        import datafusion  # noqa: F401

        has_datafusion = True
    except ImportError:
        warnings.warn(
            "datafusion not available — DataFusion benchmarks will be skipped.",
            stacklevel=1,
        )

    try:
        import pybedtools  # noqa: F401

        has_bedtools = True
    except ImportError:
        warnings.warn(
            "pybedtools not available — bedtools benchmarks will be skipped.",
            stacklevel=1,
        )

    return has_datafusion, has_bedtools


def _parse_sizes(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",")]


def _parse_ops(s: str) -> list[str]:
    ops = [x.strip() for x in s.split(",")]
    unknown = [o for o in ops if o not in ALL_OPS]
    if unknown:
        sys.exit(
            f"Unknown operation(s): {', '.join(unknown)}. "
            f"Valid: {', '.join(ALL_OPS)}"
        )
    return ops


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark GIQL DataFusion engine vs bedtools. "
            "The same DataFusionEngine instance is reused across all reps of a "
            "(size, op) pair — DataFusion's planning cache may benefit later reps. "
            "This is expected behavior, not an artifact."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--sizes",
        default="10000,100000,1000000",
        help="Comma-separated dataset sizes (ignored if --input given). Default: 10000,100000,1000000",
    )
    parser.add_argument(
        "--input",
        metavar="FILE",
        help=(
            "Path to a .bed or .parquet file to use as the primary dataset. "
            "Runs once at the file's actual size; --sizes is ignored."
        ),
    )
    parser.add_argument(
        "--reps",
        type=int,
        default=3,
        help="Number of timed repetitions per (engine, op, size). Default: 3",
    )
    parser.add_argument(
        "--output",
        metavar="FILE",
        help="Write results to this CSV file.",
    )
    parser.add_argument(
        "--ops",
        default=",".join(ALL_OPS),
        help=(
            f"Comma-separated operations to run. Default: all "
            f"({', '.join(ALL_OPS)})"
        ),
    )
    parser.add_argument(
        "--row-group-per-chrom",
        action="store_true",
        default=False,
        help=(
            "Write one Parquet row group per chromosome, enabling DataFusion "
            "row-group pruning via min/max statistics. Default: single row group "
            "(worst case)."
        ),
    )
    parser.add_argument(
        "--unsorted",
        action="store_true",
        default=False,
        help=(
            "Shuffle intervals within each chromosome (chromosomes remain "
            "contiguous). DataFusion merge/cluster queries will include an "
            "ORDER BY subquery to sort before aggregating."
        ),
    )
    args = parser.parse_args()

    has_datafusion, has_bedtools = _detect_backends()
    if not has_datafusion and not has_bedtools:
        sys.exit("No backends available. Install datafusion and/or pybedtools.")

    ops = _parse_ops(args.ops)

    from benchmarks.runner import print_results_table, run_benchmark, write_results_csv

    input_path: Path | None = None
    if args.input:
        input_path = Path(args.input)
        if not input_path.exists():
            sys.exit(f"File not found: {input_path}")
        if input_path.suffix not in (".bed", ".parquet"):
            sys.exit(f"Unknown file format: {input_path.suffix!r}. Use .bed or .parquet")

    sizes = _parse_sizes(args.sizes)

    results = run_benchmark(
        sizes=sizes,
        ops=ops,
        n_reps=args.reps,
        input_path=input_path,
        has_datafusion=has_datafusion,
        has_bedtools=has_bedtools,
        row_group_per_chrom=args.row_group_per_chrom,
        unsorted=args.unsorted,
    )

    print_results_table(results)

    if args.output:
        output_path = Path(args.output)
        write_results_csv(results, output_path)
        print(f"\nResults written to {output_path}")


if __name__ == "__main__":
    main()
