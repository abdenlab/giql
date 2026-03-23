"""CLI entry point for the bin-size optimization experiment.

Usage::

    pixi run python benchmarks/experiments/bin_size/run.py --workers 4 --reps 3
    pixi run python benchmarks/experiments/bin_size/run.py \\
        --workers 1 --reps 1 --export-csv results.csv
"""

from __future__ import annotations

import argparse
import asyncio
import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parent.parent.parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

MEDIAN_LENGTHS = [100, 500, 1000, 5000, 10000, 50000]
SIGMAS = [0.5, 1.0, 2.0]
N_INTERVALS_LIST = [100_000, 250_000, 500_000, 1_000_000]
BIN_SIZES = [1000, 5000, 10000, 50000, 100000, 500000]


def _distributions():
    """Optuna parameter distributions for the search space."""
    import optuna

    return {
        "median_length": (
            optuna.distributions.CategoricalDistribution(
                MEDIAN_LENGTHS,
            )
        ),
        "sigma": (
            optuna.distributions.CategoricalDistribution(
                SIGMAS,
            )
        ),
        "n_intervals": (
            optuna.distributions.CategoricalDistribution(
                N_INTERVALS_LIST,
            )
        ),
        "bin_size": (
            optuna.distributions.CategoricalDistribution(
                BIN_SIZES,
            )
        ),
    }


def create_study(db_path: str):
    """Create an Optuna study backed by SQLite."""
    import optuna

    storage = f"sqlite:///{db_path}"
    sampler = optuna.samplers.GridSampler(
        {
            "median_length": MEDIAN_LENGTHS,
            "sigma": SIGMAS,
            "n_intervals": N_INTERVALS_LIST,
            "bin_size": BIN_SIZES,
        }
    )
    return optuna.create_study(
        study_name="bin_size_sweep",
        storage=storage,
        sampler=sampler,
        direction="minimize",
        load_if_exists=True,
    )


async def _consume_routine(study, routine, distributions):
    """Consume results from a wool routine and add to study."""
    import optuna

    async for result in routine:
        state = (
            optuna.trial.TrialState.COMPLETE
            if result.n_pairs >= 0
            else optuna.trial.TrialState.FAIL
        )
        trial = optuna.trial.create_trial(
            state=state,
            params={
                "median_length": result.median_length,
                "sigma": result.sigma,
                "n_intervals": result.n_intervals,
                "bin_size": result.bin_size,
            },
            distributions=distributions,
            values=[result.time_s],
            user_attrs={
                "rep": result.rep,
                "n_pairs": result.n_pairs,
            },
        )
        study.add_trial(trial)
        print(
            f"  ml={result.median_length:>5} "
            f"σ={result.sigma:.1f} "
            f"n={result.n_intervals:>7} "
            f"b={result.bin_size:>6} "
            f"rep={result.rep} "
            f"→ {result.time_s:.3f}s "
            f"({result.n_pairs} pairs)"
        )


async def run_experiment(
    study, n_reps: int, n_workers: int,
):
    """Dispatch wool routines for all dataset configs."""
    import wool

    from benchmarks.experiments.bin_size.routine import bin_size_sweep

    distributions = _distributions()

    async with wool.WorkerPool(size=n_workers):
        async with asyncio.TaskGroup() as tg:
            for median_length in MEDIAN_LENGTHS:
                for sigma in SIGMAS:
                    for n_intervals in N_INTERVALS_LIST:
                        routine = bin_size_sweep(
                            median_length=median_length,
                            sigma=sigma,
                            n_intervals=n_intervals,
                            bin_sizes=BIN_SIZES,
                            n_reps=n_reps,
                        )
                        tg.create_task(
                            _consume_routine(
                                study,
                                routine,
                                distributions,
                            )
                        )


def export_results_csv(study, path: str) -> None:
    """Export study trials to CSV."""
    df = study.trials_dataframe()
    df.to_csv(path, index=False)
    print(f"Results exported to {path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Bin-size optimization experiment for "
            "INTERSECTS join."
        ),
    )
    parser.add_argument(
        "--db",
        default="bin_size_study.db",
        help=(
            "SQLite database for Optuna study. "
            "Default: bin_size_study.db"
        ),
    )
    parser.add_argument(
        "--reps",
        type=int,
        default=3,
        help=(
            "Timed repetitions per (config, bin_size). "
            "Default: 3"
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of wool workers. Default: 4",
    )
    parser.add_argument(
        "--export-csv",
        metavar="FILE",
        help="Export results to CSV after the run.",
    )
    args = parser.parse_args()

    study = create_study(args.db)

    n_configs = (
        len(MEDIAN_LENGTHS)
        * len(SIGMAS)
        * len(N_INTERVALS_LIST)
    )
    n_total = n_configs * len(BIN_SIZES) * args.reps
    print(
        f"Running bin-size experiment: "
        f"{n_configs} configs × {len(BIN_SIZES)} bin_sizes "
        f"× {args.reps} reps = {n_total} trials"
    )
    print(f"Workers: {args.workers}, DB: {args.db}")

    asyncio.run(
        run_experiment(study, args.reps, args.workers)
    )

    if args.export_csv:
        export_results_csv(study, args.export_csv)


if __name__ == "__main__":
    main()
