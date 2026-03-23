"""Analysis and visualization for bin-size experiment results.

Usage::

    pixi run python benchmarks/experiments/bin_size/analysis.py --db bin_size_study.db
    pixi run python benchmarks/experiments/bin_size/analysis.py \\
        --db bin_size_study.db --output-dir plots/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Total hg38 genome size (sum of benchmarks.data.HG38_CHROM_SIZES).
_GENOME_SIZE = 3_088_269_832


def load_study(db_path: str):
    """Load an Optuna study from SQLite."""
    import optuna

    return optuna.load_study(
        study_name="bin_size_sweep",
        storage=f"sqlite:///{db_path}",
    )


def results_dataframe(study) -> pd.DataFrame:
    """Convert study trials to a tidy DataFrame with derived metrics."""
    records = []
    for trial in study.trials:
        if trial.state.name != "COMPLETE":
            continue
        records.append({
            "median_length": trial.params["median_length"],
            "sigma": trial.params["sigma"],
            "n_intervals": trial.params["n_intervals"],
            "bin_size": trial.params["bin_size"],
            "time_s": trial.values[0],
            "rep": trial.user_attrs.get("rep", 0),
            "n_pairs": trial.user_attrs.get("n_pairs", -1),
        })

    df = pd.DataFrame(records)

    # Lognormal mean = median * exp(sigma^2 / 2)
    df["mean_length"] = df["median_length"] * np.exp(
        df["sigma"] ** 2 / 2
    )
    df["expansion_factor"] = (
        df["mean_length"] / df["bin_size"] + 1
    )
    df["density"] = (
        df["n_intervals"] * df["mean_length"] / _GENOME_SIZE
    )
    df["bin_occupancy"] = (
        df["n_intervals"]
        * df["expansion_factor"]
        / (_GENOME_SIZE / df["bin_size"])
    )

    return df


def plot_heatmaps(df: pd.DataFrame, output_dir: Path) -> None:
    """Plot time_s vs (median_length, bin_size) heatmaps.

    One heatmap per (sigma, n_intervals) combination.
    """
    agg = (
        df.groupby([
            "median_length",
            "sigma",
            "n_intervals",
            "bin_size",
        ])
        .agg(time_mean=("time_s", "mean"))
        .reset_index()
    )

    for sigma in sorted(df["sigma"].unique()):
        for n_int in sorted(df["n_intervals"].unique()):
            subset = agg[
                (agg["sigma"] == sigma)
                & (agg["n_intervals"] == n_int)
            ]
            if subset.empty:
                continue

            pivot = subset.pivot(
                index="bin_size",
                columns="median_length",
                values="time_mean",
            )

            fig, ax = plt.subplots(figsize=(10, 7))
            im = ax.pcolormesh(
                range(len(pivot.columns) + 1),
                range(len(pivot.index) + 1),
                np.log10(pivot.values),
                cmap="viridis",
                shading="flat",
            )
            ax.set_xticks(
                np.arange(len(pivot.columns)) + 0.5,
                [str(c) for c in pivot.columns],
            )
            ax.set_yticks(
                np.arange(len(pivot.index)) + 0.5,
                [str(i) for i in pivot.index],
            )
            ax.set_xlabel("median_length (bp)")
            ax.set_ylabel("bin_size (bp)")
            ax.set_title(
                f"log10(time_s) — "
                f"sigma={sigma}, n={n_int:,}"
            )
            fig.colorbar(im, ax=ax, label="log10(time_s)")

            for i, bs in enumerate(pivot.index):
                for j, ml in enumerate(pivot.columns):
                    val = pivot.loc[bs, ml]
                    if not np.isnan(val):
                        ax.text(
                            j + 0.5,
                            i + 0.5,
                            f"{val:.2f}",
                            ha="center",
                            va="center",
                            fontsize=7,
                            color="white",
                        )

            fig.tight_layout()
            fname = (
                f"heatmap_sigma{sigma}_n{n_int}.png"
            )
            fig.savefig(output_dir / fname, dpi=150)
            plt.close(fig)
            print(f"  Saved {fname}")


def plot_optimal_bin_size(
    df: pd.DataFrame, output_dir: Path,
) -> None:
    """Plot optimal bin_size vs median_length (log-log)."""
    agg = (
        df.groupby([
            "median_length",
            "sigma",
            "n_intervals",
            "bin_size",
        ])
        .agg(time_mean=("time_s", "mean"))
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(10, 7))
    markers = {0.5: "o", 1.0: "s", 2.0: "^"}
    colors = {100_000: "tab:blue", 1_000_000: "tab:red"}

    for sigma in sorted(agg["sigma"].unique()):
        for n_int in sorted(agg["n_intervals"].unique()):
            subset = agg[
                (agg["sigma"] == sigma)
                & (agg["n_intervals"] == n_int)
            ]
            if subset.empty:
                continue

            optimal = (
                subset.loc[
                    subset.groupby("median_length")[
                        "time_mean"
                    ].idxmin()
                ].sort_values("median_length")
            )
            ax.plot(
                optimal["median_length"],
                optimal["bin_size"],
                marker=markers.get(sigma, "o"),
                color=colors.get(n_int, "tab:gray"),
                label=f"σ={sigma}, n={n_int:,}",
                linewidth=1.5,
                markersize=6,
            )

    ml_range = sorted(agg["median_length"].unique())
    ax.plot(
        ml_range,
        ml_range,
        "--",
        color="gray",
        alpha=0.5,
        label="bin = median_length",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("median_length (bp)")
    ax.set_ylabel("optimal bin_size (bp)")
    ax.set_title(
        "Optimal bin size vs median interval length"
    )
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(
        output_dir / "optimal_bin_size.png", dpi=150,
    )
    plt.close(fig)
    print("  Saved optimal_bin_size.png")


def plot_ratio(
    df: pd.DataFrame, output_dir: Path,
) -> None:
    """Plot optimal_bin_size / p90(length) ratio."""
    agg = (
        df.groupby([
            "median_length",
            "sigma",
            "n_intervals",
            "bin_size",
        ])
        .agg(time_mean=("time_s", "mean"))
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(10, 7))
    markers = {0.5: "o", 1.0: "s", 2.0: "^"}
    colors = {100_000: "tab:blue", 1_000_000: "tab:red"}

    for sigma in sorted(agg["sigma"].unique()):
        for n_int in sorted(agg["n_intervals"].unique()):
            subset = agg[
                (agg["sigma"] == sigma)
                & (agg["n_intervals"] == n_int)
            ]
            if subset.empty:
                continue

            optimal = (
                subset.loc[
                    subset.groupby("median_length")[
                        "time_mean"
                    ].idxmin()
                ].sort_values("median_length")
            )
            # Lognormal p90 = median * exp(z_0.9 * sigma)
            p90 = optimal["median_length"] * np.exp(
                sigma * 1.2816
            )
            ratio = optimal["bin_size"] / p90

            ax.plot(
                optimal["median_length"],
                ratio,
                marker=markers.get(sigma, "o"),
                color=colors.get(n_int, "tab:gray"),
                label=f"σ={sigma}, n={n_int:,}",
                linewidth=1.5,
                markersize=6,
            )

    ax.axhline(
        1.0, linestyle="--", color="gray", alpha=0.5,
        label="ratio = 1",
    )
    ax.set_xscale("log")
    ax.set_xlabel("median_length (bp)")
    ax.set_ylabel("optimal_bin_size / p90(length)")
    ax.set_title(
        "Databricks heuristic test: optimal / p90 ratio"
    )
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_dir / "ratio_plot.png", dpi=150)
    plt.close(fig)
    print("  Saved ratio_plot.png")


def plot_expansion_vs_time(
    df: pd.DataFrame, output_dir: Path,
) -> None:
    """Plot time_s vs expansion_factor, colored by density."""
    agg = (
        df.groupby([
            "median_length",
            "sigma",
            "n_intervals",
            "bin_size",
        ])
        .agg(
            time_mean=("time_s", "mean"),
            expansion_factor=("expansion_factor", "first"),
            density=("density", "first"),
        )
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(10, 7))
    scatter = ax.scatter(
        agg["expansion_factor"],
        agg["time_mean"],
        c=np.log10(agg["density"]),
        cmap="coolwarm",
        alpha=0.7,
        s=20,
    )
    fig.colorbar(scatter, ax=ax, label="log10(density)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(
        "expansion_factor (mean_length / bin_size + 1)"
    )
    ax.set_ylabel("time_s")
    ax.set_title("Expansion factor vs execution time")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(
        output_dir / "expansion_vs_time.png", dpi=150,
    )
    plt.close(fig)
    print("  Saved expansion_vs_time.png")


def plot_regression(
    df: pd.DataFrame, output_dir: Path,
) -> tuple[np.ndarray, float]:
    """Fit log(optimal_b) = a*log(mu) + b*log(G/n) + c."""
    agg = (
        df.groupby([
            "median_length",
            "sigma",
            "n_intervals",
            "bin_size",
        ])
        .agg(
            time_mean=("time_s", "mean"),
            mean_length=("mean_length", "first"),
        )
        .reset_index()
    )

    optimal = agg.loc[
        agg.groupby([
            "median_length", "sigma", "n_intervals",
        ])["time_mean"].idxmin()
    ]

    log_b = np.log(
        optimal["bin_size"].values
    ).astype(float)
    log_mu = np.log(
        optimal["mean_length"].values
    ).astype(float)
    log_sparsity = np.log(
        _GENOME_SIZE / optimal["n_intervals"].values
    ).astype(float)

    A = np.column_stack([
        log_mu, log_sparsity, np.ones(len(log_b)),
    ])
    coeffs, _, _, _ = np.linalg.lstsq(A, log_b, rcond=None)

    predicted = A @ coeffs
    ss_res = np.sum((log_b - predicted) ** 2)
    ss_tot = np.sum((log_b - np.mean(log_b)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(predicted, log_b, alpha=0.7, s=30)
    lims = [
        min(predicted.min(), log_b.min()),
        max(predicted.max(), log_b.max()),
    ]
    ax.plot(lims, lims, "--", color="gray", alpha=0.5)
    ax.set_xlabel("predicted log(optimal_bin_size)")
    ax.set_ylabel("actual log(optimal_bin_size)")
    ax.set_title(
        f"Heuristic regression "
        f"(R\u00b2={r_squared:.3f})\n"
        f"log(b) = {coeffs[0]:.2f}\u00b7log(\u03bc) "
        f"+ {coeffs[1]:.2f}\u00b7log(G/n) "
        f"+ {coeffs[2]:.2f}"
    )
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_dir / "regression.png", dpi=150)
    plt.close(fig)
    print("  Saved regression.png")
    print(
        f"  Fit: log(b) = {coeffs[0]:.3f}*log(mu) "
        f"+ {coeffs[1]:.3f}*log(G/n) "
        f"+ {coeffs[2]:.3f} "
        f"(R²={r_squared:.3f})"
    )
    return coeffs, r_squared


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze bin-size experiment results.",
    )
    parser.add_argument(
        "--db",
        default="bin_size_study.db",
        help=(
            "Optuna SQLite database. "
            "Default: bin_size_study.db"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="plots",
        help="Directory for plot output. Default: plots/",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading study from {args.db}...")
    study = load_study(args.db)
    df = results_dataframe(study)
    print(f"  {len(df)} trials loaded")

    print("Plotting heatmaps...")
    plot_heatmaps(df, output_dir)

    print("Plotting optimal bin size...")
    plot_optimal_bin_size(df, output_dir)

    print("Plotting ratio (Databricks heuristic test)...")
    plot_ratio(df, output_dir)

    print("Plotting expansion factor vs time...")
    plot_expansion_vs_time(df, output_dir)

    print("Fitting regression...")
    plot_regression(df, output_dir)


if __name__ == "__main__":
    main()
