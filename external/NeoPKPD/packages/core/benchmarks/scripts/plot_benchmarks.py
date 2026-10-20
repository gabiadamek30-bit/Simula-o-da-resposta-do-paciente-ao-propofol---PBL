#!/usr/bin/env python3
"""
NeoPKPD Benchmark Plotting Script
==================================

Aggregates benchmark results from all platforms and generates
publication-quality figures for the SoftwareX paper.

Usage:
    python plot_benchmarks.py

Output:
    benchmarks/figures/
        - figure_simulation_comparison.png
        - figure_population_scaling.png
        - figure_speedup_summary.png
        - benchmark_summary_table.csv (for LaTeX)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import sys

# =============================================================================
# Configuration
# =============================================================================

RESULTS_DIR = Path(__file__).parent.parent / "results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"

# Color palette (colorblind-friendly)
COLORS = {
    "NeoPKPD (Julia)": "#2563EB",    # Blue
    "NeoPKPD (Python)": "#3B82F6",   # Light blue
    "nlmixr2": "#059669",            # Green
    "mrgsolve": "#D97706",           # Orange
}

# Plot style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# =============================================================================
# Data Loading
# =============================================================================

def load_benchmark_data():
    """Load all benchmark CSV files."""
    data = {}

    files = {
        "NeoPKPD (Julia)": "neopkpd_julia_benchmarks.csv",
        "NeoPKPD (Python)": "neopkpd_python_benchmarks.csv",
        "nlmixr2": "nlmixr2_benchmarks.csv",
        "mrgsolve": "mrgsolve_benchmarks.csv",
    }

    for platform, filename in files.items():
        filepath = RESULTS_DIR / filename
        if filepath.exists():
            df = pd.read_csv(filepath)
            df["platform"] = platform
            data[platform] = df
            print(f"Loaded: {filename} ({len(df)} rows)")
        else:
            print(f"Warning: {filename} not found")

    if not data:
        print("Error: No benchmark data found!")
        sys.exit(1)

    return pd.concat(data.values(), ignore_index=True)


# =============================================================================
# Figure 1: Single Simulation Comparison
# =============================================================================

def plot_simulation_comparison(df):
    """Bar chart comparing single simulation times across platforms."""
    print("\nGenerating: Single Simulation Comparison...")

    # Filter for single simulation benchmarks
    sim_df = df[
        (df["category"] == "simulation") &
        (df["name"] == "single_simulation")
    ].copy()

    if sim_df.empty:
        print("  No simulation data available")
        return

    # Models to compare
    models = ["OneCompIV", "TwoCompIV", "TwoCompOral", "MichaelisMenten"]
    models = [m for m in models if m in sim_df["model"].values]

    platforms = [p for p in COLORS.keys() if p in sim_df["platform"].values]

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(models))
    width = 0.2
    n_platforms = len(platforms)

    for i, platform in enumerate(platforms):
        platform_data = sim_df[sim_df["platform"] == platform]
        means = []
        errors = []

        for model in models:
            model_data = platform_data[platform_data["model"] == model]
            if not model_data.empty:
                means.append(model_data["mean_ms"].values[0])
                errors.append(model_data["std_ms"].values[0])
            else:
                means.append(0)
                errors.append(0)

        offset = (i - n_platforms/2 + 0.5) * width
        bars = ax.bar(x + offset, means, width, yerr=errors,
                      label=platform, color=COLORS[platform],
                      capsize=3, alpha=0.9)

    ax.set_xlabel("Model Type")
    ax.set_ylabel("Time (ms)")
    ax.set_title("Single Simulation Performance Comparison")
    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=15, ha="right")
    ax.legend(loc="upper left")
    ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.3)

    # Add speedup annotations
    for i, model in enumerate(models):
        neopkpd_time = sim_df[
            (sim_df["platform"] == "NeoPKPD (Julia)") &
            (sim_df["model"] == model)
        ]["mean_ms"].values

        if len(neopkpd_time) > 0:
            for platform in ["nlmixr2", "mrgsolve"]:
                other_time = sim_df[
                    (sim_df["platform"] == platform) &
                    (sim_df["model"] == model)
                ]["mean_ms"].values

                if len(other_time) > 0 and neopkpd_time[0] > 0:
                    speedup = other_time[0] / neopkpd_time[0]
                    if speedup > 1.5:
                        ax.annotate(f"{speedup:.1f}x",
                                    xy=(i, other_time[0]),
                                    xytext=(0, 5),
                                    textcoords="offset points",
                                    ha="center", fontsize=7,
                                    color=COLORS[platform])

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_simulation_comparison.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


# =============================================================================
# Figure 2: Population Scaling
# =============================================================================

def plot_population_scaling(df):
    """Line plot showing scaling with number of subjects."""
    print("\nGenerating: Population Scaling Analysis...")

    pop_df = df[df["category"] == "population"].copy()

    if pop_df.empty:
        print("  No population data available")
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    platforms = [p for p in COLORS.keys() if p in pop_df["platform"].values]

    for platform in platforms:
        platform_data = pop_df[pop_df["platform"] == platform].sort_values("n_subjects")

        if not platform_data.empty:
            ax.errorbar(
                platform_data["n_subjects"],
                platform_data["mean_ms"],
                yerr=platform_data["std_ms"],
                label=platform,
                color=COLORS[platform],
                marker="o",
                capsize=3,
                linewidth=2,
                markersize=6
            )

    ax.set_xlabel("Number of Subjects")
    ax.set_ylabel("Time (ms)")
    ax.set_title("Population Simulation Scaling")
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)

    # Add ideal scaling line (linear)
    if "NeoPKPD (Julia)" in platforms:
        julia_data = pop_df[pop_df["platform"] == "NeoPKPD (Julia)"].sort_values("n_subjects")
        if len(julia_data) >= 2:
            x_ref = julia_data["n_subjects"].values
            y_ref = julia_data["mean_ms"].values
            # Fit linear scaling from first point
            slope = y_ref[0] / x_ref[0]
            x_line = np.array([x_ref[0], x_ref[-1]])
            y_line = slope * x_line
            ax.plot(x_line, y_line, "--", color="gray", alpha=0.5,
                    label="Linear scaling")

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_population_scaling.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


# =============================================================================
# Figure 3: Overall Speedup Summary
# =============================================================================

def plot_speedup_summary(df):
    """Horizontal bar chart showing speedup factors."""
    print("\nGenerating: Speedup Summary...")

    # Calculate speedups relative to NeoPKPD Julia
    speedups = []

    categories = df["category"].unique()
    ref_platform = "NeoPKPD (Julia)"

    for category in categories:
        cat_df = df[df["category"] == category]

        for model in cat_df["model"].unique():
            ref_data = cat_df[
                (cat_df["platform"] == ref_platform) &
                (cat_df["model"] == model)
            ]

            if ref_data.empty:
                continue

            ref_time = ref_data["mean_ms"].values[0]

            for platform in cat_df["platform"].unique():
                if platform == ref_platform:
                    continue

                other_data = cat_df[
                    (cat_df["platform"] == platform) &
                    (cat_df["model"] == model)
                ]

                if not other_data.empty:
                    other_time = other_data["mean_ms"].values[0]
                    if ref_time > 0:
                        speedup = other_time / ref_time
                        speedups.append({
                            "category": category,
                            "model": model,
                            "platform": platform,
                            "speedup": speedup
                        })

    if not speedups:
        print("  No speedup data available")
        return

    speedup_df = pd.DataFrame(speedups)

    # Aggregate by platform
    agg_speedup = speedup_df.groupby("platform")["speedup"].agg(["mean", "std"]).reset_index()
    agg_speedup = agg_speedup.sort_values("mean", ascending=True)

    fig, ax = plt.subplots(figsize=(8, 5))

    y_pos = np.arange(len(agg_speedup))
    colors = [COLORS.get(p, "#888888") for p in agg_speedup["platform"]]

    bars = ax.barh(y_pos, agg_speedup["mean"], xerr=agg_speedup["std"],
                   color=colors, capsize=5, alpha=0.9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(agg_speedup["platform"])
    ax.set_xlabel("Speedup Factor (vs NeoPKPD Julia)")
    ax.set_title("Average Speedup Comparison")
    ax.axvline(x=1, color="black", linestyle="--", alpha=0.5, label="NeoPKPD Julia baseline")

    # Add value labels
    for i, (mean, std) in enumerate(zip(agg_speedup["mean"], agg_speedup["std"])):
        ax.text(mean + std + 0.1, i, f"{mean:.2f}x", va="center", fontsize=9)

    ax.set_xlim(0, max(agg_speedup["mean"] + agg_speedup["std"]) * 1.2)
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_speedup_summary.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


# =============================================================================
# Generate LaTeX Table
# =============================================================================

def generate_latex_table(df):
    """Generate LaTeX-formatted summary table for the paper."""
    print("\nGenerating: LaTeX Summary Table...")

    # Create summary table
    summary_rows = []

    # Filter for key benchmarks
    key_benchmarks = [
        ("simulation", "single_simulation", "OneCompIV"),
        ("simulation", "single_simulation", "TwoCompIV"),
        ("simulation", "single_simulation", "MichaelisMenten"),
        ("population", "population_simulation", "TwoCompIV"),
        ("pkpd", "pkpd_simulation", "DirectEmax"),
        ("pkpd", "pkpd_simulation", "IndirectResponse"),
    ]

    for category, name, model in key_benchmarks:
        row_data = {"Task": f"{category.title()} ({model})"}

        for platform in COLORS.keys():
            data = df[
                (df["category"] == category) &
                (df["name"] == name) &
                (df["model"] == model) &
                (df["platform"] == platform)
            ]

            if not data.empty:
                mean = data["mean_ms"].values[0]
                std = data["std_ms"].values[0]
                row_data[platform] = f"{mean:.2f} ± {std:.2f}"
            else:
                row_data[platform] = "—"

        summary_rows.append(row_data)

    summary_df = pd.DataFrame(summary_rows)

    # Save as CSV
    output_path = FIGURES_DIR / "benchmark_summary_table.csv"
    summary_df.to_csv(output_path, index=False)
    print(f"  Saved: {output_path}")

    # Generate LaTeX
    latex_output = FIGURES_DIR / "benchmark_summary_table.tex"
    with open(latex_output, "w") as f:
        f.write("% Auto-generated benchmark table\n")
        f.write("% Include in paper with: \\input{figures/benchmark_summary_table.tex}\n\n")
        f.write("\\begin{table}[!ht]\n")
        f.write("\\centering\n")
        f.write("\\footnotesize\n")
        f.write("\\caption{Benchmark results (mean ± SD in milliseconds, N=100 runs)}\n")
        f.write("\\label{tab:benchmarks}\n")

        # Determine columns
        cols = ["Task"] + [p for p in COLORS.keys() if p in summary_df.columns]
        col_spec = "@{}l" + "c" * (len(cols) - 1) + "@{}"

        f.write(f"\\begin{{tabular}}{{{col_spec}}}\n")
        f.write("\\toprule\n")

        # Header
        header = " & ".join([f"\\textbf{{{c}}}" for c in cols])
        f.write(f"{header} \\\\\n")
        f.write("\\midrule\n")

        # Data rows
        for _, row in summary_df.iterrows():
            values = [str(row.get(c, "—")) for c in cols]
            f.write(" & ".join(values) + " \\\\\n")

        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    print(f"  Saved: {latex_output}")


# =============================================================================
# Main
# =============================================================================

def main():
    print("\n" + "=" * 70)
    print("NeoPKPD Benchmark Plotting")
    print("=" * 70)

    # Create output directory
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    df = load_benchmark_data()

    # Generate figures
    plot_simulation_comparison(df)
    plot_population_scaling(df)
    plot_speedup_summary(df)
    generate_latex_table(df)

    print("\n" + "=" * 70)
    print("PLOTTING COMPLETE")
    print(f"Figures saved to: {FIGURES_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
