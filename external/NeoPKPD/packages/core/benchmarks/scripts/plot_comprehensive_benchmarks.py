#!/usr/bin/env python3
"""
NeoPKPD Comprehensive Benchmark Plotting
=========================================

Generates publication-quality figures comparing NeoPKPD vs mrgsolve vs nlmixr2
across ALL benchmark categories.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# =============================================================================
# Configuration
# =============================================================================

RESULTS_DIR = Path(__file__).parent.parent / "results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"

# Color palette
COLORS = {
    "NeoPKPD (Julia)": "#2563EB",
    "mrgsolve": "#D97706",
    "nlmixr2": "#059669",
}

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

def load_data():
    """Load all comprehensive benchmark CSV files."""
    data = {}

    files = {
        "NeoPKPD (Julia)": "neopkpd_comprehensive_benchmarks.csv",
        "mrgsolve": "mrgsolve_comprehensive_benchmarks.csv",
        "nlmixr2": "nlmixr2_comprehensive_benchmarks.csv",
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
# Figure 1: PK Models Comparison
# =============================================================================

def plot_pk_comparison(df):
    """Bar chart comparing PK model simulation times."""
    print("\nGenerating: PK Models Comparison...")

    pk_df = df[df["category"] == "pk"].copy()

    if pk_df.empty:
        print("  No PK data available")
        return

    # Models to compare
    models = ["OneCompIVBolus", "TwoCompIVBolus", "TwoCompOral",
              "ThreeCompIVBolus", "MichaelisMenten", "TransitAbsorption"]
    models = [m for m in models if m in pk_df["model"].values]

    platforms = [p for p in COLORS.keys() if p in pk_df["platform"].values]

    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(len(models))
    width = 0.25
    n_platforms = len(platforms)

    for i, platform in enumerate(platforms):
        platform_data = pk_df[pk_df["platform"] == platform]
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
    ax.set_title("PK Model Simulation Performance Comparison")
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace("Comp", "-Comp\n") for m in models], fontsize=8)
    ax.legend(loc="upper left")
    ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.3)

    # Add speedup annotations for NeoPKPD vs others
    if "NeoPKPD (Julia)" in platforms:
        for i, model in enumerate(models):
            neo_time = pk_df[
                (pk_df["platform"] == "NeoPKPD (Julia)") &
                (pk_df["model"] == model)
            ]["mean_ms"].values

            if len(neo_time) > 0 and neo_time[0] > 0:
                for platform in ["mrgsolve", "nlmixr2"]:
                    other_time = pk_df[
                        (pk_df["platform"] == platform) &
                        (pk_df["model"] == model)
                    ]["mean_ms"].values

                    if len(other_time) > 0:
                        speedup = other_time[0] / neo_time[0]
                        if speedup > 1.5:
                            ax.annotate(f"{speedup:.0f}x",
                                        xy=(i, other_time[0]),
                                        xytext=(0, 5),
                                        textcoords="offset points",
                                        ha="center", fontsize=7,
                                        color=COLORS[platform])

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_pk_comparison.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Figure 2: Population Scaling
# =============================================================================

def plot_population_scaling(df):
    """Line plot showing population simulation scaling."""
    print("\nGenerating: Population Scaling Analysis...")

    pop_df = df[df["category"] == "population"].copy()

    if pop_df.empty:
        print("  No population data available")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    platforms = [p for p in COLORS.keys() if p in pop_df["platform"].values]

    for platform in platforms:
        platform_data = pop_df[pop_df["platform"] == platform].sort_values("n_subjects")

        if not platform_data.empty and platform_data["n_subjects"].nunique() > 1:
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
    ax.set_title("Population Simulation Scaling Comparison")
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_population_scaling.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Figure 3: Category Summary
# =============================================================================

def plot_category_summary(df):
    """Bar chart summarizing performance by category."""
    print("\nGenerating: Category Summary...")

    # Calculate mean time per category per platform
    summary = df.groupby(["category", "platform"])["mean_ms"].mean().reset_index()

    categories = ["pk", "pd", "pkpd", "nca", "population", "sensitivity"]
    categories = [c for c in categories if c in summary["category"].values]

    platforms = [p for p in COLORS.keys() if p in summary["platform"].values]

    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(len(categories))
    width = 0.25
    n_platforms = len(platforms)

    for i, platform in enumerate(platforms):
        platform_data = summary[summary["platform"] == platform]
        means = []

        for cat in categories:
            cat_data = platform_data[platform_data["category"] == cat]
            if not cat_data.empty:
                means.append(cat_data["mean_ms"].values[0])
            else:
                means.append(0)

        offset = (i - n_platforms/2 + 0.5) * width
        ax.bar(x + offset, means, width,
               label=platform, color=COLORS[platform], alpha=0.9)

    ax.set_xlabel("Benchmark Category")
    ax.set_ylabel("Mean Time (ms)")
    ax.set_title("Average Performance by Category")
    ax.set_xticks(x)
    ax.set_xticklabels([c.upper() for c in categories])
    ax.legend()
    ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_category_summary.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Figure 4: Speedup Summary
# =============================================================================

def plot_speedup_summary(df):
    """Horizontal bar chart showing speedup factors."""
    print("\nGenerating: Speedup Summary...")

    speedups = []
    ref_platform = "NeoPKPD (Julia)"

    for category in df["category"].unique():
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

                if not other_data.empty and ref_time > 0:
                    other_time = other_data["mean_ms"].values[0]
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

    # Aggregate by category and platform
    agg = speedup_df.groupby(["category", "platform"])["speedup"].agg(["mean", "std"]).reset_index()

    # Create subplot for each category
    categories = agg["category"].unique()
    n_cats = len(categories)

    fig, axes = plt.subplots(1, n_cats, figsize=(4 * n_cats, 5), sharey=True)
    if n_cats == 1:
        axes = [axes]

    for idx, cat in enumerate(categories):
        ax = axes[idx]
        cat_data = agg[agg["category"] == cat].sort_values("mean")

        y_pos = np.arange(len(cat_data))
        colors = [COLORS.get(p, "#888888") for p in cat_data["platform"]]

        ax.barh(y_pos, cat_data["mean"], xerr=cat_data["std"],
                color=colors, capsize=5, alpha=0.9)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(cat_data["platform"])
        ax.set_xlabel("Speedup vs NeoPKPD")
        ax.set_title(cat.upper())
        ax.axvline(x=1, color="black", linestyle="--", alpha=0.5)
        ax.grid(axis="x", alpha=0.3)

        # Add value labels
        for i, (mean, std) in enumerate(zip(cat_data["mean"], cat_data["std"])):
            ax.text(mean + std + 0.1, i, f"{mean:.1f}x", va="center", fontsize=8)

    plt.suptitle("Performance Comparison (higher = slower than NeoPKPD)", fontsize=12)
    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_speedup_by_category.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Generate Summary Table
# =============================================================================

def generate_summary_table(df):
    """Generate comprehensive summary table."""
    print("\nGenerating: Summary Tables...")

    # Pivot table by category and platform
    summary = df.groupby(["category", "platform"]).agg({
        "mean_ms": "mean",
        "std_ms": "mean",
        "model": "count"
    }).reset_index()
    summary.columns = ["Category", "Platform", "Mean (ms)", "Std (ms)", "N Benchmarks"]

    # Calculate speedup
    ref = summary[summary["Platform"] == "NeoPKPD (Julia)"][["Category", "Mean (ms)"]].copy()
    ref.columns = ["Category", "ref_ms"]

    summary = summary.merge(ref, on="Category", how="left")
    summary["Speedup"] = summary["ref_ms"] / summary["Mean (ms)"]
    summary = summary.drop("ref_ms", axis=1)

    # Save CSV
    output_csv = FIGURES_DIR / "comprehensive_benchmark_summary.csv"
    summary.to_csv(output_csv, index=False)
    print(f"  Saved: {output_csv}")

    # Generate LaTeX table
    latex_output = FIGURES_DIR / "comprehensive_benchmark_table.tex"
    with open(latex_output, "w") as f:
        f.write("% Comprehensive NeoPKPD Benchmark Results\n")
        f.write("\\begin{table}[!ht]\n")
        f.write("\\centering\n")
        f.write("\\small\n")
        f.write("\\caption{Comprehensive benchmark comparison (mean time in ms)}\n")
        f.write("\\label{tab:comprehensive_benchmarks}\n")
        f.write("\\begin{tabular}{@{}llrrr@{}}\n")
        f.write("\\toprule\n")
        f.write("\\textbf{Category} & \\textbf{Platform} & \\textbf{Mean (ms)} & \\textbf{Std} & \\textbf{Speedup} \\\\\n")
        f.write("\\midrule\n")

        for cat in summary["Category"].unique():
            cat_data = summary[summary["Category"] == cat]
            for idx, row in cat_data.iterrows():
                speedup_str = f"{row['Speedup']:.2f}x" if pd.notna(row['Speedup']) else "—"
                f.write(f"{row['Category'] if idx == cat_data.index[0] else ''} & "
                        f"{row['Platform']} & {row['Mean (ms)']:.2f} & "
                        f"{row['Std (ms)']:.2f} & {speedup_str} \\\\\n")
            f.write("\\midrule\n")

        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    print(f"  Saved: {latex_output}")

# =============================================================================
# Main
# =============================================================================

def main():
    print("\n" + "=" * 70)
    print("NeoPKPD Comprehensive Benchmark Plotting")
    print("=" * 70)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    df = load_data()

    plot_pk_comparison(df)
    plot_population_scaling(df)
    plot_category_summary(df)
    plot_speedup_summary(df)
    generate_summary_table(df)

    print("\n" + "=" * 70)
    print("PLOTTING COMPLETE")
    print(f"Figures saved to: {FIGURES_DIR}")
    print("=" * 70)

if __name__ == "__main__":
    main()
