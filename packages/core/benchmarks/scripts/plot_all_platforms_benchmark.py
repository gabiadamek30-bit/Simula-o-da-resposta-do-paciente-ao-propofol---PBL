#!/usr/bin/env python3
"""
NeoPKPD Comprehensive 6-Platform Benchmark Plotting
====================================================

Generates publication-quality figures comparing:
- NeoPKPD (Julia) - Open source
- mrgsolve (R) - Open source
- nlmixr2 (R) - Open source
- NONMEM (Fortran) - Commercial [Reference data]
- Monolix (C++) - Commercial [Reference data]
- Pumas (Julia) - Commercial [Reference data]

Data sources:
- NeoPKPD, mrgsolve, nlmixr2: Measured benchmarks
- NONMEM, Monolix, Pumas: Literature reference data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path
import sys

# =============================================================================
# Configuration
# =============================================================================

RESULTS_DIR = Path(__file__).parent.parent / "results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"

# Color palette - distinct colors for 6 platforms
COLORS = {
    "NeoPKPD (Julia)": "#2563EB",   # Blue
    "Pumas": "#8B5CF6",              # Purple
    "mrgsolve": "#D97706",           # Orange
    "Monolix": "#059669",            # Green
    "nlmixr2": "#DC2626",            # Red
    "NONMEM": "#6B7280",             # Gray
}

# Platform categories
OPEN_SOURCE = ["NeoPKPD (Julia)", "mrgsolve", "nlmixr2"]
COMMERCIAL = ["Pumas", "Monolix", "NONMEM"]
JULIA_BASED = ["NeoPKPD (Julia)", "Pumas"]
R_BASED = ["mrgsolve", "nlmixr2"]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# =============================================================================
# Data Loading
# =============================================================================

def load_all_data():
    """Load benchmark data from all 6 platforms."""
    data = {}

    files = {
        "NeoPKPD (Julia)": "neopkpd_comprehensive_benchmarks.csv",
        "mrgsolve": "mrgsolve_comprehensive_benchmarks.csv",
        "nlmixr2": "nlmixr2_comprehensive_benchmarks.csv",
        "NONMEM": "nonmem_comprehensive_benchmarks.csv",
        "Monolix": "monolix_comprehensive_benchmarks.csv",
        "Pumas": "pumas_comprehensive_benchmarks.csv",
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
# Figure 1: PK Models - All 6 Platforms
# =============================================================================

def plot_pk_comparison_all(df):
    """Bar chart comparing PK model simulation times across all platforms."""
    print("\nGenerating: PK Models Comparison (6 Platforms)...")

    pk_df = df[df["category"] == "pk"].copy()

    if pk_df.empty:
        print("  No PK data available")
        return

    models = ["OneCompIVBolus", "TwoCompIVBolus", "TwoCompOral",
              "ThreeCompIVBolus", "MichaelisMenten", "TransitAbsorption"]
    models = [m for m in models if m in pk_df["model"].values]

    platforms = [p for p in COLORS.keys() if p in pk_df["platform"].values]

    fig, ax = plt.subplots(figsize=(14, 7))

    x = np.arange(len(models))
    width = 0.12
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
        hatch = '//' if platform in COMMERCIAL else None
        bars = ax.bar(x + offset, means, width, yerr=errors,
                      label=platform, color=COLORS[platform],
                      capsize=2, alpha=0.85, hatch=hatch, edgecolor='white')

    ax.set_xlabel("Model Type")
    ax.set_ylabel("Time (ms) - Log Scale")
    ax.set_title("PK Model Simulation Performance: 6 Platform Comparison")
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace("Comp", "-Cmp\n") for m in models], fontsize=8)

    # Custom legend with open source / commercial distinction
    handles = [Patch(facecolor=COLORS[p], label=p,
                    hatch='//' if p in COMMERCIAL else None,
                    edgecolor='white' if p in COMMERCIAL else COLORS[p])
               for p in platforms]
    ax.legend(handles=handles, loc="upper left", ncol=2)

    ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.3)

    # Add annotation for commercial vs open source
    ax.text(0.98, 0.02, "Hatched = Commercial (Reference Data)",
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=8, style='italic', color='gray')

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_pk_comparison_6platforms.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Figure 2: Category Summary - All Platforms
# =============================================================================

def plot_category_summary_all(df):
    """Bar chart summarizing performance by category for all platforms."""
    print("\nGenerating: Category Summary (6 Platforms)...")

    summary = df.groupby(["category", "platform"])["mean_ms"].mean().reset_index()

    categories = ["pk", "pd", "pkpd", "nca", "population", "sensitivity", "vpc", "trial", "estimation"]
    categories = [c for c in categories if c in summary["category"].values]

    platforms = [p for p in COLORS.keys() if p in summary["platform"].values]

    fig, ax = plt.subplots(figsize=(14, 7))

    x = np.arange(len(categories))
    width = 0.12
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
        hatch = '//' if platform in COMMERCIAL else None
        ax.bar(x + offset, means, width,
               label=platform, color=COLORS[platform],
               alpha=0.85, hatch=hatch, edgecolor='white')

    ax.set_xlabel("Benchmark Category")
    ax.set_ylabel("Mean Time (ms) - Log Scale")
    ax.set_title("Average Performance by Category: 6 Platform Comparison")
    ax.set_xticks(x)
    ax.set_xticklabels([c.upper() for c in categories])

    handles = [Patch(facecolor=COLORS[p], label=p,
                    hatch='//' if p in COMMERCIAL else None,
                    edgecolor='white' if p in COMMERCIAL else COLORS[p])
               for p in platforms]
    ax.legend(handles=handles, loc="upper left", ncol=2)

    ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.3)

    ax.text(0.98, 0.02, "Hatched = Commercial (Reference Data)",
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=8, style='italic', color='gray')

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_category_summary_6platforms.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Figure 3: Population Scaling - All Platforms
# =============================================================================

def plot_population_scaling_all(df):
    """Line plot showing population simulation scaling across platforms."""
    print("\nGenerating: Population Scaling (6 Platforms)...")

    pop_df = df[df["category"] == "population"].copy()

    if pop_df.empty:
        print("  No population data available")
        return

    fig, ax = plt.subplots(figsize=(12, 7))

    platforms = [p for p in COLORS.keys() if p in pop_df["platform"].values]

    for platform in platforms:
        platform_data = pop_df[pop_df["platform"] == platform].sort_values("n_subjects")

        if not platform_data.empty and platform_data["n_subjects"].nunique() > 1:
            linestyle = '--' if platform in COMMERCIAL else '-'
            marker = 's' if platform in COMMERCIAL else 'o'

            ax.errorbar(
                platform_data["n_subjects"],
                platform_data["mean_ms"],
                yerr=platform_data["std_ms"],
                label=platform,
                color=COLORS[platform],
                marker=marker,
                linestyle=linestyle,
                capsize=3,
                linewidth=2,
                markersize=6
            )

    ax.set_xlabel("Number of Subjects")
    ax.set_ylabel("Time (ms) - Log Scale")
    ax.set_title("Population Simulation Scaling: 6 Platform Comparison")
    ax.legend(loc="upper left")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)

    ax.text(0.98, 0.02, "Dashed = Commercial (Reference Data)",
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=8, style='italic', color='gray')

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_population_scaling_6platforms.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Figure 4: Speedup Heatmap
# =============================================================================

def plot_speedup_heatmap(df):
    """Heatmap showing speedup of NeoPKPD vs all other platforms."""
    print("\nGenerating: Speedup Heatmap...")

    ref_platform = "NeoPKPD (Julia)"

    # Calculate mean time per category per platform
    summary = df.groupby(["category", "platform"])["mean_ms"].mean().reset_index()

    categories = ["pk", "pd", "pkpd", "nca", "population", "sensitivity", "vpc", "trial", "estimation"]
    categories = [c for c in categories if c in summary["category"].values]

    other_platforms = [p for p in COLORS.keys() if p != ref_platform and p in summary["platform"].values]

    # Build speedup matrix
    speedup_matrix = np.zeros((len(categories), len(other_platforms)))

    for i, cat in enumerate(categories):
        ref_data = summary[(summary["category"] == cat) & (summary["platform"] == ref_platform)]
        if ref_data.empty:
            continue
        ref_time = ref_data["mean_ms"].values[0]

        for j, platform in enumerate(other_platforms):
            other_data = summary[(summary["category"] == cat) & (summary["platform"] == platform)]
            if not other_data.empty and ref_time > 0:
                speedup_matrix[i, j] = other_data["mean_ms"].values[0] / ref_time

    fig, ax = plt.subplots(figsize=(12, 8))

    # Use log scale for colors
    im = ax.imshow(np.log10(speedup_matrix + 0.01), cmap='RdYlGn_r', aspect='auto')

    # Add text annotations
    for i in range(len(categories)):
        for j in range(len(other_platforms)):
            value = speedup_matrix[i, j]
            if value >= 1:
                text = f"{value:.0f}x" if value >= 10 else f"{value:.1f}x"
                color = 'white' if value > 10 else 'black'
            else:
                text = f"{value:.2f}x"
                color = 'black'
            ax.text(j, i, text, ha='center', va='center', color=color, fontsize=9)

    ax.set_xticks(np.arange(len(other_platforms)))
    ax.set_yticks(np.arange(len(categories)))
    ax.set_xticklabels(other_platforms, rotation=45, ha='right')
    ax.set_yticklabels([c.upper() for c in categories])

    ax.set_title("Performance Ratio vs NeoPKPD (>1 = slower than NeoPKPD)")
    ax.set_xlabel("Platform")
    ax.set_ylabel("Category")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Log10(Speedup Ratio)')

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_speedup_heatmap.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Figure 5: Open Source vs Commercial Comparison
# =============================================================================

def plot_opensource_vs_commercial(df):
    """Compare open source vs commercial platforms."""
    print("\nGenerating: Open Source vs Commercial Comparison...")

    summary = df.groupby(["category", "platform"])["mean_ms"].mean().reset_index()

    categories = ["pk", "pd", "pkpd", "population", "vpc", "trial", "estimation"]
    categories = [c for c in categories if c in summary["category"].values]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Open Source subplot
    ax1 = axes[0]
    open_source_platforms = [p for p in OPEN_SOURCE if p in summary["platform"].values]

    x = np.arange(len(categories))
    width = 0.25

    for i, platform in enumerate(open_source_platforms):
        platform_data = summary[summary["platform"] == platform]
        means = []
        for cat in categories:
            cat_data = platform_data[platform_data["category"] == cat]
            means.append(cat_data["mean_ms"].values[0] if not cat_data.empty else 0)

        offset = (i - len(open_source_platforms)/2 + 0.5) * width
        ax1.bar(x + offset, means, width, label=platform, color=COLORS[platform], alpha=0.9)

    ax1.set_xlabel("Category")
    ax1.set_ylabel("Mean Time (ms)")
    ax1.set_title("Open Source Platforms")
    ax1.set_xticks(x)
    ax1.set_xticklabels([c.upper() for c in categories])
    ax1.legend()
    ax1.set_yscale("log")
    ax1.grid(axis="y", alpha=0.3)

    # Commercial subplot
    ax2 = axes[1]
    commercial_platforms = [p for p in COMMERCIAL if p in summary["platform"].values]

    for i, platform in enumerate(commercial_platforms):
        platform_data = summary[summary["platform"] == platform]
        means = []
        for cat in categories:
            cat_data = platform_data[platform_data["category"] == cat]
            means.append(cat_data["mean_ms"].values[0] if not cat_data.empty else 0)

        offset = (i - len(commercial_platforms)/2 + 0.5) * width
        ax2.bar(x + offset, means, width, label=platform, color=COLORS[platform], alpha=0.9, hatch='//')

    ax2.set_xlabel("Category")
    ax2.set_ylabel("Mean Time (ms)")
    ax2.set_title("Commercial Platforms (Reference Data)")
    ax2.set_xticks(x)
    ax2.set_xticklabels([c.upper() for c in categories])
    ax2.legend()
    ax2.set_yscale("log")
    ax2.grid(axis="y", alpha=0.3)

    plt.suptitle("Pharmacometrics Platform Performance Comparison", fontsize=14, y=1.02)
    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_opensource_vs_commercial.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Figure 6: Julia vs R vs Fortran
# =============================================================================

def plot_language_comparison(df):
    """Compare platforms by implementation language."""
    print("\nGenerating: Language Comparison...")

    summary = df.groupby(["category", "platform"])["mean_ms"].mean().reset_index()

    # Group by language
    language_groups = {
        "Julia": ["NeoPKPD (Julia)", "Pumas"],
        "R/C++": ["mrgsolve", "nlmixr2"],
        "Fortran": ["NONMEM"],
        "C++": ["Monolix"]
    }

    categories = ["pk", "pkpd", "population", "vpc", "trial"]
    categories = [c for c in categories if c in summary["category"].values]

    fig, ax = plt.subplots(figsize=(12, 7))

    # Calculate mean per language per category
    lang_data = {}
    for lang, platforms in language_groups.items():
        for cat in categories:
            cat_summary = summary[(summary["category"] == cat) & (summary["platform"].isin(platforms))]
            if not cat_summary.empty:
                if lang not in lang_data:
                    lang_data[lang] = {}
                lang_data[lang][cat] = cat_summary["mean_ms"].mean()

    # Plot
    x = np.arange(len(categories))
    width = 0.2
    lang_colors = {"Julia": "#2563EB", "R/C++": "#D97706", "Fortran": "#6B7280", "C++": "#059669"}

    for i, (lang, cat_times) in enumerate(lang_data.items()):
        means = [cat_times.get(cat, 0) for cat in categories]
        offset = (i - len(lang_data)/2 + 0.5) * width
        ax.bar(x + offset, means, width, label=lang, color=lang_colors.get(lang, '#888'), alpha=0.9)

    ax.set_xlabel("Category")
    ax.set_ylabel("Mean Time (ms) - Log Scale")
    ax.set_title("Performance by Implementation Language")
    ax.set_xticks(x)
    ax.set_xticklabels([c.upper() for c in categories])
    ax.legend()
    ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    output_path = FIGURES_DIR / "figure_language_comparison.png"
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")

# =============================================================================
# Generate Summary Tables
# =============================================================================

def generate_summary_tables(df):
    """Generate comprehensive summary tables."""
    print("\nGenerating: Summary Tables...")

    # Full summary
    summary = df.groupby(["category", "platform"]).agg({
        "mean_ms": "mean",
        "std_ms": "mean",
        "model": "count"
    }).reset_index()
    summary.columns = ["Category", "Platform", "Mean (ms)", "Std (ms)", "N Benchmarks"]

    # Calculate speedup vs NeoPKPD
    ref = summary[summary["Platform"] == "NeoPKPD (Julia)"][["Category", "Mean (ms)"]].copy()
    ref.columns = ["Category", "ref_ms"]

    summary = summary.merge(ref, on="Category", how="left")
    summary["Speedup"] = summary["ref_ms"] / summary["Mean (ms)"]
    summary = summary.drop("ref_ms", axis=1)

    # Save CSV
    output_csv = FIGURES_DIR / "comprehensive_benchmark_summary_6platforms.csv"
    summary.to_csv(output_csv, index=False)
    print(f"  Saved: {output_csv}")

    # Generate LaTeX table
    latex_output = FIGURES_DIR / "comprehensive_benchmark_table_6platforms.tex"
    with open(latex_output, "w") as f:
        f.write("% Comprehensive NeoPKPD 6-Platform Benchmark Results\n")
        f.write("% Commercial platforms marked with * use reference data from literature\n")
        f.write("\\begin{table}[!ht]\n")
        f.write("\\centering\n")
        f.write("\\small\n")
        f.write("\\caption{Comprehensive benchmark comparison across 6 pharmacometrics platforms}\n")
        f.write("\\label{tab:6platform_benchmarks}\n")
        f.write("\\begin{tabular}{@{}llrrr@{}}\n")
        f.write("\\toprule\n")
        f.write("\\textbf{Category} & \\textbf{Platform} & \\textbf{Mean (ms)} & \\textbf{Std} & \\textbf{vs NeoPKPD} \\\\\n")
        f.write("\\midrule\n")

        for cat in summary["Category"].unique():
            cat_data = summary[summary["Category"] == cat]
            for idx, row in cat_data.iterrows():
                speedup_str = f"{row['Speedup']:.2f}x" if pd.notna(row['Speedup']) else "—"
                platform_name = row['Platform']
                if platform_name in COMMERCIAL:
                    platform_name += "*"
                f.write(f"{row['Category'] if idx == cat_data.index[0] else ''} & "
                        f"{platform_name} & {row['Mean (ms)']:.2f} & "
                        f"{row['Std (ms)']:.2f} & {speedup_str} \\\\\n")
            f.write("\\midrule\n")

        f.write("\\bottomrule\n")
        f.write("\\multicolumn{5}{l}{\\footnotesize *Commercial platforms use reference data from published literature}\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    print(f"  Saved: {latex_output}")

# =============================================================================
# Main
# =============================================================================

def main():
    print("\n" + "=" * 70)
    print("NeoPKPD 6-Platform Comprehensive Benchmark Plotting")
    print("=" * 70)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    df = load_all_data()

    plot_pk_comparison_all(df)
    plot_category_summary_all(df)
    plot_population_scaling_all(df)
    plot_speedup_heatmap(df)
    plot_opensource_vs_commercial(df)
    plot_language_comparison(df)
    generate_summary_tables(df)

    print("\n" + "=" * 70)
    print("PLOTTING COMPLETE")
    print(f"Figures saved to: {FIGURES_DIR}")
    print("=" * 70)

if __name__ == "__main__":
    main()
