"""
NeoPKPD PK Visualization

Pharmacokinetic visualization functions including:
- Concentration-time plots
- Spaghetti plots
- Mean with confidence ribbon plots
- Individual fit plots
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

from .backends import _get_plotter
from .themes import get_theme_config, get_color, get_colors


# ============================================================================
# Concentration-Time Plots
# ============================================================================

def plot_conc_time(
    result: Dict[str, Any],
    observation: str = "conc",
    log_scale: bool = False,
    show_markers: bool = True,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (10, 6),
    color: Optional[str] = None,
    label: Optional[str] = None,
) -> Any:
    """
    Plot concentration-time profile from simulation result.

    Args:
        result: Simulation result dictionary
        observation: Observation key to plot (default: "conc")
        log_scale: Use log scale for y-axis
        show_markers: Show data point markers
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size (width, height) in inches
        color: Line color (default: theme primary)
        label: Legend label

    Returns:
        Figure object (matplotlib Figure or plotly Figure)

    Example:
        >>> result = neopkpd.simulate_pk_iv_bolus(...)
        >>> fig = plot_conc_time(result, log_scale=True)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    t = result["t"]
    c = result["observations"][observation]

    if color is None:
        color = get_color(0)

    fig = plotter.create_figure(figsize=figsize, title=title)

    if show_markers:
        plotter.scatter_plot(fig, t, c, color=color, size=theme["marker_size"] * 5)

    plotter.line_plot(fig, t, c, label=label, color=color,
                      linewidth=theme["line_width"])

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    if label:
        plotter.add_legend(fig)

    return plotter.finalize(fig)


def plot_multi_conc_time(
    results: List[Dict[str, Any]],
    labels: Optional[List[str]] = None,
    observation: str = "conc",
    log_scale: bool = False,
    show_markers: bool = True,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (10, 6),
    colors: Optional[List[str]] = None,
) -> Any:
    """
    Plot multiple concentration-time profiles overlaid.

    Args:
        results: List of simulation result dictionaries
        labels: Labels for each result
        observation: Observation key to plot
        log_scale: Use log scale for y-axis
        show_markers: Show data point markers
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size
        colors: List of colors for each series

    Returns:
        Figure object

    Example:
        >>> results = [result_100mg, result_200mg, result_400mg]
        >>> labels = ["100 mg", "200 mg", "400 mg"]
        >>> fig = plot_multi_conc_time(results, labels=labels)
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    n = len(results)
    if colors is None:
        colors = get_colors(n)
    if labels is None:
        labels = [f"Series {i+1}" for i in range(n)]

    fig = plotter.create_figure(figsize=figsize, title=title)

    for i, (result, label, color) in enumerate(zip(results, labels, colors)):
        t = result["t"]
        c = result["observations"][observation]

        if show_markers:
            plotter.scatter_plot(fig, t, c, color=color, size=theme["marker_size"] * 5)

        plotter.line_plot(fig, t, c, label=label, color=color,
                          linewidth=theme["line_width"])

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig)

    return plotter.finalize(fig)


# ============================================================================
# Spaghetti Plots
# ============================================================================

def plot_spaghetti(
    population_result: Dict[str, Any],
    observation: str = "conc",
    n_subjects: Optional[int] = None,
    log_scale: bool = False,
    alpha: Optional[float] = None,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (10, 6),
    color: Optional[str] = None,
    show_mean: bool = True,
    mean_color: Optional[str] = None,
) -> Any:
    """
    Create spaghetti plot showing individual profiles from population simulation.

    Args:
        population_result: Population simulation result dictionary
        observation: Observation key to plot
        n_subjects: Number of subjects to plot (None = all)
        log_scale: Use log scale for y-axis
        alpha: Transparency for individual lines
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size
        color: Color for individual lines
        show_mean: Show mean line
        mean_color: Color for mean line

    Returns:
        Figure object

    Example:
        >>> pop_result = neopkpd.simulate_population_iv_bolus(...)
        >>> fig = plot_spaghetti(pop_result, n_subjects=50, show_mean=True)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    individuals = population_result["individuals"]
    if n_subjects is not None:
        individuals = individuals[:n_subjects]

    if alpha is None:
        alpha = theme["alpha_spaghetti"]
    if color is None:
        color = get_color(0)
    if mean_color is None:
        mean_color = get_color(1)

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Plot individual profiles
    for i, ind in enumerate(individuals):
        t = ind["t"]
        c = ind["observations"][observation]
        label = None if i > 0 else "Individual"
        plotter.line_plot(fig, t, c, label=label, color=color,
                          linewidth=0.8, alpha=alpha)

    # Plot mean if available in summaries
    if show_mean and "summaries" in population_result:
        summary_key = f"pk_{observation}" if f"pk_{observation}" in population_result["summaries"] else observation
        if summary_key in population_result["summaries"]:
            summary = population_result["summaries"][summary_key]
            t_mean = individuals[0]["t"]  # Assume same time points
            c_mean = summary["mean"]
            plotter.line_plot(fig, t_mean, c_mean, label="Mean", color=mean_color,
                              linewidth=theme["line_width"] * 1.5)

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig)

    return plotter.finalize(fig)


# ============================================================================
# Mean + Ribbon Plots
# ============================================================================

def plot_mean_ribbon(
    population_result: Dict[str, Any],
    observation: str = "conc",
    ci_levels: List[float] = [0.05, 0.95],
    show_median: bool = True,
    log_scale: bool = False,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (10, 6),
    mean_color: Optional[str] = None,
    ribbon_color: Optional[str] = None,
) -> Any:
    """
    Plot mean profile with confidence interval ribbon.

    Args:
        population_result: Population simulation result dictionary
        observation: Observation key to plot
        ci_levels: Quantile levels for ribbon [lower, upper]
        show_median: Show median line
        log_scale: Use log scale for y-axis
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size
        mean_color: Color for mean line
        ribbon_color: Color for confidence ribbon

    Returns:
        Figure object

    Example:
        >>> pop_result = neopkpd.simulate_population_iv_bolus(...)
        >>> fig = plot_mean_ribbon(pop_result, ci_levels=[0.025, 0.975])
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    if mean_color is None:
        mean_color = get_color(0)
    if ribbon_color is None:
        ribbon_color = get_color(0)

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Get summary data
    individuals = population_result["individuals"]
    t = individuals[0]["t"]

    # Calculate statistics manually if summaries not available
    import numpy as np

    n_time = len(t)
    all_conc = []
    for ind in individuals:
        all_conc.append(ind["observations"][observation])

    all_conc = np.array(all_conc)
    mean_c = np.mean(all_conc, axis=0)
    median_c = np.median(all_conc, axis=0)
    lower_c = np.percentile(all_conc, ci_levels[0] * 100, axis=0)
    upper_c = np.percentile(all_conc, ci_levels[1] * 100, axis=0)

    # Plot ribbon
    ci_pct = int((ci_levels[1] - ci_levels[0]) * 100)
    plotter.fill_between(fig, list(t), list(lower_c), list(upper_c),
                         color=ribbon_color, alpha=theme["alpha_ribbon"],
                         label=f"{ci_pct}% CI")

    # Plot mean
    plotter.line_plot(fig, list(t), list(mean_c), label="Mean",
                      color=mean_color, linewidth=theme["line_width"])

    # Plot median if requested
    if show_median:
        plotter.line_plot(fig, list(t), list(median_c), label="Median",
                          color=mean_color, linestyle="--",
                          linewidth=theme["line_width"])

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig)

    return plotter.finalize(fig)


# ============================================================================
# Individual Fits
# ============================================================================

def plot_individual_fits(
    population_result: Dict[str, Any],
    observed_data: Optional[List[Dict[str, List[float]]]] = None,
    observation: str = "conc",
    n_subjects: Optional[int] = None,
    n_cols: int = 4,
    subplot_size: Tuple[float, float] = (3, 2.5),
    log_scale: bool = False,
    title: Optional[str] = None,
) -> Any:
    """
    Create panel of individual fit plots.

    Args:
        population_result: Population simulation result dictionary
        observed_data: Optional list of observed data per subject
        observation: Observation key to plot
        n_subjects: Number of subjects to plot
        n_cols: Number of columns in panel
        subplot_size: Size of each subplot
        log_scale: Use log scale for y-axis
        title: Overall title

    Returns:
        Figure object

    Example:
        >>> fig = plot_individual_fits(pop_result, n_subjects=16, n_cols=4)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    individuals = population_result["individuals"]
    if n_subjects is not None:
        individuals = individuals[:n_subjects]

    n = len(individuals)
    n_rows = (n + n_cols - 1) // n_cols
    figsize = (subplot_size[0] * n_cols, subplot_size[1] * n_rows)

    # For matplotlib, create subplots
    backend = plotter.__class__.__name__

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        axes = axes.flatten() if n > 1 else [axes]

        for i, ind in enumerate(individuals):
            ax = axes[i]
            t = ind["t"]
            c = ind["observations"][observation]

            ax.plot(t, c, color=get_color(0), linewidth=theme["line_width"],
                    label="Predicted")

            if observed_data is not None and i < len(observed_data):
                obs = observed_data[i]
                ax.scatter(obs["t"], obs["c"], color=get_color(1),
                           s=theme["marker_size"] * 10, label="Observed")

            ax.set_title(f"Subject {i+1}")
            ax.set_xlabel("Time")
            ax.set_ylabel("Conc")

            if log_scale:
                ax.set_yscale("log")

        # Hide empty subplots
        for i in range(n, len(axes)):
            axes[i].set_visible(False)

        if title:
            fig.suptitle(title)

        fig.tight_layout()
        return fig

    else:
        # Plotly - use subplots
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(rows=n_rows, cols=n_cols,
                            subplot_titles=[f"Subject {i+1}" for i in range(n)])

        for i, ind in enumerate(individuals):
            row = i // n_cols + 1
            col = i % n_cols + 1

            t = ind["t"]
            c = ind["observations"][observation]

            fig.add_trace(
                go.Scatter(x=t, y=c, mode="lines", name="Predicted",
                           line=dict(color=get_color(0)), showlegend=(i == 0)),
                row=row, col=col
            )

            if observed_data is not None and i < len(observed_data):
                obs = observed_data[i]
                fig.add_trace(
                    go.Scatter(x=obs["t"], y=obs["c"], mode="markers",
                               name="Observed", marker=dict(color=get_color(1)),
                               showlegend=(i == 0)),
                    row=row, col=col
                )

        fig.update_layout(height=figsize[1] * 100, width=figsize[0] * 100,
                          title=title)

        if log_scale:
            fig.update_yaxes(type="log")

        return fig
