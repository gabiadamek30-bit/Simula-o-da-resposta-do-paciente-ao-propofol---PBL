"""
NeoPKPD VPC Visualization

Visual Predictive Check visualization functions including:
- Standard VPC plots with percentile ribbons
- Prediction-corrected VPC (pcVPC)
- Stratified VPC with faceting
- VPC with BLQ handling visualization
- Confidence interval bands
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np

from .backends import _get_plotter, get_backend
from .themes import get_theme_config, get_color, get_colors, NEOPKPD_COLORS


# ============================================================================
# Standard VPC Plot
# ============================================================================

def plot_vpc(
    vpc_result: Any,
    log_scale: bool = False,
    show_observed: bool = True,
    show_simulated_median: bool = True,
    show_ci: bool = True,
    pi_levels: Optional[List[float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (10, 6),
    obs_color: Optional[str] = None,
    sim_color: Optional[str] = None,
    ci_alpha: float = 0.2,
    median_alpha: float = 0.5,
) -> Any:
    """
    Plot standard Visual Predictive Check.

    Creates a VPC plot showing observed percentiles (lines) overlaid on
    simulated percentile confidence intervals (ribbons). This is the
    primary diagnostic for population PK model evaluation.

    Args:
        vpc_result: VPCResult from compute_vpc or compute_vpc_python
        log_scale: Use log scale for y-axis
        show_observed: Show observed percentile lines
        show_simulated_median: Show simulated median lines
        show_ci: Show confidence interval ribbons
        pi_levels: PI levels to plot (default: all from result)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size (width, height) in inches
        obs_color: Color for observed lines (default: theme colors)
        sim_color: Color for simulated elements (default: theme colors)
        ci_alpha: Transparency for CI ribbons
        median_alpha: Transparency for simulated median lines

    Returns:
        Figure object (matplotlib Figure or plotly Figure)

    Example:
        >>> from neopkpd import vpc
        >>> result = vpc.compute_vpc(observed, pop_spec, grid, config)
        >>> from neopkpd.viz import plot_vpc
        >>> fig = plot_vpc(result, log_scale=True)
        >>> fig.savefig("vpc_plot.png")
    """
    # Import vpc module for type access
    from ..vpc import (
        VPCResult, get_bin_midpoints, get_observed_percentile,
        get_simulated_median, get_simulated_ci
    )

    plotter = _get_plotter()
    theme = get_theme_config()

    # Get time points
    times = get_bin_midpoints(vpc_result)

    # Determine PI levels to plot
    if pi_levels is None:
        pi_levels = vpc_result.config.pi_levels

    # Define colors for percentiles (5th, 50th, 95th style)
    n_levels = len(pi_levels)
    if obs_color is None:
        obs_colors = get_colors(n_levels)
    else:
        obs_colors = [obs_color] * n_levels

    if sim_color is None:
        sim_colors = get_colors(n_levels)
    else:
        sim_colors = [sim_color] * n_levels

    # Create figure
    plot_title = title or ("pcVPC" if vpc_result.prediction_corrected else "Visual Predictive Check")
    fig = plotter.create_figure(figsize=figsize, title=plot_title)

    # Plot each percentile level
    for i, level in enumerate(pi_levels):
        pct_label = f"{int(level * 100)}th"

        # Get data
        obs_values = get_observed_percentile(vpc_result, level)
        sim_median = get_simulated_median(vpc_result, level)
        sim_lower, sim_upper = get_simulated_ci(vpc_result, level)

        # Plot CI ribbon
        if show_ci and len(sim_lower) > 0:
            plotter.fill_between(
                fig, list(times), list(sim_lower), list(sim_upper),
                color=sim_colors[i], alpha=ci_alpha,
                label=f"Sim {pct_label} CI"
            )

        # Plot simulated median line
        if show_simulated_median and len(sim_median) > 0:
            plotter.line_plot(
                fig, list(times), list(sim_median),
                color=sim_colors[i], linestyle="--",
                linewidth=theme["line_width"], alpha=median_alpha,
                label=f"Sim {pct_label}"
            )

        # Plot observed percentile line
        if show_observed and len(obs_values) > 0:
            plotter.line_plot(
                fig, list(times), list(obs_values),
                color=obs_colors[i], linestyle="-",
                linewidth=theme["line_width"] * 1.2,
                label=f"Obs {pct_label}"
            )

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=plot_title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig, location="best")

    return plotter.finalize(fig)


# ============================================================================
# Prediction-Corrected VPC Plot
# ============================================================================

def plot_pcvpc(
    vpc_result: Any,
    log_scale: bool = False,
    show_observed: bool = True,
    show_ci: bool = True,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Prediction-Corrected Concentration",
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Plot prediction-corrected Visual Predictive Check.

    pcVPC removes the structural model trend, making it easier to
    assess residual variability. This is particularly useful when
    structural model predictions vary substantially over time.

    Args:
        vpc_result: VPCResult from compute_pcvpc (prediction_corrected=True)
        log_scale: Use log scale for y-axis
        show_observed: Show observed percentile lines
        show_ci: Show confidence interval ribbons
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> config = vpc.VPCConfig(prediction_corrected=True)
        >>> result = vpc.compute_pcvpc(observed, pop_spec, grid, config)
        >>> fig = plot_pcvpc(result)
    """
    # Use plot_vpc with pcVPC-specific defaults
    plot_title = title or "Prediction-Corrected VPC"

    return plot_vpc(
        vpc_result,
        log_scale=log_scale,
        show_observed=show_observed,
        show_ci=show_ci,
        title=plot_title,
        xlabel=xlabel,
        ylabel=ylabel,
        figsize=figsize,
    )


# ============================================================================
# Stratified VPC Plot
# ============================================================================

def plot_stratified_vpc(
    stratified_result: Any,
    log_scale: bool = False,
    n_cols: int = 2,
    subplot_size: Tuple[float, float] = (5, 4),
    show_observed: bool = True,
    show_ci: bool = True,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
) -> Any:
    """
    Plot stratified VPC as faceted panels.

    Creates a multi-panel plot with separate VPCs for each stratum,
    allowing visual comparison of model performance across covariate groups.

    Args:
        stratified_result: StratifiedVPCResult from compute_stratified_vpc
        log_scale: Use log scale for y-axis
        n_cols: Number of columns in panel layout
        subplot_size: Size of each subplot (width, height)
        show_observed: Show observed percentile lines
        show_ci: Show confidence interval ribbons
        title: Overall title
        xlabel: X-axis label
        ylabel: Y-axis label

    Returns:
        Figure object

    Example:
        >>> config = vpc.VPCConfig(stratify_by=["WT_GROUP"])
        >>> result = vpc.compute_stratified_vpc(observed, pop_spec, grid, config)
        >>> fig = plot_stratified_vpc(result)
    """
    from ..vpc import (
        StratifiedVPCResult, get_bin_midpoints, get_observed_percentile,
        get_simulated_median, get_simulated_ci
    )

    plotter = _get_plotter()
    theme = get_theme_config()
    backend = plotter.__class__.__name__

    results = stratified_result.results
    strata_names = stratified_result.strata_names
    n_strata = len(results)
    n_rows = (n_strata + n_cols - 1) // n_cols
    figsize = (subplot_size[0] * n_cols, subplot_size[1] * n_rows)

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n_strata > 1:
            axes = axes.flatten()
        else:
            axes = [axes]

        colors = get_colors(3)  # For 5th, 50th, 95th percentiles

        for idx, (vpc_result, stratum_name) in enumerate(zip(results, strata_names)):
            ax = axes[idx]
            times = get_bin_midpoints(vpc_result)
            pi_levels = vpc_result.config.pi_levels

            for i, level in enumerate(pi_levels):
                obs_values = get_observed_percentile(vpc_result, level)
                sim_lower, sim_upper = get_simulated_ci(vpc_result, level)

                if show_ci and len(sim_lower) > 0:
                    ax.fill_between(times, sim_lower, sim_upper,
                                    color=colors[i], alpha=0.2)

                if show_observed and len(obs_values) > 0:
                    ax.plot(times, obs_values, color=colors[i],
                            linewidth=theme["line_width"],
                            label=f"{int(level * 100)}th")

            ax.set_title(stratum_name)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            if log_scale:
                ax.set_yscale("log")
            if idx == 0:
                ax.legend(fontsize="small")

        # Hide empty subplots
        for idx in range(n_strata, len(axes)):
            axes[idx].set_visible(False)

        if title:
            fig.suptitle(title)

        fig.tight_layout()
        return fig

    else:
        # Plotly implementation
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(
            rows=n_rows, cols=n_cols,
            subplot_titles=strata_names
        )

        colors = get_colors(3)

        for idx, vpc_result in enumerate(results):
            row = idx // n_cols + 1
            col = idx % n_cols + 1

            times = get_bin_midpoints(vpc_result)
            pi_levels = vpc_result.config.pi_levels

            for i, level in enumerate(pi_levels):
                obs_values = get_observed_percentile(vpc_result, level)
                sim_lower, sim_upper = get_simulated_ci(vpc_result, level)

                if show_ci and len(sim_lower) > 0:
                    # Upper bound
                    fig.add_trace(
                        go.Scatter(x=times, y=sim_upper, mode="lines",
                                   line=dict(width=0), showlegend=False),
                        row=row, col=col
                    )
                    # Lower with fill
                    fig.add_trace(
                        go.Scatter(x=times, y=sim_lower, mode="lines",
                                   line=dict(width=0), fill="tonexty",
                                   fillcolor=f"rgba(0,0,0,0.1)",
                                   showlegend=False),
                        row=row, col=col
                    )

                if show_observed and len(obs_values) > 0:
                    fig.add_trace(
                        go.Scatter(x=times, y=obs_values, mode="lines",
                                   name=f"{int(level * 100)}th",
                                   line=dict(color=colors[i]),
                                   showlegend=(idx == 0)),
                        row=row, col=col
                    )

        fig.update_layout(
            height=figsize[1] * 100,
            width=figsize[0] * 100,
            title=title
        )

        if log_scale:
            fig.update_yaxes(type="log")

        return fig


# ============================================================================
# VPC with BLQ Plot
# ============================================================================

def plot_vpc_with_blq(
    vpc_result: Any,
    log_scale: bool = False,
    blq_panel_height: float = 0.25,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    blq_ylabel: str = "% BLQ",
    figsize: Tuple[float, float] = (10, 8),
    obs_color: Optional[str] = None,
    blq_color: Optional[str] = None,
) -> Any:
    """
    Plot VPC with BLQ (below limit of quantitation) panel.

    Creates a two-panel plot:
    - Upper: Standard VPC with percentile ribbons
    - Lower: %BLQ comparison between observed and simulated

    This is essential for drugs with significant BLQ observations.

    Args:
        vpc_result: VPCResult with blq_stats from compute_vpc_with_blq
        log_scale: Use log scale for concentration y-axis
        blq_panel_height: Relative height of BLQ panel (0-1)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label for concentration
        blq_ylabel: Y-axis label for BLQ panel
        figsize: Figure size
        obs_color: Color for observed data
        blq_color: Color for BLQ elements

    Returns:
        Figure object

    Example:
        >>> config = vpc.VPCConfig(lloq=0.1, blq_method=vpc.BLQMethod.M4)
        >>> result = vpc.compute_vpc_with_blq(observed, pop_spec, grid, config)
        >>> fig = plot_vpc_with_blq(result)
    """
    from ..vpc import (
        get_bin_midpoints, get_observed_percentile, get_simulated_ci,
        get_blq_observed, get_blq_simulated
    )

    theme = get_theme_config()

    if obs_color is None:
        obs_color = get_color(0)
    if blq_color is None:
        blq_color = get_color(1)

    times = get_bin_midpoints(vpc_result)
    pi_levels = vpc_result.config.pi_levels

    backend = get_backend()

    if backend == "matplotlib":
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        fig = plt.figure(figsize=figsize)
        gs = GridSpec(2, 1, height_ratios=[1 - blq_panel_height, blq_panel_height],
                      hspace=0.1)

        # Upper panel: VPC
        ax_vpc = fig.add_subplot(gs[0])
        colors = get_colors(len(pi_levels))

        for i, level in enumerate(pi_levels):
            obs_values = get_observed_percentile(vpc_result, level)
            sim_lower, sim_upper = get_simulated_ci(vpc_result, level)

            if len(sim_lower) > 0:
                ax_vpc.fill_between(times, sim_lower, sim_upper,
                                    color=colors[i], alpha=0.2)
            if len(obs_values) > 0:
                ax_vpc.plot(times, obs_values, color=colors[i],
                            linewidth=theme["line_width"],
                            label=f"{int(level * 100)}th pctl")

        ax_vpc.set_ylabel(ylabel)
        if log_scale:
            ax_vpc.set_yscale("log")
        ax_vpc.legend(loc="upper right", fontsize="small")
        ax_vpc.set_xticklabels([])

        if title:
            ax_vpc.set_title(title)

        # Lower panel: BLQ
        ax_blq = fig.add_subplot(gs[1], sharex=ax_vpc)

        blq_obs = get_blq_observed(vpc_result)
        blq_sim_median, blq_sim_lower, blq_sim_upper = get_blq_simulated(vpc_result)

        if len(blq_obs) > 0:
            # BLQ CI ribbon
            if len(blq_sim_lower) > 0:
                ax_blq.fill_between(times, blq_sim_lower * 100, blq_sim_upper * 100,
                                    color=blq_color, alpha=0.2, label="Sim %BLQ CI")
            # Simulated median
            if len(blq_sim_median) > 0:
                ax_blq.plot(times, blq_sim_median * 100, color=blq_color,
                            linestyle="--", label="Sim %BLQ")
            # Observed
            ax_blq.plot(times, blq_obs * 100, color=obs_color,
                        linewidth=theme["line_width"], marker='o',
                        markersize=4, label="Obs %BLQ")

            ax_blq.set_ylim(0, 100)
            ax_blq.legend(loc="upper right", fontsize="small")

        ax_blq.set_xlabel(xlabel)
        ax_blq.set_ylabel(blq_ylabel)

        fig.tight_layout()
        return fig

    else:
        # Plotly implementation
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(
            rows=2, cols=1,
            row_heights=[1 - blq_panel_height, blq_panel_height],
            shared_xaxes=True,
            vertical_spacing=0.05
        )

        colors = get_colors(len(pi_levels))

        # Upper panel: VPC
        for i, level in enumerate(pi_levels):
            obs_values = get_observed_percentile(vpc_result, level)
            sim_lower, sim_upper = get_simulated_ci(vpc_result, level)

            if len(sim_lower) > 0:
                fig.add_trace(
                    go.Scatter(x=times, y=sim_upper, mode="lines",
                               line=dict(width=0), showlegend=False),
                    row=1, col=1
                )
                fig.add_trace(
                    go.Scatter(x=times, y=sim_lower, mode="lines",
                               line=dict(width=0), fill="tonexty",
                               fillcolor=f"rgba(0,100,200,0.2)",
                               showlegend=False),
                    row=1, col=1
                )

            if len(obs_values) > 0:
                fig.add_trace(
                    go.Scatter(x=times, y=obs_values, mode="lines",
                               name=f"{int(level * 100)}th pctl",
                               line=dict(color=colors[i])),
                    row=1, col=1
                )

        # Lower panel: BLQ
        blq_obs = get_blq_observed(vpc_result)
        blq_sim_median, blq_sim_lower, blq_sim_upper = get_blq_simulated(vpc_result)

        if len(blq_obs) > 0:
            if len(blq_sim_lower) > 0:
                fig.add_trace(
                    go.Scatter(x=times, y=blq_sim_upper * 100, mode="lines",
                               line=dict(width=0), showlegend=False),
                    row=2, col=1
                )
                fig.add_trace(
                    go.Scatter(x=times, y=blq_sim_lower * 100, mode="lines",
                               line=dict(width=0), fill="tonexty",
                               name="Sim %BLQ CI"),
                    row=2, col=1
                )

            fig.add_trace(
                go.Scatter(x=times, y=blq_obs * 100, mode="lines+markers",
                           name="Obs %BLQ", line=dict(color=obs_color)),
                row=2, col=1
            )

        fig.update_layout(
            height=figsize[1] * 100,
            width=figsize[0] * 100,
            title=title
        )
        fig.update_yaxes(title_text=ylabel, row=1, col=1)
        fig.update_yaxes(title_text=blq_ylabel, range=[0, 100], row=2, col=1)
        fig.update_xaxes(title_text=xlabel, row=2, col=1)

        if log_scale:
            fig.update_yaxes(type="log", row=1, col=1)

        return fig


# ============================================================================
# VPC with CI Ribbons (Detailed)
# ============================================================================

def plot_vpc_ci(
    vpc_result: Any,
    log_scale: bool = False,
    ci_levels: List[float] = [0.05, 0.95],
    show_bin_boundaries: bool = False,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Plot VPC with detailed confidence interval ribbons.

    Shows nested CI ribbons for each percentile, making it easy to
    assess whether observed percentiles fall within expected simulation
    uncertainty bands.

    Args:
        vpc_result: VPCResult from compute_vpc
        log_scale: Use log scale for y-axis
        ci_levels: CI levels to highlight [lower, upper]
        show_bin_boundaries: Show vertical lines at bin boundaries
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> fig = plot_vpc_ci(result, ci_levels=[0.025, 0.975])
    """
    from ..vpc import (
        get_bin_midpoints, get_observed_percentile,
        get_simulated_median, get_simulated_ci
    )

    plotter = _get_plotter()
    theme = get_theme_config()

    times = get_bin_midpoints(vpc_result)
    pi_levels = vpc_result.config.pi_levels

    # Color scheme for different percentiles
    # 5th percentile - blue shades
    # 50th percentile - green shades
    # 95th percentile - red shades
    level_colors = {
        0.05: {"ci": "#3498DB", "obs": "#2E86AB", "sim": "#5DADE2"},
        0.50: {"ci": "#27AE60", "obs": "#27AE60", "sim": "#58D68D"},
        0.95: {"ci": "#E74C3C", "obs": "#E74C3C", "sim": "#EC7063"},
    }

    # Default colors for other levels
    default_colors = get_colors(len(pi_levels))

    fig = plotter.create_figure(figsize=figsize, title=title or "VPC with Confidence Intervals")

    # Plot from outer to inner percentiles for proper layering
    for i, level in enumerate(sorted(pi_levels)):
        # Get colors for this level
        if level in level_colors:
            colors = level_colors[level]
            ci_color = colors["ci"]
            obs_color = colors["obs"]
        else:
            obs_color = default_colors[i]
            ci_color = default_colors[i]

        obs_values = get_observed_percentile(vpc_result, level)
        sim_median = get_simulated_median(vpc_result, level)
        sim_lower, sim_upper = get_simulated_ci(vpc_result, level)

        pct_label = f"{int(level * 100)}th"

        # Plot CI ribbon
        if len(sim_lower) > 0:
            plotter.fill_between(
                fig, list(times), list(sim_lower), list(sim_upper),
                color=ci_color, alpha=0.4,
                label=f"{pct_label} CI"
            )

        # Plot observed line (solid, thicker)
        if len(obs_values) > 0:
            plotter.line_plot(
                fig, list(times), list(obs_values),
                color=obs_color, linestyle="-",
                linewidth=theme["line_width"] * 1.5,
                label=f"Obs {pct_label}"
            )

    # Add bin boundaries if requested
    if show_bin_boundaries:
        backend = plotter.__class__.__name__
        if "Matplotlib" in backend:
            ax = fig["ax"]
            for bin_data in vpc_result.bins:
                ax.axvline(bin_data.time_min, color="gray", linestyle=":",
                           alpha=0.3, linewidth=0.5)
                ax.axvline(bin_data.time_max, color="gray", linestyle=":",
                           alpha=0.3, linewidth=0.5)

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel,
                       title=title or "VPC with Confidence Intervals")

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig, location="best")

    return plotter.finalize(fig)
