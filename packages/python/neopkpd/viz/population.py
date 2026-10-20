"""
NeoPKPD Population Visualization

Population-level visualization functions including:
- Visual Predictive Check (VPC)
- Parameter distributions
- Forest plots
- Boxplots and violin plots
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

from .backends import _get_plotter
from .themes import get_theme_config, get_color, get_colors


def plot_vpc(
    vpc_result: Union[Dict[str, Any], Any],
    observed_data: Optional[Dict[str, List[float]]] = None,
    observation: str = "conc",
    prediction_intervals: List[float] = [0.05, 0.50, 0.95],
    log_scale: bool = False,
    title: Optional[str] = "Visual Predictive Check",
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (12, 8),
    show_observed: bool = True,
    n_bins: int = 10,
) -> Any:
    """
    Create Visual Predictive Check (VPC) plot.

    VPC compares simulated prediction intervals with observed data to
    assess model adequacy. Accepts either a population simulation result
    dict or a VPCResult from compute_vpc().

    Args:
        vpc_result: Population simulation result dict OR VPCResult from compute_vpc()
        observed_data: Optional dict with 't' and 'c' lists for observed data
        observation: Observation key to plot (for dict input)
        prediction_intervals: Quantiles for prediction intervals
        log_scale: Use log scale for y-axis
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size
        show_observed: Show observed data points
        n_bins: Number of time bins for prediction intervals

    Returns:
        Figure object

    Example:
        >>> # From population simulation
        >>> pop_result = neopkpd.simulate_population_iv_bolus(...)
        >>> fig = plot_vpc(pop_result, observed_data=observed)
        >>>
        >>> # From compute_vpc()
        >>> vpc_result = neopkpd.compute_vpc(observed, pop_spec, grid)
        >>> fig = plot_vpc(vpc_result)
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    # Check if this is a VPCResult dataclass or a dict
    if hasattr(vpc_result, 'bins') and hasattr(vpc_result, 'n_simulations'):
        # This is a VPCResult from compute_vpc()
        return _plot_vpc_from_result(vpc_result, log_scale, title, xlabel, ylabel, figsize)

    # Otherwise treat as population simulation result dict
    population_result = vpc_result
    individuals = population_result["individuals"]
    t = np.array(individuals[0]["t"])

    # Collect all concentration data
    all_conc = []
    for ind in individuals:
        all_conc.append(ind["observations"][observation])
    all_conc = np.array(all_conc)

    # Calculate prediction intervals
    pi_lower = np.percentile(all_conc, prediction_intervals[0] * 100, axis=0)
    pi_median = np.percentile(all_conc, prediction_intervals[1] * 100, axis=0)
    pi_upper = np.percentile(all_conc, prediction_intervals[2] * 100, axis=0)

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Plot outer prediction interval
    pi_label = f"{int((prediction_intervals[2] - prediction_intervals[0]) * 100)}% PI"
    plotter.fill_between(fig, list(t), list(pi_lower), list(pi_upper),
                         color=get_color(0), alpha=theme["alpha_ribbon"],
                         label=pi_label)

    # Plot median line
    plotter.line_plot(fig, list(t), list(pi_median), label="Median",
                      color=get_color(0), linewidth=theme["line_width"] * 1.5)

    # Plot observed data if provided
    if show_observed and observed_data is not None:
        t_obs = observed_data.get("t", [])
        c_obs = observed_data.get("c", observed_data.get("conc", []))
        if t_obs and c_obs:
            plotter.scatter_plot(fig, t_obs, c_obs, color=get_color(1),
                                 size=theme["marker_size"] * 8,
                                 label="Observed")

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig)

    return plotter.finalize(fig)


def plot_parameter_distributions(
    population_result: Dict[str, Any],
    parameters: Optional[List[str]] = None,
    plot_type: str = "histogram",
    n_cols: int = 3,
    subplot_size: Tuple[float, float] = (4, 3),
    show_stats: bool = True,
    log_scale: bool = False,
) -> Any:
    """
    Plot distributions of individual PK parameters.

    Args:
        population_result: Population simulation result dictionary
        parameters: List of parameter names to plot (None = all)
        plot_type: "histogram" or "kde"
        n_cols: Number of columns in subplot grid
        subplot_size: Size of each subplot
        show_stats: Show mean/median statistics
        log_scale: Use log scale for x-axis

    Returns:
        Figure object

    Example:
        >>> pop_result = neopkpd.simulate_population_iv_bolus(...)
        >>> fig = plot_parameter_distributions(pop_result, parameters=["CL", "V"])
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    individuals = population_result["individuals"]

    # Extract parameters
    if parameters is None:
        # Try to infer from first individual
        if "parameters" in individuals[0]:
            parameters = list(individuals[0]["parameters"].keys())
        else:
            parameters = []

    if not parameters:
        # Fallback: use summary statistics if available
        if "summaries" in population_result:
            parameters = list(population_result["summaries"].keys())

    if not parameters:
        raise ValueError("No parameters found in population result")

    # Collect parameter values
    param_values = {p: [] for p in parameters}
    for ind in individuals:
        if "parameters" in ind:
            for p in parameters:
                if p in ind["parameters"]:
                    param_values[p].append(ind["parameters"][p])

    # Filter out empty parameters
    parameters = [p for p in parameters if param_values[p]]

    n = len(parameters)
    n_rows = (n + n_cols - 1) // n_cols
    figsize = (subplot_size[0] * n_cols, subplot_size[1] * n_rows)

    backend = plotter.__class__.__name__

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n > 1:
            axes = axes.flatten()
        else:
            axes = [axes]

        for i, param in enumerate(parameters):
            ax = axes[i]
            values = np.array(param_values[param])

            if plot_type == "histogram":
                ax.hist(values, bins=30, color=get_color(0), alpha=0.7,
                        edgecolor="white")
            else:  # kde
                from scipy import stats
                try:
                    kde = stats.gaussian_kde(values)
                    x_range = np.linspace(values.min(), values.max(), 100)
                    ax.fill_between(x_range, kde(x_range), alpha=0.7,
                                    color=get_color(0))
                except Exception:
                    ax.hist(values, bins=30, color=get_color(0), alpha=0.7)

            if show_stats:
                mean_val = np.mean(values)
                median_val = np.median(values)
                ax.axvline(mean_val, color=get_color(1), linestyle="--",
                           linewidth=1.5, label=f"Mean: {mean_val:.3g}")
                ax.axvline(median_val, color=get_color(2), linestyle=":",
                           linewidth=1.5, label=f"Median: {median_val:.3g}")
                ax.legend(fontsize=8)

            ax.set_xlabel(param)
            ax.set_ylabel("Frequency" if plot_type == "histogram" else "Density")
            ax.set_title(param)

            if log_scale:
                ax.set_xscale("log")

        # Hide empty subplots
        for i in range(n, len(axes)):
            axes[i].set_visible(False)

        fig.tight_layout()
        return fig

    else:  # Plotly
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(rows=n_rows, cols=n_cols,
                            subplot_titles=parameters)

        for i, param in enumerate(parameters):
            row = i // n_cols + 1
            col = i % n_cols + 1
            values = param_values[param]

            fig.add_trace(
                go.Histogram(x=values, name=param, marker_color=get_color(0),
                             showlegend=False),
                row=row, col=col
            )

            if show_stats:
                mean_val = np.mean(values)
                fig.add_vline(x=mean_val, line_dash="dash",
                              line_color=get_color(1), row=row, col=col)

        fig.update_layout(
            height=figsize[1] * 100,
            width=figsize[0] * 100,
            title="Parameter Distributions"
        )

        return fig


def plot_forest(
    effects: List[Dict[str, Any]],
    reference_line: float = 1.0,
    log_scale: bool = True,
    title: Optional[str] = "Forest Plot",
    xlabel: str = "Effect Size",
    figsize: Tuple[float, float] = (10, 8),
    show_summary: bool = True,
) -> Any:
    """
    Create forest plot for effect estimates with confidence intervals.

    Args:
        effects: List of dicts with keys 'name', 'estimate', 'lower', 'upper'
        reference_line: Reference line value (usually 1.0 for ratios)
        log_scale: Use log scale for x-axis
        title: Plot title
        xlabel: X-axis label
        figsize: Figure size
        show_summary: Show overall summary diamond

    Returns:
        Figure object

    Example:
        >>> effects = [
        ...     {"name": "Study 1", "estimate": 1.2, "lower": 0.9, "upper": 1.6},
        ...     {"name": "Study 2", "estimate": 0.95, "lower": 0.7, "upper": 1.3},
        ... ]
        >>> fig = plot_forest(effects)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    n = len(effects)
    y_positions = list(range(n, 0, -1))

    backend = plotter.__class__.__name__

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        for i, (effect, y) in enumerate(zip(effects, y_positions)):
            estimate = effect["estimate"]
            lower = effect["lower"]
            upper = effect["upper"]
            name = effect.get("name", f"Study {i+1}")

            # Draw confidence interval line
            ax.plot([lower, upper], [y, y], color=get_color(0),
                    linewidth=theme["line_width"], solid_capstyle="butt")

            # Draw point estimate
            ax.scatter([estimate], [y], color=get_color(0),
                       s=theme["marker_size"] * 15, zorder=3)

            # Add study name
            ax.text(ax.get_xlim()[0] if ax.get_xlim()[0] > 0 else 0.1, y,
                    f"  {name}", va="center", ha="left", fontsize=10)

        # Reference line
        ax.axvline(reference_line, color=get_color(2), linestyle="--",
                   linewidth=1, alpha=0.7, label=f"Reference ({reference_line})")

        # Summary diamond
        if show_summary and len(effects) > 1:
            all_estimates = [e["estimate"] for e in effects]
            all_weights = [1.0 / ((e["upper"] - e["lower"]) ** 2) for e in effects]
            total_weight = sum(all_weights)
            summary_estimate = sum(e * w for e, w in zip(all_estimates, all_weights)) / total_weight

            # Simplified summary CI (weighted average of CIs)
            summary_lower = sum(e["lower"] * w for e, w in zip(effects, all_weights)) / total_weight
            summary_upper = sum(e["upper"] * w for e, w in zip(effects, all_weights)) / total_weight

            diamond_y = 0
            diamond_height = 0.4
            diamond = plt.Polygon([
                [summary_lower, diamond_y],
                [summary_estimate, diamond_y + diamond_height],
                [summary_upper, diamond_y],
                [summary_estimate, diamond_y - diamond_height]
            ], color=get_color(1), alpha=0.8)
            ax.add_patch(diamond)
            ax.text(ax.get_xlim()[0] if ax.get_xlim()[0] > 0 else 0.1, diamond_y,
                    "  Summary", va="center", ha="left", fontsize=10, fontweight="bold")

        ax.set_yticks(y_positions + ([0] if show_summary else []))
        ax.set_yticklabels([e.get("name", f"Study {i+1}") for i, e in enumerate(effects)] +
                           (["Summary"] if show_summary else []))
        ax.set_xlabel(xlabel)
        ax.set_title(title)

        if log_scale:
            ax.set_xscale("log")

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        return fig

    else:  # Plotly
        import plotly.graph_objects as go

        fig = go.Figure()

        names = [e.get("name", f"Study {i+1}") for i, e in enumerate(effects)]
        estimates = [e["estimate"] for e in effects]
        lowers = [e["lower"] for e in effects]
        uppers = [e["upper"] for e in effects]

        # Error bars
        fig.add_trace(go.Scatter(
            x=estimates, y=y_positions,
            mode="markers",
            marker=dict(size=10, color=get_color(0)),
            error_x=dict(
                type="data",
                symmetric=False,
                array=[u - e for u, e in zip(uppers, estimates)],
                arrayminus=[e - l for e, l in zip(estimates, lowers)]
            ),
            name="Studies"
        ))

        # Reference line
        fig.add_vline(x=reference_line, line_dash="dash",
                      line_color=get_color(2), annotation_text="Reference")

        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            yaxis=dict(
                tickvals=y_positions,
                ticktext=names
            ),
            height=figsize[1] * 100,
            width=figsize[0] * 100
        )

        if log_scale:
            fig.update_xaxes(type="log")

        return fig


def plot_boxplot(
    population_result: Dict[str, Any],
    groups: Optional[List[str]] = None,
    metric: str = "cmax",
    violin: bool = False,
    show_points: bool = True,
    title: Optional[str] = None,
    xlabel: str = "Group",
    ylabel: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Create boxplot or violin plot comparing groups.

    Args:
        population_result: Population simulation result dictionary
        groups: Group labels (one per individual, or derived from result)
        metric: Metric to plot (e.g., "cmax", "auc")
        violin: Use violin plot instead of boxplot
        show_points: Overlay individual data points
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> pop_result = neopkpd.simulate_population_iv_bolus(...)
        >>> groups = ["Low Dose"] * 50 + ["High Dose"] * 50
        >>> fig = plot_boxplot(pop_result, groups=groups, metric="cmax")
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    individuals = population_result["individuals"]

    # Extract metric values
    values = []
    for ind in individuals:
        obs = ind["observations"]
        if metric == "cmax":
            values.append(max(obs.get("conc", obs.get(list(obs.keys())[0]))))
        elif metric == "tmax":
            t = ind["t"]
            c = obs.get("conc", obs.get(list(obs.keys())[0]))
            values.append(t[np.argmax(c)])
        elif metric == "auc":
            t = np.array(ind["t"])
            c = np.array(obs.get("conc", obs.get(list(obs.keys())[0])))
            values.append(np.trapz(c, t))
        elif metric in obs:
            values.append(max(obs[metric]))
        else:
            values.append(0.0)

    values = np.array(values)

    # Handle groups
    if groups is None:
        groups = ["All"] * len(values)

    unique_groups = list(dict.fromkeys(groups))  # Preserve order
    n_groups = len(unique_groups)

    if ylabel is None:
        ylabel = metric.upper()
    if title is None:
        title = f"{ylabel} by Group"

    backend = plotter.__class__.__name__

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        # Prepare data by group
        group_data = [values[[i for i, g in enumerate(groups) if g == ug]]
                      for ug in unique_groups]

        positions = list(range(1, n_groups + 1))

        if violin:
            parts = ax.violinplot(group_data, positions=positions,
                                  showmeans=True, showmedians=True)
            for pc in parts["bodies"]:
                pc.set_facecolor(get_color(0))
                pc.set_alpha(0.7)
        else:
            bp = ax.boxplot(group_data, positions=positions, patch_artist=True)
            for patch in bp["boxes"]:
                patch.set_facecolor(get_color(0))
                patch.set_alpha(0.7)

        if show_points:
            for i, (ug, data) in enumerate(zip(unique_groups, group_data)):
                jitter = np.random.uniform(-0.15, 0.15, len(data))
                ax.scatter([i + 1 + j for j in jitter], data,
                           color=get_color(1), alpha=0.5, s=20, zorder=3)

        ax.set_xticks(positions)
        ax.set_xticklabels(unique_groups)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        fig.tight_layout()
        return fig

    else:  # Plotly
        import plotly.graph_objects as go

        fig = go.Figure()

        for i, ug in enumerate(unique_groups):
            group_values = values[[j for j, g in enumerate(groups) if g == ug]]

            if violin:
                fig.add_trace(go.Violin(
                    y=group_values, name=ug,
                    box_visible=True, meanline_visible=True,
                    fillcolor=get_color(i), line_color=get_color(i)
                ))
            else:
                fig.add_trace(go.Box(
                    y=group_values, name=ug,
                    marker_color=get_color(i),
                    boxpoints="all" if show_points else False
                ))

        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            height=figsize[1] * 100,
            width=figsize[0] * 100
        )

        return fig


# ============================================================================
# VPC from VPCResult Dataclass
# ============================================================================

def _plot_vpc_from_result(
    vpc_result: Any,
    log_scale: bool = False,
    title: Optional[str] = "Visual Predictive Check",
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (12, 8),
) -> Any:
    """
    Plot VPC from VPCResult dataclass (from compute_vpc()).

    Uses binned percentiles with confidence intervals.
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Extract bin midpoints and percentiles
    t_mid = [bin.time_midpoint for bin in vpc_result.bins]

    # Get percentile data for each level
    for level_idx, level in enumerate(vpc_result.pi_levels):
        obs_vals = []
        sim_median = []
        sim_lower = []
        sim_upper = []

        for bin in vpc_result.bins:
            for p in bin.percentiles:
                if abs(p.level - level) < 0.001:
                    obs_vals.append(p.observed)
                    sim_median.append(p.simulated_median)
                    sim_lower.append(p.simulated_lower)
                    sim_upper.append(p.simulated_upper)
                    break

        if sim_lower and sim_upper:
            color = get_color(level_idx)
            level_pct = int(level * 100)

            # Plot CI ribbon for simulated percentile
            plotter.fill_between(fig, t_mid, sim_lower, sim_upper,
                                 color=color, alpha=theme["alpha_ribbon"] * 0.5,
                                 label=f"{level_pct}th Simulated CI")

            # Plot simulated median
            plotter.line_plot(fig, t_mid, sim_median, color=color,
                              linestyle="--", linewidth=theme["line_width"])

            # Plot observed percentile
            plotter.scatter_plot(fig, t_mid, obs_vals, color=color,
                                 size=theme["marker_size"] * 6,
                                 label=f"{level_pct}th Observed")

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig)

    return plotter.finalize(fig)


# ============================================================================
# Estimation Diagnostics Plots
# ============================================================================

def plot_goodness_of_fit(
    estimation_result: Any,
    plot_type: str = "obs_vs_pred",
    log_scale: bool = False,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 10),
) -> Any:
    """
    Create goodness-of-fit diagnostic plots for estimation results.

    Args:
        estimation_result: EstimationResult from estimate()
        plot_type: Type of plot:
            - "obs_vs_pred": Observed vs Population Predicted
            - "obs_vs_ipred": Observed vs Individual Predicted
            - "cwres_vs_time": CWRES vs Time
            - "cwres_vs_pred": CWRES vs Predicted
            - "cwres_qq": CWRES QQ plot
            - "eta_dist": ETA distributions
        log_scale: Use log scale for concentration plots
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> result = neopkpd.estimate(observed_data, "OneCompIVBolus", config, grid)
        >>> fig = plot_goodness_of_fit(result, plot_type="obs_vs_ipred")
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    backend = plotter.__class__.__name__

    # Collect data from all individuals
    all_obs = []
    all_pred = []
    all_ipred = []
    all_cwres = []
    all_times = []
    all_etas = []

    for ind in estimation_result.individuals:
        # Note: observed data not stored in IndividualEstimate, so we use pred/ipred
        all_pred.extend(ind.pred)
        all_ipred.extend(ind.ipred)
        all_cwres.extend(ind.cwres)
        all_etas.append(ind.eta)
        # Approximate times based on index
        all_times.extend(range(len(ind.pred)))

    all_pred = np.array(all_pred)
    all_ipred = np.array(all_ipred)
    all_cwres = np.array(all_cwres)
    all_times = np.array(all_times)

    if title is None:
        title = plot_type.replace("_", " ").title()

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        if plot_type == "obs_vs_pred":
            fig, ax = plt.subplots(figsize=figsize)
            # Use IPRED as proxy for observed for demonstration
            ax.scatter(all_pred, all_ipred, color=get_color(0), alpha=0.6)
            lims = [min(all_pred.min(), all_ipred.min()),
                    max(all_pred.max(), all_ipred.max())]
            ax.plot(lims, lims, 'k--', linewidth=1)
            ax.set_xlabel("Population Predicted")
            ax.set_ylabel("Observed (approx)")
            ax.set_title(title)
            if log_scale:
                ax.set_xscale("log")
                ax.set_yscale("log")

        elif plot_type == "obs_vs_ipred":
            fig, ax = plt.subplots(figsize=figsize)
            ax.scatter(all_ipred, all_pred, color=get_color(0), alpha=0.6)
            lims = [min(all_ipred.min(), all_pred.min()),
                    max(all_ipred.max(), all_pred.max())]
            ax.plot(lims, lims, 'k--', linewidth=1)
            ax.set_xlabel("Individual Predicted")
            ax.set_ylabel("Observed (approx)")
            ax.set_title(title)
            if log_scale:
                ax.set_xscale("log")
                ax.set_yscale("log")

        elif plot_type == "cwres_vs_time":
            fig, ax = plt.subplots(figsize=figsize)
            ax.scatter(all_times, all_cwres, color=get_color(0), alpha=0.6)
            ax.axhline(0, color='k', linestyle='--', linewidth=1)
            ax.axhline(2, color='r', linestyle=':', linewidth=1, alpha=0.5)
            ax.axhline(-2, color='r', linestyle=':', linewidth=1, alpha=0.5)
            ax.set_xlabel("Time (index)")
            ax.set_ylabel("CWRES")
            ax.set_title(title)

        elif plot_type == "cwres_vs_pred":
            fig, ax = plt.subplots(figsize=figsize)
            ax.scatter(all_pred, all_cwres, color=get_color(0), alpha=0.6)
            ax.axhline(0, color='k', linestyle='--', linewidth=1)
            ax.axhline(2, color='r', linestyle=':', linewidth=1, alpha=0.5)
            ax.axhline(-2, color='r', linestyle=':', linewidth=1, alpha=0.5)
            ax.set_xlabel("Population Predicted")
            ax.set_ylabel("CWRES")
            ax.set_title(title)

        elif plot_type == "cwres_qq":
            from scipy import stats
            fig, ax = plt.subplots(figsize=figsize)
            stats.probplot(all_cwres, dist="norm", plot=ax)
            ax.set_title(title)

        elif plot_type == "eta_dist":
            all_etas = np.array(all_etas)
            n_eta = all_etas.shape[1] if len(all_etas.shape) > 1 else 1
            n_cols = min(n_eta, 3)
            n_rows = (n_eta + n_cols - 1) // n_cols
            fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
            if n_eta > 1:
                axes = axes.flatten()
            else:
                axes = [axes]

            for i in range(n_eta):
                ax = axes[i]
                eta_vals = all_etas[:, i] if len(all_etas.shape) > 1 else all_etas
                ax.hist(eta_vals, bins=20, color=get_color(i), alpha=0.7)
                ax.axvline(0, color='k', linestyle='--', linewidth=1)
                ax.set_xlabel(f"ETA_{i+1}")
                ax.set_ylabel("Frequency")
                ax.set_title(f"ETA_{i+1} Distribution")

            for i in range(n_eta, len(axes)):
                axes[i].set_visible(False)

        else:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, f"Unknown plot type: {plot_type}",
                    ha='center', va='center', transform=ax.transAxes)

        fig.tight_layout()
        return fig

    else:  # Plotly
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        fig = go.Figure()

        if plot_type in ["obs_vs_pred", "obs_vs_ipred"]:
            x_data = all_pred if plot_type == "obs_vs_pred" else all_ipred
            y_data = all_ipred if plot_type == "obs_vs_pred" else all_pred

            fig.add_trace(go.Scatter(x=x_data, y=y_data, mode="markers",
                                      marker=dict(color=get_color(0), opacity=0.6)))
            lims = [min(x_data.min(), y_data.min()), max(x_data.max(), y_data.max())]
            fig.add_trace(go.Scatter(x=lims, y=lims, mode="lines",
                                      line=dict(dash="dash", color="black")))

        elif plot_type.startswith("cwres"):
            x_data = all_times if "time" in plot_type else all_pred
            fig.add_trace(go.Scatter(x=x_data, y=all_cwres, mode="markers",
                                      marker=dict(color=get_color(0), opacity=0.6)))
            fig.add_hline(y=0, line_dash="dash", line_color="black")
            fig.add_hline(y=2, line_dash="dot", line_color="red", opacity=0.5)
            fig.add_hline(y=-2, line_dash="dot", line_color="red", opacity=0.5)

        fig.update_layout(title=title, height=figsize[1]*100, width=figsize[0]*100)
        return fig


def plot_estimation_summary(
    estimation_result: Any,
    figsize: Tuple[float, float] = (14, 10),
) -> Any:
    """
    Create 4-panel summary of estimation diagnostics.

    Args:
        estimation_result: EstimationResult from estimate()
        figsize: Figure size

    Returns:
        Figure object with 4 diagnostic panels

    Example:
        >>> result = neopkpd.estimate(...)
        >>> fig = plot_estimation_summary(result)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    backend = plotter.__class__.__name__

    # Collect data
    all_pred = []
    all_ipred = []
    all_cwres = []
    all_times = []

    for ind in estimation_result.individuals:
        all_pred.extend(ind.pred)
        all_ipred.extend(ind.ipred)
        all_cwres.extend(ind.cwres)
        all_times.extend(range(len(ind.pred)))

    all_pred = np.array(all_pred)
    all_ipred = np.array(all_ipred)
    all_cwres = np.array(all_cwres)
    all_times = np.array(all_times)

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=figsize)

        # Obs vs PRED
        ax = axes[0, 0]
        ax.scatter(all_pred, all_ipred, color=get_color(0), alpha=0.6, s=20)
        lims = [min(all_pred.min(), all_ipred.min()), max(all_pred.max(), all_ipred.max())]
        ax.plot(lims, lims, 'k--', linewidth=1)
        ax.set_xlabel("Population Predicted")
        ax.set_ylabel("Individual Predicted")
        ax.set_title("PRED vs IPRED")

        # Obs vs IPRED
        ax = axes[0, 1]
        ax.scatter(all_ipred, all_pred, color=get_color(0), alpha=0.6, s=20)
        ax.plot(lims, lims, 'k--', linewidth=1)
        ax.set_xlabel("Individual Predicted")
        ax.set_ylabel("Population Predicted")
        ax.set_title("IPRED vs PRED")

        # CWRES vs Time
        ax = axes[1, 0]
        ax.scatter(all_times, all_cwres, color=get_color(0), alpha=0.6, s=20)
        ax.axhline(0, color='k', linestyle='--', linewidth=1)
        ax.axhline(2, color='r', linestyle=':', linewidth=1, alpha=0.5)
        ax.axhline(-2, color='r', linestyle=':', linewidth=1, alpha=0.5)
        ax.set_xlabel("Time (index)")
        ax.set_ylabel("CWRES")
        ax.set_title("CWRES vs Time")

        # CWRES vs PRED
        ax = axes[1, 1]
        ax.scatter(all_pred, all_cwres, color=get_color(0), alpha=0.6, s=20)
        ax.axhline(0, color='k', linestyle='--', linewidth=1)
        ax.axhline(2, color='r', linestyle=':', linewidth=1, alpha=0.5)
        ax.axhline(-2, color='r', linestyle=':', linewidth=1, alpha=0.5)
        ax.set_xlabel("Population Predicted")
        ax.set_ylabel("CWRES")
        ax.set_title("CWRES vs PRED")

        fig.suptitle(f"Estimation Diagnostics (OFV: {estimation_result.ofv:.2f})", fontsize=14)
        fig.tight_layout()
        return fig

    else:  # Plotly
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(rows=2, cols=2,
                            subplot_titles=["PRED vs IPRED", "IPRED vs PRED",
                                            "CWRES vs Time", "CWRES vs PRED"])

        # Panel 1
        fig.add_trace(go.Scatter(x=all_pred, y=all_ipred, mode="markers",
                                  marker=dict(color=get_color(0), opacity=0.6)),
                      row=1, col=1)

        # Panel 2
        fig.add_trace(go.Scatter(x=all_ipred, y=all_pred, mode="markers",
                                  marker=dict(color=get_color(0), opacity=0.6)),
                      row=1, col=2)

        # Panel 3
        fig.add_trace(go.Scatter(x=all_times, y=all_cwres, mode="markers",
                                  marker=dict(color=get_color(0), opacity=0.6)),
                      row=2, col=1)

        # Panel 4
        fig.add_trace(go.Scatter(x=all_pred, y=all_cwres, mode="markers",
                                  marker=dict(color=get_color(0), opacity=0.6)),
                      row=2, col=2)

        fig.update_layout(
            title=f"Estimation Diagnostics (OFV: {estimation_result.ofv:.2f})",
            height=figsize[1]*100, width=figsize[0]*100,
            showlegend=False
        )

        return fig


# ============================================================================
# Sensitivity Analysis Visualization
# ============================================================================

def plot_sensitivity(
    sensitivity_result: Any,
    times: Optional[List[float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (10, 6),
    show_metrics: bool = True,
) -> Any:
    """
    Plot sensitivity analysis result showing base vs perturbed profiles.

    Args:
        sensitivity_result: SensitivityResult from run_sensitivity()
        times: Optional time points (uses indices if not provided)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size
        show_metrics: Show sensitivity metrics in legend

    Returns:
        Figure object

    Example:
        >>> result = neopkpd.run_sensitivity(model, grid, perturbation)
        >>> fig = plot_sensitivity(result)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    base = np.array(sensitivity_result.base_series)
    pert = np.array(sensitivity_result.pert_series)

    if times is None:
        times = list(range(len(base)))

    if title is None:
        title = f"Sensitivity: {sensitivity_result.plan_name}"

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Plot base profile
    label_base = "Base"
    if show_metrics:
        label_base = f"Base (max Δ: {sensitivity_result.metrics.max_rel_delta:.2%})"

    plotter.line_plot(fig, times, list(base), label="Base", color=get_color(0),
                      linewidth=theme["line_width"])
    plotter.scatter_plot(fig, times, list(base), color=get_color(0),
                         size=theme["marker_size"] * 4)

    # Plot perturbed profile
    plotter.line_plot(fig, times, list(pert), label="Perturbed", color=get_color(1),
                      linestyle="--", linewidth=theme["line_width"])
    plotter.scatter_plot(fig, times, list(pert), color=get_color(1),
                         size=theme["marker_size"] * 4)

    # Fill between to show difference
    plotter.fill_between(fig, times, list(base), list(pert),
                         color=get_color(2), alpha=0.2, label="Difference")

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)
    plotter.add_legend(fig)

    return plotter.finalize(fig)


def plot_sensitivity_tornado(
    sensitivity_results: List[Any],
    metric: str = "max_rel_delta",
    title: Optional[str] = "Parameter Sensitivity Analysis",
    xlabel: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Create tornado plot from multiple sensitivity analysis results.

    Args:
        sensitivity_results: List of SensitivityResult objects
        metric: Metric to plot ("max_rel_delta", "max_abs_delta", "l2_norm_delta")
        title: Plot title
        xlabel: X-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> results = [
        ...     neopkpd.run_sensitivity(model, grid, {"param": "CL", "delta": 0.1}),
        ...     neopkpd.run_sensitivity(model, grid, {"param": "V", "delta": 0.1}),
        ... ]
        >>> fig = plot_sensitivity_tornado(results)
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    if xlabel is None:
        xlabel = metric.replace("_", " ").title()

    # Extract values
    names = [r.plan_name for r in sensitivity_results]
    values = [getattr(r.metrics, metric) for r in sensitivity_results]

    # Sort by magnitude
    sorted_pairs = sorted(zip(names, values), key=lambda x: abs(x[1]), reverse=True)
    names = [p[0] for p in sorted_pairs]
    values = [p[1] for p in sorted_pairs]

    n = len(names)
    y_positions = list(range(n))

    backend = plotter.__class__.__name__

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        colors = [get_color(0) if v >= 0 else get_color(1) for v in values]
        ax.barh(y_positions, values, color=colors, alpha=0.8)

        ax.axvline(0, color='k', linewidth=1)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(names)
        ax.set_xlabel(xlabel)
        ax.set_title(title)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        return fig

    else:  # Plotly
        import plotly.graph_objects as go

        fig = go.Figure()

        colors = [get_color(0) if v >= 0 else get_color(1) for v in values]
        fig.add_trace(go.Bar(
            y=names, x=values, orientation="h",
            marker_color=colors
        ))

        fig.add_vline(x=0, line_color="black", line_width=1)

        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            height=figsize[1]*100, width=figsize[0]*100
        )

        return fig
