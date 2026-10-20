"""
NeoPKPD Sensitivity Analysis Visualization

Visualization functions for parameter sensitivity analysis including:
- Tornado diagrams
- Spider plots
- Sensitivity heatmaps
- Waterfall charts
- One-at-a-time (OFAT) sensitivity curves
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np

from .backends import _get_plotter, get_backend
from .themes import get_theme_config, get_color, get_colors, NEOPKPD_COLORS


# ============================================================================
# Tornado Diagram
# ============================================================================

def plot_tornado(
    sensitivities: Dict[str, Tuple[float, float]],
    metric_name: str = "AUC",
    baseline_value: Optional[float] = None,
    sort_by_magnitude: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Create tornado diagram showing parameter sensitivities.

    Tornado diagrams show the impact of parameter variations on an
    output metric, with bars extending from baseline in both directions.

    Args:
        sensitivities: Dict mapping parameter names to (low, high) output values
            when parameter is at its low/high bounds
        metric_name: Name of the output metric
        baseline_value: Baseline metric value (shown as vertical line)
        sort_by_magnitude: Sort parameters by sensitivity magnitude
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> sensitivities = {
        ...     "CL": (85.0, 115.0),   # AUC at -20% CL, +20% CL
        ...     "V": (98.0, 102.0),
        ...     "Ka": (90.0, 110.0),
        ... }
        >>> fig = plot_tornado(sensitivities, metric_name="AUC", baseline_value=100.0)
    """
    backend = get_backend()
    theme = get_theme_config()

    # Prepare data
    params = list(sensitivities.keys())
    low_vals = np.array([sensitivities[p][0] for p in params])
    high_vals = np.array([sensitivities[p][1] for p in params])

    if baseline_value is None:
        baseline_value = (np.mean(low_vals) + np.mean(high_vals)) / 2

    # Calculate ranges
    low_delta = low_vals - baseline_value
    high_delta = high_vals - baseline_value

    # Sort by magnitude if requested
    if sort_by_magnitude:
        magnitudes = np.abs(high_vals - low_vals)
        order = np.argsort(magnitudes)[::-1]
        params = [params[i] for i in order]
        low_delta = low_delta[order]
        high_delta = high_delta[order]

    n_params = len(params)
    y_pos = np.arange(n_params)

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        # Low side (typically negative if parameter decrease reduces output)
        bars_low = ax.barh(y_pos, low_delta, height=0.6,
                           color=NEOPKPD_COLORS["primary"], alpha=0.8,
                           label='-20% Parameter')

        # High side
        bars_high = ax.barh(y_pos, high_delta, height=0.6,
                            color=NEOPKPD_COLORS["secondary"], alpha=0.8,
                            label='+20% Parameter')

        # Baseline
        ax.axvline(0, color='black', linewidth=1)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(params)
        ax.set_xlabel(f"Change in {metric_name} from Baseline ({baseline_value:.2f})")
        ax.set_title(title or f"Parameter Sensitivity - {metric_name}")
        ax.legend(loc='lower right')
        ax.grid(axis='x', alpha=0.3)

        # Invert y-axis so most sensitive at top
        ax.invert_yaxis()

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        fig = go.Figure()

        fig.add_trace(go.Bar(
            y=params,
            x=low_delta,
            orientation='h',
            name='-20% Parameter',
            marker_color=NEOPKPD_COLORS["primary"]
        ))

        fig.add_trace(go.Bar(
            y=params,
            x=high_delta,
            orientation='h',
            name='+20% Parameter',
            marker_color=NEOPKPD_COLORS["secondary"]
        ))

        fig.add_vline(x=0, line_color="black", line_width=1)

        fig.update_layout(
            title=title or f"Parameter Sensitivity - {metric_name}",
            xaxis_title=f"Change in {metric_name} from Baseline ({baseline_value:.2f})",
            barmode='overlay',
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Spider Plot
# ============================================================================

def plot_spider(
    sensitivity_curves: Dict[str, Tuple[List[float], List[float]]],
    metric_name: str = "AUC",
    normalize: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Create spider (radar) plot showing parameter sensitivity curves.

    Shows how the output metric changes as each parameter is varied
    over its range, with all curves overlaid.

    Args:
        sensitivity_curves: Dict mapping parameter names to
            (parameter_values, metric_values) tuples
        metric_name: Name of the output metric
        normalize: Normalize to percent change from baseline
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> # Each entry: (param_values_as_fraction_of_nominal, metric_values)
        >>> curves = {
        ...     "CL": ([0.8, 0.9, 1.0, 1.1, 1.2], [125, 111, 100, 91, 83]),
        ...     "V":  ([0.8, 0.9, 1.0, 1.1, 1.2], [100, 100, 100, 100, 100]),
        ... }
        >>> fig = plot_spider(curves, metric_name="AUC")
    """
    backend = get_backend()
    theme = get_theme_config()

    params = list(sensitivity_curves.keys())
    colors = get_colors(len(params))

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        for i, (param, (x_vals, y_vals)) in enumerate(sensitivity_curves.items()):
            x = np.array(x_vals)
            y = np.array(y_vals)

            if normalize:
                # Find baseline (where x=1.0 or closest)
                baseline_idx = np.argmin(np.abs(x - 1.0))
                baseline = y[baseline_idx]
                y_plot = (y / baseline - 1) * 100  # Percent change
                ylabel = f"% Change in {metric_name}"
            else:
                y_plot = y
                ylabel = metric_name

            ax.plot(x * 100, y_plot, color=colors[i], linewidth=theme["line_width"],
                    marker='o', markersize=4, label=param)

        ax.axhline(0, color='black', linestyle='--', alpha=0.5)
        ax.axvline(100, color='black', linestyle='--', alpha=0.5)

        ax.set_xlabel("Parameter Value (% of Nominal)")
        ax.set_ylabel(ylabel)
        ax.set_title(title or f"Parameter Sensitivity - {metric_name}")
        ax.legend(loc='best')
        ax.grid(alpha=0.3)

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        fig = go.Figure()

        for i, (param, (x_vals, y_vals)) in enumerate(sensitivity_curves.items()):
            x = np.array(x_vals) * 100
            y = np.array(y_vals)

            if normalize:
                baseline_idx = np.argmin(np.abs(np.array(x_vals) - 1.0))
                baseline = y[baseline_idx]
                y_plot = (y / baseline - 1) * 100
            else:
                y_plot = y

            fig.add_trace(go.Scatter(
                x=x, y=y_plot,
                mode='lines+markers',
                name=param,
                line=dict(color=colors[i])
            ))

        fig.add_hline(y=0, line_dash="dash", line_color="black", opacity=0.5)
        fig.add_vline(x=100, line_dash="dash", line_color="black", opacity=0.5)

        ylabel = f"% Change in {metric_name}" if normalize else metric_name

        fig.update_layout(
            title=title or f"Parameter Sensitivity - {metric_name}",
            xaxis_title="Parameter Value (% of Nominal)",
            yaxis_title=ylabel,
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Sensitivity Heatmap
# ============================================================================

def plot_sensitivity_heatmap(
    sensitivity_matrix: np.ndarray,
    param_names: List[str],
    metric_names: List[str],
    normalize: bool = True,
    annotate: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 8),
    cmap: str = "RdYlBu_r",
) -> Any:
    """
    Create heatmap showing parameter-metric sensitivity matrix.

    Shows the sensitivity of multiple output metrics to multiple
    parameters in a matrix format.

    Args:
        sensitivity_matrix: 2D array of sensitivities [n_params x n_metrics]
        param_names: Names of parameters (rows)
        metric_names: Names of output metrics (columns)
        normalize: Normalize sensitivities to [-1, 1] range
        annotate: Show values in cells
        title: Plot title
        figsize: Figure size
        cmap: Colormap name

    Returns:
        Figure object

    Example:
        >>> # Sensitivity matrix: rows=params, cols=metrics
        >>> matrix = np.array([
        ...     [0.95, 0.85, 0.10],  # CL sensitivities for [AUC, Cmax, Tmax]
        ...     [0.02, 0.30, 0.01],  # V sensitivities
        ...     [0.15, 0.50, 0.80],  # Ka sensitivities
        ... ])
        >>> fig = plot_sensitivity_heatmap(matrix, ["CL", "V", "Ka"], ["AUC", "Cmax", "Tmax"])
    """
    backend = get_backend()

    if normalize:
        # Normalize each column independently
        matrix = sensitivity_matrix / np.abs(sensitivity_matrix).max(axis=0, keepdims=True)
        vmin, vmax = -1, 1
    else:
        matrix = sensitivity_matrix
        vmin, vmax = None, None

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
        plt.colorbar(im, ax=ax, label='Sensitivity')

        if annotate:
            for i in range(len(param_names)):
                for j in range(len(metric_names)):
                    val = matrix[i, j]
                    text = f"{val:.2f}"
                    color = "white" if abs(val) > 0.5 else "black"
                    ax.text(j, i, text, ha="center", va="center",
                            color=color, fontsize=9)

        ax.set_xticks(np.arange(len(metric_names)))
        ax.set_yticks(np.arange(len(param_names)))
        ax.set_xticklabels(metric_names)
        ax.set_yticklabels(param_names)
        ax.set_xlabel("Output Metric")
        ax.set_ylabel("Parameter")
        ax.set_title(title or "Parameter-Metric Sensitivity Matrix")

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=metric_names,
            y=param_names,
            colorscale=cmap,
            zmin=vmin, zmax=vmax,
            text=[[f"{v:.2f}" for v in row] for row in matrix] if annotate else None,
            texttemplate="%{text}" if annotate else None
        ))

        fig.update_layout(
            title=title or "Parameter-Metric Sensitivity Matrix",
            xaxis_title="Output Metric",
            yaxis_title="Parameter",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Waterfall Chart
# ============================================================================

def plot_waterfall(
    sensitivities: Dict[str, float],
    metric_name: str = "AUC",
    baseline_value: float = 100.0,
    sort_by_magnitude: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Create waterfall chart showing cumulative parameter contributions.

    Shows how each parameter contributes to the total change from
    baseline, with running totals.

    Args:
        sensitivities: Dict mapping parameter names to their contribution
            (change in metric from baseline)
        metric_name: Name of the output metric
        baseline_value: Starting baseline value
        sort_by_magnitude: Sort by absolute magnitude
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> contributions = {"CL": -15.0, "V": 2.0, "Ka": -5.0, "F": 3.0}
        >>> fig = plot_waterfall(contributions, baseline_value=100.0)
    """
    backend = get_backend()
    theme = get_theme_config()

    params = list(sensitivities.keys())
    values = np.array([sensitivities[p] for p in params])

    if sort_by_magnitude:
        order = np.argsort(np.abs(values))[::-1]
        params = [params[i] for i in order]
        values = values[order]

    n_params = len(params)

    # Calculate cumulative
    cumulative = np.zeros(n_params + 2)
    cumulative[0] = baseline_value
    cumulative[1:-1] = baseline_value + np.cumsum(values)
    cumulative[-1] = cumulative[-2]

    colors = [NEOPKPD_COLORS["success"] if v >= 0 else NEOPKPD_COLORS["error"]
              for v in values]

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        # All labels including baseline and total
        labels = ["Baseline"] + params + ["Total"]
        x_pos = np.arange(len(labels))

        # Plot bars
        for i in range(1, len(labels) - 1):
            bottom = cumulative[i - 1]
            height = cumulative[i] - cumulative[i - 1]
            color = colors[i - 1]
            ax.bar(x_pos[i], height, bottom=bottom, color=color, width=0.6)

            # Connector line
            ax.plot([x_pos[i - 1] + 0.3, x_pos[i] - 0.3],
                    [cumulative[i - 1], cumulative[i - 1]],
                    color='gray', linewidth=1, linestyle='--')

        # Baseline and Total bars
        ax.bar(x_pos[0], baseline_value, color=NEOPKPD_COLORS["info"], width=0.6)
        ax.bar(x_pos[-1], cumulative[-1], color=NEOPKPD_COLORS["dark"], width=0.6)

        # Connector to total
        ax.plot([x_pos[-2] + 0.3, x_pos[-1] - 0.3],
                [cumulative[-2], cumulative[-2]],
                color='gray', linewidth=1, linestyle='--')

        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_ylabel(metric_name)
        ax.set_title(title or f"Waterfall - {metric_name} Contributions")
        ax.grid(axis='y', alpha=0.3)

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        labels = ["Baseline"] + params + ["Total"]

        # Plotly waterfall format
        measures = ["absolute"] + ["relative"] * len(params) + ["total"]
        vals = [baseline_value] + list(values) + [None]

        fig = go.Figure(go.Waterfall(
            x=labels,
            measure=measures,
            y=vals,
            decreasing=dict(marker_color=NEOPKPD_COLORS["error"]),
            increasing=dict(marker_color=NEOPKPD_COLORS["success"]),
            totals=dict(marker_color=NEOPKPD_COLORS["dark"]),
            connector=dict(line=dict(color="gray", dash="dash"))
        ))

        fig.update_layout(
            title=title or f"Waterfall - {metric_name} Contributions",
            yaxis_title=metric_name,
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# One-at-a-Time Sensitivity Curves
# ============================================================================

def plot_one_at_a_time(
    result: Any,
    param_name: str,
    metric_names: Optional[List[str]] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Plot one-at-a-time (OFAT) sensitivity curves for a single parameter.

    Shows how multiple output metrics respond to changes in a single
    parameter while all other parameters are held at nominal values.

    Args:
        result: SensitivityResult from run_sensitivity()
        param_name: Name of parameter to plot
        metric_names: Which metrics to show (default: all)
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> from neopkpd import run_sensitivity
        >>> result = run_sensitivity(model, param_ranges, metrics)
        >>> fig = plot_one_at_a_time(result, "CL")
    """
    backend = get_backend()
    theme = get_theme_config()

    # Extract OFAT data for the parameter
    if hasattr(result, 'ofat_results'):
        ofat_data = result.ofat_results.get(param_name)
    else:
        raise ValueError("Result does not contain OFAT data")

    if ofat_data is None:
        raise ValueError(f"No OFAT data for parameter: {param_name}")

    param_values = ofat_data['param_values']
    metrics_data = ofat_data['metrics']

    if metric_names is None:
        metric_names = list(metrics_data.keys())

    colors = get_colors(len(metric_names))

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        for i, metric in enumerate(metric_names):
            if metric in metrics_data:
                values = metrics_data[metric]
                # Normalize to percent change
                baseline_idx = len(values) // 2
                baseline = values[baseline_idx]
                pct_change = (np.array(values) / baseline - 1) * 100

                ax.plot(np.array(param_values) * 100, pct_change,
                        color=colors[i], linewidth=theme["line_width"],
                        marker='o', markersize=4, label=metric)

        ax.axhline(0, color='black', linestyle='--', alpha=0.5)
        ax.axvline(100, color='black', linestyle='--', alpha=0.5)

        ax.set_xlabel(f"{param_name} (% of Nominal)")
        ax.set_ylabel("% Change in Metric")
        ax.set_title(title or f"OFAT Sensitivity - {param_name}")
        ax.legend(loc='best')
        ax.grid(alpha=0.3)

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        fig = go.Figure()

        for i, metric in enumerate(metric_names):
            if metric in metrics_data:
                values = metrics_data[metric]
                baseline_idx = len(values) // 2
                baseline = values[baseline_idx]
                pct_change = (np.array(values) / baseline - 1) * 100

                fig.add_trace(go.Scatter(
                    x=np.array(param_values) * 100,
                    y=pct_change,
                    mode='lines+markers',
                    name=metric,
                    line=dict(color=colors[i])
                ))

        fig.add_hline(y=0, line_dash="dash", line_color="black", opacity=0.5)
        fig.add_vline(x=100, line_dash="dash", line_color="black", opacity=0.5)

        fig.update_layout(
            title=title or f"OFAT Sensitivity - {param_name}",
            xaxis_title=f"{param_name} (% of Nominal)",
            yaxis_title="% Change in Metric",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig
