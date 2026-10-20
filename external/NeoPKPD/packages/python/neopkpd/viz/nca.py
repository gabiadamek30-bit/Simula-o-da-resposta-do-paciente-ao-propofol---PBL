"""
NeoPKPD NCA Visualization

Non-Compartmental Analysis visualization functions including:
- Lambda-z fit plots
- AUC visualization
- Dose proportionality plots
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from .backends import _get_plotter
from .themes import get_theme_config, get_color, get_colors


def plot_lambda_z_fit(
    nca_result: Any,
    times: List[float],
    concentrations: List[float],
    show_excluded: bool = True,
    log_scale: bool = True,
    title: Optional[str] = "Lambda-z Regression Fit",
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Plot terminal phase lambda-z regression fit.

    Shows the concentration-time data with the lambda-z regression line
    and highlights the points used in the regression.

    Args:
        nca_result: NCA result object (with lambda_z, intercept, points_used)
        times: Time points
        concentrations: Concentration values
        show_excluded: Show excluded points in different style
        log_scale: Use log scale for y-axis (recommended for lambda-z)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> result = nca.run_nca(t, c, dose=100)
        >>> fig = plot_lambda_z_fit(result, t, c)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Get lambda-z info
    if hasattr(nca_result, 'lambda_z_result'):
        # Python NCAResult
        lz_result = nca_result.lambda_z_result if hasattr(nca_result, 'lambda_z_result') else nca_result
        lambda_z = nca_result.lambda_z
    elif isinstance(nca_result, dict):
        # Dict result
        lambda_z = nca_result.get("lambda_z")
        lz_result = nca_result
    else:
        lambda_z = getattr(nca_result, 'lambda_z', None)
        lz_result = nca_result

    # Determine which points were used (assuming terminal phase after tmax)
    import numpy as np
    t = np.array(times)
    c = np.array(concentrations)

    # Find Cmax index
    tmax_idx = np.argmax(c)

    # Points after Cmax with positive concentration
    terminal_mask = (np.arange(len(t)) > tmax_idx) & (c > 0)
    used_mask = np.zeros(len(t), dtype=bool)

    if hasattr(lz_result, 'points_used'):
        for idx in lz_result.points_used:
            if 0 <= idx < len(used_mask):
                used_mask[idx] = True
    elif isinstance(lz_result, dict) and 'points_used' in lz_result:
        for idx in lz_result['points_used']:
            if 0 <= idx < len(used_mask):
                used_mask[idx] = True
    else:
        # Default: use last 3+ points after Cmax
        used_mask = terminal_mask

    # Plot all data points
    if show_excluded:
        excluded_mask = ~used_mask
        if np.any(excluded_mask):
            plotter.scatter_plot(fig, list(t[excluded_mask]), list(c[excluded_mask]),
                                 color=get_color(2), marker="o",
                                 size=theme["marker_size"] * 4, alpha=0.5,
                                 label="Excluded")

    # Plot used points
    plotter.scatter_plot(fig, list(t[used_mask]), list(c[used_mask]),
                         color=get_color(0), marker="o",
                         size=theme["marker_size"] * 6,
                         label="Used for regression")

    # Plot regression line if lambda_z available
    if lambda_z is not None and lambda_z > 0:
        intercept = getattr(lz_result, 'intercept', None)
        if intercept is None and isinstance(lz_result, dict):
            intercept = lz_result.get('intercept')

        if intercept is not None:
            # Create regression line
            t_fit = t[used_mask]
            if len(t_fit) > 0:
                t_line = np.linspace(min(t_fit), max(t_fit), 50)
                c_line = np.exp(intercept - lambda_z * t_line)
                plotter.line_plot(fig, list(t_line), list(c_line),
                                  color=get_color(1), linestyle="--",
                                  linewidth=theme["line_width"],
                                  label=f"Fit (λz={lambda_z:.4f})")

    # Connect all points with a line
    plotter.line_plot(fig, list(t), list(c), color=get_color(0),
                      linewidth=0.5, alpha=0.5)

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig)

    return plotter.finalize(fig)


def plot_auc_visualization(
    times: List[float],
    concentrations: List[float],
    nca_result: Optional[Any] = None,
    show_extrapolation: bool = True,
    log_scale: bool = False,
    title: Optional[str] = "AUC Visualization",
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Visualize AUC with shaded area under the curve.

    Args:
        times: Time points
        concentrations: Concentration values
        nca_result: Optional NCA result for extrapolation info
        show_extrapolation: Show extrapolated AUC portion
        log_scale: Use log scale for y-axis
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> fig = plot_auc_visualization(t, c, nca_result=result, show_extrapolation=True)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np
    t = np.array(times)
    c = np.array(concentrations)

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Fill AUC area
    zeros = [0.0] * len(t)
    plotter.fill_between(fig, list(t), zeros, list(c),
                         color=get_color(0), alpha=0.3,
                         label="AUC0-t")

    # Plot concentration line
    plotter.line_plot(fig, list(t), list(c), color=get_color(0),
                      linewidth=theme["line_width"])

    # Plot data points
    plotter.scatter_plot(fig, list(t), list(c), color=get_color(0),
                         size=theme["marker_size"] * 5)

    # Show extrapolation if available
    if show_extrapolation and nca_result is not None:
        lambda_z = getattr(nca_result, 'lambda_z', None)
        if lambda_z is None and isinstance(nca_result, dict):
            lambda_z = nca_result.get('lambda_z')

        if lambda_z is not None and lambda_z > 0:
            clast = c[-1] if c[-1] > 0 else c[c > 0][-1] if np.any(c > 0) else 0
            tlast = t[len(c) - 1]

            if clast > 0:
                # Extrapolation line
                t_extra = np.linspace(tlast, tlast + 3 / lambda_z, 50)
                c_extra = clast * np.exp(-lambda_z * (t_extra - tlast))

                # Extrapolation area
                zeros_extra = [0.0] * len(t_extra)
                plotter.fill_between(fig, list(t_extra), zeros_extra, list(c_extra),
                                     color=get_color(1), alpha=0.2,
                                     label="AUC extrapolated")

                plotter.line_plot(fig, list(t_extra), list(c_extra),
                                  color=get_color(1), linestyle="--",
                                  linewidth=theme["line_width"])

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig)

    return plotter.finalize(fig)


def plot_dose_proportionality(
    nca_results: List[Any],
    doses: List[float],
    metric: str = "auc_0_inf",
    normalize: bool = True,
    show_fit: bool = True,
    title: Optional[str] = "Dose Proportionality",
    xlabel: str = "Dose",
    ylabel: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Plot dose proportionality for a PK parameter.

    Args:
        nca_results: List of NCA results at different doses
        doses: List of dose values
        metric: PK parameter to plot (e.g., "auc_0_inf", "cmax")
        normalize: Plot dose-normalized values
        show_fit: Show power model fit line
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label (auto-generated if None)
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> results = [nca.run_nca(t, c, dose=d) for d in [25, 50, 100, 200]]
        >>> fig = plot_dose_proportionality(results, [25, 50, 100, 200])
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    # Extract metric values
    values = []
    for result in nca_results:
        if hasattr(result, metric):
            val = getattr(result, metric)
        elif isinstance(result, dict):
            val = result.get(metric)
        else:
            val = None

        if val is not None:
            values.append(float(val))
        else:
            values.append(np.nan)

    doses = np.array(doses)
    values = np.array(values)

    # Remove NaN
    valid = ~np.isnan(values)
    doses_valid = doses[valid]
    values_valid = values[valid]

    if normalize:
        values_plot = values_valid / doses_valid
        if ylabel is None:
            ylabel = f"{metric} / Dose"
    else:
        values_plot = values_valid
        if ylabel is None:
            ylabel = metric

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Plot data points
    plotter.scatter_plot(fig, list(doses_valid), list(values_plot),
                         color=get_color(0), size=theme["marker_size"] * 8,
                         label="Observed")

    # Add power model fit line
    if show_fit and len(doses_valid) >= 2:
        # Log-log regression
        log_dose = np.log(doses_valid)
        log_values = np.log(values_valid)

        # Simple linear regression on log scale
        n = len(log_dose)
        mean_x = np.mean(log_dose)
        mean_y = np.mean(log_values)
        slope = np.sum((log_dose - mean_x) * (log_values - mean_y)) / np.sum((log_dose - mean_x)**2)
        intercept = mean_y - slope * mean_x

        # Fitted line
        dose_line = np.linspace(min(doses_valid), max(doses_valid), 50)
        values_line = np.exp(intercept + slope * np.log(dose_line))

        if normalize:
            values_line = values_line / dose_line

        plotter.line_plot(fig, list(dose_line), list(values_line),
                          color=get_color(1), linestyle="--",
                          linewidth=theme["line_width"],
                          label=f"Power fit (β={slope:.2f})")

        # Linear proportionality line (β=1)
        values_linear = np.exp(intercept) * dose_line
        if normalize:
            values_linear = values_linear / dose_line

        plotter.line_plot(fig, list(dose_line), list(values_linear),
                          color=get_color(2), linestyle=":",
                          linewidth=theme["line_width"],
                          label="Linear (β=1)")

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)
    plotter.add_legend(fig)

    return plotter.finalize(fig)
