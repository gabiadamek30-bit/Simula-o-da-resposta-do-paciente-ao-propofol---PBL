"""
NeoPKPD PKPD Visualization

Pharmacokinetic-Pharmacodynamic visualization functions including:
- Effect-concentration plots
- Hysteresis plots
- Dose-response plots
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from .backends import _get_plotter
from .themes import get_theme_config, get_color, get_colors


def plot_effect_conc(
    result: Dict[str, Any],
    conc_observation: str = "conc",
    effect_observation: str = "effect",
    color_by_time: bool = True,
    show_line: bool = True,
    title: Optional[str] = "Effect vs Concentration",
    xlabel: str = "Concentration",
    ylabel: str = "Effect",
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Plot effect vs concentration relationship.

    Args:
        result: Simulation result dictionary with both PK and PD observations
        conc_observation: Key for concentration data
        effect_observation: Key for effect data
        color_by_time: Color points by time to show temporal progression
        show_line: Connect points with line
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> result = neopkpd.simulate_pkpd_sigmoid_emax(...)
        >>> fig = plot_effect_conc(result, color_by_time=True)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    t = result["t"]
    c = result["observations"][conc_observation]
    e = result["observations"][effect_observation]

    fig = plotter.create_figure(figsize=figsize, title=title)

    if show_line:
        plotter.line_plot(fig, c, e, color=get_color(0),
                          linewidth=theme["line_width"] * 0.8, alpha=0.5)

    if color_by_time:
        # Use matplotlib colormap directly for time coloring
        backend = plotter.__class__.__name__
        if "Matplotlib" in backend:
            import matplotlib.pyplot as plt
            import numpy as np
            ax = fig["ax"]
            scatter = ax.scatter(c, e, c=t, cmap="viridis", s=theme["marker_size"] * 10)
            plt.colorbar(scatter, ax=ax, label="Time")
        else:
            # Plotly - use color scale
            import plotly.graph_objects as go
            fig.add_trace(go.Scatter(
                x=c, y=e, mode="markers",
                marker=dict(color=t, colorscale="Viridis", showscale=True,
                            colorbar=dict(title="Time")),
            ))
            return fig
    else:
        plotter.scatter_plot(fig, c, e, color=get_color(0),
                             size=theme["marker_size"] * 6)

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    return plotter.finalize(fig)


def plot_hysteresis(
    result: Dict[str, Any],
    conc_observation: str = "conc",
    effect_observation: str = "effect",
    arrow_frequency: int = 5,
    title: Optional[str] = "Hysteresis Plot",
    xlabel: str = "Concentration",
    ylabel: str = "Effect",
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Plot hysteresis loop showing temporal lag between concentration and effect.

    Args:
        result: Simulation result dictionary
        conc_observation: Key for concentration data
        effect_observation: Key for effect data
        arrow_frequency: Add arrow every N points to show direction
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> result = neopkpd.simulate_pkpd_biophase_equilibration(...)
        >>> fig = plot_hysteresis(result, arrow_frequency=10)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    t = result["t"]
    c = result["observations"][conc_observation]
    e = result["observations"][effect_observation]

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Plot the hysteresis loop
    plotter.line_plot(fig, c, e, color=get_color(0),
                      linewidth=theme["line_width"])

    # Add arrows to show direction
    backend = plotter.__class__.__name__
    if "Matplotlib" in backend:
        ax = fig["ax"]
        for i in range(arrow_frequency, len(c) - 1, arrow_frequency):
            dx = c[i+1] - c[i]
            dy = e[i+1] - e[i]
            if abs(dx) > 1e-10 or abs(dy) > 1e-10:
                ax.annotate("", xy=(c[i+1], e[i+1]), xytext=(c[i], e[i]),
                            arrowprops=dict(arrowstyle="->", color=get_color(1),
                                            lw=1.5))

    # Mark start and end
    plotter.scatter_plot(fig, [c[0]], [e[0]], color=get_color(2),
                         marker="o", size=theme["marker_size"] * 10,
                         label="Start")
    plotter.scatter_plot(fig, [c[-1]], [e[-1]], color=get_color(3),
                         marker="s", size=theme["marker_size"] * 10,
                         label="End")

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)
    plotter.add_legend(fig)

    return plotter.finalize(fig)


def plot_dose_response(
    dose_results: List[Dict[str, Any]],
    doses: List[float],
    effect_observation: str = "effect",
    metric: str = "max",
    fit_model: Optional[str] = "emax",
    title: Optional[str] = "Dose-Response Curve",
    xlabel: str = "Dose",
    ylabel: str = "Effect",
    figsize: Tuple[float, float] = (10, 6),
    log_dose: bool = True,
) -> Any:
    """
    Plot dose-response relationship.

    Args:
        dose_results: List of simulation results at different doses
        doses: List of dose values
        effect_observation: Key for effect data
        metric: How to summarize effect ("max", "auc", "final")
        fit_model: Model to fit ("emax", "sigmoid", None)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size
        log_dose: Use log scale for dose axis

    Returns:
        Figure object

    Example:
        >>> results = [simulate_pkpd_emax(dose=d) for d in [10, 25, 50, 100, 200]]
        >>> fig = plot_dose_response(results, [10, 25, 50, 100, 200])
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    # Extract effect metric for each dose
    effects = []
    for result in dose_results:
        e = result["observations"][effect_observation]
        if metric == "max":
            effects.append(max(e))
        elif metric == "min":
            effects.append(min(e))
        elif metric == "final":
            effects.append(e[-1])
        elif metric == "auc":
            t = result["t"]
            effects.append(np.trapz(e, t))
        else:
            effects.append(max(e))

    doses = np.array(doses)
    effects = np.array(effects)

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Plot data points
    plotter.scatter_plot(fig, list(doses), list(effects), color=get_color(0),
                         size=theme["marker_size"] * 8, label="Observed")

    # Fit Emax model if requested
    if fit_model in ["emax", "sigmoid"]:
        try:
            from scipy.optimize import curve_fit

            if fit_model == "emax":
                def emax_func(d, emax, ec50):
                    return emax * d / (ec50 + d)
                p0 = [max(effects), np.median(doses)]
            else:  # sigmoid
                def emax_func(d, emax, ec50, gamma):
                    return emax * d**gamma / (ec50**gamma + d**gamma)
                p0 = [max(effects), np.median(doses), 1.0]

            popt, _ = curve_fit(emax_func, doses, effects, p0=p0, maxfev=2000)

            # Plot fit line
            dose_line = np.linspace(min(doses), max(doses), 100)
            effect_line = emax_func(dose_line, *popt)

            label = f"Emax fit (EC50={popt[1]:.1f})"
            if fit_model == "sigmoid" and len(popt) > 2:
                label = f"Sigmoid fit (EC50={popt[1]:.1f}, γ={popt[2]:.2f})"

            plotter.line_plot(fig, list(dose_line), list(effect_line),
                              color=get_color(1), linestyle="--",
                              linewidth=theme["line_width"], label=label)

        except (ImportError, RuntimeError):
            pass  # Skip fitting if scipy not available or fit fails

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_dose:
        plotter.set_log_scale(fig, x=True)

    plotter.add_legend(fig)

    return plotter.finalize(fig)
