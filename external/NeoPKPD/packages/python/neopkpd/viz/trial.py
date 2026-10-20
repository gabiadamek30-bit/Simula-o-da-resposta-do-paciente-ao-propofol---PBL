"""
NeoPKPD Trial Visualization

Clinical trial visualization functions including:
- Power curves
- Tornado plots (sensitivity analysis)
- Kaplan-Meier survival curves
- Endpoint distributions
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

from .backends import _get_plotter
from .themes import get_theme_config, get_color, get_colors


def plot_power_curve(
    sample_sizes: List[int],
    powers: List[float],
    target_power: float = 0.80,
    achieved_n: Optional[int] = None,
    title: Optional[str] = "Power Curve",
    xlabel: str = "Sample Size (N)",
    ylabel: str = "Power",
    figsize: Tuple[float, float] = (10, 6),
    show_target: bool = True,
) -> Any:
    """
    Plot power as a function of sample size.

    Args:
        sample_sizes: List of sample sizes evaluated
        powers: Corresponding power values
        target_power: Target power level (usually 0.80)
        achieved_n: Sample size achieving target power (will be highlighted)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size
        show_target: Show target power horizontal line

    Returns:
        Figure object

    Example:
        >>> sample_sizes = [20, 40, 60, 80, 100, 120]
        >>> powers = [0.35, 0.55, 0.70, 0.82, 0.89, 0.94]
        >>> fig = plot_power_curve(sample_sizes, powers, target_power=0.80)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Main power curve
    plotter.line_plot(fig, sample_sizes, powers, label="Power",
                      color=get_color(0), linewidth=theme["line_width"] * 1.5)
    plotter.scatter_plot(fig, sample_sizes, powers, color=get_color(0),
                         size=theme["marker_size"] * 6)

    # Target power line
    if show_target:
        backend = plotter.__class__.__name__
        if "Matplotlib" in backend:
            ax = fig["ax"]
            ax.axhline(target_power, color=get_color(2), linestyle="--",
                       linewidth=1.5, alpha=0.7,
                       label=f"Target ({int(target_power*100)}%)")
        else:
            import plotly.graph_objects as go
            fig.add_hline(y=target_power, line_dash="dash",
                          line_color=get_color(2),
                          annotation_text=f"Target ({int(target_power*100)}%)")

    # Highlight achieved sample size
    if achieved_n is not None:
        # Find closest power value
        import numpy as np
        idx = np.argmin(np.abs(np.array(sample_sizes) - achieved_n))
        achieved_power = powers[idx]

        plotter.scatter_plot(fig, [achieved_n], [achieved_power],
                             color=get_color(1), size=theme["marker_size"] * 12,
                             marker="s", label=f"N={achieved_n}")

        backend = plotter.__class__.__name__
        if "Matplotlib" in backend:
            ax = fig["ax"]
            ax.axvline(achieved_n, color=get_color(1), linestyle=":",
                       linewidth=1, alpha=0.5)

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    # Set y-axis limits
    backend = plotter.__class__.__name__
    if "Matplotlib" in backend:
        ax = fig["ax"]
        ax.set_ylim(0, 1.05)
    else:
        fig.update_yaxes(range=[0, 1.05])

    plotter.add_legend(fig)

    return plotter.finalize(fig)


def plot_tornado(
    sensitivity_results: List[Dict[str, Any]],
    baseline_value: float = 0.0,
    sort_by_impact: bool = True,
    title: Optional[str] = "Tornado Plot - Sensitivity Analysis",
    xlabel: str = "Effect on Outcome",
    figsize: Tuple[float, float] = (10, 8),
    show_values: bool = True,
) -> Any:
    """
    Create tornado plot for sensitivity analysis.

    Shows how changes in input parameters affect the outcome.

    Args:
        sensitivity_results: List of dicts with keys:
            - 'parameter': Parameter name
            - 'low_value': Outcome when parameter is low
            - 'high_value': Outcome when parameter is high
            - 'low_label': Optional label for low value
            - 'high_label': Optional label for high value
        baseline_value: Baseline outcome value (center line)
        sort_by_impact: Sort bars by total impact
        title: Plot title
        xlabel: X-axis label
        figsize: Figure size
        show_values: Show numeric values on bars

    Returns:
        Figure object

    Example:
        >>> sensitivity = [
        ...     {"parameter": "CL", "low_value": 8.5, "high_value": 11.2},
        ...     {"parameter": "V", "low_value": 9.2, "high_value": 10.5},
        ...     {"parameter": "ka", "low_value": 9.8, "high_value": 10.1},
        ... ]
        >>> fig = plot_tornado(sensitivity, baseline_value=10.0)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    # Calculate impacts
    results = []
    for r in sensitivity_results:
        low_impact = r["low_value"] - baseline_value
        high_impact = r["high_value"] - baseline_value
        total_impact = abs(high_impact - low_impact)
        results.append({
            **r,
            "low_impact": low_impact,
            "high_impact": high_impact,
            "total_impact": total_impact
        })

    # Sort by impact if requested
    if sort_by_impact:
        results = sorted(results, key=lambda x: x["total_impact"], reverse=True)

    n = len(results)
    y_positions = list(range(n))
    parameters = [r["parameter"] for r in results]

    backend = plotter.__class__.__name__

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        for i, r in enumerate(results):
            low_impact = r["low_impact"]
            high_impact = r["high_impact"]

            # Draw bar from baseline to low value (extends left or right)
            if low_impact < 0:
                ax.barh(i, low_impact, left=baseline_value, height=0.6,
                        color=get_color(0), alpha=0.8, label="Low" if i == 0 else "")
            else:
                ax.barh(i, low_impact, left=baseline_value, height=0.6,
                        color=get_color(0), alpha=0.8, label="Low" if i == 0 else "")

            # Draw bar from baseline to high value
            if high_impact > 0:
                ax.barh(i, high_impact, left=baseline_value, height=0.6,
                        color=get_color(1), alpha=0.8, label="High" if i == 0 else "")
            else:
                ax.barh(i, high_impact, left=baseline_value, height=0.6,
                        color=get_color(1), alpha=0.8, label="High" if i == 0 else "")

            # Show values
            if show_values:
                low_label = r.get("low_label", f"{r['low_value']:.2f}")
                high_label = r.get("high_label", f"{r['high_value']:.2f}")

                # Position labels at ends of bars
                if low_impact < 0:
                    ax.text(baseline_value + low_impact - 0.02, i, low_label,
                            va="center", ha="right", fontsize=8)
                else:
                    ax.text(baseline_value + low_impact + 0.02, i, low_label,
                            va="center", ha="left", fontsize=8)

                if high_impact > 0:
                    ax.text(baseline_value + high_impact + 0.02, i, high_label,
                            va="center", ha="left", fontsize=8)
                else:
                    ax.text(baseline_value + high_impact - 0.02, i, high_label,
                            va="center", ha="right", fontsize=8)

        # Baseline line
        ax.axvline(baseline_value, color="black", linewidth=1.5, zorder=3)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(parameters)
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.legend(loc="lower right")

        # Add some padding
        xlim = ax.get_xlim()
        padding = (xlim[1] - xlim[0]) * 0.15
        ax.set_xlim(xlim[0] - padding, xlim[1] + padding)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        return fig

    else:  # Plotly
        import plotly.graph_objects as go

        fig = go.Figure()

        # Low values (left bars)
        fig.add_trace(go.Bar(
            y=parameters,
            x=[r["low_impact"] for r in results],
            orientation="h",
            name="Low",
            marker_color=get_color(0),
            base=baseline_value
        ))

        # High values (right bars)
        fig.add_trace(go.Bar(
            y=parameters,
            x=[r["high_impact"] for r in results],
            orientation="h",
            name="High",
            marker_color=get_color(1),
            base=baseline_value
        ))

        # Baseline line
        fig.add_vline(x=baseline_value, line_width=2, line_color="black")

        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            barmode="overlay",
            height=figsize[1] * 100,
            width=figsize[0] * 100
        )

        return fig


def plot_kaplan_meier(
    time_to_event: List[float],
    event_occurred: List[bool],
    groups: Optional[List[str]] = None,
    title: Optional[str] = "Kaplan-Meier Survival Curve",
    xlabel: str = "Time",
    ylabel: str = "Survival Probability",
    figsize: Tuple[float, float] = (10, 6),
    show_censored: bool = True,
    show_ci: bool = True,
) -> Any:
    """
    Create Kaplan-Meier survival curve.

    Args:
        time_to_event: Time to event or censoring for each subject
        event_occurred: Whether the event occurred (True) or was censored (False)
        groups: Optional group labels for stratified analysis
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size
        show_censored: Show censored observations as tick marks
        show_ci: Show confidence intervals

    Returns:
        Figure object

    Example:
        >>> times = [10, 15, 20, 25, 30, 35, 40, 45]
        >>> events = [True, False, True, True, False, True, False, True]
        >>> fig = plot_kaplan_meier(times, events)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    def compute_km_curve(times, events):
        """Compute Kaplan-Meier survival estimates."""
        # Sort by time
        sorted_indices = np.argsort(times)
        times = np.array(times)[sorted_indices]
        events = np.array(events)[sorted_indices]

        unique_times = np.unique(times[events])
        n_at_risk = []
        n_events = []
        survival = []

        n = len(times)
        current_survival = 1.0

        # Add time 0
        km_times = [0.0]
        km_survival = [1.0]

        for t in unique_times:
            at_risk = np.sum(times >= t)
            d = np.sum((times == t) & events)

            if at_risk > 0:
                current_survival *= (1 - d / at_risk)

            km_times.append(t)
            km_survival.append(current_survival)

        # Extend to max time
        km_times.append(max(times))
        km_survival.append(km_survival[-1])

        return km_times, km_survival

    def compute_km_ci(times, events, km_times, km_survival, alpha=0.05):
        """Compute Greenwood confidence intervals."""
        import numpy as np
        from scipy import stats

        n = len(times)
        times = np.array(times)
        events = np.array(events)

        variance = []
        current_var = 0.0
        prev_survival = 1.0

        lower = [1.0]
        upper = [1.0]

        for i, t in enumerate(km_times[1:-1]):
            at_risk = np.sum(times >= t)
            d = np.sum((times == t) & events)

            if at_risk > 0 and at_risk - d > 0:
                current_var += d / (at_risk * (at_risk - d))

            s = km_survival[i + 1]
            if s > 0 and s < 1:
                se = s * np.sqrt(current_var)
                z = stats.norm.ppf(1 - alpha / 2)
                lower.append(max(0, s - z * se))
                upper.append(min(1, s + z * se))
            else:
                lower.append(s)
                upper.append(s)

        lower.append(lower[-1])
        upper.append(upper[-1])

        return lower, upper

    fig = plotter.create_figure(figsize=figsize, title=title)

    if groups is None:
        # Single group
        km_times, km_survival = compute_km_curve(time_to_event, event_occurred)

        # Step function
        backend = plotter.__class__.__name__
        if "Matplotlib" in backend:
            ax = fig["ax"]
            ax.step(km_times, km_survival, where="post", color=get_color(0),
                    linewidth=theme["line_width"], label="Survival")

            if show_ci:
                lower, upper = compute_km_ci(time_to_event, event_occurred,
                                             km_times, km_survival)
                ax.fill_between(km_times, lower, upper, step="post",
                                alpha=theme["alpha_ribbon"], color=get_color(0))

            if show_censored:
                censored_times = [t for t, e in zip(time_to_event, event_occurred) if not e]
                # Find survival at censored times
                censored_surv = []
                for ct in censored_times:
                    for i, kt in enumerate(km_times):
                        if kt >= ct:
                            censored_surv.append(km_survival[max(0, i-1)])
                            break
                    else:
                        censored_surv.append(km_survival[-1])

                ax.scatter(censored_times, censored_surv, marker="|", s=50,
                           color=get_color(0), zorder=3, label="Censored")
        else:
            import plotly.graph_objects as go
            fig.add_trace(go.Scatter(
                x=km_times, y=km_survival, mode="lines",
                line=dict(shape="hv", color=get_color(0)),
                name="Survival"
            ))

    else:
        # Multiple groups
        unique_groups = list(dict.fromkeys(groups))
        colors = get_colors(len(unique_groups))

        for gi, group in enumerate(unique_groups):
            group_mask = [g == group for g in groups]
            group_times = [t for t, m in zip(time_to_event, group_mask) if m]
            group_events = [e for e, m in zip(event_occurred, group_mask) if m]

            km_times, km_survival = compute_km_curve(group_times, group_events)

            backend = plotter.__class__.__name__
            if "Matplotlib" in backend:
                ax = fig["ax"]
                ax.step(km_times, km_survival, where="post", color=colors[gi],
                        linewidth=theme["line_width"], label=group)
            else:
                import plotly.graph_objects as go
                fig.add_trace(go.Scatter(
                    x=km_times, y=km_survival, mode="lines",
                    line=dict(shape="hv", color=colors[gi]),
                    name=group
                ))

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    # Set y-axis limits
    backend = plotter.__class__.__name__
    if "Matplotlib" in backend:
        ax = fig["ax"]
        ax.set_ylim(0, 1.05)
    else:
        fig.update_yaxes(range=[0, 1.05])

    plotter.add_legend(fig)

    return plotter.finalize(fig)


def plot_endpoint_distribution(
    trial_result: Dict[str, Any],
    endpoint: str,
    by_arm: bool = True,
    plot_type: str = "histogram",
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: str = "Frequency",
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Plot distribution of endpoint values from trial simulation.

    Args:
        trial_result: Trial simulation result dictionary
        endpoint: Endpoint name to plot
        by_arm: Separate distributions by treatment arm
        plot_type: "histogram", "kde", or "boxplot"
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> trial_result = simulate_trial(...)
        >>> fig = plot_endpoint_distribution(trial_result, "auc", by_arm=True)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    if title is None:
        title = f"Distribution of {endpoint}"
    if xlabel is None:
        xlabel = endpoint

    backend = plotter.__class__.__name__

    if by_arm and "arms" in trial_result:
        arms = trial_result["arms"]
        arm_names = list(arms.keys())
        n_arms = len(arm_names)
        colors = get_colors(n_arms)

        if "Matplotlib" in backend:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=figsize)

            if plot_type == "boxplot":
                data = []
                for arm_name in arm_names:
                    arm_data = arms[arm_name]
                    if "endpoint_values" in arm_data and endpoint in arm_data["endpoint_values"]:
                        data.append(arm_data["endpoint_values"][endpoint])
                    else:
                        data.append([])

                bp = ax.boxplot(data, labels=arm_names, patch_artist=True)
                for patch, color in zip(bp["boxes"], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
            else:
                for i, arm_name in enumerate(arm_names):
                    arm_data = arms[arm_name]
                    if "endpoint_values" in arm_data and endpoint in arm_data["endpoint_values"]:
                        values = arm_data["endpoint_values"][endpoint]
                        ax.hist(values, bins=30, alpha=0.6, color=colors[i],
                                label=arm_name, edgecolor="white")

                ax.legend()

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            fig.tight_layout()
            return fig

        else:  # Plotly
            import plotly.graph_objects as go

            fig = go.Figure()

            for i, arm_name in enumerate(arm_names):
                arm_data = arms[arm_name]
                if "endpoint_values" in arm_data and endpoint in arm_data["endpoint_values"]:
                    values = arm_data["endpoint_values"][endpoint]

                    if plot_type == "boxplot":
                        fig.add_trace(go.Box(y=values, name=arm_name,
                                             marker_color=colors[i]))
                    else:
                        fig.add_trace(go.Histogram(x=values, name=arm_name,
                                                   marker_color=colors[i],
                                                   opacity=0.6))

            fig.update_layout(
                title=title,
                xaxis_title=xlabel,
                yaxis_title=ylabel,
                barmode="overlay",
                height=figsize[1] * 100,
                width=figsize[0] * 100
            )

            return fig

    else:
        # Single distribution
        if "endpoint_values" in trial_result and endpoint in trial_result["endpoint_values"]:
            values = trial_result["endpoint_values"][endpoint]
        else:
            values = []

        fig = plotter.create_figure(figsize=figsize, title=title)

        if "Matplotlib" in backend:
            ax = fig["ax"]
            ax.hist(values, bins=30, color=get_color(0), alpha=0.7,
                    edgecolor="white")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            fig["fig"].tight_layout()
            return fig["fig"]

        else:
            import plotly.graph_objects as go
            fig.add_trace(go.Histogram(x=values, marker_color=get_color(0)))
            fig.update_layout(xaxis_title=xlabel, yaxis_title=ylabel)
            return fig
