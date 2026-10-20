"""
NeoPKPD Bootstrap Visualization

Bootstrap analysis visualization including:
- Parameter distribution histograms
- Confidence interval comparisons
- Parameter stability plots
- Correlation matrices from bootstrap
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np

from .backends import _get_plotter, get_backend
from .themes import get_theme_config, get_color, get_colors, NEOPKPD_COLORS


# ============================================================================
# Bootstrap Parameter Distributions
# ============================================================================

def plot_bootstrap_distributions(
    result: Any,
    param_names: Optional[List[str]] = None,
    n_cols: int = 3,
    show_original: bool = True,
    show_ci: bool = True,
    title: Optional[str] = None,
    subplot_size: Tuple[float, float] = (4, 3),
) -> Any:
    """
    Plot histograms of bootstrap parameter estimates.

    Shows the distribution of parameter estimates across bootstrap
    replicates with optional reference lines for original estimate
    and confidence intervals.

    Args:
        result: BootstrapResult from run_bootstrap()
        param_names: Which parameters to plot (default: all)
        n_cols: Number of columns in panel
        show_original: Show original estimate as vertical line
        show_ci: Show CI bounds as vertical lines
        title: Plot title
        subplot_size: Size of each subplot

    Returns:
        Figure object

    Example:
        >>> from neopkpd import run_bootstrap
        >>> boot_result = run_bootstrap(result, data, config)
        >>> fig = plot_bootstrap_distributions(boot_result)
    """
    backend = get_backend()
    theme = get_theme_config()

    estimates = result.theta_estimates  # n_bootstrap x n_params
    n_params = estimates.shape[1]

    all_names = getattr(result, 'theta_names', None) or [f"θ{i+1}" for i in range(n_params)]

    if param_names is not None:
        indices = [all_names.index(n) for n in param_names if n in all_names]
        names = [all_names[i] for i in indices]
    else:
        indices = list(range(n_params))
        names = all_names

    n_plots = len(names)
    n_rows = (n_plots + n_cols - 1) // n_cols
    figsize = (subplot_size[0] * n_cols, subplot_size[1] * n_rows)
    colors = get_colors(n_plots)

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n_plots > 1:
            axes = axes.flatten()
        else:
            axes = [axes]

        for i, (idx, name) in enumerate(zip(indices, names)):
            ax = axes[i]
            data = estimates[:, idx]

            ax.hist(data, bins=30, density=True, alpha=0.7,
                    color=colors[i], edgecolor='white')

            # Original estimate
            if show_original and result.original_estimate:
                orig = result.original_estimate[idx]
                ax.axvline(orig, color='red', linewidth=2, linestyle='-',
                           label=f'Original: {orig:.3g}')

            # CI bounds
            if show_ci:
                ci_lower = result.theta_ci_lower[idx]
                ci_upper = result.theta_ci_upper[idx]
                ax.axvline(ci_lower, color='gray', linewidth=1.5, linestyle='--',
                           label=f'{int(result.ci_level*100)}% CI')
                ax.axvline(ci_upper, color='gray', linewidth=1.5, linestyle='--')

            ax.set_title(name)
            ax.set_xlabel("Estimate")
            ax.set_ylabel("Density")

            if i == 0:
                ax.legend(fontsize='small')

        for i in range(n_plots, len(axes)):
            axes[i].set_visible(False)

        if title:
            fig.suptitle(title)

        fig.tight_layout()
        return fig

    else:
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(rows=n_rows, cols=n_cols, subplot_titles=names)

        for i, (idx, name) in enumerate(zip(indices, names)):
            row = i // n_cols + 1
            col = i % n_cols + 1
            data = estimates[:, idx]

            fig.add_trace(
                go.Histogram(x=data, nbinsx=30, histnorm='probability density',
                             marker_color=colors[i], name=name, showlegend=False),
                row=row, col=col
            )

            if show_original and result.original_estimate:
                fig.add_vline(x=result.original_estimate[idx], line_color="red",
                              line_width=2, row=row, col=col)

            if show_ci:
                fig.add_vline(x=result.theta_ci_lower[idx], line_dash="dash",
                              line_color="gray", row=row, col=col)
                fig.add_vline(x=result.theta_ci_upper[idx], line_dash="dash",
                              line_color="gray", row=row, col=col)

        fig.update_layout(
            title=title or "Bootstrap Parameter Distributions",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Bootstrap CI Comparison
# ============================================================================

def plot_bootstrap_ci(
    result: Any,
    show_bias_corrected: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Plot confidence interval comparison between methods.

    Compares CI estimates from different methods (percentile, BCa)
    and shows original estimate with SE-based CI.

    Args:
        result: BootstrapResult from run_bootstrap()
        show_bias_corrected: Show bias-corrected estimate
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> fig = plot_bootstrap_ci(boot_result, show_bias_corrected=True)
    """
    backend = get_backend()
    theme = get_theme_config()

    n_params = len(result.theta_mean)
    param_names = getattr(result, 'theta_names', None) or [f"θ{i+1}" for i in range(n_params)]

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        y_pos = np.arange(n_params)
        y_offset = 0.15

        # Original estimate with SE-based CI
        ax.scatter(result.original_estimate, y_pos - y_offset, color=get_color(0),
                   s=100, marker='o', label='Original', zorder=3)
        for i, (est, se) in enumerate(zip(result.original_estimate, result.theta_se)):
            ax.plot([est - 1.96 * se, est + 1.96 * se], [i - y_offset, i - y_offset],
                    color=get_color(0), linewidth=2)

        # Bootstrap percentile CI
        ax.scatter(result.theta_mean, y_pos, color=get_color(1), s=100,
                   marker='s', label=f'Bootstrap ({result.ci_method})', zorder=3)
        for i in range(n_params):
            ax.plot([result.theta_ci_lower[i], result.theta_ci_upper[i]], [i, i],
                    color=get_color(1), linewidth=2)

        # Bias-corrected estimate
        if show_bias_corrected and result.bias_corrected is not None:
            ax.scatter(result.bias_corrected, y_pos + y_offset, color=get_color(2),
                       s=100, marker='^', label='Bias-corrected', zorder=3)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(param_names)
        ax.set_xlabel("Estimate")
        ax.set_title(title or "Bootstrap Confidence Intervals")
        ax.legend(loc='upper right')
        ax.grid(axis='x', alpha=0.3)

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        fig = go.Figure()

        y_pos = list(range(n_params))

        # Original estimate
        fig.add_trace(go.Scatter(
            x=result.original_estimate, y=[y - 0.1 for y in y_pos],
            mode='markers',
            marker=dict(size=12, color=get_color(0), symbol='circle'),
            name='Original',
            error_x=dict(type='data', symmetric=True,
                         array=[1.96 * se for se in result.theta_se])
        ))

        # Bootstrap estimate
        fig.add_trace(go.Scatter(
            x=result.theta_mean, y=y_pos,
            mode='markers',
            marker=dict(size=12, color=get_color(1), symbol='square'),
            name=f'Bootstrap ({result.ci_method})',
            error_x=dict(type='data', symmetric=False,
                         array=[u - m for u, m in zip(result.theta_ci_upper, result.theta_mean)],
                         arrayminus=[m - l for m, l in zip(result.theta_mean, result.theta_ci_lower)])
        ))

        # Bias-corrected
        if show_bias_corrected and result.bias_corrected is not None:
            fig.add_trace(go.Scatter(
                x=result.bias_corrected, y=[y + 0.1 for y in y_pos],
                mode='markers',
                marker=dict(size=12, color=get_color(2), symbol='triangle-up'),
                name='Bias-corrected'
            ))

        fig.update_layout(
            title=title or "Bootstrap Confidence Intervals",
            xaxis_title="Estimate",
            yaxis=dict(ticktext=param_names, tickvals=y_pos),
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Bootstrap Stability
# ============================================================================

def plot_bootstrap_stability(
    result: Any,
    param_names: Optional[List[str]] = None,
    n_points: int = 50,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 6),
) -> Any:
    """
    Plot parameter estimate stability across bootstrap replicates.

    Shows cumulative mean and CI as function of number of bootstrap
    samples to assess whether enough replicates were run.

    Args:
        result: BootstrapResult from run_bootstrap()
        param_names: Which parameters to plot
        n_points: Number of evaluation points
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> fig = plot_bootstrap_stability(boot_result)
    """
    backend = get_backend()
    theme = get_theme_config()

    estimates = result.theta_estimates
    n_boot, n_params = estimates.shape

    all_names = getattr(result, 'theta_names', None) or [f"θ{i+1}" for i in range(n_params)]

    if param_names is not None:
        indices = [all_names.index(n) for n in param_names if n in all_names]
        names = [all_names[i] for i in indices]
    else:
        indices = list(range(n_params))
        names = all_names

    # Compute cumulative statistics at n_points evaluation points
    eval_points = np.linspace(10, n_boot, n_points, dtype=int)
    cumulative_mean = np.zeros((len(eval_points), len(indices)))
    cumulative_lower = np.zeros((len(eval_points), len(indices)))
    cumulative_upper = np.zeros((len(eval_points), len(indices)))

    for i, n in enumerate(eval_points):
        for j, idx in enumerate(indices):
            data = estimates[:n, idx]
            cumulative_mean[i, j] = np.mean(data)
            cumulative_lower[i, j] = np.percentile(data, 2.5)
            cumulative_upper[i, j] = np.percentile(data, 97.5)

    colors = get_colors(len(names))

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, len(names), figsize=figsize, squeeze=False)
        axes = axes.flatten()

        for i, (idx, name) in enumerate(zip(indices, names)):
            ax = axes[i]

            ax.fill_between(eval_points, cumulative_lower[:, i], cumulative_upper[:, i],
                            color=colors[i], alpha=0.2)
            ax.plot(eval_points, cumulative_mean[:, i], color=colors[i],
                    linewidth=theme["line_width"])

            # Reference: final estimate
            ax.axhline(result.theta_mean[idx], color='red', linestyle='--',
                       alpha=0.7, label='Final')

            ax.set_xlabel("N Bootstrap")
            ax.set_ylabel("Estimate")
            ax.set_title(name)
            ax.grid(alpha=0.3)

        fig.suptitle(title or "Bootstrap Stability")
        fig.tight_layout()
        return fig

    else:
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(rows=1, cols=len(names), subplot_titles=names)

        for i, (idx, name) in enumerate(zip(indices, names)):
            # CI ribbon
            fig.add_trace(
                go.Scatter(x=list(eval_points), y=list(cumulative_upper[:, i]),
                           mode='lines', line=dict(width=0), showlegend=False),
                row=1, col=i + 1
            )
            fig.add_trace(
                go.Scatter(x=list(eval_points), y=list(cumulative_lower[:, i]),
                           mode='lines', line=dict(width=0),
                           fill='tonexty', fillcolor=f'rgba(0,100,200,0.2)',
                           showlegend=False),
                row=1, col=i + 1
            )

            # Mean line
            fig.add_trace(
                go.Scatter(x=list(eval_points), y=list(cumulative_mean[:, i]),
                           mode='lines', line=dict(color=colors[i]),
                           name=name, showlegend=False),
                row=1, col=i + 1
            )

            # Reference line
            fig.add_hline(y=result.theta_mean[idx], line_dash="dash",
                          line_color="red", row=1, col=i + 1)

        fig.update_layout(
            title=title or "Bootstrap Stability",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Bootstrap Correlation
# ============================================================================

def plot_bootstrap_correlation(
    result: Any,
    annotate: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 8),
    cmap: str = "RdBu_r",
) -> Any:
    """
    Plot correlation matrix of bootstrap parameter estimates.

    Shows inter-parameter correlations computed from bootstrap samples.
    High correlations indicate parameters are not independently estimated.

    Args:
        result: BootstrapResult from run_bootstrap()
        annotate: Show correlation values in cells
        title: Plot title
        figsize: Figure size
        cmap: Colormap name

    Returns:
        Figure object

    Example:
        >>> fig = plot_bootstrap_correlation(boot_result)
    """
    backend = get_backend()

    estimates = result.theta_estimates
    corr = np.corrcoef(estimates.T)

    n_params = estimates.shape[1]
    param_names = getattr(result, 'theta_names', None) or [f"θ{i+1}" for i in range(n_params)]

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        im = ax.imshow(corr, cmap=cmap, vmin=-1, vmax=1)
        plt.colorbar(im, ax=ax)

        if annotate:
            for i in range(n_params):
                for j in range(n_params):
                    text = f"{corr[i, j]:.2f}"
                    color = "white" if abs(corr[i, j]) > 0.5 else "black"
                    ax.text(j, i, text, ha="center", va="center",
                            color=color, fontsize=9)

        ax.set_xticks(np.arange(n_params))
        ax.set_yticks(np.arange(n_params))
        ax.set_xticklabels(param_names, rotation=45, ha="right")
        ax.set_yticklabels(param_names)
        ax.set_title(title or "Bootstrap Parameter Correlations")

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        fig = go.Figure(data=go.Heatmap(
            z=corr,
            x=param_names,
            y=param_names,
            colorscale="RdBu_r",
            zmin=-1, zmax=1,
            text=[[f"{v:.2f}" for v in row] for row in corr] if annotate else None,
            texttemplate="%{text}" if annotate else None
        ))

        fig.update_layout(
            title=title or "Bootstrap Parameter Correlations",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig
