"""
NeoPKPD Estimation Visualization

Parameter estimation diagnostics visualization including:
- Parameter estimates forest plots
- Omega correlation heatmaps
- Convergence diagnostics
- Eta shrinkage plots
- Residual distributions
- Individual EBE plots
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np

from .backends import _get_plotter, get_backend
from .themes import get_theme_config, get_color, get_colors, NEOPKPD_COLORS


# ============================================================================
# Parameter Estimates Forest Plot
# ============================================================================

def plot_parameter_estimates(
    result: Any,
    include_omega: bool = False,
    show_ci: bool = True,
    reference_values: Optional[Dict[str, float]] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Create forest plot of parameter estimates with confidence intervals.

    Displays fixed effect (theta) estimates with their CIs. Optionally
    includes omega (random effect variance) estimates.

    Args:
        result: EstimationResult from estimate()
        include_omega: Include omega diagonal elements
        show_ci: Show confidence intervals
        reference_values: Optional dict of reference values to mark
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> result = estimate(data, "OneCompIVBolus", config, grid)
        >>> fig = plot_parameter_estimates(result, include_omega=True)
    """
    theme = get_theme_config()
    backend = get_backend()

    # Extract parameter data
    theta = result.theta
    theta_se = result.theta_se
    theta_names = getattr(result, 'theta_names', None) or [f"θ{i+1}" for i in range(len(theta))]

    if include_omega:
        omega_diag = [result.omega[i][i] for i in range(len(result.omega))]
        omega_names = getattr(result, 'omega_names', None) or [f"ω²{i+1}" for i in range(len(omega_diag))]
        omega_se = None
        if result.omega_se:
            omega_se = [result.omega_se[i][i] for i in range(len(result.omega_se))]

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        # Combine all parameters
        params = list(theta)
        names = list(theta_names)
        ses = list(theta_se) if theta_se else [None] * len(theta)

        if include_omega:
            params.extend(omega_diag)
            names.extend(omega_names)
            if omega_se:
                ses.extend(omega_se)
            else:
                ses.extend([None] * len(omega_diag))

        n_params = len(params)
        y_pos = np.arange(n_params)

        # Plot point estimates
        ax.scatter(params, y_pos, color=get_color(0), s=100, zorder=3)

        # Plot CIs if available
        if show_ci and any(s is not None for s in ses):
            for i, (param, se) in enumerate(zip(params, ses)):
                if se is not None:
                    ci_lower = param - 1.96 * se
                    ci_upper = param + 1.96 * se
                    ax.plot([ci_lower, ci_upper], [i, i], color=get_color(0),
                            linewidth=2, zorder=2)
                    ax.plot([ci_lower, ci_lower], [i - 0.1, i + 0.1],
                            color=get_color(0), linewidth=2)
                    ax.plot([ci_upper, ci_upper], [i - 0.1, i + 0.1],
                            color=get_color(0), linewidth=2)

        # Plot reference values if provided
        if reference_values:
            for name, ref in reference_values.items():
                if name in names:
                    idx = names.index(name)
                    ax.scatter([ref], [idx], color=get_color(1), s=80,
                               marker='D', zorder=4, label='Reference' if idx == 0 else None)

        # Add vertical line at zero for omega
        if include_omega:
            ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(names)
        ax.set_xlabel("Estimate")
        ax.set_title(title or "Parameter Estimates")
        ax.grid(axis='x', alpha=0.3)

        if reference_values:
            ax.legend()

        fig.tight_layout()
        return fig

    else:
        # Plotly implementation
        import plotly.graph_objects as go

        params = list(theta)
        names = list(theta_names)
        ses = list(theta_se) if theta_se else [None] * len(theta)

        if include_omega:
            params.extend(omega_diag)
            names.extend(omega_names)
            if omega_se:
                ses.extend(omega_se)
            else:
                ses.extend([None] * len(omega_diag))

        fig = go.Figure()

        # Error bars for CIs
        error_x = dict(type='data', symmetric=True,
                       array=[1.96 * s if s else 0 for s in ses])

        fig.add_trace(go.Scatter(
            x=params, y=names,
            mode='markers',
            marker=dict(size=12, color=get_color(0)),
            error_x=error_x if show_ci else None,
            name='Estimate'
        ))

        # Reference values
        if reference_values:
            ref_x = []
            ref_y = []
            for name, ref in reference_values.items():
                if name in names:
                    ref_x.append(ref)
                    ref_y.append(name)
            if ref_x:
                fig.add_trace(go.Scatter(
                    x=ref_x, y=ref_y,
                    mode='markers',
                    marker=dict(size=10, color=get_color(1), symbol='diamond'),
                    name='Reference'
                ))

        fig.update_layout(
            title=title or "Parameter Estimates",
            xaxis_title="Estimate",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Omega Correlation Matrix Heatmap
# ============================================================================

def plot_omega_matrix(
    result: Any,
    show_correlation: bool = True,
    annotate: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 8),
    cmap: str = "RdBu_r",
) -> Any:
    """
    Plot heatmap of omega (random effect covariance/correlation) matrix.

    Args:
        result: EstimationResult from estimate()
        show_correlation: Show correlation matrix (vs covariance)
        annotate: Show values in cells
        title: Plot title
        figsize: Figure size
        cmap: Colormap name

    Returns:
        Figure object

    Example:
        >>> fig = plot_omega_matrix(result, show_correlation=True)
    """
    backend = get_backend()

    if show_correlation:
        matrix = np.array(result.omega_corr)
        mat_title = "Random Effects Correlation Matrix"
        vmin, vmax = -1, 1
    else:
        matrix = np.array(result.omega)
        mat_title = "Random Effects Covariance Matrix"
        vmin, vmax = None, None

    omega_names = getattr(result, 'omega_names', None) or [f"η{i+1}" for i in range(len(matrix))]

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax)

        # Add colorbar
        plt.colorbar(im, ax=ax)

        # Add annotations
        if annotate:
            for i in range(len(matrix)):
                for j in range(len(matrix)):
                    text = f"{matrix[i, j]:.3f}"
                    color = "white" if abs(matrix[i, j]) > 0.5 else "black"
                    ax.text(j, i, text, ha="center", va="center", color=color)

        ax.set_xticks(np.arange(len(omega_names)))
        ax.set_yticks(np.arange(len(omega_names)))
        ax.set_xticklabels(omega_names, rotation=45, ha="right")
        ax.set_yticklabels(omega_names)
        ax.set_title(title or mat_title)

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=omega_names,
            y=omega_names,
            colorscale="RdBu_r",
            zmin=vmin,
            zmax=vmax,
            text=[[f"{v:.3f}" for v in row] for row in matrix] if annotate else None,
            texttemplate="%{text}" if annotate else None
        ))

        fig.update_layout(
            title=title or mat_title,
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Convergence Plot
# ============================================================================

def plot_convergence(
    result: Any,
    show_gradient: bool = False,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 5),
) -> Any:
    """
    Plot convergence diagnostics (OFV vs iteration).

    Shows how the objective function value evolved during estimation.
    For SAEM, shows burn-in and estimation phases.

    Args:
        result: EstimationResult with convergence history
        show_gradient: Show gradient norm subplot
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> fig = plot_convergence(result)
    """
    backend = get_backend()
    theme = get_theme_config()

    # Check if convergence history available
    if not hasattr(result, 'convergence_history') or result.convergence_history is None:
        # Create a simple summary plot
        if backend == "matplotlib":
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=figsize)

            # Show final OFV
            ax.bar(['Final OFV'], [result.ofv], color=get_color(0))
            ax.set_title(f"Final OFV: {result.ofv:.2f} (after {result.n_iterations} iterations)")
            ax.set_ylabel("Objective Function Value")

            # Add convergence status
            status = "Converged" if result.convergence else "Not Converged"
            color = "green" if result.convergence else "red"
            ax.text(0, result.ofv / 2, status, ha='center', fontsize=14,
                    color=color, fontweight='bold')

            fig.tight_layout()
            return fig
        else:
            import plotly.graph_objects as go
            fig = go.Figure(go.Bar(x=['Final OFV'], y=[result.ofv],
                                   marker_color=get_color(0)))
            fig.update_layout(
                title=f"Final OFV: {result.ofv:.2f} ({result.n_iterations} iterations)",
                yaxis_title="Objective Function Value"
            )
            return fig

    # With convergence history
    history = result.convergence_history
    iterations = list(range(len(history['ofv'])))
    ofv_values = history['ofv']

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        if show_gradient and 'gradient_norm' in history:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
        else:
            fig, ax1 = plt.subplots(figsize=figsize)

        ax1.plot(iterations, ofv_values, color=get_color(0),
                 linewidth=theme["line_width"])
        ax1.set_ylabel("OFV")
        ax1.set_title(title or f"Convergence ({result.method.upper()})")
        ax1.grid(alpha=0.3)

        # Mark SAEM phases if applicable
        if result.method.lower() == 'saem' and hasattr(result, 'saem_n_burn'):
            n_burn = result.saem_n_burn
            ax1.axvline(n_burn, color='red', linestyle='--', label='End burn-in')
            ax1.legend()

        if show_gradient and 'gradient_norm' in history:
            ax2.semilogy(iterations, history['gradient_norm'],
                         color=get_color(1), linewidth=theme["line_width"])
            ax2.set_xlabel("Iteration")
            ax2.set_ylabel("Gradient Norm")
            ax2.grid(alpha=0.3)
        else:
            ax1.set_xlabel("Iteration")

        fig.tight_layout()
        return fig

    else:
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        if show_gradient and 'gradient_norm' in history:
            fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                                subplot_titles=["OFV", "Gradient Norm"])
            fig.add_trace(go.Scatter(x=iterations, y=ofv_values,
                                     mode='lines', name='OFV'), row=1, col=1)
            fig.add_trace(go.Scatter(x=iterations, y=history['gradient_norm'],
                                     mode='lines', name='Gradient'), row=2, col=1)
            fig.update_yaxes(type="log", row=2, col=1)
        else:
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=iterations, y=ofv_values,
                                     mode='lines', name='OFV'))

        fig.update_layout(
            title=title or f"Convergence ({result.method.upper()})",
            xaxis_title="Iteration",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Parameter Convergence (Individual Parameters)
# ============================================================================

def plot_parameter_convergence(
    result: Any,
    param_names: Optional[List[str]] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 8),
) -> Any:
    """
    Plot individual parameter traces during estimation.

    Shows how each parameter evolved during the estimation process.
    Useful for diagnosing convergence issues.

    Args:
        result: EstimationResult with parameter history
        param_names: Names of parameters to plot (None = all)
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object
    """
    backend = get_backend()
    theme = get_theme_config()

    if not hasattr(result, 'convergence_history') or result.convergence_history is None:
        raise ValueError("Result does not contain convergence history")

    history = result.convergence_history
    if 'theta_history' not in history:
        raise ValueError("No parameter history available")

    theta_history = np.array(history['theta_history'])
    n_iter, n_params = theta_history.shape

    all_names = getattr(result, 'theta_names', None) or [f"θ{i+1}" for i in range(n_params)]

    if param_names is not None:
        indices = [all_names.index(n) for n in param_names if n in all_names]
        names = [all_names[i] for i in indices]
    else:
        indices = list(range(n_params))
        names = all_names

    colors = get_colors(len(names))
    iterations = list(range(n_iter))

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        n_plots = len(names)
        n_cols = min(3, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n_plots > 1:
            axes = axes.flatten()
        else:
            axes = [axes]

        for i, (idx, name) in enumerate(zip(indices, names)):
            ax = axes[i]
            ax.plot(iterations, theta_history[:, idx], color=colors[i],
                    linewidth=theme["line_width"])
            ax.axhline(result.theta[idx], color='red', linestyle='--',
                       alpha=0.7, label='Final')
            ax.set_title(name)
            ax.set_xlabel("Iteration")
            ax.grid(alpha=0.3)

        for i in range(len(names), len(axes)):
            axes[i].set_visible(False)

        if title:
            fig.suptitle(title)

        fig.tight_layout()
        return fig

    else:
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        n_plots = len(names)
        n_cols = min(3, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols

        fig = make_subplots(rows=n_rows, cols=n_cols, subplot_titles=names)

        for i, (idx, name) in enumerate(zip(indices, names)):
            row = i // n_cols + 1
            col = i % n_cols + 1
            fig.add_trace(
                go.Scatter(x=iterations, y=theta_history[:, idx],
                           mode='lines', name=name, line=dict(color=colors[i])),
                row=row, col=col
            )
            fig.add_hline(y=result.theta[idx], line_dash="dash",
                          line_color="red", row=row, col=col)

        fig.update_layout(
            title=title or "Parameter Convergence",
            width=figsize[0] * 100,
            height=figsize[1] * 100,
            showlegend=False
        )

        return fig


# ============================================================================
# Shrinkage Plot
# ============================================================================

def plot_shrinkage(
    result: Any,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 5),
    warning_threshold: float = 30.0,
) -> Any:
    """
    Plot eta (random effect) shrinkage as bar chart.

    Shrinkage >30% indicates limited information in the data
    to estimate individual parameters.

    Args:
        result: EstimationResult with diagnostics
        title: Plot title
        figsize: Figure size
        warning_threshold: Shrinkage level above which to show warning color

    Returns:
        Figure object

    Example:
        >>> fig = plot_shrinkage(result, warning_threshold=25.0)
    """
    backend = get_backend()
    theme = get_theme_config()

    if result.diagnostics is None:
        raise ValueError("No diagnostics available in result")

    shrinkage = np.array(result.diagnostics.eta_shrinkage) * 100  # Convert to %
    omega_names = getattr(result, 'omega_names', None) or [f"η{i+1}" for i in range(len(shrinkage))]

    # Color based on threshold
    colors = [NEOPKPD_COLORS["error"] if s > warning_threshold
              else NEOPKPD_COLORS["success"] for s in shrinkage]

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        bars = ax.bar(omega_names, shrinkage, color=colors)

        # Add threshold line
        ax.axhline(warning_threshold, color='red', linestyle='--',
                   alpha=0.7, label=f'Warning ({warning_threshold}%)')

        # Add value labels
        for bar, val in zip(bars, shrinkage):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f'{val:.1f}%', ha='center', va='bottom', fontsize=9)

        ax.set_ylabel("Shrinkage (%)")
        ax.set_title(title or "Eta Shrinkage")
        ax.set_ylim(0, max(100, max(shrinkage) * 1.1))
        ax.legend()
        ax.grid(axis='y', alpha=0.3)

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        fig = go.Figure(go.Bar(
            x=omega_names,
            y=shrinkage,
            marker_color=colors,
            text=[f"{s:.1f}%" for s in shrinkage],
            textposition='outside'
        ))

        fig.add_hline(y=warning_threshold, line_dash="dash",
                      line_color="red", annotation_text=f"Warning ({warning_threshold}%)")

        fig.update_layout(
            title=title or "Eta Shrinkage",
            yaxis_title="Shrinkage (%)",
            yaxis_range=[0, max(100, max(shrinkage) * 1.1)],
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Eta Distributions
# ============================================================================

def plot_eta_distributions(
    result: Any,
    n_cols: int = 3,
    show_reference: bool = True,
    title: Optional[str] = None,
    subplot_size: Tuple[float, float] = (4, 3),
) -> Any:
    """
    Plot histograms of individual random effect (eta) estimates.

    Shows the distribution of EBEs (Empirical Bayes Estimates) for
    each random effect. Reference line shows expected N(0, omega).

    Args:
        result: EstimationResult with individuals
        n_cols: Number of columns in panel
        show_reference: Show expected normal distribution
        title: Plot title
        subplot_size: Size of each subplot

    Returns:
        Figure object

    Example:
        >>> fig = plot_eta_distributions(result)
    """
    backend = get_backend()
    theme = get_theme_config()

    # Extract etas from individuals
    etas = []
    for ind in result.individuals:
        etas.append(ind.eta)

    etas = np.array(etas)
    n_etas = etas.shape[1]

    omega_names = getattr(result, 'omega_names', None) or [f"η{i+1}" for i in range(n_etas)]
    omega_var = [result.omega[i][i] for i in range(n_etas)]

    n_rows = (n_etas + n_cols - 1) // n_cols
    figsize = (subplot_size[0] * n_cols, subplot_size[1] * n_rows)

    colors = get_colors(n_etas)

    if backend == "matplotlib":
        import matplotlib.pyplot as plt
        from scipy import stats

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n_etas > 1:
            axes = axes.flatten()
        else:
            axes = [axes]

        for i in range(n_etas):
            ax = axes[i]
            eta_i = etas[:, i]

            ax.hist(eta_i, bins=20, density=True, alpha=0.7,
                    color=colors[i], edgecolor='white')

            if show_reference:
                x = np.linspace(eta_i.min() - 0.5, eta_i.max() + 0.5, 100)
                y = stats.norm.pdf(x, 0, np.sqrt(omega_var[i]))
                ax.plot(x, y, color='red', linewidth=2, linestyle='--',
                        label='Expected')

            ax.axvline(0, color='black', linestyle='-', alpha=0.5)
            ax.set_title(omega_names[i])
            ax.set_xlabel("EBE")
            ax.set_ylabel("Density")

            if i == 0 and show_reference:
                ax.legend()

        for i in range(n_etas, len(axes)):
            axes[i].set_visible(False)

        if title:
            fig.suptitle(title)

        fig.tight_layout()
        return fig

    else:
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(rows=n_rows, cols=n_cols, subplot_titles=omega_names)

        for i in range(n_etas):
            row = i // n_cols + 1
            col = i % n_cols + 1

            eta_i = etas[:, i]

            fig.add_trace(
                go.Histogram(x=eta_i, nbinsx=20, histnorm='probability density',
                             marker_color=colors[i], name=omega_names[i],
                             showlegend=False),
                row=row, col=col
            )

            fig.add_vline(x=0, line_dash="solid", line_color="black",
                          opacity=0.5, row=row, col=col)

        fig.update_layout(
            title=title or "Eta Distributions",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Individual Parameters (EBE)
# ============================================================================

def plot_individual_parameters(
    result: Any,
    param_names: Optional[List[str]] = None,
    sort_by: Optional[str] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 6),
) -> Any:
    """
    Plot individual parameter estimates (EBEs) vs population values.

    Shows the spread of individual parameter estimates around
    the population typical values.

    Args:
        result: EstimationResult with individuals
        param_names: Which parameters to plot (default: all)
        sort_by: Sort individuals by this parameter
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object
    """
    backend = get_backend()
    theme = get_theme_config()

    # Get individual etas and convert to parameters
    theta = np.array(result.theta)
    n_theta = len(theta)

    # Extract etas
    etas = []
    subject_ids = []
    for ind in result.individuals:
        etas.append(ind.eta)
        subject_ids.append(ind.subject_id)

    etas = np.array(etas)
    n_subjects = len(subject_ids)
    n_etas = etas.shape[1] if len(etas.shape) > 1 else 1

    # Calculate individual parameters (assuming log-normal)
    # Individual = theta * exp(eta)
    # Note: etas may have fewer columns than theta (not all params have IIV)
    n_params_with_iiv = min(n_theta, n_etas)
    individual_params = theta[:n_params_with_iiv] * np.exp(etas[:, :n_params_with_iiv])

    theta_names = getattr(result, 'theta_names', None) or [f"θ{i+1}" for i in range(n_theta)]
    # Limit to params with IIV
    theta_names = theta_names[:n_params_with_iiv]
    n_theta = n_params_with_iiv

    if param_names is not None:
        indices = [theta_names.index(n) for n in param_names if n in theta_names]
        names = [theta_names[i] for i in indices]
    else:
        indices = list(range(n_theta))
        names = theta_names

    colors = get_colors(len(names))

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, len(names), figsize=figsize)
        if len(names) == 1:
            axes = [axes]

        for i, (idx, name) in enumerate(zip(indices, names)):
            ax = axes[i]

            # Sort if requested
            ind_vals = individual_params[:, idx]
            order = np.argsort(ind_vals) if sort_by == name else np.arange(n_subjects)

            x = np.arange(n_subjects)
            ax.scatter(x, ind_vals[order], color=colors[i], s=20, alpha=0.6)
            ax.axhline(theta[idx], color='red', linewidth=2, linestyle='--',
                       label='Population')

            ax.set_xlabel("Subject")
            ax.set_ylabel(name)
            ax.set_title(name)
            ax.grid(alpha=0.3)

            if i == 0:
                ax.legend()

        if title:
            fig.suptitle(title)

        fig.tight_layout()
        return fig

    else:
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(rows=1, cols=len(names), subplot_titles=names)

        for i, (idx, name) in enumerate(zip(indices, names)):
            ind_vals = individual_params[:, idx]
            x = list(range(len(ind_vals)))

            fig.add_trace(
                go.Scatter(x=x, y=ind_vals, mode='markers',
                           marker=dict(size=6, color=colors[i]),
                           name=name, showlegend=False),
                row=1, col=i + 1
            )

            fig.add_hline(y=theta[idx], line_dash="dash",
                          line_color="red", row=1, col=i + 1)

        fig.update_layout(
            title=title or "Individual Parameters",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# OFV Model Comparison
# ============================================================================

def plot_ofv_comparison(
    comparison: Any,
    show_delta: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 5),
) -> Any:
    """
    Plot model comparison metrics (OFV, AIC, BIC).

    Args:
        comparison: ModelComparisonResult from compare_models()
        show_delta: Show delta values from best model
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> comparison = compare_models([result1, result2], ["Base", "Full"])
        >>> fig = plot_ofv_comparison(comparison)
    """
    backend = get_backend()
    theme = get_theme_config()

    model_names = comparison.model_names
    ofv = comparison.ofv_values
    aic = comparison.aic_values
    bic = comparison.bic_values

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=figsize)

        metrics = [("OFV", ofv), ("AIC", aic), ("BIC", bic)]
        colors_list = [get_color(0), get_color(1), get_color(2)]

        for ax, (metric_name, values), color in zip(axes, metrics, colors_list):
            bars = ax.bar(model_names, values, color=color, alpha=0.8)
            ax.set_ylabel(metric_name)
            ax.set_title(metric_name)
            ax.grid(axis='y', alpha=0.3)

            # Rotate labels if needed
            if len(model_names) > 3:
                ax.set_xticklabels(model_names, rotation=45, ha='right')

            # Add value labels
            for bar, val in zip(bars, values):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                        f'{val:.1f}', ha='center', va='bottom', fontsize=8)

        if title:
            fig.suptitle(title)

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        fig = make_subplots(rows=1, cols=3, subplot_titles=["OFV", "AIC", "BIC"])

        metrics = [ofv, aic, bic]
        colors_list = [get_color(0), get_color(1), get_color(2)]

        for i, (values, color) in enumerate(zip(metrics, colors_list)):
            fig.add_trace(
                go.Bar(x=model_names, y=values, marker_color=color,
                       text=[f"{v:.1f}" for v in values], textposition='outside'),
                row=1, col=i + 1
            )

        fig.update_layout(
            title=title or "Model Comparison",
            width=figsize[0] * 100,
            height=figsize[1] * 100,
            showlegend=False
        )

        return fig


# ============================================================================
# Correlation Matrix
# ============================================================================

def plot_correlation_matrix(
    result: Any,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 8),
    cmap: str = "RdBu_r",
) -> Any:
    """
    Plot parameter correlation matrix from estimation.

    Shows correlations between fixed effect estimates. High correlations
    (>0.9) may indicate overparameterization.

    Args:
        result: EstimationResult
        title: Plot title
        figsize: Figure size
        cmap: Colormap name

    Returns:
        Figure object
    """
    backend = get_backend()

    # Check if correlation info available
    if not hasattr(result, 'theta_correlation') or result.theta_correlation is None:
        # Compute from covariance if available
        if hasattr(result, 'theta_covariance') and result.theta_covariance is not None:
            cov = np.array(result.theta_covariance)
            std = np.sqrt(np.diag(cov))
            corr = cov / np.outer(std, std)
        else:
            raise ValueError("No correlation or covariance information available")
    else:
        corr = np.array(result.theta_correlation)

    theta_names = getattr(result, 'theta_names', None) or [f"θ{i+1}" for i in range(len(corr))]

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        im = ax.imshow(corr, cmap=cmap, vmin=-1, vmax=1)
        plt.colorbar(im, ax=ax)

        # Annotations
        for i in range(len(corr)):
            for j in range(len(corr)):
                text = f"{corr[i, j]:.2f}"
                color = "white" if abs(corr[i, j]) > 0.5 else "black"
                ax.text(j, i, text, ha="center", va="center", color=color, fontsize=9)

        ax.set_xticks(np.arange(len(theta_names)))
        ax.set_yticks(np.arange(len(theta_names)))
        ax.set_xticklabels(theta_names, rotation=45, ha="right")
        ax.set_yticklabels(theta_names)
        ax.set_title(title or "Parameter Correlation Matrix")

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        fig = go.Figure(data=go.Heatmap(
            z=corr,
            x=theta_names,
            y=theta_names,
            colorscale="RdBu_r",
            zmin=-1, zmax=1,
            text=[[f"{v:.2f}" for v in row] for row in corr],
            texttemplate="%{text}"
        ))

        fig.update_layout(
            title=title or "Parameter Correlation Matrix",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig


# ============================================================================
# Sigma Residuals
# ============================================================================

def plot_sigma_residuals(
    result: Any,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 5),
) -> Any:
    """
    Plot sigma (residual error) estimates visualization.

    Args:
        result: EstimationResult
        title: Plot title
        figsize: Figure size

    Returns:
        Figure object
    """
    backend = get_backend()
    theme = get_theme_config()

    sigma = result.sigma
    sigma_se = result.sigma_se

    # Handle different sigma structures
    if isinstance(sigma, dict):
        sigma_names = list(sigma.keys())
        sigma_vals = list(sigma.values())
        sigma_se_vals = [sigma_se.get(k, None) for k in sigma_names] if sigma_se else [None] * len(sigma_names)
    else:
        sigma_names = ["σ"]
        sigma_vals = [sigma]
        sigma_se_vals = [sigma_se] if sigma_se else [None]

    if backend == "matplotlib":
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        x = np.arange(len(sigma_names))
        bars = ax.bar(x, sigma_vals, color=get_color(0), alpha=0.8)

        # Error bars for SE
        if any(se is not None for se in sigma_se_vals):
            for i, (val, se) in enumerate(zip(sigma_vals, sigma_se_vals)):
                if se is not None:
                    ax.errorbar(i, val, yerr=1.96 * se, color='black',
                                capsize=5, capthick=2)

        ax.set_xticks(x)
        ax.set_xticklabels(sigma_names)
        ax.set_ylabel("Estimate")
        ax.set_title(title or "Residual Error Estimates")
        ax.grid(axis='y', alpha=0.3)

        # Add CV annotation for proportional error
        for i, (bar, val, name) in enumerate(zip(bars, sigma_vals, sigma_names)):
            if 'prop' in name.lower():
                cv = np.sqrt(val) * 100 if val > 0 else 0
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                        f'CV={cv:.1f}%', ha='center', va='bottom', fontsize=9)

        fig.tight_layout()
        return fig

    else:
        import plotly.graph_objects as go

        error_y = None
        if any(se is not None for se in sigma_se_vals):
            error_y = dict(type='data', array=[1.96 * se if se else 0 for se in sigma_se_vals])

        fig = go.Figure(go.Bar(
            x=sigma_names,
            y=sigma_vals,
            marker_color=get_color(0),
            error_y=error_y
        ))

        fig.update_layout(
            title=title or "Residual Error Estimates",
            yaxis_title="Estimate",
            width=figsize[0] * 100,
            height=figsize[1] * 100
        )

        return fig
