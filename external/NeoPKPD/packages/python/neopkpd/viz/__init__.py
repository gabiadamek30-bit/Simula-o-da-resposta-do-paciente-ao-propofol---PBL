"""
NeoPKPD Visualization Module

Professional visualization tools for PK/PD data, supporting both
static (matplotlib) and interactive (plotly) backends.

Features:
- PK plots: concentration-time, spaghetti, mean+ribbon
- VPC plots: standard VPC, pcVPC, stratified, BLQ
- NCA plots: lambda_z fit, AUC visualization, dose proportionality, BE
- PKPD plots: effect-concentration, hysteresis, dose-response
- Population plots: forest plots, parameter distributions, GOF
- Estimation diagnostics: convergence, shrinkage, correlations
- Bootstrap plots: distributions, CIs, stability
- Sensitivity plots: tornado, spider, heatmap, waterfall
- Trial plots: power curves, tornado plots, Kaplan-Meier

Example:
    >>> import neopkpd
    >>> from neopkpd import viz
    >>>
    >>> # Set backend (matplotlib or plotly)
    >>> viz.set_backend("matplotlib")
    >>>
    >>> # Plot simulation result
    >>> result = neopkpd.simulate_pk_iv_bolus(...)
    >>> fig = viz.plot_conc_time(result)
    >>> fig.show()
    >>>
    >>> # VPC plot
    >>> vpc_result = neopkpd.vpc.compute_vpc(...)
    >>> fig = viz.plot_vpc_detailed(vpc_result)
"""

from .backends import (
    get_backend,
    set_backend,
    available_backends,
    PlotBackend,
)

from .themes import (
    get_theme,
    set_theme,
    available_themes,
    NEOPKPD_COLORS,
)

from .pk import (
    plot_conc_time,
    plot_multi_conc_time,
    plot_spaghetti,
    plot_mean_ribbon,
    plot_individual_fits,
)

from .nca import (
    plot_lambda_z_fit,
    plot_auc_visualization,
    plot_dose_proportionality,
)

from .pkpd import (
    plot_effect_conc,
    plot_hysteresis,
    plot_dose_response,
)

from .population import (
    plot_vpc,
    plot_parameter_distributions,
    plot_forest,
    plot_boxplot,
    plot_goodness_of_fit,
    plot_estimation_summary,
    plot_sensitivity,
    plot_sensitivity_tornado,
)

from .trial import (
    plot_power_curve,
    plot_tornado as plot_trial_tornado,
    plot_kaplan_meier,
    plot_endpoint_distribution,
)

# New VPC visualization module
from .vpc import (
    plot_vpc as plot_vpc_detailed,
    plot_pcvpc,
    plot_stratified_vpc,
    plot_vpc_with_blq,
    plot_vpc_ci,
)

# New estimation diagnostics module
from .estimation import (
    plot_parameter_estimates,
    plot_omega_matrix,
    plot_convergence,
    plot_parameter_convergence,
    plot_shrinkage,
    plot_eta_distributions,
    plot_individual_parameters,
    plot_ofv_comparison,
    plot_correlation_matrix,
    plot_sigma_residuals,
)

# New bootstrap visualization module
from .bootstrap import (
    plot_bootstrap_distributions,
    plot_bootstrap_ci,
    plot_bootstrap_stability,
    plot_bootstrap_correlation,
)

# New sensitivity visualization module
from .sensitivity import (
    plot_tornado,
    plot_spider,
    plot_sensitivity_heatmap,
    plot_waterfall,
    plot_one_at_a_time,
)


__all__ = [
    # Backends
    "get_backend",
    "set_backend",
    "available_backends",
    "PlotBackend",
    # Themes
    "get_theme",
    "set_theme",
    "available_themes",
    "NEOPKPD_COLORS",
    # PK plots
    "plot_conc_time",
    "plot_multi_conc_time",
    "plot_spaghetti",
    "plot_mean_ribbon",
    "plot_individual_fits",
    # NCA plots
    "plot_lambda_z_fit",
    "plot_auc_visualization",
    "plot_dose_proportionality",
    # PKPD plots
    "plot_effect_conc",
    "plot_hysteresis",
    "plot_dose_response",
    # Population plots (basic)
    "plot_vpc",
    "plot_parameter_distributions",
    "plot_forest",
    "plot_boxplot",
    # Estimation diagnostics (basic)
    "plot_goodness_of_fit",
    "plot_estimation_summary",
    # Sensitivity plots (basic)
    "plot_sensitivity",
    "plot_sensitivity_tornado",
    # Trial plots
    "plot_power_curve",
    "plot_trial_tornado",
    "plot_kaplan_meier",
    "plot_endpoint_distribution",

    # ==========================================
    # NEW VISUALIZATION FUNCTIONS
    # ==========================================

    # VPC visualization (detailed)
    "plot_vpc_detailed",
    "plot_pcvpc",
    "plot_stratified_vpc",
    "plot_vpc_with_blq",
    "plot_vpc_ci",

    # Estimation diagnostics (comprehensive)
    "plot_parameter_estimates",
    "plot_omega_matrix",
    "plot_convergence",
    "plot_parameter_convergence",
    "plot_shrinkage",
    "plot_eta_distributions",
    "plot_individual_parameters",
    "plot_ofv_comparison",
    "plot_correlation_matrix",
    "plot_sigma_residuals",

    # Bootstrap visualization
    "plot_bootstrap_distributions",
    "plot_bootstrap_ci",
    "plot_bootstrap_stability",
    "plot_bootstrap_correlation",

    # Sensitivity analysis visualization
    "plot_tornado",
    "plot_spider",
    "plot_sensitivity_heatmap",
    "plot_waterfall",
    "plot_one_at_a_time",
]
