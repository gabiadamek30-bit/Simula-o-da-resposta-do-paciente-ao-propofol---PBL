"""
NeoPKPD NCA Module - Non-Compartmental Analysis

FDA/EMA compliant Non-Compartmental Analysis for pharmacokinetic data.

Features:
- Primary exposure metrics (Cmax, Tmax, Cmin, AUC)
- Lambda-z terminal slope estimation
- AUC calculations (0-t, 0-inf, 0-tau, partial)
- PK parameters (t1/2, MRT, CL/F, Vz/F, Vss)
- Multiple dose metrics (accumulation, PTF, swing)
- Dose linearity assessment
- Bioequivalence analysis (90% CI, TOST)
- Reference-scaled BE for highly variable drugs (RSABE/ABEL)
- Study design support (2x2 crossover, replicate designs)

Example:
    >>> import neopkpd
    >>> from neopkpd.nca import run_nca, NCAConfig
    >>> neopkpd.init_julia()
    >>>
    >>> # Run NCA on concentration-time data
    >>> t = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    >>> c = [0.0, 1.2, 2.0, 1.8, 1.2, 0.6, 0.3, 0.075]
    >>> result = run_nca(t, c, dose=100.0)
    >>> print(f"Cmax: {result.cmax}")
    >>> print(f"AUC0-inf: {result.auc_0_inf}")
    >>> print(f"t1/2: {result.t_half}")
"""

from .metrics import (
    # Primary exposure
    nca_cmax,
    nca_tmax,
    nca_cmin,
    nca_clast,
    nca_tlast,
    nca_cavg,
    # Additional exposure
    nca_ctrough,
    nca_c_at_time,
    time_above_concentration,
    # Lambda-z
    estimate_lambda_z,
    nca_half_life,
    # AUC
    auc_0_t,
    auc_0_inf,
    auc_0_tau,
    auc_partial,
    aumc_0_t,
    aumc_0_inf,
    # PK parameters
    nca_mrt,
    nca_cl_f,
    nca_vz_f,
    nca_vss,
    nca_bioavailability,
    nca_cl,
    nca_vz,
    nca_mrt_iv,
    nca_cl_ss,
    nca_vss_from_aumc,
    nca_vc,
    nca_mean_absorption_time,
    nca_c0_backextrap,
    nca_c0_from_regression,
    # Multiple dose
    nca_accumulation_index,
    nca_ptf,
    nca_swing,
    nca_linearity_index,
    nca_time_to_steady_state,
    nca_accumulation_predicted,
    nca_accumulation_cmax,
    nca_accumulation_cmin,
    nca_dose_normalized_auc,
    nca_dose_normalized_cmax,
    nca_time_to_steady_state_doses,
    nca_effective_half_life,
)

from .analysis import (
    run_nca,
    run_population_nca,
    summarize_population_nca,
    NCAConfig,
    NCAResult,
)

from .bioequivalence import (
    # Basic BE
    bioequivalence_90ci,
    tost_analysis,
    be_conclusion,
    geometric_mean_ratio,
    geometric_mean,
    within_subject_cv,
    run_bioequivalence,
    BioequivalenceResult,
    # Study designs
    BEStudyDesign,
    ReplicateDesign,
    RegulatoryGuidance,
    # Design-aware CI
    bioequivalence_90ci_design,
    compute_be_degrees_of_freedom,
    # Reference-scaled BE
    RSABEConfig,
    RSABEResult,
    ABELResult,
    compute_swr,
    compute_within_subject_variance_reference,
    rsabe_analysis,
    abel_analysis,
    reference_scaled_be,
    abel_scaled_limits,
)

__all__ = [
    # Config and Results
    "NCAConfig",
    "NCAResult",
    # Full workflow
    "run_nca",
    "run_population_nca",
    "summarize_population_nca",
    # Primary exposure
    "nca_cmax",
    "nca_tmax",
    "nca_cmin",
    "nca_clast",
    "nca_tlast",
    "nca_cavg",
    # Additional exposure
    "nca_ctrough",
    "nca_c_at_time",
    "time_above_concentration",
    # Lambda-z
    "estimate_lambda_z",
    "nca_half_life",
    # AUC
    "auc_0_t",
    "auc_0_inf",
    "auc_0_tau",
    "auc_partial",
    "aumc_0_t",
    "aumc_0_inf",
    # PK parameters
    "nca_mrt",
    "nca_cl_f",
    "nca_vz_f",
    "nca_vss",
    "nca_bioavailability",
    "nca_cl",
    "nca_vz",
    "nca_mrt_iv",
    "nca_cl_ss",
    "nca_vss_from_aumc",
    "nca_vc",
    "nca_mean_absorption_time",
    "nca_c0_backextrap",
    "nca_c0_from_regression",
    # Multiple dose
    "nca_accumulation_index",
    "nca_ptf",
    "nca_swing",
    "nca_linearity_index",
    "nca_time_to_steady_state",
    "nca_accumulation_predicted",
    "nca_accumulation_cmax",
    "nca_accumulation_cmin",
    "nca_dose_normalized_auc",
    "nca_dose_normalized_cmax",
    "nca_time_to_steady_state_doses",
    "nca_effective_half_life",
    # Basic Bioequivalence
    "bioequivalence_90ci",
    "tost_analysis",
    "be_conclusion",
    "geometric_mean_ratio",
    "geometric_mean",
    "within_subject_cv",
    "run_bioequivalence",
    "BioequivalenceResult",
    # Study designs
    "BEStudyDesign",
    "ReplicateDesign",
    "RegulatoryGuidance",
    # Design-aware CI
    "bioequivalence_90ci_design",
    "compute_be_degrees_of_freedom",
    # Reference-scaled BE
    "RSABEConfig",
    "RSABEResult",
    "ABELResult",
    "compute_swr",
    "compute_within_subject_variance_reference",
    "rsabe_analysis",
    "abel_analysis",
    "reference_scaled_be",
    "abel_scaled_limits",
]
