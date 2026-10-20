"""
NeoPKPD - Professional PK/PD Simulation Platform

A Python interface to the NeoPKPD Julia core for pharmacokinetic and
pharmacodynamic simulations.

Features:
- Single PK simulations:
  - One-compartment: IV bolus, IV infusion, oral first-order
  - Two-compartment: IV bolus, oral first-order
  - Three-compartment: IV bolus
  - Advanced absorption: Transit compartment model
  - Nonlinear elimination: Michaelis-Menten kinetics
- Coupled PK-PD simulations:
  - Direct Emax model
  - Sigmoid Emax model (Hill equation)
  - Biophase equilibration (effect compartment)
  - Indirect response turnover model
- Population simulations with IIV, IOV, and covariates
- Parameter estimation (FOCE-I, SAEM, Laplacian)
- Visual Predictive Check (VPC, pcVPC)
- Non-Compartmental Analysis (NCA)
- CDISC/SDTM data import
- NONMEM/Monolix model import
- Residual error models
- Sensitivity analysis
- Artifact replay for reproducibility
- PK/PD metrics (Cmax, AUC, Tmax, Emin, time below threshold)

Quick Start:
    >>> import neopkpd
    >>> neopkpd.init_julia()  # Initialize once per session

    >>> # Run a simple PK simulation
    >>> result = neopkpd.simulate_pk_iv_bolus(
    ...     cl=1.0, v=10.0,
    ...     doses=[{"time": 0.0, "amount": 100.0}],
    ...     t0=0.0, t1=24.0,
    ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    ... )

    >>> # Run IV infusion
    >>> result = neopkpd.simulate_pk_iv_bolus(
    ...     cl=1.0, v=10.0,
    ...     doses=[{"time": 0.0, "amount": 100.0, "duration": 1.0}],  # 1-hour infusion
    ...     t0=0.0, t1=24.0,
    ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    ... )

    >>> # Run parameter estimation
    >>> from neopkpd.estimation import estimate, EstimationConfig
    >>> config = EstimationConfig(method="foce", theta_init=[1.0, 10.0])
    >>> result = estimate(observed_data, "OneCompIVBolus", config, grid)

    >>> # Run VPC
    >>> from neopkpd.vpc import compute_vpc, VPCConfig
    >>> vpc_result = compute_vpc(observed, pop_spec, grid)

For more information, see the docstrings for individual functions.
"""

# Core initialization and utilities
from ._core import (
    init_julia,
    version,
    SensitivityMetrics,
    SensitivityResult,
    PopulationSensitivityResult,
    # Global Sensitivity Analysis types
    SobolIndex,
    SobolResult,
    MorrisIndex,
    MorrisResult,
)

# Artifact operations
from .artifacts import (
    replay_artifact,
    write_single_artifact,
)

# PK/PD metrics
from .metrics import (
    cmax,
    auc_trapezoid,
    emin,
    time_below,
    auc_above_baseline,
    tmax,
    half_life,
)

# Simulations - One-compartment PK
from .simulations.pk_onecomp import (
    simulate_pk_iv_bolus,
    simulate_pk_oral_first_order,
)

# Simulations - Two-compartment PK
from .simulations.pk_twocomp import (
    simulate_pk_twocomp_iv_bolus,
    simulate_pk_twocomp_oral,
)

# Simulations - Three-compartment PK
from .simulations.pk_threecomp import (
    simulate_pk_threecomp_iv_bolus,
)

# Simulations - Advanced PK
from .simulations.pk_advanced import (
    simulate_pk_transit_absorption,
    simulate_pk_michaelis_menten,
)

# Simulations - Custom PK models
from .simulations.pk_custom import (
    simulate_pk_tmdd_custom,
    simulate_pk_parallel_absorption,
    simulate_pk_enterohepatic_recirculation,
    simulate_pk_autoinduction,
)

# Simulations - Coupled PKPD
from .simulations.pkpd import (
    # Basic effect models
    simulate_pkpd_direct_emax,
    simulate_pkpd_sigmoid_emax,
    simulate_pkpd_biophase_equilibration,
    # Indirect response models (IRM)
    simulate_pkpd_indirect_response,  # IRM-III (inhibition of Kout)
    simulate_pkpd_irm1,               # IRM-I (inhibition of Kin)
    simulate_pkpd_irm2,               # IRM-II (stimulation of Kin)
    simulate_pkpd_irm4,               # IRM-IV (stimulation of Kout)
    # Advanced PD models
    simulate_pkpd_transit_compartment,
    simulate_pkpd_disease_progression,
    # Tolerance models
    simulate_pkpd_tolerance_counter_regulation,
    simulate_pkpd_receptor_regulation,
    # Combination effect models
    simulate_pkpd_bliss_independence,
    simulate_pkpd_competitive_inhibition,
    simulate_pkpd_drug_interaction,
)

# Simulations - Population
from .simulations.population import (
    simulate_population_iv_bolus,
    simulate_population_oral,
)

# Simulations - Sensitivity (Local)
from .simulations.sensitivity import (
    run_sensitivity,
)

# Simulations - Global Sensitivity Analysis
from .simulations.gsa import (
    run_sobol_sensitivity,
    run_morris_sensitivity,
)

# NCA - Non-Compartmental Analysis
from .nca import (
    # Full workflow
    run_nca,
    run_population_nca,
    summarize_population_nca,
    NCAConfig,
    NCAResult,
    # Primary metrics
    nca_cmax,
    nca_tmax,
    nca_cmin,
    nca_clast,
    nca_tlast,
    nca_cavg,
    # AUC calculations
    auc_0_t,
    auc_0_inf,
    auc_0_tau,
    auc_partial,
    aumc_0_t,
    aumc_0_inf,
    # Lambda-z and half-life
    estimate_lambda_z,
    nca_half_life,
    # PK parameters
    nca_mrt,
    nca_cl_f,
    nca_vz_f,
    nca_vss,
    nca_cl,
    nca_vz,
    # Multiple dose metrics
    nca_accumulation_index,
    nca_ptf,
    nca_swing,
    nca_linearity_index,
    nca_time_to_steady_state,
    # Bioequivalence
    bioequivalence_90ci,
    tost_analysis,
    be_conclusion,
    geometric_mean_ratio,
    geometric_mean,
    within_subject_cv,
    # Study designs and RSABE
    BEStudyDesign,
    ReplicateDesign,
    RegulatoryGuidance,
    RSABEConfig,
)

# Parameter Estimation (NLME)
from .estimation import (
    estimate,
    run_bootstrap,
    compare_models,
    likelihood_ratio_test,
    compute_diagnostics,
    get_individual_predictions,
    EstimationConfig,
    BootstrapConfig,
    BLQConfig,
    IOVSpec,
    CovariateOnIIV,
    EstimationResult,
    BootstrapResult,
    IndividualEstimate,
    DiagnosticsSummary,
    BLQSummary,
    ModelComparisonResult,
    BLQMethod,
    OmegaStructure,
    BootstrapCIMethod,
)

# Visual Predictive Check
from .vpc import (
    # VPC computation
    compute_vpc,
    compute_pcvpc,
    compute_vpc_python,
    compute_stratified_vpc,
    compute_vpc_with_blq,
    # Configuration and types
    VPCConfig,
    VPCResult,
    VPCBin,
    VPCPercentileData,
    StratifiedVPCResult,
    BLQBinStats,
    # Binning strategies
    BinningStrategy,
    BinDefinition,
    QuantileBinning,
    EqualWidthBinning,
    KMeansBinning,
    # BLQ handling
    BLQMethod as VPCBLQMethod,
    handle_blq,
    # Result extraction helpers
    get_bin_midpoints,
    get_observed_percentile,
    get_simulated_median,
    get_simulated_ci,
    get_blq_observed,
    get_blq_simulated,
)

# CDISC Data Import
from .data import (
    read_cdisc_pc,
    read_cdisc_ex,
    read_cdisc_dm,
    read_cdisc_dataset,
    cdisc_to_observed_data,
    validate_cdisc_dataset,
    SubjectData,
    ObservedData,
)

# Model Import (NONMEM/Monolix)
from .import_ import (
    import_nonmem,
    import_monolix,
    import_model,
    parse_nonmem_control,
    parse_monolix_project,
    ImportedModel,
    NONMEMControlFile,
    MonolixProject,
)


__all__ = [
    # Core
    "init_julia",
    "version",

    # Data classes - Local Sensitivity
    "SensitivityMetrics",
    "SensitivityResult",
    "PopulationSensitivityResult",
    # Data classes - Global Sensitivity Analysis
    "SobolIndex",
    "SobolResult",
    "MorrisIndex",
    "MorrisResult",

    # Artifacts
    "replay_artifact",
    "write_single_artifact",

    # Single PK - One-compartment
    "simulate_pk_iv_bolus",
    "simulate_pk_oral_first_order",

    # Single PK - Two-compartment
    "simulate_pk_twocomp_iv_bolus",
    "simulate_pk_twocomp_oral",

    # Single PK - Three-compartment
    "simulate_pk_threecomp_iv_bolus",

    # Single PK - Advanced models
    "simulate_pk_transit_absorption",
    "simulate_pk_michaelis_menten",

    # Coupled PKPD - Direct PD models
    "simulate_pkpd_direct_emax",
    "simulate_pkpd_sigmoid_emax",
    "simulate_pkpd_biophase_equilibration",

    # Coupled PKPD - ODE-based PD models
    "simulate_pkpd_indirect_response",

    # Population
    "simulate_population_iv_bolus",
    "simulate_population_oral",

    # Sensitivity - Local
    "run_sensitivity",
    # Sensitivity - Global (Sobol', Morris)
    "run_sobol_sensitivity",
    "run_morris_sensitivity",

    # Metrics
    "cmax",
    "auc_trapezoid",
    "emin",
    "time_below",
    "auc_above_baseline",
    "tmax",
    "half_life",

    # NCA - Full workflow
    "run_nca",
    "run_population_nca",
    "summarize_population_nca",
    "NCAConfig",
    "NCAResult",
    # NCA - Primary metrics
    "nca_cmax",
    "nca_tmax",
    "nca_cmin",
    "nca_clast",
    "nca_tlast",
    "nca_cavg",
    # NCA - AUC
    "auc_0_t",
    "auc_0_inf",
    "auc_0_tau",
    "auc_partial",
    "aumc_0_t",
    "aumc_0_inf",
    # NCA - Lambda-z
    "estimate_lambda_z",
    "nca_half_life",
    # NCA - PK parameters
    "nca_mrt",
    "nca_cl_f",
    "nca_vz_f",
    "nca_vss",
    "nca_cl",
    "nca_vz",
    # NCA - Multiple dose
    "nca_accumulation_index",
    "nca_ptf",
    "nca_swing",
    "nca_linearity_index",
    "nca_time_to_steady_state",
    # NCA - Bioequivalence
    "bioequivalence_90ci",
    "tost_analysis",
    "be_conclusion",
    "geometric_mean_ratio",
    "geometric_mean",
    "within_subject_cv",
    # NCA - Study designs and RSABE
    "BEStudyDesign",
    "ReplicateDesign",
    "RegulatoryGuidance",
    "RSABEConfig",

    # Estimation (NLME)
    "estimate",
    "run_bootstrap",
    "compare_models",
    "likelihood_ratio_test",
    "compute_diagnostics",
    "get_individual_predictions",
    "EstimationConfig",
    "BootstrapConfig",
    "BLQConfig",
    "IOVSpec",
    "CovariateOnIIV",
    "EstimationResult",
    "BootstrapResult",
    "IndividualEstimate",
    "DiagnosticsSummary",
    "BLQSummary",
    "ModelComparisonResult",
    "BLQMethod",
    "OmegaStructure",
    "BootstrapCIMethod",

    # VPC - Computation
    "compute_vpc",
    "compute_pcvpc",
    "compute_vpc_python",
    "compute_stratified_vpc",
    "compute_vpc_with_blq",
    # VPC - Types
    "VPCConfig",
    "VPCResult",
    "VPCBin",
    "VPCPercentileData",
    "StratifiedVPCResult",
    "BLQBinStats",
    # VPC - Binning
    "BinningStrategy",
    "BinDefinition",
    "QuantileBinning",
    "EqualWidthBinning",
    "KMeansBinning",
    # VPC - BLQ
    "VPCBLQMethod",
    "handle_blq",
    # VPC - Result extraction
    "get_bin_midpoints",
    "get_observed_percentile",
    "get_simulated_median",
    "get_simulated_ci",
    "get_blq_observed",
    "get_blq_simulated",

    # CDISC Data
    "read_cdisc_pc",
    "read_cdisc_ex",
    "read_cdisc_dm",
    "read_cdisc_dataset",
    "cdisc_to_observed_data",
    "validate_cdisc_dataset",
    "SubjectData",
    "ObservedData",

    # Model Import
    "import_nonmem",
    "import_monolix",
    "import_model",
    "parse_nonmem_control",
    "parse_monolix_project",
    "ImportedModel",
    "NONMEMControlFile",
    "MonolixProject",

    # Compliance (FDA 21 CFR Part 11)
    "compliance",
]


# FDA 21 CFR Part 11 Compliance Module
from . import compliance
