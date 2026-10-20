"""
NeoPKPD Parameter Estimation Module

This module provides Python bindings for NLME (Nonlinear Mixed Effects)
parameter estimation methods including FOCE-I, SAEM, and Laplacian.

Features:
- Multiple estimation methods: FOCE-I, SAEM, Laplacian
- BLQ (Below Limit of Quantification) handling: M1, M2, M3
- Inter-Occasion Variability (IOV) support
- Covariate effects on IIV
- Bootstrap confidence intervals
- Diagnostics: CWRES, IWRES, shrinkage
- Model comparison: LRT, AIC, BIC
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union
from enum import Enum
import numpy as np

from .._core import _require_julia


# ============================================================================
# Enumerations
# ============================================================================

class BLQMethod(Enum):
    """Method for handling Below Limit of Quantification observations."""
    M1_DISCARD = "m1_discard"      # Discard all BLQ observations
    M2_IMPUTE_HALF = "m2_half"     # Replace BLQ with LLOQ/2
    M2_IMPUTE_ZERO = "m2_zero"     # Replace BLQ with 0
    M3_LIKELIHOOD = "m3_likelihood"  # Censored likelihood (gold standard)


class OmegaStructure(Enum):
    """Structure of omega (IIV variance) matrix."""
    DIAGONAL = "diagonal"  # No correlations between random effects
    BLOCK = "block"        # Block diagonal (correlations within blocks)
    FULL = "full"          # Full matrix (all correlations)


class BootstrapCIMethod(Enum):
    """Method for computing bootstrap confidence intervals."""
    PERCENTILE = "percentile"  # Standard percentile method
    BCA = "bca"                # Bias-corrected and accelerated
    BASIC = "basic"            # Basic (pivotal) method


# ============================================================================
# Configuration Types
# ============================================================================

@dataclass
class BLQConfig:
    """Configuration for BLQ/censoring handling in estimation.

    Attributes:
        method: BLQ handling method (M1, M2, or M3)
        lloq: Lower limit of quantification
        max_consecutive_blq: Maximum consecutive BLQ before warning
    """
    method: BLQMethod = BLQMethod.M3_LIKELIHOOD
    lloq: float = 0.0
    max_consecutive_blq: int = 5


@dataclass
class IOVSpec:
    """Specification for inter-occasion variability.

    Attributes:
        eta_name: Name of the random effect with IOV
        occasion_names: Names/labels for each occasion
        omega_iov: Variance of IOV (separate from IIV variance)
    """
    eta_name: str
    occasion_names: List[str]
    omega_iov: float = 0.04  # Default 20% CV


@dataclass
class CovariateOnIIV:
    """Specification for covariate effects on inter-individual variability.

    Allows variance parameters to depend on covariates:
        Var(η_i) = ω² * exp(θ_cov * (COV_i - COV_ref))

    Attributes:
        eta_name: Name of the random effect (e.g., "CL", "V")
        covariate_name: Name of the covariate (e.g., "WT", "CRCL")
        reference_value: Reference/centering value for the covariate
        effect_type: "exponential" or "linear"
    """
    eta_name: str
    covariate_name: str
    reference_value: float = 0.0
    effect_type: str = "exponential"


@dataclass
class EstimationConfig:
    """Configuration for parameter estimation.

    Attributes:
        method: Estimation method ("foce", "saem", or "laplacian")
        theta_init: Initial values for fixed effects
        theta_lower: Lower bounds for theta
        theta_upper: Upper bounds for theta
        theta_names: Names of theta parameters
        omega_init: Initial omega matrix (variance of random effects)
        omega_names: Names of eta parameters
        omega_structure: Structure of omega matrix
        sigma_type: Residual error type
        sigma_init: Initial sigma value
        max_iter: Maximum iterations
        tol: Convergence tolerance
        compute_se: Compute standard errors
        compute_ci: Compute confidence intervals
        ci_level: Confidence level
        verbose: Print progress
        seed: Random seed
        blq_config: BLQ handling configuration
        iov_specs: Inter-occasion variability specifications
        covariate_on_iiv: Covariate effects on IIV
    """
    method: str = "foce"
    theta_init: Optional[List[float]] = None
    theta_lower: Optional[List[float]] = None
    theta_upper: Optional[List[float]] = None
    theta_names: Optional[List[str]] = None
    omega_init: Optional[List[List[float]]] = None
    omega_names: Optional[List[str]] = None
    omega_structure: OmegaStructure = OmegaStructure.DIAGONAL
    sigma_type: str = "proportional"
    sigma_init: float = 0.1
    max_iter: int = 500
    tol: float = 1e-6
    compute_se: bool = True
    compute_ci: bool = True
    ci_level: float = 0.95
    verbose: bool = False
    seed: int = 12345
    # Advanced options
    blq_config: Optional[BLQConfig] = None
    iov_specs: Optional[List[IOVSpec]] = None
    covariate_on_iiv: Optional[List[CovariateOnIIV]] = None
    # FOCE-specific options
    foce_max_inner_iter: int = 50
    foce_inner_tol: float = 1e-6
    foce_centered: bool = False
    foce_compute_robust_se: bool = True
    # SAEM-specific options
    saem_n_burn: int = 300
    saem_n_iter: int = 200
    saem_n_chains: int = 3
    saem_n_mcmc_steps: int = 100


@dataclass
class BootstrapConfig:
    """Configuration for bootstrap analysis.

    Attributes:
        n_bootstrap: Number of bootstrap replicates (FDA recommends ≥500)
        stratify_by: Variables to stratify by (e.g., ["study", "dose"])
        ci_level: Confidence interval level
        ci_method: Method for computing CIs
        parallel: Run in parallel
        seed: Random seed
        min_success_rate: Minimum required success rate
    """
    n_bootstrap: int = 1000
    stratify_by: Optional[List[str]] = None
    ci_level: float = 0.95
    ci_method: BootstrapCIMethod = BootstrapCIMethod.PERCENTILE
    parallel: bool = False
    seed: int = 12345
    min_success_rate: float = 0.8


# ============================================================================
# Result Types
# ============================================================================

@dataclass
class IndividualEstimate:
    """Individual-level estimation results."""
    subject_id: str
    eta: List[float]
    eta_se: Optional[List[float]]
    ipred: List[float]
    pred: List[float]
    cwres: List[float]
    iwres: List[float]
    wres: List[float]
    ofv_contribution: float


@dataclass
class BLQSummary:
    """Summary of BLQ observations."""
    total_observations: int
    blq_observations: int
    blq_percentage: float
    blq_by_subject: Dict[str, int]
    method_used: str


@dataclass
class DiagnosticsSummary:
    """Summary of estimation diagnostics."""
    cwres_mean: float
    cwres_std: float
    iwres_mean: float
    iwres_std: float
    eta_shrinkage: List[float]
    epsilon_shrinkage: float
    condition_number: float
    eigenvalue_ratio: float


@dataclass
class EstimationResult:
    """Result from NLME parameter estimation."""
    method: str
    # Fixed effects
    theta: List[float]
    theta_se: Optional[List[float]]
    theta_se_robust: Optional[List[float]]
    theta_rse: Optional[List[float]]
    theta_ci_lower: Optional[List[float]]
    theta_ci_upper: Optional[List[float]]
    # Random effects
    omega: List[List[float]]
    omega_se: Optional[List[List[float]]]
    omega_corr: List[List[float]]
    # Residual error
    sigma: Dict[str, Any]
    sigma_se: Optional[Dict[str, Any]]
    # Fit statistics
    ofv: float
    aic: float
    bic: float
    # Convergence info
    convergence: bool
    n_iterations: int
    gradient_norm: float
    condition_number: float
    eigenvalue_ratio: float
    covariance_successful: bool
    # Individual results
    individuals: List[IndividualEstimate]
    # Diagnostics
    diagnostics: Optional[DiagnosticsSummary]
    blq_summary: Optional[BLQSummary]
    # Runtime
    runtime_seconds: float
    warnings: List[str]


@dataclass
class BootstrapResult:
    """Result from bootstrap analysis."""
    # Theta estimates
    theta_estimates: np.ndarray  # n_bootstrap x n_params
    theta_mean: List[float]
    theta_se: List[float]
    theta_rse: List[float]
    theta_ci_lower: List[float]
    theta_ci_upper: List[float]
    # Original estimate and bias
    original_estimate: List[float]
    bias: List[float]
    bias_corrected: List[float]
    # Omega and sigma summaries
    omega_mean: Optional[List[List[float]]]
    omega_se: Optional[List[List[float]]]
    omega_ci_lower: Optional[List[List[float]]]
    omega_ci_upper: Optional[List[List[float]]]
    sigma_mean: Optional[float]
    sigma_se: Optional[float]
    sigma_ci_lower: Optional[float]
    sigma_ci_upper: Optional[float]
    # Diagnostics
    n_successful: int
    n_failed: int
    convergence_rate: float
    ci_level: float
    ci_method: str


@dataclass
class ModelComparisonResult:
    """Result from model comparison."""
    model_names: List[str]
    ofv_values: List[float]
    aic_values: List[float]
    bic_values: List[float]
    n_params: List[int]
    best_model_aic: str
    best_model_bic: str
    lrt_statistic: Optional[float]
    lrt_pvalue: Optional[float]
    lrt_df: Optional[int]


# ============================================================================
# Mixture Model Types
# ============================================================================

class MixtureParameterization(Enum):
    """Type of mixture parameterization."""
    FULL = "full"                # Both theta and omega differ across components
    THETA_ONLY = "theta_only"    # Only theta differs (common for metabolizers)
    OMEGA_ONLY = "omega_only"    # Only omega differs


@dataclass
class MixtureComponent:
    """Specification for a single mixture component (subpopulation)."""
    name: str
    theta: Optional[List[float]] = None
    omega: Optional[List[List[float]]] = None
    theta_indices: List[int] = field(default_factory=list)
    omega_indices: List[int] = field(default_factory=list)


@dataclass
class MixtureConfig:
    """Configuration for mixture model estimation."""
    n_components: int
    components: List[MixtureComponent]
    mixing_probabilities: List[float]
    parameterization: MixtureParameterization = MixtureParameterization.THETA_ONLY
    estimate_probabilities: bool = True
    probability_bounds: tuple = (0.01, 0.99)
    # Base parameters
    base_theta: List[float] = field(default_factory=list)
    base_omega: List[List[float]] = field(default_factory=list)
    theta_lower: Optional[List[float]] = None
    theta_upper: Optional[List[float]] = None
    # EM settings
    max_iter: int = 100
    tol: float = 1e-4
    n_init: int = 5
    verbose: bool = False
    # Sigma settings
    sigma_type: str = "proportional"
    sigma_init: float = 0.1
    # Classification
    classification_threshold: float = 0.5
    seed: int = 12345


@dataclass
class MixtureSubjectResult:
    """Result for a single subject in mixture model."""
    subject_id: str
    posterior_probabilities: List[float]
    most_likely_component: int
    classification_confidence: float
    component_likelihoods: List[float]
    eta: List[float]
    ipred: List[float]


@dataclass
class MixtureResult:
    """Result from mixture model estimation."""
    # Convergence info
    converged: bool
    n_iterations: int
    # Fit statistics
    log_likelihood: float
    aic: float
    bic: float
    # Mixing probabilities
    mixing_probabilities: List[float]
    # Component-specific parameters
    component_theta: List[List[float]]
    component_omega: List[List[List[float]]]
    # Subject-level results
    subjects: List[MixtureSubjectResult]
    # Classification summary
    classification_summary: Dict[int, int]
    classification_entropy: float
    # Runtime
    runtime_seconds: float
    messages: List[str]


# ============================================================================
# Main Estimation Function
# ============================================================================

def estimate(
    observed_data: Dict[str, Any],
    model_kind: str,
    config: EstimationConfig,
    grid: Dict[str, Any],
    solver: Optional[Dict[str, Any]] = None,
) -> EstimationResult:
    """
    Run NLME parameter estimation.

    Supports three estimation methods:
    - FOCE-I: First-Order Conditional Estimation with Interaction (gold standard)
    - SAEM: Stochastic Approximation Expectation Maximization
    - Laplacian: Simple Laplacian approximation

    Args:
        observed_data: Dict with keys:
            - subjects: List of subject data dicts
            - Each subject has: subject_id, times, observations, doses
        model_kind: PK model type (e.g., "OneCompIVBolus", "TwoCompOral")
        config: EstimationConfig with initial values and method settings
        grid: Simulation grid dict with t0, t1, saveat
        solver: Optional solver settings

    Returns:
        EstimationResult containing estimated parameters, SEs, and diagnostics

    Example:
        >>> config = EstimationConfig(
        ...     method="foce",
        ...     theta_init=[10.0, 50.0],
        ...     theta_names=["CL", "V"],
        ...     omega_init=[[0.09, 0.0], [0.0, 0.04]],
        ...     sigma_type="proportional",
        ...     sigma_init=0.1
        ... )
        >>> result = estimate(observed_data, "OneCompIVBolus", config, grid)
        >>> print(f"CL = {result.theta[0]} +/- {result.theta_se[0]}")
    """
    jl = _require_julia()

    # Build observed data structure
    SubjectData = jl.NeoPKPD.SubjectData
    ObservedData = jl.NeoPKPD.ObservedData
    DoseEvent = jl.NeoPKPD.DoseEvent

    # Create Julia vector to hold subjects
    subjects_vec = jl.seval("SubjectData[]")
    for subj in observed_data["subjects"]:
        # Create dose vector
        doses_vec = jl.seval("DoseEvent[]")
        for d in subj.get("doses", []):
            dose = DoseEvent(float(d["time"]), float(d["amount"]))
            jl.seval("push!")(doses_vec, dose)

        # Convert times and observations to Julia Vector{Float64}
        times_vec = jl.seval("Float64[]")
        for t in subj["times"]:
            jl.seval("push!")(times_vec, float(t))

        obs_vec = jl.seval("Float64[]")
        for o in subj["observations"]:
            jl.seval("push!")(obs_vec, float(o))

        subj_data = SubjectData(
            subj["subject_id"],
            times_vec,
            obs_vec,
            doses_vec
        )
        jl.seval("push!")(subjects_vec, subj_data)

    obs = ObservedData(subjects_vec)

    # Build model spec (just for type inference - estimation creates individual specs)
    model_spec = _build_base_model_spec(jl, model_kind, config.theta_init)

    # Build estimation method
    if config.method == "foce":
        method = jl.NeoPKPD.FOCEIMethod(
            max_inner_iter=config.foce_max_inner_iter,
            inner_tol=config.foce_inner_tol,
            centered=config.foce_centered,
            compute_robust_se=config.foce_compute_robust_se
        )
    elif config.method == "saem":
        # CRITICAL: parallel_chains=false prevents Bus error in Python-Julia bridge
        # Julia's @threads with closures conflicts with Python's GIL via JuliaCall
        method = jl.NeoPKPD.SAEMMethod(
            n_burn=config.saem_n_burn,
            n_iter=config.saem_n_iter,
            n_chains=config.saem_n_chains,
            n_mcmc_steps=config.saem_n_mcmc_steps,
            parallel_chains=False
        )
    elif config.method == "laplacian":
        method = jl.NeoPKPD.LaplacianMethod()
    else:
        raise ValueError(f"Unknown estimation method: {config.method}")

    # Build sigma specification
    sigma_spec = _build_sigma_spec(jl, config.sigma_type, config.sigma_init)

    # Build omega matrix
    omega_init = np.array(config.omega_init) if config.omega_init else np.eye(len(config.theta_init)) * 0.09

    # Build theta bounds
    n_theta = len(config.theta_init)
    theta_lower = config.theta_lower if config.theta_lower else [1e-6] * n_theta
    theta_upper = config.theta_upper if config.theta_upper else [1e6] * n_theta

    # Build omega structure
    if config.omega_structure == OmegaStructure.DIAGONAL:
        omega_structure = jl.NeoPKPD.DiagonalOmega()
    elif config.omega_structure == OmegaStructure.FULL:
        omega_structure = jl.NeoPKPD.FullOmega()
    else:
        omega_structure = jl.NeoPKPD.DiagonalOmega()

    # Build BLQ config if provided
    blq_config_jl = None
    if config.blq_config is not None:
        blq_method_map = {
            BLQMethod.M1_DISCARD: jl.NeoPKPD.BLQ_M1_DISCARD,
            BLQMethod.M2_IMPUTE_HALF: jl.NeoPKPD.BLQ_M2_IMPUTE_HALF,
            BLQMethod.M2_IMPUTE_ZERO: jl.NeoPKPD.BLQ_M2_IMPUTE_ZERO,
            BLQMethod.M3_LIKELIHOOD: jl.NeoPKPD.BLQ_M3_LIKELIHOOD,
        }
        blq_config_jl = jl.NeoPKPD.BLQConfig(
            blq_method_map[config.blq_config.method],
            float(config.blq_config.lloq),
            max_consecutive_blq=config.blq_config.max_consecutive_blq
        )

    # Build estimation config - convert all Python lists to Julia vectors
    theta_init_jl = _to_julia_float_vec(jl, config.theta_init)
    theta_lower_jl = _to_julia_float_vec(jl, theta_lower)
    theta_upper_jl = _to_julia_float_vec(jl, theta_upper)
    theta_names_jl = _to_julia_symbol_vec(jl, config.theta_names or [f"theta_{i}" for i in range(n_theta)])
    omega_init_jl = _to_julia_matrix(jl, omega_init.tolist())
    omega_names_jl = _to_julia_symbol_vec(jl, config.omega_names or [f"eta_{i}" for i in range(omega_init.shape[0])])

    # Build kwargs dict for EstimationConfig - only include non-default seed if needed
    est_kwargs = dict(
        theta_init=theta_init_jl,
        theta_lower=theta_lower_jl,
        theta_upper=theta_upper_jl,
        theta_names=theta_names_jl,
        omega_init=omega_init_jl,
        omega_names=omega_names_jl,
        omega_structure=omega_structure,
        sigma_init=sigma_spec,
        max_iter=config.max_iter,
        tol=config.tol,
        compute_se=config.compute_se,
        compute_ci=config.compute_ci,
        ci_level=config.ci_level,
        verbose=config.verbose,
    )
    if blq_config_jl is not None:
        est_kwargs["blq_config"] = blq_config_jl

    est_config = jl.NeoPKPD.EstimationConfig(method, **est_kwargs)

    # Build grid
    saveat_jl = _to_julia_float_vec(jl, [float(x) for x in grid["saveat"]])
    grid_jl = jl.NeoPKPD.SimGrid(
        float(grid["t0"]),
        float(grid["t1"]),
        saveat_jl,
    )

    # Build solver
    if solver is None:
        solver_jl = jl.NeoPKPD.SolverSpec(jl.Symbol("Tsit5"), 1e-10, 1e-12, 10**7)
    else:
        solver_jl = jl.NeoPKPD.SolverSpec(
            jl.Symbol(solver.get("alg", "Tsit5")),
            float(solver.get("reltol", 1e-10)),
            float(solver.get("abstol", 1e-12)),
            int(solver.get("maxiters", 10**7)),
        )

    # Run estimation
    result = jl.NeoPKPD.estimate(obs, model_spec, est_config, grid=grid_jl, solver=solver_jl)

    # Convert result to Python
    return _convert_estimation_result(result, config.method)


# ============================================================================
# Bootstrap Function
# ============================================================================

def run_bootstrap(
    observed_data: Dict[str, Any],
    model_kind: str,
    estimation_config: EstimationConfig,
    bootstrap_config: BootstrapConfig,
    grid: Dict[str, Any],
    solver: Optional[Dict[str, Any]] = None,
) -> BootstrapResult:
    """
    Run bootstrap analysis for parameter uncertainty estimation.

    Implements FDA/EMA-recommended non-parametric case bootstrap.

    Args:
        observed_data: Observed data dictionary
        model_kind: PK model type
        estimation_config: Configuration for each estimation run
        bootstrap_config: Bootstrap analysis configuration
        grid: Simulation time grid
        solver: Optional solver settings

    Returns:
        BootstrapResult with bootstrap estimates and confidence intervals

    Example:
        >>> boot_config = BootstrapConfig(n_bootstrap=500, ci_level=0.95)
        >>> result = run_bootstrap(data, "OneCompIVBolus", est_config, boot_config, grid)
        >>> print(f"CL 95% CI: [{result.theta_ci_lower[0]}, {result.theta_ci_upper[0]}]")
    """
    jl = _require_julia()

    # Build observed data
    SubjectData = jl.NeoPKPD.SubjectData
    ObservedData = jl.NeoPKPD.ObservedData
    DoseEvent = jl.NeoPKPD.DoseEvent

    # Create Julia vector to hold subjects
    subjects_vec = jl.seval("SubjectData[]")
    for subj in observed_data["subjects"]:
        # Create dose vector
        doses_vec = jl.seval("DoseEvent[]")
        for d in subj.get("doses", []):
            dose = DoseEvent(float(d["time"]), float(d["amount"]))
            jl.seval("push!")(doses_vec, dose)

        # Convert times and observations to Julia Vector{Float64}
        times_vec = jl.seval("Float64[]")
        for t in subj["times"]:
            jl.seval("push!")(times_vec, float(t))

        obs_vec = jl.seval("Float64[]")
        for o in subj["observations"]:
            jl.seval("push!")(obs_vec, float(o))

        subj_data = SubjectData(
            subj["subject_id"],
            times_vec,
            obs_vec,
            doses_vec
        )
        jl.seval("push!")(subjects_vec, subj_data)

    obs = ObservedData(subjects_vec)

    # Build estimation config (reuse from main estimation)
    model_spec = _build_base_model_spec(jl, model_kind, estimation_config.theta_init)

    # Build estimation method
    if estimation_config.method == "foce":
        method = jl.NeoPKPD.FOCEIMethod()
    elif estimation_config.method == "saem":
        # parallel_chains=False prevents Bus error in Python-Julia bridge
        method = jl.NeoPKPD.SAEMMethod(n_burn=100, n_iter=100, parallel_chains=False)
    else:
        method = jl.NeoPKPD.LaplacianMethod()

    sigma_spec = _build_sigma_spec(jl, estimation_config.sigma_type, estimation_config.sigma_init)
    omega_init = np.array(estimation_config.omega_init) if estimation_config.omega_init else np.eye(len(estimation_config.theta_init)) * 0.09
    n_theta = len(estimation_config.theta_init)

    # Convert to Julia types using helper functions
    theta_init_jl = _to_julia_float_vec(jl, [float(x) for x in estimation_config.theta_init])
    theta_lower_jl = _to_julia_float_vec(jl, [float(x) for x in (estimation_config.theta_lower or [1e-6] * n_theta)])
    theta_upper_jl = _to_julia_float_vec(jl, [float(x) for x in (estimation_config.theta_upper or [1e6] * n_theta)])
    theta_names_jl = _to_julia_symbol_vec(jl, estimation_config.theta_names or [f"theta_{i}" for i in range(n_theta)])
    omega_init_jl = _to_julia_matrix(jl, omega_init.tolist())
    omega_names_jl = _to_julia_symbol_vec(jl, estimation_config.omega_names or [f"eta_{i}" for i in range(omega_init.shape[0])])

    est_config = jl.NeoPKPD.EstimationConfig(
        method,
        theta_init=theta_init_jl,
        theta_lower=theta_lower_jl,
        theta_upper=theta_upper_jl,
        theta_names=theta_names_jl,
        omega_init=omega_init_jl,
        omega_names=omega_names_jl,
        omega_structure=jl.NeoPKPD.DiagonalOmega(),
        sigma_init=sigma_spec,
        max_iter=200,  # Reduced for bootstrap
        tol=1e-4,
        compute_se=False,  # Skip SE computation for speed
        compute_ci=False,
        ci_level=0.95,
        verbose=False,
        # Note: seed is not passed here - Julia uses default; bootstrap seeds are in BootstrapSpec
    )

    # Build grid and solver
    saveat_jl = _to_julia_float_vec(jl, [float(x) for x in grid["saveat"]])
    grid_jl = jl.NeoPKPD.SimGrid(
        float(grid["t0"]),
        float(grid["t1"]),
        saveat_jl,
    )

    if solver is None:
        solver_jl = jl.NeoPKPD.SolverSpec(jl.Symbol("Tsit5"), 1e-8, 1e-10, 10**6)
    else:
        solver_jl = jl.NeoPKPD.SolverSpec(
            jl.Symbol(solver.get("alg", "Tsit5")),
            float(solver.get("reltol", 1e-8)),
            float(solver.get("abstol", 1e-10)),
            int(solver.get("maxiters", 10**6)),
        )

    # Build bootstrap spec
    ci_method_map = {
        BootstrapCIMethod.PERCENTILE: jl.NeoPKPD.PercentileCI(),
        BootstrapCIMethod.BCA: jl.NeoPKPD.BCCI(),
        BootstrapCIMethod.BASIC: jl.NeoPKPD.BasicCI(),
    }

    boot_spec = jl.NeoPKPD.BootstrapSpec(
        n_bootstrap=bootstrap_config.n_bootstrap,
        # Note: seed not passed - Julia uses default UInt64(12345) to avoid Python->Julia Int64/UInt64 issues
        parallel=False,  # Disable parallel to avoid threading issues with JuliaCall
        ci_level=bootstrap_config.ci_level,
        ci_method=ci_method_map[bootstrap_config.ci_method],
        min_success_rate=bootstrap_config.min_success_rate,
    )

    # Use run_bootstrap_native which doesn't require a callback
    # This avoids Python-Julia callback interop issues that cause crashes
    result = jl.NeoPKPD.run_bootstrap_native(
        obs,
        model_spec,
        est_config,
        boot_spec,
        grid=grid_jl,
        solver=solver_jl,
        verbose=False
    )

    return _convert_bootstrap_result(result)


# ============================================================================
# Model Comparison Functions
# ============================================================================

def compare_models(
    results: List[EstimationResult],
    model_names: List[str],
    n_obs: int,
) -> ModelComparisonResult:
    """
    Compare multiple estimation results using AIC, BIC, and LRT.

    Args:
        results: List of EstimationResult objects
        model_names: Names for each model
        n_obs: Total number of observations

    Returns:
        ModelComparisonResult with comparison metrics

    Example:
        >>> comparison = compare_models([result1, result2], ["Full", "Reduced"], n_obs=100)
        >>> print(f"Best model (AIC): {comparison.best_model_aic}")
    """
    ofv_values = [r.ofv for r in results]
    aic_values = [r.aic for r in results]
    bic_values = [r.bic for r in results]
    n_params = [len(r.theta) + _count_omega_params(r.omega) + 1 for r in results]

    best_aic_idx = np.argmin(aic_values)
    best_bic_idx = np.argmin(bic_values)

    # Compute LRT if exactly 2 models (assumes first is full, second is reduced)
    lrt_statistic = None
    lrt_pvalue = None
    lrt_df = None
    if len(results) == 2:
        lrt_df = abs(n_params[0] - n_params[1])
        if lrt_df > 0:
            lrt_statistic = abs(ofv_values[1] - ofv_values[0])
            from scipy.stats import chi2
            lrt_pvalue = 1 - chi2.cdf(lrt_statistic, lrt_df)

    return ModelComparisonResult(
        model_names=model_names,
        ofv_values=ofv_values,
        aic_values=aic_values,
        bic_values=bic_values,
        n_params=n_params,
        best_model_aic=model_names[best_aic_idx],
        best_model_bic=model_names[best_bic_idx],
        lrt_statistic=lrt_statistic,
        lrt_pvalue=lrt_pvalue,
        lrt_df=lrt_df,
    )


def likelihood_ratio_test(
    ofv_full: float,
    ofv_reduced: float,
    df: int,
) -> Dict[str, float]:
    """
    Perform likelihood ratio test between nested models.

    Args:
        ofv_full: OFV of full (more complex) model
        ofv_reduced: OFV of reduced (simpler) model
        df: Degrees of freedom (difference in parameters)

    Returns:
        Dict with chi_sq statistic and p_value

    Example:
        >>> result = likelihood_ratio_test(450.0, 460.0, df=2)
        >>> if result["p_value"] < 0.05:
        ...     print("Full model significantly better")
    """
    chi_sq = max(0.0, ofv_reduced - ofv_full)
    from scipy.stats import chi2
    p_value = 1 - chi2.cdf(chi_sq, df)
    return {"chi_sq": chi_sq, "p_value": p_value, "df": df}


# ============================================================================
# Diagnostic Functions
# ============================================================================

def compute_diagnostics(result: EstimationResult) -> DiagnosticsSummary:
    """
    Compute diagnostic summary from estimation result.

    Args:
        result: EstimationResult from estimate()

    Returns:
        DiagnosticsSummary with residual and shrinkage metrics

    Example:
        >>> diag = compute_diagnostics(result)
        >>> print(f"CWRES mean: {diag.cwres_mean:.3f} (should be ~0)")
        >>> print(f"Eta shrinkage: {[f'{s:.1f}%' for s in diag.eta_shrinkage]}")
    """
    # Collect residuals from all individuals
    all_cwres = []
    all_iwres = []
    all_etas = []

    for ind in result.individuals:
        all_cwres.extend(ind.cwres)
        all_iwres.extend(ind.iwres)
        all_etas.append(ind.eta)

    cwres_arr = np.array(all_cwres)
    iwres_arr = np.array(all_iwres)
    etas_arr = np.array(all_etas)

    # CWRES and IWRES statistics
    cwres_mean = float(np.nanmean(cwres_arr))
    cwres_std = float(np.nanstd(cwres_arr))
    iwres_mean = float(np.nanmean(iwres_arr))
    iwres_std = float(np.nanstd(iwres_arr))

    # Eta shrinkage: 1 - SD(eta) / sqrt(omega)
    n_eta = len(result.omega)
    eta_shrinkage = []
    for j in range(n_eta):
        eta_j = etas_arr[:, j]
        sd_empirical = float(np.nanstd(eta_j))
        sd_theoretical = np.sqrt(result.omega[j][j])
        if sd_theoretical > 0:
            shrinkage = (1.0 - sd_empirical / sd_theoretical) * 100
        else:
            shrinkage = float('nan')
        eta_shrinkage.append(shrinkage)

    # Epsilon shrinkage: 1 - SD(IWRES)
    epsilon_shrinkage = (1.0 - iwres_std) * 100

    return DiagnosticsSummary(
        cwres_mean=cwres_mean,
        cwres_std=cwres_std,
        iwres_mean=iwres_mean,
        iwres_std=iwres_std,
        eta_shrinkage=eta_shrinkage,
        epsilon_shrinkage=epsilon_shrinkage,
        condition_number=result.condition_number,
        eigenvalue_ratio=result.eigenvalue_ratio,
    )


def get_individual_predictions(result: EstimationResult) -> Dict[str, Dict[str, List[float]]]:
    """
    Extract individual predictions from estimation result.

    Args:
        result: EstimationResult from estimate()

    Returns:
        Dict mapping subject_id to dict with ipred, pred, cwres, iwres

    Example:
        >>> preds = get_individual_predictions(result)
        >>> for subj_id, data in preds.items():
        ...     print(f"Subject {subj_id}: IPRED = {data['ipred'][:3]}...")
    """
    predictions = {}
    for ind in result.individuals:
        predictions[ind.subject_id] = {
            "ipred": ind.ipred,
            "pred": ind.pred,
            "cwres": ind.cwres,
            "iwres": ind.iwres,
            "wres": ind.wres,
            "eta": ind.eta,
            "ofv_contribution": ind.ofv_contribution,
        }
    return predictions


# ============================================================================
# NPDE (Normalized Prediction Distribution Errors)
# ============================================================================

@dataclass
class NPDEResult:
    """
    Result of NPDE computation.

    Attributes:
        npde: All NPDE values (pooled across subjects)
        pde: All pde values (prediction distribution errors)
        npde_by_subject: NPDE values per subject
        pde_by_subject: pde values per subject
        n_observations: Total number of observations
        n_simulations: Number of Monte Carlo simulations used
        mean_npde: Mean of NPDE (should be ~0 for good model)
        std_npde: Std of NPDE (should be ~1 for good model)
    """
    npde: List[float]
    pde: List[float]
    npde_by_subject: List[List[float]]
    pde_by_subject: List[List[float]]
    n_observations: int
    n_simulations: int
    mean_npde: float
    std_npde: float


def compute_npde(
    observed_data: Dict[str, Any],
    model_kind: str,
    result: EstimationResult,
    grid: Dict[str, Any],
    solver: Optional[Dict[str, Any]] = None,
    n_sim: int = 1000,
    seed: int = 12345,
) -> NPDEResult:
    """
    Compute NPDE (Normalized Prediction Distribution Errors) using Monte Carlo simulation.

    NPDE is an industry-standard diagnostic for assessing model adequacy.
    If the model correctly describes the data, NPDE should follow a standard
    normal distribution N(0, 1).

    Reference: Brendel et al. (2006), Comets et al. (2008)

    Args:
        observed_data: Observed data dictionary with subjects
        model_kind: PK model type (e.g., "OneCompIVBolus")
        result: EstimationResult from parameter estimation
        grid: Simulation time grid
        solver: Optional solver settings
        n_sim: Number of Monte Carlo simulations (default: 1000, FDA recommends ≥500)
        seed: Random seed for reproducibility

    Returns:
        NPDEResult with NPDE values and diagnostics

    Example:
        >>> npde_result = compute_npde(data, "OneCompIVBolus", est_result, grid)
        >>> print(f"NPDE mean: {npde_result.mean_npde:.3f} (should be ~0)")
        >>> print(f"NPDE std: {npde_result.std_npde:.3f} (should be ~1)")
        >>>
        >>> # Check model adequacy
        >>> if abs(npde_result.mean_npde) < 0.1 and abs(npde_result.std_npde - 1) < 0.2:
        ...     print("Model appears adequate")
    """
    jl = _require_julia()

    # Build observed data
    SubjectData = jl.NeoPKPD.SubjectData
    ObservedData = jl.NeoPKPD.ObservedData
    DoseEvent = jl.NeoPKPD.DoseEvent

    subjects_vec = jl.seval("SubjectData[]")
    for subj in observed_data["subjects"]:
        doses_vec = jl.seval("DoseEvent[]")
        for d in subj.get("doses", []):
            dose = DoseEvent(float(d["time"]), float(d["amount"]))
            jl.seval("push!")(doses_vec, dose)

        times_vec = jl.seval("Float64[]")
        for t in subj["times"]:
            jl.seval("push!")(times_vec, float(t))

        obs_vec = jl.seval("Float64[]")
        for o in subj["observations"]:
            jl.seval("push!")(obs_vec, float(o))

        subj_data = SubjectData(
            subj["subject_id"],
            times_vec,
            obs_vec,
            doses_vec
        )
        jl.seval("push!")(subjects_vec, subj_data)

    obs = ObservedData(subjects_vec)

    # Build model spec
    model_spec = _build_base_model_spec(jl, model_kind, result.theta)

    # Convert theta, omega to Julia
    theta_jl = _to_julia_float_vec(jl, [float(x) for x in result.theta])
    omega_jl = _to_julia_matrix(jl, result.omega)

    # Build sigma spec from result
    sigma_jl = jl.NeoPKPD.ResidualErrorSpec(
        jl.NeoPKPD.ProportionalError(),
        jl.NeoPKPD.ProportionalErrorParams(float(result.sigma))
    )

    # Build grid and solver
    saveat_jl = _to_julia_float_vec(jl, [float(x) for x in grid["saveat"]])
    grid_jl = jl.NeoPKPD.SimGrid(
        float(grid["t0"]),
        float(grid["t1"]),
        saveat_jl,
    )

    if solver is None:
        solver_jl = jl.NeoPKPD.SolverSpec(jl.Symbol("Tsit5"), 1e-8, 1e-10, 10**6)
    else:
        solver_jl = jl.NeoPKPD.SolverSpec(
            jl.Symbol(solver.get("alg", "Tsit5")),
            float(solver.get("reltol", 1e-8)),
            float(solver.get("abstol", 1e-10)),
            int(solver.get("maxiters", 10**6)),
        )

    # Call Julia NPDE function
    npde_jl = jl.NeoPKPD.compute_npde_monte_carlo(
        obs, model_spec, theta_jl, omega_jl, sigma_jl, grid_jl, solver_jl,
        n_sim=n_sim, seed=jl.UInt64(seed)
    )

    # Convert results to Python
    npde = [float(x) for x in npde_jl[jl.Symbol("npde")]]
    pde = [float(x) for x in npde_jl[jl.Symbol("pde")]]

    npde_by_subject = []
    for subj_npde in npde_jl[jl.Symbol("npde_by_subject")]:
        npde_by_subject.append([float(x) for x in subj_npde])

    pde_by_subject = []
    for subj_pde in npde_jl[jl.Symbol("pde_by_subject")]:
        pde_by_subject.append([float(x) for x in subj_pde])

    n_observations = int(npde_jl[jl.Symbol("n_observations")])
    n_simulations = int(npde_jl[jl.Symbol("n_simulations")])

    mean_npde = float(np.mean(npde))
    std_npde = float(np.std(npde))

    return NPDEResult(
        npde=npde,
        pde=pde,
        npde_by_subject=npde_by_subject,
        pde_by_subject=pde_by_subject,
        n_observations=n_observations,
        n_simulations=n_simulations,
        mean_npde=mean_npde,
        std_npde=std_npde,
    )


# ============================================================================
# Helper Functions
# ============================================================================

def _to_julia_float_vec(jl, py_list: List[float]):
    """Convert Python list to Julia Vector{Float64}."""
    vec = jl.seval("Float64[]")
    for x in py_list:
        jl.seval("push!")(vec, float(x))
    return vec


def _to_julia_symbol_vec(jl, py_list: List[str]):
    """Convert Python list of strings to Julia Vector{Symbol}."""
    vec = jl.seval("Symbol[]")
    for s in py_list:
        jl.seval("push!")(vec, jl.Symbol(s))
    return vec


def _to_julia_matrix(jl, py_matrix: List[List[float]]):
    """Convert Python 2D list to Julia Matrix{Float64}."""
    n_rows = len(py_matrix)
    n_cols = len(py_matrix[0]) if n_rows > 0 else 0

    # Create Julia matrix
    mat = jl.seval(f"zeros({n_rows}, {n_cols})")
    for i, row in enumerate(py_matrix):
        for j, val in enumerate(row):
            jl.seval("setindex!")(mat, float(val), i + 1, j + 1)  # Julia is 1-indexed
    return mat


def _build_base_model_spec(jl, model_kind: str, theta_init: List[float]):
    """Build a model spec for the specified model kind."""
    ModelSpec = jl.NeoPKPD.ModelSpec
    DoseEvent = jl.NeoPKPD.DoseEvent

    # Create Julia dose vector
    doses = jl.seval("DoseEvent[]")
    dose = DoseEvent(0.0, 100.0)  # Placeholder
    jl.seval("push!")(doses, dose)

    if model_kind == "OneCompIVBolus":
        Kind = jl.NeoPKPD.OneCompIVBolus
        Params = jl.NeoPKPD.OneCompIVBolusParams
        params = Params(float(theta_init[0]), float(theta_init[1]))
    elif model_kind == "OneCompOralFirstOrder":
        Kind = jl.NeoPKPD.OneCompOralFirstOrder
        Params = jl.NeoPKPD.OneCompOralFirstOrderParams
        params = Params(float(theta_init[0]), float(theta_init[1]), float(theta_init[2]))
    elif model_kind == "TwoCompIVBolus":
        Kind = jl.NeoPKPD.TwoCompIVBolus
        Params = jl.NeoPKPD.TwoCompIVBolusParams
        params = Params(float(theta_init[0]), float(theta_init[1]), float(theta_init[2]), float(theta_init[3]))
    elif model_kind == "TwoCompOral":
        Kind = jl.NeoPKPD.TwoCompOral
        Params = jl.NeoPKPD.TwoCompOralParams
        params = Params(*[float(x) for x in theta_init[:5]])
    else:
        raise ValueError(f"Unsupported model kind: {model_kind}")

    return ModelSpec(Kind(), "estimation", params, doses)


def _build_sigma_spec(jl, sigma_type: str, sigma_init: float):
    """Build residual error specification."""
    if sigma_type == "additive":
        kind = jl.NeoPKPD.AdditiveError()
        params = jl.NeoPKPD.AdditiveErrorParams(float(sigma_init))
    elif sigma_type == "proportional":
        kind = jl.NeoPKPD.ProportionalError()
        params = jl.NeoPKPD.ProportionalErrorParams(float(sigma_init))
    elif sigma_type == "combined":
        kind = jl.NeoPKPD.CombinedError()
        params = jl.NeoPKPD.CombinedErrorParams(float(sigma_init), float(sigma_init))
    elif sigma_type == "exponential":
        kind = jl.NeoPKPD.ExponentialError()
        params = jl.NeoPKPD.ExponentialErrorParams(float(sigma_init))
    else:
        raise ValueError(f"Unknown sigma type: {sigma_type}")

    return jl.NeoPKPD.ResidualErrorSpec(kind, params, jl.Symbol("conc"), jl.UInt64(1))


def _convert_estimation_result(result, method: str) -> EstimationResult:
    """Convert Julia EstimationResult to Python dataclass."""
    # Convert individuals
    individuals = []
    for ind in result.individuals:
        individuals.append(IndividualEstimate(
            subject_id=str(ind.subject_id),
            eta=list(ind.eta),
            eta_se=list(ind.eta_se) if ind.eta_se is not None else None,
            ipred=list(ind.ipred),
            pred=list(ind.pred),
            cwres=list(ind.cwres),
            iwres=list(ind.iwres),
            wres=list(ind.wres),
            ofv_contribution=float(ind.ofv_contribution),
        ))

    # Convert omega to list of lists
    omega = [[float(result.omega[i, j]) for j in range(result.omega.shape[1])]
             for i in range(result.omega.shape[0])]
    omega_corr = [[float(result.omega_corr[i, j]) for j in range(result.omega_corr.shape[1])]
                  for i in range(result.omega_corr.shape[0])]

    omega_se = None
    if result.omega_se is not None:
        omega_se = [[float(result.omega_se[i, j]) for j in range(result.omega_se.shape[1])]
                    for i in range(result.omega_se.shape[0])]

    # Get robust SE if available
    theta_se_robust = None
    if hasattr(result, 'theta_se_robust') and result.theta_se_robust is not None:
        theta_se_robust = list(result.theta_se_robust)

    # Get BLQ summary if available
    blq_summary = None
    if hasattr(result, 'blq_summary') and result.blq_summary is not None:
        blq = result.blq_summary
        blq_summary = BLQSummary(
            total_observations=int(blq.total_observations),
            blq_observations=int(blq.blq_observations),
            blq_percentage=float(blq.blq_percentage),
            blq_by_subject={str(k): int(v) for k, v in blq.blq_by_subject.items()},
            method_used=str(blq.method_used),
        )

    return EstimationResult(
        method=method,
        theta=list(result.theta),
        theta_se=list(result.theta_se) if result.theta_se is not None else None,
        theta_se_robust=theta_se_robust,
        theta_rse=list(result.theta_rse) if result.theta_rse is not None else None,
        theta_ci_lower=list(result.theta_ci_lower) if result.theta_ci_lower is not None else None,
        theta_ci_upper=list(result.theta_ci_upper) if result.theta_ci_upper is not None else None,
        omega=omega,
        omega_se=omega_se,
        omega_corr=omega_corr,
        sigma={"type": "estimated"},
        sigma_se=None,
        ofv=float(result.ofv),
        aic=float(result.aic),
        bic=float(result.bic),
        convergence=bool(result.convergence),
        n_iterations=int(result.n_iterations),
        gradient_norm=float(result.gradient_norm),
        condition_number=float(result.condition_number) if hasattr(result, 'condition_number') else float('nan'),
        eigenvalue_ratio=float(result.eigenvalue_ratio) if hasattr(result, 'eigenvalue_ratio') else float('nan'),
        covariance_successful=bool(result.covariance_step_successful) if hasattr(result, 'covariance_step_successful') else bool(result.covariance_successful),
        individuals=individuals,
        diagnostics=None,  # Computed separately
        blq_summary=blq_summary,
        runtime_seconds=float(result.runtime_seconds),
        warnings=list(result.messages) if hasattr(result, 'messages') else [],
    )


def _convert_bootstrap_result(result) -> BootstrapResult:
    """Convert Julia BootstrapResult to Python dataclass."""
    # Convert theta estimates matrix
    theta_estimates = np.array(result.theta_estimates)

    # Omega summary
    omega_mean = None
    omega_se = None
    omega_ci_lower = None
    omega_ci_upper = None
    if result.omega_summary is not None:
        os = result.omega_summary
        omega_mean = [[float(os.mean[i, j]) for j in range(os.mean.shape[1])]
                      for i in range(os.mean.shape[0])]
        omega_se = [[float(os.se[i, j]) for j in range(os.se.shape[1])]
                    for i in range(os.se.shape[0])]
        omega_ci_lower = [[float(os.ci_lower[i, j]) for j in range(os.ci_lower.shape[1])]
                          for i in range(os.ci_lower.shape[0])]
        omega_ci_upper = [[float(os.ci_upper[i, j]) for j in range(os.ci_upper.shape[1])]
                          for i in range(os.ci_upper.shape[0])]

    # Sigma summary
    sigma_mean = None
    sigma_se = None
    sigma_ci_lower = None
    sigma_ci_upper = None
    if result.sigma_summary is not None:
        ss = result.sigma_summary
        sigma_mean = float(ss.mean)
        sigma_se = float(ss.se)
        sigma_ci_lower = float(ss.ci_lower)
        sigma_ci_upper = float(ss.ci_upper)

    return BootstrapResult(
        theta_estimates=theta_estimates,
        theta_mean=list(result.theta_mean),
        theta_se=list(result.theta_se),
        theta_rse=list(result.theta_rse),
        theta_ci_lower=list(result.theta_ci_lower),
        theta_ci_upper=list(result.theta_ci_upper),
        original_estimate=list(result.original_estimate),
        bias=list(result.bias),
        bias_corrected=list(result.bias_corrected),
        omega_mean=omega_mean,
        omega_se=omega_se,
        omega_ci_lower=omega_ci_lower,
        omega_ci_upper=omega_ci_upper,
        sigma_mean=sigma_mean,
        sigma_se=sigma_se,
        sigma_ci_lower=sigma_ci_lower,
        sigma_ci_upper=sigma_ci_upper,
        n_successful=int(result.diagnostics.n_successful),
        n_failed=int(result.diagnostics.n_failed),
        convergence_rate=float(result.diagnostics.convergence_rate),
        ci_level=float(result.ci_level),
        ci_method=str(result.ci_method),
    )


def _count_omega_params(omega: List[List[float]]) -> int:
    """Count number of parameters in omega matrix (diagonal only)."""
    return len(omega)


# ============================================================================
# Stepwise Covariate Modeling (SCM)
# ============================================================================

class CovariateRelationshipType(str, Enum):
    """Type of covariate-parameter relationship."""
    LINEAR = "linear"
    POWER = "power"
    EXPONENTIAL = "exponential"


@dataclass
class CovariateRelationship:
    """
    Defines a potential covariate-parameter relationship to test in SCM.

    Attributes:
        param: Parameter affected by covariate (e.g., "CL", "V")
        covariate: Covariate name (e.g., "WT", "AGE")
        relationship: Type of relationship (linear, power, exponential)
        reference: Reference value for centering (e.g., median weight)
    """
    param: str
    covariate: str
    relationship: CovariateRelationshipType = CovariateRelationshipType.POWER
    reference: float = 1.0


@dataclass
class SCMConfig:
    """
    Configuration for stepwise covariate modeling.

    Attributes:
        relationships: List of covariate relationships to test
        forward_p_value: P-value threshold for forward inclusion (default: 0.05)
        backward_p_value: P-value threshold for backward elimination (default: 0.01)
        max_iterations: Maximum iterations for forward/backward steps
    """
    relationships: List[CovariateRelationship]
    forward_p_value: float = 0.05
    backward_p_value: float = 0.01
    max_iterations: int = 100


@dataclass
class SCMStepResult:
    """Result of a single SCM step."""
    param: str
    covariate: str
    relationship: str
    ofv_change: float
    p_value: float
    coefficient: float
    coefficient_se: float
    included: bool
    step_type: str  # "forward" or "backward"


@dataclass
class SCMResult:
    """
    Complete result of stepwise covariate modeling.

    Attributes:
        final_model: List of included covariate relationships
        forward_steps: Results from forward selection
        backward_steps: Results from backward elimination
        base_ofv: OFV of base model (no covariates)
        final_ofv: OFV of final model
        n_forward_iterations: Number of forward selection iterations
        n_backward_iterations: Number of backward elimination iterations
    """
    final_model: List[CovariateRelationship]
    forward_steps: List[SCMStepResult]
    backward_steps: List[SCMStepResult]
    base_ofv: float
    final_ofv: float
    n_forward_iterations: int
    n_backward_iterations: int


def run_scm(
    observed_data: Dict[str, Any],
    model_kind: str,
    estimation_config: EstimationConfig,
    scm_config: SCMConfig,
    grid: Dict[str, Any],
    solver: Optional[Dict[str, Any]] = None,
) -> SCMResult:
    """
    Run stepwise covariate modeling (SCM) analysis.

    Performs forward selection and backward elimination to identify
    significant covariate effects on PK parameters.

    Args:
        observed_data: Observed data dictionary with subjects containing covariates
        model_kind: PK model type (e.g., "OneCompIVBolus")
        estimation_config: Configuration for base estimation
        scm_config: SCM configuration with relationships to test
        grid: Simulation time grid
        solver: Optional solver settings

    Returns:
        SCMResult with selected covariates, p-values, and OFV changes

    Example:
        >>> relationships = [
        ...     CovariateRelationship("CL", "WT", CovariateRelationshipType.POWER, 70.0),
        ...     CovariateRelationship("V", "WT", CovariateRelationshipType.POWER, 70.0),
        ...     CovariateRelationship("CL", "AGE", CovariateRelationshipType.LINEAR, 40.0),
        ... ]
        >>> scm_config = SCMConfig(relationships, forward_p_value=0.05)
        >>> result = run_scm(data, "OneCompIVBolus", est_config, scm_config, grid)
        >>> print(f"Selected covariates: {len(result.final_model)}")
    """
    jl = _require_julia()

    # Build observed data with covariates
    SubjectData = jl.NeoPKPD.SubjectData
    ObservedData = jl.NeoPKPD.ObservedData
    DoseEvent = jl.NeoPKPD.DoseEvent

    # Create Julia vector to hold subjects
    subjects_vec = jl.seval("SubjectData[]")
    for subj in observed_data["subjects"]:
        # Create dose vector
        doses_vec = jl.seval("DoseEvent[]")
        for d in subj.get("doses", []):
            dose = DoseEvent(float(d["time"]), float(d["amount"]))
            jl.seval("push!")(doses_vec, dose)

        # Convert times and observations to Julia Vector{Float64}
        times_vec = jl.seval("Float64[]")
        for t in subj["times"]:
            jl.seval("push!")(times_vec, float(t))

        obs_vec = jl.seval("Float64[]")
        for o in subj["observations"]:
            jl.seval("push!")(obs_vec, float(o))

        # Convert covariates to Julia Dict{Symbol, Any}
        cov_dict = jl.seval("Dict{Symbol, Any}()")
        for key, val in subj.get("covariates", {}).items():
            jl.seval("setindex!")(cov_dict, float(val), jl.Symbol(key))

        subj_data = SubjectData(
            subj["subject_id"],
            times_vec,
            obs_vec,
            doses_vec,
            covariates=cov_dict
        )
        jl.seval("push!")(subjects_vec, subj_data)

    obs = ObservedData(subjects_vec)

    # Build model spec
    model_spec = _build_base_model_spec(jl, model_kind, estimation_config.theta_init)

    # Build estimation config
    n_theta = len(estimation_config.theta_init)
    omega_init = np.array(estimation_config.omega_init) if estimation_config.omega_init else np.eye(n_theta) * 0.09

    if estimation_config.method == "foce":
        method = jl.NeoPKPD.FOCEIMethod()
    elif estimation_config.method == "saem":
        method = jl.NeoPKPD.SAEMMethod(n_burn=100, n_iter=100, parallel_chains=False)
    else:
        method = jl.NeoPKPD.LaplacianMethod()

    sigma_spec = _build_sigma_spec(jl, estimation_config.sigma_type, estimation_config.sigma_init)

    theta_init_jl = _to_julia_float_vec(jl, [float(x) for x in estimation_config.theta_init])
    theta_lower_jl = _to_julia_float_vec(jl, [float(x) for x in (estimation_config.theta_lower or [1e-6] * n_theta)])
    theta_upper_jl = _to_julia_float_vec(jl, [float(x) for x in (estimation_config.theta_upper or [1e6] * n_theta)])
    theta_names_jl = _to_julia_symbol_vec(jl, estimation_config.theta_names or [f"theta_{i}" for i in range(n_theta)])
    omega_init_jl = _to_julia_matrix(jl, omega_init.tolist())
    omega_names_jl = _to_julia_symbol_vec(jl, estimation_config.omega_names or [f"eta_{i}" for i in range(omega_init.shape[0])])

    est_config = jl.NeoPKPD.EstimationConfig(
        method,
        theta_init=theta_init_jl,
        theta_lower=theta_lower_jl,
        theta_upper=theta_upper_jl,
        theta_names=theta_names_jl,
        omega_init=omega_init_jl,
        omega_names=omega_names_jl,
        omega_structure=jl.NeoPKPD.DiagonalOmega(),
        sigma_init=sigma_spec,
        max_iter=estimation_config.max_iter,
        tol=estimation_config.tol,
        compute_se=True,
        compute_ci=False,
        verbose=False,
    )

    # Build grid and solver
    saveat_jl = _to_julia_float_vec(jl, [float(x) for x in grid["saveat"]])
    grid_jl = jl.NeoPKPD.SimGrid(
        float(grid["t0"]),
        float(grid["t1"]),
        saveat_jl,
    )

    if solver is None:
        solver_jl = jl.NeoPKPD.SolverSpec(jl.Symbol("Tsit5"), 1e-8, 1e-10, 10**6)
    else:
        solver_jl = jl.NeoPKPD.SolverSpec(
            jl.Symbol(solver.get("alg", "Tsit5")),
            float(solver.get("reltol", 1e-8)),
            float(solver.get("abstol", 1e-10)),
            int(solver.get("maxiters", 10**6)),
        )

    # Build covariate relationships
    relationships_vec = jl.seval("CovariateRelationship[]")
    for rel in scm_config.relationships:
        jl_rel = jl.NeoPKPD.CovariateRelationship(
            jl.Symbol(rel.param),
            jl.Symbol(rel.covariate),
            relationship=jl.Symbol(rel.relationship.value),
            reference=float(rel.reference)
        )
        jl.seval("push!")(relationships_vec, jl_rel)

    # Build SCM spec
    scm_spec = jl.NeoPKPD.SCMSpec(
        relationships_vec,
        forward_p_value=scm_config.forward_p_value,
        backward_p_value=scm_config.backward_p_value,
        max_iterations=scm_config.max_iterations,
    )

    # Run SCM
    result = jl.NeoPKPD.run_scm_native(
        obs,
        model_spec,
        est_config,
        scm_spec,
        grid=grid_jl,
        solver=solver_jl,
        verbose=False
    )

    # Convert result to Python types
    final_model = []
    for rel in result.final_model:
        final_model.append(CovariateRelationship(
            param=str(rel.param),
            covariate=str(rel.covariate),
            relationship=CovariateRelationshipType(str(rel.relationship)),
            reference=float(rel.reference)
        ))

    forward_steps = []
    for step in result.forward_steps:
        forward_steps.append(SCMStepResult(
            param=str(step.relationship.param),
            covariate=str(step.relationship.covariate),
            relationship=str(step.relationship.relationship),
            ofv_change=float(step.ofv_change),
            p_value=float(step.p_value),
            coefficient=float(step.coefficient),
            coefficient_se=float(step.coefficient_se),
            included=bool(step.included),
            step_type=str(step.step_type)
        ))

    backward_steps = []
    for step in result.backward_steps:
        backward_steps.append(SCMStepResult(
            param=str(step.relationship.param),
            covariate=str(step.relationship.covariate),
            relationship=str(step.relationship.relationship),
            ofv_change=float(step.ofv_change),
            p_value=float(step.p_value),
            coefficient=float(step.coefficient),
            coefficient_se=float(step.coefficient_se),
            included=bool(step.included),
            step_type=str(step.step_type)
        ))

    return SCMResult(
        final_model=final_model,
        forward_steps=forward_steps,
        backward_steps=backward_steps,
        base_ofv=float(result.base_ofv),
        final_ofv=float(result.final_ofv),
        n_forward_iterations=int(result.n_forward_iterations),
        n_backward_iterations=int(result.n_backward_iterations)
    )


# ============================================================================
# Mixture Model Estimation
# ============================================================================

def estimate_mixture(
    observed_data: Dict[str, Any],
    model_kind: str,
    config: MixtureConfig,
    grid: Dict[str, Any],
    solver: Optional[Dict[str, Any]] = None,
) -> MixtureResult:
    """
    Run mixture model estimation for subpopulation analysis.

    Industry-standard implementation following NONMEM $MIX and Monolix methodology.
    Uses EM algorithm with proper individual eta estimation and omega updates.

    Args:
        observed_data: Dict with keys:
            - subjects: List of subject data dicts
            - Each subject has: subject_id, times, observations, doses
        model_kind: PK model type (e.g., "OneCompIVBolus", "TwoCompOral")
        config: MixtureConfig with mixture specification and estimation settings
        grid: Simulation grid dict with t0, t1, saveat
        solver: Optional solver settings

    Returns:
        MixtureResult containing:
            - converged: Whether estimation converged
            - mixing_probabilities: Estimated π_k for each component
            - component_theta: Component-specific theta values
            - component_omega: Component-specific omega matrices
            - subjects: Per-subject posterior probabilities and classifications
            - classification_entropy: Measure of separation quality

    Example:
        >>> # Two-component mixture for metabolizer phenotypes
        >>> config = MixtureConfig(
        ...     n_components=2,
        ...     components=[
        ...         MixtureComponent(name="Slow", theta=[2.0, 50.0]),
        ...         MixtureComponent(name="Fast", theta=[6.0, 50.0])
        ...     ],
        ...     mixing_probabilities=[0.7, 0.3],
        ...     base_theta=[4.0, 50.0],
        ...     base_omega=[[0.09, 0.0], [0.0, 0.04]],
        ... )
        >>> result = estimate_mixture(data, "OneCompIVBolus", config, grid)
        >>> for subj in result.subjects:
        ...     print(f"{subj.subject_id}: Component {subj.most_likely_component}")
    """
    jl = _require_julia()

    # Build observed data structure
    SubjectData = jl.NeoPKPD.SubjectData
    ObservedData = jl.NeoPKPD.ObservedData
    DoseEvent = jl.NeoPKPD.DoseEvent

    # Create Julia vector to hold subjects
    subjects_vec = jl.seval("SubjectData[]")
    for subj in observed_data["subjects"]:
        doses_vec = jl.seval("DoseEvent[]")
        for d in subj.get("doses", []):
            dose = DoseEvent(float(d["time"]), float(d["amount"]))
            jl.seval("push!")(doses_vec, dose)

        times_vec = jl.seval("Float64[]")
        for t in subj["times"]:
            jl.seval("push!")(times_vec, float(t))

        obs_vec = jl.seval("Float64[]")
        for o in subj["observations"]:
            jl.seval("push!")(obs_vec, float(o))

        subj_data = SubjectData(
            subj["subject_id"],
            times_vec,
            obs_vec,
            doses_vec
        )
        jl.seval("push!")(subjects_vec, subj_data)

    obs = ObservedData(subjects_vec)

    # Build model spec
    model_spec = _build_base_model_spec(jl, model_kind, config.base_theta)

    # Build mixture components
    MixtureComponent_jl = jl.NeoPKPD.MixtureComponent
    components_vec = jl.seval("MixtureComponent[]")

    for comp in config.components:
        theta_jl = jl.seval("nothing")
        if comp.theta is not None:
            theta_jl = jl.seval("Float64[]")
            for v in comp.theta:
                jl.seval("push!")(theta_jl, float(v))

        omega_jl = jl.seval("nothing")
        if comp.omega is not None:
            omega_jl = _to_julia_matrix(jl, comp.omega)

        theta_indices_jl = jl.seval("Int[]")
        for idx in comp.theta_indices:
            jl.seval("push!")(theta_indices_jl, int(idx))

        omega_indices_jl = jl.seval("Int[]")
        for idx in comp.omega_indices:
            jl.seval("push!")(omega_indices_jl, int(idx))

        jl_comp = MixtureComponent_jl(
            name=comp.name,
            theta=theta_jl,
            omega=omega_jl,
            theta_indices=theta_indices_jl,
            omega_indices=omega_indices_jl
        )
        jl.seval("push!")(components_vec, jl_comp)

    # Build mixing probabilities
    mixing_probs_jl = jl.seval("Float64[]")
    for p in config.mixing_probabilities:
        jl.seval("push!")(mixing_probs_jl, float(p))

    # Build parameterization type
    if config.parameterization == MixtureParameterization.FULL:
        param_type = jl.NeoPKPD.FullMixture()
    elif config.parameterization == MixtureParameterization.THETA_ONLY:
        param_type = jl.NeoPKPD.ThetaOnlyMixture()
    else:
        param_type = jl.NeoPKPD.OmegaOnlyMixture()

    # Build mixture spec
    MixtureSpec = jl.NeoPKPD.MixtureSpec
    mixture_spec = MixtureSpec(
        n_components=config.n_components,
        components=components_vec,
        mixing_probabilities=mixing_probs_jl,
        parameterization=param_type,
        estimate_probabilities=config.estimate_probabilities,
        probability_bounds=(float(config.probability_bounds[0]), float(config.probability_bounds[1]))
    )

    # Build EM method
    MixtureEM = jl.NeoPKPD.MixtureEM
    FOCEIMethod = jl.NeoPKPD.FOCEIMethod
    em_method = MixtureEM(
        max_iter=config.max_iter,
        tol=config.tol,
        inner_method=FOCEIMethod(),
        n_init=config.n_init,
        verbose=config.verbose
    )

    # Build base theta and omega
    base_theta_jl = jl.seval("Float64[]")
    for v in config.base_theta:
        jl.seval("push!")(base_theta_jl, float(v))

    base_omega_jl = _to_julia_matrix(jl, config.base_omega)

    # Build sigma
    sigma_jl = _build_sigma_spec(jl, config.sigma_type, config.sigma_init)

    # Build theta bounds
    n_theta = len(config.base_theta)
    theta_lower_jl = jl.seval("Float64[]")
    theta_upper_jl = jl.seval("Float64[]")

    lower_vals = config.theta_lower if config.theta_lower else [0.0] * n_theta
    upper_vals = config.theta_upper if config.theta_upper else [float("inf")] * n_theta

    for v in lower_vals:
        jl.seval("push!")(theta_lower_jl, float(v))
    for v in upper_vals:
        jl.seval("push!")(theta_upper_jl, float(v))

    # Build mixture config
    MixtureConfig_jl = jl.NeoPKPD.MixtureConfig
    mix_config = MixtureConfig_jl(
        mixture_spec,
        em_method,
        base_theta=base_theta_jl,
        base_omega=base_omega_jl,
        sigma=sigma_jl,
        theta_lower=theta_lower_jl,
        theta_upper=theta_upper_jl,
        compute_posteriors=True,
        classification_threshold=config.classification_threshold,
        seed=jl.UInt64(config.seed)
    )

    # Build grid and solver
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    saveat_jl = jl.seval("Float64[]")
    for t in grid.get("saveat", []):
        jl.seval("push!")(saveat_jl, float(t))

    grid_jl = SimGrid(
        float(grid.get("t0", 0.0)),
        float(grid.get("t1", 24.0)),
        saveat_jl
    )

    solver_settings = solver or {}
    solver_jl = SolverSpec(
        jl.Symbol(solver_settings.get("algorithm", "Tsit5")),
        float(solver_settings.get("abstol", 1e-8)),
        float(solver_settings.get("reltol", 1e-8)),
        int(solver_settings.get("maxiters", 100000))
    )

    # Run mixture estimation
    result = jl.NeoPKPD.mixture_estimate(
        obs, model_spec, mix_config,
        grid=grid_jl, solver=solver_jl
    )

    # Convert results to Python
    subjects = []
    for subj in result.subjects:
        subjects.append(MixtureSubjectResult(
            subject_id=str(subj.subject_id),
            posterior_probabilities=list(subj.posterior_probabilities),
            most_likely_component=int(subj.most_likely_component),
            classification_confidence=float(subj.classification_confidence),
            component_likelihoods=list(subj.component_likelihoods),
            eta=list(subj.eta),
            ipred=list(subj.ipred)
        ))

    component_theta = []
    for theta in result.component_theta:
        component_theta.append(list(theta))

    component_omega = []
    for omega in result.component_omega:
        omega_list = []
        n = jl.seval("size")(omega, 1)
        for i in range(1, n + 1):
            row = []
            for j in range(1, n + 1):
                row.append(float(omega[i, j]))
            omega_list.append(row)
        component_omega.append(omega_list)

    # Convert classification summary
    classification_summary = {}
    for k in result.classification_summary:
        classification_summary[int(k)] = int(result.classification_summary[k])

    return MixtureResult(
        converged=bool(result.converged),
        n_iterations=int(result.n_iterations),
        log_likelihood=float(result.log_likelihood),
        aic=float(result.aic),
        bic=float(result.bic),
        mixing_probabilities=list(result.mixing_probabilities),
        component_theta=component_theta,
        component_omega=component_omega,
        subjects=subjects,
        classification_summary=classification_summary,
        classification_entropy=float(result.entropy),
        runtime_seconds=float(result.runtime_seconds),
        messages=list(result.messages)
    )


__all__ = [
    # Main functions
    "estimate",
    "estimate_mixture",
    "run_bootstrap",
    "run_scm",
    "compare_models",
    "likelihood_ratio_test",
    "compute_diagnostics",
    "get_individual_predictions",
    "compute_npde",
    # Configuration types
    "EstimationConfig",
    "BootstrapConfig",
    "BLQConfig",
    "IOVSpec",
    "CovariateOnIIV",
    "SCMConfig",
    "CovariateRelationship",
    "MixtureConfig",
    "MixtureComponent",
    # Result types
    "EstimationResult",
    "BootstrapResult",
    "IndividualEstimate",
    "DiagnosticsSummary",
    "BLQSummary",
    "ModelComparisonResult",
    "SCMResult",
    "SCMStepResult",
    "NPDEResult",
    "MixtureResult",
    "MixtureSubjectResult",
    # Enums
    "BLQMethod",
    "OmegaStructure",
    "BootstrapCIMethod",
    "CovariateRelationshipType",
    "MixtureParameterization",
]
