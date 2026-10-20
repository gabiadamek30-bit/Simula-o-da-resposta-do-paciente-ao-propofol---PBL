"""
NeoPKPD NCA Bioequivalence - FDA/EMA compliant bioequivalence analysis.

This module provides bioequivalence assessment tools including:
- 90% confidence interval calculation
- TOST (Two One-Sided Tests) analysis
- Geometric mean ratio calculation
- Within-subject CV estimation
- Reference-scaled average bioequivalence (RSABE) for highly variable drugs
- Study design types (Crossover2x2, Crossover2x4, ParallelGroup)
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union

from .._core import _require_julia, _to_julia_float_vector


# ============================================================================
# Study Design Types
# ============================================================================

class BEStudyDesign(str, Enum):
    """Bioequivalence study design types."""
    CROSSOVER_2X2 = "crossover_2x2"
    CROSSOVER_2X4 = "crossover_2x4"
    PARALLEL_GROUP = "parallel_group"


class ReplicateDesign(str, Enum):
    """Replicate design types for reference-scaled BE."""
    PARTIAL_REPLICATE_3X3 = "partial_replicate_3x3"
    FULL_REPLICATE_2X4 = "full_replicate_2x4"
    FULL_REPLICATE_2X3 = "full_replicate_2x3"


class RegulatoryGuidance(str, Enum):
    """Regulatory guidance for reference-scaled BE."""
    FDA = "fda"
    EMA = "ema"
    HEALTH_CANADA = "health_canada"


# ============================================================================
# Geometric Mean and Ratio
# ============================================================================

def geometric_mean(values: List[float]) -> float:
    """
    Calculate geometric mean.

    GM = exp(mean(log(values)))

    Args:
        values: List of positive values

    Returns:
        Geometric mean

    Example:
        >>> gm = geometric_mean([10, 20, 30])
        >>> print(f"Geometric mean: {gm:.2f}")
    """
    jl = _require_julia()
    v = _to_julia_float_vector(jl, values)
    return float(jl.NeoPKPD.geometric_mean(v))


def geometric_mean_ratio(
    test_values: List[float],
    reference_values: List[float]
) -> float:
    """
    Calculate geometric mean ratio (GMR) of test to reference.

    GMR = exp(mean(log(test)) - mean(log(reference)))

    Args:
        test_values: Test formulation values (Cmax or AUC)
        reference_values: Reference formulation values

    Returns:
        Geometric mean ratio (test/reference)

    Example:
        >>> gmr = geometric_mean_ratio(test_cmax, reference_cmax)
        >>> print(f"GMR: {gmr:.4f}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    return float(jl.NeoPKPD.geometric_mean_ratio(t, r))


def within_subject_cv(
    test_values: List[float],
    reference_values: List[float]
) -> float:
    """
    Estimate within-subject coefficient of variation from crossover data.

    CV_intra = sqrt(exp(MSE) - 1) * 100%

    Where MSE is from ANOVA on log-transformed data.

    Args:
        test_values: Test values (paired with reference)
        reference_values: Reference values (paired)

    Returns:
        Within-subject CV (%)

    Example:
        >>> cv = within_subject_cv(test_cmax, reference_cmax)
        >>> print(f"Intra-subject CV: {cv:.1f}%")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    return float(jl.NeoPKPD.within_subject_cv(t, r))


# ============================================================================
# 90% Confidence Interval
# ============================================================================

def bioequivalence_90ci(
    test_values: List[float],
    reference_values: List[float],
    log_transform: bool = True
) -> Dict[str, float]:
    """
    Calculate 90% confidence interval for geometric mean ratio.

    Uses the standard two-sequence, two-period crossover analysis.

    Args:
        test_values: Test formulation values (paired with reference)
        reference_values: Reference formulation values
        log_transform: Apply log transformation (default: True)

    Returns:
        Dict with keys:
        - gmr: Geometric mean ratio
        - ci_lower: Lower bound of 90% CI
        - ci_upper: Upper bound of 90% CI
        - cv_intra: Intra-subject CV (%)
        - n: Number of subjects

    Example:
        >>> result = bioequivalence_90ci(test_cmax, reference_cmax)
        >>> print(f"GMR: {result['gmr']:.4f}")
        >>> print(f"90% CI: ({result['ci_lower']:.4f}, {result['ci_upper']:.4f})")
        >>> if result['ci_lower'] >= 0.80 and result['ci_upper'] <= 1.25:
        ...     print("Bioequivalent!")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    result = jl.NeoPKPD.bioequivalence_90ci(t, r, log_transform=log_transform)

    return {
        "gmr": float(result.gmr),
        "ci_lower": float(result.ci_lower),
        "ci_upper": float(result.ci_upper),
        "cv_intra": float(result.cv_intra),
        "n": int(result.n),
    }


# ============================================================================
# TOST Analysis
# ============================================================================

def tost_analysis(
    test_values: List[float],
    reference_values: List[float],
    theta_lower: float = 0.80,
    theta_upper: float = 1.25,
    alpha: float = 0.05
) -> Dict[str, Any]:
    """
    Perform Two One-Sided Tests (TOST) procedure for bioequivalence.

    Tests:
    - H01: muT/muR <= theta_lower (lower bound)
    - H02: muT/muR >= theta_upper (upper bound)

    BE is concluded if both null hypotheses are rejected.

    Args:
        test_values: Test formulation values
        reference_values: Reference formulation values
        theta_lower: Lower equivalence bound (default: 0.80)
        theta_upper: Upper equivalence bound (default: 1.25)
        alpha: Significance level (default: 0.05)

    Returns:
        Dict with keys:
        - parameter: Parameter analyzed
        - t_lower: T-statistic for lower bound test
        - t_upper: T-statistic for upper bound test
        - p_lower: P-value for lower bound test
        - p_upper: P-value for upper bound test
        - reject_lower: Whether lower bound H0 is rejected
        - reject_upper: Whether upper bound H0 is rejected
        - conclusion: "bioequivalent" or "not_bioequivalent"

    Example:
        >>> result = tost_analysis(test_auc, reference_auc)
        >>> print(f"Conclusion: {result['conclusion']}")
        >>> print(f"p-values: lower={result['p_lower']:.4f}, upper={result['p_upper']:.4f}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    result = jl.NeoPKPD.tost_analysis(
        t, r,
        theta_lower=theta_lower,
        theta_upper=theta_upper,
        alpha=alpha
    )

    return {
        "parameter": str(result.parameter),
        "t_lower": float(result.t_lower),
        "t_upper": float(result.t_upper),
        "p_lower": float(result.p_lower),
        "p_upper": float(result.p_upper),
        "reject_lower": bool(result.reject_lower),
        "reject_upper": bool(result.reject_upper),
        "conclusion": str(result.be_conclusion),
    }


# ============================================================================
# BE Conclusion
# ============================================================================

def be_conclusion(
    ci_lower: float,
    ci_upper: float,
    theta_lower: float = 0.80,
    theta_upper: float = 1.25
) -> str:
    """
    Determine bioequivalence conclusion from confidence interval.

    Args:
        ci_lower: Lower bound of 90% CI
        ci_upper: Upper bound of 90% CI
        theta_lower: Lower equivalence bound (default: 0.80)
        theta_upper: Upper equivalence bound (default: 1.25)

    Returns:
        "bioequivalent", "not_bioequivalent", or "inconclusive"

    Example:
        >>> result = bioequivalence_90ci(test, reference)
        >>> conclusion = be_conclusion(result['ci_lower'], result['ci_upper'])
        >>> print(f"BE conclusion: {conclusion}")
    """
    jl = _require_julia()
    result = jl.NeoPKPD.be_conclusion(
        ci_lower, ci_upper,
        theta_lower=theta_lower,
        theta_upper=theta_upper
    )
    return str(result)


# ============================================================================
# Complete BE Analysis
# ============================================================================

@dataclass
class BioequivalenceResult:
    """
    Complete bioequivalence analysis result.

    Attributes:
        parameter: Parameter analyzed (e.g., "cmax", "auc")
        n_test: Number of test subjects
        n_reference: Number of reference subjects
        gmr: Geometric mean ratio
        ci_lower: Lower bound of 90% CI
        ci_upper: Upper bound of 90% CI
        cv_intra: Intra-subject CV (%)
        conclusion: BE conclusion
        be_limits: BE acceptance limits used
    """
    parameter: str
    n_test: int
    n_reference: int
    gmr: float
    ci_lower: float
    ci_upper: float
    cv_intra: float
    conclusion: str
    be_limits: Tuple[float, float]


def run_bioequivalence(
    parameter: str,
    test_values: List[float],
    reference_values: List[float],
    be_limits: Tuple[float, float] = (0.80, 1.25)
) -> BioequivalenceResult:
    """
    Run complete bioequivalence analysis for a parameter.

    Args:
        parameter: Parameter name (e.g., "cmax", "auc")
        test_values: Test formulation values
        reference_values: Reference formulation values
        be_limits: BE acceptance limits (default: (0.80, 1.25))

    Returns:
        BioequivalenceResult with complete analysis

    Example:
        >>> result = run_bioequivalence("cmax", test_cmax, ref_cmax)
        >>> print(f"Parameter: {result.parameter}")
        >>> print(f"GMR: {result.gmr:.4f} ({result.ci_lower:.4f}, {result.ci_upper:.4f})")
        >>> print(f"Conclusion: {result.conclusion}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    param_sym = jl.Symbol(parameter)

    result = jl.NeoPKPD.create_be_result(
        param_sym, t, r,
        be_limits=be_limits
    )

    return BioequivalenceResult(
        parameter=str(result.parameter),
        n_test=int(result.n_test),
        n_reference=int(result.n_reference),
        gmr=float(result.gmr),
        ci_lower=float(result.ci_lower),
        ci_upper=float(result.ci_upper),
        cv_intra=float(result.cv_intra),
        conclusion=str(result.be_conclusion),
        be_limits=be_limits,
    )


# ============================================================================
# Design-Aware 90% Confidence Interval
# ============================================================================

def bioequivalence_90ci_design(
    test_values: List[float],
    reference_values: List[float],
    design: BEStudyDesign,
    log_transform: bool = True,
    alpha: float = 0.10
) -> Dict[str, Any]:
    """
    Calculate 90% CI for GMR with proper degrees of freedom based on study design.

    This function provides regulatory-compliant CI calculation using the
    exact t-distribution with design-appropriate degrees of freedom.

    Args:
        test_values: Test formulation values (paired with reference)
        reference_values: Reference formulation values
        design: Study design type
        log_transform: Apply log transformation (default: True)
        alpha: Significance level (default: 0.10 for 90% CI)

    Returns:
        Dict with keys:
        - gmr: Geometric mean ratio
        - ci_lower: Lower bound of 90% CI
        - ci_upper: Upper bound of 90% CI
        - cv_intra: Intra-subject CV (%)
        - n: Number of subjects
        - df: Degrees of freedom used
        - t_critical: Critical t-value
        - design: Study design name

    Example:
        >>> result = bioequivalence_90ci_design(test, ref, BEStudyDesign.CROSSOVER_2X2)
        >>> print(f"GMR: {result['gmr']:.4f}")
        >>> print(f"df: {result['df']}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)

    # Create Julia design object
    if design == BEStudyDesign.CROSSOVER_2X2:
        jl_design = jl.NeoPKPD.Crossover2x2()
    elif design == BEStudyDesign.CROSSOVER_2X4:
        jl_design = jl.NeoPKPD.Crossover2x4()
    else:
        jl_design = jl.NeoPKPD.ParallelGroup()

    result = jl.NeoPKPD.bioequivalence_90ci_design(
        t, r, jl_design,
        log_transform=log_transform,
        alpha=alpha
    )

    return {
        "gmr": float(result.gmr),
        "ci_lower": float(result.ci_lower),
        "ci_upper": float(result.ci_upper),
        "cv_intra": float(result.cv_intra),
        "n": int(result.n),
        "df": int(result.df),
        "t_critical": float(result.t_critical),
        "design": str(result.design),
    }


def compute_be_degrees_of_freedom(
    design: BEStudyDesign,
    n: int
) -> int:
    """
    Compute correct degrees of freedom for bioequivalence CI based on study design.

    Args:
        design: Study design type
        n: Number of subjects

    Returns:
        Degrees of freedom for t-distribution
    """
    jl = _require_julia()

    # Create Julia design object
    if design == BEStudyDesign.CROSSOVER_2X2:
        jl_design = jl.NeoPKPD.Crossover2x2()
    elif design == BEStudyDesign.CROSSOVER_2X4:
        jl_design = jl.NeoPKPD.Crossover2x4()
    else:
        jl_design = jl.NeoPKPD.ParallelGroup()

    return int(jl.NeoPKPD.compute_be_degrees_of_freedom(jl_design, n))


# ============================================================================
# Reference-Scaled Average Bioequivalence (RSABE)
# ============================================================================

@dataclass
class RSABEConfig:
    """
    Configuration for reference-scaled average bioequivalence analysis.

    Attributes:
        guidance: Regulatory guidance (FDA, EMA, Health Canada)
        design: Replicate design type
        parameter: Parameter to analyze (e.g., "cmax", "auc_0_inf")
        alpha: Significance level (default: 0.05)
    """
    guidance: RegulatoryGuidance = RegulatoryGuidance.FDA
    design: ReplicateDesign = ReplicateDesign.FULL_REPLICATE_2X4
    parameter: str = "cmax"
    alpha: float = 0.05


@dataclass
class RSABEResult:
    """
    Result of FDA Reference-Scaled Average Bioequivalence analysis.

    Attributes:
        parameter: Parameter analyzed
        n_subjects: Number of subjects
        gmr: Geometric mean ratio (point estimate)
        log_diff: Mean log difference
        se_diff: Standard error of log difference
        swr: Within-subject reference SD
        swr_squared: Within-subject reference variance
        cv_wr: Within-subject reference CV (%)
        use_scaled: Whether scaled criterion was used
        rsabe_criterion: Value of RSABE criterion (if scaled)
        rsabe_upper_bound: Upper 95% CI for criterion
        point_estimate_pass: GMR within 80-125%
        scaled_criterion_pass: RSABE criterion satisfied
        be_conclusion: Overall BE conclusion
        scaled_limit_lower: Lower BE limit used
        scaled_limit_upper: Upper BE limit used
    """
    parameter: str
    n_subjects: int
    gmr: float
    log_diff: float
    se_diff: float
    swr: float
    swr_squared: float
    cv_wr: float
    use_scaled: bool
    rsabe_criterion: Optional[float]
    rsabe_upper_bound: Optional[float]
    point_estimate_pass: bool
    scaled_criterion_pass: bool
    be_conclusion: str
    scaled_limit_lower: float
    scaled_limit_upper: float


@dataclass
class ABELResult:
    """
    Result of EMA Average Bioequivalence with Expanding Limits analysis.

    Attributes:
        parameter: Parameter analyzed
        n_subjects: Number of subjects
        gmr: Geometric mean ratio
        ci_lower: 90% CI lower bound
        ci_upper: 90% CI upper bound
        swr: Within-subject reference SD
        cv_wr: Within-subject reference CV (%)
        use_scaled: Whether expanded limits were used
        be_limit_lower: Lower BE limit (possibly expanded)
        be_limit_upper: Upper BE limit (possibly expanded)
        be_conclusion: Overall BE conclusion
    """
    parameter: str
    n_subjects: int
    gmr: float
    ci_lower: float
    ci_upper: float
    swr: float
    cv_wr: float
    use_scaled: bool
    be_limit_lower: float
    be_limit_upper: float
    be_conclusion: str


def compute_swr(
    reference_values_1: List[float],
    reference_values_2: List[float]
) -> float:
    """
    Compute within-subject SD for reference (sigma_WR) from paired reference values.

    Args:
        reference_values_1: First reference measurement per subject
        reference_values_2: Second reference measurement per subject

    Returns:
        Within-subject reference SD (sigma_WR)
    """
    jl = _require_julia()
    r1 = _to_julia_float_vector(jl, reference_values_1)
    r2 = _to_julia_float_vector(jl, reference_values_2)
    return float(jl.NeoPKPD.compute_swr(r1, r2))


def compute_within_subject_variance_reference(
    reference_observations: List[List[float]]
) -> Dict[str, float]:
    """
    Compute within-subject variance for reference formulation from replicate data.

    This is the key estimate for reference-scaled bioequivalence.

    Args:
        reference_observations: List of reference observations per subject
            (each subject should have 2+ observations)

    Returns:
        Dict with keys:
        - swr_squared: Within-subject variance on log scale
        - swr: Within-subject standard deviation
        - cv_wr: Within-subject CV for reference (%)
        - df: Degrees of freedom
    """
    jl = _require_julia()

    # Convert to Julia Vector{Vector{Float64}}
    ref_obs_jl = jl.seval("Vector{Float64}[]")
    for obs in reference_observations:
        obs_vec = _to_julia_float_vector(jl, obs)
        jl.seval("push!")(ref_obs_jl, obs_vec)

    result = jl.NeoPKPD.compute_within_subject_variance_reference(ref_obs_jl)

    return {
        "swr_squared": float(result.swr_squared),
        "swr": float(result.swr),
        "cv_wr": float(result.cv_wr),
        "df": int(result.df),
    }


def rsabe_analysis(
    test_observations: List[List[float]],
    reference_observations: List[List[float]],
    config: Optional[RSABEConfig] = None
) -> RSABEResult:
    """
    Perform FDA Reference-Scaled Average Bioequivalence analysis.

    For highly variable drugs (CV > 30%), FDA allows scaling of BE limits
    based on within-subject reference variability.

    Args:
        test_observations: Test observations per subject (list of lists)
        reference_observations: Reference observations per subject
        config: RSABE configuration

    Returns:
        RSABEResult with complete analysis

    Example:
        >>> # Full replicate 2x4 design: each subject gets T, R, T, R
        >>> test_obs = [[100, 105], [95, 100], [110, 115], ...]  # per subject
        >>> ref_obs = [[98, 102], [92, 97], [108, 112], ...]
        >>> result = rsabe_analysis(test_obs, ref_obs)
        >>> print(f"Scaled: {result.use_scaled}")
        >>> print(f"Conclusion: {result.be_conclusion}")
    """
    jl = _require_julia()

    if config is None:
        config = RSABEConfig()

    # Convert observations to Julia Vector{Vector{Float64}}
    test_obs_jl = jl.seval("Vector{Float64}[]")
    for obs in test_observations:
        obs_vec = _to_julia_float_vector(jl, obs)
        jl.seval("push!")(test_obs_jl, obs_vec)

    ref_obs_jl = jl.seval("Vector{Float64}[]")
    for obs in reference_observations:
        obs_vec = _to_julia_float_vector(jl, obs)
        jl.seval("push!")(ref_obs_jl, obs_vec)

    # Create subject IDs and sequences (simple defaults)
    n_subjects = len(test_observations)
    subjects = [f"S{i+1}" for i in range(n_subjects)]
    sequences = [""] * n_subjects

    # Create ReplicateData
    subjects_jl = jl.seval("String[]")
    for s in subjects:
        jl.seval("push!")(subjects_jl, s)

    sequences_jl = jl.seval("String[]")
    for s in sequences:
        jl.seval("push!")(sequences_jl, s)

    n_test = len(test_observations[0]) if test_observations else 1
    n_ref = len(reference_observations[0]) if reference_observations else 2

    data = jl.NeoPKPD.ReplicateData(
        subjects_jl, sequences_jl, test_obs_jl, ref_obs_jl,
        n_subjects, n_test, n_ref
    )

    # Create guidance
    if config.guidance == RegulatoryGuidance.FDA:
        guidance = jl.NeoPKPD.FDAGuidance()
    elif config.guidance == RegulatoryGuidance.EMA:
        guidance = jl.NeoPKPD.EMAGuidance()
    else:
        guidance = jl.NeoPKPD.HealthCanadaGuidance()

    # Create design
    if config.design == ReplicateDesign.PARTIAL_REPLICATE_3X3:
        design = jl.NeoPKPD.PartialReplicate3x3()
    elif config.design == ReplicateDesign.FULL_REPLICATE_2X4:
        design = jl.NeoPKPD.FullReplicate2x4()
    else:
        design = jl.NeoPKPD.FullReplicate2x3()

    # Create config
    rsabe_config = jl.NeoPKPD.RSABEConfig(
        guidance=guidance,
        design=design,
        parameter=jl.Symbol(config.parameter),
        alpha=config.alpha
    )

    result = jl.NeoPKPD.rsabe_analysis(data, rsabe_config)

    return RSABEResult(
        parameter=str(result.parameter),
        n_subjects=int(result.n_subjects),
        gmr=float(result.gmr),
        log_diff=float(result.log_diff),
        se_diff=float(result.se_diff),
        swr=float(result.swr),
        swr_squared=float(result.swr_squared),
        cv_wr=float(result.cv_wr),
        use_scaled=bool(result.use_scaled),
        rsabe_criterion=_maybe_float(result.rsabe_criterion),
        rsabe_upper_bound=_maybe_float(result.rsabe_upper_bound),
        point_estimate_pass=bool(result.point_estimate_pass),
        scaled_criterion_pass=bool(result.scaled_criterion_pass),
        be_conclusion=str(result.be_conclusion),
        scaled_limit_lower=float(result.scaled_limit_lower),
        scaled_limit_upper=float(result.scaled_limit_upper),
    )


def abel_analysis(
    test_observations: List[List[float]],
    reference_observations: List[List[float]],
    config: Optional[RSABEConfig] = None
) -> ABELResult:
    """
    Perform EMA Average Bioequivalence with Expanding Limits analysis.

    For highly variable drugs, EMA allows expansion of BE limits based on
    within-subject reference variability.

    Args:
        test_observations: Test observations per subject
        reference_observations: Reference observations per subject
        config: Configuration (guidance will be set to EMA)

    Returns:
        ABELResult with complete analysis

    Example:
        >>> result = abel_analysis(test_obs, ref_obs)
        >>> print(f"Expanded limits: {result.be_limit_lower:.4f} - {result.be_limit_upper:.4f}")
    """
    jl = _require_julia()

    if config is None:
        config = RSABEConfig(guidance=RegulatoryGuidance.EMA)
    else:
        config = RSABEConfig(
            guidance=RegulatoryGuidance.EMA,
            design=config.design,
            parameter=config.parameter,
            alpha=config.alpha
        )

    # Convert observations to Julia
    test_obs_jl = jl.seval("Vector{Float64}[]")
    for obs in test_observations:
        obs_vec = _to_julia_float_vector(jl, obs)
        jl.seval("push!")(test_obs_jl, obs_vec)

    ref_obs_jl = jl.seval("Vector{Float64}[]")
    for obs in reference_observations:
        obs_vec = _to_julia_float_vector(jl, obs)
        jl.seval("push!")(ref_obs_jl, obs_vec)

    n_subjects = len(test_observations)
    subjects_jl = jl.seval("String[]")
    for i in range(n_subjects):
        jl.seval("push!")(subjects_jl, f"S{i+1}")

    sequences_jl = jl.seval("String[]")
    for _ in range(n_subjects):
        jl.seval("push!")(sequences_jl, "")

    n_test = len(test_observations[0]) if test_observations else 1
    n_ref = len(reference_observations[0]) if reference_observations else 2

    data = jl.NeoPKPD.ReplicateData(
        subjects_jl, sequences_jl, test_obs_jl, ref_obs_jl,
        n_subjects, n_test, n_ref
    )

    if config.design == ReplicateDesign.PARTIAL_REPLICATE_3X3:
        design = jl.NeoPKPD.PartialReplicate3x3()
    elif config.design == ReplicateDesign.FULL_REPLICATE_2X4:
        design = jl.NeoPKPD.FullReplicate2x4()
    else:
        design = jl.NeoPKPD.FullReplicate2x3()

    rsabe_config = jl.NeoPKPD.RSABEConfig(
        guidance=jl.NeoPKPD.EMAGuidance(),
        design=design,
        parameter=jl.Symbol(config.parameter),
        alpha=config.alpha
    )

    result = jl.NeoPKPD.abel_analysis(data, rsabe_config)

    return ABELResult(
        parameter=str(result.parameter),
        n_subjects=int(result.n_subjects),
        gmr=float(result.gmr),
        ci_lower=float(result.ci_lower),
        ci_upper=float(result.ci_upper),
        swr=float(result.swr),
        cv_wr=float(result.cv_wr),
        use_scaled=bool(result.use_scaled),
        be_limit_lower=float(result.be_limit_lower),
        be_limit_upper=float(result.be_limit_upper),
        be_conclusion=str(result.be_conclusion),
    )


def reference_scaled_be(
    test_observations: List[List[float]],
    reference_observations: List[List[float]],
    config: Optional[RSABEConfig] = None
) -> Union[RSABEResult, ABELResult]:
    """
    Unified entry point for reference-scaled bioequivalence analysis.

    Automatically dispatches to FDA RSABE or EMA ABEL based on guidance.

    Args:
        test_observations: Test observations per subject
        reference_observations: Reference observations per subject
        config: Configuration specifying guidance and design

    Returns:
        RSABEResult (FDA/Health Canada) or ABELResult (EMA)
    """
    if config is None:
        config = RSABEConfig()

    if config.guidance == RegulatoryGuidance.EMA:
        return abel_analysis(test_observations, reference_observations, config)
    else:
        return rsabe_analysis(test_observations, reference_observations, config)


def abel_scaled_limits(swr: float) -> Dict[str, Any]:
    """
    Calculate EMA ABEL expanded bioequivalence limits.

    When CV > 30%, EMA allows expanded limits based on sigma_WR.

    Args:
        swr: Within-subject reference SD

    Returns:
        Dict with keys:
        - lower: Lower BE limit
        - upper: Upper BE limit
        - is_scaled: Whether scaling was applied
    """
    jl = _require_julia()
    guidance = jl.NeoPKPD.EMAGuidance()
    result = jl.NeoPKPD.abel_scaled_limits(swr, guidance)

    return {
        "lower": float(result.lower),
        "upper": float(result.upper),
        "is_scaled": bool(result.is_scaled),
    }


# ============================================================================
# Helper
# ============================================================================

def _maybe_float(val: Any) -> Optional[float]:
    """Convert value to float, returning None for NaN or None."""
    if val is None:
        return None
    try:
        f = float(val)
        import math
        if math.isnan(f):
            return None
        return f
    except (TypeError, ValueError):
        return None
