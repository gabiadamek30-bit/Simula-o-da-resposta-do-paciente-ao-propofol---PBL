"""
NeoPKPD Trial Analysis

Statistical analysis for clinical trials including:
- Power analysis
- Sample size estimation
- Alpha spending functions
- Arm comparisons
- Responder analysis
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import math


@dataclass
class PowerResult:
    """
    Power analysis result.

    Attributes:
        power: Estimated power
        alpha: Significance level
        effect_size: Effect size
        n_per_arm: Sample size per arm
        n_simulations: Number of simulations (if simulation-based)
        method: Analysis method ('analytical' or 'simulation')
    """
    power: float
    alpha: float
    effect_size: float
    n_per_arm: int
    n_simulations: int = 0
    method: str = "analytical"


@dataclass
class SampleSizeResult:
    """
    Sample size estimation result.

    Attributes:
        n_per_arm: Required sample size per arm
        total_n: Total required sample size
        target_power: Target power
        achieved_power: Achieved power at this sample size
        alpha: Significance level
        effect_size: Effect size
    """
    n_per_arm: int
    total_n: int
    target_power: float
    achieved_power: float
    alpha: float
    effect_size: float


@dataclass
class ComparisonResult:
    """
    Arm comparison result.

    Attributes:
        arm1: First arm name
        arm2: Second arm name
        difference: Mean difference
        ci_lower: Lower confidence interval
        ci_upper: Upper confidence interval
        p_value: P-value
        significant: Whether difference is significant
    """
    arm1: str
    arm2: str
    difference: float
    ci_lower: float
    ci_upper: float
    p_value: float
    significant: bool


@dataclass
class ResponderResult:
    """
    Responder analysis result.

    Attributes:
        arm: Arm name
        n_total: Total subjects
        n_responders: Number of responders
        response_rate: Response rate
        ci_lower: Lower confidence interval
        ci_upper: Upper confidence interval
    """
    arm: str
    n_total: int
    n_responders: int
    response_rate: float
    ci_lower: float
    ci_upper: float


def _norm_cdf(x: float) -> float:
    """Standard normal CDF approximation."""
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def _norm_ppf(p: float) -> float:
    """Standard normal quantile function (inverse CDF) approximation."""
    if p <= 0:
        return float('-inf')
    if p >= 1:
        return float('inf')

    # Rational approximation for normal inverse CDF
    # Abramowitz and Stegun approximation
    a = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ]
    b = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ]
    c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ]
    d = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ]

    p_low = 0.02425
    p_high = 1 - p_low

    if p < p_low:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    elif p <= p_high:
        q = p - 0.5
        r = q * q
        return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q / \
               (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)
    else:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)


def estimate_power_analytical(
    n_per_arm: int,
    effect_size: float,
    sd: float = 1.0,
    alpha: float = 0.05,
    alternative: str = "two-sided",
) -> PowerResult:
    """
    Estimate power analytically for a two-sample t-test.

    Args:
        n_per_arm: Sample size per arm
        effect_size: Standardized effect size (Cohen's d)
        sd: Standard deviation (default: 1.0 for standardized)
        alpha: Significance level
        alternative: 'two-sided', 'greater', or 'less'

    Returns:
        PowerResult object

    Example:
        >>> result = estimate_power_analytical(n_per_arm=50, effect_size=0.5)
        >>> print(f"Power: {result.power:.2%}")
    """
    # Calculate non-centrality parameter
    ncp = effect_size * math.sqrt(n_per_arm / 2)

    # Calculate critical value
    if alternative == "two-sided":
        z_alpha = _norm_ppf(1 - alpha / 2)
        # Power = P(|Z| > z_alpha | H1)
        power = 1 - _norm_cdf(z_alpha - ncp) + _norm_cdf(-z_alpha - ncp)
    elif alternative == "greater":
        z_alpha = _norm_ppf(1 - alpha)
        power = 1 - _norm_cdf(z_alpha - ncp)
    else:  # less
        z_alpha = _norm_ppf(alpha)
        power = _norm_cdf(z_alpha - ncp)

    return PowerResult(
        power=power,
        alpha=alpha,
        effect_size=effect_size,
        n_per_arm=n_per_arm,
        method="analytical",
    )


def estimate_sample_size(
    target_power: float,
    effect_size: float,
    sd: float = 1.0,
    alpha: float = 0.05,
    alternative: str = "two-sided",
    n_arms: int = 2,
) -> SampleSizeResult:
    """
    Estimate required sample size for target power.

    Args:
        target_power: Target power (e.g., 0.80)
        effect_size: Standardized effect size
        sd: Standard deviation
        alpha: Significance level
        alternative: 'two-sided', 'greater', or 'less'
        n_arms: Number of arms (default: 2)

    Returns:
        SampleSizeResult object

    Example:
        >>> result = estimate_sample_size(target_power=0.80, effect_size=0.5)
        >>> print(f"Required n per arm: {result.n_per_arm}")
    """
    # Binary search for sample size
    n_low = 2
    n_high = 10000

    while n_high - n_low > 1:
        n_mid = (n_low + n_high) // 2
        power_result = estimate_power_analytical(
            n_per_arm=n_mid,
            effect_size=effect_size,
            sd=sd,
            alpha=alpha,
            alternative=alternative,
        )
        if power_result.power < target_power:
            n_low = n_mid
        else:
            n_high = n_mid

    # Use the higher value to ensure target power is met
    n_per_arm = n_high
    achieved = estimate_power_analytical(
        n_per_arm=n_per_arm,
        effect_size=effect_size,
        sd=sd,
        alpha=alpha,
        alternative=alternative,
    )

    return SampleSizeResult(
        n_per_arm=n_per_arm,
        total_n=n_per_arm * n_arms,
        target_power=target_power,
        achieved_power=achieved.power,
        alpha=alpha,
        effect_size=effect_size,
    )


def alpha_spending_function(
    information_fraction: float,
    total_alpha: float = 0.05,
    spending_type: str = "obrien_fleming",
) -> float:
    """
    Calculate cumulative alpha spent at a given information fraction.

    Args:
        information_fraction: Information fraction (0 to 1)
        total_alpha: Total alpha to spend
        spending_type: 'obrien_fleming', 'pocock', or 'haybittle_peto'

    Returns:
        Cumulative alpha spent

    Example:
        >>> # At 50% information
        >>> alpha = alpha_spending_function(0.5, 0.05, 'obrien_fleming')
    """
    t = information_fraction

    if spending_type == "obrien_fleming":
        # O'Brien-Fleming spending function
        if t <= 0:
            return 0.0
        if t >= 1:
            return total_alpha
        z = _norm_ppf(1 - total_alpha / 2)
        return 2 * (1 - _norm_cdf(z / math.sqrt(t)))

    elif spending_type == "pocock":
        # Pocock spending function
        return total_alpha * math.log(1 + (math.e - 1) * t)

    elif spending_type == "haybittle_peto":
        # Haybittle-Peto: fixed boundary until final
        if t < 1:
            return 0.001  # Use very small alpha for interim
        return total_alpha

    else:
        # Linear spending
        return total_alpha * t


def incremental_alpha(
    information_fractions: List[float],
    total_alpha: float = 0.05,
    spending_type: str = "obrien_fleming",
) -> List[float]:
    """
    Calculate incremental alpha to spend at each analysis.

    Args:
        information_fractions: List of information fractions
        total_alpha: Total alpha to spend
        spending_type: Alpha spending function type

    Returns:
        List of incremental alphas for each analysis

    Example:
        >>> # Two interim analyses at 50% and 75%, final at 100%
        >>> alphas = incremental_alpha([0.5, 0.75, 1.0])
    """
    cumulative = [
        alpha_spending_function(t, total_alpha, spending_type)
        for t in information_fractions
    ]

    incremental = [cumulative[0]]
    for i in range(1, len(cumulative)):
        incremental.append(cumulative[i] - cumulative[i-1])

    return incremental


def compare_arms(
    values1: List[float],
    values2: List[float],
    arm1_name: str = "Arm 1",
    arm2_name: str = "Arm 2",
    alpha: float = 0.05,
) -> ComparisonResult:
    """
    Compare two treatment arms using a two-sample t-test approximation.

    Args:
        values1: Values from first arm
        values2: Values from second arm
        arm1_name: First arm name
        arm2_name: Second arm name
        alpha: Significance level

    Returns:
        ComparisonResult object

    Example:
        >>> arm1_values = [10.2, 11.5, 9.8, 12.1, 10.5]
        >>> arm2_values = [8.1, 7.9, 8.5, 7.2, 8.8]
        >>> result = compare_arms(arm1_values, arm2_values, "Active", "Placebo")
    """
    n1 = len(values1)
    n2 = len(values2)

    if n1 < 2 or n2 < 2:
        return ComparisonResult(
            arm1=arm1_name,
            arm2=arm2_name,
            difference=0.0,
            ci_lower=0.0,
            ci_upper=0.0,
            p_value=1.0,
            significant=False,
        )

    mean1 = sum(values1) / n1
    mean2 = sum(values2) / n2
    diff = mean1 - mean2

    var1 = sum((x - mean1) ** 2 for x in values1) / (n1 - 1)
    var2 = sum((x - mean2) ** 2 for x in values2) / (n2 - 1)

    # Pooled standard error
    se = math.sqrt(var1 / n1 + var2 / n2)

    if se == 0:
        return ComparisonResult(
            arm1=arm1_name,
            arm2=arm2_name,
            difference=diff,
            ci_lower=diff,
            ci_upper=diff,
            p_value=1.0 if diff == 0 else 0.0,
            significant=diff != 0,
        )

    # Z-score (using normal approximation for large n)
    z = diff / se

    # Two-sided p-value
    p_value = 2 * (1 - _norm_cdf(abs(z)))

    # Confidence interval
    z_alpha = _norm_ppf(1 - alpha / 2)
    ci_lower = diff - z_alpha * se
    ci_upper = diff + z_alpha * se

    return ComparisonResult(
        arm1=arm1_name,
        arm2=arm2_name,
        difference=diff,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        p_value=p_value,
        significant=p_value < alpha,
    )


def responder_analysis(
    values: List[float],
    threshold: float,
    arm_name: str = "Arm",
    direction: str = "greater",
    confidence: float = 0.95,
) -> ResponderResult:
    """
    Perform responder analysis.

    Args:
        values: Endpoint values
        threshold: Response threshold
        arm_name: Arm name
        direction: 'greater' or 'less' (response if value > or < threshold)
        confidence: Confidence level for interval

    Returns:
        ResponderResult object

    Example:
        >>> values = [10.2, 15.5, 8.1, 12.3, 9.5, 14.2, 11.8]
        >>> result = responder_analysis(values, threshold=10.0, direction='greater')
        >>> print(f"Response rate: {result.response_rate:.1%}")
    """
    n_total = len(values)

    if direction == "greater":
        n_responders = sum(1 for v in values if v > threshold)
    else:
        n_responders = sum(1 for v in values if v < threshold)

    response_rate = n_responders / n_total if n_total > 0 else 0.0

    # Wilson score interval for proportion
    alpha = 1 - confidence
    z = _norm_ppf(1 - alpha / 2)

    if n_total == 0:
        ci_lower = 0.0
        ci_upper = 1.0
    else:
        p = response_rate
        denominator = 1 + z**2 / n_total
        center = (p + z**2 / (2 * n_total)) / denominator
        margin = z * math.sqrt((p * (1 - p) + z**2 / (4 * n_total)) / n_total) / denominator

        ci_lower = max(0.0, center - margin)
        ci_upper = min(1.0, center + margin)

    return ResponderResult(
        arm=arm_name,
        n_total=n_total,
        n_responders=n_responders,
        response_rate=response_rate,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
    )


def bioequivalence_90ci(
    test_values: List[float],
    reference_values: List[float],
    log_transform: bool = True,
) -> Tuple[float, float]:
    """
    Calculate 90% confidence interval for bioequivalence assessment.

    Args:
        test_values: Test formulation values
        reference_values: Reference formulation values
        log_transform: Whether to log-transform values

    Returns:
        Tuple of (lower bound, upper bound) of 90% CI for ratio

    Example:
        >>> test = [95.2, 102.1, 98.5, 105.3, 97.8]
        >>> ref = [100.0, 98.5, 101.2, 99.8, 100.5]
        >>> lower, upper = bioequivalence_90ci(test, ref)
        >>> if 0.80 <= lower and upper <= 1.25:
        ...     print("Bioequivalent")
    """
    if log_transform:
        test_values = [math.log(v) for v in test_values if v > 0]
        reference_values = [math.log(v) for v in reference_values if v > 0]

    n_test = len(test_values)
    n_ref = len(reference_values)

    if n_test < 2 or n_ref < 2:
        return (0.0, float('inf'))

    mean_test = sum(test_values) / n_test
    mean_ref = sum(reference_values) / n_ref
    diff = mean_test - mean_ref

    var_test = sum((x - mean_test) ** 2 for x in test_values) / (n_test - 1)
    var_ref = sum((x - mean_ref) ** 2 for x in reference_values) / (n_ref - 1)

    se = math.sqrt(var_test / n_test + var_ref / n_ref)

    # 90% CI (two one-sided tests at alpha=0.05)
    t_alpha = _norm_ppf(0.95)  # Normal approximation
    ci_lower = diff - t_alpha * se
    ci_upper = diff + t_alpha * se

    if log_transform:
        return (math.exp(ci_lower), math.exp(ci_upper))
    else:
        return (ci_lower, ci_upper)


def assess_bioequivalence(
    test_values: List[float],
    reference_values: List[float],
    lower_limit: float = 0.80,
    upper_limit: float = 1.25,
    log_transform: bool = True,
) -> Dict[str, Any]:
    """
    Perform complete bioequivalence assessment.

    Args:
        test_values: Test formulation values
        reference_values: Reference formulation values
        lower_limit: Lower BE limit (default: 0.80)
        upper_limit: Upper BE limit (default: 1.25)
        log_transform: Whether to log-transform values

    Returns:
        Dictionary with BE assessment results

    Example:
        >>> result = assess_bioequivalence(test, reference)
        >>> print(f"BE conclusion: {'Pass' if result['bioequivalent'] else 'Fail'}")
    """
    ci_lower, ci_upper = bioequivalence_90ci(test_values, reference_values, log_transform)

    # Calculate point estimate
    if log_transform:
        log_test = [math.log(v) for v in test_values if v > 0]
        log_ref = [math.log(v) for v in reference_values if v > 0]
        point_estimate = math.exp(sum(log_test) / len(log_test) - sum(log_ref) / len(log_ref))
    else:
        point_estimate = sum(test_values) / len(test_values) / (sum(reference_values) / len(reference_values))

    bioequivalent = ci_lower >= lower_limit and ci_upper <= upper_limit

    return {
        "point_estimate": point_estimate,
        "ci_90_lower": ci_lower,
        "ci_90_upper": ci_upper,
        "be_lower_limit": lower_limit,
        "be_upper_limit": upper_limit,
        "bioequivalent": bioequivalent,
        "n_test": len(test_values),
        "n_reference": len(reference_values),
    }
