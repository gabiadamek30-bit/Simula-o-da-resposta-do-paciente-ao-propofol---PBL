"""
NeoPKPD NCA Metrics - Individual NCA metric functions.

This module provides Python bindings to the Julia NCA metric functions.
All functions call the corresponding Julia implementations for FDA/EMA
compliant calculations.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from .._core import _require_julia, _to_julia_float_vector


# ============================================================================
# Primary Exposure Metrics
# ============================================================================

def nca_cmax(concentrations: List[float]) -> float:
    """
    Find maximum observed concentration (Cmax).

    Args:
        concentrations: List of concentration values

    Returns:
        Maximum concentration

    Example:
        >>> cmax = nca_cmax([0.0, 1.2, 2.0, 1.5, 0.8])
        >>> print(f"Cmax: {cmax}")
    """
    jl = _require_julia()
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.NeoPKPD.nca_cmax(c))


def nca_tmax(times: List[float], concentrations: List[float]) -> float:
    """
    Find time of maximum observed concentration (Tmax).

    Args:
        times: List of time points
        concentrations: List of concentration values

    Returns:
        Time of maximum concentration

    Example:
        >>> tmax = nca_tmax([0.0, 1.0, 2.0, 4.0], [0.0, 1.5, 2.0, 1.0])
        >>> print(f"Tmax: {tmax}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.NeoPKPD.nca_tmax(t, c))


def nca_cmin(concentrations: List[float]) -> float:
    """
    Find minimum observed concentration (Cmin).

    Args:
        concentrations: List of concentration values

    Returns:
        Minimum concentration

    Example:
        >>> cmin = nca_cmin([2.5, 4.2, 3.8, 2.5])
        >>> print(f"Cmin: {cmin}")
    """
    jl = _require_julia()
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.NeoPKPD.nca_cmin(c))


def nca_clast(
    times: List[float],
    concentrations: List[float],
    lloq: float = 0.0
) -> float:
    """
    Find last measurable concentration (Clast).

    Args:
        times: List of time points
        concentrations: List of concentration values
        lloq: Lower limit of quantification

    Returns:
        Last measurable concentration
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.NeoPKPD.nca_clast(t, c, lloq=lloq))


def nca_tlast(
    times: List[float],
    concentrations: List[float],
    lloq: float = 0.0
) -> float:
    """
    Find time of last measurable concentration (Tlast).

    Args:
        times: List of time points
        concentrations: List of concentration values
        lloq: Lower limit of quantification

    Returns:
        Time of last measurable concentration
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.NeoPKPD.nca_tlast(t, c, lloq=lloq))


def nca_cavg(
    times: List[float],
    concentrations: List[float],
    tau: float,
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate average concentration over a dosing interval (Cavg).

    Cavg = AUC0-tau / tau

    Args:
        times: List of time points
        concentrations: List of concentration values
        tau: Dosing interval
        method: AUC calculation method ("linear", "log_linear", "lin_log_mixed")

    Returns:
        Average concentration over dosing interval
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.NeoPKPD.nca_cavg(t, c, tau, config))


# ============================================================================
# Lambda-z Estimation
# ============================================================================

def estimate_lambda_z(
    times: List[float],
    concentrations: List[float],
    min_points: int = 3,
    r2_threshold: float = 0.9,
    method: str = "lin_log_mixed"
) -> Dict[str, Any]:
    """
    Estimate terminal elimination rate constant (lambda_z).

    Uses log-linear regression on the terminal phase to estimate lambda_z.

    Args:
        times: List of time points
        concentrations: List of concentration values
        min_points: Minimum points for regression (default: 3)
        r2_threshold: Minimum R-squared threshold (default: 0.9)
        method: AUC calculation method

    Returns:
        Dict with keys:
        - lambda_z: Terminal elimination rate constant (or None)
        - t_half: Terminal half-life (or None)
        - r_squared: R-squared of regression (or None)
        - n_points: Number of points used
        - quality_flag: Quality assessment ("good", "warning", "insufficient")
        - warnings: List of quality warnings

    Example:
        >>> result = estimate_lambda_z([0, 1, 2, 4, 8, 12, 24], [0, 2, 1.8, 1.2, 0.6, 0.3, 0.075])
        >>> print(f"Lambda_z: {result['lambda_z']}")
        >>> print(f"t1/2: {result['t_half']}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method, lambda_z_min_points=min_points,
                              lambda_z_r2_threshold=r2_threshold)

    result = jl.NeoPKPD.estimate_lambda_z(t, c, config)

    return {
        "lambda_z": _maybe_float(result.lambda_z),
        "t_half": _maybe_float(result.t_half),
        "r_squared": _maybe_float(result.r_squared),
        "adjusted_r_squared": _maybe_float(result.adjusted_r_squared),
        "intercept": _maybe_float(result.intercept),
        "n_points": int(result.n_points),
        "start_time": float(result.start_time),
        "end_time": float(result.end_time),
        "points_used": list(result.points_used),
        "quality_flag": str(result.quality_flag),
        "warnings": list(result.warnings),
    }


def nca_half_life(lambda_z: float) -> float:
    """
    Calculate terminal half-life from lambda_z.

    t1/2 = ln(2) / lambda_z

    Args:
        lambda_z: Terminal elimination rate constant

    Returns:
        Terminal half-life

    Example:
        >>> t_half = nca_half_life(0.1)
        >>> print(f"Half-life: {t_half:.2f} hours")
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_half_life(lambda_z))


# ============================================================================
# AUC Calculations
# ============================================================================

def auc_0_t(
    times: List[float],
    concentrations: List[float],
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate AUC from time 0 to last measurable concentration (AUC0-t).

    Args:
        times: List of time points
        concentrations: List of concentration values
        method: Calculation method ("linear", "log_linear", "lin_log_mixed")

    Returns:
        AUC from 0 to last measurable concentration

    Example:
        >>> auc = auc_0_t([0, 1, 2, 4, 8], [0, 2, 1.5, 1.0, 0.5])
        >>> print(f"AUC0-t: {auc:.2f}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.NeoPKPD.auc_0_t(t, c, config))


def auc_0_inf(
    times: List[float],
    concentrations: List[float],
    lambda_z: float,
    clast: float,
    method: str = "lin_log_mixed"
) -> Tuple[float, float]:
    """
    Calculate AUC extrapolated to infinity (AUC0-inf).

    AUC0-inf = AUC0-t + Clast/lambda_z

    Args:
        times: List of time points
        concentrations: List of concentration values
        lambda_z: Terminal elimination rate constant
        clast: Last measurable concentration
        method: Calculation method

    Returns:
        Tuple of (AUC0-inf, extrapolation percentage)

    Example:
        >>> auc_inf, extra_pct = auc_0_inf(t, c, lambda_z=0.1, clast=0.5)
        >>> print(f"AUC0-inf: {auc_inf:.2f}, Extrapolation: {extra_pct:.1f}%")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    result = jl.NeoPKPD.auc_0_inf(t, c, lambda_z, clast, config)
    return (float(result[0]), float(result[1]))


def auc_0_tau(
    times: List[float],
    concentrations: List[float],
    tau: float,
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate AUC over a dosing interval (AUC0-tau).

    Args:
        times: List of time points
        concentrations: List of concentration values
        tau: Dosing interval
        method: Calculation method

    Returns:
        AUC over the dosing interval
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.NeoPKPD.auc_0_tau(t, c, tau, config))


def auc_partial(
    times: List[float],
    concentrations: List[float],
    t_start: float,
    t_end: float,
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate partial AUC between two time points.

    Args:
        times: List of time points
        concentrations: List of concentration values
        t_start: Start time
        t_end: End time
        method: Calculation method

    Returns:
        Partial AUC
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.NeoPKPD.auc_partial(t, c, t_start, t_end, config))


def aumc_0_t(
    times: List[float],
    concentrations: List[float],
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate AUMC from time 0 to last measurable concentration.

    AUMC = Area Under the (first) Moment Curve = integral of t*C(t)

    Args:
        times: List of time points
        concentrations: List of concentration values
        method: Calculation method

    Returns:
        AUMC from 0 to last measurable concentration
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.NeoPKPD.aumc_0_t(t, c, config))


def aumc_0_inf(
    times: List[float],
    concentrations: List[float],
    lambda_z: float,
    clast: float,
    tlast: float,
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate AUMC extrapolated to infinity.

    Args:
        times: List of time points
        concentrations: List of concentration values
        lambda_z: Terminal elimination rate constant
        clast: Last measurable concentration
        tlast: Time of last measurable concentration
        method: Calculation method

    Returns:
        AUMC extrapolated to infinity
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.NeoPKPD.aumc_0_inf(t, c, lambda_z, clast, tlast, config))


# ============================================================================
# PK Parameters
# ============================================================================

def nca_mrt(
    aumc_0_inf: float,
    auc_0_inf: float,
    route: str = "extravascular",
    t_inf: float = 0.0
) -> float:
    """
    Calculate Mean Residence Time (MRT).

    For extravascular: MRT = AUMC0-inf / AUC0-inf
    For IV infusion: MRT = AUMC0-inf / AUC0-inf - Tinf/2

    Args:
        aumc_0_inf: AUMC extrapolated to infinity
        auc_0_inf: AUC extrapolated to infinity
        route: Administration route ("extravascular", "iv_bolus", "iv_infusion")
        t_inf: Infusion duration (for iv_infusion)

    Returns:
        Mean residence time
    """
    jl = _require_julia()
    route_sym = jl.Symbol(route)
    return float(jl.NeoPKPD.nca_mrt(aumc_0_inf, auc_0_inf, route=route_sym, t_inf=t_inf))


def nca_cl_f(dose: float, auc_0_inf: float) -> float:
    """
    Calculate apparent clearance (CL/F).

    CL/F = Dose / AUC0-inf

    Args:
        dose: Administered dose
        auc_0_inf: AUC extrapolated to infinity

    Returns:
        Apparent clearance
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_cl_f(dose, auc_0_inf))


def nca_vz_f(dose: float, lambda_z: float, auc_0_inf: float) -> float:
    """
    Calculate apparent volume of distribution (Vz/F).

    Vz/F = Dose / (lambda_z * AUC0-inf)

    Args:
        dose: Administered dose
        lambda_z: Terminal elimination rate constant
        auc_0_inf: AUC extrapolated to infinity

    Returns:
        Apparent volume of distribution
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_vz_f(dose, lambda_z, auc_0_inf))


def nca_vss(cl: float, mrt: float) -> float:
    """
    Calculate volume of distribution at steady state (Vss).

    Vss = CL * MRT

    Args:
        cl: Clearance
        mrt: Mean residence time

    Returns:
        Volume at steady state
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_vss(cl, mrt))


def nca_bioavailability(
    auc_test: float,
    dose_test: float,
    auc_reference: float,
    dose_reference: float
) -> float:
    """
    Calculate relative bioavailability (F).

    F = (AUC_test / Dose_test) / (AUC_ref / Dose_ref)

    Args:
        auc_test: AUC of test formulation
        dose_test: Dose of test formulation
        auc_reference: AUC of reference formulation
        dose_reference: Dose of reference formulation

    Returns:
        Relative bioavailability
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_bioavailability(auc_test, dose_test, auc_reference, dose_reference))


# ============================================================================
# Additional Exposure Metrics
# ============================================================================

def nca_ctrough(
    times: List[float],
    concentrations: List[float],
    tau: float
) -> float:
    """
    Find trough concentration at end of dosing interval.

    Args:
        times: List of time points
        concentrations: List of concentration values
        tau: Dosing interval

    Returns:
        Trough concentration
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.NeoPKPD.nca_ctrough(t, c, tau))


def nca_c_at_time(
    times: List[float],
    concentrations: List[float],
    target_time: float
) -> float:
    """
    Get concentration at a specific time point (with interpolation if needed).

    Args:
        times: List of time points
        concentrations: List of concentration values
        target_time: Target time for concentration

    Returns:
        Concentration at target time
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.NeoPKPD.nca_c_at_time(t, c, target_time))


def time_above_concentration(
    times: List[float],
    concentrations: List[float],
    threshold: float
) -> float:
    """
    Calculate time above a concentration threshold.

    Args:
        times: List of time points
        concentrations: List of concentration values
        threshold: Concentration threshold

    Returns:
        Total time where concentration > threshold
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.NeoPKPD.time_above_concentration(t, c, threshold))


# ============================================================================
# Additional PK Parameters
# ============================================================================

def nca_cl(dose: float, auc_0_inf: float) -> float:
    """
    Calculate systemic clearance (CL) for IV administration.

    CL = Dose / AUC0-inf

    Args:
        dose: Administered dose
        auc_0_inf: AUC extrapolated to infinity

    Returns:
        Systemic clearance
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_cl(dose, auc_0_inf))


def nca_vz(dose: float, lambda_z: float, auc_0_inf: float) -> float:
    """
    Calculate terminal volume of distribution (Vz) for IV administration.

    Vz = Dose / (lambda_z * AUC0-inf)

    Args:
        dose: Administered dose
        lambda_z: Terminal elimination rate constant
        auc_0_inf: AUC extrapolated to infinity

    Returns:
        Terminal volume of distribution
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_vz(dose, lambda_z, auc_0_inf))


def nca_mrt_iv(mrt_extravascular: float, absorption_time: float) -> float:
    """
    Estimate MRT for IV administration from extravascular MRT.

    MRTiv = MRText - MAT

    Args:
        mrt_extravascular: MRT from extravascular administration
        absorption_time: Mean absorption time estimate

    Returns:
        Estimated IV MRT
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_mrt_iv(mrt_extravascular, absorption_time))


def nca_cl_ss(dose: float, auc_0_tau: float) -> float:
    """
    Calculate clearance at steady state.

    CLss = Dose / AUC0-tau

    Args:
        dose: Administered dose per interval
        auc_0_tau: AUC over dosing interval

    Returns:
        Steady-state clearance
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_cl_ss(dose, auc_0_tau))


def nca_vss_from_aumc(dose: float, auc_0_inf: float, aumc_0_inf: float) -> float:
    """
    Calculate Vss directly from moment curves.

    Vss = Dose * AUMC0-inf / (AUC0-inf)^2

    Args:
        dose: Administered dose
        auc_0_inf: AUC from 0 to infinity
        aumc_0_inf: AUMC from 0 to infinity

    Returns:
        Volume at steady state
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_vss_from_aumc(dose, auc_0_inf, aumc_0_inf))


def nca_vc(dose: float, c0: float) -> float:
    """
    Calculate central volume of distribution from back-extrapolated C0.

    Vc = Dose / C0

    For IV bolus administration.

    Args:
        dose: Administered dose
        c0: Back-extrapolated concentration at time 0

    Returns:
        Central volume of distribution
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_vc(dose, c0))


def nca_mean_absorption_time(mrt_po: float, mrt_iv: float) -> float:
    """
    Calculate Mean Absorption Time (MAT).

    MAT = MRTpo - MRTiv

    Args:
        mrt_po: MRT from oral/extravascular administration
        mrt_iv: MRT from IV administration

    Returns:
        Mean absorption time
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_mean_absorption_time(mrt_po, mrt_iv))


def nca_c0_backextrap(
    times: List[float],
    concentrations: List[float],
    lambda_z: float
) -> float:
    """
    Back-extrapolate C0 for IV bolus administration.

    Uses the lambda_z regression parameters.

    Args:
        times: Time points (from terminal phase regression)
        concentrations: Concentration values
        lambda_z: Terminal elimination rate constant

    Returns:
        Back-extrapolated C0
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.NeoPKPD.nca_c0_backextrap(t, c, lambda_z))


def nca_c0_from_regression(intercept: float) -> float:
    """
    Get C0 from lambda_z regression intercept.

    For log-linear regression: ln(C) = intercept - lambda_z * t
    C0 = exp(intercept)

    Args:
        intercept: Y-intercept from log-linear regression

    Returns:
        Back-extrapolated C0
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_c0_from_regression(intercept))


# ============================================================================
# Multiple Dose Metrics
# ============================================================================

def nca_accumulation_index(auc_ss: float, auc_sd: float) -> float:
    """
    Calculate accumulation index (Rac).

    Rac = AUC_ss / AUC_sd

    Args:
        auc_ss: AUC0-tau at steady state
        auc_sd: AUC0-inf from single dose

    Returns:
        Accumulation index
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_accumulation_index(auc_ss, auc_sd))


def nca_ptf(cmax: float, cmin: float, cavg: float) -> float:
    """
    Calculate Peak-Trough Fluctuation (PTF) percentage.

    PTF = 100 * (Cmax - Cmin) / Cavg

    Args:
        cmax: Maximum concentration in dosing interval
        cmin: Minimum (trough) concentration
        cavg: Average concentration over dosing interval

    Returns:
        Peak-trough fluctuation (%)
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_ptf(cmax, cmin, cavg))


def nca_swing(cmax: float, cmin: float) -> float:
    """
    Calculate Swing percentage.

    Swing = 100 * (Cmax - Cmin) / Cmin

    Args:
        cmax: Maximum concentration in dosing interval
        cmin: Minimum (trough) concentration

    Returns:
        Swing (%)
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_swing(cmax, cmin))


def nca_linearity_index(doses: List[float], aucs: List[float]) -> Dict[str, Any]:
    """
    Assess dose proportionality using power model.

    AUC = alpha * Dose^beta
    Linear if beta is approximately 1.0.

    Args:
        doses: List of dose levels
        aucs: List of corresponding AUC values

    Returns:
        Dict with keys: beta, r_squared, is_linear
    """
    jl = _require_julia()
    d = _to_julia_float_vector(jl, doses)
    a = _to_julia_float_vector(jl, aucs)
    result = jl.NeoPKPD.nca_linearity_index(d, a)
    return {
        "beta": float(result.beta),
        "r_squared": float(result.r_squared),
        "is_linear": bool(result.is_linear),
    }


def nca_time_to_steady_state(lambda_z: float, fraction: float = 0.90) -> float:
    """
    Estimate time to reach a fraction of steady state.

    t_ss = -ln(1 - fraction) / lambda_z

    Args:
        lambda_z: Terminal elimination rate constant
        fraction: Fraction of steady state (default: 0.90)

    Returns:
        Time to reach specified fraction of steady state
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_time_to_steady_state(lambda_z, fraction=fraction))


def nca_accumulation_predicted(lambda_z: float, tau: float) -> float:
    """
    Predict accumulation index from lambda_z and dosing interval.

    Rac_pred = 1 / (1 - exp(-lambda_z * tau))

    This is the theoretical accumulation for a one-compartment model.

    Args:
        lambda_z: Terminal elimination rate constant
        tau: Dosing interval

    Returns:
        Predicted accumulation index
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_accumulation_predicted(lambda_z, tau))


def nca_accumulation_cmax(cmax_ss: float, cmax_sd: float) -> float:
    """
    Calculate Cmax accumulation ratio.

    Rac_Cmax = Cmax_ss / Cmax_sd

    Args:
        cmax_ss: Cmax at steady state
        cmax_sd: Cmax from single dose

    Returns:
        Cmax accumulation ratio
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_accumulation_cmax(cmax_ss, cmax_sd))


def nca_accumulation_cmin(cmin_ss: float, c_at_tau_sd: float) -> float:
    """
    Calculate Cmin accumulation ratio.

    Rac_Cmin = Cmin_ss / C(tau)_sd

    Args:
        cmin_ss: Cmin (trough) at steady state
        c_at_tau_sd: Concentration at time tau after single dose

    Returns:
        Cmin accumulation ratio
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_accumulation_cmin(cmin_ss, c_at_tau_sd))


def nca_dose_normalized_auc(auc: float, dose: float) -> float:
    """
    Calculate dose-normalized AUC.

    AUC/D = AUC / Dose

    Args:
        auc: AUC value
        dose: Dose administered

    Returns:
        Dose-normalized AUC
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_dose_normalized_auc(auc, dose))


def nca_dose_normalized_cmax(cmax: float, dose: float) -> float:
    """
    Calculate dose-normalized Cmax.

    Cmax/D = Cmax / Dose

    Args:
        cmax: Maximum concentration
        dose: Dose administered

    Returns:
        Dose-normalized Cmax
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_dose_normalized_cmax(cmax, dose))


def nca_time_to_steady_state_doses(
    lambda_z: float,
    tau: float,
    fraction: float = 0.90
) -> int:
    """
    Estimate number of doses to reach a fraction of steady state.

    n_doses = ceil(t_ss / tau)

    Args:
        lambda_z: Terminal elimination rate constant
        tau: Dosing interval
        fraction: Fraction of steady state (default: 0.90)

    Returns:
        Number of doses to reach steady state
    """
    jl = _require_julia()
    return int(jl.NeoPKPD.nca_time_to_steady_state_doses(lambda_z, tau, fraction=fraction))


def nca_effective_half_life(auc_ss: float, auc_inf_sd: float) -> float:
    """
    Calculate effective half-life from steady-state data.

    Args:
        auc_ss: AUC0-tau at steady state
        auc_inf_sd: AUC0-inf from single dose

    Returns:
        Effective half-life
    """
    jl = _require_julia()
    return float(jl.NeoPKPD.nca_effective_half_life(auc_ss, auc_inf_sd))


# ============================================================================
# Helper Functions
# ============================================================================

def _make_nca_config(
    jl: Any,
    method: str = "lin_log_mixed",
    lambda_z_min_points: int = 3,
    lambda_z_r2_threshold: float = 0.9,
    extrapolation_max_pct: float = 20.0,
    lloq: Optional[float] = None,
) -> Any:
    """Create Julia NCAConfig object."""
    # Select method type
    if method == "linear":
        method_obj = jl.NeoPKPD.LinearMethod()
    elif method == "log_linear":
        method_obj = jl.NeoPKPD.LogLinearMethod()
    else:
        method_obj = jl.NeoPKPD.LinLogMixedMethod()

    if lloq is not None:
        return jl.NeoPKPD.NCAConfig(
            method=method_obj,
            lambda_z_min_points=lambda_z_min_points,
            lambda_z_r2_threshold=lambda_z_r2_threshold,
            extrapolation_max_pct=extrapolation_max_pct,
            lloq=lloq,
        )
    else:
        return jl.NeoPKPD.NCAConfig(
            method=method_obj,
            lambda_z_min_points=lambda_z_min_points,
            lambda_z_r2_threshold=lambda_z_r2_threshold,
            extrapolation_max_pct=extrapolation_max_pct,
        )


def _maybe_float(val: Any) -> Optional[float]:
    """Convert Julia value to float, returning None for Julia nothing."""
    if val is None:
        return None
    try:
        return float(val)
    except (TypeError, ValueError):
        return None
