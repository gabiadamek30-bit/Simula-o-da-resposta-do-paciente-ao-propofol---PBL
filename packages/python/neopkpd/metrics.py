"""
PK/PD Metrics

This module provides functions for computing pharmacokinetic and
pharmacodynamic metrics from simulation results.
"""

import math
from typing import Any, Dict

from ._core import _require_julia, _to_julia_float_vector


def cmax(result: Dict[str, Any], observation: str = "conc") -> float:
    """
    Compute maximum concentration (Cmax) from simulation result.

    Args:
        result: Simulation result dict (from simulate_* functions)
        observation: Name of observation to analyze (default: "conc")

    Returns:
        float: Maximum value of the specified observation

    Example:
        >>> result = neopkpd.simulate_pk_iv_bolus(cl=1.0, v=10.0, ...)
        >>> print(f"Cmax: {neopkpd.cmax(result)}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, result["t"])
    y = _to_julia_float_vector(jl, result["observations"][observation])
    return float(jl.NeoPKPD.cmax(t, y))


def auc_trapezoid(result: Dict[str, Any], observation: str = "conc") -> float:
    """
    Compute area under the curve (AUC) using trapezoidal rule.

    Args:
        result: Simulation result dict (from simulate_* functions)
        observation: Name of observation to analyze (default: "conc")

    Returns:
        float: AUC computed over the simulation time grid

    Example:
        >>> result = neopkpd.simulate_pk_iv_bolus(cl=1.0, v=10.0, ...)
        >>> print(f"AUC: {neopkpd.auc_trapezoid(result)}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, result["t"])
    y = _to_julia_float_vector(jl, result["observations"][observation])
    return float(jl.NeoPKPD.auc_trapezoid(t, y))


def emin(result: Dict[str, Any], observation: str = "effect") -> float:
    """
    Compute minimum value of a response observation (Emin).

    Typically used for PD endpoints where the drug causes suppression
    below baseline.

    Args:
        result: Simulation result dict (from simulate_* functions)
        observation: Name of observation to analyze (default: "effect")

    Returns:
        float: Minimum value of the specified observation

    Example:
        >>> result = neopkpd.simulate_pkpd_indirect_response(...)
        >>> print(f"Emin: {neopkpd.emin(result, 'response')}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, result["t"])
    y = _to_julia_float_vector(jl, result["observations"][observation])
    return float(jl.NeoPKPD.emin(t, y))


def time_below(result: Dict[str, Any], threshold: float, observation: str = "conc") -> float:
    """
    Compute total time where observation is below a threshold.

    Uses left-constant interpolation rule: for interval [t[i-1], t[i]],
    uses y[i-1] to determine if interval is below threshold.

    Args:
        result: Simulation result dict (from simulate_* functions)
        threshold: Threshold value
        observation: Name of observation to analyze (default: "conc")

    Returns:
        float: Total time spent below the threshold

    Example:
        >>> result = neopkpd.simulate_pk_iv_bolus(cl=1.0, v=10.0, ...)
        >>> time_subtherapeutic = neopkpd.time_below(result, threshold=1.0)
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, result["t"])
    y = _to_julia_float_vector(jl, result["observations"][observation])
    return float(jl.NeoPKPD.time_below(t, y, float(threshold)))


def auc_above_baseline(result: Dict[str, Any], baseline: float, observation: str = "effect") -> float:
    """
    Compute AUC of the area where baseline exceeds the observation.

    Measures the "suppression burden" - the integrated area where
    the response is below baseline. Useful for indirect response models.

    Uses: AUC of max(0, baseline - y) over the time grid.

    Args:
        result: Simulation result dict (from simulate_* functions)
        baseline: Baseline value to compare against
        observation: Name of observation to analyze (default: "effect")

    Returns:
        float: Integrated suppression area

    Example:
        >>> result = neopkpd.simulate_pkpd_indirect_response(..., r0=100.0, ...)
        >>> suppression = neopkpd.auc_above_baseline(result, baseline=100.0, observation='response')
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, result["t"])
    y = _to_julia_float_vector(jl, result["observations"][observation])
    return float(jl.NeoPKPD.auc_above_baseline(t, y, float(baseline)))


def tmax(result: Dict[str, Any], observation: str = "conc") -> float:
    """
    Compute time of maximum concentration (Tmax).

    Args:
        result: Simulation result dict (from simulate_* functions)
        observation: Name of observation to analyze (default: "conc")

    Returns:
        float: Time at which maximum concentration occurs

    Example:
        >>> result = neopkpd.simulate_pk_oral_first_order(ka=0.5, cl=1.0, v=10.0, ...)
        >>> print(f"Tmax: {neopkpd.tmax(result)}")
    """
    t = result["t"]
    y = result["observations"][observation]
    max_idx = y.index(max(y))
    return float(t[max_idx])


def half_life(cl: float, v: float) -> float:
    """
    Compute elimination half-life from clearance and volume.

    Formula: t1/2 = ln(2) * V / CL

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution

    Returns:
        float: Elimination half-life (same time units as CL)

    Example:
        >>> t_half = neopkpd.half_life(cl=1.0, v=10.0)
        >>> print(f"Half-life: {t_half} hours")
    """
    return math.log(2) * v / cl
