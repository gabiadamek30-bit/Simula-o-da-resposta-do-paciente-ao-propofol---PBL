"""
NeoPKPD NCA Analysis - Full NCA workflow functions.

This module provides the main NCA entry points including full NCA analysis
and population NCA.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Union

from .._core import _require_julia, _to_julia_float_vector
from .metrics import _make_nca_config, _maybe_float


# ============================================================================
# Configuration and Result Types
# ============================================================================

@dataclass
class NCAConfig:
    """
    Configuration for NCA analysis.

    Attributes:
        method: AUC calculation method ("linear", "log_linear", "lin_log_mixed")
        lambda_z_min_points: Minimum points for lambda_z regression
        lambda_z_r2_threshold: Minimum R-squared for acceptable lambda_z
        extrapolation_max_pct: Maximum allowed AUC extrapolation %
        significant_digits: Significant digits for reporting
        blq_handling: How to handle BLQ values ("zero", "missing", "lloq_half")
        lloq: Lower limit of quantification
    """
    method: str = "lin_log_mixed"
    lambda_z_min_points: int = 3
    lambda_z_r2_threshold: float = 0.9
    extrapolation_max_pct: float = 20.0
    significant_digits: int = 3
    blq_handling: str = "zero"
    lloq: Optional[float] = None

    def _to_julia(self, jl: Any) -> Any:
        """Convert to Julia NCAConfig object."""
        # Select method type
        if self.method == "linear":
            method_obj = jl.NeoPKPD.LinearMethod()
        elif self.method == "log_linear":
            method_obj = jl.NeoPKPD.LogLinearMethod()
        else:
            method_obj = jl.NeoPKPD.LinLogMixedMethod()

        # Select BLQ handling
        if self.blq_handling == "missing":
            blq_obj = jl.NeoPKPD.BLQMissing()
        elif self.blq_handling == "lloq_half":
            blq_obj = jl.NeoPKPD.BLQLLOQHalf()
        else:
            blq_obj = jl.NeoPKPD.BLQZero()

        if self.lloq is not None:
            return jl.NeoPKPD.NCAConfig(
                method=method_obj,
                lambda_z_min_points=self.lambda_z_min_points,
                lambda_z_r2_threshold=self.lambda_z_r2_threshold,
                extrapolation_max_pct=self.extrapolation_max_pct,
                significant_digits=self.significant_digits,
                blq_handling=blq_obj,
                lloq=self.lloq,
            )
        else:
            return jl.NeoPKPD.NCAConfig(
                method=method_obj,
                lambda_z_min_points=self.lambda_z_min_points,
                lambda_z_r2_threshold=self.lambda_z_r2_threshold,
                extrapolation_max_pct=self.extrapolation_max_pct,
                significant_digits=self.significant_digits,
                blq_handling=blq_obj,
            )


@dataclass
class NCAResult:
    """
    Complete NCA analysis results.

    Attributes:
        cmax: Maximum observed concentration
        tmax: Time of Cmax
        cmin: Minimum concentration (multiple dose only)
        clast: Last measurable concentration
        tlast: Time of last measurable concentration
        cavg: Average concentration over dosing interval
        auc_0_t: AUC from 0 to last measurable concentration
        auc_0_inf: AUC extrapolated to infinity
        auc_extra_pct: Percent of AUC extrapolated
        auc_0_tau: AUC over dosing interval
        aumc_0_t: AUMC from 0 to last concentration
        aumc_0_inf: AUMC extrapolated to infinity
        lambda_z: Terminal elimination rate constant
        t_half: Terminal half-life
        mrt: Mean residence time
        cl_f: Apparent clearance
        vz_f: Apparent volume of distribution
        vss: Volume at steady state
        accumulation_index: AUC ratio at steady state
        ptf: Peak-trough fluctuation (%)
        swing: Swing (%)
        cmax_dn: Dose-normalized Cmax
        auc_dn: Dose-normalized AUC
        quality_flags: Quality assessment flags
        warnings: Analysis warnings
        metadata: Additional metadata
    """
    cmax: float
    tmax: float
    cmin: Optional[float]
    clast: float
    tlast: float
    cavg: Optional[float]
    auc_0_t: float
    auc_0_inf: Optional[float]
    auc_extra_pct: Optional[float]
    auc_0_tau: Optional[float]
    aumc_0_t: float
    aumc_0_inf: Optional[float]
    lambda_z: Optional[float]
    t_half: Optional[float]
    mrt: Optional[float]
    cl_f: Optional[float]
    vz_f: Optional[float]
    vss: Optional[float]
    accumulation_index: Optional[float]
    ptf: Optional[float]
    swing: Optional[float]
    cmax_dn: Optional[float]
    auc_dn: Optional[float]
    quality_flags: List[str]
    warnings: List[str]
    metadata: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "cmax": self.cmax,
            "tmax": self.tmax,
            "cmin": self.cmin,
            "clast": self.clast,
            "tlast": self.tlast,
            "cavg": self.cavg,
            "auc_0_t": self.auc_0_t,
            "auc_0_inf": self.auc_0_inf,
            "auc_extra_pct": self.auc_extra_pct,
            "auc_0_tau": self.auc_0_tau,
            "aumc_0_t": self.aumc_0_t,
            "aumc_0_inf": self.aumc_0_inf,
            "lambda_z": self.lambda_z,
            "t_half": self.t_half,
            "mrt": self.mrt,
            "cl_f": self.cl_f,
            "vz_f": self.vz_f,
            "vss": self.vss,
            "accumulation_index": self.accumulation_index,
            "ptf": self.ptf,
            "swing": self.swing,
            "cmax_dn": self.cmax_dn,
            "auc_dn": self.auc_dn,
            "quality_flags": self.quality_flags,
            "warnings": self.warnings,
            "metadata": self.metadata,
        }


# ============================================================================
# Main NCA Functions
# ============================================================================

def run_nca(
    times: List[float],
    concentrations: List[float],
    dose: float,
    config: Optional[NCAConfig] = None,
    dosing_type: str = "single",
    tau: Optional[float] = None,
    route: str = "extravascular",
    t_inf: float = 0.0,
) -> NCAResult:
    """
    Perform complete Non-Compartmental Analysis on concentration-time data.

    This is the main entry point for NCA analysis, following FDA/EMA guidance.

    Args:
        times: Time points (sorted, ascending)
        concentrations: Concentration values
        dose: Administered dose
        config: NCA configuration (default: standard FDA/EMA settings)
        dosing_type: "single", "multiple", or "steady_state"
        tau: Dosing interval (required for multiple dose)
        route: Administration route ("extravascular", "iv_bolus", "iv_infusion")
        t_inf: Infusion duration (for iv_infusion)

    Returns:
        NCAResult with all NCA metrics

    Example:
        >>> t = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        >>> c = [0.0, 1.2, 2.0, 1.8, 1.2, 0.6, 0.3, 0.075]
        >>> result = run_nca(t, c, dose=100.0)
        >>> print(f"Cmax: {result.cmax}")
        >>> print(f"AUC0-inf: {result.auc_0_inf}")
        >>> print(f"t1/2: {result.t_half}")

    For multiple dose analysis:
        >>> result = run_nca(t, c, dose=100.0, dosing_type="multiple", tau=12.0)
        >>> print(f"Cmin: {result.cmin}")
        >>> print(f"Cavg: {result.cavg}")
        >>> print(f"PTF: {result.ptf}%")
    """
    jl = _require_julia()

    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)

    if config is None:
        config = NCAConfig()

    jl_config = config._to_julia(jl)
    dosing_sym = jl.Symbol(dosing_type)
    route_sym = jl.Symbol(route)

    if tau is not None:
        result = jl.NeoPKPD.run_nca(
            t, c, dose,
            config=jl_config,
            dosing_type=dosing_sym,
            tau=tau,
            route=route_sym,
            t_inf=t_inf,
        )
    else:
        result = jl.NeoPKPD.run_nca(
            t, c, dose,
            config=jl_config,
            dosing_type=dosing_sym,
            route=route_sym,
            t_inf=t_inf,
        )

    return _nca_result_to_py(result)


def run_population_nca(
    population_result: Dict[str, Any],
    dose: float,
    config: Optional[NCAConfig] = None,
    observation: str = "conc",
    dosing_type: str = "single",
    tau: Optional[float] = None,
    route: str = "extravascular",
) -> List[NCAResult]:
    """
    Perform NCA analysis on each individual in a population simulation result.

    Args:
        population_result: Population simulation result dictionary
        dose: Administered dose
        config: NCA configuration
        observation: Observation to analyze (default: "conc")
        dosing_type: Dosing type
        tau: Dosing interval for multiple dose
        route: Administration route

    Returns:
        List of NCAResult for each individual

    Example:
        >>> # Run population simulation first
        >>> pop_result = neopkpd.simulate_population_iv_bolus(...)
        >>>
        >>> # Then run NCA on each individual
        >>> nca_results = run_population_nca(pop_result, dose=100.0)
        >>> for i, result in enumerate(nca_results):
        ...     print(f"Subject {i+1}: Cmax={result.cmax}, AUC={result.auc_0_inf}")
    """
    results = []

    for individual in population_result["individuals"]:
        t = individual["t"]
        c = individual["observations"][observation]

        result = run_nca(
            t, c, dose,
            config=config,
            dosing_type=dosing_type,
            tau=tau,
            route=route,
        )
        results.append(result)

    return results


def summarize_population_nca(
    nca_results: List[NCAResult],
    parameters: Optional[List[str]] = None,
) -> Dict[str, Dict[str, float]]:
    """
    Summarize NCA results across a population.

    Args:
        nca_results: List of NCA results from population
        parameters: Parameters to summarize (default: cmax, auc_0_inf, t_half, cl_f)

    Returns:
        Dict mapping parameter name to summary statistics

    Example:
        >>> nca_results = run_population_nca(pop_result, dose=100.0)
        >>> summary = summarize_population_nca(nca_results)
        >>> print(f"Cmax mean: {summary['cmax']['mean']}")
        >>> print(f"Cmax CV%: {summary['cmax']['cv_pct']}")
    """
    if parameters is None:
        parameters = ["cmax", "auc_0_inf", "t_half", "cl_f"]

    summaries = {}

    for param in parameters:
        values = []
        for result in nca_results:
            val = getattr(result, param, None)
            if val is not None:
                values.append(val)

        if values:
            n = len(values)
            mean_val = sum(values) / n
            sorted_vals = sorted(values)
            median_val = sorted_vals[n // 2] if n % 2 == 1 else (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2
            variance = sum((v - mean_val) ** 2 for v in values) / (n - 1) if n > 1 else 0.0
            sd_val = variance ** 0.5
            cv_val = (sd_val / mean_val * 100.0) if mean_val != 0 else 0.0

            import math
            gm_val = math.exp(sum(math.log(v) for v in values if v > 0) / len([v for v in values if v > 0])) if all(v > 0 for v in values) else None

            summaries[param] = {
                "n": n,
                "mean": mean_val,
                "median": median_val,
                "sd": sd_val,
                "cv_pct": cv_val,
                "geometric_mean": gm_val,
                "min": min(values),
                "max": max(values),
            }

    return summaries


# ============================================================================
# Helper Functions
# ============================================================================

def _nca_result_to_py(result: Any) -> NCAResult:
    """Convert Julia NCAResult to Python NCAResult."""
    # Extract lambda_z from lambda_z_result
    lambda_z = None
    if result.lambda_z_result is not None:
        lambda_z = _maybe_float(result.lambda_z_result.lambda_z)

    return NCAResult(
        cmax=float(result.cmax),
        tmax=float(result.tmax),
        cmin=_maybe_float(result.cmin),
        clast=float(result.clast),
        tlast=float(result.tlast),
        cavg=_maybe_float(result.cavg),
        auc_0_t=float(result.auc_0_t),
        auc_0_inf=_maybe_float(result.auc_0_inf),
        auc_extra_pct=_maybe_float(result.auc_extra_pct),
        auc_0_tau=_maybe_float(result.auc_0_tau),
        aumc_0_t=float(result.aumc_0_t),
        aumc_0_inf=_maybe_float(result.aumc_0_inf),
        lambda_z=lambda_z,
        t_half=_maybe_float(result.t_half),
        mrt=_maybe_float(result.mrt),
        cl_f=_maybe_float(result.cl_f),
        vz_f=_maybe_float(result.vz_f),
        vss=_maybe_float(result.vss),
        accumulation_index=_maybe_float(result.accumulation_index),
        ptf=_maybe_float(result.ptf),
        swing=_maybe_float(result.swing),
        cmax_dn=_maybe_float(result.cmax_dn),
        auc_dn=_maybe_float(result.auc_dn),
        quality_flags=[str(f) for f in result.quality_flags],
        warnings=list(result.warnings),
        metadata=dict(result.metadata),
    )
