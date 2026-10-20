"""
Global Sensitivity Analysis (GSA) for NeoPKPD.

This module provides Sobol' and Morris methods for comprehensive
global sensitivity analysis of PK/PD models.
"""

from typing import Any, Callable, Dict, List, Optional, Union

from neopkpd._core import (
    MorrisIndex,
    MorrisResult,
    SobolIndex,
    SobolResult,
    _create_dose_events,
    _require_julia,
    _to_julia_float_vector,
)

__all__ = [
    "run_sobol_sensitivity",
    "run_morris_sensitivity",
    "SobolResult",
    "SobolIndex",
    "MorrisResult",
    "MorrisIndex",
]


def run_sobol_sensitivity(
    model: Dict[str, Any],
    grid: Dict[str, Any],
    bounds: Dict[str, tuple],
    *,
    base_sample_size: int = 1024,
    compute_second_order: bool = False,
    bootstrap_samples: int = 1000,
    bootstrap_ci_level: float = 0.95,
    observation: str = "conc",
    seed: int = 12345,
    aggregator: str = "mean",
) -> SobolResult:
    """
    Run Sobol' variance-based global sensitivity analysis.

    Computes first-order (Si) and total-order (STi) sensitivity indices
    for each parameter using the Saltelli sampling scheme. Optionally
    computes second-order indices (Sij) for parameter interactions.

    Args:
        model: Model specification dict with keys:
            - kind: Model type (e.g., "OneCompIVBolus")
            - params: Dict of parameter values
            - doses: List of dose dicts with 'time' and 'amount'
        grid: Simulation grid dict with keys:
            - t0: Start time
            - t1: End time
            - saveat: List of time points or step size
        bounds: Parameter bounds as Dict[str, (lower, upper)]
            e.g., {"CL": (0.5, 2.0), "V": (5.0, 20.0)}
        base_sample_size: N in Saltelli scheme. Total evals = N*(d+2)
            With second_order=True: N*(d+2 + d*(d-1)/2)
        compute_second_order: If True, compute Sij interaction indices.
            This requires additional model evaluations but reveals
            important parameter interactions for model reduction.
        bootstrap_samples: Number of bootstrap samples for CIs (0 = no CI)
        bootstrap_ci_level: Confidence level (default 0.95)
        observation: Model output to analyze (default "conc")
        seed: Random seed for reproducibility
        aggregator: How to aggregate time series ("mean", "max", "auc")

    Returns:
        SobolResult with sensitivity indices for each parameter.
        If compute_second_order=True, result.second_order contains
        Dict[Tuple[str, str], float] with Sij indices.

    Example:
        >>> # First-order only
        >>> result = run_sobol_sensitivity(
        ...     model={"kind": "OneCompIVBolus", "params": {"CL": 1.0, "V": 10.0},
        ...            "doses": [{"time": 0.0, "amount": 100.0}]},
        ...     grid={"t0": 0.0, "t1": 24.0, "saveat": list(range(25))},
        ...     bounds={"CL": (0.5, 2.0), "V": (5.0, 20.0)},
        ...     base_sample_size=512
        ... )
        >>> print(f"CL importance: Si={result.indices['CL'].Si:.3f}")

        >>> # With second-order interactions
        >>> result = run_sobol_sensitivity(
        ...     model={"kind": "OneCompIVBolus", "params": {"CL": 1.0, "V": 10.0},
        ...            "doses": [{"time": 0.0, "amount": 100.0}]},
        ...     grid={"t0": 0.0, "t1": 24.0, "saveat": list(range(25))},
        ...     bounds={"CL": (0.5, 2.0), "V": (5.0, 20.0)},
        ...     base_sample_size=512,
        ...     compute_second_order=True
        ... )
        >>> print(f"CL-V interaction: Sij={result.second_order[('CL', 'V')]:.3f}")
    """
    jl = _require_julia()
    NC = jl.NeoPKPD

    # Build model spec
    model_kind = getattr(NC, model["kind"])()
    params_type = getattr(NC, f"{model['kind']}Params")
    params_obj = params_type(*[float(v) for v in model["params"].values()])
    doses = _create_dose_events(jl, model["doses"])
    spec = NC.ModelSpec(model_kind, "gsa_model", params_obj, doses)

    # Build grid
    saveat = grid.get("saveat", 0.5)
    if isinstance(saveat, list):
        saveat_vec = _to_julia_float_vector(jl, saveat)
        grid_obj = NC.SimGrid(float(grid["t0"]), float(grid["t1"]), saveat_vec)
    else:
        grid_obj = NC.SimGrid(float(grid["t0"]), float(grid["t1"]), float(saveat))

    # Build solver
    solver_obj = NC.SolverSpec()

    # Build parameter bounds
    params_list = list(bounds.keys())
    params_vec = jl.Vector[jl.Symbol]([jl.Symbol(p) for p in params_list])
    lower_vec = _to_julia_float_vector(jl, [bounds[p][0] for p in params_list])
    upper_vec = _to_julia_float_vector(jl, [bounds[p][1] for p in params_list])
    bounds_obj = NC.ParameterBounds(params_vec, lower_vec, upper_vec)

    # Build method
    method = NC.SobolMethod(
        base_sample_size=base_sample_size,
        compute_second_order=compute_second_order,
        bootstrap_samples=bootstrap_samples,
        bootstrap_ci_level=bootstrap_ci_level,
    )

    # Build GSA spec
    gsa_spec = NC.GlobalSensitivitySpec(
        method, bounds_obj,
        observation=jl.Symbol(observation),
        seed=jl.UInt64(seed),
    )

    # Build aggregator
    if aggregator == "mean":
        agg_fn = jl.Statistics.mean
    elif aggregator == "max":
        agg_fn = jl.maximum
    elif aggregator == "auc":
        # Use trapezoidal integration
        agg_fn = jl.seval("x -> sum(diff(x) .* (x[1:end-1] .+ x[2:end]) ./ 2)")
    else:
        agg_fn = jl.Statistics.mean

    # Run analysis
    result = NC.run_sobol_sensitivity(spec, grid_obj, solver_obj, gsa_spec, aggregator=agg_fn)

    # Convert to Python
    indices = {}
    for param in params_list:
        idx = result.indices[jl.Symbol(param)]
        indices[param] = SobolIndex(
            Si=float(idx.Si),
            Si_ci_lower=float(idx.Si_ci_lower),
            Si_ci_upper=float(idx.Si_ci_upper),
            STi=float(idx.STi),
            STi_ci_lower=float(idx.STi_ci_lower),
            STi_ci_upper=float(idx.STi_ci_upper),
        )

    # Convert second-order indices if computed
    second_order = None
    if result.second_order is not None and compute_second_order:
        second_order = {}
        # Julia uses Symbol tuples as keys, convert to Python string tuples
        for key in result.second_order.keys():
            # Key is a Julia Tuple{Symbol, Symbol}
            param_i = str(key[0])
            param_j = str(key[1])
            sij = float(result.second_order[key])
            second_order[(param_i, param_j)] = sij

    metadata = {str(k): v for k, v in result.metadata.items()}

    return SobolResult(
        params=params_list,
        indices=indices,
        second_order=second_order,
        n_evaluations=int(result.n_evaluations),
        convergence_metric=float(result.convergence_metric),
        output_variance=float(result.output_variance),
        computation_time=float(result.computation_time),
        metadata=metadata,
    )


def run_morris_sensitivity(
    model: Dict[str, Any],
    grid: Dict[str, Any],
    bounds: Dict[str, tuple],
    *,
    n_trajectories: int = 20,
    n_levels: int = 4,
    observation: str = "conc",
    seed: int = 12345,
    aggregator: str = "mean",
) -> MorrisResult:
    """
    Run Morris Elementary Effects screening analysis.

    A computationally efficient screening method to identify important
    parameters before more expensive variance-based analysis.

    Args:
        model: Model specification dict with keys:
            - kind: Model type (e.g., "OneCompIVBolus")
            - params: Dict of parameter values
            - doses: List of dose dicts with 'time' and 'amount'
        grid: Simulation grid dict with keys:
            - t0: Start time
            - t1: End time
            - saveat: List of time points or step size
        bounds: Parameter bounds as Dict[str, (lower, upper)]
        n_trajectories: Number of Morris trajectories. Total evals = r*(d+1)
        n_levels: Number of grid levels (typically 4)
        observation: Model output to analyze (default "conc")
        seed: Random seed for reproducibility
        aggregator: How to aggregate time series ("mean", "max", "auc")

    Returns:
        MorrisResult with mu, mu_star, sigma for each parameter

    Example:
        >>> result = run_morris_sensitivity(
        ...     model={"kind": "OneCompIVBolus", "params": {"CL": 1.0, "V": 10.0},
        ...            "doses": [{"time": 0.0, "amount": 100.0}]},
        ...     grid={"t0": 0.0, "t1": 24.0, "saveat": list(range(25))},
        ...     bounds={"CL": (0.5, 2.0), "V": (5.0, 20.0)},
        ...     n_trajectories=20
        ... )
        >>> print(f"CL importance: mu*={result.indices['CL'].mu_star:.3f}")
    """
    jl = _require_julia()
    NC = jl.NeoPKPD

    # Build model spec
    model_kind = getattr(NC, model["kind"])()
    params_type = getattr(NC, f"{model['kind']}Params")
    params_obj = params_type(*[float(v) for v in model["params"].values()])
    doses = _create_dose_events(jl, model["doses"])
    spec = NC.ModelSpec(model_kind, "gsa_model", params_obj, doses)

    # Build grid
    saveat = grid.get("saveat", 0.5)
    if isinstance(saveat, list):
        saveat_vec = _to_julia_float_vector(jl, saveat)
        grid_obj = NC.SimGrid(float(grid["t0"]), float(grid["t1"]), saveat_vec)
    else:
        grid_obj = NC.SimGrid(float(grid["t0"]), float(grid["t1"]), float(saveat))

    # Build solver
    solver_obj = NC.SolverSpec()

    # Build parameter bounds
    params_list = list(bounds.keys())
    params_vec = jl.Vector[jl.Symbol]([jl.Symbol(p) for p in params_list])
    lower_vec = _to_julia_float_vector(jl, [bounds[p][0] for p in params_list])
    upper_vec = _to_julia_float_vector(jl, [bounds[p][1] for p in params_list])
    bounds_obj = NC.ParameterBounds(params_vec, lower_vec, upper_vec)

    # Build method
    method = NC.MorrisMethod(
        n_trajectories=n_trajectories,
        n_levels=n_levels,
    )

    # Build GSA spec
    gsa_spec = NC.GlobalSensitivitySpec(
        method, bounds_obj,
        observation=jl.Symbol(observation),
        seed=jl.UInt64(seed),
    )

    # Build aggregator
    if aggregator == "mean":
        agg_fn = jl.Statistics.mean
    elif aggregator == "max":
        agg_fn = jl.maximum
    elif aggregator == "auc":
        agg_fn = jl.seval("x -> sum(diff(x) .* (x[1:end-1] .+ x[2:end]) ./ 2)")
    else:
        agg_fn = jl.Statistics.mean

    # Run analysis
    result = NC.run_morris_sensitivity(spec, grid_obj, solver_obj, gsa_spec, aggregator=agg_fn)

    # Convert to Python
    indices = {}
    elementary_effects = {}
    for param in params_list:
        idx = result.indices[jl.Symbol(param)]
        indices[param] = MorrisIndex(
            mu=float(idx.mu),
            mu_star=float(idx.mu_star),
            sigma=float(idx.sigma),
        )
        ees = result.elementary_effects[jl.Symbol(param)]
        elementary_effects[param] = [float(e) for e in ees]

    metadata = {str(k): v for k, v in result.metadata.items()}

    return MorrisResult(
        params=params_list,
        indices=indices,
        elementary_effects=elementary_effects,
        n_evaluations=int(result.n_evaluations),
        computation_time=float(result.computation_time),
        metadata=metadata,
    )
