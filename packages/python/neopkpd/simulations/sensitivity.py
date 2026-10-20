"""
Sensitivity Analysis

This module provides functions for parameter sensitivity analysis.
"""

from typing import Any, Dict, List, Optional

from .._core import (
    _require_julia,
    _to_julia_vector,
    SensitivityMetrics,
    SensitivityResult,
)


def run_sensitivity(
    model: Dict[str, Any],
    grid: Dict[str, Any],
    perturbation: Dict[str, Any],
    observation: str = "conc",
    solver: Optional[Dict[str, Any]] = None,
) -> SensitivityResult:
    """
    Run sensitivity analysis by perturbing a single parameter.

    Compares base simulation to perturbed simulation and computes
    metrics quantifying the impact of the parameter change.

    Args:
        model: Model specification dict with keys:
            - kind: "OneCompIVBolus" or "OneCompOralFirstOrder"
            - params: Dict of model parameters
            - doses: List of {time, amount} dicts
        grid: Simulation grid dict with keys:
            - t0: Start time
            - t1: End time
            - saveat: List of output time points
        perturbation: Perturbation specification dict with keys:
            - name: Name for this sensitivity analysis
            - param: Parameter to perturb (e.g., "CL")
            - delta: Relative perturbation (e.g., 0.01 for 1%)
        observation: Name of observation to analyze (default: "conc")
        solver: Optional solver settings

    Returns:
        SensitivityResult containing base/perturbed series and metrics

    Example:
        >>> result = neopkpd.run_sensitivity(
        ...     model={
        ...         "kind": "OneCompIVBolus",
        ...         "params": {"CL": 1.0, "V": 10.0},
        ...         "doses": [{"time": 0.0, "amount": 100.0}]
        ...     },
        ...     grid={"t0": 0.0, "t1": 24.0, "saveat": [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]},
        ...     perturbation={"name": "cl_sens", "param": "CL", "delta": 0.01}
        ... )
        >>> print(f"Max relative delta: {result.metrics.max_rel_delta}")
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PerturbationPlan = jl.NeoPKPD.PerturbationPlan

    # Parse model
    kind_str = model["kind"]
    params_dict = model["params"]
    doses_raw = model.get("doses", [])
    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses_raw]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    if kind_str == "OneCompIVBolus":
        OneCompIVBolus = jl.NeoPKPD.OneCompIVBolus
        OneCompIVBolusParams = jl.NeoPKPD.OneCompIVBolusParams
        params = OneCompIVBolusParams(float(params_dict["CL"]), float(params_dict["V"]))
        model_spec = ModelSpec(OneCompIVBolus(), "py_sens", params, doses_vec)
    elif kind_str == "OneCompOralFirstOrder":
        OneCompOralFirstOrder = jl.NeoPKPD.OneCompOralFirstOrder
        OneCompOralFirstOrderParams = jl.NeoPKPD.OneCompOralFirstOrderParams
        params = OneCompOralFirstOrderParams(
            float(params_dict["Ka"]),
            float(params_dict["CL"]),
            float(params_dict["V"]),
        )
        model_spec = ModelSpec(OneCompOralFirstOrder(), "py_sens", params, doses_vec)
    else:
        raise ValueError(f"Unsupported model kind: {kind_str}")

    # Parse grid
    grid_jl = SimGrid(
        float(grid["t0"]),
        float(grid["t1"]),
        [float(x) for x in grid["saveat"]],
    )

    # Parse solver
    if solver is None:
        solver_jl = SolverSpec(jl.Symbol("Tsit5"), 1e-10, 1e-12, 10**7)
    else:
        solver_jl = SolverSpec(
            jl.Symbol(solver.get("alg", "Tsit5")),
            float(solver.get("reltol", 1e-10)),
            float(solver.get("abstol", 1e-12)),
            int(solver.get("maxiters", 10**7)),
        )

    # Create perturbation plan
    plan = PerturbationPlan(
        perturbation["name"],
        jl.Symbol(perturbation["param"]),
        float(perturbation["delta"]),
    )

    # Run sensitivity analysis
    res = jl.NeoPKPD.run_sensitivity(model_spec, grid_jl, solver_jl, plan, jl.Symbol(observation))

    return SensitivityResult(
        plan_name=str(res.plan.name),
        observation=str(res.observation),
        base_series=list(res.base_metric_series),
        pert_series=list(res.pert_metric_series),
        metrics=SensitivityMetrics(
            max_abs_delta=float(res.metrics.max_abs_delta),
            max_rel_delta=float(res.metrics.max_rel_delta),
            l2_norm_delta=float(res.metrics.l2_norm_delta),
        ),
        metadata=dict(res.metadata),
    )
