"""
Artifact Operations

This module provides functions for saving and replaying simulation artifacts.
Artifacts contain all information needed to reproduce a simulation exactly.
"""

from pathlib import Path
from typing import Any, Dict, Optional, Union

from ._core import (
    _require_julia,
    _simresult_to_py,
    _popresult_to_py,
    _to_julia_vector,
)


def replay_artifact(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Replay a simulation from a saved artifact file.

    This function re-executes a simulation using the exact parameters stored
    in an artifact file, validating reproducibility across platforms.

    Supports:
    - Single PK/PD artifacts
    - Population artifacts
    - Sensitivity analysis artifacts (single and population)

    Args:
        path: Path to the artifact JSON file

    Returns:
        Dict containing the simulation results. Structure depends on artifact type:
        - Single: {t, states, observations, metadata}
        - Population: {individuals, params, summaries, metadata}
        - Sensitivity: {plan, observation, base_series, pert_series, metrics, metadata}

    Example:
        >>> result = neopkpd.replay_artifact("validation/golden/pk_iv_bolus.json")
        >>> print(f"Time points: {len(result['t'])}")
    """
    jl = _require_julia()
    artifact = jl.NeoPKPD.read_execution_json(str(Path(path).resolve()))

    atype = "single"
    if "artifact_type" in artifact:
        atype = str(artifact["artifact_type"])

    if atype == "population":
        res = jl.NeoPKPD.replay_population_execution(artifact)
        return _popresult_to_py(res)

    if atype == "sensitivity_single":
        res = jl.NeoPKPD.replay_sensitivity_execution(artifact)
        return {
            "plan": {"name": str(res.plan.name)},
            "observation": str(res.observation),
            "base_series": list(res.base_metric_series),
            "pert_series": list(res.pert_metric_series),
            "metrics": {
                "max_abs_delta": float(res.metrics.max_abs_delta),
                "max_rel_delta": float(res.metrics.max_rel_delta),
                "l2_norm_delta": float(res.metrics.l2_norm_delta),
            },
            "metadata": dict(res.metadata),
        }

    if atype == "sensitivity_population":
        res = jl.NeoPKPD.replay_population_sensitivity_execution(artifact)
        return {
            "plan": {"name": str(res.plan.name)},
            "observation": str(res.observation),
            "probs": list(res.probs),
            "base_mean": list(res.base_summary_mean),
            "pert_mean": list(res.pert_summary_mean),
            "metrics_mean": {
                "max_abs_delta": float(res.metrics_mean.max_abs_delta),
                "max_rel_delta": float(res.metrics_mean.max_rel_delta),
                "l2_norm_delta": float(res.metrics_mean.l2_norm_delta),
            },
            "metadata": dict(res.metadata),
        }

    res = jl.NeoPKPD.replay_execution(artifact)
    return _simresult_to_py(res)


def write_single_artifact(
    path: Union[str, Path],
    *,
    model: Dict[str, Any],
    grid: Dict[str, Any],
    solver: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Run a simulation and write a complete artifact file.

    The artifact file contains all information needed to reproduce the
    simulation, including the model specification, grid, solver settings,
    and results.

    Args:
        path: Output path for the artifact JSON file
        model: Model specification dict with keys:
            - kind: "OneCompIVBolus" or "OneCompOralFirstOrder"
            - params: Dict of model parameters (CL, V, Ka)
            - doses: List of {time, amount} dicts
        grid: Simulation grid dict with keys:
            - t0: Start time
            - t1: End time
            - saveat: List of output time points
        solver: Optional solver settings dict with keys:
            - alg: Algorithm name (default: "Tsit5")
            - reltol: Relative tolerance (default: 1e-10)
            - abstol: Absolute tolerance (default: 1e-12)
            - maxiters: Maximum iterations (default: 10^7)

    Example:
        >>> neopkpd.write_single_artifact(
        ...     "my_simulation.json",
        ...     model={
        ...         "kind": "OneCompIVBolus",
        ...         "params": {"CL": 1.0, "V": 10.0},
        ...         "doses": [{"time": 0.0, "amount": 100.0}]
        ...     },
        ...     grid={"t0": 0.0, "t1": 24.0, "saveat": [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]}
        ... )
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    kind = model["kind"]
    dose_list = [DoseEvent(float(d["time"]), float(d["amount"])) for d in model["doses"]]
    doses = _to_julia_vector(jl, dose_list, DoseEvent)

    if kind == "OneCompIVBolus":
        OneCompIVBolus = jl.NeoPKPD.OneCompIVBolus
        Params = jl.NeoPKPD.OneCompIVBolusParams
        params = Params(float(model["params"]["CL"]), float(model["params"]["V"]))
        spec = ModelSpec(OneCompIVBolus(), "py_artifact", params, doses)
    elif kind == "OneCompOralFirstOrder":
        OneCompOralFirstOrder = jl.NeoPKPD.OneCompOralFirstOrder
        Params = jl.NeoPKPD.OneCompOralFirstOrderParams
        params = Params(
            float(model["params"]["KA"]),
            float(model["params"]["CL"]),
            float(model["params"]["V"]),
        )
        spec = ModelSpec(OneCompOralFirstOrder(), "py_artifact", params, doses)
    else:
        raise ValueError(f"Unsupported model kind: {kind}")

    grid_jl = SimGrid(
        float(grid["t0"]),
        float(grid["t1"]),
        [float(x) for x in grid["saveat"]],
    )

    if solver is None:
        solver_jl = SolverSpec(jl.Symbol("Tsit5"), 1e-10, 1e-12, 10**7)
    else:
        solver_jl = SolverSpec(
            jl.Symbol(solver.get("alg", "Tsit5")),
            float(solver.get("reltol", 1e-10)),
            float(solver.get("abstol", 1e-12)),
            int(solver.get("maxiters", 10**7)),
        )

    res = jl.NeoPKPD.simulate(spec, grid_jl, solver_jl)

    jl.NeoPKPD.write_execution_json(
        str(Path(path).resolve()),
        model_spec=spec,
        grid=grid_jl,
        solver=solver_jl,
        result=res,
    )
