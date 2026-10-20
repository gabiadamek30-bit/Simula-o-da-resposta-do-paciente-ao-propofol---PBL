"""
NeoPKPD Python Bridge

This module provides Python bindings for the NeoPKPD Julia core,
enabling full PK/PD simulation capabilities from Python.

Features:
- Single PK simulations (IV bolus, oral first-order)
- Coupled PK-PD simulations (direct Emax, indirect response turnover)
- Population simulations with IIV, IOV, and covariates
- Sensitivity analysis
- Artifact replay for reproducibility
- PK/PD metrics (Cmax, AUC, Emin, time below threshold)

Example:
    >>> import neopkpd
    >>> neopkpd.init_julia()
    >>> result = neopkpd.simulate_pk_iv_bolus(
    ...     cl=1.0, v=10.0,
    ...     doses=[{"time": 0.0, "amount": 100.0}],
    ...     t0=0.0, t1=24.0,
    ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    ... )
    >>> print(f"Cmax: {neopkpd.cmax(result)}")
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

_JL = None


# ============================================================================
# Data Classes
# ============================================================================

@dataclass(frozen=True)
class RepoPaths:
    """Repository path configuration."""
    repo_root: Path
    core_project: Path


@dataclass(frozen=True)
class SensitivityMetrics:
    """Results from sensitivity analysis."""
    max_abs_delta: float
    max_rel_delta: float
    l2_norm_delta: float


@dataclass(frozen=True)
class SensitivityResult:
    """Full sensitivity analysis result."""
    plan_name: str
    observation: str
    base_series: List[float]
    pert_series: List[float]
    metrics: SensitivityMetrics
    metadata: Dict[str, Any]


@dataclass(frozen=True)
class PopulationSensitivityResult:
    """Population sensitivity analysis result."""
    plan_name: str
    observation: str
    probs: List[float]
    base_mean: List[float]
    pert_mean: List[float]
    metrics_mean: SensitivityMetrics
    metadata: Dict[str, Any]


# ============================================================================
# Initialization
# ============================================================================

def version() -> str:
    """
    Get the NeoPKPD version string.

    Returns:
        str: The version string (e.g., "0.1.0")

    Example:
        >>> neopkpd.version()
        '0.1.0'
    """
    jl = _require_julia()
    return str(jl.NeoPKPD.NEOPKPD_VERSION)


def _detect_repo_root(start: Optional[Path] = None) -> Path:
    """Detect the NeoPKPD repository root directory."""
    here = (start or Path(__file__)).resolve()
    for p in [here] + list(here.parents):
        # Check for packages/core/Project.toml (current structure)
        if (p / "packages" / "core" / "Project.toml").exists():
            return p
        # Legacy: check for core/NeoPKPD/Project.toml
        if (p / "core" / "NeoPKPD" / "Project.toml").exists():
            return p
    raise RuntimeError("Could not locate repo root (packages/core/Project.toml not found).")


def init_julia(repo_root: Optional[Union[str, Path]] = None) -> None:
    """
    Initialize the Julia runtime and load NeoPKPD.

    This function activates the Julia project and loads the core simulation
    engine. It's safe to call multiple times - subsequent calls are no-ops.

    Args:
        repo_root: Optional path to the NeoPKPD repository root.
                   If not provided, auto-detection is used.

    Example:
        >>> import neopkpd
        >>> neopkpd.init_julia()
        >>> # Now ready to run simulations
    """
    global _JL
    if _JL is not None:
        return

    root = Path(repo_root).resolve() if repo_root else _detect_repo_root()

    # Determine core project path - try new structure first, then legacy
    if (root / "packages" / "core" / "Project.toml").exists():
        core_project = root / "packages" / "core"
    else:
        core_project = root / "core" / "NeoPKPD"

    from juliacall import Main as jl  # type: ignore

    # Activate and instantiate the exact Julia project in this repo
    jl.seval("import Pkg")
    jl.Pkg.activate(str(core_project))
    jl.Pkg.instantiate()

    jl.seval("using NeoPKPD")

    _JL = jl


def _require_julia() -> Any:
    """Ensure Julia is initialized and return the Julia main module."""
    if _JL is None:
        init_julia()
    return _JL


def get_julia() -> Any:
    """
    Get the Julia main module, initializing if necessary.

    This is the public API for getting the Julia runtime.
    Use this in modules that need access to Julia.

    Returns:
        The Julia main module with NeoPKPD loaded.
    """
    return _require_julia()


# ============================================================================
# Result Conversion Utilities
# ============================================================================

def _simresult_to_py(res: Any) -> Dict[str, Any]:
    """
    Convert NeoPKPD.SimResult to a Python dictionary.

    Args:
        res: Julia SimResult object

    Returns:
        Dict with keys: t, states, observations, metadata
    """
    t = list(res.t)
    states = {str(k): list(v) for (k, v) in res.states.items()}
    obs = {str(k): list(v) for (k, v) in res.observations.items()}
    meta = dict(res.metadata)
    return {"t": t, "states": states, "observations": obs, "metadata": meta}


def _popresult_to_py(popres: Any) -> Dict[str, Any]:
    """
    Convert NeoPKPD.PopulationResult to a Python dictionary.

    Args:
        popres: Julia PopulationResult object

    Returns:
        Dict with keys: individuals, params, summaries, metadata
    """
    individuals = [_simresult_to_py(r) for r in popres.individuals]
    params = [{str(k): float(v) for (k, v) in d.items()} for d in popres.params]
    summaries = {}

    for (k, s) in popres.summaries.items():
        summaries[str(k)] = {
            "observation": str(s.observation),
            "probs": list(s.probs),
            "mean": list(s.mean),
            "median": list(s.median),
            "quantiles": {str(p): list(v) for (p, v) in s.quantiles.items()},
        }

    meta = dict(popres.metadata)
    return {"individuals": individuals, "params": params, "summaries": summaries, "metadata": meta}


def _to_julia_vector(jl: Any, items: list, item_type: Any) -> Any:
    """Convert a Python list to a Julia Vector of the specified type."""
    vec = jl.Vector[item_type](jl.undef, len(items))
    for i, item in enumerate(items):
        vec[i] = item  # PythonCall uses 0-based indexing from Python
    return vec


def _to_julia_float_vector(jl: Any, items: list) -> Any:
    """Convert a Python list to a Julia Vector{Float64}."""
    vec = jl.Vector[jl.Float64](jl.undef, len(items))
    for i, item in enumerate(items):
        vec[i] = float(item)
    return vec


def _create_dose_event(jl: Any, dose_dict: Dict[str, Any]) -> Any:
    """
    Create a Julia DoseEvent from a Python dict.

    Supports both bolus and infusion administration:
    - Bolus: {"time": 0.0, "amount": 100.0}  (duration defaults to 0.0)
    - Infusion: {"time": 0.0, "amount": 100.0, "duration": 1.0}  (100 mg over 1 hour)

    Args:
        jl: Julia module
        dose_dict: Dict with 'time', 'amount', and optional 'duration'

    Returns:
        Julia DoseEvent object
    """
    DoseEvent = jl.NeoPKPD.DoseEvent
    return DoseEvent(
        float(dose_dict["time"]),
        float(dose_dict["amount"]),
        float(dose_dict.get("duration", 0.0))
    )


def _create_dose_events(jl: Any, doses: List[Dict[str, Any]]) -> Any:
    """
    Create a Julia Vector of DoseEvents from a list of Python dicts.

    Args:
        jl: Julia module
        doses: List of dose dicts with 'time', 'amount', and optional 'duration'

    Returns:
        Julia Vector{DoseEvent}
    """
    DoseEvent = jl.NeoPKPD.DoseEvent
    dose_objs = [_create_dose_event(jl, d) for d in doses]
    return _to_julia_vector(jl, dose_objs, DoseEvent)


# ============================================================================
# Artifact Replay
# ============================================================================

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


# ============================================================================
# Artifact Writing
# ============================================================================

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


# ============================================================================
# Single PK Simulation
# ============================================================================

def simulate_pk_iv_bolus(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a one-compartment IV bolus PK simulation.

    Model equations:
        dA/dt = -k * A
        C = A / V
        k = CL / V

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events, each a dict with 'time' and 'amount'
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (A)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pk_iv_bolus(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     t0=0.0, t1=24.0,
        ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        ... )
        >>> print(f"Cmax: {max(result['observations']['conc'])}")
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    OneCompIVBolus = jl.NeoPKPD.OneCompIVBolus
    OneCompIVBolusParams = jl.NeoPKPD.OneCompIVBolusParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(OneCompIVBolus(), "py_iv_bolus", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_oral_first_order(
    ka: float,
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a one-compartment oral first-order absorption PK simulation.

    Model equations:
        dA_gut/dt = -Ka * A_gut
        dA/dt = Ka * A_gut - k * A
        C = A / V
        k = CL / V

    Args:
        ka: Absorption rate constant (1/time)
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events, each a dict with 'time' and 'amount'
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (A_gut, A)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pk_oral_first_order(
        ...     ka=0.5, cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     t0=0.0, t1=24.0,
        ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        ... )
        >>> print(f"Tmax at approximately: {result['t'][result['observations']['conc'].index(max(result['observations']['conc']))]}")
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    OneCompOralFirstOrder = jl.NeoPKPD.OneCompOralFirstOrder
    OneCompOralFirstOrderParams = jl.NeoPKPD.OneCompOralFirstOrderParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "py_oral_first_order",
        OneCompOralFirstOrderParams(float(ka), float(cl), float(v)),
        doses_vec,
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


# ============================================================================
# Two-Compartment PK Simulations
# ============================================================================

def simulate_pk_twocomp_iv_bolus(
    cl: float,
    v1: float,
    q: float,
    v2: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a two-compartment IV bolus PK simulation.

    Model structure:
    - Central compartment: receives IV bolus, connected to peripheral
    - Peripheral compartment: equilibrates with central
    - Elimination from central compartment

    Model equations:
        dA_central/dt = -k10*A_central - k12*A_central + k21*A_peripheral
        dA_peripheral/dt = k12*A_central - k21*A_peripheral
        C = A_central / V1

    Micro-rate constants:
        k10 = CL/V1 (elimination)
        k12 = Q/V1 (central to peripheral)
        k21 = Q/V2 (peripheral to central)

    Args:
        cl: Clearance from central compartment (volume/time)
        v1: Volume of central compartment
        q: Inter-compartmental clearance
        v2: Volume of peripheral compartment
        doses: List of dose events, each a dict with 'time' and 'amount'
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (A_central, A_peripheral)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pk_twocomp_iv_bolus(
        ...     cl=10.0, v1=50.0, q=5.0, v2=100.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=48.0,
        ...     saveat=list(range(49))
        ... )
        >>> print(f"Cmax: {max(result['observations']['conc'])}")
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    TwoCompIVBolus = jl.NeoPKPD.TwoCompIVBolus
    TwoCompIVBolusParams = jl.NeoPKPD.TwoCompIVBolusParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        TwoCompIVBolus(), "py_twocomp_iv",
        TwoCompIVBolusParams(float(cl), float(v1), float(q), float(v2)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_twocomp_oral(
    ka: float,
    cl: float,
    v1: float,
    q: float,
    v2: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a two-compartment oral first-order absorption PK simulation.

    Model structure:
    - Gut compartment: first-order absorption into central
    - Central compartment: connected to peripheral, elimination
    - Peripheral compartment: equilibrates with central

    Args:
        ka: Absorption rate constant (1/time)
        cl: Clearance from central compartment (volume/time)
        v1: Volume of central compartment
        q: Inter-compartmental clearance
        v2: Volume of peripheral compartment
        doses: List of dose events
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (A_gut, A_central, A_peripheral)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pk_twocomp_oral(
        ...     ka=1.0, cl=10.0, v1=50.0, q=5.0, v2=100.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=48.0,
        ...     saveat=list(range(49))
        ... )
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    TwoCompOral = jl.NeoPKPD.TwoCompOral
    TwoCompOralParams = jl.NeoPKPD.TwoCompOralParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        TwoCompOral(), "py_twocomp_oral",
        TwoCompOralParams(float(ka), float(cl), float(v1), float(q), float(v2)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_threecomp_iv_bolus(
    cl: float,
    v1: float,
    q2: float,
    v2: float,
    q3: float,
    v3: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a three-compartment IV bolus PK simulation (mammillary model).

    Model structure:
    - Central compartment: receives IV bolus, connected to two peripheral compartments
    - Shallow peripheral (periph1): rapid equilibration with central
    - Deep peripheral (periph2): slow equilibration with central
    - Elimination from central compartment

    The concentration-time profile shows tri-exponential decay:
    - Alpha phase: rapid initial decline (distribution to shallow peripheral)
    - Beta phase: intermediate decline (distribution to deep peripheral)
    - Gamma phase: terminal elimination phase

    Args:
        cl: Clearance from central compartment (volume/time)
        v1: Volume of central compartment
        q2: Inter-compartmental clearance to shallow peripheral
        v2: Volume of shallow peripheral compartment
        q3: Inter-compartmental clearance to deep peripheral
        v3: Volume of deep peripheral compartment
        doses: List of dose events
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (A_central, A_periph1, A_periph2)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pk_threecomp_iv_bolus(
        ...     cl=10.0, v1=50.0, q2=10.0, v2=80.0, q3=2.0, v3=200.0,
        ...     doses=[{"time": 0.0, "amount": 1000.0}],
        ...     t0=0.0, t1=72.0,
        ...     saveat=list(range(73))
        ... )
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    ThreeCompIVBolus = jl.NeoPKPD.ThreeCompIVBolus
    ThreeCompIVBolusParams = jl.NeoPKPD.ThreeCompIVBolusParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        ThreeCompIVBolus(), "py_threecomp_iv",
        ThreeCompIVBolusParams(float(cl), float(v1), float(q2), float(v2), float(q3), float(v3)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_transit_absorption(
    n: int,
    ktr: float,
    ka: float,
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a transit compartment absorption PK simulation.

    Model structure:
    - Chain of N transit compartments modeling GI transit
    - Final absorption from last transit into central compartment
    - First-order elimination from central

    This model captures delayed absorption profiles commonly seen
    with oral drugs, providing a gamma-like input function.

    Mean transit time (MTT) ≈ (N+1) / Ktr

    Args:
        n: Number of transit compartments (1-20)
        ktr: Transit rate constant between compartments (1/time)
        ka: Absorption rate constant from last transit to central (1/time)
        cl: Clearance from central compartment (volume/time)
        v: Volume of central compartment
        doses: List of dose events
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (Transit_1, ..., Transit_N, A_central)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata including N_transit

    Example:
        >>> result = neopkpd.simulate_pk_transit_absorption(
        ...     n=5, ktr=2.0, ka=1.0, cl=10.0, v=50.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=24.0,
        ...     saveat=[x * 0.25 for x in range(97)]  # every 15 min
        ... )
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    TransitAbsorption = jl.NeoPKPD.TransitAbsorption
    TransitAbsorptionParams = jl.NeoPKPD.TransitAbsorptionParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        TransitAbsorption(), "py_transit",
        TransitAbsorptionParams(int(n), float(ktr), float(ka), float(cl), float(v)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_michaelis_menten(
    vmax: float,
    km: float,
    v: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a one-compartment PK simulation with Michaelis-Menten (saturable) elimination.

    Model equation:
        dA/dt = -Vmax * C / (Km + C)
        C = A / V

    This model describes nonlinear pharmacokinetics where elimination
    becomes saturated at high concentrations. Key characteristics:
    - At low C (C << Km): Approximately linear elimination, CL ≈ Vmax/Km
    - At high C (C >> Km): Zero-order elimination, rate ≈ Vmax
    - Dose-dependent half-life (increases with dose)
    - Disproportionate increase in AUC with dose

    Common examples: Phenytoin, Ethanol, high-dose Aspirin

    Args:
        vmax: Maximum elimination rate (amount/time)
        km: Michaelis constant (concentration at half-Vmax)
        v: Volume of distribution
        doses: List of dose events
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (A_central)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pk_michaelis_menten(
        ...     vmax=100.0, km=5.0, v=50.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=48.0,
        ...     saveat=list(range(49))
        ... )
        >>> # Note: half-life varies with concentration
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    MichaelisMentenElimination = jl.NeoPKPD.MichaelisMentenElimination
    MichaelisMentenEliminationParams = jl.NeoPKPD.MichaelisMentenEliminationParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        MichaelisMentenElimination(), "py_mm",
        MichaelisMentenEliminationParams(float(vmax), float(km), float(v)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


# ============================================================================
# Coupled PKPD Simulation
# ============================================================================

def simulate_pkpd_direct_emax(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    e0: float,
    emax: float,
    ec50: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with direct Emax effect model.

    PD equation:
        Effect = E0 + (Emax * C) / (EC50 + C)

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events
        e0: Baseline effect
        emax: Maximum effect
        ec50: Concentration at 50% of Emax
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type ("OneCompIVBolus" or "OneCompOralFirstOrder")
        ka: Absorption rate constant (required if pk_kind is oral)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables
        - observations: Dict of observables (conc, effect)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_direct_emax(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     e0=0.0, emax=100.0, ec50=5.0,
        ...     t0=0.0, t1=24.0,
        ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        ... )
        >>> print(f"Max effect: {max(result['observations']['effect'])}")
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    DirectEmax = jl.NeoPKPD.DirectEmax
    DirectEmaxParams = jl.NeoPKPD.DirectEmaxParams

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    if pk_kind == "OneCompIVBolus":
        OneCompIVBolus = jl.NeoPKPD.OneCompIVBolus
        OneCompIVBolusParams = jl.NeoPKPD.OneCompIVBolusParams
        pk_spec = ModelSpec(OneCompIVBolus(), "py_pkpd_pk", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    elif pk_kind == "OneCompOralFirstOrder":
        if ka is None:
            raise ValueError("ka required for OneCompOralFirstOrder PK model")
        OneCompOralFirstOrder = jl.NeoPKPD.OneCompOralFirstOrder
        OneCompOralFirstOrderParams = jl.NeoPKPD.OneCompOralFirstOrderParams
        pk_spec = ModelSpec(
            OneCompOralFirstOrder(), "py_pkpd_pk",
            OneCompOralFirstOrderParams(float(ka), float(cl), float(v)), doses_vec
        )
    else:
        raise ValueError(f"Unsupported pk_kind: {pk_kind}")

    pd_spec = PDSpec(
        DirectEmax(), "py_pkpd_pd",
        DirectEmaxParams(float(e0), float(emax), float(ec50)),
        jl.Symbol("conc"), jl.Symbol("effect")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pkpd_indirect_response(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    kin: float,
    kout: float,
    r0: float,
    imax: float,
    ic50: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with indirect response turnover model.

    PD equations:
        I(C) = (Imax * C) / (IC50 + C)
        dR/dt = Kin - Kout * (1 - I(C)) * R

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events
        kin: Zero-order production rate
        kout: First-order elimination rate
        r0: Baseline response (should equal Kin/Kout for steady state)
        imax: Maximum inhibition (0 to 1)
        ic50: Concentration at 50% of Imax
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type ("OneCompIVBolus" or "OneCompOralFirstOrder")
        ka: Absorption rate constant (required if pk_kind is oral)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (including R for response)
        - observations: Dict of observables (conc, response)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_indirect_response(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     kin=10.0, kout=0.1, r0=100.0, imax=0.9, ic50=5.0,
        ...     t0=0.0, t1=120.0,
        ...     saveat=list(range(121))
        ... )
        >>> print(f"Min response: {min(result['observations']['response'])}")
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    IndirectResponseTurnover = jl.NeoPKPD.IndirectResponseTurnover
    IndirectResponseTurnoverParams = jl.NeoPKPD.IndirectResponseTurnoverParams

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    if pk_kind == "OneCompIVBolus":
        OneCompIVBolus = jl.NeoPKPD.OneCompIVBolus
        OneCompIVBolusParams = jl.NeoPKPD.OneCompIVBolusParams
        pk_spec = ModelSpec(OneCompIVBolus(), "py_pkpd_pk", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    elif pk_kind == "OneCompOralFirstOrder":
        if ka is None:
            raise ValueError("ka required for OneCompOralFirstOrder PK model")
        OneCompOralFirstOrder = jl.NeoPKPD.OneCompOralFirstOrder
        OneCompOralFirstOrderParams = jl.NeoPKPD.OneCompOralFirstOrderParams
        pk_spec = ModelSpec(
            OneCompOralFirstOrder(), "py_pkpd_pk",
            OneCompOralFirstOrderParams(float(ka), float(cl), float(v)), doses_vec
        )
    else:
        raise ValueError(f"Unsupported pk_kind: {pk_kind}")

    pd_spec = PDSpec(
        IndirectResponseTurnover(), "py_pkpd_pd",
        IndirectResponseTurnoverParams(float(kin), float(kout), float(r0), float(imax), float(ic50)),
        jl.Symbol("conc"), jl.Symbol("response")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pkpd_sigmoid_emax(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with Sigmoid Emax (Hill equation) effect model.

    PD equation (Hill equation):
        Effect = E0 + (Emax * C^gamma) / (EC50^gamma + C^gamma)

    The gamma (Hill coefficient) controls the steepness of the response:
    - gamma = 1: Standard hyperbolic Emax model
    - gamma > 1: Steeper, more "switch-like" response (threshold effect)
    - gamma < 1: More gradual response

    Common gamma values: 0.5-5 for most drugs, 3-6 for neuromuscular blockers

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events
        e0: Baseline effect
        emax: Maximum effect (can be negative for inhibitory effects)
        ec50: Concentration at 50% of Emax
        gamma: Hill coefficient (steepness parameter)
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type ("OneCompIVBolus" or "OneCompOralFirstOrder")
        ka: Absorption rate constant (required if pk_kind is oral)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables
        - observations: Dict of observables (conc, effect)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_sigmoid_emax(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     e0=0.0, emax=100.0, ec50=5.0, gamma=2.0,
        ...     t0=0.0, t1=24.0,
        ...     saveat=list(range(25))
        ... )
        >>> # Effect curve will be steeper than standard Emax
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    SigmoidEmax = jl.NeoPKPD.SigmoidEmax
    SigmoidEmaxParams = jl.NeoPKPD.SigmoidEmaxParams

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    if pk_kind == "OneCompIVBolus":
        OneCompIVBolus = jl.NeoPKPD.OneCompIVBolus
        OneCompIVBolusParams = jl.NeoPKPD.OneCompIVBolusParams
        pk_spec = ModelSpec(OneCompIVBolus(), "py_pkpd_pk", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    elif pk_kind == "OneCompOralFirstOrder":
        if ka is None:
            raise ValueError("ka required for OneCompOralFirstOrder PK model")
        OneCompOralFirstOrder = jl.NeoPKPD.OneCompOralFirstOrder
        OneCompOralFirstOrderParams = jl.NeoPKPD.OneCompOralFirstOrderParams
        pk_spec = ModelSpec(
            OneCompOralFirstOrder(), "py_pkpd_pk",
            OneCompOralFirstOrderParams(float(ka), float(cl), float(v)), doses_vec
        )
    else:
        raise ValueError(f"Unsupported pk_kind: {pk_kind}")

    pd_spec = PDSpec(
        SigmoidEmax(), "py_pkpd_pd",
        SigmoidEmaxParams(float(e0), float(emax), float(ec50), float(gamma)),
        jl.Symbol("conc"), jl.Symbol("effect")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pkpd_biophase_equilibration(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    ke0: float,
    e0: float,
    emax: float,
    ec50: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with biophase equilibration (effect compartment) model.

    This model introduces a hypothetical effect site compartment to account for
    temporal delay between plasma concentration changes and pharmacodynamic effects.

    PD equations:
        dCe/dt = ke0 * (Cp - Ce)  # Effect compartment equilibration
        Effect = E0 + (Emax * Ce) / (EC50 + Ce)  # Emax effect from effect site

    The equilibration half-life t1/2,ke0 = ln(2)/ke0 indicates:
    - Small ke0 (long t1/2): Slow equilibration, significant hysteresis
    - Large ke0 (short t1/2): Fast equilibration, near-direct effect

    Common applications: Anesthetics, CNS-active drugs, neuromuscular blockers

    Note: This function uses quasi-steady state approximation (Ce ≈ Cp).
    For true effect compartment dynamics with hysteresis, use the full ODE approach.

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events
        ke0: Effect site equilibration rate constant (1/time)
        e0: Baseline effect
        emax: Maximum effect
        ec50: Effect site concentration at 50% of Emax
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type ("OneCompIVBolus" or "OneCompOralFirstOrder")
        ka: Absorption rate constant (required if pk_kind is oral)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables
        - observations: Dict of observables (conc, effect)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_biophase_equilibration(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     ke0=0.5, e0=0.0, emax=100.0, ec50=5.0,
        ...     t0=0.0, t1=24.0,
        ...     saveat=list(range(25))
        ... )
        >>> # Effect will lag behind concentration due to equilibration delay
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    BiophaseEquilibration = jl.NeoPKPD.BiophaseEquilibration
    BiophaseEquilibrationParams = jl.NeoPKPD.BiophaseEquilibrationParams

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    if pk_kind == "OneCompIVBolus":
        OneCompIVBolus = jl.NeoPKPD.OneCompIVBolus
        OneCompIVBolusParams = jl.NeoPKPD.OneCompIVBolusParams
        pk_spec = ModelSpec(OneCompIVBolus(), "py_pkpd_pk", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    elif pk_kind == "OneCompOralFirstOrder":
        if ka is None:
            raise ValueError("ka required for OneCompOralFirstOrder PK model")
        OneCompOralFirstOrder = jl.NeoPKPD.OneCompOralFirstOrder
        OneCompOralFirstOrderParams = jl.NeoPKPD.OneCompOralFirstOrderParams
        pk_spec = ModelSpec(
            OneCompOralFirstOrder(), "py_pkpd_pk",
            OneCompOralFirstOrderParams(float(ka), float(cl), float(v)), doses_vec
        )
    else:
        raise ValueError(f"Unsupported pk_kind: {pk_kind}")

    pd_spec = PDSpec(
        BiophaseEquilibration(), "py_pkpd_pd",
        BiophaseEquilibrationParams(float(ke0), float(e0), float(emax), float(ec50)),
        jl.Symbol("conc"), jl.Symbol("effect")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


# ============================================================================
# Population Simulation
# ============================================================================

def simulate_population_iv_bolus(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    n: int,
    seed: int,
    omegas: Dict[str, float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a population PK simulation with inter-individual variability (IIV).

    Uses log-normal IIV:
        theta_i = theta_pop * exp(eta_i)
        eta_i ~ Normal(0, omega^2)

    Args:
        cl: Population typical clearance
        v: Population typical volume
        doses: List of dose events
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        n: Number of individuals
        seed: Random seed for reproducibility
        omegas: Dict of omega values for each parameter (e.g., {"CL": 0.3, "V": 0.2})
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - individuals: List of individual simulation results
        - params: List of realized parameters for each individual
        - summaries: Dict of summary statistics (mean, median, quantiles)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_population_iv_bolus(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     t0=0.0, t1=24.0,
        ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
        ...     n=100, seed=12345,
        ...     omegas={"CL": 0.3, "V": 0.2}
        ... )
        >>> print(f"Number of individuals: {len(result['individuals'])}")
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    OneCompIVBolus = jl.NeoPKPD.OneCompIVBolus
    OneCompIVBolusParams = jl.NeoPKPD.OneCompIVBolusParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    IIVSpec = jl.NeoPKPD.IIVSpec
    LogNormalIIV = jl.NeoPKPD.LogNormalIIV
    PopulationSpec = jl.NeoPKPD.PopulationSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    base = ModelSpec(OneCompIVBolus(), "py_pop_iv", OneCompIVBolusParams(float(cl), float(v)), doses_vec)

    omega_j = {jl.Symbol(k): float(val) for (k, val) in omegas.items()}
    iiv = IIVSpec(LogNormalIIV(), omega_j, jl.UInt64(int(seed)), int(n))

    # No IOV, no covariate model, no covariates
    pop = PopulationSpec(base, iiv, None, None, [])

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_population(pop, grid, solver)
    return _popresult_to_py(res)


def simulate_population_oral(
    ka: float,
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    n: int,
    seed: int,
    omegas: Dict[str, float],
    iov_pis: Optional[Dict[str, float]] = None,
    iov_seed: Optional[int] = None,
    covariates: Optional[List[Dict[str, float]]] = None,
    covariate_effects: Optional[List[Dict[str, Any]]] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a population oral PK simulation with full IIV, IOV, and covariate support.

    Args:
        ka: Population typical absorption rate
        cl: Population typical clearance
        v: Population typical volume
        doses: List of dose events
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        n: Number of individuals
        seed: Random seed for IIV
        omegas: Dict of omega values for IIV (e.g., {"CL": 0.3, "V": 0.2})
        iov_pis: Optional dict of pi values for IOV (e.g., {"CL": 0.1})
        iov_seed: Optional separate seed for IOV (required if iov_pis provided)
        covariates: Optional list of covariate dicts, one per individual
                    e.g., [{"WT": 70.0}, {"WT": 85.0}, ...]
        covariate_effects: Optional list of covariate effect specifications
                    e.g., [{"kind": "PowerCovariate", "param": "CL", "covariate": "WT", "beta": 0.75, "ref": 70.0}]
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - individuals: List of individual simulation results
        - params: List of realized parameters for each individual
        - summaries: Dict of summary statistics
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_population_oral(
        ...     ka=0.5, cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     t0=0.0, t1=24.0,
        ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
        ...     n=100, seed=12345,
        ...     omegas={"CL": 0.3, "V": 0.2, "Ka": 0.3}
        ... )
    """
    jl = _require_julia()

    DoseEvent = jl.NeoPKPD.DoseEvent
    ModelSpec = jl.NeoPKPD.ModelSpec
    OneCompOralFirstOrder = jl.NeoPKPD.OneCompOralFirstOrder
    OneCompOralFirstOrderParams = jl.NeoPKPD.OneCompOralFirstOrderParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    IIVSpec = jl.NeoPKPD.IIVSpec
    IOVSpec = jl.NeoPKPD.IOVSpec
    LogNormalIIV = jl.NeoPKPD.LogNormalIIV
    OccasionDefinition = jl.NeoPKPD.OccasionDefinition
    PopulationSpec = jl.NeoPKPD.PopulationSpec
    IndividualCovariates = jl.NeoPKPD.IndividualCovariates
    CovariateModel = jl.NeoPKPD.CovariateModel
    CovariateEffect = jl.NeoPKPD.CovariateEffect
    LinearCovariate = jl.NeoPKPD.LinearCovariate
    PowerCovariate = jl.NeoPKPD.PowerCovariate
    ExpCovariate = jl.NeoPKPD.ExpCovariate

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    base = ModelSpec(
        OneCompOralFirstOrder(), "py_pop_oral",
        OneCompOralFirstOrderParams(float(ka), float(cl), float(v)),
        doses_vec
    )

    omega_j = {jl.Symbol(k): float(val) for (k, val) in omegas.items()}
    iiv = IIVSpec(LogNormalIIV(), omega_j, jl.UInt64(int(seed)), int(n))

    # IOV
    iov = None
    if iov_pis is not None:
        if iov_seed is None:
            raise ValueError("iov_seed required when iov_pis is provided")
        pi_j = {jl.Symbol(k): float(val) for (k, val) in iov_pis.items()}
        iov = IOVSpec(LogNormalIIV(), pi_j, jl.UInt64(int(iov_seed)), OccasionDefinition(jl.Symbol("dose_times")))

    # Covariate model
    cov_model = None
    if covariate_effects is not None and len(covariate_effects) > 0:
        effects = []
        for eff in covariate_effects:
            kind_str = eff["kind"]
            if kind_str == "LinearCovariate":
                kind = LinearCovariate()
            elif kind_str == "PowerCovariate":
                kind = PowerCovariate()
            elif kind_str == "ExpCovariate":
                kind = ExpCovariate()
            else:
                raise ValueError(f"Unknown covariate kind: {kind_str}")
            effects.append(CovariateEffect(
                kind,
                jl.Symbol(eff["param"]),
                jl.Symbol(eff["covariate"]),
                float(eff["beta"]),
                float(eff["ref"]),
            ))
        cov_model = CovariateModel("py_cov_model", effects)

    # Individual covariates
    indiv_covs = []
    if covariates is not None:
        for cov_dict in covariates:
            vals = {jl.Symbol(k): float(v) for (k, v) in cov_dict.items()}
            indiv_covs.append(IndividualCovariates(vals, None))

    pop = PopulationSpec(base, iiv, iov, cov_model, indiv_covs)

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_population(pop, grid, solver)
    return _popresult_to_py(res)


# ============================================================================
# Sensitivity Analysis
# ============================================================================

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


# ============================================================================
# PK/PD Metrics
# ============================================================================

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
    import math
    return math.log(2) * v / cl
