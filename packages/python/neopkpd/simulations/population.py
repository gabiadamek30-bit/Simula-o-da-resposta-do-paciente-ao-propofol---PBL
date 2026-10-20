"""
Population Simulations

This module provides simulation functions for population PK simulations
with inter-individual variability (IIV), inter-occasion variability (IOV),
and covariate effects.
"""

from typing import Any, Dict, List, Optional

from .._core import _require_julia, _popresult_to_py, _to_julia_vector


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
