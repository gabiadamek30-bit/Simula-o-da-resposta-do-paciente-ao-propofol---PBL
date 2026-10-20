"""
Advanced PK Simulations

This module provides simulation functions for advanced PK models:
- Transit compartment absorption (TransitAbsorption)
- Michaelis-Menten elimination (MichaelisMentenElimination)
"""

from typing import Any, Dict, List, Optional, Union

from .._core import _require_julia, _simresult_to_py, _create_dose_events


def simulate_pk_transit_absorption(
    n: int,
    ktr: float,
    ka: float,
    cl: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
    alag: Optional[float] = None,
    bioavailability: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Run a transit compartment absorption PK simulation.

    Model structure:
    - Chain of N transit compartments modeling GI transit
    - Final absorption from last transit into central compartment
    - First-order elimination from central

    This model captures delayed absorption profiles commonly seen
    with oral drugs, providing a gamma-like input function.

    Mean transit time (MTT) = (N+1) / Ktr

    Args:
        n: Number of transit compartments (1-20)
        ktr: Transit rate constant between compartments (1/time)
        ka: Absorption rate constant from last transit to central (1/time)
        cl: Clearance from central compartment (volume/time)
        v: Volume of central compartment
        doses: List of dose events, each a dict with:
            - 'time': Dose administration time
            - 'amount': Total drug amount
            - 'duration' (optional): For extended-release formulations
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)
        alag: Absorption lag time - delays dose by this amount (default: None)
        bioavailability: Fraction of dose absorbed (0-1, default: None = 1.0)

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

    ModelSpec = jl.NeoPKPD.ModelSpec
    TransitAbsorption = jl.NeoPKPD.TransitAbsorption
    TransitAbsorptionParams = jl.NeoPKPD.TransitAbsorptionParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    # Create dose events with optional duration support
    doses_vec = _create_dose_events(jl, doses)

    spec = ModelSpec(
        TransitAbsorption(), "py_transit",
        TransitAbsorptionParams(int(n), float(ktr), float(ka), float(cl), float(v)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    # Apply dose modifiers if specified
    if alag is not None or bioavailability is not None:
        AbsorptionModifiers = jl.NeoPKPD.AbsorptionModifiers
        ModelSpecWithModifiers = jl.NeoPKPD.ModelSpecWithModifiers
        modifiers = AbsorptionModifiers(
            alag=float(alag) if alag is not None else 0.0,
            bioavailability=float(bioavailability) if bioavailability is not None else 1.0
        )
        spec_with_mod = ModelSpecWithModifiers(
            TransitAbsorption(), "py_transit",
            TransitAbsorptionParams(int(n), float(ktr), float(ka), float(cl), float(v)),
            doses_vec, modifiers
        )
        res = jl.NeoPKPD.simulate(spec_with_mod, grid, solver)
    else:
        res = jl.NeoPKPD.simulate(spec, grid, solver)

    return _simresult_to_py(res)


def simulate_pk_michaelis_menten(
    vmax: float,
    km: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
    alag: Optional[float] = None,
    bioavailability: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Run a one-compartment PK simulation with Michaelis-Menten (saturable) elimination.

    Model equation:
        dA/dt = -Vmax * C / (Km + C)
        C = A / V

    This model describes nonlinear pharmacokinetics where elimination
    becomes saturated at high concentrations. Key characteristics:
    - At low C (C << Km): Approximately linear elimination, CL = Vmax/Km
    - At high C (C >> Km): Zero-order elimination, rate = Vmax
    - Dose-dependent half-life (increases with dose)
    - Disproportionate increase in AUC with dose

    Common examples: Phenytoin, Ethanol, high-dose Aspirin

    Args:
        vmax: Maximum elimination rate (amount/time)
        km: Michaelis constant (concentration at half-Vmax)
        v: Volume of distribution
        doses: List of dose events, each a dict with:
            - 'time': Dose administration time
            - 'amount': Total drug amount
            - 'duration' (optional): Infusion duration (0 = bolus, default)
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)
        alag: Absorption lag time - delays dose by this amount (default: None)
        bioavailability: Fraction of dose absorbed (0-1, default: None = 1.0)

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

    ModelSpec = jl.NeoPKPD.ModelSpec
    MichaelisMentenElimination = jl.NeoPKPD.MichaelisMentenElimination
    MichaelisMentenEliminationParams = jl.NeoPKPD.MichaelisMentenEliminationParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    # Create dose events with optional duration support
    doses_vec = _create_dose_events(jl, doses)

    spec = ModelSpec(
        MichaelisMentenElimination(), "py_mm",
        MichaelisMentenEliminationParams(float(vmax), float(km), float(v)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    # Apply dose modifiers if specified
    if alag is not None or bioavailability is not None:
        AbsorptionModifiers = jl.NeoPKPD.AbsorptionModifiers
        ModelSpecWithModifiers = jl.NeoPKPD.ModelSpecWithModifiers
        modifiers = AbsorptionModifiers(
            alag=float(alag) if alag is not None else 0.0,
            bioavailability=float(bioavailability) if bioavailability is not None else 1.0
        )
        spec_with_mod = ModelSpecWithModifiers(
            MichaelisMentenElimination(), "py_mm",
            MichaelisMentenEliminationParams(float(vmax), float(km), float(v)),
            doses_vec, modifiers
        )
        res = jl.NeoPKPD.simulate(spec_with_mod, grid, solver)
    else:
        res = jl.NeoPKPD.simulate(spec, grid, solver)

    return _simresult_to_py(res)
