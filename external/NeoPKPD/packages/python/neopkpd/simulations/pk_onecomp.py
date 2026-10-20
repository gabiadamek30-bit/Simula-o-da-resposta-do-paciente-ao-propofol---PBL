"""
One-Compartment PK Simulations

This module provides simulation functions for one-compartment PK models:
- IV bolus (OneCompIVBolus)
- IV infusion (OneCompIVBolus with duration > 0)
- Oral first-order absorption (OneCompOralFirstOrder)
"""

from typing import Any, Dict, List, Optional, Union

from .._core import _require_julia, _simresult_to_py, _create_dose_events


def simulate_pk_iv_bolus(
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
    Run a one-compartment IV bolus or infusion PK simulation.

    Supports both IV bolus (instantaneous) and IV infusion (zero-order input)
    administration by specifying optional 'duration' in dose events.

    Model equations:
        For bolus: dA/dt = -k * A (with instantaneous input at dose times)
        For infusion: dA/dt = R - k * A (where R = amount/duration during infusion)
        C = A / V
        k = CL / V

    Args:
        cl: Clearance (volume/time)
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
        - states: Dict of state variables (A)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> # IV bolus
        >>> result = neopkpd.simulate_pk_iv_bolus(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     t0=0.0, t1=24.0,
        ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        ... )
        >>>
        >>> # IV infusion (100 mg over 1 hour)
        >>> result = neopkpd.simulate_pk_iv_bolus(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0, "duration": 1.0}],
        ...     t0=0.0, t1=24.0,
        ...     saveat=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        ... )
    """
    jl = _require_julia()

    ModelSpec = jl.NeoPKPD.ModelSpec
    OneCompIVBolus = jl.NeoPKPD.OneCompIVBolus
    OneCompIVBolusParams = jl.NeoPKPD.OneCompIVBolusParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    # Create dose events with optional duration support
    doses_vec = _create_dose_events(jl, doses)

    spec = ModelSpec(OneCompIVBolus(), "py_iv_bolus", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
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
            OneCompIVBolus(), "py_iv_bolus",
            OneCompIVBolusParams(float(cl), float(v)),
            doses_vec, modifiers
        )
        res = jl.NeoPKPD.simulate(spec_with_mod, grid, solver)
    else:
        res = jl.NeoPKPD.simulate(spec, grid, solver)

    return _simresult_to_py(res)


def simulate_pk_oral_first_order(
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
        >>>
        >>> # With absorption lag time (30 min) and 70% bioavailability
        >>> result = neopkpd.simulate_pk_oral_first_order(
        ...     ka=0.5, cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     t0=0.0, t1=24.0,
        ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
        ...     alag=0.5, bioavailability=0.7
        ... )
    """
    jl = _require_julia()

    ModelSpec = jl.NeoPKPD.ModelSpec
    OneCompOralFirstOrder = jl.NeoPKPD.OneCompOralFirstOrder
    OneCompOralFirstOrderParams = jl.NeoPKPD.OneCompOralFirstOrderParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    # Create dose events with optional duration support
    doses_vec = _create_dose_events(jl, doses)

    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "py_oral_first_order",
        OneCompOralFirstOrderParams(float(ka), float(cl), float(v)),
        doses_vec,
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
            OneCompOralFirstOrder(), "py_oral_first_order",
            OneCompOralFirstOrderParams(float(ka), float(cl), float(v)),
            doses_vec, modifiers
        )
        res = jl.NeoPKPD.simulate(spec_with_mod, grid, solver)
    else:
        res = jl.NeoPKPD.simulate(spec, grid, solver)

    return _simresult_to_py(res)
