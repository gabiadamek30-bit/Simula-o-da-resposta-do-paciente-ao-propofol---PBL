"""
Two-Compartment PK Simulations

This module provides simulation functions for two-compartment PK models:
- IV bolus (TwoCompIVBolus)
- IV infusion (TwoCompIVBolus with duration > 0)
- Oral first-order absorption (TwoCompOral)
"""

from typing import Any, Dict, List, Optional, Union

from .._core import _require_julia, _simresult_to_py, _create_dose_events


def simulate_pk_twocomp_iv_bolus(
    cl: float,
    v1: float,
    q: float,
    v2: float,
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
    Run a two-compartment IV bolus or infusion PK simulation.

    Supports both IV bolus (instantaneous) and IV infusion (zero-order input)
    administration by specifying optional 'duration' in dose events.

    Model structure:
    - Central compartment: receives IV bolus/infusion, connected to peripheral
    - Peripheral compartment: equilibrates with central
    - Elimination from central compartment

    Model equations:
        For bolus: dA_central/dt = -k10*A_central - k12*A_central + k21*A_peripheral
        For infusion: dA_central/dt = R - k10*A_central - k12*A_central + k21*A_peripheral
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
        - states: Dict of state variables (A_central, A_peripheral)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> # IV bolus
        >>> result = neopkpd.simulate_pk_twocomp_iv_bolus(
        ...     cl=10.0, v1=50.0, q=5.0, v2=100.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=48.0,
        ...     saveat=list(range(49))
        ... )
        >>>
        >>> # IV infusion (500 mg over 2 hours)
        >>> result = neopkpd.simulate_pk_twocomp_iv_bolus(
        ...     cl=10.0, v1=50.0, q=5.0, v2=100.0,
        ...     doses=[{"time": 0.0, "amount": 500.0, "duration": 2.0}],
        ...     t0=0.0, t1=48.0,
        ...     saveat=list(range(49))
        ... )
    """
    jl = _require_julia()

    ModelSpec = jl.NeoPKPD.ModelSpec
    TwoCompIVBolus = jl.NeoPKPD.TwoCompIVBolus
    TwoCompIVBolusParams = jl.NeoPKPD.TwoCompIVBolusParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    # Create dose events with optional duration support
    doses_vec = _create_dose_events(jl, doses)

    spec = ModelSpec(
        TwoCompIVBolus(), "py_twocomp_iv",
        TwoCompIVBolusParams(float(cl), float(v1), float(q), float(v2)),
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
            TwoCompIVBolus(), "py_twocomp_iv",
            TwoCompIVBolusParams(float(cl), float(v1), float(q), float(v2)),
            doses_vec, modifiers
        )
        res = jl.NeoPKPD.simulate(spec_with_mod, grid, solver)
    else:
        res = jl.NeoPKPD.simulate(spec, grid, solver)

    return _simresult_to_py(res)


def simulate_pk_twocomp_oral(
    ka: float,
    cl: float,
    v1: float,
    q: float,
    v2: float,
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
    Run a two-compartment oral first-order absorption PK simulation.

    Model structure:
    - Gut compartment: first-order absorption into central
    - Central compartment: connected to peripheral, elimination
    - Peripheral compartment: equilibrates with central

    Model equations:
        dA_gut/dt = -Ka * A_gut
        dA_central/dt = Ka * A_gut - (CL/V1)*A_central - (Q/V1)*A_central + (Q/V2)*A_peripheral
        dA_peripheral/dt = (Q/V1)*A_central - (Q/V2)*A_peripheral
        C = A_central / V1

    Args:
        ka: Absorption rate constant (1/time)
        cl: Clearance from central compartment (volume/time)
        v1: Volume of central compartment
        q: Inter-compartmental clearance
        v2: Volume of peripheral compartment
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
        >>>
        >>> # With absorption lag time (30 min) and 80% bioavailability
        >>> result = neopkpd.simulate_pk_twocomp_oral(
        ...     ka=1.0, cl=10.0, v1=50.0, q=5.0, v2=100.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=48.0,
        ...     saveat=list(range(49)),
        ...     alag=0.5, bioavailability=0.8
        ... )
    """
    jl = _require_julia()

    ModelSpec = jl.NeoPKPD.ModelSpec
    TwoCompOral = jl.NeoPKPD.TwoCompOral
    TwoCompOralParams = jl.NeoPKPD.TwoCompOralParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    # Create dose events with optional duration support
    doses_vec = _create_dose_events(jl, doses)

    spec = ModelSpec(
        TwoCompOral(), "py_twocomp_oral",
        TwoCompOralParams(float(ka), float(cl), float(v1), float(q), float(v2)),
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
            TwoCompOral(), "py_twocomp_oral",
            TwoCompOralParams(float(ka), float(cl), float(v1), float(q), float(v2)),
            doses_vec, modifiers
        )
        res = jl.NeoPKPD.simulate(spec_with_mod, grid, solver)
    else:
        res = jl.NeoPKPD.simulate(spec, grid, solver)

    return _simresult_to_py(res)
