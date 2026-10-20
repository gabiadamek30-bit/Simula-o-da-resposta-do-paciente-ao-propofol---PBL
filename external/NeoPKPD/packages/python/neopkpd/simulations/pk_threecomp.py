"""
Three-Compartment PK Simulations

This module provides simulation functions for three-compartment PK models:
- IV bolus (ThreeCompIVBolus)
- IV infusion (ThreeCompIVBolus with duration > 0)
"""

from typing import Any, Dict, List, Optional, Union

from .._core import _require_julia, _simresult_to_py, _create_dose_events


def simulate_pk_threecomp_iv_bolus(
    cl: float,
    v1: float,
    q2: float,
    v2: float,
    q3: float,
    v3: float,
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
    Run a three-compartment IV bolus or infusion PK simulation (mammillary model).

    Supports both IV bolus (instantaneous) and IV infusion (zero-order input)
    administration by specifying optional 'duration' in dose events.

    Model structure:
    - Central compartment: receives IV bolus/infusion, connected to two peripheral compartments
    - Shallow peripheral (periph1): rapid equilibration with central
    - Deep peripheral (periph2): slow equilibration with central
    - Elimination from central compartment

    Model equations:
        dA_central/dt = [R] - k10*A_central - k12*A_central + k21*A_periph1 - k13*A_central + k31*A_periph2
        dA_periph1/dt = k12*A_central - k21*A_periph1
        dA_periph2/dt = k13*A_central - k31*A_periph2
        C = A_central / V1

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
        - states: Dict of state variables (A_central, A_periph1, A_periph2)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> # IV bolus
        >>> result = neopkpd.simulate_pk_threecomp_iv_bolus(
        ...     cl=10.0, v1=50.0, q2=10.0, v2=80.0, q3=2.0, v3=200.0,
        ...     doses=[{"time": 0.0, "amount": 1000.0}],
        ...     t0=0.0, t1=72.0,
        ...     saveat=list(range(73))
        ... )
        >>>
        >>> # IV infusion (1000 mg over 4 hours)
        >>> result = neopkpd.simulate_pk_threecomp_iv_bolus(
        ...     cl=10.0, v1=50.0, q2=10.0, v2=80.0, q3=2.0, v3=200.0,
        ...     doses=[{"time": 0.0, "amount": 1000.0, "duration": 4.0}],
        ...     t0=0.0, t1=72.0,
        ...     saveat=list(range(73))
        ... )
    """
    jl = _require_julia()

    ModelSpec = jl.NeoPKPD.ModelSpec
    ThreeCompIVBolus = jl.NeoPKPD.ThreeCompIVBolus
    ThreeCompIVBolusParams = jl.NeoPKPD.ThreeCompIVBolusParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    # Create dose events with optional duration support
    doses_vec = _create_dose_events(jl, doses)

    spec = ModelSpec(
        ThreeCompIVBolus(), "py_threecomp_iv",
        ThreeCompIVBolusParams(float(cl), float(v1), float(q2), float(v2), float(q3), float(v3)),
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
            ThreeCompIVBolus(), "py_threecomp_iv",
            ThreeCompIVBolusParams(float(cl), float(v1), float(q2), float(v2), float(q3), float(v3)),
            doses_vec, modifiers
        )
        res = jl.NeoPKPD.simulate(spec_with_mod, grid, solver)
    else:
        res = jl.NeoPKPD.simulate(spec, grid, solver)

    return _simresult_to_py(res)
