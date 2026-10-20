"""
Custom ODE PK Simulations

This module provides simulation functions for pre-built custom PK models:
- Target-Mediated Drug Disposition (TMDD)
- Parallel first-order absorption
- Enterohepatic recirculation
- Autoinduction

These models demonstrate advanced PK mechanisms that go beyond standard
compartmental models.
"""

from typing import Any, Dict, List, Optional, Union

from .._core import _require_julia, _simresult_to_py, _create_dose_events


def simulate_pk_tmdd_custom(
    kel: float,
    kon: float,
    koff: float,
    ksyn: float,
    kdeg: float,
    kint: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a Target-Mediated Drug Disposition (TMDD) PK simulation.

    This model describes drugs that bind to their pharmacological target,
    forming a drug-target complex that affects both PK and PD.

    States:
    - L_free: Free drug (ligand) concentration
    - R_free: Free receptor concentration
    - RL_complex: Drug-receptor complex concentration

    Model equations:
        dL/dt = -kel*L - kon*L*R + koff*RL
        dR/dt = ksyn - kdeg*R - kon*L*R + koff*RL
        dRL/dt = kon*L*R - koff*RL - kint*RL

    Initial conditions:
        L(0) = dose/V (after bolus)
        R(0) = ksyn/kdeg (steady-state receptor)
        RL(0) = 0

    Args:
        kel: Drug elimination rate constant (1/time)
        kon: Association rate constant (1/(concentration*time))
        koff: Dissociation rate constant (1/time)
        ksyn: Receptor synthesis rate (concentration/time)
        kdeg: Receptor degradation rate (1/time)
        kint: Complex internalization rate (1/time)
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
        - states: Dict with L_free, R_free, RL_complex
        - observations: Dict with conc (free drug concentration)
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pk_tmdd_custom(
        ...     kel=0.1, kon=0.01, koff=0.001,
        ...     ksyn=1.0, kdeg=0.1, kint=0.05, v=50.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=72.0,
        ...     saveat=list(range(73))
        ... )
    """
    jl = _require_julia()

    ModelSpec = jl.NeoPKPD.ModelSpec
    CustomODE = jl.NeoPKPD.target_mediated_drug_disposition()
    CustomODEParams = jl.NeoPKPD.CustomODEParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    doses_vec = _create_dose_events(jl, doses)

    # Create parameter dict for CustomODEParams
    params = CustomODEParams(
        kel=float(kel),
        kon=float(kon),
        koff=float(koff),
        ksyn=float(ksyn),
        kdeg=float(kdeg),
        kint=float(kint),
        V=float(v),
    )

    spec = ModelSpec(CustomODE, "py_tmdd_custom", params, doses_vec)
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_parallel_absorption(
    ka1: float,
    ka2: float,
    f1: float,
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
) -> Dict[str, Any]:
    """
    Run a parallel first-order absorption PK simulation.

    This model describes drugs with multiple absorption sites or mechanisms,
    each with its own absorption rate constant.

    States:
    - A_depot1: Amount in first depot
    - A_depot2: Amount in second depot
    - A_central: Amount in central compartment

    Model equations:
        dA1/dt = -Ka1 * A1
        dA2/dt = -Ka2 * A2
        dAc/dt = Ka1*A1 + Ka2*A2 - (CL/V)*Ac
        C = Ac / V

    Dose splitting: F1 fraction goes to depot1, (1-F1) to depot2

    Args:
        ka1: First absorption rate constant (1/time)
        ka2: Second absorption rate constant (1/time)
        f1: Fraction to first depot (0-1, F2 = 1-F1)
        cl: Clearance (volume/time)
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
        - states: Dict with A_depot1, A_depot2, A_central
        - observations: Dict with conc
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pk_parallel_absorption(
        ...     ka1=2.0, ka2=0.5, f1=0.6, cl=10.0, v=50.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=24.0,
        ...     saveat=list(range(25))
        ... )
    """
    jl = _require_julia()

    ModelSpec = jl.NeoPKPD.ModelSpec
    CustomODE = jl.NeoPKPD.parallel_first_order_absorption()
    CustomODEParams = jl.NeoPKPD.CustomODEParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    doses_vec = _create_dose_events(jl, doses)

    params = CustomODEParams(
        Ka1=float(ka1),
        Ka2=float(ka2),
        F1=float(f1),
        CL=float(cl),
        V=float(v),
    )

    spec = ModelSpec(CustomODE, "py_parallel_absorption", params, doses_vec)
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_enterohepatic_recirculation(
    ka: float,
    cl: float,
    v: float,
    kbile: float,
    kreab: float,
    f_reab: float,
    doses: List[Dict[str, Union[float, int]]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run an enterohepatic recirculation (EHR) PK simulation.

    This model describes drugs that undergo biliary excretion and
    reabsorption from the GI tract, leading to secondary peaks.

    States:
    - A_gut: Amount in GI tract
    - A_central: Amount in central compartment
    - A_bile: Amount in biliary compartment

    Model equations:
        dAgut/dt = -Ka*Agut + F_reab*Kreab*Abile
        dAc/dt = Ka*Agut - (CL/V)*Ac - Kbile*Ac
        dAbile/dt = Kbile*Ac - Kreab*Abile
        C = Ac / V

    Args:
        ka: Absorption rate constant (1/time)
        cl: Clearance (volume/time)
        v: Volume of distribution
        kbile: Biliary excretion rate constant (1/time)
        kreab: Reabsorption rate constant from bile (1/time)
        f_reab: Fraction reabsorbed (0-1)
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
        - states: Dict with A_gut, A_central, A_bile
        - observations: Dict with conc
        - metadata: Dict of run metadata

    Example:
        >>> result = neopkpd.simulate_pk_enterohepatic_recirculation(
        ...     ka=1.0, cl=10.0, v=50.0,
        ...     kbile=0.5, kreab=0.3, f_reab=0.7,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=48.0,
        ...     saveat=list(range(49))
        ... )
        >>> # Look for secondary peaks in concentration
    """
    jl = _require_julia()

    ModelSpec = jl.NeoPKPD.ModelSpec
    CustomODE = jl.NeoPKPD.enterohepatic_recirculation()
    CustomODEParams = jl.NeoPKPD.CustomODEParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    doses_vec = _create_dose_events(jl, doses)

    params = CustomODEParams(
        Ka=float(ka),
        CL=float(cl),
        V=float(v),
        Kbile=float(kbile),
        Kreab=float(kreab),
        F_reab=float(f_reab),
    )

    spec = ModelSpec(CustomODE, "py_ehr", params, doses_vec)
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_autoinduction(
    cl0: float,
    v: float,
    emax: float,
    ec50: float,
    kenz: float,
    doses: List[Dict[str, Union[float, int]]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run an autoinduction PK simulation.

    This model describes drugs that induce their own metabolism,
    leading to time-varying clearance.

    States:
    - A_central: Amount in central compartment
    - E_enzyme: Enzyme level (relative to baseline, starts at 1.0)

    Model equations:
        E_induced = 1 + Emax*C/(EC50 + C)
        dE/dt = kenz*(E_induced - E)
        dAc/dt = -(CL0*E/V)*Ac
        C = Ac / V

    The enzyme turnover model causes clearance to gradually increase
    during chronic dosing, then return to baseline after drug washout.

    Args:
        cl0: Baseline clearance (volume/time)
        v: Volume of distribution
        emax: Maximum enzyme induction (fold increase)
        ec50: Concentration for 50% induction
        kenz: Enzyme turnover rate (1/time)
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
        - states: Dict with A_central, E_enzyme
        - observations: Dict with conc
        - metadata: Dict of run metadata

    Example:
        >>> # Chronic dosing with autoinduction
        >>> doses = [{"time": i*24.0, "amount": 200.0} for i in range(7)]
        >>> result = neopkpd.simulate_pk_autoinduction(
        ...     cl0=10.0, v=50.0, emax=2.0, ec50=5.0, kenz=0.1,
        ...     doses=doses,
        ...     t0=0.0, t1=168.0,
        ...     saveat=[x*0.5 for x in range(337)]
        ... )
        >>> # Note: Clearance increases over time due to enzyme induction
    """
    jl = _require_julia()

    ModelSpec = jl.NeoPKPD.ModelSpec
    CustomODE = jl.NeoPKPD.autoinduction()
    CustomODEParams = jl.NeoPKPD.CustomODEParams
    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec

    doses_vec = _create_dose_events(jl, doses)

    params = CustomODEParams(
        CL0=float(cl0),
        V=float(v),
        Emax=float(emax),
        EC50=float(ec50),
        kenz=float(kenz),
    )

    spec = ModelSpec(CustomODE, "py_autoinduction", params, doses_vec)
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate(spec, grid, solver)
    return _simresult_to_py(res)
