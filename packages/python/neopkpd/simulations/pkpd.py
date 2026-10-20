"""
PK-PD Simulations

This module provides simulation functions for coupled PK-PD models.

Basic Effect Models:
- Direct Emax (DirectEmax)
- Sigmoid Emax / Hill equation (SigmoidEmax)
- Biophase equilibration / effect compartment (BiophaseEquilibration)

Indirect Response Models:
- IRM-I: Inhibition of Kin (IndirectResponseIRM1)
- IRM-II: Stimulation of Kin (IndirectResponseIRM2)
- IRM-III: Inhibition of Kout (IndirectResponseTurnover)
- IRM-IV: Stimulation of Kout (IndirectResponseIRM4)

Advanced PD Models:
- Transit Compartment PD (TransitCompartmentPD)
- Disease Progression (DiseaseProgressionPD)

Combination Effect Models:
- Bliss Independence (BlissIndependence)
- Competitive Inhibition (CompetitiveInhibition)
- Drug Interaction / Greco model (DrugInteraction)

Tolerance Models:
- Counter-regulation tolerance (ToleranceCounterRegulation)
- Receptor regulation (ReceptorRegulation)
"""

from typing import Any, Dict, List, Optional, Union

from .._core import _require_julia, _simresult_to_py, _to_julia_vector, _create_dose_events


def _create_pk_spec(jl, pk_kind: str, cl: float, v: float, ka: Optional[float],
                    doses: List[Dict[str, Union[float, int]]],
                    q: Optional[float] = None, v2: Optional[float] = None):
    """Helper to create PK model spec based on pk_kind string."""
    ModelSpec = jl.NeoPKPD.ModelSpec

    doses_vec = _create_dose_events(jl, doses)

    if pk_kind == "OneCompIVBolus":
        Kind = jl.NeoPKPD.OneCompIVBolus
        Params = jl.NeoPKPD.OneCompIVBolusParams
        return ModelSpec(Kind(), "py_pkpd_pk", Params(float(cl), float(v)), doses_vec)

    elif pk_kind == "OneCompOralFirstOrder":
        if ka is None:
            raise ValueError("ka required for OneCompOralFirstOrder PK model")
        Kind = jl.NeoPKPD.OneCompOralFirstOrder
        Params = jl.NeoPKPD.OneCompOralFirstOrderParams
        return ModelSpec(Kind(), "py_pkpd_pk", Params(float(ka), float(cl), float(v)), doses_vec)

    elif pk_kind == "TwoCompIVBolus":
        if q is None or v2 is None:
            raise ValueError("q and v2 required for TwoCompIVBolus PK model")
        Kind = jl.NeoPKPD.TwoCompIVBolus
        Params = jl.NeoPKPD.TwoCompIVBolusParams
        return ModelSpec(Kind(), "py_pkpd_pk", Params(float(cl), float(v), float(q), float(v2)), doses_vec)

    elif pk_kind == "TwoCompOral":
        if ka is None or q is None or v2 is None:
            raise ValueError("ka, q, and v2 required for TwoCompOral PK model")
        Kind = jl.NeoPKPD.TwoCompOral
        Params = jl.NeoPKPD.TwoCompOralParams
        return ModelSpec(Kind(), "py_pkpd_pk", Params(float(ka), float(cl), float(v), float(q), float(v2)), doses_vec)

    else:
        raise ValueError(f"Unsupported pk_kind: {pk_kind}. Supported: OneCompIVBolus, OneCompOralFirstOrder, TwoCompIVBolus, TwoCompOral")


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

    # DirectEmax is algebraic, use simulate_pkpd (not coupled)
    res = jl.NeoPKPD.simulate_pkpd(pk_spec, pd_spec, grid, solver)
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

    Note: This function uses quasi-steady state approximation (Ce = Cp).
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


# =============================================================================
# Indirect Response Models (IRM-I, IRM-II, IRM-IV)
# =============================================================================

def simulate_pkpd_irm1(
    cl: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
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
    q: Optional[float] = None,
    v2: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with Indirect Response Model Type I (IRM-I).

    IRM-I: Drug INHIBITS PRODUCTION (Kin)

    PD equations:
        I(C) = (Imax * C) / (IC50 + C)  # Inhibition function
        dR/dt = Kin * (1 - I(C)) - Kout * R  # Response dynamics

    Effect: Drug decreases production rate -> Response decreases below baseline.

    Clinical Applications:
    - Corticosteroids on cortisol production
    - Statins on cholesterol synthesis
    - Immunosuppressants on cytokine production

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events with 'time', 'amount', optional 'duration'
        kin: Zero-order production rate (amount/time)
        kout: First-order elimination rate constant (1/time)
        r0: Baseline response (= Kin/Kout at steady state)
        imax: Maximum inhibition [0, 1]
        ic50: Concentration at 50% of maximum inhibition
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type (see supported models)
        ka: Absorption rate constant (for oral models)
        q: Inter-compartmental clearance (for two-compartment models)
        v2: Peripheral volume (for two-compartment models)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys: t, states (including R), observations (conc, response), metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_irm1(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     kin=10.0, kout=0.1, r0=100.0, imax=0.9, ic50=5.0,
        ...     t0=0.0, t1=120.0,
        ...     saveat=list(range(121))
        ... )
    """
    jl = _require_julia()

    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    IndirectResponseIRM1 = jl.NeoPKPD.IndirectResponseIRM1
    IndirectResponseIRM1Params = jl.NeoPKPD.IndirectResponseIRM1Params

    pk_spec = _create_pk_spec(jl, pk_kind, cl, v, ka, doses, q, v2)

    pd_spec = PDSpec(
        IndirectResponseIRM1(), "py_pkpd_pd",
        IndirectResponseIRM1Params(float(kin), float(kout), float(r0), float(imax), float(ic50)),
        jl.Symbol("conc"), jl.Symbol("response")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pkpd_irm2(
    cl: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
    kin: float,
    kout: float,
    r0: float,
    smax: float,
    sc50: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    q: Optional[float] = None,
    v2: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with Indirect Response Model Type II (IRM-II).

    IRM-II: Drug STIMULATES PRODUCTION (Kin)

    PD equations:
        S(C) = (Smax * C) / (SC50 + C)  # Stimulation function
        dR/dt = Kin * (1 + S(C)) - Kout * R  # Response dynamics

    Effect: Drug increases production rate -> Response increases above baseline.

    Clinical Applications:
    - EPO on red blood cell production
    - G-CSF on neutrophil production
    - Growth factors on tissue growth

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events with 'time', 'amount', optional 'duration'
        kin: Zero-order production rate (amount/time)
        kout: First-order elimination rate constant (1/time)
        r0: Baseline response (= Kin/Kout at steady state)
        smax: Maximum stimulation (can exceed 1)
        sc50: Concentration at 50% of maximum stimulation
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type (see supported models)
        ka: Absorption rate constant (for oral models)
        q: Inter-compartmental clearance (for two-compartment models)
        v2: Peripheral volume (for two-compartment models)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys: t, states (including R), observations (conc, response), metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_irm2(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     kin=10.0, kout=0.1, r0=100.0, smax=2.0, sc50=5.0,
        ...     t0=0.0, t1=120.0,
        ...     saveat=list(range(121))
        ... )
    """
    jl = _require_julia()

    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    IndirectResponseIRM2 = jl.NeoPKPD.IndirectResponseIRM2
    IndirectResponseIRM2Params = jl.NeoPKPD.IndirectResponseIRM2Params

    pk_spec = _create_pk_spec(jl, pk_kind, cl, v, ka, doses, q, v2)

    pd_spec = PDSpec(
        IndirectResponseIRM2(), "py_pkpd_pd",
        IndirectResponseIRM2Params(float(kin), float(kout), float(r0), float(smax), float(sc50)),
        jl.Symbol("conc"), jl.Symbol("response")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pkpd_irm4(
    cl: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
    kin: float,
    kout: float,
    r0: float,
    smax: float,
    sc50: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    q: Optional[float] = None,
    v2: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with Indirect Response Model Type IV (IRM-IV).

    IRM-IV: Drug STIMULATES ELIMINATION (Kout)

    PD equations:
        S(C) = (Smax * C) / (SC50 + C)  # Stimulation function
        dR/dt = Kin - Kout * (1 + S(C)) * R  # Response dynamics

    Effect: Drug increases elimination rate -> Response decreases below baseline.

    Clinical Applications:
    - Diuretics on sodium excretion
    - Thyroid hormone on metabolic rate
    - Laxatives on bowel motility

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events with 'time', 'amount', optional 'duration'
        kin: Zero-order production rate (amount/time)
        kout: First-order elimination rate constant (1/time)
        r0: Baseline response (= Kin/Kout at steady state)
        smax: Maximum stimulation (can exceed 1)
        sc50: Concentration at 50% of maximum stimulation
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type (see supported models)
        ka: Absorption rate constant (for oral models)
        q: Inter-compartmental clearance (for two-compartment models)
        v2: Peripheral volume (for two-compartment models)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys: t, states (including R), observations (conc, response), metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_irm4(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     kin=10.0, kout=0.1, r0=100.0, smax=2.0, sc50=5.0,
        ...     t0=0.0, t1=120.0,
        ...     saveat=list(range(121))
        ... )
    """
    jl = _require_julia()

    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    IndirectResponseIRM4 = jl.NeoPKPD.IndirectResponseIRM4
    IndirectResponseIRM4Params = jl.NeoPKPD.IndirectResponseIRM4Params

    pk_spec = _create_pk_spec(jl, pk_kind, cl, v, ka, doses, q, v2)

    pd_spec = PDSpec(
        IndirectResponseIRM4(), "py_pkpd_pd",
        IndirectResponseIRM4Params(float(kin), float(kout), float(r0), float(smax), float(sc50)),
        jl.Symbol("conc"), jl.Symbol("response")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


# =============================================================================
# Transit Compartment PD
# =============================================================================

def simulate_pkpd_transit_compartment(
    cl: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
    n_transit: int,
    ktr: float,
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    q: Optional[float] = None,
    v2: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with Transit Compartment PD model.

    This model introduces a chain of transit compartments to model delayed
    drug effects with signal transduction.

    PD equations:
        Signal(C) = E0 + Emax * C^gamma / (EC50^gamma + C^gamma)
        dA1/dt = ktr * (Signal(C) - A1)
        dAi/dt = ktr * (A(i-1) - Ai)  for i = 2..N
        Effect = AN

    Mean Transit Time (MTT) = (N + 1) / ktr

    Clinical Applications:
    - Delayed myelosuppression (neutropenia, thrombocytopenia)
    - Delayed biomarker responses
    - Signal transduction cascades

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events with 'time', 'amount', optional 'duration'
        n_transit: Number of transit compartments (1-20 typical)
        ktr: Transit rate constant (1/time)
        e0: Baseline effect/signal
        emax: Maximum effect above baseline
        ec50: Concentration at 50% of maximum effect
        gamma: Hill coefficient (steepness)
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type (see supported models)
        ka: Absorption rate constant (for oral models)
        q: Inter-compartmental clearance (for two-compartment models)
        v2: Peripheral volume (for two-compartment models)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys: t, states (transit compartments), observations (conc, effect), metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_transit_compartment(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     n_transit=5, ktr=0.5, e0=100.0, emax=-50.0, ec50=5.0, gamma=1.0,
        ...     t0=0.0, t1=168.0,  # One week
        ...     saveat=list(range(169))
        ... )
    """
    jl = _require_julia()

    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    TransitCompartmentPD = jl.NeoPKPD.TransitCompartmentPD
    TransitCompartmentPDParams = jl.NeoPKPD.TransitCompartmentPDParams

    pk_spec = _create_pk_spec(jl, pk_kind, cl, v, ka, doses, q, v2)

    pd_spec = PDSpec(
        TransitCompartmentPD(), "py_pkpd_pd",
        TransitCompartmentPDParams(int(n_transit), float(ktr), float(e0), float(emax), float(ec50), float(gamma)),
        jl.Symbol("conc"), jl.Symbol("effect")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


# =============================================================================
# Disease Progression PD
# =============================================================================

def simulate_pkpd_disease_progression(
    cl: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
    growth_model: str,
    s0: float,
    kgrow: float,
    smax: float,
    alpha: float,
    kdrug: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    q: Optional[float] = None,
    v2: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with Disease Progression model.

    Models tumor growth dynamics with drug-induced cell kill.

    Growth Models:
    - "linear": dS/dt = alpha - kdrug * C * S
    - "asymptotic": dS/dt = kgrow * (Smax - S) - kdrug * C * S
    - "gompertz": dS/dt = kgrow * S * log(Smax/S) - kdrug * C * S
    - "logistic": dS/dt = kgrow * S * (1 - S/Smax) - kdrug * C * S
    - "exponential": dS/dt = kgrow * S - kdrug * C * S

    Clinical Applications:
    - Tumor growth modeling
    - Oncology dose-response
    - Survival analysis

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events with 'time', 'amount', optional 'duration'
        growth_model: One of "linear", "asymptotic", "gompertz", "logistic", "exponential"
        s0: Initial tumor size
        kgrow: Growth rate constant
        smax: Maximum size (carrying capacity, for bounded models)
        alpha: Linear growth rate (for linear model only)
        kdrug: Drug-induced cell kill rate constant
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type (see supported models)
        ka: Absorption rate constant (for oral models)
        q: Inter-compartmental clearance (for two-compartment models)
        v2: Peripheral volume (for two-compartment models)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys: t, states (S), observations (conc, tumor_size), metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_disease_progression(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     growth_model="gompertz",
        ...     s0=100.0, kgrow=0.05, smax=1000.0, alpha=0.0, kdrug=0.01,
        ...     t0=0.0, t1=168.0,
        ...     saveat=list(range(169))
        ... )
    """
    jl = _require_julia()

    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    DiseaseProgressionPD = jl.NeoPKPD.DiseaseProgressionPD
    DiseaseProgressionPDParams = jl.NeoPKPD.DiseaseProgressionPDParams

    # Map string to Julia enum
    growth_map = {
        "linear": jl.NeoPKPD.LinearGrowth,
        "asymptotic": jl.NeoPKPD.AsymptoticGrowth,
        "gompertz": jl.NeoPKPD.GompertzGrowth,
        "logistic": jl.NeoPKPD.LogisticGrowth,
        "exponential": jl.NeoPKPD.ExponentialGrowth,
    }

    if growth_model.lower() not in growth_map:
        raise ValueError(f"Unsupported growth_model: {growth_model}. "
                        f"Supported: {list(growth_map.keys())}")

    growth_type = growth_map[growth_model.lower()]

    pk_spec = _create_pk_spec(jl, pk_kind, cl, v, ka, doses, q, v2)

    pd_spec = PDSpec(
        DiseaseProgressionPD(growth_type), "py_pkpd_pd",
        DiseaseProgressionPDParams(float(s0), float(kgrow), float(smax), float(alpha), float(kdrug)),
        jl.Symbol("conc"), jl.Symbol("tumor_size")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


# =============================================================================
# Tolerance Models
# =============================================================================

def simulate_pkpd_tolerance_counter_regulation(
    cl: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
    kin_mod: float,
    kout_mod: float,
    alpha_feedback: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    q: Optional[float] = None,
    v2: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with Tolerance Counter-Regulation model.

    Models tolerance development through a feedback mechanism that opposes
    the drug effect over time.

    PD equations:
        E_drug = Emax * C^gamma / (EC50^gamma + C^gamma)
        dM/dt = kin_mod * E_drug - kout_mod * M  # Moderator dynamics
        E_net = E0 + E_drug - alpha * M  # Net effect with tolerance

    Clinical Applications:
    - Opioid tolerance development
    - Beta-blocker tolerance
    - Benzodiazepine adaptation

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events with 'time', 'amount', optional 'duration'
        e0: Baseline effect
        emax: Maximum drug effect
        ec50: EC50 for drug effect
        gamma: Hill coefficient
        kin_mod: Moderator production rate constant
        kout_mod: Moderator elimination rate constant
        alpha_feedback: Feedback strength (moderator effect coefficient)
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type (see supported models)
        ka: Absorption rate constant (for oral models)
        q: Inter-compartmental clearance (for two-compartment models)
        v2: Peripheral volume (for two-compartment models)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys: t, states (M), observations (conc, effect), metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_tolerance_counter_regulation(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": i*8.0, "amount": 50.0} for i in range(21)],  # TID for 7 days
        ...     e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
        ...     kin_mod=0.1, kout_mod=0.05, alpha_feedback=1.0,
        ...     t0=0.0, t1=168.0,
        ...     saveat=[i*0.5 for i in range(337)]
        ... )
    """
    jl = _require_julia()

    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    ToleranceCounterRegulation = jl.NeoPKPD.ToleranceCounterRegulation
    ToleranceCounterRegulationParams = jl.NeoPKPD.ToleranceCounterRegulationParams

    pk_spec = _create_pk_spec(jl, pk_kind, cl, v, ka, doses, q, v2)

    pd_spec = PDSpec(
        ToleranceCounterRegulation(), "py_pkpd_pd",
        ToleranceCounterRegulationParams(
            float(e0), float(emax), float(ec50), float(gamma),
            float(kin_mod), float(kout_mod), float(alpha_feedback)
        ),
        jl.Symbol("conc"), jl.Symbol("effect")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pkpd_receptor_regulation(
    cl: float,
    v: float,
    doses: List[Dict[str, Union[float, int]]],
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
    r_baseline: float,
    kreg: float,
    rmax: float,
    kchange: float,
    direction: str,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    q: Optional[float] = None,
    v2: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with Receptor Regulation model.

    Models tolerance through changes in receptor density in response to
    sustained drug exposure.

    PD equations:
        E_drug = Emax * C^gamma / (EC50^gamma + C^gamma)
        dR/dt = kreg * (R_baseline - R) + regulation_effect
        E_net = E0 + R * E_drug  # Receptor amplifies/attenuates drug effect

    Direction:
    - "down": Down-regulation (regulation_effect = -kchange * E_drug * R)
    - "up": Up-regulation (regulation_effect = +kchange * E_drug * (Rmax - R))

    Clinical Applications:
    - Beta-receptor down-regulation with chronic agonist exposure
    - Opioid receptor adaptation
    - Hormone receptor regulation

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events with 'time', 'amount', optional 'duration'
        e0: Baseline effect
        emax: Maximum drug effect on receptors
        ec50: EC50 for drug effect
        gamma: Hill coefficient
        r_baseline: Baseline receptor density (normalized, typically 1.0)
        kreg: Receptor return-to-baseline rate constant
        rmax: Maximum receptor density (for up-regulation)
        kchange: Rate of receptor change
        direction: "down" for down-regulation, "up" for up-regulation
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type (see supported models)
        ka: Absorption rate constant (for oral models)
        q: Inter-compartmental clearance (for two-compartment models)
        v2: Peripheral volume (for two-compartment models)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys: t, states (R), observations (conc, effect), metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_receptor_regulation(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": i*8.0, "amount": 50.0} for i in range(21)],
        ...     e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
        ...     r_baseline=1.0, kreg=0.1, rmax=2.0, kchange=0.05, direction="down",
        ...     t0=0.0, t1=168.0,
        ...     saveat=[i*0.5 for i in range(337)]
        ... )
    """
    jl = _require_julia()

    SimGrid = jl.NeoPKPD.SimGrid
    SolverSpec = jl.NeoPKPD.SolverSpec
    PDSpec = jl.NeoPKPD.PDSpec
    ReceptorRegulation = jl.NeoPKPD.ReceptorRegulation
    ReceptorRegulationParams = jl.NeoPKPD.ReceptorRegulationParams

    if direction.lower() not in ["down", "up"]:
        raise ValueError(f"direction must be 'down' or 'up', got '{direction}'")

    pk_spec = _create_pk_spec(jl, pk_kind, cl, v, ka, doses, q, v2)

    pd_spec = PDSpec(
        ReceptorRegulation(), "py_pkpd_pd",
        ReceptorRegulationParams(
            float(e0), float(emax), float(ec50), float(gamma),
            float(r_baseline), float(kreg), float(rmax), float(kchange),
            jl.Symbol(direction.lower())
        ),
        jl.Symbol("conc"), jl.Symbol("effect")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.NeoPKPD.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


# =============================================================================
# Combination Effect Models
# =============================================================================

def simulate_pkpd_bliss_independence(
    # Drug A PK parameters
    cl_a: float,
    v_a: float,
    doses_a: List[Dict[str, Union[float, int]]],
    # Drug B PK parameters
    cl_b: float,
    v_b: float,
    doses_b: List[Dict[str, Union[float, int]]],
    # PD parameters
    e0: float,
    emax_a: float,
    ec50_a: float,
    gamma_a: float,
    emax_b: float,
    ec50_b: float,
    gamma_b: float,
    # Simulation settings
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind_a: str = "OneCompIVBolus",
    pk_kind_b: str = "OneCompIVBolus",
    ka_a: Optional[float] = None,
    ka_b: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a combination PK-PD simulation with Bliss Independence model.

    Bliss Independence assumes two drugs act independently through separate
    mechanisms. The combined effect is:
        f_combined = f_A + f_B - f_A × f_B

    where f_A and f_B are fractional effects (0-1 scale).

    PD equations:
        f_A = C_A^γ_A / (EC50_A^γ_A + C_A^γ_A)
        f_B = C_B^γ_B / (EC50_B^γ_B + C_B^γ_B)
        E = E0 + f_combined × max(Emax_A, Emax_B)

    Clinical Applications:
        - Combination chemotherapy
        - Antibacterial combinations
        - Multi-target therapies

    Args:
        cl_a: Clearance for drug A
        v_a: Volume of distribution for drug A
        doses_a: Dose events for drug A
        cl_b: Clearance for drug B
        v_b: Volume of distribution for drug B
        doses_b: Dose events for drug B
        e0: Baseline effect
        emax_a: Maximum effect for drug A
        ec50_a: EC50 for drug A
        gamma_a: Hill coefficient for drug A
        emax_b: Maximum effect for drug B
        ec50_b: EC50 for drug B
        gamma_b: Hill coefficient for drug B
        t0: Simulation start time
        t1: Simulation end time
        saveat: Time points for output
        pk_kind_a: PK model type for drug A
        pk_kind_b: PK model type for drug B
        ka_a: Absorption rate for drug A (if oral)
        ka_b: Absorption rate for drug B (if oral)
        alg: ODE solver algorithm
        reltol: Relative tolerance
        abstol: Absolute tolerance
        maxiters: Maximum solver iterations

    Returns:
        Dict with keys:
        - times: Time points
        - observations: Dict with conc_a, conc_b, effect
        - metadata: Run metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_bliss_independence(
        ...     cl_a=1.0, v_a=10.0, doses_a=[{"time": 0.0, "amount": 100.0}],
        ...     cl_b=2.0, v_b=20.0, doses_b=[{"time": 0.0, "amount": 50.0}],
        ...     e0=0.0, emax_a=80.0, ec50_a=5.0, gamma_a=1.0,
        ...     emax_b=60.0, ec50_b=3.0, gamma_b=1.0,
        ...     t0=0.0, t1=24.0, saveat=list(range(25))
        ... )
    """
    jl = _require_julia()

    import numpy as np

    # Simulate PK for drug A
    pk_spec_a = _create_pk_spec(jl, pk_kind_a, cl_a, v_a, ka_a, doses_a)
    grid = jl.NeoPKPD.SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = jl.NeoPKPD.SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))
    res_a = jl.NeoPKPD.simulate(pk_spec_a, grid, solver)
    conc_a = np.array([float(x) for x in res_a.observations[jl.Symbol("conc")]])

    # Simulate PK for drug B
    pk_spec_b = _create_pk_spec(jl, pk_kind_b, cl_b, v_b, ka_b, doses_b)
    res_b = jl.NeoPKPD.simulate(pk_spec_b, grid, solver)
    conc_b = np.array([float(x) for x in res_b.observations[jl.Symbol("conc")]])

    # Calculate Bliss Independence combined effect
    times = np.array([float(x) for x in res_a.t])
    effect = np.zeros(len(times))

    for i in range(len(times)):
        c_a = max(conc_a[i], 0.0)
        c_b = max(conc_b[i], 0.0)

        # Fractional effects
        if c_a <= 0.0:
            f_a = 0.0
        else:
            ratio_a = (ec50_a / c_a) ** gamma_a
            f_a = 1.0 / (1.0 + ratio_a)

        if c_b <= 0.0:
            f_b = 0.0
        else:
            ratio_b = (ec50_b / c_b) ** gamma_b
            f_b = 1.0 / (1.0 + ratio_b)

        # Bliss Independence
        f_combined = f_a + f_b - f_a * f_b
        emax_combined = max(emax_a, emax_b)
        effect[i] = e0 + f_combined * emax_combined

    return {
        "times": times.tolist(),
        "observations": {
            "conc_a": conc_a.tolist(),
            "conc_b": conc_b.tolist(),
            "effect": effect.tolist(),
        },
        "metadata": {
            "model": "BlissIndependence",
            "pk_kind_a": pk_kind_a,
            "pk_kind_b": pk_kind_b,
        }
    }


def simulate_pkpd_competitive_inhibition(
    # Drug PK parameters
    cl_drug: float,
    v_drug: float,
    doses_drug: List[Dict[str, Union[float, int]]],
    # Inhibitor PK parameters
    cl_inhib: float,
    v_inhib: float,
    doses_inhib: List[Dict[str, Union[float, int]]],
    # PD parameters
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
    ki: float,
    # Simulation settings
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind_drug: str = "OneCompIVBolus",
    pk_kind_inhib: str = "OneCompIVBolus",
    ka_drug: Optional[float] = None,
    ka_inhib: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a PK-PD simulation with Competitive Inhibition model.

    Models the effect of a competitive inhibitor on drug action.
    The inhibitor competes for the same binding site, shifting the EC50.

    PD equation:
        EC50_apparent = EC50 × (1 + I/Ki)
        E = E0 + Emax × C^γ / (EC50_apparent^γ + C^γ)

    Clinical Applications:
        - Drug-drug interactions
        - Receptor antagonism
        - Enzyme inhibition

    Args:
        cl_drug: Clearance for the drug (agonist)
        v_drug: Volume of distribution for the drug
        doses_drug: Dose events for the drug
        cl_inhib: Clearance for the inhibitor
        v_inhib: Volume of distribution for the inhibitor
        doses_inhib: Dose events for the inhibitor
        e0: Baseline effect
        emax: Maximum effect
        ec50: Intrinsic EC50 (without inhibitor)
        gamma: Hill coefficient
        ki: Inhibitor binding constant (lower = more potent inhibitor)
        t0: Simulation start time
        t1: Simulation end time
        saveat: Time points for output
        pk_kind_drug: PK model type for drug
        pk_kind_inhib: PK model type for inhibitor
        ka_drug: Absorption rate for drug (if oral)
        ka_inhib: Absorption rate for inhibitor (if oral)
        alg: ODE solver algorithm
        reltol: Relative tolerance
        abstol: Absolute tolerance
        maxiters: Maximum solver iterations

    Returns:
        Dict with keys:
        - times: Time points
        - observations: Dict with conc_drug, conc_inhibitor, effect
        - metadata: Run metadata

    Example:
        >>> result = neopkpd.simulate_pkpd_competitive_inhibition(
        ...     cl_drug=1.0, v_drug=10.0, doses_drug=[{"time": 0.0, "amount": 100.0}],
        ...     cl_inhib=0.5, v_inhib=15.0, doses_inhib=[{"time": 0.0, "amount": 50.0}],
        ...     e0=0.0, emax=100.0, ec50=5.0, gamma=1.0, ki=2.0,
        ...     t0=0.0, t1=24.0, saveat=list(range(25))
        ... )
    """
    jl = _require_julia()

    import numpy as np

    # Simulate PK for drug
    pk_spec_drug = _create_pk_spec(jl, pk_kind_drug, cl_drug, v_drug, ka_drug, doses_drug)
    grid = jl.NeoPKPD.SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = jl.NeoPKPD.SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))
    res_drug = jl.NeoPKPD.simulate(pk_spec_drug, grid, solver)
    conc_drug = np.array([float(x) for x in res_drug.observations[jl.Symbol("conc")]])

    # Simulate PK for inhibitor
    pk_spec_inhib = _create_pk_spec(jl, pk_kind_inhib, cl_inhib, v_inhib, ka_inhib, doses_inhib)
    res_inhib = jl.NeoPKPD.simulate(pk_spec_inhib, grid, solver)
    conc_inhib = np.array([float(x) for x in res_inhib.observations[jl.Symbol("conc")]])

    # Calculate competitive inhibition effect
    times = np.array([float(x) for x in res_drug.t])
    effect = np.zeros(len(times))

    for i in range(len(times)):
        c = max(conc_drug[i], 0.0)
        inhib = max(conc_inhib[i], 0.0)

        # Apparent EC50 with competitive inhibition
        ec50_apparent = ec50 * (1.0 + inhib / ki)

        if c <= 0.0:
            effect[i] = e0
        else:
            c_gamma = c ** gamma
            ec50_app_gamma = ec50_apparent ** gamma
            effect[i] = e0 + emax * c_gamma / (ec50_app_gamma + c_gamma)

    return {
        "times": times.tolist(),
        "observations": {
            "conc_drug": conc_drug.tolist(),
            "conc_inhibitor": conc_inhib.tolist(),
            "effect": effect.tolist(),
        },
        "metadata": {
            "model": "CompetitiveInhibition",
            "pk_kind_drug": pk_kind_drug,
            "pk_kind_inhib": pk_kind_inhib,
        }
    }


def simulate_pkpd_drug_interaction(
    # Drug A PK parameters
    cl_a: float,
    v_a: float,
    doses_a: List[Dict[str, Union[float, int]]],
    # Drug B PK parameters
    cl_b: float,
    v_b: float,
    doses_b: List[Dict[str, Union[float, int]]],
    # PD parameters
    e0: float,
    emax: float,
    ec50_a: float,
    ec50_b: float,
    psi: float,
    # Simulation settings
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind_a: str = "OneCompIVBolus",
    pk_kind_b: str = "OneCompIVBolus",
    ka_a: Optional[float] = None,
    ka_b: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a combination PK-PD simulation with Drug Interaction (Greco) model.

    The Greco model quantifies synergy or antagonism between two drugs
    using an interaction parameter (psi/ψ).

    PD equation (Greco model, γ=1):
        a = C_A / EC50_A
        b = C_B / EC50_B
        interaction = ψ × a × b
        E = E0 + Emax × (a + b + interaction) / (1 + a + b + interaction)

    Interaction parameter ψ:
        - ψ > 0: Synergy (greater than additive effect)
        - ψ = 0: Additivity (Loewe additivity)
        - ψ < 0: Antagonism (less than additive effect)

    Clinical Applications:
        - Synergy detection in combination therapy
        - Drug interaction studies
        - Isobologram analysis

    Args:
        cl_a: Clearance for drug A
        v_a: Volume of distribution for drug A
        doses_a: Dose events for drug A
        cl_b: Clearance for drug B
        v_b: Volume of distribution for drug B
        doses_b: Dose events for drug B
        e0: Baseline effect
        emax: Maximum combined effect
        ec50_a: EC50 for drug A
        ec50_b: EC50 for drug B
        psi: Interaction parameter (>0 synergy, <0 antagonism, =0 additive)
        t0: Simulation start time
        t1: Simulation end time
        saveat: Time points for output
        pk_kind_a: PK model type for drug A
        pk_kind_b: PK model type for drug B
        ka_a: Absorption rate for drug A (if oral)
        ka_b: Absorption rate for drug B (if oral)
        alg: ODE solver algorithm
        reltol: Relative tolerance
        abstol: Absolute tolerance
        maxiters: Maximum solver iterations

    Returns:
        Dict with keys:
        - times: Time points
        - observations: Dict with conc_a, conc_b, effect
        - metadata: Run metadata including synergy classification

    Example:
        >>> # Synergistic combination (psi > 0)
        >>> result = neopkpd.simulate_pkpd_drug_interaction(
        ...     cl_a=1.0, v_a=10.0, doses_a=[{"time": 0.0, "amount": 50.0}],
        ...     cl_b=2.0, v_b=20.0, doses_b=[{"time": 0.0, "amount": 30.0}],
        ...     e0=0.0, emax=100.0, ec50_a=5.0, ec50_b=3.0, psi=2.0,
        ...     t0=0.0, t1=24.0, saveat=list(range(25))
        ... )
        >>> print(f"Max effect: {max(result['observations']['effect']):.1f}")
    """
    jl = _require_julia()

    import numpy as np

    # Simulate PK for drug A
    pk_spec_a = _create_pk_spec(jl, pk_kind_a, cl_a, v_a, ka_a, doses_a)
    grid = jl.NeoPKPD.SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = jl.NeoPKPD.SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))
    res_a = jl.NeoPKPD.simulate(pk_spec_a, grid, solver)
    conc_a = np.array([float(x) for x in res_a.observations[jl.Symbol("conc")]])

    # Simulate PK for drug B
    pk_spec_b = _create_pk_spec(jl, pk_kind_b, cl_b, v_b, ka_b, doses_b)
    res_b = jl.NeoPKPD.simulate(pk_spec_b, grid, solver)
    conc_b = np.array([float(x) for x in res_b.observations[jl.Symbol("conc")]])

    # Calculate Greco drug interaction effect
    times = np.array([float(x) for x in res_a.t])
    effect = np.zeros(len(times))

    for i in range(len(times)):
        c_a = max(conc_a[i], 0.0)
        c_b = max(conc_b[i], 0.0)

        # Normalized concentrations
        a = c_a / ec50_a
        b = c_b / ec50_b

        # Interaction term (Greco model)
        interaction = psi * a * b

        # Combined effect
        numerator = a + b + interaction
        denominator = 1.0 + a + b + interaction

        if denominator <= 0.0:
            # Handle edge case (strong antagonism)
            effect[i] = e0
        else:
            effect[i] = e0 + emax * numerator / denominator

    # Classify interaction type
    if psi > 0.1:
        interaction_type = "synergy"
    elif psi < -0.1:
        interaction_type = "antagonism"
    else:
        interaction_type = "additive"

    return {
        "times": times.tolist(),
        "observations": {
            "conc_a": conc_a.tolist(),
            "conc_b": conc_b.tolist(),
            "effect": effect.tolist(),
        },
        "metadata": {
            "model": "DrugInteraction_Greco",
            "pk_kind_a": pk_kind_a,
            "pk_kind_b": pk_kind_b,
            "psi": psi,
            "interaction_type": interaction_type,
        }
    }
