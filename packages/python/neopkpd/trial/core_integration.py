"""
NeoPKPD Trial Core Integration

Julia-connected clinical trial simulation using actual PK/PD models.
Provides:
- Model-based subject exposure simulation
- Dose escalation algorithms (3+3, mTPI, CRM)
- Crossover analysis with period/sequence effects
- Adaptive trial simulation with interim analyses
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np

from .designs import (
    ParallelDesign,
    CrossoverDesign,
    DoseEscalationDesign,
    AdaptiveDesign,
    BioequivalenceDesign,
)
from .regimens import DosingRegimen, TitrationRegimen


# ============================================================================
# Result Types
# ============================================================================


@dataclass
class SubjectExposure:
    """
    Individual subject PK/PD exposure from model-based simulation.

    Attributes:
        subject_id: Subject identifier
        times: Observation times (hours)
        concentrations: Drug concentrations
        pd_response: PD response values (if applicable)
        pk_metrics: Derived PK metrics (Cmax, Tmax, AUC, etc.)
    """
    subject_id: str
    times: np.ndarray
    concentrations: np.ndarray
    pd_response: Optional[np.ndarray] = None
    pk_metrics: Dict[str, float] = field(default_factory=dict)


@dataclass
class CohortResult:
    """
    Dose escalation cohort result.

    Attributes:
        cohort_number: Cohort number
        dose_level: Dose level index
        dose_amount: Dose amount
        n_subjects: Number of subjects
        n_dlt: Number of DLTs
        subject_ids: Subject identifiers
        dlt_flags: DLT flags for each subject
        pk_exposures: PK exposures for each subject
        decision: Escalation decision
    """
    cohort_number: int
    dose_level: int
    dose_amount: float
    n_subjects: int
    n_dlt: int
    subject_ids: List[str]
    dlt_flags: List[bool]
    pk_exposures: List[float]
    decision: str


@dataclass
class EscalationResult:
    """
    Dose escalation simulation result.

    Attributes:
        design_type: Escalation design type
        cohorts: List of cohort results
        mtd_level: MTD dose level index
        mtd_dose: MTD dose amount
        total_subjects: Total subjects enrolled
        total_dlt: Total DLTs observed
        dlt_rate_by_dose: DLT rate at each dose level
        n_subjects_by_dose: Number of subjects at each dose level
        completed: Whether escalation completed
        termination_reason: Reason for termination
        seed: Random seed used
    """
    design_type: str
    cohorts: List[CohortResult]
    mtd_level: Optional[int]
    mtd_dose: Optional[float]
    total_subjects: int
    total_dlt: int
    dlt_rate_by_dose: Dict[int, float]
    n_subjects_by_dose: Dict[int, int]
    completed: bool
    termination_reason: str
    seed: int


@dataclass
class CrossoverAnalysis:
    """
    Crossover study analysis result.

    Attributes:
        treatment_effect: Treatment effect estimate
        treatment_se: Standard error of treatment effect
        treatment_ci: 90% CI for treatment effect
        period_effect: Period effect (p-value)
        sequence_effect: Sequence/carryover effect (p-value)
        within_subject_cv: Within-subject coefficient of variation
        be_assessment: Bioequivalence assessment result
    """
    treatment_effect: float
    treatment_se: float
    treatment_ci: Tuple[float, float]
    period_effect: Dict[str, float]
    sequence_effect: Dict[str, float]
    within_subject_cv: float
    be_assessment: Optional[Dict[str, Any]] = None


@dataclass
class InterimResult:
    """
    Interim analysis result.

    Attributes:
        analysis_number: Interim analysis number
        information_fraction: Information fraction
        n_enrolled: Subjects enrolled at interim
        effect_estimate: Effect size estimate
        effect_se: Standard error
        z_stat: Z-statistic
        p_value: P-value
        efficacy_boundary: Efficacy stopping boundary
        futility_boundary: Futility stopping boundary
        conditional_power: Conditional power
        stop_for_efficacy: Whether to stop for efficacy
        stop_for_futility: Whether to stop for futility
    """
    analysis_number: int
    information_fraction: float
    n_enrolled: int
    effect_estimate: float
    effect_se: float
    z_stat: float
    p_value: float
    efficacy_boundary: float
    futility_boundary: float
    conditional_power: float
    stop_for_efficacy: bool
    stop_for_futility: bool


@dataclass
class AdaptiveTrialResult:
    """
    Adaptive trial simulation result.

    Attributes:
        final_effect: Final effect estimate
        final_se: Final standard error
        final_p_value: Final p-value
        interim_results: List of interim analysis results
        stopped_early: Whether trial stopped early
        stop_reason: Reason for early stopping
        final_n: Final sample size
        initial_n: Initial planned sample size
        seed: Random seed
    """
    final_effect: float
    final_se: float
    final_p_value: float
    interim_results: List[InterimResult]
    stopped_early: bool
    stop_reason: Optional[str]
    final_n: int
    initial_n: int
    seed: int


# ============================================================================
# Julia Connection Utilities
# ============================================================================


def _get_julia():
    """Get the Julia connection from neopkpd."""
    try:
        from neopkpd.bridge import get_julia
        return get_julia()
    except ImportError:
        raise RuntimeError(
            "Julia bridge not available. Please ensure Julia is installed "
            "and the neopkpd package is properly configured."
        )


def _to_julia_float_vec(jl, py_list: List[float]):
    """Convert Python list to Julia Vector{Float64}."""
    vec = jl.seval("Float64[]")
    for x in py_list:
        jl.seval("push!")(vec, float(x))
    return vec


def _to_julia_vec(jl, py_list: list, elem_type) -> Any:
    """Convert Python list to Julia Vector of the specified element type."""
    vec = jl.Vector[elem_type](jl.undef, len(py_list))
    for i, item in enumerate(py_list):
        vec[i] = item  # PythonCall uses 0-based indexing from Python
    return vec


def _to_julia_int_vec(jl, py_list: List[int]):
    """Convert Python list to Julia Vector{Int}."""
    vec = jl.seval("Int[]")
    for x in py_list:
        jl.seval("push!")(vec, int(x))
    return vec


# ============================================================================
# Subject Exposure Simulation
# ============================================================================


def simulate_subject_exposure(
    model_kind: str,
    dose_events: List[Dict[str, float]],
    observation_times: List[float],
    subject: Optional[Dict[str, Any]] = None,
    pk_params: Optional[Dict[str, float]] = None,
    grid: Optional[Dict[str, Any]] = None,
    solver: Optional[Dict[str, Any]] = None,
    include_iiv: bool = True,
    omega: Optional[List[List[float]]] = None,
    seed: int = 12345,
) -> SubjectExposure:
    """
    Simulate PK/PD exposure for a single subject using Julia models.

    Args:
        model_kind: PK model type (e.g., 'OneCompIVBolus', 'TwoCompOral')
        dose_events: List of dose events with 'time', 'amount', 'rate' keys
        observation_times: Times to observe concentrations (hours)
        subject: Subject covariates (id, weight, age, etc.)
        pk_params: Base PK parameters (CL, V, etc.)
        grid: Simulation grid specification
        solver: ODE solver specification
        include_iiv: Whether to include inter-individual variability
        omega: Omega matrix for IIV (if include_iiv=True)
        seed: Random seed

    Returns:
        SubjectExposure with times, concentrations, and PK metrics

    Example:
        >>> exposure = simulate_subject_exposure(
        ...     model_kind='OneCompIVBolus',
        ...     dose_events=[{'time': 0.0, 'amount': 100.0, 'rate': 0.0}],
        ...     observation_times=[0, 0.5, 1, 2, 4, 8, 12, 24],
        ...     pk_params={'CL': 10.0, 'V': 50.0}
        ... )
        >>> print(f"Cmax: {exposure.pk_metrics['cmax']:.2f}")
    """
    jl = _get_julia()

    # Default parameters if not provided
    if pk_params is None:
        pk_params = {'CL': 10.0, 'V': 50.0}

    # Build Julia model spec - use correct Julia type names
    if model_kind == 'OneCompIVBolus':
        pk_model = jl.NeoPKPD.OneCompIVBolus()
        pk_params_jl = jl.NeoPKPD.OneCompIVBolusParams(
            float(pk_params.get('CL', 10.0)),
            float(pk_params.get('V', 50.0))
        )
    elif model_kind in ['OneCompOral', 'OneCompOralFirstOrder']:
        pk_model = jl.NeoPKPD.OneCompOralFirstOrder()
        pk_params_jl = jl.NeoPKPD.OneCompOralFirstOrderParams(
            float(pk_params.get('Ka', 1.0)),
            float(pk_params.get('CL', 10.0)),
            float(pk_params.get('V', 50.0))
        )
    elif model_kind == 'TwoCompIVBolus':
        pk_model = jl.NeoPKPD.TwoCompIVBolus()
        pk_params_jl = jl.NeoPKPD.TwoCompIVBolusParams(
            float(pk_params.get('CL', 10.0)),
            float(pk_params.get('V1', 50.0)),
            float(pk_params.get('Q', 5.0)),
            float(pk_params.get('V2', 100.0))
        )
    elif model_kind == 'TwoCompOral':
        pk_model = jl.NeoPKPD.TwoCompOral()
        pk_params_jl = jl.NeoPKPD.TwoCompOralParams(
            float(pk_params.get('Ka', 1.0)),
            float(pk_params.get('CL', 10.0)),
            float(pk_params.get('V1', 50.0)),
            float(pk_params.get('Q', 5.0)),
            float(pk_params.get('V2', 100.0))
        )
    else:
        raise ValueError(f"Unknown model kind: {model_kind}")

    # Build dosing events as Julia vector
    doses_list = []
    for dose in dose_events:
        dose_event = jl.NeoPKPD.DoseEvent(
            float(dose['time']),
            float(dose['amount']),
            float(dose.get('duration', dose.get('rate', 0.0)))  # duration for infusion
        )
        doses_list.append(dose_event)

    # Convert to Julia Vector{DoseEvent}
    DoseEvent = jl.NeoPKPD.DoseEvent
    doses = _to_julia_vec(jl, doses_list, DoseEvent)

    # Build ModelSpec with doses
    subject_id = subject.get('id', '1') if subject else '1'
    model_spec = jl.NeoPKPD.ModelSpec(
        pk_model,
        f"subject_{subject_id}",
        pk_params_jl,
        doses
    )

    # Build observation times
    obs_times = _to_julia_float_vec(jl, observation_times)

    # Build grid
    if grid is None:
        t_max = max(observation_times) if observation_times else 24.0
        grid_jl = jl.NeoPKPD.SimGrid(0.0, t_max, obs_times)
    else:
        saveat = _to_julia_float_vec(jl, [float(x) for x in grid.get('saveat', observation_times)])
        grid_jl = jl.NeoPKPD.SimGrid(
            float(grid.get('t0', 0.0)),
            float(grid.get('t1', max(observation_times))),
            saveat
        )

    # Build solver
    if solver is None:
        solver_jl = jl.NeoPKPD.SolverSpec(jl.Symbol("Tsit5"), 1e-8, 1e-10, 10**6)
    else:
        solver_jl = jl.NeoPKPD.SolverSpec(
            jl.Symbol(solver.get('alg', 'Tsit5')),
            float(solver.get('reltol', 1e-8)),
            float(solver.get('abstol', 1e-10)),
            int(solver.get('maxiters', 10**6))
        )

    # Run simulation
    result = jl.NeoPKPD.simulate(model_spec, grid_jl, solver_jl)

    # Extract results - SimResult has t, observations, and states
    times = np.array([float(t) for t in result.t])
    # Get concentration from observations (typically 'conc')
    conc_obs = result.observations[jl.Symbol("conc")]
    concentrations = np.array([float(c) for c in conc_obs])

    # Calculate PK metrics
    pk_metrics = {}
    if len(concentrations) > 0:
        pk_metrics['cmax'] = float(np.max(concentrations))
        pk_metrics['tmax'] = float(times[np.argmax(concentrations)])
        pk_metrics['auc_0_last'] = float(np.trapezoid(concentrations, times))

        # Calculate half-life from terminal phase if enough points
        if len(concentrations) >= 3:
            try:
                # Find terminal phase (last few points where conc is declining)
                for i in range(len(concentrations) - 1, 0, -1):
                    if concentrations[i] > 0 and concentrations[i-1] > 0:
                        lambda_z = (np.log(concentrations[i-1]) - np.log(concentrations[i])) / (times[i] - times[i-1])
                        if lambda_z > 0:
                            pk_metrics['lambda_z'] = float(lambda_z)
                            pk_metrics['t_half'] = float(np.log(2) / lambda_z)
                            break
            except Exception:
                pass

    subject_id = subject.get('id', '1') if subject else '1'

    return SubjectExposure(
        subject_id=str(subject_id),
        times=times,
        concentrations=concentrations,
        pd_response=None,
        pk_metrics=pk_metrics
    )


# ============================================================================
# Dose Escalation Simulation
# ============================================================================


def simulate_dose_escalation_3plus3(
    dose_levels: List[float],
    model_kind: str,
    pk_params: Dict[str, float],
    dlt_threshold: float = 100.0,
    grid: Optional[Dict[str, Any]] = None,
    solver: Optional[Dict[str, Any]] = None,
    seed: int = 12345,
    verbose: bool = False,
) -> EscalationResult:
    """
    Simulate 3+3 dose escalation study.

    The 3+3 design is the most common Phase I design:
    - Start with 3 subjects at starting dose
    - If 0/3 DLTs: escalate
    - If 1/3 DLTs: expand to 6
    - If 2+/3 or 2+/6 DLTs: MTD is previous dose
    - If 1/6 DLTs: escalate

    Args:
        dose_levels: List of dose levels to test
        model_kind: PK model type
        pk_params: Base PK parameters
        dlt_threshold: Concentration threshold for DLT
        grid: Simulation grid
        solver: ODE solver
        seed: Random seed
        verbose: Print progress

    Returns:
        EscalationResult with cohorts and MTD determination

    Example:
        >>> result = simulate_dose_escalation_3plus3(
        ...     dose_levels=[10, 25, 50, 100, 200],
        ...     model_kind='OneCompIVBolus',
        ...     pk_params={'CL': 10.0, 'V': 50.0},
        ...     dlt_threshold=100.0
        ... )
        >>> print(f"MTD: {result.mtd_dose}")
    """
    np.random.seed(seed)

    cohorts: List[CohortResult] = []
    current_level = 0
    total_subjects = 0
    total_dlt = 0
    n_by_dose: Dict[int, int] = {}
    dlt_by_dose: Dict[int, int] = {}
    cohort_num = 0

    observation_times = [0, 0.5, 1, 2, 4, 8, 12, 24]

    while current_level < len(dose_levels):
        cohort_num += 1
        dose = dose_levels[current_level]
        n_at_level = n_by_dose.get(current_level, 0)

        # Determine cohort size
        cohort_size = 3 if n_at_level == 0 else 3  # Always add 3

        if verbose:
            print(f"Cohort {cohort_num}: Dose level {current_level} ({dose}), n={cohort_size}")

        # Simulate cohort
        subject_ids = []
        dlt_flags = []
        pk_exposures = []

        for i in range(cohort_size):
            subject_id = f"S{total_subjects + 1:03d}"
            subject_ids.append(subject_id)
            total_subjects += 1

            # Simulate exposure
            try:
                exposure = simulate_subject_exposure(
                    model_kind=model_kind,
                    dose_events=[{'time': 0.0, 'amount': dose, 'rate': 0.0}],
                    observation_times=observation_times,
                    pk_params=pk_params,
                    grid=grid,
                    solver=solver,
                    seed=seed + total_subjects
                )
                cmax = exposure.pk_metrics.get('cmax', 0.0)
                pk_exposures.append(cmax)

                # DLT based on threshold
                has_dlt = cmax > dlt_threshold
                dlt_flags.append(has_dlt)
            except Exception as e:
                if verbose:
                    print(f"  Subject {subject_id} simulation failed: {e}")
                pk_exposures.append(0.0)
                dlt_flags.append(False)

        n_dlt = sum(dlt_flags)
        total_dlt += n_dlt

        # Update counts
        n_by_dose[current_level] = n_by_dose.get(current_level, 0) + cohort_size
        dlt_by_dose[current_level] = dlt_by_dose.get(current_level, 0) + n_dlt

        # Decision logic
        total_at_level = n_by_dose[current_level]
        total_dlt_at_level = dlt_by_dose[current_level]

        if total_at_level == 3:
            if total_dlt_at_level == 0:
                decision = "escalate"
            elif total_dlt_at_level == 1:
                decision = "expand"
            else:  # 2 or more
                decision = "de-escalate"
        else:  # total_at_level == 6
            if total_dlt_at_level <= 1:
                decision = "escalate"
            else:
                decision = "de-escalate"

        cohorts.append(CohortResult(
            cohort_number=cohort_num,
            dose_level=current_level,
            dose_amount=dose,
            n_subjects=cohort_size,
            n_dlt=n_dlt,
            subject_ids=subject_ids,
            dlt_flags=dlt_flags,
            pk_exposures=pk_exposures,
            decision=decision
        ))

        if verbose:
            print(f"  DLTs: {n_dlt}/{cohort_size}, Total: {total_dlt_at_level}/{total_at_level}, Decision: {decision}")

        # Apply decision
        if decision == "escalate":
            current_level += 1
        elif decision == "expand":
            pass  # Stay at same level, will add 3 more
        else:  # de-escalate
            break

    # Determine MTD
    mtd_level = None
    mtd_dose = None
    termination_reason = "completed"

    if cohorts:
        last_cohort = cohorts[-1]
        if last_cohort.decision == "de-escalate":
            if current_level > 0:
                mtd_level = current_level - 1
                mtd_dose = dose_levels[mtd_level]
                termination_reason = "mtd_found"
            else:
                mtd_level = None
                mtd_dose = None
                termination_reason = "too_toxic_at_lowest"
        elif current_level >= len(dose_levels):
            mtd_level = len(dose_levels) - 1
            mtd_dose = dose_levels[-1]
            termination_reason = "max_dose_reached"

    # Calculate DLT rates
    dlt_rates = {
        level: dlt_by_dose.get(level, 0) / n_by_dose.get(level, 1)
        for level in n_by_dose.keys()
    }

    return EscalationResult(
        design_type="3+3",
        cohorts=cohorts,
        mtd_level=mtd_level,
        mtd_dose=mtd_dose,
        total_subjects=total_subjects,
        total_dlt=total_dlt,
        dlt_rate_by_dose=dlt_rates,
        n_subjects_by_dose=n_by_dose,
        completed=True,
        termination_reason=termination_reason,
        seed=seed
    )


def simulate_dose_escalation_mtpi(
    dose_levels: List[float],
    model_kind: str,
    pk_params: Dict[str, float],
    dlt_threshold: float = 100.0,
    target_dlt_rate: float = 0.25,
    equivalence_interval: Tuple[float, float] = (0.20, 0.30),
    cohort_size: int = 3,
    max_subjects: int = 30,
    grid: Optional[Dict[str, Any]] = None,
    solver: Optional[Dict[str, Any]] = None,
    seed: int = 12345,
    verbose: bool = False,
) -> EscalationResult:
    """
    Simulate mTPI (modified Toxicity Probability Interval) dose escalation.

    mTPI uses a Bayesian framework with predefined probability intervals
    to make escalation decisions.

    Args:
        dose_levels: List of dose levels
        model_kind: PK model type
        pk_params: Base PK parameters
        dlt_threshold: Concentration threshold for DLT
        target_dlt_rate: Target DLT rate
        equivalence_interval: Equivalence interval around target
        cohort_size: Subjects per cohort
        max_subjects: Maximum total subjects
        grid: Simulation grid
        solver: ODE solver
        seed: Random seed
        verbose: Print progress

    Returns:
        EscalationResult with mTPI-based MTD

    Example:
        >>> result = simulate_dose_escalation_mtpi(
        ...     dose_levels=[10, 25, 50, 100, 200],
        ...     model_kind='OneCompIVBolus',
        ...     pk_params={'CL': 10.0, 'V': 50.0},
        ...     target_dlt_rate=0.25
        ... )
    """
    from scipy import stats

    np.random.seed(seed)

    cohorts: List[CohortResult] = []
    current_level = 0
    total_subjects = 0
    total_dlt = 0
    n_by_dose: Dict[int, int] = {i: 0 for i in range(len(dose_levels))}
    dlt_by_dose: Dict[int, int] = {i: 0 for i in range(len(dose_levels))}
    cohort_num = 0

    observation_times = [0, 0.5, 1, 2, 4, 8, 12, 24]
    ei_low, ei_high = equivalence_interval

    while total_subjects < max_subjects and current_level < len(dose_levels):
        cohort_num += 1
        dose = dose_levels[current_level]

        if verbose:
            print(f"Cohort {cohort_num}: Dose level {current_level} ({dose})")

        # Simulate cohort
        subject_ids = []
        dlt_flags = []
        pk_exposures = []

        for i in range(cohort_size):
            subject_id = f"S{total_subjects + 1:03d}"
            subject_ids.append(subject_id)
            total_subjects += 1

            try:
                exposure = simulate_subject_exposure(
                    model_kind=model_kind,
                    dose_events=[{'time': 0.0, 'amount': dose, 'rate': 0.0}],
                    observation_times=observation_times,
                    pk_params=pk_params,
                    grid=grid,
                    solver=solver,
                    seed=seed + total_subjects
                )
                cmax = exposure.pk_metrics.get('cmax', 0.0)
                pk_exposures.append(cmax)
                has_dlt = cmax > dlt_threshold
                dlt_flags.append(has_dlt)
            except Exception:
                pk_exposures.append(0.0)
                dlt_flags.append(False)

        n_dlt = sum(dlt_flags)
        total_dlt += n_dlt

        # Update counts
        n_by_dose[current_level] += cohort_size
        dlt_by_dose[current_level] += n_dlt

        # mTPI decision using Beta-Binomial
        n_current = n_by_dose[current_level]
        y_current = dlt_by_dose[current_level]

        # Calculate probability intervals using Beta(1,1) prior
        alpha_post = 1 + y_current
        beta_post = 1 + n_current - y_current

        prob_under = stats.beta.cdf(ei_low, alpha_post, beta_post)
        prob_eq = stats.beta.cdf(ei_high, alpha_post, beta_post) - prob_under
        prob_over = 1 - stats.beta.cdf(ei_high, alpha_post, beta_post)

        # Decision based on highest probability
        if prob_under >= prob_eq and prob_under >= prob_over:
            decision = "escalate"
        elif prob_over >= prob_eq and prob_over >= prob_under:
            decision = "de-escalate"
        else:
            decision = "stay"

        cohorts.append(CohortResult(
            cohort_number=cohort_num,
            dose_level=current_level,
            dose_amount=dose,
            n_subjects=cohort_size,
            n_dlt=n_dlt,
            subject_ids=subject_ids,
            dlt_flags=dlt_flags,
            pk_exposures=pk_exposures,
            decision=decision
        ))

        if verbose:
            print(f"  P(under)={prob_under:.3f}, P(eq)={prob_eq:.3f}, P(over)={prob_over:.3f}")
            print(f"  Decision: {decision}")

        # Apply decision
        if decision == "escalate" and current_level < len(dose_levels) - 1:
            current_level += 1
        elif decision == "de-escalate":
            if current_level > 0:
                current_level -= 1
            else:
                break  # Too toxic at lowest dose

    # Determine MTD as dose with posterior closest to target
    best_level = None
    best_diff = float('inf')

    for level, n in n_by_dose.items():
        if n > 0:
            y = dlt_by_dose[level]
            post_mean = (1 + y) / (2 + n)
            diff = abs(post_mean - target_dlt_rate)
            if diff < best_diff:
                best_diff = diff
                best_level = level

    mtd_level = best_level
    mtd_dose = dose_levels[best_level] if best_level is not None else None

    dlt_rates = {
        level: dlt_by_dose[level] / max(n_by_dose[level], 1)
        for level in n_by_dose.keys()
    }

    return EscalationResult(
        design_type="mTPI",
        cohorts=cohorts,
        mtd_level=mtd_level,
        mtd_dose=mtd_dose,
        total_subjects=total_subjects,
        total_dlt=total_dlt,
        dlt_rate_by_dose=dlt_rates,
        n_subjects_by_dose=n_by_dose,
        completed=total_subjects >= max_subjects or current_level >= len(dose_levels),
        termination_reason="max_subjects" if total_subjects >= max_subjects else "completed",
        seed=seed
    )


def simulate_dose_escalation_crm(
    dose_levels: List[float],
    model_kind: str,
    pk_params: Dict[str, float],
    dlt_threshold: float = 100.0,
    target_dlt_rate: float = 0.25,
    skeleton: Optional[List[float]] = None,
    cohort_size: int = 1,
    max_subjects: int = 30,
    grid: Optional[Dict[str, Any]] = None,
    solver: Optional[Dict[str, Any]] = None,
    seed: int = 12345,
    verbose: bool = False,
) -> EscalationResult:
    """
    Simulate CRM (Continual Reassessment Method) dose escalation.

    CRM uses a Bayesian model to estimate the dose-toxicity relationship
    and selects doses closest to the target DLT rate.

    Args:
        dose_levels: List of dose levels
        model_kind: PK model type
        pk_params: Base PK parameters
        dlt_threshold: Concentration threshold for DLT
        target_dlt_rate: Target DLT probability
        skeleton: Prior DLT probabilities for each dose (defaults to logistic model)
        cohort_size: Subjects per cohort (typically 1 for CRM)
        max_subjects: Maximum total subjects
        grid: Simulation grid
        solver: ODE solver
        seed: Random seed
        verbose: Print progress

    Returns:
        EscalationResult with CRM-based MTD

    Example:
        >>> result = simulate_dose_escalation_crm(
        ...     dose_levels=[10, 25, 50, 100, 200],
        ...     model_kind='OneCompIVBolus',
        ...     pk_params={'CL': 10.0, 'V': 50.0},
        ...     target_dlt_rate=0.25
        ... )
    """
    from scipy import stats
    from scipy.optimize import minimize_scalar

    np.random.seed(seed)

    n_doses = len(dose_levels)

    # Default skeleton using logistic model
    if skeleton is None:
        # Create skeleton assuming dose-proportional toxicity
        max_dose = max(dose_levels)
        skeleton = [0.05 + 0.20 * (d / max_dose) for d in dose_levels]

    cohorts: List[CohortResult] = []
    current_level = 0  # Start at lowest dose
    total_subjects = 0
    total_dlt = 0
    n_by_dose: Dict[int, int] = {i: 0 for i in range(n_doses)}
    dlt_by_dose: Dict[int, int] = {i: 0 for i in range(n_doses)}
    cohort_num = 0

    # CRM model parameter (power model: p(d) = skeleton[d]^exp(a))
    a_current = 0.0  # Start with prior

    observation_times = [0, 0.5, 1, 2, 4, 8, 12, 24]

    while total_subjects < max_subjects:
        cohort_num += 1
        dose = dose_levels[current_level]

        if verbose:
            print(f"Cohort {cohort_num}: Dose level {current_level} ({dose})")

        # Simulate cohort
        subject_ids = []
        dlt_flags = []
        pk_exposures = []

        for i in range(cohort_size):
            subject_id = f"S{total_subjects + 1:03d}"
            subject_ids.append(subject_id)
            total_subjects += 1

            try:
                exposure = simulate_subject_exposure(
                    model_kind=model_kind,
                    dose_events=[{'time': 0.0, 'amount': dose, 'rate': 0.0}],
                    observation_times=observation_times,
                    pk_params=pk_params,
                    grid=grid,
                    solver=solver,
                    seed=seed + total_subjects
                )
                cmax = exposure.pk_metrics.get('cmax', 0.0)
                pk_exposures.append(cmax)
                has_dlt = cmax > dlt_threshold
                dlt_flags.append(has_dlt)
            except Exception:
                pk_exposures.append(0.0)
                dlt_flags.append(False)

        n_dlt = sum(dlt_flags)
        total_dlt += n_dlt

        # Update counts
        n_by_dose[current_level] += cohort_size
        dlt_by_dose[current_level] += n_dlt

        # Update CRM model using all data
        # Likelihood: product of Bernoulli(p_i^exp(a)) for each observation
        def neg_log_likelihood(a):
            ll = 0.0
            for level in range(n_doses):
                n = n_by_dose[level]
                y = dlt_by_dose[level]
                if n > 0:
                    p = skeleton[level] ** np.exp(a)
                    p = max(1e-10, min(1 - 1e-10, p))  # Bound for numerical stability
                    ll += y * np.log(p) + (n - y) * np.log(1 - p)
            # Add normal prior on a
            ll -= 0.5 * a ** 2
            return -ll

        # Optimize
        result = minimize_scalar(neg_log_likelihood, bounds=(-3, 3), method='bounded')
        a_current = result.x

        # Calculate current DLT probability estimates
        prob_estimates = [skeleton[i] ** np.exp(a_current) for i in range(n_doses)]

        # Select next dose closest to target
        best_level = 0
        best_diff = float('inf')
        for i in range(n_doses):
            diff = abs(prob_estimates[i] - target_dlt_rate)
            if diff < best_diff:
                best_diff = diff
                best_level = i

        # Safety constraint: don't skip more than one dose level up
        if best_level > current_level + 1:
            best_level = current_level + 1

        decision = "escalate" if best_level > current_level else (
            "de-escalate" if best_level < current_level else "stay"
        )

        cohorts.append(CohortResult(
            cohort_number=cohort_num,
            dose_level=current_level,
            dose_amount=dose,
            n_subjects=cohort_size,
            n_dlt=n_dlt,
            subject_ids=subject_ids,
            dlt_flags=dlt_flags,
            pk_exposures=pk_exposures,
            decision=decision
        ))

        if verbose:
            print(f"  DLT: {n_dlt}/{cohort_size}")
            print(f"  Prob estimates: {[f'{p:.3f}' for p in prob_estimates]}")
            print(f"  Next dose level: {best_level}")

        current_level = best_level

    # Final MTD estimate
    prob_estimates = [skeleton[i] ** np.exp(a_current) for i in range(n_doses)]
    mtd_level = min(range(n_doses), key=lambda i: abs(prob_estimates[i] - target_dlt_rate))
    mtd_dose = dose_levels[mtd_level]

    dlt_rates = {
        level: dlt_by_dose[level] / max(n_by_dose[level], 1)
        for level in n_by_dose.keys()
    }

    return EscalationResult(
        design_type="CRM",
        cohorts=cohorts,
        mtd_level=mtd_level,
        mtd_dose=mtd_dose,
        total_subjects=total_subjects,
        total_dlt=total_dlt,
        dlt_rate_by_dose=dlt_rates,
        n_subjects_by_dose=n_by_dose,
        completed=True,
        termination_reason="max_subjects",
        seed=seed
    )


# ============================================================================
# Model-Connected Trial Simulation
# ============================================================================


@dataclass
class ModelTrialSubjectResult:
    """
    Subject result from model-connected trial simulation.

    Attributes:
        subject_id: Subject identifier
        arm_name: Treatment arm name
        completed: Whether subject completed trial
        dropout_day: Day of dropout (if applicable)
        times: Observation times
        concentrations: Drug concentrations
        pk_metrics: Derived PK metrics
        endpoint_values: Endpoint values
    """
    subject_id: str
    arm_name: str
    completed: bool
    dropout_day: Optional[float]
    times: np.ndarray
    concentrations: np.ndarray
    pk_metrics: Dict[str, float]
    endpoint_values: Dict[str, float]


@dataclass
class ModelTrialArmResult:
    """
    Treatment arm result from model-connected trial.

    Attributes:
        name: Arm name
        n_enrolled: Number enrolled
        n_completed: Number completed
        subjects: Individual subject results
        pk_summary: Summary of PK metrics
        endpoint_summary: Summary of endpoints
    """
    name: str
    n_enrolled: int
    n_completed: int
    subjects: List[ModelTrialSubjectResult]
    pk_summary: Dict[str, Dict[str, float]]
    endpoint_summary: Dict[str, Dict[str, float]]


@dataclass
class ModelTrialResult:
    """
    Complete model-connected trial result.

    Attributes:
        trial_name: Trial name
        arms: Results by arm
        comparisons: Statistical comparisons between arms
        seed: Random seed used
    """
    trial_name: str
    arms: Dict[str, ModelTrialArmResult]
    comparisons: Dict[str, Dict[str, float]]
    seed: int


def simulate_trial_with_model(
    trial_name: str,
    arms: List[Dict[str, Any]],
    model_kind: str,
    pk_params: Dict[str, float],
    observation_times: List[float],
    duration_days: float = 28.0,
    omega: Optional[List[List[float]]] = None,
    dropout_rate: float = 0.0,
    grid: Optional[Dict[str, Any]] = None,
    solver: Optional[Dict[str, Any]] = None,
    seed: int = 12345,
) -> ModelTrialResult:
    """
    Simulate a clinical trial using actual PK/PD models.

    This function connects to Julia PK models to simulate realistic
    drug exposures for each subject with inter-individual variability.

    Args:
        trial_name: Name of the trial
        arms: List of arm specifications, each with:
            - name: Arm name
            - dose: Dose amount
            - n_subjects: Number of subjects
            - placebo: Whether this is placebo arm (optional)
        model_kind: PK model type
        pk_params: Base PK parameters
        observation_times: PK sampling times (hours)
        duration_days: Trial duration
        omega: Omega matrix for IIV (optional)
        dropout_rate: Daily dropout rate
        grid: Simulation grid
        solver: ODE solver
        seed: Random seed

    Returns:
        ModelTrialResult with subject-level PK data and summaries

    Example:
        >>> result = simulate_trial_with_model(
        ...     trial_name="Phase 2 PK Study",
        ...     arms=[
        ...         {'name': 'Placebo', 'dose': 0.0, 'n_subjects': 10, 'placebo': True},
        ...         {'name': 'Low Dose', 'dose': 50.0, 'n_subjects': 10},
        ...         {'name': 'High Dose', 'dose': 100.0, 'n_subjects': 10},
        ...     ],
        ...     model_kind='OneCompIVBolus',
        ...     pk_params={'CL': 10.0, 'V': 50.0},
        ...     observation_times=[0, 1, 2, 4, 8, 12, 24]
        ... )
    """
    np.random.seed(seed)

    arm_results: Dict[str, ModelTrialArmResult] = {}
    subject_counter = 0

    for arm_spec in arms:
        arm_name = arm_spec['name']
        dose = arm_spec['dose']
        n_subjects = arm_spec['n_subjects']
        is_placebo = arm_spec.get('placebo', False)

        subjects: List[ModelTrialSubjectResult] = []

        for i in range(n_subjects):
            subject_counter += 1
            subject_id = f"S{subject_counter:03d}"

            # Simulate dropout
            dropout_day = None
            completed = True
            if dropout_rate > 0:
                for day in range(int(duration_days)):
                    if np.random.random() < dropout_rate:
                        dropout_day = float(day)
                        completed = False
                        break

            # Simulate PK
            if is_placebo or dose == 0:
                # Placebo - just noise
                times = np.array(observation_times)
                concentrations = np.abs(np.random.normal(0, 0.01, len(observation_times)))
                pk_metrics = {'cmax': 0.0, 'tmax': 0.0, 'auc_0_last': 0.0}
            else:
                try:
                    # Add IIV if omega provided
                    subject_params = pk_params.copy()
                    if omega is not None:
                        omega_mat = np.array(omega)
                        n_params = omega_mat.shape[0]
                        eta = np.random.multivariate_normal(np.zeros(n_params), omega_mat)

                        # Apply exponential IIV to CL and V
                        param_names = list(pk_params.keys())
                        for j, pname in enumerate(param_names[:n_params]):
                            if pname in subject_params:
                                subject_params[pname] = pk_params[pname] * np.exp(eta[j])

                    exposure = simulate_subject_exposure(
                        model_kind=model_kind,
                        dose_events=[{'time': 0.0, 'amount': dose, 'rate': 0.0}],
                        observation_times=observation_times,
                        pk_params=subject_params,
                        grid=grid,
                        solver=solver,
                        seed=seed + subject_counter
                    )
                    times = exposure.times
                    concentrations = exposure.concentrations
                    pk_metrics = exposure.pk_metrics
                except Exception as e:
                    times = np.array(observation_times)
                    concentrations = np.zeros(len(observation_times))
                    pk_metrics = {'cmax': 0.0, 'error': str(e)}

            # Calculate endpoints
            endpoint_values = {
                'cmax': pk_metrics.get('cmax', 0.0),
                'auc': pk_metrics.get('auc_0_last', 0.0),
            }

            subjects.append(ModelTrialSubjectResult(
                subject_id=subject_id,
                arm_name=arm_name,
                completed=completed,
                dropout_day=dropout_day,
                times=times,
                concentrations=concentrations,
                pk_metrics=pk_metrics,
                endpoint_values=endpoint_values
            ))

        # Summarize arm
        completed_subjects = [s for s in subjects if s.completed]
        n_completed = len(completed_subjects)

        # PK summary
        pk_summary = {}
        for metric in ['cmax', 'tmax', 'auc_0_last', 't_half']:
            values = [s.pk_metrics.get(metric, np.nan) for s in completed_subjects
                     if metric in s.pk_metrics and not np.isnan(s.pk_metrics.get(metric, np.nan))]
            if values:
                pk_summary[metric] = {
                    'n': len(values),
                    'mean': float(np.mean(values)),
                    'sd': float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
                    'median': float(np.median(values)),
                    'min': float(np.min(values)),
                    'max': float(np.max(values)),
                }

        # Endpoint summary
        endpoint_summary = {}
        for endpoint in ['cmax', 'auc']:
            values = [s.endpoint_values.get(endpoint, 0.0) for s in completed_subjects]
            if values:
                endpoint_summary[endpoint] = {
                    'n': len(values),
                    'mean': float(np.mean(values)),
                    'sd': float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
                }

        arm_results[arm_name] = ModelTrialArmResult(
            name=arm_name,
            n_enrolled=n_subjects,
            n_completed=n_completed,
            subjects=subjects,
            pk_summary=pk_summary,
            endpoint_summary=endpoint_summary
        )

    # Statistical comparisons between arms
    comparisons = {}
    arm_names = list(arm_results.keys())

    if len(arm_names) >= 2:
        from scipy import stats

        # Compare each arm to the first (reference)
        ref_arm = arm_names[0]
        ref_values = [s.endpoint_values.get('cmax', 0.0)
                     for s in arm_results[ref_arm].subjects if s.completed]

        for arm_name in arm_names[1:]:
            test_values = [s.endpoint_values.get('cmax', 0.0)
                          for s in arm_results[arm_name].subjects if s.completed]

            if len(ref_values) > 1 and len(test_values) > 1:
                t_stat, p_value = stats.ttest_ind(test_values, ref_values)
                diff = np.mean(test_values) - np.mean(ref_values)

                comparisons[f"{arm_name}_vs_{ref_arm}"] = {
                    'difference': float(diff),
                    't_statistic': float(t_stat),
                    'p_value': float(p_value),
                    'significant': p_value < 0.05
                }

    return ModelTrialResult(
        trial_name=trial_name,
        arms=arm_results,
        comparisons=comparisons,
        seed=seed
    )


# ============================================================================
# Crossover Analysis
# ============================================================================


def analyze_crossover(
    period1_values: List[float],
    period2_values: List[float],
    sequences: List[str],
    design: Optional[CrossoverDesign] = None,
    log_transform: bool = True,
) -> CrossoverAnalysis:
    """
    Analyze 2x2 crossover study data.

    Performs analysis of variance for crossover design including:
    - Treatment effect estimation
    - Period effect testing
    - Sequence (carryover) effect testing
    - Within-subject CV calculation
    - Bioequivalence assessment

    Args:
        period1_values: Values from period 1 for all subjects
        period2_values: Values from period 2 for all subjects
        sequences: Sequence assignment ('AB' or 'BA') for each subject
        design: Crossover design specification
        log_transform: Whether to log-transform for BE analysis

    Returns:
        CrossoverAnalysis with effects and BE assessment

    Example:
        >>> p1 = [100, 110, 95, 105, 102, 98]
        >>> p2 = [98, 108, 92, 100, 105, 95]
        >>> seq = ['AB', 'AB', 'AB', 'BA', 'BA', 'BA']
        >>> result = analyze_crossover(p1, p2, seq)
        >>> print(f"Treatment effect: {result.treatment_effect:.3f}")
    """
    import math

    p1 = np.array(period1_values)
    p2 = np.array(period2_values)
    seqs = np.array(sequences)

    n = len(p1)

    # Log transform if requested
    if log_transform:
        p1_log = np.log(p1)
        p2_log = np.log(p2)
    else:
        p1_log = p1
        p2_log = p2

    # Separate by sequence
    ab_mask = seqs == 'AB'
    ba_mask = seqs == 'BA'

    # For AB sequence: Period 1 = T, Period 2 = R
    # For BA sequence: Period 1 = R, Period 2 = T

    # Calculate subject differences (within-subject)
    diff = p1_log - p2_log  # For AB: T - R, For BA: R - T

    # Treatment effect = mean of (T - R) differences
    # For AB: diff = T - R (correct sign)
    # For BA: diff = R - T (need to negate)
    treatment_diffs = np.where(ab_mask, diff, -diff)
    treatment_effect = np.mean(treatment_diffs)
    treatment_se = np.std(treatment_diffs, ddof=1) / np.sqrt(n)

    # 90% CI for BE
    from scipy import stats
    t_crit = stats.t.ppf(0.95, df=n-2)
    ci_lower = treatment_effect - t_crit * treatment_se
    ci_upper = treatment_effect + t_crit * treatment_se

    # Period effect
    period_effect_value = np.mean(p1_log - p2_log)  # Overall period difference
    period_se = np.std(p1_log - p2_log, ddof=1) / np.sqrt(n)
    period_t = period_effect_value / period_se if period_se > 0 else 0
    period_p = 2 * (1 - stats.t.cdf(abs(period_t), df=n-1))

    # Sequence effect (carryover)
    ab_sums = p1_log[ab_mask] + p2_log[ab_mask]
    ba_sums = p1_log[ba_mask] + p2_log[ba_mask]
    if len(ab_sums) > 0 and len(ba_sums) > 0:
        seq_t, seq_p = stats.ttest_ind(ab_sums, ba_sums)
    else:
        seq_t, seq_p = 0, 1.0

    # Within-subject CV
    # CV = sqrt(exp(MSE) - 1) where MSE is residual variance
    residual_var = np.var(treatment_diffs, ddof=1) / 2  # Divide by 2 for within-subject
    if log_transform and residual_var > 0:
        within_cv = np.sqrt(np.exp(residual_var) - 1) * 100
    else:
        within_cv = np.std(treatment_diffs, ddof=1) / np.mean(np.exp(p1_log)) * 100

    # BE assessment
    be_assessment = None
    if log_transform:
        gmr = np.exp(treatment_effect)
        ci_lower_ratio = np.exp(ci_lower)
        ci_upper_ratio = np.exp(ci_upper)

        be_assessment = {
            'geometric_mean_ratio': gmr,
            'ci_90_lower': ci_lower_ratio,
            'ci_90_upper': ci_upper_ratio,
            'is_bioequivalent': 0.80 <= ci_lower_ratio and ci_upper_ratio <= 1.25,
            'theta_limits': (0.80, 1.25)
        }

    return CrossoverAnalysis(
        treatment_effect=float(treatment_effect),
        treatment_se=float(treatment_se),
        treatment_ci=(float(ci_lower), float(ci_upper)),
        period_effect={'estimate': float(period_effect_value), 'p_value': float(period_p)},
        sequence_effect={'t_statistic': float(seq_t), 'p_value': float(seq_p)},
        within_subject_cv=float(within_cv),
        be_assessment=be_assessment
    )


def test_period_effect(
    period1_values: List[float],
    period2_values: List[float],
) -> Dict[str, float]:
    """
    Test for period effect in crossover study.

    Args:
        period1_values: Values from period 1
        period2_values: Values from period 2

    Returns:
        Dict with effect estimate, t-statistic, and p-value
    """
    from scipy import stats

    p1 = np.array(period1_values)
    p2 = np.array(period2_values)

    diff = p1 - p2
    effect = np.mean(diff)
    se = np.std(diff, ddof=1) / np.sqrt(len(diff))
    t_stat = effect / se if se > 0 else 0
    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=len(diff)-1))

    return {
        'effect': float(effect),
        'se': float(se),
        't_statistic': float(t_stat),
        'p_value': float(p_value)
    }


def test_sequence_effect(
    period1_values: List[float],
    period2_values: List[float],
    sequences: List[str],
) -> Dict[str, float]:
    """
    Test for sequence (carryover) effect in crossover study.

    Args:
        period1_values: Values from period 1
        period2_values: Values from period 2
        sequences: Sequence assignments

    Returns:
        Dict with effect estimate, t-statistic, and p-value
    """
    from scipy import stats

    p1 = np.array(period1_values)
    p2 = np.array(period2_values)
    seqs = np.array(sequences)

    # Sum of both periods for each subject
    sums = p1 + p2

    ab_sums = sums[seqs == 'AB']
    ba_sums = sums[seqs == 'BA']

    if len(ab_sums) == 0 or len(ba_sums) == 0:
        return {'effect': 0.0, 't_statistic': 0.0, 'p_value': 1.0}

    effect = np.mean(ab_sums) - np.mean(ba_sums)
    t_stat, p_value = stats.ttest_ind(ab_sums, ba_sums)

    return {
        'effect': float(effect),
        't_statistic': float(t_stat),
        'p_value': float(p_value)
    }


def compute_within_subject_cv(
    values1: List[float],
    values2: List[float],
    log_scale: bool = True,
) -> float:
    """
    Compute within-subject coefficient of variation.

    Args:
        values1: First measurements
        values2: Second measurements (same subjects)
        log_scale: Whether data is on log scale

    Returns:
        Within-subject CV as percentage
    """
    v1 = np.array(values1)
    v2 = np.array(values2)

    if log_scale:
        diff = np.log(v1) - np.log(v2)
        mse = np.var(diff, ddof=1) / 2
        cv = np.sqrt(np.exp(mse) - 1) * 100
    else:
        diff = v1 - v2
        mean_val = (np.mean(v1) + np.mean(v2)) / 2
        cv = np.std(diff, ddof=1) / np.sqrt(2) / mean_val * 100

    return float(cv)


# ============================================================================
# Adaptive Trial Simulation
# ============================================================================


def simulate_adaptive_trial(
    n_per_arm: int,
    effect_size: float,
    sd: float,
    interim_fractions: List[float] = [0.5],
    alpha: float = 0.05,
    spending_type: str = 'obrien_fleming',
    futility_threshold: float = 0.10,
    seed: int = 12345,
    n_replicates: int = 1,
) -> AdaptiveTrialResult:
    """
    Simulate adaptive trial with interim analyses.

    Supports O'Brien-Fleming, Pocock, and other alpha spending functions
    for group sequential designs.

    Args:
        n_per_arm: Planned subjects per arm
        effect_size: True effect size
        sd: Standard deviation
        interim_fractions: Information fractions for interim analyses
        alpha: Overall alpha level
        spending_type: Alpha spending function type
        futility_threshold: Conditional power threshold for futility
        seed: Random seed
        n_replicates: Number of replicates (for power estimation)

    Returns:
        AdaptiveTrialResult with interim analyses and final outcome

    Example:
        >>> result = simulate_adaptive_trial(
        ...     n_per_arm=50,
        ...     effect_size=0.5,
        ...     sd=1.0,
        ...     interim_fractions=[0.5],
        ...     alpha=0.05
        ... )
        >>> print(f"Stopped early: {result.stopped_early}")
    """
    from scipy import stats
    from .analysis import alpha_spending_function

    np.random.seed(seed)

    # All analysis points (interim + final)
    all_fractions = sorted(interim_fractions + [1.0])

    # Calculate alpha at each analysis
    alphas = []
    cumulative_alpha = 0.0
    for frac in all_fractions:
        spent = alpha_spending_function(frac, alpha, spending_type)
        incremental = spent - cumulative_alpha
        alphas.append(incremental)
        cumulative_alpha = spent

    interim_results: List[InterimResult] = []
    stopped_early = False
    stop_reason = None
    final_n = 0

    # Generate all data upfront
    control = np.random.normal(0, sd, n_per_arm)
    treatment = np.random.normal(effect_size, sd, n_per_arm)

    for i, frac in enumerate(all_fractions):
        n_at_analysis = int(n_per_arm * frac)
        final_n = n_at_analysis * 2  # Both arms

        # Get data up to this point
        ctrl_data = control[:n_at_analysis]
        trt_data = treatment[:n_at_analysis]

        # Calculate statistics
        effect_est = np.mean(trt_data) - np.mean(ctrl_data)
        pooled_se = np.sqrt(np.var(ctrl_data, ddof=1)/n_at_analysis +
                           np.var(trt_data, ddof=1)/n_at_analysis)
        z_stat = effect_est / pooled_se if pooled_se > 0 else 0
        p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))

        # Efficacy boundary (two-sided)
        cumulative_alpha_spent = sum(alphas[:i+1])
        efficacy_z = stats.norm.ppf(1 - cumulative_alpha_spent/2)

        # Conditional power for futility
        if frac < 1.0:
            remaining_info = 1.0 - frac
            # Assume effect size is true for conditional power
            projected_z = z_stat * np.sqrt(1/frac) + effect_size/sd * np.sqrt(n_per_arm * remaining_info)
            cond_power = 1 - stats.norm.cdf(stats.norm.ppf(1-alpha/2) - projected_z)
        else:
            cond_power = float(p_value < alpha)

        futility_z = 0.0 if frac < 1.0 else efficacy_z

        stop_efficacy = abs(z_stat) > efficacy_z
        stop_futility = frac < 1.0 and cond_power < futility_threshold

        interim_results.append(InterimResult(
            analysis_number=i + 1,
            information_fraction=frac,
            n_enrolled=final_n,
            effect_estimate=float(effect_est),
            effect_se=float(pooled_se),
            z_stat=float(z_stat),
            p_value=float(p_value),
            efficacy_boundary=float(efficacy_z),
            futility_boundary=float(futility_z),
            conditional_power=float(cond_power),
            stop_for_efficacy=stop_efficacy,
            stop_for_futility=stop_futility
        ))

        if stop_efficacy:
            stopped_early = True
            stop_reason = "efficacy"
            break
        elif stop_futility:
            stopped_early = True
            stop_reason = "futility"
            break

    # Final analysis stats
    final_result = interim_results[-1]

    return AdaptiveTrialResult(
        final_effect=final_result.effect_estimate,
        final_se=final_result.effect_se,
        final_p_value=final_result.p_value,
        interim_results=interim_results,
        stopped_early=stopped_early,
        stop_reason=stop_reason,
        final_n=final_n,
        initial_n=n_per_arm * 2,
        seed=seed
    )


# ============================================================================
# PK Metrics Calculation
# ============================================================================


def calculate_pk_metrics(
    times: List[float],
    concentrations: List[float],
) -> Dict[str, float]:
    """
    Calculate PK metrics from concentration-time data.

    Args:
        times: Time points (hours)
        concentrations: Concentration values

    Returns:
        Dict with Cmax, Tmax, AUC_0_last, Lambda_z, T_half

    Example:
        >>> times = [0, 0.5, 1, 2, 4, 8, 12, 24]
        >>> concs = [0, 1.5, 2.0, 1.8, 1.2, 0.6, 0.3, 0.1]
        >>> metrics = calculate_pk_metrics(times, concs)
        >>> print(f"Cmax: {metrics['cmax']:.2f}, Tmax: {metrics['tmax']:.1f}h")
    """
    t = np.array(times)
    c = np.array(concentrations)

    metrics = {}

    if len(c) == 0:
        return metrics

    # Cmax and Tmax
    max_idx = np.argmax(c)
    metrics['cmax'] = float(c[max_idx])
    metrics['tmax'] = float(t[max_idx])

    # AUC by trapezoidal rule (handle both old and new numpy API)
    try:
        metrics['auc_0_last'] = float(np.trapezoid(c, t))
    except AttributeError:
        metrics['auc_0_last'] = float(np.trapezoid(c, t))

    # Terminal elimination rate constant and half-life
    # Use last 3 points in terminal phase
    terminal_points = 3
    if len(c) >= terminal_points:
        # Find terminal phase (declining portion after Cmax)
        terminal_start = max_idx + 1
        if terminal_start < len(c) - 1:
            t_term = t[terminal_start:]
            c_term = c[terminal_start:]

            # Only use positive concentrations
            mask = c_term > 0
            if np.sum(mask) >= 2:
                t_term = t_term[mask]
                c_term = c_term[mask]

                # Linear regression on log-concentration
                log_c = np.log(c_term)
                slope, intercept = np.polyfit(t_term, log_c, 1)

                if slope < 0:
                    metrics['lambda_z'] = float(-slope)
                    metrics['t_half'] = float(np.log(2) / (-slope))

    # AUC extrapolated to infinity
    if 'lambda_z' in metrics and c[-1] > 0:
        auc_extrap = c[-1] / metrics['lambda_z']
        metrics['auc_0_inf'] = metrics['auc_0_last'] + auc_extrap

    return metrics
