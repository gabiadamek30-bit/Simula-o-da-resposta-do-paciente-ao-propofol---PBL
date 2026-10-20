"""
NeoPKPD Trial Simulation

Trial simulation including:
- Trial specification and execution
- Dropout simulation
- Compliance modeling
- Treatment arm management
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union
import random
import math

from .designs import (
    ParallelDesign,
    CrossoverDesign,
    DoseEscalationDesign,
    AdaptiveDesign,
    BioequivalenceDesign,
)
from .regimens import DosingRegimen, TitrationRegimen
from .population import VirtualSubject, VirtualPopulationSpec, generate_virtual_population


@dataclass
class DropoutSpec:
    """
    Dropout specification for trial simulation.

    Attributes:
        random_rate_per_day: Random dropout rate per day
        ae_threshold: Adverse event threshold for dropout
        ae_dropout_prob: Probability of dropout given AE threshold exceeded

    Example:
        >>> spec = DropoutSpec(random_rate_per_day=0.005)
    """
    random_rate_per_day: float = 0.005
    ae_threshold: Optional[float] = None
    ae_dropout_prob: float = 0.5


@dataclass
class ComplianceSpec:
    """
    Compliance specification for trial simulation.

    Attributes:
        mean_compliance: Mean compliance rate (0-1)
        compliance_sd: Compliance standard deviation
        pattern: Compliance pattern ('random', 'weekend_miss', 'decay', 'early_good')

    Example:
        >>> spec = ComplianceSpec(mean_compliance=0.85, pattern='decay')
    """
    mean_compliance: float = 0.90
    compliance_sd: float = 0.10
    pattern: str = "random"


@dataclass
class TreatmentArm:
    """
    Treatment arm specification.

    Attributes:
        name: Arm name
        regimen: Dosing regimen
        n_subjects: Number of subjects in arm
        model_spec: Optional PK/PD model specification
        placebo: Whether this is a placebo arm

    Example:
        >>> from neopkpd.trial import dosing_qd
        >>> arm = TreatmentArm(
        ...     name="Active",
        ...     regimen=dosing_qd(100.0, 28),
        ...     n_subjects=50
        ... )
    """
    name: str
    regimen: Union[DosingRegimen, TitrationRegimen]
    n_subjects: int
    model_spec: Optional[Dict[str, Any]] = None
    placebo: bool = False


@dataclass
class TrialSpec:
    """
    Complete trial specification.

    Attributes:
        name: Trial name
        design: Study design
        arms: Treatment arms
        population_spec: Virtual population specification
        duration_days: Trial duration in days
        enrollment_rate: Subjects enrolled per day
        dropout: Dropout specification
        compliance: Compliance specification
        pk_sampling_times: PK sampling times (hours)
        endpoints: Endpoint names
        n_replicates: Number of simulation replicates
        seed: Random seed

    Example:
        >>> from neopkpd.trial import parallel_design, dosing_qd, TreatmentArm
        >>> spec = TrialSpec(
        ...     name="Phase 2 Study",
        ...     design=parallel_design(2),
        ...     arms=[
        ...         TreatmentArm("Placebo", dosing_qd(0.0, 28), 25, placebo=True),
        ...         TreatmentArm("Active", dosing_qd(100.0, 28), 25),
        ...     ],
        ...     duration_days=28
        ... )
    """
    name: str
    design: Union[ParallelDesign, CrossoverDesign, DoseEscalationDesign, AdaptiveDesign, BioequivalenceDesign]
    arms: List[TreatmentArm]
    population_spec: Optional[VirtualPopulationSpec] = None
    duration_days: float = 28.0
    enrollment_rate: float = 5.0
    dropout: Optional[DropoutSpec] = None
    compliance: Optional[ComplianceSpec] = None
    pk_sampling_times: List[float] = field(default_factory=lambda: [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0])
    endpoints: List[str] = field(default_factory=lambda: ["pk_exposure"])
    n_replicates: int = 1
    seed: int = 12345


@dataclass
class SubjectResult:
    """
    Individual subject result from trial simulation.

    Attributes:
        subject_id: Subject ID
        arm_name: Treatment arm name
        completed: Whether subject completed the trial
        dropout_day: Day of dropout (if applicable)
        compliance_rate: Actual compliance rate
        endpoint_values: Endpoint values
    """
    subject_id: int
    arm_name: str
    completed: bool
    dropout_day: Optional[float] = None
    compliance_rate: float = 1.0
    endpoint_values: Dict[str, float] = field(default_factory=dict)


@dataclass
class ArmResult:
    """
    Treatment arm result from trial simulation.

    Attributes:
        name: Arm name
        n_enrolled: Number enrolled
        n_completed: Number completed
        completion_rate: Completion rate
        mean_compliance: Mean compliance
        subjects: Individual subject results
        endpoint_summaries: Endpoint summary statistics
    """
    name: str
    n_enrolled: int
    n_completed: int
    completion_rate: float
    mean_compliance: float
    subjects: List[SubjectResult]
    endpoint_summaries: Dict[str, Dict[str, float]] = field(default_factory=dict)


@dataclass
class TrialResult:
    """
    Complete trial simulation result.

    Attributes:
        trial_name: Trial name
        arms: Results by arm
        overall_completion_rate: Overall completion rate
        overall_compliance: Overall mean compliance
        endpoint_comparisons: Endpoint comparison statistics
        replicate: Replicate number
    """
    trial_name: str
    arms: Dict[str, ArmResult]
    overall_completion_rate: float
    overall_compliance: float
    endpoint_comparisons: Dict[str, Dict[str, float]] = field(default_factory=dict)
    replicate: int = 1


def simulate_dropout(
    n_subjects: int,
    duration_days: float,
    spec: Optional[DropoutSpec] = None,
    seed: Optional[int] = None,
) -> List[Optional[float]]:
    """
    Simulate dropout times for a group of subjects.

    Args:
        n_subjects: Number of subjects
        duration_days: Trial duration in days
        spec: Dropout specification
        seed: Random seed

    Returns:
        List of dropout days (None if completed)

    Example:
        >>> dropout_days = simulate_dropout(100, 28.0, DropoutSpec(random_rate_per_day=0.01))
        >>> n_dropouts = sum(1 for d in dropout_days if d is not None)
    """
    if seed is not None:
        random.seed(seed)

    if spec is None:
        return [None] * n_subjects

    dropout_times: List[Optional[float]] = []

    for _ in range(n_subjects):
        # Simulate daily dropout
        dropped = False
        for day in range(int(duration_days)):
            if random.random() < spec.random_rate_per_day:
                dropout_times.append(float(day))
                dropped = True
                break

        if not dropped:
            dropout_times.append(None)

    return dropout_times


def apply_compliance(
    dose_times: List[float],
    dose_amounts: List[float],
    spec: Optional[ComplianceSpec] = None,
    seed: Optional[int] = None,
) -> List[float]:
    """
    Apply compliance model to dose schedule.

    Args:
        dose_times: Scheduled dose times
        dose_amounts: Scheduled dose amounts
        spec: Compliance specification
        seed: Random seed

    Returns:
        Actual dose amounts (0 for missed doses)

    Example:
        >>> times = [0, 24, 48, 72]
        >>> amounts = [100.0, 100.0, 100.0, 100.0]
        >>> actual = apply_compliance(times, amounts, ComplianceSpec(mean_compliance=0.8))
    """
    if seed is not None:
        random.seed(seed)

    if spec is None:
        return list(dose_amounts)

    actual_doses = []
    n_doses = len(dose_amounts)

    for i, dose in enumerate(dose_amounts):
        if spec.pattern == "random":
            # Random compliance
            if random.random() <= spec.mean_compliance:
                actual_doses.append(dose)
            else:
                actual_doses.append(0.0)

        elif spec.pattern == "weekend_miss":
            # Higher miss rate on weekends
            day = int(dose_times[i] / 24.0) % 7
            miss_prob = (1 - spec.mean_compliance) * 2 if day >= 5 else (1 - spec.mean_compliance)
            if random.random() > miss_prob:
                actual_doses.append(dose)
            else:
                actual_doses.append(0.0)

        elif spec.pattern == "decay":
            # Compliance decays over time
            relative_time = i / max(n_doses - 1, 1)
            current_compliance = spec.mean_compliance * (1 - 0.3 * relative_time)
            if random.random() <= current_compliance:
                actual_doses.append(dose)
            else:
                actual_doses.append(0.0)

        elif spec.pattern == "early_good":
            # Better compliance early
            relative_time = i / max(n_doses - 1, 1)
            if relative_time < 0.5:
                current_compliance = min(spec.mean_compliance + 0.1, 1.0)
            else:
                current_compliance = max(spec.mean_compliance - 0.1, 0.0)
            if random.random() <= current_compliance:
                actual_doses.append(dose)
            else:
                actual_doses.append(0.0)

        else:
            actual_doses.append(dose)

    return actual_doses


def _simulate_subject(
    subject: VirtualSubject,
    arm: TreatmentArm,
    trial_spec: TrialSpec,
    rng: random.Random,
) -> SubjectResult:
    """Simulate a single subject."""
    # Simulate dropout
    dropout_day = None
    if trial_spec.dropout is not None:
        for day in range(int(trial_spec.duration_days)):
            if rng.random() < trial_spec.dropout.random_rate_per_day:
                dropout_day = float(day)
                break

    completed = dropout_day is None

    # Calculate compliance
    compliance_rate = 1.0
    if trial_spec.compliance is not None:
        compliance_rate = max(0.0, min(1.0,
            rng.gauss(trial_spec.compliance.mean_compliance, trial_spec.compliance.compliance_sd)
        ))

    # Generate endpoint values (simplified)
    endpoint_values = {}
    for endpoint in trial_spec.endpoints:
        if endpoint == "pk_exposure":
            # Simplified exposure calculation
            base_exposure = arm.regimen.dose_amount * compliance_rate
            if not arm.placebo:
                endpoint_values[endpoint] = base_exposure * (1 + rng.gauss(0, 0.2))
            else:
                endpoint_values[endpoint] = rng.gauss(0, 0.1)
        else:
            endpoint_values[endpoint] = rng.gauss(0, 1)

    return SubjectResult(
        subject_id=subject.id,
        arm_name=arm.name,
        completed=completed,
        dropout_day=dropout_day,
        compliance_rate=compliance_rate,
        endpoint_values=endpoint_values,
    )


def _summarize_endpoint(values: List[float]) -> Dict[str, float]:
    """Calculate summary statistics for an endpoint."""
    if not values:
        return {"n": 0, "mean": 0.0, "sd": 0.0, "min": 0.0, "max": 0.0}

    n = len(values)
    mean = sum(values) / n
    sd = math.sqrt(sum((x - mean) ** 2 for x in values) / (n - 1)) if n > 1 else 0.0

    return {
        "n": float(n),
        "mean": mean,
        "sd": sd,
        "min": min(values),
        "max": max(values),
    }


def simulate_trial(
    spec: TrialSpec,
    seed: Optional[int] = None,
) -> TrialResult:
    """
    Simulate a clinical trial.

    Args:
        spec: Trial specification
        seed: Random seed (overrides spec.seed if provided)

    Returns:
        TrialResult object

    Example:
        >>> from neopkpd.trial import (
        ...     TrialSpec, TreatmentArm, parallel_design, dosing_qd
        ... )
        >>> spec = TrialSpec(
        ...     name="Test Trial",
        ...     design=parallel_design(2),
        ...     arms=[
        ...         TreatmentArm("Placebo", dosing_qd(0.0, 28), 25, placebo=True),
        ...         TreatmentArm("Active", dosing_qd(100.0, 28), 25),
        ...     ],
        ...     duration_days=28
        ... )
        >>> result = simulate_trial(spec, seed=42)
        >>> print(f"Completion rate: {result.overall_completion_rate:.2%}")
    """
    if seed is None:
        seed = spec.seed

    rng = random.Random(seed)

    # Generate virtual population
    total_subjects = sum(arm.n_subjects for arm in spec.arms)
    population = generate_virtual_population(
        total_subjects,
        spec=spec.population_spec,
        seed=rng.randint(0, 2**31),
    )

    # Assign subjects to arms
    subject_idx = 0
    arm_results: Dict[str, ArmResult] = {}

    for arm in spec.arms:
        subjects: List[SubjectResult] = []

        for _ in range(arm.n_subjects):
            if subject_idx >= len(population):
                break

            subject = population[subject_idx]
            subject_idx += 1

            result = _simulate_subject(subject, arm, spec, rng)
            subjects.append(result)

        # Calculate arm summaries
        n_completed = sum(1 for s in subjects if s.completed)
        mean_compliance = sum(s.compliance_rate for s in subjects) / len(subjects) if subjects else 0.0

        # Endpoint summaries
        endpoint_summaries: Dict[str, Dict[str, float]] = {}
        for endpoint in spec.endpoints:
            values = [s.endpoint_values.get(endpoint, 0.0) for s in subjects if s.completed]
            endpoint_summaries[endpoint] = _summarize_endpoint(values)

        arm_results[arm.name] = ArmResult(
            name=arm.name,
            n_enrolled=len(subjects),
            n_completed=n_completed,
            completion_rate=n_completed / len(subjects) if subjects else 0.0,
            mean_compliance=mean_compliance,
            subjects=subjects,
            endpoint_summaries=endpoint_summaries,
        )

    # Overall statistics
    total_enrolled = sum(ar.n_enrolled for ar in arm_results.values())
    total_completed = sum(ar.n_completed for ar in arm_results.values())
    overall_completion = total_completed / total_enrolled if total_enrolled > 0 else 0.0

    all_compliance = []
    for ar in arm_results.values():
        all_compliance.extend(s.compliance_rate for s in ar.subjects)
    overall_compliance = sum(all_compliance) / len(all_compliance) if all_compliance else 0.0

    return TrialResult(
        trial_name=spec.name,
        arms=arm_results,
        overall_completion_rate=overall_completion,
        overall_compliance=overall_compliance,
        endpoint_comparisons={},
        replicate=1,
    )


def simulate_trial_replicates(
    spec: TrialSpec,
    n_replicates: Optional[int] = None,
    seed: Optional[int] = None,
) -> List[TrialResult]:
    """
    Simulate multiple replicates of a clinical trial.

    Args:
        spec: Trial specification
        n_replicates: Number of replicates (overrides spec.n_replicates)
        seed: Random seed

    Returns:
        List of TrialResult objects

    Example:
        >>> results = simulate_trial_replicates(spec, n_replicates=100)
        >>> power = sum(1 for r in results if r.endpoint_comparisons.get("significant", False)) / len(results)
    """
    if n_replicates is None:
        n_replicates = spec.n_replicates

    if seed is None:
        seed = spec.seed

    rng = random.Random(seed)
    results = []

    for rep in range(n_replicates):
        rep_seed = rng.randint(0, 2**31)
        result = simulate_trial(spec, seed=rep_seed)
        result.replicate = rep + 1
        results.append(result)

    return results
