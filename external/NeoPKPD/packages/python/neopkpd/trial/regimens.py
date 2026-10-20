"""
NeoPKPD Dosing Regimens

Dosing regimen definitions including:
- Standard regimens (QD, BID, TID, QID)
- Custom timing regimens
- Titration regimens
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, List, Optional, Union
from enum import Enum


class DosingFrequency(str, Enum):
    """Dosing frequency types."""
    QD = "QD"      # Once daily
    BID = "BID"    # Twice daily
    TID = "TID"    # Three times daily
    QID = "QID"    # Four times daily
    CUSTOM = "custom"


# Default dose times (hours after midnight)
DOSE_TIMES = {
    DosingFrequency.QD: [8.0],  # 8 AM
    DosingFrequency.BID: [8.0, 20.0],  # 8 AM, 8 PM
    DosingFrequency.TID: [8.0, 14.0, 20.0],  # 8 AM, 2 PM, 8 PM
    DosingFrequency.QID: [8.0, 12.0, 16.0, 20.0],  # 8 AM, 12 PM, 4 PM, 8 PM
}


@dataclass
class DosingRegimen:
    """
    Dosing regimen specification.

    Attributes:
        frequency: Dosing frequency (QD, BID, TID, QID, or custom)
        dose_amount: Dose amount per administration
        duration_days: Duration of treatment in days
        loading_dose: Optional loading dose for first administration
        dose_times: Times of administration (hours after midnight)

    Example:
        >>> regimen = DosingRegimen(
        ...     frequency=DosingFrequency.BID,
        ...     dose_amount=50.0,
        ...     duration_days=14
        ... )
    """
    frequency: DosingFrequency
    dose_amount: float
    duration_days: int
    loading_dose: Optional[float] = None
    dose_times: List[float] = field(default_factory=list)

    def __post_init__(self):
        if not self.dose_times:
            self.dose_times = DOSE_TIMES.get(self.frequency, [8.0])


@dataclass
class TitrationStep:
    """
    Single step in a titration regimen.

    Attributes:
        dose: Dose amount for this step
        duration_days: Duration at this dose level
    """
    dose: float
    duration_days: int


@dataclass
class TitrationRegimen:
    """
    Titration dosing regimen with gradual dose changes.

    Attributes:
        steps: List of titration steps
        frequency: Dosing frequency
        loading_dose: Optional loading dose

    Example:
        >>> regimen = TitrationRegimen(
        ...     steps=[
        ...         TitrationStep(25.0, 7),
        ...         TitrationStep(50.0, 7),
        ...         TitrationStep(100.0, 14)
        ...     ],
        ...     frequency=DosingFrequency.QD
        ... )
    """
    steps: List[TitrationStep]
    frequency: DosingFrequency = DosingFrequency.QD
    loading_dose: Optional[float] = None


def dosing_qd(
    dose: float,
    duration_days: int,
    loading_dose: Optional[float] = None,
) -> DosingRegimen:
    """
    Create a once-daily dosing regimen.

    Args:
        dose: Daily dose amount
        duration_days: Treatment duration in days
        loading_dose: Optional loading dose on day 1

    Returns:
        DosingRegimen object

    Example:
        >>> regimen = dosing_qd(100.0, 28)
        >>> regimen = dosing_qd(100.0, 28, loading_dose=200.0)
    """
    return DosingRegimen(
        frequency=DosingFrequency.QD,
        dose_amount=dose,
        duration_days=duration_days,
        loading_dose=loading_dose,
    )


def dosing_bid(
    dose: float,
    duration_days: int,
    loading_dose: Optional[float] = None,
) -> DosingRegimen:
    """
    Create a twice-daily dosing regimen.

    Args:
        dose: Dose amount per administration
        duration_days: Treatment duration in days
        loading_dose: Optional loading dose

    Returns:
        DosingRegimen object

    Example:
        >>> regimen = dosing_bid(50.0, 14)  # 50mg twice daily for 14 days
    """
    return DosingRegimen(
        frequency=DosingFrequency.BID,
        dose_amount=dose,
        duration_days=duration_days,
        loading_dose=loading_dose,
    )


def dosing_tid(
    dose: float,
    duration_days: int,
    loading_dose: Optional[float] = None,
) -> DosingRegimen:
    """
    Create a three-times-daily dosing regimen.

    Args:
        dose: Dose amount per administration
        duration_days: Treatment duration in days
        loading_dose: Optional loading dose

    Returns:
        DosingRegimen object

    Example:
        >>> regimen = dosing_tid(25.0, 7)  # 25mg three times daily for 7 days
    """
    return DosingRegimen(
        frequency=DosingFrequency.TID,
        dose_amount=dose,
        duration_days=duration_days,
        loading_dose=loading_dose,
    )


def dosing_qid(
    dose: float,
    duration_days: int,
    loading_dose: Optional[float] = None,
) -> DosingRegimen:
    """
    Create a four-times-daily dosing regimen.

    Args:
        dose: Dose amount per administration
        duration_days: Treatment duration in days
        loading_dose: Optional loading dose

    Returns:
        DosingRegimen object

    Example:
        >>> regimen = dosing_qid(25.0, 5)  # 25mg four times daily for 5 days
    """
    return DosingRegimen(
        frequency=DosingFrequency.QID,
        dose_amount=dose,
        duration_days=duration_days,
        loading_dose=loading_dose,
    )


def dosing_custom(
    dose: float,
    duration_days: int,
    dose_times: List[float],
    loading_dose: Optional[float] = None,
) -> DosingRegimen:
    """
    Create a custom dosing regimen with specified dose times.

    Args:
        dose: Dose amount per administration
        duration_days: Treatment duration in days
        dose_times: Times of administration (hours after midnight)
        loading_dose: Optional loading dose

    Returns:
        DosingRegimen object

    Example:
        >>> # Dose at 6 AM and 10 PM
        >>> regimen = dosing_custom(75.0, 21, [6.0, 22.0])
    """
    return DosingRegimen(
        frequency=DosingFrequency.CUSTOM,
        dose_amount=dose,
        duration_days=duration_days,
        loading_dose=loading_dose,
        dose_times=dose_times,
    )


def titration_regimen(
    start_dose: float,
    target_dose: float,
    n_steps: int,
    days_per_step: int,
    frequency: DosingFrequency = DosingFrequency.QD,
    loading_dose: Optional[float] = None,
    maintenance_days: int = 0,
) -> TitrationRegimen:
    """
    Create a titration dosing regimen.

    Args:
        start_dose: Starting dose
        target_dose: Target maintenance dose
        n_steps: Number of titration steps
        days_per_step: Days at each dose level
        frequency: Dosing frequency
        loading_dose: Optional loading dose
        maintenance_days: Days at maintenance dose (default: 0)

    Returns:
        TitrationRegimen object

    Example:
        >>> # Titrate from 25mg to 100mg in 4 steps, 7 days each
        >>> regimen = titration_regimen(25.0, 100.0, 4, 7)
        >>> # With 28 days maintenance at target dose
        >>> regimen = titration_regimen(25.0, 100.0, 4, 7, maintenance_days=28)
    """
    if n_steps < 2:
        raise ValueError("Titration requires at least 2 steps")

    dose_increment = (target_dose - start_dose) / (n_steps - 1)

    steps = []
    for i in range(n_steps):
        dose = start_dose + i * dose_increment
        steps.append(TitrationStep(dose=dose, duration_days=days_per_step))

    # Add maintenance phase if specified
    if maintenance_days > 0:
        steps.append(TitrationStep(dose=target_dose, duration_days=maintenance_days))

    return TitrationRegimen(
        steps=steps,
        frequency=frequency,
        loading_dose=loading_dose,
    )


def dose_event_times(regimen: Union[DosingRegimen, TitrationRegimen]) -> List[float]:
    """
    Calculate all dose event times for a regimen.

    Args:
        regimen: Dosing regimen

    Returns:
        List of dose times in hours from start

    Example:
        >>> regimen = dosing_bid(50.0, 3)
        >>> times = dose_event_times(regimen)  # [8, 20, 32, 44, 56, 68]
    """
    times = []

    if isinstance(regimen, DosingRegimen):
        dose_times = regimen.dose_times
        for day in range(regimen.duration_days):
            for t in dose_times:
                times.append(day * 24.0 + t)

    elif isinstance(regimen, TitrationRegimen):
        dose_times = DOSE_TIMES.get(regimen.frequency, [8.0])
        current_time = 0.0

        for step in regimen.steps:
            for day in range(step.duration_days):
                for t in dose_times:
                    times.append(current_time + day * 24.0 + t)
            current_time += step.duration_days * 24.0

    return times


def total_regimen_duration(regimen: Union[DosingRegimen, TitrationRegimen]) -> int:
    """
    Calculate total duration of a dosing regimen in days.

    Args:
        regimen: Dosing regimen

    Returns:
        Total duration in days

    Example:
        >>> regimen = dosing_qd(100.0, 7)
        >>> total_regimen_duration(regimen)  # 7
    """
    if isinstance(regimen, DosingRegimen):
        return regimen.duration_days
    elif isinstance(regimen, TitrationRegimen):
        return sum(step.duration_days for step in regimen.steps)
    else:
        raise TypeError(f"Unknown regimen type: {type(regimen)}")


def generate_doses(
    regimen: Union[DosingRegimen, TitrationRegimen],
    compliance_rate: float = 1.0,
    seed: Optional[int] = None,
) -> List[float]:
    """
    Generate dose amounts for each dose event, accounting for compliance.

    Args:
        regimen: Dosing regimen
        compliance_rate: Compliance rate (0-1)
        seed: Random seed for reproducibility

    Returns:
        List of dose amounts (0 for missed doses)

    Example:
        >>> regimen = dosing_qd(100.0, 7)
        >>> doses = generate_doses(regimen)  # [100, 100, 100, 100, 100, 100, 100]
        >>> doses = generate_doses(regimen, compliance_rate=0.8, seed=42)
    """
    import random

    if seed is not None:
        random.seed(seed)

    times = dose_event_times(regimen)
    doses = []

    if isinstance(regimen, DosingRegimen):
        for i, _ in enumerate(times):
            if i == 0 and regimen.loading_dose is not None:
                dose = regimen.loading_dose
            else:
                dose = regimen.dose_amount

            # Apply compliance
            if random.random() <= compliance_rate:
                doses.append(dose)
            else:
                doses.append(0.0)

    elif isinstance(regimen, TitrationRegimen):
        dose_times = DOSE_TIMES.get(regimen.frequency, [8.0])
        doses_per_day = len(dose_times)

        idx = 0
        for step in regimen.steps:
            for _ in range(step.duration_days):
                for _ in range(doses_per_day):
                    dose = step.dose
                    if idx == 0 and regimen.loading_dose is not None:
                        dose = regimen.loading_dose

                    if random.random() <= compliance_rate:
                        doses.append(dose)
                    else:
                        doses.append(0.0)
                    idx += 1

    return doses
