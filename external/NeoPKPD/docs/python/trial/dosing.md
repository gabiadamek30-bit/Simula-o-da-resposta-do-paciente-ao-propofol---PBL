# Dosing Regimens

Comprehensive guide for configuring dosing schedules in NeoPKPD Python.

---

## Overview

The dosing module provides flexible specification of drug administration schedules.

```python
from neopkpd import trial

# Once daily dosing
regimen = trial.dosing_qd(dose=100.0, duration_days=28)

# Custom schedule
regimen = trial.dosing_custom(
    dose=100.0,
    times_per_day=[0.0, 8.0, 16.0],
    duration_days=14
)
```

---

## Standard Regimens

### Once Daily (QD)

```python
from neopkpd.trial import dosing_qd

# Basic QD dosing
regimen = dosing_qd(
    dose=100.0,           # Dose amount
    duration_days=28,     # Treatment duration
    dose_time=8.0         # Time of day (hours, default 8 AM)
)

# With loading dose
regimen = dosing_qd(
    dose=100.0,
    duration_days=28,
    loading_dose=200.0    # Double first dose
)

# Access schedule
print(f"Doses per day: {regimen.doses_per_day}")
print(f"Total doses: {regimen.total_doses}")
print(f"Dose times: {regimen.dose_times}")
```

### Twice Daily (BID)

```python
from neopkpd.trial import dosing_bid

# Standard BID
regimen = dosing_bid(
    dose=50.0,
    duration_days=14,
    dose_times=[8.0, 20.0]  # 8 AM and 8 PM
)

# With morning loading
regimen = dosing_bid(
    dose=50.0,
    duration_days=14,
    loading_dose=100.0,
    loading_applies_to="first_only"  # or "morning_only"
)
```

### Three Times Daily (TID)

```python
from neopkpd.trial import dosing_tid

# Standard TID (every 8 hours)
regimen = dosing_tid(
    dose=25.0,
    duration_days=7,
    dose_times=[6.0, 14.0, 22.0]  # 6 AM, 2 PM, 10 PM
)
```

### Four Times Daily (QID)

```python
from neopkpd.trial import dosing_qid

# Standard QID (every 6 hours)
regimen = dosing_qid(
    dose=20.0,
    duration_days=7,
    dose_times=[6.0, 12.0, 18.0, 24.0]
)
```

---

## DosingRegimen Class

### Class Definition

```python
@dataclass
class DosingRegimen:
    """Complete dosing regimen specification."""

    dose: float                      # Standard dose amount
    duration_days: int               # Total treatment duration
    dose_times: list[float]          # Times of day for doses (hours)
    loading_dose: float | None = None  # Optional loading dose
    route: str = "oral"              # Administration route
    formulation: str = "tablet"      # Drug formulation

    @property
    def doses_per_day(self) -> int:
        """Number of doses per day."""
        return len(self.dose_times)

    @property
    def total_doses(self) -> int:
        """Total doses over regimen."""
        return self.doses_per_day * self.duration_days

    def get_dose_schedule(self) -> list[dict]:
        """Return complete dose schedule."""
        ...
```

### Creating Custom Regimens

```python
from neopkpd.trial import DosingRegimen

# Manual specification
regimen = DosingRegimen(
    dose=100.0,
    duration_days=28,
    dose_times=[8.0],           # Once at 8 AM
    loading_dose=200.0,
    route="oral",
    formulation="tablet"
)

# Get schedule as list
schedule = regimen.get_dose_schedule()
for event in schedule[:5]:
    print(f"Day {event['day']}, Time {event['time']}: {event['dose']} mg")
```

---

## Custom Schedules

### Irregular Dosing

```python
from neopkpd.trial import dosing_custom

# Custom times
regimen = dosing_custom(
    dose=100.0,
    times_per_day=[7.0, 13.0, 19.0],  # 7 AM, 1 PM, 7 PM
    duration_days=14
)

# Variable doses per day
regimen = dosing_custom(
    doses=[100.0, 50.0, 50.0],  # Different amounts
    times_per_day=[8.0, 14.0, 20.0],
    duration_days=14
)
```

### Weekly Dosing

```python
# Once weekly
regimen = trial.dosing_weekly(
    dose=500.0,
    duration_weeks=12,
    dose_day=1  # Monday (1=Mon, 7=Sun)
)

# Twice weekly
regimen = trial.dosing_twice_weekly(
    dose=250.0,
    duration_weeks=12,
    dose_days=[1, 4]  # Monday and Thursday
)
```

### PRN (As Needed)

```python
# PRN dosing with max daily dose
regimen = trial.dosing_prn(
    dose=50.0,
    max_doses_per_day=4,
    duration_days=14,
    min_interval_hours=4.0
)
```

---

## Titration Regimens

### Linear Titration

```python
from neopkpd.trial import titration_regimen

# Gradual dose increase
regimen = titration_regimen(
    start_dose=25.0,
    target_dose=100.0,
    steps=[25, 50, 75, 100],
    days_per_step=7,
    frequency="qd",
    maintenance_days=28
)

# Schedule:
# Days 1-7:   25 mg QD
# Days 8-14:  50 mg QD
# Days 15-21: 75 mg QD
# Days 22+:   100 mg QD (maintenance)
```

### TitrationRegimen Class

```python
@dataclass
class TitrationRegimen:
    """Dose titration regimen."""

    start_dose: float               # Initial dose
    target_dose: float              # Target maintenance dose
    steps: list[float]              # Dose levels
    days_per_step: int              # Days at each level
    frequency: str = "qd"           # Dosing frequency
    maintenance_days: int = 0       # Days at target dose
    back_titration_allowed: bool = False  # Can decrease?

    @property
    def total_titration_days(self) -> int:
        return len(self.steps) * self.days_per_step

    @property
    def total_duration(self) -> int:
        return self.total_titration_days + self.maintenance_days
```

### Flexible Titration

```python
# With tolerability-based adjustment
regimen = trial.flexible_titration(
    start_dose=25.0,
    target_dose=100.0,
    step_size=25.0,
    min_days_per_step=3,
    max_days_per_step=14,
    tolerability_criterion="ae_grade < 2"
)
```

---

## Infusion Regimens

### IV Infusion

```python
from neopkpd.trial import dosing_infusion

# Short infusion
regimen = dosing_infusion(
    dose=500.0,              # mg
    infusion_duration=1.0,   # hours
    frequency="qd",
    duration_days=5
)

# Long infusion
regimen = dosing_infusion(
    dose=1000.0,
    infusion_duration=24.0,  # Continuous over 24h
    frequency="qd",
    duration_days=7
)
```

### IV Bolus

```python
# Bolus injection
regimen = trial.dosing_iv_bolus(
    dose=100.0,
    frequency="qd",
    duration_days=3
)
```

### Loading + Maintenance Infusion

```python
# Loading bolus followed by infusion
regimen = trial.dosing_loading_infusion(
    loading_dose=500.0,
    loading_duration=0.5,     # 30-min loading
    maintenance_rate=50.0,    # mg/hour
    maintenance_duration=24.0,
    duration_days=5
)
```

---

## Multiple Formulations

### Formulation Specification

```python
from neopkpd.trial import FormulationSpec

# Define formulations
tablet = FormulationSpec(
    name="tablet",
    route="oral",
    bioavailability=0.80,
    absorption_rate=1.5  # Ka
)

solution = FormulationSpec(
    name="solution",
    route="oral",
    bioavailability=0.95,
    absorption_rate=2.5
)

# Use in regimen
regimen = trial.dosing_qd(
    dose=100.0,
    duration_days=28,
    formulation=tablet
)
```

### Switching Formulations

```python
# Switch from IV to oral
regimen = trial.sequential_formulation(
    phase1=trial.dosing_infusion(dose=500.0, duration_days=3),
    phase2=trial.dosing_qd(dose=250.0, duration_days=25)
)
```

---

## Compliance Modeling

### ComplianceSpec Class

```python
@dataclass
class ComplianceSpec:
    """Patient compliance specification."""

    mean: float = 0.90          # Mean compliance rate
    sd: float = 0.10            # Standard deviation
    pattern: str = "random"     # "random", "decay", "weekend_miss"
    min_compliance: float = 0.50  # Minimum allowed

    def sample_compliance(self, n: int, seed: int = None) -> list[float]:
        """Generate individual compliance rates."""
        ...
```

### Compliance Patterns

```python
# Random missing doses
compliance = trial.ComplianceSpec(
    mean=0.85,
    sd=0.10,
    pattern="random"
)

# Weekend-miss pattern
compliance = trial.ComplianceSpec(
    mean=0.90,
    pattern="weekend_miss",
    weekday_rate=0.95,
    weekend_rate=0.70
)

# Declining compliance
compliance = trial.ComplianceSpec(
    mean=0.85,
    pattern="decay",
    initial_rate=0.95,
    final_rate=0.75,
    half_time_days=14
)
```

### Apply Compliance

```python
# Apply to regimen
actual_doses = trial.apply_compliance(
    regimen=regimen,
    compliance_spec=compliance,
    n_subjects=100,
    seed=42
)

# Each subject gets individual dose schedule
for subject_id, doses in actual_doses.items():
    taken = sum(1 for d in doses if d["taken"])
    total = len(doses)
    print(f"Subject {subject_id}: {taken}/{total} doses taken")
```

---

## Dose Modifications

### Dose Reduction Rules

```python
from neopkpd.trial import DoseModificationRule

# Reduce for toxicity
rule = DoseModificationRule(
    trigger="ae_grade >= 3",
    action="reduce_by_25%",
    min_dose=25.0,
    max_reductions=2
)

# Skip dose for lab value
rule = DoseModificationRule(
    trigger="neutrophil_count < 1000",
    action="hold_until_recovery",
    recovery_criterion="neutrophil_count >= 1500"
)
```

### Applying Modifications

```python
# Add rules to regimen
regimen = trial.dosing_qd(
    dose=100.0,
    duration_days=28,
    modification_rules=[
        DoseModificationRule(trigger="ae_grade >= 3", action="reduce_by_25%"),
        DoseModificationRule(trigger="ae_grade >= 4", action="discontinue")
    ]
)
```

---

## Regimen Validation

### Check Regimen

```python
# Validate regimen
validation = trial.validate_regimen(regimen)

print(f"Valid: {validation.is_valid}")
print(f"Warnings: {validation.warnings}")
print(f"Daily dose: {validation.daily_dose}")
print(f"Total exposure: {validation.total_dose}")
```

### Comparison

```python
# Compare regimens
comparison = trial.compare_regimens(regimen1, regimen2)

print(f"Dose ratio: {comparison.dose_ratio}")
print(f"Frequency difference: {comparison.frequency_diff}")
print(f"Duration difference: {comparison.duration_diff}")
```

---

## Complete Example

```python
from neopkpd import trial

# ================================
# Complex Dosing Schedule Setup
# ================================

# 1. Define titration for new patients
titration = trial.titration_regimen(
    start_dose=25.0,
    target_dose=100.0,
    steps=[25, 50, 75, 100],
    days_per_step=7,
    frequency="qd"
)

# 2. Define maintenance regimen
maintenance = trial.dosing_bid(
    dose=50.0,
    duration_days=56,
    dose_times=[8.0, 20.0]
)

# 3. Combine into full treatment
full_regimen = trial.sequential_regimen([
    titration,
    maintenance
])

# 4. Define compliance
compliance = trial.ComplianceSpec(
    mean=0.90,
    sd=0.08,
    pattern="decay",
    initial_rate=0.95,
    final_rate=0.80,
    half_time_days=28
)

# 5. Add modification rules
full_regimen = trial.add_modification_rules(
    full_regimen,
    rules=[
        trial.DoseModificationRule(
            trigger="ae_grade >= 3",
            action="reduce_by_25%"
        ),
        trial.DoseModificationRule(
            trigger="discontinuation_criterion",
            action="discontinue"
        )
    ]
)

# 6. Print schedule summary
print("=== Dosing Schedule Summary ===")
print(f"Titration phase: {titration.total_duration} days")
print(f"Maintenance phase: {maintenance.duration_days} days")
print(f"Total duration: {full_regimen.total_duration} days")
print(f"Total doses: {full_regimen.total_doses}")

# 7. Generate individual schedules
schedules = trial.generate_individual_schedules(
    regimen=full_regimen,
    compliance=compliance,
    n_subjects=50,
    seed=42
)

# 8. Calculate actual exposure
for i, schedule in enumerate(schedules[:5]):
    doses_taken = sum(1 for d in schedule if d["taken"])
    total_dose = sum(d["dose"] for d in schedule if d["taken"])
    print(f"Subject {i+1}: {doses_taken} doses, {total_dose:.0f} mg total")
```

---

## Regimen Functions Reference

| Function | Description |
|----------|-------------|
| `dosing_qd` | Once daily dosing |
| `dosing_bid` | Twice daily dosing |
| `dosing_tid` | Three times daily |
| `dosing_qid` | Four times daily |
| `dosing_weekly` | Once weekly |
| `dosing_custom` | Custom schedule |
| `dosing_infusion` | IV infusion |
| `titration_regimen` | Dose escalation |
| `sequential_regimen` | Combine regimens |
| `apply_compliance` | Add missed doses |
| `validate_regimen` | Check regimen validity |

---

## See Also

- [Study Designs](designs.md) - Trial design types
- [Virtual Population](population.md) - Subject generation
- [Power Analysis](power.md) - Sample size calculation
- [Julia Dosing](../../julia/trial/parallel.md) - Julia interface
