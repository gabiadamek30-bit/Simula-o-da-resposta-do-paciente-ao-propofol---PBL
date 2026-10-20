# Study Designs

Comprehensive guide for clinical trial study designs in NeoPKPD Python.

---

## Overview

NeoPKPD supports multiple trial design types for different clinical development phases.

```python
from neopkpd import trial

# Parallel design
design = trial.parallel_design(n_arms=3, randomization_ratio=[1, 1, 1])

# Crossover design
design = trial.crossover_2x2(washout_duration=14.0)

# Dose escalation
design = trial.dose_escalation_3plus3(dose_levels=[10, 25, 50, 100])
```

---

## Parallel Designs

### Basic Parallel Design

```python
from neopkpd.trial import ParallelDesign, parallel_design

# Create 2-arm parallel design
design = parallel_design(
    n_arms=2,
    randomization_ratio=[1, 1]
)

# Access design properties
print(f"Arms: {design.n_arms}")
print(f"Ratio: {design.randomization_ratio}")
```

### Multi-Arm Designs

```python
# 4-arm dose-finding study
design = ParallelDesign(
    n_arms=4,
    arm_names=["Placebo", "Low", "Medium", "High"],
    randomization_ratio=[1, 1, 1, 1],
    stratification_factors=["sex", "age_group"]
)

# Unequal randomization (2:1:1:1)
design = ParallelDesign(
    n_arms=4,
    arm_names=["Placebo", "25mg", "50mg", "100mg"],
    randomization_ratio=[2, 1, 1, 1]
)
```

### ParallelDesign Class

```python
@dataclass
class ParallelDesign:
    """Parallel group trial design."""

    n_arms: int                            # Number of treatment arms
    arm_names: list[str] | None = None     # Names for each arm
    randomization_ratio: list[int] | None = None  # Allocation ratio
    stratification_factors: list[str] | None = None  # Stratification
    block_size: int = 4                    # Block randomization size
```

### Stratified Randomization

```python
# Stratify by sex and age
design = ParallelDesign(
    n_arms=2,
    arm_names=["Control", "Treatment"],
    randomization_ratio=[1, 1],
    stratification_factors=["sex", "age_group"],
    block_size=4
)

# Generate balanced randomization
assignments = trial.generate_randomization(
    design=design,
    n_subjects=100,
    strata_distribution={
        "sex": {"male": 0.5, "female": 0.5},
        "age_group": {"young": 0.3, "middle": 0.4, "elderly": 0.3}
    },
    seed=42
)
```

---

## Crossover Designs

### 2x2 Crossover

```python
from neopkpd.trial import crossover_2x2, CrossoverDesign

# Standard AB/BA crossover
design = crossover_2x2(
    treatments=["Test", "Reference"],
    washout_duration=14.0,
    period_duration=1.0  # 1 day per period
)

print(f"Sequences: {design.sequences}")
# [['Test', 'Reference'], ['Reference', 'Test']]
```

### CrossoverDesign Class

```python
@dataclass
class CrossoverDesign:
    """Crossover trial design."""

    n_periods: int                    # Number of periods
    n_sequences: int                  # Number of sequences
    sequences: list[list[str]]        # Treatment sequences
    washout_duration: float           # Washout between periods (days)
    period_duration: float = 1.0      # Duration of each period (days)
    treatments: list[str] | None = None  # Treatment names
```

### Williams Design

```python
# Balanced for first-order carryover
# For 3 treatments: 6 sequences
design = trial.williams_design(
    treatments=["A", "B", "C"],
    washout_duration=7.0
)

print(f"Sequences: {design.sequences}")
# [['A', 'B', 'C'], ['B', 'C', 'A'], ['C', 'A', 'B'],
#  ['A', 'C', 'B'], ['C', 'B', 'A'], ['B', 'A', 'C']]
```

### Latin Square Design

```python
# 3×3 Latin square
design = trial.latin_square_design(
    treatments=["A", "B", "C"],
    washout_duration=7.0
)

# Sequences: ABC, BCA, CAB
```

### Replicate Designs

```python
# Full replicate (TRTR/RTRT)
design = trial.replicate_crossover_2x4(
    treatments=["Test", "Reference"],
    washout_duration=7.0
)
# Sequences: TRTR, RTRT

# Partial replicate (TRR/RTR/RRT)
design = trial.partial_replicate_3x3(
    treatments=["Test", "Reference"],
    washout_duration=7.0
)
```

---

## Dose Escalation Designs

### 3+3 Design

```python
from neopkpd.trial import DoseEscalation3plus3

# Standard 3+3
design = DoseEscalation3plus3(
    dose_levels=[10.0, 25.0, 50.0, 100.0, 200.0, 400.0],
    starting_dose_index=0,
    target_dlt_rate=0.33
)

# Convenience function
design = trial.dose_escalation_3plus3(
    dose_levels=[10, 25, 50, 100, 200],
    starting_dose=10.0
)
```

### DoseEscalation3plus3 Class

```python
@dataclass
class DoseEscalation3plus3:
    """Rule-based 3+3 dose escalation design."""

    dose_levels: list[float]         # Available doses
    starting_dose_index: int = 0     # Starting dose (0-indexed)
    target_dlt_rate: float = 0.33    # Target toxicity rate
    max_patients_per_dose: int = 6   # Maximum per cohort

    def next_dose_decision(
        self,
        n_dlt: int,
        n_patients: int
    ) -> str:
        """Return decision: 'escalate', 'expand', 'deescalate', 'stop'"""
        ...
```

### 3+3 Decision Rules

```python
# Decision logic
def get_3plus3_decision(n_dlt: int, n_patients: int) -> str:
    """
    3+3 decision rules:
    - 0/3 DLT → Escalate
    - 1/3 DLT → Expand to 6
    - 2-3/3 DLT → De-escalate
    - 0-1/6 DLT → Escalate
    - 2+/6 DLT → MTD = previous dose
    """
    if n_patients == 3:
        if n_dlt == 0:
            return "escalate"
        elif n_dlt == 1:
            return "expand"
        else:
            return "deescalate"
    elif n_patients == 6:
        if n_dlt <= 1:
            return "escalate"
        else:
            return "deescalate"
    return "stay"
```

### mTPI Design

```python
from neopkpd.trial import MTPI

design = MTPI(
    dose_levels=[10.0, 25.0, 50.0, 100.0, 200.0],
    target_toxicity=0.25,
    epsilon1=0.05,              # Lower equivalence margin
    epsilon2=0.05,              # Upper equivalence margin
    prior_alpha=1.0,            # Beta prior
    prior_beta=1.0,
    starting_dose_index=0,
    max_sample_size=36
)
```

### CRM Design

```python
from neopkpd.trial import CRM

design = CRM(
    dose_levels=[10.0, 25.0, 50.0, 100.0, 200.0],
    target_toxicity=0.25,
    skeleton=[0.05, 0.10, 0.25, 0.40, 0.55],  # Prior p(DLT)
    model="power",              # or "logistic"
    prior_mean=0.0,
    prior_sd=1.34,
    cohort_size=1,
    max_sample_size=30
)

# Get next recommended dose
next_dose_idx = design.recommend_dose(
    dose_history=[0, 0, 0, 1, 1, 1],  # Dose indices
    dlt_history=[False, False, False, False, True, False]
)
```

### BOIN Design

```python
from neopkpd.trial import BOIN

design = BOIN(
    dose_levels=[10.0, 25.0, 50.0, 100.0, 200.0],
    target_toxicity=0.30,
    p1=0.25,                    # Highest acceptable rate
    p2=0.35,                    # Lowest unacceptable rate
    cohort_size=3,
    max_sample_size=36
)

# BOIN boundaries
print(f"Escalation boundary: {design.lambda_e:.3f}")
print(f"De-escalation boundary: {design.lambda_d:.3f}")
```

---

## Bioequivalence Designs

### Standard BE Design

```python
from neopkpd.trial import BioequivalenceDesign

design = BioequivalenceDesign(
    n_periods=2,
    test_formulation="Test",
    reference_formulation="Reference",
    washout_duration=14.0,
    theta1=0.80,                # Lower BE limit
    theta2=1.25,                # Upper BE limit
    alpha=0.05
)
```

### BE Design Types

```python
# 2×2 crossover (standard)
design = trial.be_design_2x2(washout_duration=14.0)

# 2×4 replicate (for HVD)
design = trial.be_design_2x4(washout_duration=7.0)

# 3×3 partial replicate
design = trial.be_design_3x3_partial(washout_duration=7.0)
```

### Highly Variable Drug Designs

```python
# RSABE design (FDA)
design = trial.rsabe_design(
    reference_cv=0.40,          # Reference CV > 30%
    washout_duration=7.0,
    theta_s=0.8928              # Scaling factor
)

# ABEL design (EMA)
design = trial.abel_design(
    reference_cv=0.40,
    washout_duration=7.0,
    max_widening=0.50           # Max limit widening
)
```

---

## Adaptive Designs

### Group Sequential Design

```python
from neopkpd.trial import GroupSequentialDesign

design = GroupSequentialDesign(
    n_looks=3,                   # Number of interim analyses
    information_fractions=[0.33, 0.67, 1.0],
    alpha=0.05,
    power=0.80,
    spending_function="obrien_fleming"
)

# Get boundaries
boundaries = design.get_efficacy_boundaries()
print(f"Stage 1 boundary: {boundaries[0]:.3f}")
print(f"Stage 2 boundary: {boundaries[1]:.3f}")
print(f"Final boundary: {boundaries[2]:.3f}")
```

### Sample Size Re-estimation

```python
from neopkpd.trial import AdaptiveSampleSize

design = AdaptiveSampleSize(
    initial_n=50,
    interim_fraction=0.5,
    conditional_power_target=0.80,
    max_n=150,
    min_increase=10
)

# Re-estimate at interim
new_n = design.reestimate_sample_size(
    observed_effect=0.35,
    observed_se=0.15,
    n_current=50
)
```

---

## Design Validation

### Check Design Properties

```python
# Validate crossover design
design = trial.crossover_2x2(washout_duration=14.0)

validation = trial.validate_design(design)
print(f"Valid: {validation.is_valid}")
print(f"Balanced: {validation.is_balanced}")
print(f"Warnings: {validation.warnings}")
```

### Design Summary

```python
# Print design summary
summary = trial.design_summary(design)
print(summary)

# Example output:
# Study Design: 2×2 Crossover
# Treatments: Test, Reference
# Periods: 2
# Sequences: TR, RT
# Washout: 14.0 days
# Balanced: Yes
```

---

## Complete Example

```python
from neopkpd import trial

# =====================================
# Phase III Parallel Design Setup
# =====================================

# 1. Create design
design = trial.ParallelDesign(
    n_arms=3,
    arm_names=["Placebo", "50mg", "100mg"],
    randomization_ratio=[1, 2, 2],  # More on active arms
    stratification_factors=["sex", "age_group"],
    block_size=5
)

# 2. Define dosing regimens for each arm
regimens = {
    "Placebo": trial.dosing_qd(dose=0.0, duration_days=84),
    "50mg": trial.dosing_qd(dose=50.0, duration_days=84),
    "100mg": trial.dosing_qd(dose=100.0, duration_days=84)
}

# 3. Generate population
pop = trial.generate_virtual_population(
    n=250,
    spec=trial.patient_population_spec(),
    seed=42
)

# 4. Define trial specification
spec = trial.TrialSpec(
    name="Phase 3 Efficacy Trial",
    design=design,
    arms=[
        trial.TreatmentArm(name, regimen=reg)
        for name, reg in regimens.items()
    ],
    population=pop,
    dropout=trial.DropoutSpec(rate=0.15, pattern="exponential"),
    compliance=trial.ComplianceSpec(mean=0.85, sd=0.10),
    endpoints=["auc_ss", "cmax_ss", "cmin_ss"]
)

# 5. Print summary
print("=== Trial Design Summary ===")
print(f"Design: {design.n_arms}-arm parallel")
print(f"Randomization: {design.randomization_ratio}")
print(f"Stratification: {design.stratification_factors}")
print(f"Total subjects: {len(pop)}")
print(f"Duration: 84 days")
print(f"Expected completers: {int(len(pop) * 0.85)}")
```

---

## See Also

- [Dosing Regimens](dosing.md) - Dosing schedule configuration
- [Virtual Population](population.md) - Population generation
- [Power Analysis](power.md) - Sample size calculation
- [Julia Trial Designs](../../julia/trial/index.md) - Julia interface
