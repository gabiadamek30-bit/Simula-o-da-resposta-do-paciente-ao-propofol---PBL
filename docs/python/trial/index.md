# Clinical Trials

The `neopkpd.trial` module provides comprehensive clinical trial simulation and analysis capabilities.

---

## Overview

<div class="grid cards" markdown>

-   :material-clipboard-list:{ .lg .middle } **Study Designs**

    ---

    Parallel, crossover, escalation

    [:octicons-arrow-right-24: Designs](designs.md)

-   :material-pill:{ .lg .middle } **Dosing Regimens**

    ---

    QD, BID, custom schedules

    [:octicons-arrow-right-24: Dosing](dosing.md)

-   :material-account-group:{ .lg .middle } **Virtual Population**

    ---

    Generate realistic subjects

    [:octicons-arrow-right-24: Population](population.md)

-   :material-chart-line:{ .lg .middle } **Power Analysis**

    ---

    Sample size determination

    [:octicons-arrow-right-24: Power](power.md)

</div>

---

## Quick Start

### Power Analysis

```python
from neopkpd import trial

# Calculate power
power = trial.estimate_power_analytical(
    n_per_arm=50,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05
)
print(f"Power: {power.power:.1%}")

# Calculate required sample size
result = trial.estimate_sample_size(
    target_power=0.80,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05
)
print(f"Required n per arm: {result.n_per_arm}")
```

### Generate Virtual Population

```python
# Default healthy volunteers
pop = trial.generate_virtual_population(
    n=100,
    spec=trial.healthy_volunteer_spec(),
    seed=42
)

# Summarize demographics
summary = trial.summarize_population(pop)
print(f"Age: {summary['age']['mean']:.1f} years")
print(f"Weight: {summary['weight']['mean']:.1f} kg")
print(f"Female: {summary['female_fraction']*100:.0f}%")
```

---

## Study Designs

### Parallel Design

```python
design = trial.parallel_design(
    n_arms=3,
    randomization_ratio=[1, 1, 1]
)
```

### Crossover Design

```python
# 2×2 crossover
design = trial.crossover_2x2(washout_duration=14.0)

# Williams design (3-period)
design = trial.williams_design(washout_duration=7.0)
```

### Dose Escalation

```python
# 3+3 design
design = trial.dose_escalation_3plus3(
    starting_dose=10.0,
    dose_levels=[10, 25, 50, 100, 200]
)
```

### Bioequivalence

```python
design = trial.bioequivalence_design(
    n_periods=2,
    washout_duration=14.0,
    reference_formulation="tablet",
    test_formulation="capsule"
)
```

---

## Dosing Regimens

```python
# Once daily
regimen = trial.dosing_qd(dose=100.0, duration_days=28)

# Twice daily
regimen = trial.dosing_bid(dose=50.0, duration_days=14)

# Three times daily
regimen = trial.dosing_tid(dose=25.0, duration_days=7)

# Custom schedule
regimen = trial.dosing_custom(
    dose=100.0,
    times_per_day=[0.0, 8.0, 16.0],
    duration_days=14
)

# Titration
regimen = trial.titration_regimen(
    start_dose=25.0,
    target_dose=100.0,
    steps=[25, 50, 75, 100],
    days_per_step=7
)
```

---

## Virtual Population

### Demographics Specification

```python
spec = trial.DemographicSpec(
    age_mean=55.0,
    age_sd=12.0,
    age_min=18.0,
    age_max=80.0,
    weight_mean=85.0,
    weight_sd=18.0,
    weight_min=50.0,
    weight_max=150.0,
    female_fraction=0.45,
    race_distribution={
        "white": 0.6,
        "black": 0.2,
        "asian": 0.15,
        "other": 0.05
    }
)

pop = trial.generate_virtual_population(n=200, spec=spec, seed=42)
```

### Built-in Populations

```python
# Healthy volunteers
spec = trial.healthy_volunteer_spec()

# Elderly patients
spec = trial.elderly_patient_spec()

# Pediatric
spec = trial.pediatric_spec()
```

---

## Trial Simulation

```python
# Define trial
spec = trial.TrialSpec(
    name="Phase 2 Dose Finding",
    design=trial.parallel_design(3),
    arms=[
        trial.TreatmentArm("Placebo", regimen=trial.dosing_qd(0.0, 28)),
        trial.TreatmentArm("Low Dose", regimen=trial.dosing_qd(50.0, 28)),
        trial.TreatmentArm("High Dose", regimen=trial.dosing_qd(100.0, 28)),
    ],
    population=trial.generate_virtual_population(150, seed=42),
    dropout=trial.DropoutSpec(rate=0.05, pattern="exponential"),
    compliance=trial.ComplianceSpec(mean=0.90, sd=0.10)
)

# Run simulation
result = trial.simulate_trial(spec, seed=12345)

# Analyze results
for arm_name, arm_result in result.arms.items():
    print(f"{arm_name}: n={arm_result.n_completed}")
```

---

## Statistical Analysis

### Arm Comparison

```python
comparison = trial.compare_arms(
    treatment_values=[1.2, 1.5, 1.1, 1.8, 1.4],
    control_values=[0.9, 1.0, 0.8, 1.1, 0.95],
    test="ttest"
)

print(f"Difference: {comparison.difference:.2f}")
print(f"95% CI: ({comparison.ci_lower:.2f}, {comparison.ci_upper:.2f})")
print(f"p-value: {comparison.p_value:.4f}")
```

### Responder Analysis

```python
result = trial.responder_analysis(
    values=[1.2, 0.8, 1.5, 0.6, 1.1, 2.0],
    threshold=1.0
)
print(f"Responder rate: {result.rate:.1%}")
```

---

## Power Functions

| Function | Description |
|----------|-------------|
| `estimate_power_analytical` | Analytical power calculation |
| `estimate_sample_size` | Sample size for target power |
| `alpha_spending_function` | Interim analysis alpha spending |

### Power Calculation

```python
power = trial.estimate_power_analytical(
    n_per_arm=50,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05,
    test="two_sample_t"
)
```

### Alpha Spending

```python
alpha = trial.alpha_spending_function(
    information_fraction=0.5,
    total_alpha=0.05,
    method="obrien_fleming"
)
```

---

## Next Steps

- [Study Designs](designs.md) - Detailed design options
- [Power Analysis](power.md) - Sample size determination
- [Trial Visualization](../viz/trial.md) - Power curves, Kaplan-Meier
