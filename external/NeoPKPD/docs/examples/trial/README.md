# Clinical Trial Simulation Examples

Examples demonstrating clinical trial design and simulation capabilities.

## Examples

| Example | Description | Directory |
|---------|-------------|-----------|
| Parallel Design | Two-arm parallel group study | [01_parallel_design](01_parallel_design/README.md) |
| Crossover Design | 2x2 crossover study | [02_crossover_design](02_crossover_design/README.md) |
| 3+3 Dose Escalation | Phase I dose escalation | [03_dose_escalation_3plus3](03_dose_escalation_3plus3/README.md) |
| Bioequivalence | BE study with TOST analysis | [04_bioequivalence_study](04_bioequivalence_study/README.md) |
| Power Analysis | Sample size and power | [05_power_analysis](05_power_analysis/README.md) |

## Study Design Types

| Design | Use Case | Key Features |
|--------|----------|--------------|
| Parallel | Phase II/III efficacy | Independent arms, randomization |
| Crossover | BE, Phase I | Within-subject comparison, washout |
| Dose Escalation | Phase I FIH | Safety-driven, cohort-based |
| Adaptive | Phase II/III | Interim analyses, alpha spending |

## Key Concepts

### Study Components

1. **Design**: Structure (parallel, crossover, adaptive)
2. **Arms**: Treatment groups with dosing regimens
3. **Population**: Virtual subjects with demographics
4. **Endpoints**: PK exposure, efficacy, safety

### Simulation Features

- Enrollment and randomization
- Dropout modeling
- Compliance patterns
- Multiple replicates

## Usage

```python
from neopkpd import trial

# Simple parallel study
spec = trial.TrialSpec(
    name="Phase 2 Study",
    design=trial.parallel_design(2),
    arms=[
        trial.TreatmentArm("Placebo", trial.dosing_qd(0.0, 28), 50, placebo=True),
        trial.TreatmentArm("Active", trial.dosing_qd(100.0, 28), 50),
    ],
    duration_days=28,
    seed=42
)

result = trial.simulate_trial(spec)
print(f"Completion rate: {result.overall_completion_rate:.1%}")
```

## See Also

- [Population Examples](../population/README.md) - Virtual population generation
- [NCA Examples](../nca/README.md) - Non-compartmental analysis
- [VPC Examples](../vpc/README.md) - Model validation
