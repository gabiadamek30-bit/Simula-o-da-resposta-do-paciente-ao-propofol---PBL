# Power Analysis

Comprehensive guide for sample size calculation and power analysis in NeoPKPD Python.

---

## Overview

Power analysis determines the sample size needed to detect effects with adequate statistical power.

```python
from neopkpd import trial

# Calculate required sample size
result = trial.estimate_sample_size(
    target_power=0.80,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05
)

print(f"Required n per arm: {result.n_per_arm}")
```

---

## Analytical Power

### Two-Sample t-Test

```python
from neopkpd.trial import estimate_power_analytical

# Calculate power for given sample size
power = estimate_power_analytical(
    n_per_arm=50,
    effect_size=0.5,        # Cohen's d
    sd=1.0,
    alpha=0.05,
    test="two_sample_t"
)

print(f"Power: {power.power:.1%}")
print(f"Effect size: {power.effect_size}")
```

### One-Sample t-Test

```python
power = estimate_power_analytical(
    n=30,
    effect_size=0.6,
    sd=1.0,
    alpha=0.05,
    test="one_sample_t"
)
```

### Paired t-Test

```python
power = estimate_power_analytical(
    n_pairs=25,
    effect_size=0.5,
    sd=1.0,                 # SD of differences
    alpha=0.05,
    test="paired_t"
)
```

---

## Sample Size Estimation

### Basic Estimation

```python
from neopkpd.trial import estimate_sample_size

# Find n for target power
result = estimate_sample_size(
    target_power=0.80,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05,
    test="two_sample_t"
)

print(f"Required n per arm: {result.n_per_arm}")
print(f"Total N: {result.n_total}")
print(f"Achieved power: {result.achieved_power:.1%}")
```

### SampleSizeResult Class

```python
@dataclass
class SampleSizeResult:
    """Result of sample size calculation."""

    n_per_arm: int              # Sample size per arm
    n_total: int                # Total sample size
    achieved_power: float       # Actual power achieved
    effect_size: float          # Effect size used
    alpha: float                # Alpha level
    test: str                   # Statistical test used
```

### With Dropout Adjustment

```python
# Account for anticipated dropout
result = estimate_sample_size(
    target_power=0.80,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05,
    dropout_rate=0.15           # 15% dropout
)

print(f"Enrollment needed: {result.n_to_enroll}")
print(f"Expected completers: {result.n_per_arm}")
```

---

## Bioequivalence Power

### Standard BE (2×2 Crossover)

```python
from neopkpd.trial import be_sample_size

# Calculate sample size for BE study
result = be_sample_size(
    gmr=0.95,                   # Expected geometric mean ratio
    cv_within=0.25,             # Within-subject CV
    theta1=0.80,                # Lower BE limit
    theta2=1.25,                # Upper BE limit
    alpha=0.05,
    power=0.80
)

print(f"Required N: {result.n_total}")
print(f"Power: {result.power:.1%}")
```

### BE Power Calculation

```python
from neopkpd.trial import be_power

# Calculate power for given sample size
power = be_power(
    n=24,
    gmr=0.95,
    cv_within=0.25,
    theta1=0.80,
    theta2=1.25,
    alpha=0.05
)

print(f"Power: {power:.1%}")
```

### Replicate Design Sample Size

```python
# 2×4 replicate design (for HVD)
result = be_sample_size(
    gmr=0.95,
    cv_within=0.35,             # Higher CV
    theta1=0.80,
    theta2=1.25,
    alpha=0.05,
    power=0.80,
    design="replicate_2x4"
)

# Partial replicate (3×3)
result = be_sample_size(
    gmr=0.95,
    cv_within=0.35,
    power=0.80,
    design="partial_replicate"
)
```

---

## Highly Variable Drug Power

### RSABE Sample Size (FDA)

```python
from neopkpd.trial import rsabe_sample_size

# Reference-scaled average bioequivalence
result = rsabe_sample_size(
    gmr=0.95,
    cv_reference=0.40,          # Reference CV > 30%
    theta_s=0.8928,             # Scaling factor
    sigma_w0=0.25,              # Regulatory cutoff
    alpha=0.05,
    power=0.80
)

print(f"Required N: {result.n_total}")
print(f"Scaling applied: {result.scaling_applied}")
print(f"Effective limits: [{result.lower_limit:.2f}, {result.upper_limit:.2f}]")
```

### ABEL Sample Size (EMA)

```python
from neopkpd.trial import abel_sample_size

# Average bioequivalence with expanding limits
result = abel_sample_size(
    gmr=0.95,
    cv_reference=0.40,
    cv_cutoff=0.30,
    max_widening=0.50,          # Maximum expansion
    alpha=0.05,
    power=0.80
)

print(f"Required N: {result.n_total}")
print(f"Widened limits: [{result.lower_limit:.2%}, {result.upper_limit:.2%}]")
```

---

## Power Curves

### Generate Power Curve

```python
from neopkpd.trial import power_curve

# Power vs sample size
curve = power_curve(
    effect_size=0.5,
    sd=1.0,
    alpha=0.05,
    n_range=range(10, 101, 5),
    test="two_sample_t"
)

# Find minimum N for 80% power
min_n = next(n for n, p in zip(curve.n, curve.power) if p >= 0.80)
print(f"Minimum N for 80% power: {min_n}")
```

### Effect Size Sensitivity

```python
# Power for different effect sizes at fixed N
sensitivity = trial.effect_sensitivity(
    n=50,
    alpha=0.05,
    effect_sizes=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    test="two_sample_t"
)

for es, pow in zip(sensitivity.effect_sizes, sensitivity.powers):
    print(f"Effect size {es:.1f}: Power = {pow:.1%}")
```

### Minimum Detectable Effect

```python
# What effect can we detect with given N and power?
mde = trial.minimum_detectable_effect(
    n_per_arm=50,
    alpha=0.05,
    power=0.80,
    sd=1.0
)

print(f"Minimum detectable effect: d = {mde:.3f}")
```

---

## Simulation-Based Power

### Monte Carlo Power Estimation

```python
from neopkpd.trial import simulate_power

# Power via simulation
result = simulate_power(
    n_per_arm=50,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05,
    n_simulations=1000,
    seed=42
)

print(f"Simulated power: {result.power:.1%}")
print(f"95% CI: ({result.ci_lower:.1%}, {result.ci_upper:.1%})")
```

### BE Study Simulation

```python
# Simulate BE study power
result = trial.simulate_be_power(
    n=24,
    gmr=0.95,
    cv_within=0.25,
    n_simulations=1000,
    seed=42
)

print(f"BE power: {result.power:.1%}")
print(f"Pass rate: {result.pass_rate:.1%}")
```

### Trial Simulation Power

```python
# Full trial simulation
power = trial.simulate_trial_power(
    trial_spec=trial_spec,
    n_per_arm=50,
    endpoint="auc_comparison",
    success_criterion=lambda r: r.p_value < 0.05 and r.effect > 0,
    n_simulations=1000,
    seed=42
)
```

---

## Multi-Arm Trials

### Dunnett's Test

```python
from neopkpd.trial import multiarm_sample_size

# Multiple treatment arms vs placebo
result = multiarm_sample_size(
    n_arms=4,                   # 1 placebo + 3 active
    effect_size=0.5,
    sd=1.0,
    alpha=0.05,
    power=0.80,
    comparison="dunnett"        # Many-to-one
)

print(f"N per arm: {result.n_per_arm}")
print(f"Total N: {result.n_total}")
print(f"Adjusted alpha: {result.adjusted_alpha:.4f}")
```

### ANOVA

```python
# Overall group difference
result = multiarm_sample_size(
    n_arms=4,
    effect_size=0.3,            # f-statistic
    alpha=0.05,
    power=0.80,
    comparison="anova"
)
```

### Pairwise Comparisons

```python
# All pairwise with Bonferroni correction
result = multiarm_sample_size(
    n_arms=4,
    effect_size=0.5,
    alpha=0.05,
    power=0.80,
    comparison="bonferroni"
)
```

---

## Adaptive Designs

### Group Sequential Power

```python
from neopkpd.trial import group_sequential_power

# O'Brien-Fleming boundaries
result = group_sequential_power(
    n_per_arm=100,
    effect_size=0.4,
    sd=1.0,
    alpha=0.05,
    n_looks=3,
    information_fractions=[0.33, 0.67, 1.0],
    spending_function="obrien_fleming"
)

print(f"Overall power: {result.power:.1%}")
print(f"P(stop at look 1): {result.stop_probs[0]:.1%}")
print(f"P(stop at look 2): {result.stop_probs[1]:.1%}")
```

### Sample Size Re-estimation

```python
from neopkpd.trial import ssr_sample_size

# Initial sample size with interim re-estimation
result = ssr_sample_size(
    initial_n=50,
    target_power=0.80,
    effect_size=0.5,
    sd=1.0,
    interim_fraction=0.5,
    max_n=150
)

print(f"Initial N: {result.initial_n}")
print(f"Expected final N: {result.expected_n}")
print(f"Max N: {result.max_n}")
```

---

## Special Populations

### Pediatric Studies

```python
# Pediatric with higher variability
result = estimate_sample_size(
    target_power=0.80,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05,
    variability_inflation=1.2   # 20% higher variability
)
```

### Rare Disease

```python
# Small population constraints
result = trial.rare_disease_power(
    available_patients=50,
    effect_size=0.8,            # Expect large effect
    sd=1.0,
    alpha=0.10                  # Relaxed alpha
)

print(f"Achievable power: {result.power:.1%}")
print(f"Recommended design: {result.recommended_design}")
```

---

## Alpha Spending

### O'Brien-Fleming

```python
from neopkpd.trial import alpha_spending

# O'Brien-Fleming spending
alpha = alpha_spending(
    information_fraction=0.5,
    total_alpha=0.05,
    method="obrien_fleming"
)

print(f"Alpha spent at 50%: {alpha:.4f}")
```

### Pocock

```python
alpha = alpha_spending(
    information_fraction=0.5,
    total_alpha=0.05,
    method="pocock"
)
```

### Lan-DeMets

```python
# Flexible spending function
alpha = alpha_spending(
    information_fraction=0.5,
    total_alpha=0.05,
    method="lan_demets",
    rho=1.0                     # Parameter
)
```

---

## Complete Example

```python
from neopkpd import trial

# ==========================================
# Comprehensive Power Analysis
# ==========================================

print("=== Power Analysis for Phase III Trial ===\n")

# Study parameters
effect_size = 0.40              # Expected treatment effect
sd = 1.0                        # Population SD
alpha = 0.05                    # Two-sided alpha
target_power = 0.80             # Target power
dropout_rate = 0.15             # Expected dropout

# 1. Calculate base sample size
print("--- Sample Size Calculation ---")
result = trial.estimate_sample_size(
    target_power=target_power,
    effect_size=effect_size,
    sd=sd,
    alpha=alpha,
    dropout_rate=dropout_rate
)

print(f"Required completers: {result.n_per_arm} per arm")
print(f"Enrollment needed: {result.n_to_enroll} per arm")
print(f"Total enrollment: {result.n_to_enroll * 2}")
print(f"Achieved power: {result.achieved_power:.1%}")

# 2. Power curve
print("\n--- Power by Sample Size ---")
curve = trial.power_curve(
    effect_size=effect_size,
    sd=sd,
    alpha=alpha,
    n_range=range(30, 101, 10)
)

print("N/arm    Power")
print("-" * 20)
for n, power in zip(curve.n, curve.power):
    marker = "*" if power >= 0.80 else ""
    print(f"{n:5d}    {power:.1%}{marker}")

# 3. Effect size sensitivity
print("\n--- Detectable Effect Sizes ---")
print(f"With N = {result.n_per_arm} per arm:")

for power_target in [0.70, 0.80, 0.90]:
    mde = trial.minimum_detectable_effect(
        n_per_arm=result.n_per_arm,
        alpha=alpha,
        power=power_target,
        sd=sd
    )
    print(f"  {int(power_target*100)}% power: d = {mde:.3f}")

# 4. Simulation verification
print("\n--- Simulation Verification ---")
sim_result = trial.simulate_power(
    n_per_arm=result.n_per_arm,
    effect_size=effect_size,
    sd=sd,
    alpha=alpha,
    n_simulations=1000,
    seed=42
)

print(f"Simulated power: {sim_result.power:.1%}")
print(f"95% CI: ({sim_result.ci_lower:.1%}, {sim_result.ci_upper:.1%})")

# 5. Group sequential design
print("\n--- Group Sequential Design ---")
gs_result = trial.group_sequential_power(
    n_per_arm=result.n_per_arm,
    effect_size=effect_size,
    sd=sd,
    alpha=alpha,
    n_looks=2,
    information_fractions=[0.5, 1.0],
    spending_function="obrien_fleming"
)

print(f"Overall power: {gs_result.power:.1%}")
print(f"Efficacy boundaries: {gs_result.boundaries}")
print(f"P(stop at interim): {gs_result.stop_probs[0]:.1%}")
print(f"Expected sample size: {gs_result.expected_n:.0f}")

# 6. Summary
print("\n" + "=" * 50)
print("RECOMMENDATION")
print("=" * 50)
print(f"\nEnroll {result.n_to_enroll * 2} subjects total")
print(f"({result.n_to_enroll} per arm)")
print(f"\nWith {int(dropout_rate*100)}% dropout:")
print(f"  Expected completers: {result.n_per_arm * 2}")
print(f"  Power: {result.achieved_power:.1%}")
print(f"\nTo detect effect size d = {effect_size}")
print(f"At alpha = {alpha} (two-sided)")
```

---

## Power Tables

### Parallel Design Reference

| Effect Size | CV | N per Arm (80%) | N per Arm (90%) |
|-------------|-----|-----------------|-----------------|
| 0.3 | 30% | 176 | 235 |
| 0.4 | 30% | 99 | 132 |
| 0.5 | 30% | 64 | 85 |
| 0.6 | 30% | 44 | 59 |
| 0.7 | 30% | 33 | 44 |

### BE Crossover Reference

| CV Within | GMR | N (80%) | N (90%) |
|-----------|-----|---------|---------|
| 15% | 0.95 | 10 | 14 |
| 20% | 0.95 | 16 | 22 |
| 25% | 0.95 | 24 | 32 |
| 30% | 0.95 | 36 | 48 |
| 25% | 1.00 | 18 | 24 |

---

## Function Reference

| Function | Description |
|----------|-------------|
| `estimate_power_analytical` | Analytical power calculation |
| `estimate_sample_size` | Sample size for target power |
| `be_sample_size` | BE study sample size |
| `be_power` | BE study power |
| `power_curve` | Power vs sample size curve |
| `simulate_power` | Monte Carlo power |
| `multiarm_sample_size` | Multi-arm trial sample size |
| `alpha_spending` | Alpha spending functions |
| `minimum_detectable_effect` | MDE calculation |

---

## See Also

- [Study Designs](designs.md) - Trial design types
- [Dosing Regimens](dosing.md) - Dosing schedules
- [Virtual Population](population.md) - Subject generation
- [Julia Power Analysis](../../julia/trial/power.md) - Julia interface
