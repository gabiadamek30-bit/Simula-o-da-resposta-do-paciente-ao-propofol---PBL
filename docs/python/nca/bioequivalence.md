# Bioequivalence Analysis

Bioequivalence (BE) assessment following FDA, EMA, and Health Canada guidance.

---

## Overview

```python
from neopkpd.nca import bioequivalence_90ci, tost_analysis, be_conclusion

# Calculate 90% CI for GMR
ci = bioequivalence_90ci(test_values, reference_values)
print(f"90% CI: ({ci.lower:.3f}, {ci.upper:.3f})")
print(f"BE: {ci.be_met}")
```

---

## Quick Start

```python
from neopkpd.nca import bioequivalence_90ci

# AUC values from crossover study
test_auc = [45.2, 52.1, 38.7, 61.3, 49.8, 55.4, 42.1, 58.9]
ref_auc = [48.1, 55.3, 41.2, 58.9, 52.4, 53.2, 44.8, 60.1]

# Calculate 90% confidence interval
ci = bioequivalence_90ci(test_auc, ref_auc)

print(f"GMR: {ci.gmr * 100:.2f}%")
print(f"90% CI: [{ci.lower * 100:.2f}%, {ci.upper * 100:.2f}%]")
print(f"Within-subject CV: {ci.cv_within * 100:.1f}%")
print(f"BE criteria met: {ci.be_met}")
```

---

## 90% Confidence Interval

### Basic Calculation

```python
from neopkpd.nca import bioequivalence_90ci

ci = bioequivalence_90ci(test_values, ref_values)

# Results
print(f"Geometric Mean Ratio: {ci.gmr * 100:.2f}%")
print(f"90% CI Lower: {ci.lower * 100:.2f}%")
print(f"90% CI Upper: {ci.upper * 100:.2f}%")
print(f"BE Met (80-125%): {ci.be_met}")
```

### Result Attributes

```python
class BioequivalenceResult:
    gmr: float           # Geometric mean ratio (Test/Reference)
    lower: float         # 90% CI lower bound
    upper: float         # 90% CI upper bound
    be_met: bool         # True if within 80-125%
    cv_within: float     # Within-subject CV
    n_subjects: int      # Number of subjects
    se: float            # Standard error of log GMR
```

### Custom Acceptance Limits

```python
# Custom limits for narrow therapeutic index drugs
ci = bioequivalence_90ci(
    test_values, ref_values,
    lower_limit=0.90,  # 90%
    upper_limit=1.11   # 111%
)
```

---

## TOST Analysis

Two One-Sided Tests for equivalence:

```python
from neopkpd.nca import tost_analysis

result = tost_analysis(
    test_values,
    ref_values,
    theta_lower=0.80,
    theta_upper=1.25,
    alpha=0.05
)

print(f"Lower test p-value: {result.p_lower:.4f}")
print(f"Upper test p-value: {result.p_upper:.4f}")
print(f"BE concluded: {result.be_concluded}")
```

### TOST Result

```python
class TOSTResult:
    t_lower: float       # t-statistic for lower bound
    t_upper: float       # t-statistic for upper bound
    p_lower: float       # p-value for lower test
    p_upper: float       # p-value for upper test
    be_concluded: bool   # True if both p < alpha
```

---

## Geometric Mean Ratio

```python
from neopkpd.nca import geometric_mean_ratio

gmr = geometric_mean_ratio(test_values, ref_values)
print(f"GMR: {gmr * 100:.2f}%")

# Point estimate only (no CI)
```

---

## Within-Subject CV

```python
from neopkpd.nca import within_subject_cv

# From crossover data
cv = within_subject_cv(
    data,
    value_col='auc',
    subject_col='subject_id',
    period_col='period'
)

print(f"Within-subject CV: {cv * 100:.1f}%")

# Classification
if cv < 0.30:
    print("Standard variability drug")
elif cv < 0.40:
    print("Moderately variable drug")
else:
    print("Highly variable drug (HVD)")
```

---

## Study Designs

### 2x2 Crossover (Standard)

```python
from neopkpd.nca import analyze_crossover_be

result = analyze_crossover_be(
    data,
    treatment_col='formulation',
    subject_col='subject_id',
    period_col='period',
    sequence_col='sequence',
    response_col='auc',
    design='2x2'
)

print(f"GMR: {result.gmr * 100:.2f}%")
print(f"90% CI: [{result.lower * 100:.2f}%, {result.upper * 100:.2f}%]")
print(f"Period effect p-value: {result.period_pvalue:.4f}")
print(f"Sequence effect p-value: {result.sequence_pvalue:.4f}")
```

### Replicate Designs

```python
# 2x4 Replicate (for HVD)
result = analyze_crossover_be(
    data,
    design='2x4_replicate',
    response_col='auc'
)

# Partial replicate 3x3
result = analyze_crossover_be(
    data,
    design='3x3_partial',
    response_col='auc'
)
```

### Parallel Group

```python
result = analyze_parallel_be(
    data,
    treatment_col='formulation',
    subject_col='subject_id',
    response_col='auc'
)
```

---

## Reference-Scaled BE (HVD)

### RSABE (FDA)

For highly variable drugs (CV > 30%):

```python
from neopkpd.nca import rsabe_analysis

result = rsabe_analysis(
    test_values,
    ref_values,
    design='2x4_replicate'
)

print(f"Reference CV: {result.cv_reference * 100:.1f}%")
print(f"Scaling applied: {result.scaling_applied}")
print(f"Scaled criterion met: {result.scaled_criterion_met}")
print(f"GMR constraint (80-125%): {result.gmr_constraint_met}")
print(f"RSABE conclusion: {result.be_met}")
```

### ABEL (EMA)

Average Bioequivalence with Expanding Limits:

```python
from neopkpd.nca import abel_analysis

result = abel_analysis(
    test_values,
    ref_values,
    design='2x4_replicate'
)

print(f"Reference CV: {result.cv_reference * 100:.1f}%")
print(f"Widened limits: [{result.lower_limit * 100:.2f}%, {result.upper_limit * 100:.2f}%]")
print(f"GMR within limits: {result.within_widened}")
print(f"GMR constraint (80-125%): {result.gmr_constraint_met}")
print(f"ABEL conclusion: {result.be_met}")
```

---

## Sample Size Calculation

### For BE Study

```python
from neopkpd.nca import be_sample_size

n = be_sample_size(
    cv=0.25,           # Expected within-subject CV (25%)
    gmr=0.95,          # Expected GMR
    power=0.80,        # Target power
    alpha=0.05,        # Significance level
    design='2x2'       # Study design
)

print(f"Required subjects per sequence: {n}")
print(f"Total subjects: {2 * n}")
```

### Power Calculation

```python
from neopkpd.nca import be_power

power = be_power(
    n=24,
    cv=0.25,
    gmr=0.95,
    design='2x2'
)

print(f"Expected power: {power * 100:.1f}%")
```

### Sample Size Table

| CV% | GMR=0.95 | GMR=0.90 | GMR=1.00 |
|-----|----------|----------|----------|
| 15 | 10 | 14 | 8 |
| 20 | 16 | 24 | 12 |
| 25 | 24 | 36 | 18 |
| 30 | 34 | 50 | 26 |
| 35 | 46 | 68 | 36 |

---

## Complete BE Analysis

```python
import pandas as pd
from neopkpd.nca import (
    run_population_nca, bioequivalence_90ci,
    tost_analysis, within_subject_cv, NCAConfig
)

# Load crossover study data
data = pd.read_csv('be_study_data.csv')

# Run NCA for all subjects
config = NCAConfig(
    method="lin_log_mixed",
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9
)

pop_result = run_population_nca(
    data,
    subject_col='SUBJID',
    time_col='TIME',
    conc_col='CONC',
    dose_col='DOSE',
    config=config
)

# Get NCA results by formulation
results_df = pop_result.to_dataframe()
results_df = results_df.merge(
    data[['SUBJID', 'FORMULATION']].drop_duplicates(),
    on='SUBJID'
)

test_auc = results_df[results_df['FORMULATION'] == 'Test']['auc_0_inf'].values
ref_auc = results_df[results_df['FORMULATION'] == 'Reference']['auc_0_inf'].values
test_cmax = results_df[results_df['FORMULATION'] == 'Test']['cmax'].values
ref_cmax = results_df[results_df['FORMULATION'] == 'Reference']['cmax'].values

# BE Analysis for AUC
print("=" * 60)
print("BIOEQUIVALENCE ANALYSIS")
print("=" * 60)

print("\n--- AUC0-inf ---")
auc_be = bioequivalence_90ci(test_auc, ref_auc)
print(f"N: {auc_be.n_subjects}")
print(f"GMR: {auc_be.gmr * 100:.2f}%")
print(f"90% CI: [{auc_be.lower * 100:.2f}%, {auc_be.upper * 100:.2f}%]")
print(f"Within-subject CV: {auc_be.cv_within * 100:.1f}%")
print(f"BE criteria met: {'YES' if auc_be.be_met else 'NO'}")

# BE Analysis for Cmax
print("\n--- Cmax ---")
cmax_be = bioequivalence_90ci(test_cmax, ref_cmax)
print(f"N: {cmax_be.n_subjects}")
print(f"GMR: {cmax_be.gmr * 100:.2f}%")
print(f"90% CI: [{cmax_be.lower * 100:.2f}%, {cmax_be.upper * 100:.2f}%]")
print(f"Within-subject CV: {cmax_be.cv_within * 100:.1f}%")
print(f"BE criteria met: {'YES' if cmax_be.be_met else 'NO'}")

# Overall Conclusion
print("\n" + "=" * 60)
overall_be = auc_be.be_met and cmax_be.be_met
print(f"OVERALL BIOEQUIVALENCE: {'ESTABLISHED' if overall_be else 'NOT ESTABLISHED'}")
print("=" * 60)

# TOST confirmation
print("\n--- TOST Analysis ---")
auc_tost = tost_analysis(test_auc, ref_auc)
print(f"AUC: Lower p={auc_tost.p_lower:.4f}, Upper p={auc_tost.p_upper:.4f}")
print(f"     BE concluded: {auc_tost.be_concluded}")

cmax_tost = tost_analysis(test_cmax, ref_cmax)
print(f"Cmax: Lower p={cmax_tost.p_lower:.4f}, Upper p={cmax_tost.p_upper:.4f}")
print(f"      BE concluded: {cmax_tost.be_concluded}")
```

---

## Regulatory Reports

### FDA-Style Report

```python
from neopkpd.nca import generate_be_report

report = generate_be_report(
    test_auc, ref_auc, test_cmax, ref_cmax,
    regulatory='FDA',
    study_design='2x2_crossover'
)

print(report)
```

### EMA-Style Report

```python
report = generate_be_report(
    test_auc, ref_auc, test_cmax, ref_cmax,
    regulatory='EMA',
    study_design='2x2_crossover',
    include_tmax=True  # EMA includes Tmax
)

print(report)
```

---

## Highly Variable Drug Analysis

```python
from neopkpd.nca import (
    within_subject_cv, rsabe_analysis, abel_analysis
)

# Check if HVD criteria met
cv = within_subject_cv(data, 'auc', 'subject_id', 'period')
print(f"Within-subject CV: {cv * 100:.1f}%")

if cv > 0.30:
    print("HVD criteria met - applying reference scaling")

    # FDA RSABE
    rsabe = rsabe_analysis(
        test_values, ref_values,
        design='2x4_replicate'
    )
    print(f"\nRSABE (FDA): {rsabe.be_met}")

    # EMA ABEL
    abel = abel_analysis(
        test_values, ref_values,
        design='2x4_replicate'
    )
    print(f"ABEL (EMA): {abel.be_met}")
else:
    print("Standard ABE analysis applicable")
    be = bioequivalence_90ci(test_values, ref_values)
    print(f"ABE: {be.be_met}")
```

---

## Acceptance Criteria

### Standard BE (80-125%)

| Parameter | Acceptance Range |
|-----------|-----------------|
| AUC | 80.00% - 125.00% |
| Cmax | 80.00% - 125.00% |

### Narrow Therapeutic Index (90-111%)

```python
# NTI drugs (warfarin, phenytoin, etc.)
be = bioequivalence_90ci(
    test_values, ref_values,
    lower_limit=0.90,
    upper_limit=1.11
)
```

### Highly Variable Drugs (Widened)

| Reference CV | EMA Widened Range |
|--------------|-------------------|
| 30-50% | Calculated per CV |
| >50% | 69.84% - 143.19% |

---

## See Also

- [run_nca Function](run-nca.md) - NCA calculation
- [NCA Configuration](config.md) - Configuration options
- [Population NCA](population.md) - Multi-subject analysis

