# Population NCA

Multi-subject non-compartmental analysis with summary statistics.

---

## Overview

```python
from neopkpd.nca import run_population_nca, summarize_population_nca

pop_result = run_population_nca(data, dose=100.0)
summary = summarize_population_nca(pop_result)
```

---

## Data Format

### Required Structure

```python
import pandas as pd

# DataFrame with required columns
data = pd.DataFrame({
    'subject_id': [1, 1, 1, 2, 2, 2, ...],      # Subject identifier
    'time': [0.0, 1.0, 4.0, 0.0, 1.0, 4.0, ...],  # Time points
    'conc': [0.0, 5.2, 2.1, 0.0, 4.8, 1.9, ...]   # Concentrations
})

# Dose can be column or single value
data['dose'] = [100.0, 100.0, 100.0, 150.0, ...]  # Per-subject dose
# OR
dose = 100.0  # Same dose for all
```

### With Covariates

```python
data = pd.DataFrame({
    'subject_id': [...],
    'time': [...],
    'conc': [...],
    'dose': [...],
    'weight': [70.0, 85.0, ...],
    'sex': ['M', 'F', ...],
    'formulation': ['Test', 'Reference', ...]
})
```

---

## Running Population NCA

### Basic Usage

```python
from neopkpd.nca import run_population_nca

pop_result = run_population_nca(
    data,
    subject_col='subject_id',
    time_col='time',
    conc_col='conc',
    dose=100.0  # OR dose_col='dose'
)

# Access individual results
for subject_id, result in pop_result.individual_results.items():
    print(f"Subject {subject_id}: AUC={result.auc_0_inf:.2f}")
```

### With Configuration

```python
from neopkpd.nca import run_population_nca, NCAConfig

config = NCAConfig(
    method="lin_log_mixed",
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9,
    lloq=0.05,
    blq_handling="missing"
)

pop_result = run_population_nca(
    data,
    config=config,
    route="extravascular",
    dose=100.0
)
```

### Variable Dose per Subject

```python
# Dose in data column
pop_result = run_population_nca(
    data,
    dose_col='dose'  # Use dose column
)
```

---

## Summary Statistics

### Generate Summary

```python
from neopkpd.nca import summarize_population_nca

summary = summarize_population_nca(pop_result)

# Access parameter statistics
print("=== AUC0-inf Summary ===")
print(f"N:        {summary['auc_0_inf']['n']}")
print(f"Mean:     {summary['auc_0_inf']['mean']:.2f}")
print(f"SD:       {summary['auc_0_inf']['sd']:.2f}")
print(f"CV%:      {summary['auc_0_inf']['cv_pct']:.1f}")
print(f"Median:   {summary['auc_0_inf']['median']:.2f}")
print(f"Min:      {summary['auc_0_inf']['min']:.2f}")
print(f"Max:      {summary['auc_0_inf']['max']:.2f}")
print(f"GeoMean:  {summary['auc_0_inf']['geomean']:.2f}")
print(f"GeoCV%:   {summary['auc_0_inf']['geocv_pct']:.1f}")
```

### Summary Table

```python
# Get formatted table
table = pop_result.summary_table()
print(table)

# Output:
# Parameter     N    Mean    SD     CV%    Median   GeoMean  GeoCV%
# -------------------------------------------------------------------
# Cmax          24   5.23    1.12   21.4   5.15     5.12     22.1
# Tmax          24   1.25    0.35   28.0   1.00     -        -
# AUC0-t        24   45.2    8.7    19.3   44.1     44.5     19.8
# AUC0-inf      24   48.1    9.2    19.1   47.2     47.4     19.5
# t1/2          24   4.52    0.85   18.8   4.35     4.45     19.2
# CL/F          24   2.12    0.42   19.8   2.08     2.08     20.1
```

### Specific Parameters

```python
# Summary for specific parameters only
summary = summarize_population_nca(
    pop_result,
    parameters=['cmax', 'auc_0_inf', 't_half', 'cl_f']
)
```

---

## Stratified Analysis

### By Single Variable

```python
from neopkpd.nca import stratified_population_nca

# Stratify by formulation
stratified = stratified_population_nca(
    data,
    strat_col='formulation',
    dose=100.0
)

for stratum, result in stratified.items():
    summary = summarize_population_nca(result)
    print(f"\n{stratum}:")
    print(f"  AUC mean: {summary['auc_0_inf']['mean']:.2f}")
    print(f"  Cmax mean: {summary['cmax']['mean']:.2f}")
```

### By Multiple Variables

```python
# Stratify by formulation and fed state
stratified = stratified_population_nca(
    data,
    strat_cols=['formulation', 'fed_state'],
    dose_col='dose'
)

# Access specific stratum
test_fed = stratified[('Test', 'Fed')]
ref_fasted = stratified[('Reference', 'Fasted')]
```

### Geometric Mean Ratio Between Strata

```python
test_summary = summarize_population_nca(stratified['Test'])
ref_summary = summarize_population_nca(stratified['Reference'])

# Calculate GMR
gmr_auc = test_summary['auc_0_inf']['geomean'] / ref_summary['auc_0_inf']['geomean']
gmr_cmax = test_summary['cmax']['geomean'] / ref_summary['cmax']['geomean']

print(f"AUC GMR: {gmr_auc * 100:.2f}%")
print(f"Cmax GMR: {gmr_cmax * 100:.2f}%")
```

---

## Export Results

### To DataFrame

```python
# Convert to DataFrame
results_df = pop_result.to_dataframe()

# Columns: subject_id, cmax, tmax, auc_0_t, auc_0_inf, t_half, cl_f, ...
print(results_df.head())
```

### To CSV

```python
# Export to CSV
pop_result.to_csv('nca_results.csv')
```

### To JSON

```python
# Export to JSON
pop_result.to_json('nca_results.json')
```

---

## Quality Assessment

### Check Individual Results

```python
# Identify subjects with quality issues
for subject_id, result in pop_result.individual_results.items():
    issues = []

    # Check λz R²
    if result.lambda_z_result.r_squared < 0.9:
        issues.append(f"Low λz R²: {result.lambda_z_result.r_squared:.3f}")

    # Check extrapolation
    if result.auc_extra_pct > 20:
        issues.append(f"High extrapolation: {result.auc_extra_pct:.1f}%")

    # Check λz estimation
    if math.isnan(result.lambda_z_result.lambda_z):
        issues.append("λz not estimable")

    if issues:
        print(f"Subject {subject_id}: {', '.join(issues)}")
```

### Summary Quality Metrics

```python
import math

# Count quality issues
n_total = len(pop_result.individual_results)
n_good_lz = sum(
    1 for r in pop_result.individual_results.values()
    if r.lambda_z_result.r_squared >= 0.9
)
n_low_extrap = sum(
    1 for r in pop_result.individual_results.values()
    if r.auc_extra_pct <= 20
)

print(f"λz R² ≥ 0.9: {n_good_lz}/{n_total} ({100*n_good_lz/n_total:.1f}%)")
print(f"AUC extrap ≤ 20%: {n_low_extrap}/{n_total} ({100*n_low_extrap/n_total:.1f}%)")
```

### Failed Subjects

```python
# Check for subjects with failed NCA
if pop_result.failed_subjects:
    print("Subjects with NCA failures:")
    for subject_id, reason in pop_result.failed_subjects.items():
        print(f"  {subject_id}: {reason}")
```

---

## Example: Complete Population Analysis

```python
import pandas as pd
from neopkpd.nca import (
    run_population_nca, summarize_population_nca,
    stratified_population_nca, NCAConfig
)

# Load study data
data = pd.read_csv('pk_study_data.csv')

# Configure NCA
config = NCAConfig(
    method="lin_log_mixed",
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9,
    lloq=0.1,
    blq_handling="missing"
)

# Run population NCA
pop_result = run_population_nca(
    data,
    subject_col='SUBJID',
    time_col='TIME',
    conc_col='CONC',
    dose_col='DOSE',
    config=config,
    route='extravascular'
)

# Overall summary
print("=" * 60)
print("POPULATION NCA SUMMARY")
print("=" * 60)

summary = summarize_population_nca(pop_result)
print(f"\nTotal subjects: {len(pop_result.individual_results)}")
print(f"Failed subjects: {len(pop_result.failed_subjects)}")

print("\n{:<12} {:>6} {:>10} {:>10} {:>8} {:>10} {:>8}".format(
    "Parameter", "N", "Mean", "SD", "CV%", "GeoMean", "GeoCV%"
))
print("-" * 60)

for param in ['cmax', 'tmax', 'auc_0_t', 'auc_0_inf', 't_half', 'cl_f']:
    s = summary[param]
    if param == 'tmax':
        print(f"{param:<12} {s['n']:>6} {s['mean']:>10.2f} {s['sd']:>10.2f} {s['cv_pct']:>8.1f} {'N/A':>10} {'N/A':>8}")
    else:
        print(f"{param:<12} {s['n']:>6} {s['mean']:>10.2f} {s['sd']:>10.2f} {s['cv_pct']:>8.1f} {s['geomean']:>10.2f} {s['geocv_pct']:>8.1f}")

# Stratified by formulation
if 'FORMULATION' in data.columns:
    print("\n" + "=" * 60)
    print("STRATIFIED BY FORMULATION")
    print("=" * 60)

    stratified = stratified_population_nca(
        data,
        strat_col='FORMULATION',
        config=config,
        dose_col='DOSE'
    )

    for form in stratified:
        s = summarize_population_nca(stratified[form])
        print(f"\n{form}:")
        print(f"  N:        {s['auc_0_inf']['n']}")
        print(f"  AUC0-inf: {s['auc_0_inf']['geomean']:.2f} (GeoCV: {s['auc_0_inf']['geocv_pct']:.1f}%)")
        print(f"  Cmax:     {s['cmax']['geomean']:.2f} (GeoCV: {s['cmax']['geocv_pct']:.1f}%)")

# Export results
results_df = pop_result.to_dataframe()
results_df.to_csv('nca_individual_results.csv', index=False)
print(f"\nResults exported to nca_individual_results.csv")
```

---

## Multiple Dose Population NCA

```python
# Steady-state population analysis
pop_result = run_population_nca(
    data,
    subject_col='SUBJID',
    time_col='TIME',
    conc_col='CONC',
    dose=100.0,
    dosing_type='steady_state',
    tau=12.0,
    config=config
)

# Access steady-state metrics
summary = summarize_population_nca(pop_result)
print(f"Cmax,ss:    {summary['cmax']['geomean']:.2f}")
print(f"Cmin,ss:    {summary['cmin']['geomean']:.2f}")
print(f"Cavg,ss:    {summary['cavg']['geomean']:.2f}")
print(f"AUC0-tau:   {summary['auc_0_tau']['geomean']:.2f}")
print(f"PTF:        {summary['fluctuation']['mean']:.1f}%")
```

---

## See Also

- [run_nca Function](run-nca.md) - Individual NCA
- [NCA Configuration](config.md) - Configuration options
- [Bioequivalence](bioequivalence.md) - BE analysis from population NCA

