# Covariates

Comprehensive guide for incorporating patient characteristics into population models in NeoPKPD Python.

---

## Overview

Covariate models explain inter-individual variability by relating PK parameters to measurable patient characteristics.

```python
from neopkpd import simulate_population_oral

# Weight effect on CL and V
covariate_effects = [
    {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0}
]

result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omegas={"Ka": 0.16, "CL": 0.04, "V": 0.02},
    covariates=[{"WT": w} for w in weights],
    covariate_effects=covariate_effects,
    seed=42
)
```

---

## Covariate Effect Types

### Power Model

Most common for body size effects:

$$\theta_i = \theta_{pop} \cdot \left(\frac{cov_i}{ref}\right)^{\beta}$$

```python
# Allometric scaling
effect = {
    "param": "CL",
    "cov": "WT",
    "ref": 70.0,           # Reference weight
    "kind": "PowerCovariate",
    "beta": 0.75           # Allometric exponent
}

# Common allometric exponents:
# Clearance: 0.75 (3/4 power law)
# Volume: 1.0 (proportional to body size)
# Half-life: 0.25 (1/4 power)
```

### Linear Model

Effect proportional to deviation from reference:

$$\theta_i = \theta_{pop} \cdot (1 + \beta \cdot (cov_i - ref))$$

```python
# Age effect
effect = {
    "param": "CL",
    "cov": "AGE",
    "ref": 45.0,           # Reference age
    "kind": "LinearCovariate",
    "beta": -0.006         # -0.6% per year
}

# For a 65-year-old:
# CL = CL_pop * (1 + (-0.006) * (65 - 45))
# CL = CL_pop * 0.88 (12% lower)
```

### Exponential Model

Multiplicative effect:

$$\theta_i = \theta_{pop} \cdot e^{\beta \cdot (cov_i - ref)}$$

```python
# Renal function effect
effect = {
    "param": "CL",
    "cov": "CRCL",
    "ref": 100.0,          # Normal CRCL
    "kind": "ExpCovariate",
    "beta": 0.005
}

# For CRCL = 50:
# CL = CL_pop * exp(0.005 * (50 - 100))
# CL = CL_pop * 0.78 (22% lower)
```

---

## Effect Specification Format

### Dictionary Format

```python
covariate_effects = [
    {
        "param": str,      # Target parameter: "CL", "V", "Ka", etc.
        "cov": str,        # Covariate name: "WT", "AGE", "CRCL", etc.
        "ref": float,      # Reference value for centering
        "kind": str,       # "PowerCovariate", "LinearCovariate", "ExpCovariate"
        "beta": float      # Effect coefficient
    }
]
```

### Multiple Effects

```python
# Comprehensive covariate model
covariate_effects = [
    # Body size effects
    {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0},

    # Age effect
    {"param": "CL", "cov": "AGE", "ref": 45.0, "kind": "LinearCovariate", "beta": -0.005},

    # Renal function
    {"param": "CL", "cov": "CRCL", "ref": 100.0, "kind": "LinearCovariate", "beta": 0.004},

    # Sex effect (0=male, 1=female)
    {"param": "CL", "cov": "SEX", "ref": 0.0, "kind": "LinearCovariate", "beta": -0.15}
]
```

---

## Subject Covariates

### Providing Covariate Data

```python
import numpy as np

# Generate realistic population
np.random.seed(42)
n = 200

covariates = []
for _ in range(n):
    sex = 1 if np.random.random() < 0.5 else 0
    age = np.clip(np.random.normal(50, 12), 18, 85)
    # Weight depends on sex
    wt = np.random.normal(75 if sex == 0 else 65, 12 if sex == 0 else 10)
    wt = np.clip(wt, 45, 150)
    # CRCL depends on age
    crcl = np.clip(130 - age * 0.8 + np.random.normal(0, 15), 15, 150)

    covariates.append({
        "WT": wt,
        "AGE": age,
        "SEX": sex,
        "CRCL": crcl
    })
```

### Covariate Summary

```python
def summarize_covariates(covariates):
    """Summarize population covariates."""
    wts = [c["WT"] for c in covariates]
    ages = [c["AGE"] for c in covariates]
    crcls = [c["CRCL"] for c in covariates]
    sexes = [c["SEX"] for c in covariates]

    print("Covariate Summary:")
    print(f"  WT: {np.mean(wts):.1f} ± {np.std(wts):.1f} kg")
    print(f"  AGE: {np.mean(ages):.1f} ± {np.std(ages):.1f} years")
    print(f"  CRCL: {np.mean(crcls):.1f} ± {np.std(crcls):.1f} mL/min")
    print(f"  Female: {np.mean(sexes)*100:.0f}%")

summarize_covariates(covariates)
```

---

## Common Covariate Models

### Allometric Scaling

```python
# Standard allometry for all structural parameters
allometric_effects = [
    {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V1", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0},
    {"param": "Q", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V2", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0}
]
```

### Renal Function

```python
# Creatinine clearance effect
renal_effects = [
    {"param": "CL", "cov": "CRCL", "ref": 100.0, "kind": "LinearCovariate", "beta": 0.005}
]

# Or using eGFR
egfr_effects = [
    {"param": "CL", "cov": "EGFR", "ref": 90.0, "kind": "PowerCovariate", "beta": 0.5}
]
```

### Age Effects

```python
# Linear decline with age
age_linear = {"param": "CL", "cov": "AGE", "ref": 40.0, "kind": "LinearCovariate", "beta": -0.006}

# Maturation (pediatric)
maturation = {"param": "CL", "cov": "PMA", "ref": 40.0, "kind": "PowerCovariate", "beta": 0.75}
```

### Sex/Gender

```python
# Female vs male (reference)
sex_effect = {"param": "CL", "cov": "SEX", "ref": 0.0, "kind": "LinearCovariate", "beta": -0.12}
# Female CL = CL_pop * (1 - 0.12) = 88% of male CL
```

### Genetic Polymorphisms

```python
# CYP2D6 status as numeric score
# PM=0, IM=0.5, EM=1.0 (ref), UM=1.5
cyp_effect = {"param": "CL", "cov": "CYP2D6", "ref": 1.0, "kind": "PowerCovariate", "beta": 1.0}
```

---

## Simulation with Covariates

### Complete Example

```python
from neopkpd import simulate_population_oral
import numpy as np

# 1. Generate population
np.random.seed(42)
n = 200

covariates = []
for _ in range(n):
    sex = 1 if np.random.random() < 0.48 else 0
    age = np.clip(np.random.normal(55, 14), 18, 85)
    wt = np.clip(np.random.normal(70 + sex * -8, 15), 40, 150)
    crcl = np.clip(130 - age * 0.9 + np.random.normal(0, 20), 15, 150)

    covariates.append({"WT": wt, "AGE": age, "SEX": sex, "CRCL": crcl})

# 2. Define covariate model
covariate_effects = [
    {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0},
    {"param": "CL", "cov": "AGE", "ref": 45.0, "kind": "LinearCovariate", "beta": -0.005},
    {"param": "CL", "cov": "CRCL", "ref": 100.0, "kind": "LinearCovariate", "beta": 0.003},
    {"param": "CL", "cov": "SEX", "ref": 0.0, "kind": "LinearCovariate", "beta": -0.10}
]

# 3. Reduced IIV (covariates explain variability)
omegas = {"Ka": 0.16, "CL": 0.0225, "V": 0.01}

# 4. Simulate
result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=n,
    omegas=omegas,
    covariates=covariates,
    covariate_effects=covariate_effects,
    seed=12345
)

# 5. Analyze covariate impact
cl_vals = np.array([p["CL"] for p in result.individual_params])
wts = np.array([c["WT"] for c in covariates])
ages = np.array([c["AGE"] for c in covariates])
crcls = np.array([c["CRCL"] for c in covariates])
sexes = np.array([c["SEX"] for c in covariates])

print("Covariate-CL Correlations:")
print(f"  WT-CL:   r = {np.corrcoef(wts, cl_vals)[0,1]:.3f}")
print(f"  AGE-CL:  r = {np.corrcoef(ages, cl_vals)[0,1]:.3f}")
print(f"  CRCL-CL: r = {np.corrcoef(crcls, cl_vals)[0,1]:.3f}")
```

---

## Covariate Impact Analysis

### Group Comparisons

```python
# CL by sex
male_cl = cl_vals[sexes == 0]
female_cl = cl_vals[sexes == 1]
print(f"Male CL: {np.mean(male_cl):.2f} L/hr")
print(f"Female CL: {np.mean(female_cl):.2f} L/hr")
print(f"Ratio (F/M): {np.mean(female_cl)/np.mean(male_cl):.2f}")

# AUC by renal function
auc_vals = [np.trapz(ind.concentrations, result.times) for ind in result.individuals]

normal_crcl = [auc_vals[i] for i in range(n) if crcls[i] >= 90]
mild_ri = [auc_vals[i] for i in range(n) if 60 <= crcls[i] < 90]
moderate_ri = [auc_vals[i] for i in range(n) if 30 <= crcls[i] < 60]

print(f"\nAUC by Renal Function:")
print(f"  Normal (CRCL ≥90): {np.mean(normal_crcl):.1f} mg·hr/L")
print(f"  Mild RI (60-89):   {np.mean(mild_ri):.1f} mg·hr/L")
print(f"  Moderate RI (30-59): {np.mean(moderate_ri):.1f} mg·hr/L")
```

### Quartile Analysis

```python
# CL by weight quartiles
wt_quartiles = np.percentile(wts, [25, 50, 75])

q1_cl = cl_vals[wts <= wt_quartiles[0]]
q2_cl = cl_vals[(wts > wt_quartiles[0]) & (wts <= wt_quartiles[1])]
q3_cl = cl_vals[(wts > wt_quartiles[1]) & (wts <= wt_quartiles[2])]
q4_cl = cl_vals[wts > wt_quartiles[2]]

print("CL by Weight Quartile:")
print(f"  Q1 (≤{wt_quartiles[0]:.0f} kg): {np.mean(q1_cl):.2f} L/hr")
print(f"  Q2 ({wt_quartiles[0]:.0f}-{wt_quartiles[1]:.0f} kg): {np.mean(q2_cl):.2f} L/hr")
print(f"  Q3 ({wt_quartiles[1]:.0f}-{wt_quartiles[2]:.0f} kg): {np.mean(q3_cl):.2f} L/hr")
print(f"  Q4 (>{wt_quartiles[2]:.0f} kg): {np.mean(q4_cl):.2f} L/hr")
```

---

## Exposure Predictions

### Subgroup Exposure

```python
def predict_exposure_subgroup(result, covariates, subgroup_filter):
    """Calculate exposure metrics for a subgroup."""
    indices = [i for i, c in enumerate(covariates) if subgroup_filter(c)]

    if not indices:
        return None

    cmax_vals = [max(result.individuals[i].concentrations) for i in indices]
    auc_vals = [np.trapz(result.individuals[i].concentrations, result.times)
                for i in indices]

    return {
        "n": len(indices),
        "Cmax_mean": np.mean(cmax_vals),
        "Cmax_cv": np.std(cmax_vals) / np.mean(cmax_vals) * 100,
        "AUC_mean": np.mean(auc_vals),
        "AUC_cv": np.std(auc_vals) / np.mean(auc_vals) * 100
    }

# Compare subgroups
elderly = predict_exposure_subgroup(result, covariates, lambda c: c["AGE"] >= 65)
young = predict_exposure_subgroup(result, covariates, lambda c: c["AGE"] < 40)

print("Exposure by Age Group:")
print(f"  Young (<40y): AUC = {young['AUC_mean']:.1f} mg·hr/L")
print(f"  Elderly (≥65y): AUC = {elderly['AUC_mean']:.1f} mg·hr/L")
print(f"  Ratio (elderly/young): {elderly['AUC_mean']/young['AUC_mean']:.2f}")
```

---

## Visualization

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Weight vs CL
ax = axes[0, 0]
ax.scatter(wts, cl_vals, alpha=0.5)
# Add regression line
m, b = np.polyfit(wts, cl_vals, 1)
ax.plot([min(wts), max(wts)], [m*min(wts)+b, m*max(wts)+b], 'r--')
ax.set_xlabel('Weight (kg)')
ax.set_ylabel('CL (L/hr)')
ax.set_title(f'Weight vs CL (r={np.corrcoef(wts, cl_vals)[0,1]:.3f})')

# Age vs CL
ax = axes[0, 1]
ax.scatter(ages, cl_vals, alpha=0.5)
m, b = np.polyfit(ages, cl_vals, 1)
ax.plot([min(ages), max(ages)], [m*min(ages)+b, m*max(ages)+b], 'r--')
ax.set_xlabel('Age (years)')
ax.set_ylabel('CL (L/hr)')
ax.set_title(f'Age vs CL (r={np.corrcoef(ages, cl_vals)[0,1]:.3f})')

# CRCL vs CL
ax = axes[1, 0]
ax.scatter(crcls, cl_vals, alpha=0.5)
m, b = np.polyfit(crcls, cl_vals, 1)
ax.plot([min(crcls), max(crcls)], [m*min(crcls)+b, m*max(crcls)+b], 'r--')
ax.set_xlabel('CRCL (mL/min)')
ax.set_ylabel('CL (L/hr)')
ax.set_title(f'CRCL vs CL (r={np.corrcoef(crcls, cl_vals)[0,1]:.3f})')

# CL by sex (box plot)
ax = axes[1, 1]
ax.boxplot([male_cl, female_cl], labels=['Male', 'Female'])
ax.set_ylabel('CL (L/hr)')
ax.set_title('CL by Sex')

plt.tight_layout()
plt.savefig('covariate_effects.png', dpi=150)
plt.show()
```

---

## Complete Example

```python
from neopkpd import simulate_population_oral
import numpy as np
import matplotlib.pyplot as plt

# =========================================
# Comprehensive Covariate Analysis
# =========================================

print("=== Covariate Modeling ===\n")

# 1. Generate realistic patient population
np.random.seed(42)
n = 300

covariates = []
for _ in range(n):
    sex = 1 if np.random.random() < 0.48 else 0
    age = np.clip(np.random.normal(58, 15), 18, 90)
    bmi = np.clip(np.random.normal(27, 5), 18, 45)
    ht = np.random.normal(175 - sex * 12, 8)
    wt = bmi * (ht / 100) ** 2
    wt = np.clip(wt, 40, 180)
    crcl = max(15, 140 - age + (1 - sex) * 15 + np.random.normal(0, 20))

    covariates.append({
        "WT": wt, "AGE": age, "SEX": sex,
        "HT": ht, "BMI": bmi, "CRCL": crcl
    })

# 2. Covariate summary
print("--- Population Covariates ---")
for cov in ["WT", "AGE", "BMI", "CRCL"]:
    vals = [c[cov] for c in covariates]
    print(f"  {cov}: {np.mean(vals):.1f} ± {np.std(vals):.1f}")
sexes = [c["SEX"] for c in covariates]
print(f"  Female: {np.mean(sexes)*100:.0f}%")

# 3. Covariate model
covariate_effects = [
    # Allometry
    {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0},
    # Demographics
    {"param": "CL", "cov": "AGE", "ref": 50.0, "kind": "LinearCovariate", "beta": -0.006},
    {"param": "CL", "cov": "SEX", "ref": 0.0, "kind": "LinearCovariate", "beta": -0.08},
    # Renal
    {"param": "CL", "cov": "CRCL", "ref": 100.0, "kind": "LinearCovariate", "beta": 0.004}
]

print("\n--- Covariate Effects ---")
for eff in covariate_effects:
    print(f"  {eff['param']} ~ {eff['cov']}: {eff['kind']} (β={eff['beta']})")

# 4. Simulate
print("\n--- Simulation ---")
result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=48.0, saveat=0.5,
    n=n,
    omegas={"Ka": 0.16, "CL": 0.02, "V": 0.01},  # Low residual IIV
    covariates=covariates,
    covariate_effects=covariate_effects,
    seed=12345
)
print(f"Simulated {result.n_subjects} subjects")

# 5. Parameter distributions
print("\n--- Realized Parameters ---")
cl_vals = np.array([p["CL"] for p in result.individual_params])
v_vals = np.array([p["V"] for p in result.individual_params])
ka_vals = np.array([p["Ka"] for p in result.individual_params])

print(f"CL: {np.mean(cl_vals):.2f} ± {np.std(cl_vals):.2f} L/hr (CV={np.std(cl_vals)/np.mean(cl_vals)*100:.0f}%)")
print(f"V:  {np.mean(v_vals):.1f} ± {np.std(v_vals):.1f} L (CV={np.std(v_vals)/np.mean(v_vals)*100:.0f}%)")

# 6. Covariate correlations
print("\n--- Covariate Correlations with CL ---")
wts = np.array([c["WT"] for c in covariates])
ages = np.array([c["AGE"] for c in covariates])
crcls = np.array([c["CRCL"] for c in covariates])

for cov, vals in [("WT", wts), ("AGE", ages), ("CRCL", crcls)]:
    r = np.corrcoef(vals, cl_vals)[0, 1]
    print(f"  {cov}: r = {r:.3f}")

# 7. Subgroup analysis
print("\n--- Subgroup Exposure Analysis ---")
auc_vals = [np.trapz(ind.concentrations, result.times) for ind in result.individuals]

# By renal function
normal = [auc_vals[i] for i in range(n) if covariates[i]["CRCL"] >= 90]
mild = [auc_vals[i] for i in range(n) if 60 <= covariates[i]["CRCL"] < 90]
moderate = [auc_vals[i] for i in range(n) if 30 <= covariates[i]["CRCL"] < 60]
severe = [auc_vals[i] for i in range(n) if covariates[i]["CRCL"] < 30]

print("AUC by Renal Function:")
print(f"  Normal (n={len(normal)}): {np.mean(normal):.1f} mg·hr/L")
if mild: print(f"  Mild RI (n={len(mild)}): {np.mean(mild):.1f} mg·hr/L ({np.mean(mild)/np.mean(normal):.2f}x)")
if moderate: print(f"  Moderate RI (n={len(moderate)}): {np.mean(moderate):.1f} mg·hr/L ({np.mean(moderate)/np.mean(normal):.2f}x)")
if severe: print(f"  Severe RI (n={len(severe)}): {np.mean(severe):.1f} mg·hr/L ({np.mean(severe)/np.mean(normal):.2f}x)")

# By age
young = [auc_vals[i] for i in range(n) if covariates[i]["AGE"] < 40]
middle = [auc_vals[i] for i in range(n) if 40 <= covariates[i]["AGE"] < 65]
elderly = [auc_vals[i] for i in range(n) if covariates[i]["AGE"] >= 65]

print("\nAUC by Age:")
print(f"  <40 years (n={len(young)}): {np.mean(young):.1f} mg·hr/L")
print(f"  40-64 years (n={len(middle)}): {np.mean(middle):.1f} mg·hr/L")
print(f"  ≥65 years (n={len(elderly)}): {np.mean(elderly):.1f} mg·hr/L ({np.mean(elderly)/np.mean(young):.2f}x)")

# 8. Dosing implications
print("\n--- Dosing Implications ---")
# Target AUC of 50 mg·hr/L
target_auc = 50.0
dose_adjustments = []

for i in range(n):
    current_auc = auc_vals[i]
    adj_factor = target_auc / current_auc
    adjusted_dose = 500.0 * adj_factor
    dose_adjustments.append(adjusted_dose)

print(f"Target AUC: {target_auc} mg·hr/L")
print(f"Dose adjustments needed:")
print(f"  Range: {np.min(dose_adjustments):.0f} - {np.max(dose_adjustments):.0f} mg")
print(f"  Mean: {np.mean(dose_adjustments):.0f} mg")
print(f"  Subjects needing <75% dose: {sum(1 for d in dose_adjustments if d < 375)/n*100:.0f}%")
print(f"  Subjects needing >125% dose: {sum(1 for d in dose_adjustments if d > 625)/n*100:.0f}%")
```

---

## See Also

- [Population IV Bolus](iv-bolus.md) - IV bolus with covariates
- [Population Oral](oral.md) - Oral models with covariates
- [Julia Covariates](../../julia/population/covariates.md) - Detailed Julia documentation
- [Estimation](../estimation/index.md) - Estimating covariate effects
