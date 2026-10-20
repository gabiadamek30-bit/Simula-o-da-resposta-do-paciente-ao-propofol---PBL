# Population Modeling

Comprehensive guide for population pharmacokinetic and pharmacodynamic modeling in NeoPKPD Python.

---

## Overview

Population modeling accounts for variability between and within individuals, enabling realistic simulations and personalized predictions.

```python
from neopkpd import simulate_population_oral

# Simulate 100 subjects with IIV
result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omegas={"CL": 0.09, "V": 0.04, "Ka": 0.16},
    seed=42
)
```

---

## Documentation

<div class="grid cards" markdown>

-   :material-needle:{ .lg .middle } **Population IV Bolus**

    ---

    IV bolus simulation with IIV

    [:octicons-arrow-right-24: IV Bolus](iv-bolus.md)

-   :material-pill:{ .lg .middle } **Population Oral**

    ---

    Oral absorption with full variability

    [:octicons-arrow-right-24: Oral](oral.md)

-   :material-tune:{ .lg .middle } **Covariates**

    ---

    Weight, age, renal function effects

    [:octicons-arrow-right-24: Covariates](covariates.md)

</div>

---

## Key Concepts

### Inter-Individual Variability (IIV)

Individual parameters differ from population values:

$$P_i = \theta \cdot e^{\eta_i}$$

Where:
- $P_i$ = Individual parameter
- $\theta$ = Population typical value
- $\eta_i$ = Random effect, $\eta_i \sim N(0, \omega^2)$

```python
# Omega values define variability
# ω² = 0.09 corresponds to ~30% CV
omegas = {
    "CL": 0.09,   # 30% CV on clearance
    "V": 0.04,    # 20% CV on volume
    "Ka": 0.16    # 40% CV on absorption
}
```

### Omega to CV Relationship

| ω² | ω | CV (%) |
|----|---|--------|
| 0.01 | 0.10 | 10% |
| 0.04 | 0.20 | 20% |
| 0.09 | 0.30 | 30% |
| 0.16 | 0.40 | 40% |
| 0.25 | 0.50 | 50% |

For small ω: $CV \approx \omega$

For larger ω: $CV = \sqrt{e^{\omega^2} - 1}$

---

## Quick Start

### Basic Population Simulation

```python
from neopkpd import simulate_population_oral

# Define parameters
ka = 1.5        # Absorption rate (/hr)
cl = 10.0       # Clearance (L/hr)
v = 50.0        # Volume (L)

# Define dosing
doses = [{"time": 0.0, "amount": 500.0}]

# Define IIV
omegas = {"CL": 0.09, "V": 0.04, "Ka": 0.16}

# Simulate population
result = simulate_population_oral(
    ka=ka, cl=cl, v=v,
    doses=doses,
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,                    # 100 subjects
    omegas=omegas,
    seed=42                   # For reproducibility
)

# Access results
print(f"Subjects simulated: {len(result.individuals)}")
print(f"Time points: {len(result.times)}")
```

### Accessing Individual Results

```python
# Individual concentration profiles
for i, individual in enumerate(result.individuals[:5]):
    cmax = max(individual.concentrations)
    tmax_idx = individual.concentrations.index(cmax)
    tmax = result.times[tmax_idx]
    print(f"Subject {i+1}: Cmax={cmax:.2f} mg/L at Tmax={tmax:.1f} hr")

# Individual parameters
for i, params in enumerate(result.individual_params[:5]):
    print(f"Subject {i+1}: CL={params['CL']:.2f}, V={params['V']:.2f}, Ka={params['Ka']:.2f}")
```

### Population Summaries

```python
# Summary statistics at each time point
print(f"Mean concentrations: {result.mean}")
print(f"Median concentrations: {result.median}")
print(f"5th percentile: {result.percentiles[5]}")
print(f"95th percentile: {result.percentiles[95]}")

# Population metrics
import numpy as np
all_cmax = [max(ind.concentrations) for ind in result.individuals]
print(f"Mean Cmax: {np.mean(all_cmax):.2f} ± {np.std(all_cmax):.2f} mg/L")
print(f"CV(Cmax): {np.std(all_cmax)/np.mean(all_cmax)*100:.1f}%")
```

---

## Simulation Functions

### simulate_population_iv_bolus

```python
from neopkpd import simulate_population_iv_bolus

result = simulate_population_iv_bolus(
    cl=10.0,                    # Clearance (L/hr)
    v=50.0,                     # Volume (L)
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,                     # Start time
    t1=24.0,                    # End time
    saveat=0.5,                 # Time step
    n=100,                      # Number of subjects
    omegas={"CL": 0.09, "V": 0.04},
    seed=42
)
```

### simulate_population_oral

```python
from neopkpd import simulate_population_oral

result = simulate_population_oral(
    ka=1.5,                     # Absorption rate (/hr)
    cl=10.0,                    # Clearance (L/hr)
    v=50.0,                     # Volume (L)
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omegas={"CL": 0.09, "V": 0.04, "Ka": 0.16},
    seed=42,

    # Optional: IOV
    iov_pis={"CL": 0.0225, "Ka": 0.01},
    iov_seed=12345,

    # Optional: Covariates
    covariates=[{"WT": 80.0, "AGE": 45.0} for _ in range(100)],
    covariate_effects=[
        {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75}
    ]
)
```

---

## PopulationResult Structure

```python
@dataclass
class PopulationResult:
    """Result of population simulation."""

    # Time grid
    times: list[float]

    # Individual results
    individuals: list[IndividualResult]
    individual_params: list[dict[str, float]]
    etas: dict[str, list[float]]

    # Summary statistics
    mean: list[float]
    median: list[float]
    sd: list[float]
    percentiles: dict[int, list[float]]

    # Metadata
    n_subjects: int
    seed: int
```

---

## Correlated Random Effects

### Full Omega Matrix

```python
import numpy as np
from neopkpd import simulate_population_oral

# Define correlation between CL and V
# Variance: CL=0.09, V=0.04
# Covariance: 0.03 (correlation ≈ 0.5)
omega_matrix = np.array([
    [0.09, 0.03],   # CL
    [0.03, 0.04]    # V
])

result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omega_matrix=omega_matrix,
    omega_params=["CL", "V"],
    omegas={"Ka": 0.16},  # Ka uncorrelated
    seed=42
)

# Check correlation in realized parameters
cl_vals = [p["CL"] for p in result.individual_params]
v_vals = [p["V"] for p in result.individual_params]
print(f"Correlation CL-V: {np.corrcoef(cl_vals, v_vals)[0,1]:.3f}")
```

---

## Multiple Dose Simulations

### Steady-State Simulation

```python
# Multiple dose regimen (7 days QD)
doses = [
    {"time": i * 24.0, "amount": 500.0}
    for i in range(7)
]

result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=doses,
    t0=0.0, t1=168.0, saveat=1.0,
    n=100,
    omegas={"CL": 0.09, "V": 0.04, "Ka": 0.16},
    seed=42
)

# Extract steady-state metrics (last dosing interval)
ss_start_idx = result.times.index(144.0)  # Day 7
for ind in result.individuals[:5]:
    ss_conc = ind.concentrations[ss_start_idx:]
    cmax_ss = max(ss_conc)
    cmin_ss = min(ss_conc)
    print(f"Cmax,ss={cmax_ss:.2f}, Cmin,ss={cmin_ss:.2f}")
```

---

## Complete Example

```python
from neopkpd import simulate_population_oral
import numpy as np

# =========================================
# Comprehensive Population Simulation
# =========================================

print("=== Population PK Simulation ===\n")

# 1. Model parameters
print("--- Model Parameters ---")
ka = 1.5       # /hr
cl = 10.0      # L/hr
v = 50.0       # L
print(f"Ka = {ka} /hr, CL = {cl} L/hr, V = {v} L")
print(f"t1/2 = {0.693 * v / cl:.1f} hr")

# 2. IIV specification
print("\n--- Inter-Individual Variability ---")
omegas = {"CL": 0.09, "V": 0.04, "Ka": 0.16}
for param, omega in omegas.items():
    cv = np.sqrt(np.exp(omega) - 1) * 100
    print(f"  {param}: ω² = {omega}, CV ≈ {cv:.0f}%")

# 3. Covariate effects
print("\n--- Covariate Model ---")
covariate_effects = [
    {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0},
    {"param": "CL", "cov": "AGE", "ref": 45.0, "kind": "LinearCovariate", "beta": -0.005}
]
for eff in covariate_effects:
    print(f"  {eff['param']} ~ {eff['cov']}: {eff['kind']} (β={eff['beta']})")

# 4. Generate population covariates
np.random.seed(42)
n = 200
covariates = []
for _ in range(n):
    wt = np.random.normal(75, 15)
    age = np.random.normal(50, 12)
    covariates.append({"WT": max(45, min(150, wt)), "AGE": max(18, min(85, age))})

wts = [c["WT"] for c in covariates]
ages = [c["AGE"] for c in covariates]
print(f"\n--- Population Covariates ---")
print(f"  WT: {np.mean(wts):.1f} ± {np.std(wts):.1f} kg")
print(f"  AGE: {np.mean(ages):.1f} ± {np.std(ages):.1f} years")

# 5. Dosing
doses = [{"time": 0.0, "amount": 500.0}]

# 6. Simulate
print("\n--- Simulation ---")
result = simulate_population_oral(
    ka=ka, cl=cl, v=v,
    doses=doses,
    t0=0.0, t1=48.0, saveat=0.5,
    n=n,
    omegas=omegas,
    covariates=covariates,
    covariate_effects=covariate_effects,
    seed=12345
)
print(f"Simulated {result.n_subjects} subjects")

# 7. PK metrics
print("\n--- PK Metrics ---")
all_cmax = [max(ind.concentrations) for ind in result.individuals]
all_tmax = [result.times[ind.concentrations.index(max(ind.concentrations))]
            for ind in result.individuals]

# Approximate AUC
def trapezoidal_auc(times, conc):
    auc = 0.0
    for i in range(1, len(times)):
        auc += 0.5 * (conc[i] + conc[i-1]) * (times[i] - times[i-1])
    return auc

all_auc = [trapezoidal_auc(result.times, ind.concentrations)
           for ind in result.individuals]

print(f"Cmax: {np.mean(all_cmax):.2f} ± {np.std(all_cmax):.2f} mg/L")
print(f"Tmax: {np.mean(all_tmax):.2f} ± {np.std(all_tmax):.2f} hr")
print(f"AUC:  {np.mean(all_auc):.1f} ± {np.std(all_auc):.1f} mg·hr/L")

print(f"\nCV(Cmax): {np.std(all_cmax)/np.mean(all_cmax)*100:.1f}%")
print(f"CV(AUC):  {np.std(all_auc)/np.mean(all_auc)*100:.1f}%")

# 8. Percentiles
print("\n--- Population Percentiles ---")
for pct in [5, 25, 50, 75, 95]:
    cmax_pct = np.percentile(all_cmax, pct)
    auc_pct = np.percentile(all_auc, pct)
    print(f"  {pct}th: Cmax={cmax_pct:.2f}, AUC={auc_pct:.1f}")

# 9. Covariate impact
print("\n--- Covariate Impact ---")
# By weight
low_wt = [all_auc[i] for i in range(n) if covariates[i]["WT"] < 60]
high_wt = [all_auc[i] for i in range(n) if covariates[i]["WT"] > 90]
print(f"AUC (WT < 60 kg): {np.mean(low_wt):.1f} mg·hr/L")
print(f"AUC (WT > 90 kg): {np.mean(high_wt):.1f} mg·hr/L")

# By age
young = [all_auc[i] for i in range(n) if covariates[i]["AGE"] < 35]
old = [all_auc[i] for i in range(n) if covariates[i]["AGE"] > 65]
print(f"AUC (AGE < 35): {np.mean(young):.1f} mg·hr/L")
print(f"AUC (AGE > 65): {np.mean(old):.1f} mg·hr/L")
```

---

## Next Steps

- [Population IV Bolus](iv-bolus.md) - IV bolus with IIV
- [Population Oral](oral.md) - Oral absorption models
- [Covariates](covariates.md) - Covariate modeling details
- [Estimation](../estimation/index.md) - Fitting population models
- [VPC](../viz/vpc.md) - Visual predictive checks
