# Population Oral

Comprehensive guide for simulating oral absorption with inter-individual and inter-occasion variability.

---

## Overview

Population oral simulation models the variability in absorption, clearance, and volume across subjects following oral administration.

```python
from neopkpd import simulate_population_oral

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

## Function Signature

```python
def simulate_population_oral(
    ka: float,
    cl: float,
    v: float,
    doses: list[dict],
    t0: float,
    t1: float,
    saveat: float,
    n: int,
    omegas: dict[str, float] | None = None,
    omega_matrix: np.ndarray | None = None,
    omega_params: list[str] | None = None,
    seed: int | None = None,
    iov_pis: dict[str, float] | None = None,
    iov_seed: int | None = None,
    covariates: list[dict] | None = None,
    covariate_effects: list[dict] | None = None
) -> PopulationResult:
    """
    Simulate population oral pharmacokinetics.

    Parameters
    ----------
    ka : float
        Population absorption rate constant (/hr)
    cl : float
        Population clearance (L/hr)
    v : float
        Population volume of distribution (L)
    doses : list[dict]
        Dosing events [{"time": t, "amount": amt}, ...]
    t0, t1 : float
        Simulation time range
    saveat : float
        Output time step
    n : int
        Number of subjects
    omegas : dict, optional
        Diagonal omega values {"Ka": ω², "CL": ω², "V": ω²}
    omega_matrix : np.ndarray, optional
        Full omega covariance matrix
    omega_params : list[str], optional
        Parameter names for omega_matrix
    seed : int, optional
        Random seed for IIV
    iov_pis : dict, optional
        IOV variance {"Ka": π², "CL": π²}
    iov_seed : int, optional
        Random seed for IOV
    covariates : list[dict], optional
        Per-subject covariate values
    covariate_effects : list[dict], optional
        Covariate effect specifications

    Returns
    -------
    PopulationResult
        Simulation results
    """
```

---

## Basic Usage

### Single Dose

```python
from neopkpd import simulate_population_oral

result = simulate_population_oral(
    ka=1.5,         # Absorption rate (/hr)
    cl=10.0,        # Clearance (L/hr)
    v=50.0,         # Volume (L)
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=24.0,
    saveat=0.5,
    n=100,
    omegas={"Ka": 0.16, "CL": 0.09, "V": 0.04},
    seed=42
)

# Summary metrics
all_cmax = [max(ind.concentrations) for ind in result.individuals]
all_tmax = [result.times[ind.concentrations.index(max(ind.concentrations))]
            for ind in result.individuals]

print(f"Cmax: {np.mean(all_cmax):.2f} ± {np.std(all_cmax):.2f} mg/L")
print(f"Tmax: {np.mean(all_tmax):.2f} ± {np.std(all_tmax):.2f} hr")
```

### Multiple Doses

```python
# QD dosing for 7 days
doses = [{"time": i * 24.0, "amount": 500.0} for i in range(7)]

result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=doses,
    t0=0.0, t1=168.0, saveat=1.0,
    n=100,
    omegas={"Ka": 0.16, "CL": 0.09, "V": 0.04},
    seed=42
)

# Steady-state metrics (day 7)
ss_start = 144.0  # Day 7 start
for ind in result.individuals[:5]:
    ss_idx = [i for i, t in enumerate(result.times) if t >= ss_start]
    ss_conc = [ind.concentrations[i] for i in ss_idx]
    print(f"Cmax,ss={max(ss_conc):.2f}, Cmin,ss={min(ss_conc):.2f}")
```

---

## IIV on Absorption

### High Variability on Ka

Absorption rate often shows the highest variability:

```python
# Ka typically 30-60% CV
omegas = {
    "Ka": 0.25,    # ~50% CV - absorption highly variable
    "CL": 0.09,    # ~30% CV
    "V": 0.04      # ~20% CV
}

result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=200,
    omegas=omegas,
    seed=42
)

# Ka distribution
ka_vals = [p["Ka"] for p in result.individual_params]
print(f"Ka: {np.mean(ka_vals):.2f} ± {np.std(ka_vals):.2f} /hr")
print(f"Range: [{np.min(ka_vals):.2f}, {np.max(ka_vals):.2f}]")
```

### Impact on Cmax and Tmax

```python
# Correlate Ka with Cmax and Tmax
ka_vals = np.array([p["Ka"] for p in result.individual_params])
cmax_vals = np.array([max(ind.concentrations) for ind in result.individuals])
tmax_vals = np.array([result.times[ind.concentrations.index(max(ind.concentrations))]
                      for ind in result.individuals])

print(f"Correlation Ka-Cmax: {np.corrcoef(ka_vals, cmax_vals)[0,1]:.3f}")
print(f"Correlation Ka-Tmax: {np.corrcoef(ka_vals, tmax_vals)[0,1]:.3f}")
# Higher Ka → Higher Cmax, Lower Tmax
```

---

## Inter-Occasion Variability (IOV)

### Adding IOV

```python
# Multiple dose with IOV
doses = [{"time": i * 24.0, "amount": 500.0} for i in range(4)]

result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=doses,
    t0=0.0, t1=96.0, saveat=0.5,
    n=50,
    omegas={"Ka": 0.16, "CL": 0.09, "V": 0.04},   # IIV
    seed=42,
    iov_pis={"Ka": 0.04, "CL": 0.0225},           # IOV: 20% on Ka, 15% on CL
    iov_seed=12345
)
```

### IOV Effects

```python
# IOV causes occasion-to-occasion variation within subjects
# Extract Cmax for each occasion
for subj_idx in range(5):
    ind = result.individuals[subj_idx]
    print(f"Subject {subj_idx + 1}:")

    for occ in range(4):
        t_start = occ * 24.0
        t_end = (occ + 1) * 24.0
        occ_idx = [i for i, t in enumerate(result.times) if t_start <= t < t_end]
        occ_conc = [ind.concentrations[i] for i in occ_idx]
        cmax_occ = max(occ_conc)
        print(f"  Occasion {occ+1}: Cmax = {cmax_occ:.2f}")
```

---

## Two-Compartment Oral

```python
from neopkpd import simulate_population_twocomp_oral

result = simulate_population_twocomp_oral(
    ka=1.5,         # Absorption rate
    cl=10.0,        # Central clearance
    v1=50.0,        # Central volume
    q=5.0,          # Inter-compartmental clearance
    v2=100.0,       # Peripheral volume
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=72.0, saveat=0.5,
    n=100,
    omegas={"Ka": 0.16, "CL": 0.09, "V1": 0.04, "Q": 0.04, "V2": 0.04},
    seed=42
)
```

---

## With Covariates

### Allometric Scaling

```python
import numpy as np

# Generate population with weight distribution
np.random.seed(42)
n = 100
weights = np.random.normal(75, 15, n)
weights = np.clip(weights, 45, 150)

covariates = [{"WT": w} for w in weights]

covariate_effects = [
    {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0}
]

result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=n,
    omegas={"Ka": 0.16, "CL": 0.04, "V": 0.02},  # Reduced IIV (covariates explain part)
    covariates=covariates,
    covariate_effects=covariate_effects,
    seed=42
)

# Verify weight effect on CL
cl_vals = [p["CL"] for p in result.individual_params]
print(f"Correlation WT-CL: {np.corrcoef(weights, cl_vals)[0,1]:.3f}")
```

### Multiple Covariates

```python
# Weight, age, and renal function
np.random.seed(42)
n = 200

covariates = []
for _ in range(n):
    wt = np.random.normal(75, 15)
    age = np.random.normal(50, 12)
    crcl = max(30, 120 - age * 0.7 + np.random.normal(0, 15))
    covariates.append({
        "WT": np.clip(wt, 45, 150),
        "AGE": np.clip(age, 18, 85),
        "CRCL": np.clip(crcl, 15, 150)
    })

covariate_effects = [
    # Allometric scaling
    {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0},
    # Age effect
    {"param": "CL", "cov": "AGE", "ref": 45.0, "kind": "LinearCovariate", "beta": -0.005},
    # Renal function
    {"param": "CL", "cov": "CRCL", "ref": 100.0, "kind": "LinearCovariate", "beta": 0.004}
]

result = simulate_population_oral(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=n,
    omegas={"Ka": 0.16, "CL": 0.0225, "V": 0.01},
    covariates=covariates,
    covariate_effects=covariate_effects,
    seed=42
)
```

---

## PK Metric Calculation

### Comprehensive Metrics

```python
import numpy as np

def calculate_pk_metrics(times, concentrations, dose):
    """Calculate NCA-like PK metrics."""
    conc = np.array(concentrations)
    t = np.array(times)

    # Cmax and Tmax
    cmax = np.max(conc)
    tmax_idx = np.argmax(conc)
    tmax = t[tmax_idx]

    # AUC (trapezoidal)
    auc = np.trapz(conc, t)

    # AUC last dosing interval (for multiple dose)
    if len(t) > 24:
        last_24_idx = t >= (t[-1] - 24)
        auc_tau = np.trapz(conc[last_24_idx], t[last_24_idx])
    else:
        auc_tau = auc

    # Terminal half-life (from last 3 half-lives)
    # Use points after Tmax where conc is declining
    terminal_idx = (t > tmax) & (conc > 0.05 * cmax)
    if sum(terminal_idx) > 3:
        log_conc = np.log(conc[terminal_idx])
        slope, _ = np.polyfit(t[terminal_idx], log_conc, 1)
        t_half = -np.log(2) / slope
    else:
        t_half = np.nan

    return {
        "Cmax": cmax,
        "Tmax": tmax,
        "AUC": auc,
        "AUC_tau": auc_tau,
        "t_half": t_half,
        "CL_F": dose / auc if auc > 0 else np.nan
    }

# Calculate for population
all_metrics = [
    calculate_pk_metrics(result.times, ind.concentrations, 500.0)
    for ind in result.individuals
]

# Summary
print("PK Metrics Summary:")
for metric in ["Cmax", "Tmax", "AUC", "t_half"]:
    vals = [m[metric] for m in all_metrics if not np.isnan(m[metric])]
    print(f"  {metric}: {np.mean(vals):.2f} ± {np.std(vals):.2f}")
```

---

## Complete Example

```python
from neopkpd import simulate_population_oral
import numpy as np
import matplotlib.pyplot as plt

# =========================================
# Comprehensive Population Oral Simulation
# =========================================

print("=== Population Oral PK ===\n")

# 1. Parameters
ka, cl, v = 1.5, 10.0, 50.0
print(f"--- Typical Parameters ---")
print(f"Ka = {ka} /hr, CL = {cl} L/hr, V = {v} L")
print(f"Tmax (typical) ≈ {np.log(ka*v/cl) / (ka - cl/v):.1f} hr")

# 2. IIV
omegas = {"Ka": 0.16, "CL": 0.09, "V": 0.04}
print(f"\n--- IIV ---")
for p, w in omegas.items():
    cv = np.sqrt(np.exp(w) - 1) * 100
    print(f"  {p}: CV ≈ {cv:.0f}%")

# 3. IOV
iov_pis = {"Ka": 0.04, "CL": 0.0225}
print(f"\n--- IOV ---")
for p, pi in iov_pis.items():
    cv = np.sqrt(np.exp(pi) - 1) * 100
    print(f"  {p}: CV ≈ {cv:.0f}%")

# 4. Multiple dose regimen
doses = [{"time": i * 24.0, "amount": 500.0} for i in range(5)]
print(f"\n--- Dosing ---")
print(f"  500 mg QD × 5 days")

# 5. Covariates
np.random.seed(42)
n = 150
covariates = []
for _ in range(n):
    wt = np.clip(np.random.normal(75, 15), 45, 150)
    age = np.clip(np.random.normal(50, 12), 18, 85)
    covariates.append({"WT": wt, "AGE": age})

covariate_effects = [
    {"param": "CL", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 0.75},
    {"param": "V", "cov": "WT", "ref": 70.0, "kind": "PowerCovariate", "beta": 1.0}
]

print(f"\n--- Covariates ---")
wts = [c["WT"] for c in covariates]
ages = [c["AGE"] for c in covariates]
print(f"  WT: {np.mean(wts):.1f} ± {np.std(wts):.1f} kg")
print(f"  AGE: {np.mean(ages):.1f} ± {np.std(ages):.1f} years")

# 6. Simulate
print(f"\n--- Simulation ---")
result = simulate_population_oral(
    ka=ka, cl=cl, v=v,
    doses=doses,
    t0=0.0, t1=120.0, saveat=0.5,
    n=n,
    omegas=omegas,
    iov_pis=iov_pis,
    iov_seed=12345,
    covariates=covariates,
    covariate_effects=covariate_effects,
    seed=42
)
print(f"Simulated {result.n_subjects} subjects over 5 days")

# 7. First dose metrics
print(f"\n--- First Dose Metrics ---")
first_dose_idx = [i for i, t in enumerate(result.times) if t <= 24]

first_cmax = []
first_tmax = []
for ind in result.individuals:
    conc = [ind.concentrations[i] for i in first_dose_idx]
    first_cmax.append(max(conc))
    first_tmax.append(result.times[first_dose_idx[conc.index(max(conc))]])

print(f"Cmax: {np.mean(first_cmax):.2f} ± {np.std(first_cmax):.2f} mg/L")
print(f"Tmax: {np.mean(first_tmax):.2f} ± {np.std(first_tmax):.2f} hr")
print(f"CV(Cmax): {np.std(first_cmax)/np.mean(first_cmax)*100:.1f}%")

# 8. Steady-state metrics (Day 5)
print(f"\n--- Steady-State Metrics (Day 5) ---")
ss_idx = [i for i, t in enumerate(result.times) if 96 <= t <= 120]

ss_cmax = []
ss_cmin = []
ss_cavg = []

for ind in result.individuals:
    conc = [ind.concentrations[i] for i in ss_idx]
    ss_cmax.append(max(conc))
    ss_cmin.append(min(conc))
    ss_cavg.append(np.mean(conc))

print(f"Cmax,ss: {np.mean(ss_cmax):.2f} ± {np.std(ss_cmax):.2f} mg/L")
print(f"Cmin,ss: {np.mean(ss_cmin):.2f} ± {np.std(ss_cmin):.2f} mg/L")
print(f"Cavg,ss: {np.mean(ss_cavg):.2f} ± {np.std(ss_cavg):.2f} mg/L")
print(f"Fluctuation: {(np.mean(ss_cmax) - np.mean(ss_cmin))/np.mean(ss_cavg)*100:.1f}%")

# 9. Within-subject variability (IOV effect)
print(f"\n--- Within-Subject Variability ---")
cv_within = []
for ind in result.individuals:
    occasion_cmax = []
    for occ in range(5):
        t_start, t_end = occ * 24, (occ + 1) * 24
        occ_idx = [i for i, t in enumerate(result.times) if t_start <= t < t_end]
        if occ_idx:
            occ_conc = [ind.concentrations[i] for i in occ_idx]
            occasion_cmax.append(max(occ_conc))
    if len(occasion_cmax) > 1:
        cv_within.append(np.std(occasion_cmax) / np.mean(occasion_cmax) * 100)

print(f"Mean within-subject CV(Cmax): {np.mean(cv_within):.1f}%")

# 10. Visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 10a. All profiles
ax = axes[0, 0]
for ind in result.individuals[:30]:
    ax.plot(result.times, ind.concentrations, 'b-', alpha=0.2)
ax.plot(result.times, result.median, 'r-', linewidth=2, label='Median')
ax.set_xlabel('Time (hr)')
ax.set_ylabel('Concentration (mg/L)')
ax.set_title('Individual Profiles (n=30)')
ax.legend()

# 10b. Population summary
ax = axes[0, 1]
ax.fill_between(result.times, result.percentiles[5], result.percentiles[95],
                alpha=0.2, label='90% PI')
ax.fill_between(result.times, result.percentiles[25], result.percentiles[75],
                alpha=0.4, label='50% PI')
ax.plot(result.times, result.median, 'b-', linewidth=2, label='Median')
ax.set_xlabel('Time (hr)')
ax.set_ylabel('Concentration (mg/L)')
ax.set_title('Population Summary')
ax.legend()

# 10c. Cmax distribution
ax = axes[1, 0]
ax.hist(ss_cmax, bins=25, density=True, alpha=0.7)
ax.axvline(np.mean(ss_cmax), color='r', linestyle='--', label='Mean')
ax.axvline(np.median(ss_cmax), color='g', linestyle='--', label='Median')
ax.set_xlabel('Cmax,ss (mg/L)')
ax.set_ylabel('Density')
ax.set_title('Steady-State Cmax Distribution')
ax.legend()

# 10d. Weight vs CL
ax = axes[1, 1]
cl_vals = [p["CL"] for p in result.individual_params]
ax.scatter(wts, cl_vals, alpha=0.5)
ax.set_xlabel('Weight (kg)')
ax.set_ylabel('CL (L/hr)')
ax.set_title(f'Weight-CL Relationship (r={np.corrcoef(wts, cl_vals)[0,1]:.3f})')

plt.tight_layout()
plt.savefig('population_oral.png', dpi=150)
plt.show()

print("\n✓ Simulation complete")
```

---

## See Also

- [Population IV Bolus](iv-bolus.md) - IV bolus models
- [Covariates](covariates.md) - Covariate modeling
- [Julia IOV](../../julia/population/iov.md) - Detailed IOV documentation
- [NCA](../nca/index.md) - Non-compartmental analysis
