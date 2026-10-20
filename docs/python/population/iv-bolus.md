# Population IV Bolus

Comprehensive guide for simulating IV bolus administration with inter-individual variability.

---

## Overview

Population IV bolus simulation models the variability in clearance and volume across subjects following intravenous bolus administration.

```python
from neopkpd import simulate_population_iv_bolus

result = simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omegas={"CL": 0.09, "V": 0.04},
    seed=42
)
```

---

## Function Signature

```python
def simulate_population_iv_bolus(
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
    seed: int | None = None
) -> PopulationResult:
    """
    Simulate population IV bolus pharmacokinetics.

    Parameters
    ----------
    cl : float
        Population clearance (L/hr)
    v : float
        Population volume of distribution (L)
    doses : list[dict]
        Dosing events [{"time": t, "amount": amt}, ...]
    t0 : float
        Simulation start time
    t1 : float
        Simulation end time
    saveat : float
        Output time step
    n : int
        Number of subjects to simulate
    omegas : dict, optional
        Diagonal omega values {"CL": ω²_CL, "V": ω²_V}
    omega_matrix : np.ndarray, optional
        Full omega covariance matrix
    omega_params : list[str], optional
        Parameter names for omega_matrix rows/columns
    seed : int, optional
        Random seed for reproducibility

    Returns
    -------
    PopulationResult
        Simulation results with individual profiles and summaries
    """
```

---

## Basic Usage

### Single Dose

```python
from neopkpd import simulate_population_iv_bolus

# Simulate 100 subjects
result = simulate_population_iv_bolus(
    cl=5.0,                     # 5 L/hr clearance
    v=50.0,                     # 50 L volume
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=0.5,
    n=100,
    omegas={"CL": 0.09, "V": 0.04},  # 30% CV on CL, 20% on V
    seed=42
)

# Basic output
print(f"Subjects: {result.n_subjects}")
print(f"Time points: {len(result.times)}")
```

### Multiple Doses

```python
# Multiple IV bolus doses
doses = [
    {"time": 0.0, "amount": 100.0},
    {"time": 12.0, "amount": 100.0},
    {"time": 24.0, "amount": 100.0},
    {"time": 36.0, "amount": 100.0}
]

result = simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=doses,
    t0=0.0, t1=48.0, saveat=0.5,
    n=100,
    omegas={"CL": 0.09, "V": 0.04},
    seed=42
)
```

---

## IIV Specification

### Diagonal Omega (Uncorrelated)

```python
# Independent variability on CL and V
omegas = {
    "CL": 0.09,    # ω²_CL (30% CV)
    "V": 0.04      # ω²_V (20% CV)
}

result = simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omegas=omegas,
    seed=42
)
```

### Full Omega Matrix (Correlated)

```python
import numpy as np

# Correlated CL and V
# Correlation = 0.5, ω²_CL = 0.09, ω²_V = 0.04
# Covariance = 0.5 * sqrt(0.09) * sqrt(0.04) = 0.03
omega_matrix = np.array([
    [0.09, 0.03],   # CL row
    [0.03, 0.04]    # V row
])

result = simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omega_matrix=omega_matrix,
    omega_params=["CL", "V"],
    seed=42
)

# Verify correlation
cl_vals = [p["CL"] for p in result.individual_params]
v_vals = [p["V"] for p in result.individual_params]
print(f"Correlation: {np.corrcoef(cl_vals, v_vals)[0,1]:.3f}")
```

---

## Accessing Results

### Individual Profiles

```python
# Each subject's concentration-time profile
for i in range(5):
    ind = result.individuals[i]
    params = result.individual_params[i]

    cmax = max(ind.concentrations)
    auc = np.trapz(ind.concentrations, result.times)

    print(f"Subject {i+1}:")
    print(f"  CL = {params['CL']:.2f} L/hr")
    print(f"  V = {params['V']:.2f} L")
    print(f"  Cmax = {cmax:.2f} mg/L")
    print(f"  AUC = {auc:.1f} mg·hr/L")
```

### Eta Values

```python
# Random effects
print("Eta distributions:")
for param, etas in result.etas.items():
    print(f"  η_{param}: mean={np.mean(etas):.4f}, SD={np.std(etas):.4f}")

# Should be approximately N(0, ω)
```

### Population Summaries

```python
import matplotlib.pyplot as plt

# Plot population profile with percentiles
plt.figure(figsize=(10, 6))
plt.fill_between(result.times, result.percentiles[5], result.percentiles[95],
                 alpha=0.3, label='5-95th percentile')
plt.fill_between(result.times, result.percentiles[25], result.percentiles[75],
                 alpha=0.5, label='25-75th percentile')
plt.plot(result.times, result.median, 'b-', linewidth=2, label='Median')
plt.xlabel('Time (hr)')
plt.ylabel('Concentration (mg/L)')
plt.legend()
plt.yscale('log')
plt.title('Population IV Bolus Profile')
plt.show()
```

---

## Two-Compartment IV Bolus

```python
from neopkpd import simulate_population_twocomp_iv

# Two-compartment model with IIV
result = simulate_population_twocomp_iv(
    cl=10.0,        # Central clearance
    v1=50.0,        # Central volume
    q=5.0,          # Inter-compartmental clearance
    v2=100.0,       # Peripheral volume
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=72.0, saveat=0.5,
    n=100,
    omegas={"CL": 0.09, "V1": 0.04, "Q": 0.04, "V2": 0.04},
    seed=42
)
```

---

## PK Metric Calculation

### Basic Metrics

```python
import numpy as np

# Calculate PK metrics for each subject
pk_metrics = []

for i, ind in enumerate(result.individuals):
    conc = np.array(ind.concentrations)
    times = np.array(result.times)

    # C0 (initial concentration)
    c0 = conc[0]

    # Half-life (from slope of log-linear phase)
    # Use last portion of curve
    late_idx = times > times[-1] * 0.5
    if sum(late_idx) > 2:
        log_conc = np.log(conc[late_idx])
        slope = np.polyfit(times[late_idx], log_conc, 1)[0]
        half_life = -np.log(2) / slope
    else:
        half_life = np.nan

    # AUC (trapezoidal)
    auc = np.trapz(conc, times)

    # Clearance (from dose/AUC)
    dose = 100.0
    cl_calc = dose / auc

    pk_metrics.append({
        "C0": c0,
        "t_half": half_life,
        "AUC": auc,
        "CL_calc": cl_calc
    })

# Summary
c0_vals = [m["C0"] for m in pk_metrics]
thalf_vals = [m["t_half"] for m in pk_metrics if not np.isnan(m["t_half"])]
auc_vals = [m["AUC"] for m in pk_metrics]

print("PK Metrics Summary:")
print(f"C0: {np.mean(c0_vals):.2f} ± {np.std(c0_vals):.2f} mg/L")
print(f"t1/2: {np.mean(thalf_vals):.2f} ± {np.std(thalf_vals):.2f} hr")
print(f"AUC: {np.mean(auc_vals):.1f} ± {np.std(auc_vals):.1f} mg·hr/L")
```

---

## Complete Example

```python
from neopkpd import simulate_population_iv_bolus
import numpy as np
import matplotlib.pyplot as plt

# =========================================
# Population IV Bolus Simulation
# =========================================

print("=== Population IV Bolus ===\n")

# 1. Model parameters
cl = 5.0        # L/hr
v = 50.0        # L
dose = 100.0    # mg

print(f"--- Model Parameters ---")
print(f"CL = {cl} L/hr")
print(f"V = {v} L")
print(f"t1/2 = {0.693 * v / cl:.1f} hr")
print(f"Dose = {dose} mg")

# 2. IIV
omegas = {"CL": 0.09, "V": 0.04}
print(f"\n--- IIV ---")
for p, w in omegas.items():
    cv = np.sqrt(np.exp(w) - 1) * 100
    print(f"{p}: ω² = {w}, CV ≈ {cv:.0f}%")

# 3. Simulate
print(f"\n--- Simulation ---")
result = simulate_population_iv_bolus(
    cl=cl, v=v,
    doses=[{"time": 0.0, "amount": dose}],
    t0=0.0, t1=24.0, saveat=0.25,
    n=200,
    omegas=omegas,
    seed=42
)
print(f"Simulated {result.n_subjects} subjects")

# 4. Parameter distribution
print(f"\n--- Realized Parameters ---")
cl_vals = np.array([p["CL"] for p in result.individual_params])
v_vals = np.array([p["V"] for p in result.individual_params])

print(f"CL: {np.mean(cl_vals):.2f} ± {np.std(cl_vals):.2f} L/hr")
print(f"    Range: [{np.min(cl_vals):.2f}, {np.max(cl_vals):.2f}]")
print(f"V:  {np.mean(v_vals):.2f} ± {np.std(v_vals):.2f} L")
print(f"    Range: [{np.min(v_vals):.2f}, {np.max(v_vals):.2f}]")

# 5. PK metrics
print(f"\n--- PK Metrics ---")
c0_vals = [ind.concentrations[0] for ind in result.individuals]
auc_vals = [np.trapz(ind.concentrations, result.times) for ind in result.individuals]

print(f"C0: {np.mean(c0_vals):.2f} ± {np.std(c0_vals):.2f} mg/L")
print(f"CV(C0): {np.std(c0_vals)/np.mean(c0_vals)*100:.1f}%")
print(f"AUC: {np.mean(auc_vals):.1f} ± {np.std(auc_vals):.1f} mg·hr/L")
print(f"CV(AUC): {np.std(auc_vals)/np.mean(auc_vals)*100:.1f}%")

# 6. Percentiles
print(f"\n--- Population Percentiles (AUC) ---")
for pct in [5, 10, 25, 50, 75, 90, 95]:
    val = np.percentile(auc_vals, pct)
    print(f"  {pct}th: {val:.1f} mg·hr/L")

# 7. Visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 7a. Concentration-time profiles
ax = axes[0, 0]
for ind in result.individuals[:20]:
    ax.plot(result.times, ind.concentrations, 'b-', alpha=0.3)
ax.plot(result.times, result.median, 'r-', linewidth=2, label='Median')
ax.set_xlabel('Time (hr)')
ax.set_ylabel('Concentration (mg/L)')
ax.set_title('Individual Profiles')
ax.set_yscale('log')
ax.legend()

# 7b. Population summary with CI
ax = axes[0, 1]
ax.fill_between(result.times, result.percentiles[5], result.percentiles[95],
                alpha=0.3, color='blue', label='90% PI')
ax.plot(result.times, result.median, 'b-', linewidth=2, label='Median')
ax.set_xlabel('Time (hr)')
ax.set_ylabel('Concentration (mg/L)')
ax.set_title('Population Summary')
ax.legend()

# 7c. Parameter distributions
ax = axes[1, 0]
ax.hist(cl_vals, bins=20, density=True, alpha=0.7, label='CL')
ax.axvline(cl, color='r', linestyle='--', label=f'Typical = {cl}')
ax.set_xlabel('CL (L/hr)')
ax.set_ylabel('Density')
ax.set_title('CL Distribution')
ax.legend()

# 7d. AUC distribution
ax = axes[1, 1]
ax.hist(auc_vals, bins=20, density=True, alpha=0.7)
theoretical_auc = dose / cl
ax.axvline(theoretical_auc, color='r', linestyle='--',
           label=f'Typical = {theoretical_auc:.0f}')
ax.set_xlabel('AUC (mg·hr/L)')
ax.set_ylabel('Density')
ax.set_title('AUC Distribution')
ax.legend()

plt.tight_layout()
plt.savefig('population_iv_bolus.png', dpi=150)
plt.show()

print("\n✓ Simulation complete")
```

---

## Dose-Proportionality Analysis

```python
# Simulate multiple dose levels
dose_levels = [25, 50, 100, 200, 400]
results_by_dose = {}

for dose in dose_levels:
    result = simulate_population_iv_bolus(
        cl=5.0, v=50.0,
        doses=[{"time": 0.0, "amount": dose}],
        t0=0.0, t1=24.0, saveat=0.5,
        n=50,
        omegas={"CL": 0.09, "V": 0.04},
        seed=42
    )

    auc_vals = [np.trapz(ind.concentrations, result.times)
                for ind in result.individuals]
    results_by_dose[dose] = {
        "mean_auc": np.mean(auc_vals),
        "sd_auc": np.std(auc_vals)
    }

# Check linearity
print("Dose-Proportionality:")
print("Dose (mg)  |  Mean AUC  |  AUC/Dose")
print("-" * 40)
for dose in dose_levels:
    mean_auc = results_by_dose[dose]["mean_auc"]
    ratio = mean_auc / dose
    print(f"{dose:8}   | {mean_auc:9.1f}  | {ratio:9.2f}")
```

---

## See Also

- [Population Oral](oral.md) - Oral absorption models
- [Covariates](covariates.md) - Covariate effects
- [Julia IIV](../../julia/population/iiv.md) - Detailed IIV documentation
- [Estimation](../estimation/index.md) - Parameter estimation
