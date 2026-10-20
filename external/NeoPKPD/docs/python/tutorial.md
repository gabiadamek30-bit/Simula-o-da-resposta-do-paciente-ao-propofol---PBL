# Tutorial: Getting Started with Python

This tutorial provides a comprehensive introduction to using NeoPKPD with Python. By the end, you'll be able to run simulations, analyze data, and create professional visualizations.

---

## Prerequisites

Ensure you have:

- Python 3.10 or later
- Julia 1.10 or later installed
- NeoPKPD Python package installed

```bash
cd packages/python
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[all]"
```

---

## Part 1: Your First Simulation

### Initialization

```python
import neopkpd

# Initialize Julia runtime (required once per session)
neopkpd.init_julia()

# Check version
print(neopkpd.version())  # "0.1.0"
```

!!! note "First-time initialization"
    The first call to `init_julia()` takes 2-5 seconds due to Julia JIT compilation. Subsequent calls are instant.

### One-Compartment IV Bolus

```python
# Simulate IV bolus: 100 mg at time 0
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0,              # Clearance (L/h)
    v=50.0,              # Volume of distribution (L)
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,              # Start time
    t1=24.0,             # End time
    saveat=[float(t) for t in range(25)]  # Output every hour
)

# Access results
print("Time points:", result["t"])
print("Concentrations:", result["observations"]["conc"][:5])
```

**Expected output:**
```
Time points: [0.0, 1.0, 2.0, ..., 24.0]
Concentrations: [2.0, 1.809..., 1.637..., 1.481..., 1.340...]
```

### Understanding the Result

```python
# Result structure
result.keys()  # dict_keys(['t', 'states', 'observations', 'metadata'])

# Time points
result["t"]  # [0.0, 1.0, 2.0, ..., 24.0]

# Observations (concentrations)
result["observations"]["conc"]  # [2.0, 1.81, ...]

# State variables (amount in compartment)
result["states"]["A_central"]  # [100.0, 90.5, ...]

# Metadata
result["metadata"]  # {"model": "OneCompIVBolus", ...}
```

### Calculate PK Metrics

```python
import numpy as np

conc = np.array(result["observations"]["conc"])
t = np.array(result["t"])

# Cmax and Tmax
cmax = np.max(conc)
tmax = t[np.argmax(conc)]
print(f"Cmax: {cmax:.2f} mg/L at Tmax: {tmax:.1f} h")

# AUC by trapezoidal rule
auc = np.trapz(conc, t)
print(f"AUC0-24: {auc:.2f} mg·h/L")

# Half-life
t_half = 0.693 * 50.0 / 5.0  # t½ = 0.693 * V / CL
print(f"Half-life: {t_half:.2f} h")
```

---

## Part 2: Multiple Doses

### Repeated Dosing

```python
# 100 mg every 12 hours for 3 days
doses = [
    {"time": 0.0, "amount": 100.0},
    {"time": 12.0, "amount": 100.0},
    {"time": 24.0, "amount": 100.0},
    {"time": 36.0, "amount": 100.0},
    {"time": 48.0, "amount": 100.0},
    {"time": 60.0, "amount": 100.0},
]

result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=doses,
    t0=0.0, t1=72.0,
    saveat=[t * 0.5 for t in range(145)]  # Every 30 min
)

# Find steady-state trough (just before 6th dose)
conc = np.array(result["observations"]["conc"])
t = np.array(result["t"])
idx_trough = np.argmin(np.abs(t - 59.5))
print(f"Steady-state trough: {conc[idx_trough]:.2f} mg/L")
```

### IV Infusion

```python
# 100 mg infused over 1 hour
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0, "duration": 1.0}],
    t0=0.0, t1=24.0,
    saveat=[t * 0.25 for t in range(97)]
)

# Cmax occurs at end of infusion
conc = np.array(result["observations"]["conc"])
t = np.array(result["t"])
idx_1h = np.argmin(np.abs(t - 1.0))
print(f"Cmax at end of infusion: {conc[idx_1h]:.2f} mg/L")
```

---

## Part 3: Different Models

### Oral First-Order Absorption

```python
result = neopkpd.simulate_pk_oral_first_order(
    ka=1.5,              # Absorption rate constant (1/h)
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 200.0}],
    t0=0.0, t1=24.0,
    saveat=[t * 0.25 for t in range(97)]
)

conc = np.array(result["observations"]["conc"])
t = np.array(result["t"])

# Tmax is delayed for oral
tmax = t[np.argmax(conc)]
cmax = np.max(conc)
print(f"Oral: Cmax = {cmax:.2f} mg/L at Tmax = {tmax:.2f} h")
```

### Two-Compartment IV

```python
result = neopkpd.simulate_pk_twocomp_iv_bolus(
    cl=10.0,             # Clearance (L/h)
    v1=20.0,             # Central volume (L)
    q=15.0,              # Inter-compartmental clearance (L/h)
    v2=50.0,             # Peripheral volume (L)
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=48.0,
    saveat=[t * 0.1 for t in range(481)]
)

# Bi-exponential decline
conc = np.array(result["observations"]["conc"])
print(f"C at 0.5h: {conc[5]:.2f} (distribution)")
print(f"C at 24h: {conc[240]:.2f} (terminal)")
```

### Transit Absorption

```python
result = neopkpd.simulate_pk_transit_absorption(
    n_transit=5,         # 5 transit compartments
    ktr=0.5,             # Transit rate (1/h)
    ka=2.0,              # Final absorption rate (1/h)
    cl=10.0,
    v=70.0,
    doses=[{"time": 0.0, "amount": 300.0}],
    t0=0.0, t1=24.0,
    saveat=[t * 0.1 for t in range(241)]
)

# Delayed, broader peak
conc = np.array(result["observations"]["conc"])
t = np.array(result["t"])
tmax = t[np.argmax(conc)]
print(f"Transit absorption Tmax: {tmax:.2f} h")
```

---

## Part 4: PK-PD Models

### Direct Emax

```python
result = neopkpd.simulate_pkpd_direct_emax(
    # PK parameters
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    # PD parameters
    e0=0.0,              # Baseline effect
    emax=100.0,          # Maximum effect
    ec50=2.0,            # EC50 (mg/L)
    # Simulation grid
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)]
)

# Both concentration and effect
print("Concentration:", result["observations"]["conc"][:5])
print("Effect:", result["observations"]["effect"][:5])
```

### Indirect Response

```python
result = neopkpd.simulate_pkpd_indirect_response(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    kin=10.0,            # Production rate
    kout=0.5,            # Elimination rate
    ic50=2.0,            # IC50
    imax=0.9,            # Max inhibition
    baseline=None,       # Auto-calculate: kin/kout
    t0=0.0, t1=72.0,
    saveat=[float(t) for t in range(73)]
)

# Delayed response
response = result["observations"]["response"]
print(f"Baseline: {response[0]:.1f}")
print(f"Max suppression (at ~12h): {min(response):.1f}")
print(f"Return to baseline: {response[-1]:.1f}")
```

---

## Part 5: Population Simulation

### Inter-Individual Variability

```python
result = neopkpd.simulate_population_iv_bolus(
    cl=5.0,              # Typical CL
    v=50.0,              # Typical V
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=100,               # 100 subjects
    seed=12345,          # Reproducibility
    omegas={             # IIV (coefficient of variation)
        "CL": 0.3,       # 30% CV on CL
        "V": 0.2         # 20% CV on V
    }
)

# Access individual results
print(f"Number of subjects: {len(result['individuals'])}")

for i in range(3):
    cmax = max(result["individuals"][i]["observations"]["conc"])
    print(f"Subject {i+1} Cmax: {cmax:.2f} mg/L")

# Access realized parameters
for i in range(3):
    cl = result["params"][i]["CL"]
    v = result["params"][i]["V"]
    print(f"Subject {i+1}: CL={cl:.2f}, V={v:.2f}")

# Population summary
summary = result["summaries"]["conc"]
print(f"\nPopulation median Cmax: {max(summary['median']):.2f}")
print(f"90% PI: {max(summary['quantiles']['0.05']):.2f} - "
      f"{max(summary['quantiles']['0.95']):.2f}")
```

### Convert to DataFrame

```python
import pandas as pd

# Create long-format concentration DataFrame
data = []
for i, ind in enumerate(result["individuals"]):
    cl = result["params"][i]["CL"]
    v = result["params"][i]["V"]
    for t, c in zip(ind["t"], ind["observations"]["conc"]):
        data.append({
            "id": i + 1,
            "time": t,
            "conc": c,
            "CL": cl,
            "V": v
        })

df = pd.DataFrame(data)
print(df.head(10))
print(df.groupby("time")["conc"].describe())
```

---

## Part 6: Visualization

### Setup

```python
from neopkpd import viz

# Set backend
viz.set_backend("matplotlib")

# Check available backends
print(viz.available_backends())  # ["matplotlib", "plotly"]
```

### Concentration-Time Plot

```python
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[t * 0.5 for t in range(49)]
)

fig = viz.plot_conc_time(result, title="One-Compartment IV Bolus")
fig.savefig("conc_time.png", dpi=300)
```

### Population Spaghetti Plot

```python
pop_result = neopkpd.simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=50, seed=42,
    omegas={"CL": 0.3, "V": 0.2}
)

fig = viz.plot_spaghetti(pop_result, alpha=0.3)
fig.savefig("spaghetti.png", dpi=300)
```

### Mean with Confidence Ribbon

```python
fig = viz.plot_mean_ribbon(
    pop_result,
    ci_levels=[0.05, 0.95],
    show_median=True
)
fig.savefig("mean_ribbon.png", dpi=300)
```

### Interactive Plotly

```python
viz.set_backend("plotly")

fig = viz.plot_conc_time(result)
fig.write_html("conc_time_interactive.html")
fig.show()  # Opens in browser
```

---

## Part 7: Non-Compartmental Analysis

### Basic NCA

```python
from neopkpd.nca import run_nca, NCAConfig

# Concentration-time data
times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
conc = [0.0, 1.8, 2.0, 1.5, 1.0, 0.5, 0.25, 0.06]
dose = 100.0

# Run NCA
result = run_nca(times, conc, dose)

print(f"Cmax: {result.cmax:.2f} mg/L")
print(f"Tmax: {result.tmax:.2f} h")
print(f"AUC0-t: {result.auc_0_t:.2f} mg·h/L")
print(f"AUC0-inf: {result.auc_0_inf:.2f} mg·h/L")
print(f"t½: {result.t_half:.2f} h")
print(f"CL/F: {result.cl_f:.2f} L/h")
```

### NCA with Configuration

```python
config = NCAConfig(
    method="log_linear",
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9,
    extrapolation_max_pct=20.0
)

result = run_nca(times, conc, dose, config=config)
print(f"λz R²: {result.lambda_z_r_squared:.4f}")
```

---

## Part 8: Clinical Trial Simulation

### Power Analysis

```python
from neopkpd import trial

# Calculate power for given sample size
power = trial.estimate_power_analytical(
    n_per_arm=50,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05
)
print(f"Power: {power.power:.1%}")

# Calculate required sample size
sample = trial.estimate_sample_size(
    target_power=0.80,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05
)
print(f"Required n per arm: {sample.n_per_arm}")
```

### Generate Virtual Population

```python
pop = trial.generate_virtual_population(
    n=100,
    spec=trial.healthy_volunteer_spec(),
    seed=42
)

summary = trial.summarize_population(pop)
print(f"Age: {summary['age']['mean']:.1f} years")
print(f"Weight: {summary['weight']['mean']:.1f} kg")
```

---

## Summary

In this tutorial, you learned:

1. **Initialization**: `neopkpd.init_julia()`
2. **Simulation Functions**: IV bolus, oral, two-compartment, transit
3. **PK-PD Models**: Direct Emax, indirect response
4. **Population**: IIV with omega, summary statistics
5. **Visualization**: Matplotlib/Plotly dual backend
6. **NCA**: Exposure metrics, configuration
7. **Trials**: Power analysis, virtual populations

---

## Next Steps

- [Models Reference](models/index.md) - All simulation functions
- [Visualization Guide](viz/index.md) - 55+ plotting functions
- [NCA Reference](nca/index.md) - Detailed NCA documentation
- [Trial Module](trial/index.md) - Clinical trial simulation
