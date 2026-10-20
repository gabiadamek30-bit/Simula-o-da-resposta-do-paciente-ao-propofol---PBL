# One-Compartment IV Bolus

Single-compartment pharmacokinetic model with instantaneous IV bolus or zero-order infusion administration.

---

## Function Signature

```python
neopkpd.simulate_pk_iv_bolus(
    cl: float,
    v: float,
    doses: list[dict],
    t0: float,
    t1: float,
    saveat: list[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
    alag: float | None = None,
    bioavailability: float | None = None,
) -> dict
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `cl` | float | Clearance (volume/time, e.g., L/h) |
| `v` | float | Volume of distribution (L) |
| `doses` | list[dict] | Dose events (see below) |
| `t0` | float | Simulation start time |
| `t1` | float | Simulation end time |
| `saveat` | list[float] | Time points for output |
| `alg` | str | ODE solver algorithm (default: "Tsit5") |
| `reltol` | float | Relative tolerance (default: 1e-10) |
| `abstol` | float | Absolute tolerance (default: 1e-12) |
| `maxiters` | int | Maximum solver iterations |
| `alag` | float | Absorption lag time (optional) |
| `bioavailability` | float | Fraction absorbed, 0-1 (optional) |

### Dose Event Format

```python
doses = [
    {"time": 0.0, "amount": 100.0},                    # IV bolus
    {"time": 8.0, "amount": 100.0, "duration": 1.0},   # 1-hour infusion
]
```

---

## Returns

```python
{
    "t": [0.0, 1.0, 2.0, ...],           # Time points
    "states": {
        "A": [100.0, 90.5, ...]          # Amount in compartment
    },
    "observations": {
        "conc": [10.0, 9.05, ...]        # Concentration (A/V)
    },
    "metadata": {...}
}
```

---

## Model Equations

$$\frac{dA}{dt} = -k \cdot A$$

$$C = \frac{A}{V}$$

$$k = \frac{CL}{V}$$

For infusion (duration > 0):
$$\frac{dA}{dt} = R_0 - k \cdot A$$

where $R_0 = \text{amount} / \text{duration}$

---

## Basic Examples

### Single IV Bolus

```python
import neopkpd

result = neopkpd.simulate_pk_iv_bolus(
    cl=10.0,      # 10 L/h clearance
    v=50.0,       # 50 L volume
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=24.0,
    saveat=[0, 1, 2, 4, 6, 8, 12, 24]
)

# Access results
print(f"C0 = {result['observations']['conc'][0]:.1f} mg/L")  # 10.0 mg/L
print(f"Half-life = {0.693 * 50 / 10:.1f} h")                # 3.47 h
```

### IV Infusion

```python
# 500 mg over 2 hours
result = neopkpd.simulate_pk_iv_bolus(
    cl=10.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 500.0, "duration": 2.0}],
    t0=0.0,
    t1=24.0,
    saveat=[i * 0.5 for i in range(49)]
)

# Peak occurs at end of infusion
conc = result['observations']['conc']
t = result['t']
cmax_idx = max(range(len(conc)), key=lambda i: conc[i])
print(f"Cmax = {conc[cmax_idx]:.2f} mg/L at t = {t[cmax_idx]} h")
```

---

## Multiple Dosing

### Repeated IV Boluses

```python
# 500 mg every 8 hours for 5 doses
doses = [{"time": i * 8.0, "amount": 500.0} for i in range(5)]

result = neopkpd.simulate_pk_iv_bolus(
    cl=10.0,
    v=50.0,
    doses=doses,
    t0=0.0,
    t1=48.0,
    saveat=[i * 0.5 for i in range(97)]
)

# Find Cmax and Cmin in last dosing interval
conc = result['observations']['conc']
t = result['t']

# Last interval: 32-40 h
last_start = next(i for i, x in enumerate(t) if x >= 32)
last_end = next(i for i, x in enumerate(t) if x >= 40)

cmax = max(conc[last_start:last_end])
cmin = min(conc[last_start:last_end])

print(f"Steady-state Cmax: {cmax:.2f} mg/L")
print(f"Steady-state Cmin: {cmin:.2f} mg/L")
```

### Continuous Infusion

```python
# Continuous infusion: 20 mg/h for 24 hours
result = neopkpd.simulate_pk_iv_bolus(
    cl=10.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 480.0, "duration": 24.0}],
    t0=0.0,
    t1=48.0,
    saveat=[i * 0.5 for i in range(97)]
)

# Steady-state concentration
css = 480.0 / 24.0 / 10.0  # Rate / CL = 2 mg/L
print(f"Theoretical Css: {css:.1f} mg/L")
print(f"Simulated C at 24h: {result['observations']['conc'][48]:.2f} mg/L")
```

---

## Loading Dose Strategy

```python
# Loading dose followed by maintenance infusion
cl = 10.0
v = 50.0
target_css = 5.0  # mg/L

# Loading dose to immediately reach target
loading_dose = target_css * v  # 250 mg

# Maintenance rate to sustain target
maint_rate = target_css * cl  # 50 mg/h

doses = [
    {"time": 0.0, "amount": loading_dose},                     # Loading bolus
    {"time": 0.0, "amount": maint_rate * 24, "duration": 24.0} # 24h infusion
]

result = neopkpd.simulate_pk_iv_bolus(
    cl=cl, v=v,
    doses=doses,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Concentration should stay near 5 mg/L throughout
print(f"C at 0h: {result['observations']['conc'][0]:.2f} mg/L")
print(f"C at 12h: {result['observations']['conc'][48]:.2f} mg/L")
print(f"C at 24h: {result['observations']['conc'][-1]:.2f} mg/L")
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pk_profile

result = neopkpd.simulate_pk_iv_bolus(
    cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Plot concentration-time profile
fig = plot_pk_profile(
    result['t'],
    result['observations']['conc'],
    title="One-Compartment IV Bolus",
    xlabel="Time (h)",
    ylabel="Concentration (mg/L)"
)
```

---

## Derived Parameters

### From Simulation Results

```python
import numpy as np

result = neopkpd.simulate_pk_iv_bolus(
    cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

conc = np.array(result['observations']['conc'])
t = np.array(result['t'])

# AUC (trapezoidal rule)
auc = np.trapz(conc, t)
print(f"AUC: {auc:.1f} mg*h/L")
print(f"Theoretical AUC (Dose/CL): {500/10:.1f} mg*h/L")

# Half-life from log-linear regression
log_conc = np.log(conc[conc > 0])
t_valid = t[conc > 0]
slope = np.polyfit(t_valid, log_conc, 1)[0]
half_life = -0.693 / slope
print(f"Half-life: {half_life:.2f} h")
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Elimination rate constant | $k = CL/V$ |
| Half-life | $t_{1/2} = 0.693/k$ |
| C(t) after bolus | $C_0 \cdot e^{-kt}$ |
| Css (infusion) | $R_0/CL$ |
| AUC (single dose) | $Dose/CL$ |

---

## See Also

- [One-Compartment Oral](oral.md) - First-order absorption
- [Two-Compartment IV](twocomp-iv.md) - With distribution
- [NCA Analysis](../../nca/index.md) - Non-compartmental analysis
- [Population Simulation](../../population/index.md) - Adding variability
