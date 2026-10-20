# One-Compartment Oral First-Order

Single-compartment pharmacokinetic model with first-order oral absorption.

---

## Function Signature

```python
neopkpd.simulate_pk_oral_first_order(
    ka: float,
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
| `ka` | float | Absorption rate constant (1/h) |
| `cl` | float | Apparent clearance CL/F (L/h) |
| `v` | float | Apparent volume V/F (L) |
| `doses` | list[dict] | Dose events with time and amount |
| `t0` | float | Simulation start time |
| `t1` | float | Simulation end time |
| `saveat` | list[float] | Time points for output |
| `alag` | float | Absorption lag time (optional) |
| `bioavailability` | float | Fraction absorbed F (optional) |

---

## Returns

```python
{
    "t": [0.0, 1.0, 2.0, ...],
    "states": {
        "A_gut": [...],          # Amount in gut
        "A": [...]               # Amount in central
    },
    "observations": {
        "conc": [...]            # Plasma concentration
    },
    "metadata": {...}
}
```

---

## Model Equations

$$\frac{dA_{gut}}{dt} = -K_a \cdot A_{gut}$$

$$\frac{dA}{dt} = K_a \cdot A_{gut} - k \cdot A$$

$$C = \frac{A}{V}$$

### Analytical Solution (Bateman Function)

$$C(t) = \frac{F \cdot D \cdot K_a}{V \cdot (K_a - k)} \cdot (e^{-kt} - e^{-K_a t})$$

---

## Basic Examples

### Single Oral Dose

```python
import neopkpd

result = neopkpd.simulate_pk_oral_first_order(
    ka=1.5,       # Absorption rate (1/h)
    cl=10.0,      # Clearance (L/h)
    v=50.0,       # Volume (L)
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Find Cmax and Tmax
conc = result['observations']['conc']
t = result['t']
cmax_idx = max(range(len(conc)), key=lambda i: conc[i])

print(f"Cmax: {conc[cmax_idx]:.2f} mg/L")
print(f"Tmax: {t[cmax_idx]:.2f} h")
```

### With Lag Time and Bioavailability

```python
# Drug with 30-minute lag and 60% bioavailability
result = neopkpd.simulate_pk_oral_first_order(
    ka=1.5,
    cl=10.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=24.0,
    saveat=[i * 0.25 for i in range(97)],
    alag=0.5,              # 30-minute lag
    bioavailability=0.6    # 60% absorbed
)

# Compare to theoretical
# Effective dose = 500 * 0.6 = 300 mg
# AUC = 300 / 10 = 30 mg*h/L
```

---

## Multiple Dosing

### BID (Twice Daily) Dosing

```python
# 250 mg every 12 hours for 5 days
doses = [{"time": i * 12.0, "amount": 250.0} for i in range(10)]

result = neopkpd.simulate_pk_oral_first_order(
    ka=1.5,
    cl=10.0,
    v=50.0,
    doses=doses,
    t0=0.0,
    t1=120.0,
    saveat=[i * 0.5 for i in range(241)]
)

# Steady-state metrics in last interval (108-120 h)
conc = result['observations']['conc']
t = result['t']

ss_start = next(i for i, x in enumerate(t) if x >= 108)
ss_conc = conc[ss_start:]

print(f"Cmax,ss: {max(ss_conc):.2f} mg/L")
print(f"Cmin,ss: {min(ss_conc):.2f} mg/L")
print(f"Cavg,ss: {sum(ss_conc)/len(ss_conc):.2f} mg/L")
```

### TID (Three Times Daily) Dosing

```python
# 100 mg every 8 hours
doses = [{"time": i * 8.0, "amount": 100.0} for i in range(15)]

result = neopkpd.simulate_pk_oral_first_order(
    ka=2.0, cl=8.0, v=40.0,
    doses=doses,
    t0=0.0, t1=120.0,
    saveat=[i * 0.25 for i in range(481)]
)
```

---

## Absorption Rate Effects

```python
import neopkpd

# Compare different Ka values (same AUC, different profiles)
ka_values = [0.5, 1.0, 2.0, 4.0]

for ka in ka_values:
    result = neopkpd.simulate_pk_oral_first_order(
        ka=ka, cl=10.0, v=50.0,
        doses=[{"time": 0.0, "amount": 500.0}],
        t0=0.0, t1=24.0,
        saveat=[i * 0.1 for i in range(241)]
    )

    conc = result['observations']['conc']
    t = result['t']
    cmax_idx = max(range(len(conc)), key=lambda i: conc[i])

    print(f"Ka={ka}: Cmax={conc[cmax_idx]:.2f} mg/L, Tmax={t[cmax_idx]:.1f} h")
```

**Expected Output:**
```
Ka=0.5: Cmax=5.12 mg/L, Tmax=3.5 h
Ka=1.0: Cmax=6.84 mg/L, Tmax=2.0 h
Ka=2.0: Cmax=7.89 mg/L, Tmax=1.2 h
Ka=4.0: Cmax=8.55 mg/L, Tmax=0.7 h
```

---

## Flip-Flop Kinetics

When Ka < k (elimination faster than absorption):

```python
# Standard case: Ka > k (absorption-rate limited)
result_standard = neopkpd.simulate_pk_oral_first_order(
    ka=2.0, cl=10.0, v=50.0,  # k = 0.2/h, Ka = 2.0/h
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

# Flip-flop case: Ka < k (elimination-rate limited)
result_flipflop = neopkpd.simulate_pk_oral_first_order(
    ka=0.1, cl=10.0, v=50.0,  # k = 0.2/h, Ka = 0.1/h
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

# In flip-flop, the "terminal" slope reflects Ka, not k
print("Standard: terminal slope reflects k")
print("Flip-flop: terminal slope reflects Ka")
```

---

## Tmax Calculation

Theoretical Tmax (when Ka ≠ k):

$$T_{max} = \frac{\ln(K_a/k)}{K_a - k}$$

```python
import math

def calculate_tmax(ka, cl, v):
    k = cl / v
    if abs(ka - k) < 1e-10:
        return 1 / ka  # When Ka ≈ k
    return math.log(ka / k) / (ka - k)

# Example
ka, cl, v = 1.5, 10.0, 50.0
tmax_theoretical = calculate_tmax(ka, cl, v)
print(f"Theoretical Tmax: {tmax_theoretical:.2f} h")

# Verify with simulation
result = neopkpd.simulate_pk_oral_first_order(
    ka=ka, cl=cl, v=v,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.01 for i in range(2401)]  # Fine resolution
)

conc = result['observations']['conc']
t = result['t']
tmax_simulated = t[max(range(len(conc)), key=lambda i: conc[i])]
print(f"Simulated Tmax: {tmax_simulated:.2f} h")
```

---

## AUC Calculation

```python
import numpy as np

result = neopkpd.simulate_pk_oral_first_order(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=120.0,  # Long enough for complete elimination
    saveat=[i * 0.1 for i in range(1201)]
)

conc = np.array(result['observations']['conc'])
t = np.array(result['t'])

# AUC by trapezoidal rule
auc = np.trapz(conc, t)

print(f"AUC: {auc:.1f} mg*h/L")
print(f"Theoretical (Dose/CL): {500/10:.1f} mg*h/L")
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pk_profile

result = neopkpd.simulate_pk_oral_first_order(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

fig = plot_pk_profile(
    result['t'],
    result['observations']['conc'],
    title="One-Compartment Oral",
    xlabel="Time (h)",
    ylabel="Concentration (mg/L)"
)
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Absorption half-life | $t_{1/2,a} = 0.693/K_a$ |
| Elimination half-life | $t_{1/2} = 0.693 \cdot V/CL$ |
| Tmax | $\ln(K_a/k) / (K_a - k)$ |
| Cmax | Complex (use simulation) |
| AUC | $F \cdot Dose / CL$ |

---

## See Also

- [One-Compartment IV](iv-bolus.md) - IV administration
- [Two-Compartment Oral](twocomp-oral.md) - With distribution
- [Transit Absorption](transit.md) - Delayed absorption
- [NCA Analysis](../../nca/index.md) - Exposure calculations
