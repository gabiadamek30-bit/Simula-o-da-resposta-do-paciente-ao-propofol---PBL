# Three-Compartment IV Bolus

Three-compartment mammillary model with central, shallow peripheral, and deep peripheral compartments.

---

## Function Signature

```python
neopkpd.simulate_pk_threecomp_iv_bolus(
    cl: float,
    v1: float,
    q2: float,
    v2: float,
    q3: float,
    v3: float,
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
| `cl` | float | Clearance from central (L/h) |
| `v1` | float | Central compartment volume (L) |
| `q2` | float | Clearance to shallow peripheral (L/h) |
| `v2` | float | Shallow peripheral volume (L) |
| `q3` | float | Clearance to deep peripheral (L/h) |
| `v3` | float | Deep peripheral volume (L) |
| `doses` | list[dict] | Dose events |

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "A_central": [...],        # Amount in central
        "A_periph1": [...],        # Amount in shallow peripheral
        "A_periph2": [...]         # Amount in deep peripheral
    },
    "observations": {
        "conc": [...]              # Central concentration
    },
    "metadata": {...}
}
```

---

## Model Structure

- **Central (V1)**: Receives dose, site of elimination
- **Shallow Peripheral (V2)**: Rapid equilibration with central (Q2)
- **Deep Peripheral (V3)**: Slow equilibration with central (Q3)

The concentration-time profile shows tri-exponential decay:
- **Alpha phase**: Rapid initial decline
- **Beta phase**: Intermediate decline
- **Gamma phase**: Terminal elimination

---

## Basic Example

```python
import neopkpd

result = neopkpd.simulate_pk_threecomp_iv_bolus(
    cl=5.0,       # Clearance (L/h)
    v1=10.0,      # Central volume (L)
    q2=20.0,      # Rapid distribution clearance (L/h)
    v2=30.0,      # Shallow peripheral (L)
    q3=5.0,       # Slow distribution clearance (L/h)
    v3=100.0,     # Deep peripheral (L)
    doses=[{"time": 0.0, "amount": 1000.0}],
    t0=0.0,
    t1=168.0,     # 7 days
    saveat=[i * 0.5 for i in range(337)]
)

# Three phases visible in concentration decline
conc = result['observations']['conc']
t = result['t']

print(f"C at 0h: {conc[0]:.2f} mg/L (initial)")
print(f"C at 1h: {conc[2]:.2f} mg/L (alpha phase)")
print(f"C at 12h: {conc[24]:.2f} mg/L (beta phase)")
print(f"C at 168h: {conc[-1]:.4f} mg/L (gamma phase)")
```

---

## Clinical Example: Propofol

```python
import neopkpd

# Propofol typical parameters (per-minute converted to per-hour)
result = neopkpd.simulate_pk_threecomp_iv_bolus(
    cl=1.6 * 60,      # 96 L/h
    v1=4.3,           # L (central - blood)
    q2=2.3 * 60,      # 138 L/h (rapid - muscle)
    v2=22.0,          # L
    q3=0.8 * 60,      # 48 L/h (slow - fat)
    v3=200.0,         # L
    doses=[{"time": 0.0, "amount": 200.0}],  # 200 mg bolus
    t0=0.0,
    t1=6.0,           # 6 hours
    saveat=[i * 0.05 for i in range(121)]  # 3-min resolution
)

# Rapid redistribution from brain (central) to muscle
conc = result['observations']['conc']
t = result['t']

# Find time to 50% of initial concentration
c0 = conc[0]
idx_50 = next(i for i, c in enumerate(conc) if c < c0 * 0.5)
print(f"Time to 50% of C0: {t[idx_50]:.2f} h ({t[idx_50]*60:.1f} min)")
```

---

## Drug Distribution Over Time

```python
import neopkpd

result = neopkpd.simulate_pk_threecomp_iv_bolus(
    cl=5.0, v1=10.0, q2=20.0, v2=30.0, q3=5.0, v3=100.0,
    doses=[{"time": 0.0, "amount": 1000.0}],
    t0=0.0, t1=168.0,
    saveat=[0, 1, 4, 12, 24, 48, 96, 168]
)

a_central = result['states']['A_central']
a_shallow = result['states']['A_periph1']
a_deep = result['states']['A_periph2']
t = result['t']

print("Time (h) | Central | Shallow | Deep   | Total")
print("-" * 55)
for i, time in enumerate(t):
    total = a_central[i] + a_shallow[i] + a_deep[i]
    print(f"{time:8.0f} | {a_central[i]:7.1f} | {a_shallow[i]:7.1f} | {a_deep[i]:6.1f} | {total:.1f}")
```

**Expected Pattern:**
- t=0: All drug in central
- t=1-4h: Redistribution to shallow peripheral
- t=12-48h: Accumulation in deep peripheral
- t>48h: Deep peripheral contains most drug

---

## Half-Lives

```python
import numpy as np

# Parameters
cl, v1, q2, v2, q3, v3 = 5.0, 10.0, 20.0, 30.0, 5.0, 100.0

# Micro-rate constants
k10 = cl / v1      # Elimination
k12 = q2 / v1      # To shallow
k21 = q2 / v2      # From shallow
k13 = q3 / v1      # To deep
k31 = q3 / v3      # From deep

print(f"k10 (elimination): {k10:.3f} /h")
print(f"k12 (to shallow): {k12:.3f} /h")
print(f"k21 (from shallow): {k21:.3f} /h")
print(f"k13 (to deep): {k13:.3f} /h")
print(f"k31 (from deep): {k31:.3f} /h")

# Volume at steady state
vss = v1 + v2 + v3
mrt = vss / cl
print(f"\nVss: {vss} L")
print(f"MRT: {mrt:.1f} h")
```

---

## IV Infusion

```python
# Continuous infusion to steady state
result = neopkpd.simulate_pk_threecomp_iv_bolus(
    cl=5.0, v1=10.0, q2=20.0, v2=30.0, q3=5.0, v3=100.0,
    doses=[{"time": 0.0, "amount": 480.0, "duration": 24.0}],  # 20 mg/h
    t0=0.0, t1=72.0,
    saveat=[i * 0.5 for i in range(145)]
)

# Steady-state concentration
css_theoretical = (480.0 / 24.0) / 5.0  # Rate / CL
print(f"Theoretical Css: {css_theoretical:.1f} mg/L")
print(f"Simulated C at 24h: {result['observations']['conc'][48]:.2f} mg/L")
```

---

## Multiple Dosing

```python
# 1000 mg every 24 hours for 7 days
doses = [{"time": i * 24.0, "amount": 1000.0} for i in range(7)]

result = neopkpd.simulate_pk_threecomp_iv_bolus(
    cl=5.0, v1=10.0, q2=20.0, v2=30.0, q3=5.0, v3=100.0,
    doses=doses,
    t0=0.0, t1=168.0,
    saveat=[i * 1.0 for i in range(169)]
)

# With deep peripheral, accumulation continues longer
conc = result['observations']['conc']
t = result['t']

# Compare first and last trough
first_trough = conc[24]   # 24h
last_trough = conc[168]   # 168h

print(f"First trough (24h): {first_trough:.2f} mg/L")
print(f"Last trough (168h): {last_trough:.2f} mg/L")
print(f"Accumulation ratio: {last_trough/first_trough:.2f}")
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pk_profile

result = neopkpd.simulate_pk_threecomp_iv_bolus(
    cl=5.0, v1=10.0, q2=20.0, v2=30.0, q3=5.0, v3=100.0,
    doses=[{"time": 0.0, "amount": 1000.0}],
    t0=0.0, t1=168.0,
    saveat=[i * 0.5 for i in range(337)]
)

# Semi-log plot shows three phases
fig = plot_pk_profile(
    result['t'],
    result['observations']['conc'],
    title="Three-Compartment IV Bolus",
    xlabel="Time (h)",
    ylabel="Concentration (mg/L)",
    log_y=True
)
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| k10 | CL / V1 |
| k12, k21 | Q2/V1, Q2/V2 |
| k13, k31 | Q3/V1, Q3/V3 |
| C(t) | $Ae^{-\alpha t} + Be^{-\beta t} + Ce^{-\gamma t}$ |
| Vss | V1 + V2 + V3 |
| MRT | Vss / CL |
| Terminal t1/2 | 0.693 / gamma |

---

## See Also

- [Two-Compartment IV](twocomp-iv.md) - Simpler model
- [One-Compartment IV](iv-bolus.md) - Simplest model
- [Population Simulation](../../population/index.md)
