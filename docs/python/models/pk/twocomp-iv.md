# Two-Compartment IV Bolus

Two-compartment mammillary model with central and peripheral compartments, exhibiting bi-exponential concentration decline.

---

## Function Signature

```python
neopkpd.simulate_pk_twocomp_iv_bolus(
    cl: float,
    v1: float,
    q: float,
    v2: float,
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
| `q` | float | Inter-compartmental clearance (L/h) |
| `v2` | float | Peripheral compartment volume (L) |
| `doses` | list[dict] | Dose events |
| `t0`, `t1` | float | Simulation time range |
| `saveat` | list[float] | Output time points |

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "A_central": [...],        # Amount in central
        "A_peripheral": [...]      # Amount in peripheral
    },
    "observations": {
        "conc": [...]              # Central concentration (A_central/V1)
    },
    "metadata": {...}
}
```

---

## Model Equations

$$\frac{dA_1}{dt} = -(k_{10} + k_{12}) \cdot A_1 + k_{21} \cdot A_2$$

$$\frac{dA_2}{dt} = k_{12} \cdot A_1 - k_{21} \cdot A_2$$

### Micro-rate Constants

| Constant | Formula |
|----------|---------|
| k10 | CL / V1 |
| k12 | Q / V1 |
| k21 | Q / V2 |

### Bi-exponential Solution

$$C(t) = A \cdot e^{-\alpha t} + B \cdot e^{-\beta t}$$

where $\alpha > \beta$ (distribution and terminal phases).

---

## Basic Examples

### Single IV Bolus

```python
import neopkpd

result = neopkpd.simulate_pk_twocomp_iv_bolus(
    cl=10.0,      # Clearance (L/h)
    v1=20.0,      # Central volume (L)
    q=15.0,       # Inter-compartmental clearance (L/h)
    v2=50.0,      # Peripheral volume (L)
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=48.0,
    saveat=[i * 0.5 for i in range(97)]
)

# Initial concentration (C0)
print(f"C0: {result['observations']['conc'][0]:.1f} mg/L")  # 500/20 = 25 mg/L
```

### IV Infusion

```python
# 500 mg over 2 hours
result = neopkpd.simulate_pk_twocomp_iv_bolus(
    cl=10.0, v1=20.0, q=15.0, v2=50.0,
    doses=[{"time": 0.0, "amount": 500.0, "duration": 2.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

conc = result['observations']['conc']
cmax = max(conc)
print(f"Cmax: {cmax:.2f} mg/L (at end of infusion)")
```

---

## Distribution Phases

```python
import neopkpd
import numpy as np

result = neopkpd.simulate_pk_twocomp_iv_bolus(
    cl=10.0, v1=20.0, q=15.0, v2=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.1 for i in range(241)]
)

t = np.array(result['t'])
conc = np.array(result['observations']['conc'])

# Log-linear plot shows two phases
log_conc = np.log(conc)

# Distribution phase (early - fast decline)
print(f"C at 0h: {conc[0]:.2f} mg/L")
print(f"C at 1h: {conc[10]:.2f} mg/L (distribution phase)")

# Terminal phase (late - slow decline)
print(f"C at 6h: {conc[60]:.2f} mg/L")
print(f"C at 12h: {conc[120]:.2f} mg/L (terminal phase)")
```

---

## Calculate Rate Constants

```python
import numpy as np

# Parameters
cl, v1, q, v2 = 10.0, 20.0, 15.0, 50.0

# Micro-rate constants
k10 = cl / v1      # 0.5 /h
k12 = q / v1       # 0.75 /h
k21 = q / v2       # 0.3 /h

# Hybrid constants (alpha and beta)
sum_k = k10 + k12 + k21
prod_k = k10 * k21
discriminant = np.sqrt(sum_k**2 - 4*prod_k)

alpha = (sum_k + discriminant) / 2
beta = (sum_k - discriminant) / 2

# Half-lives
t_half_alpha = 0.693 / alpha
t_half_beta = 0.693 / beta

print(f"Distribution t1/2 (alpha): {t_half_alpha:.2f} h")
print(f"Terminal t1/2 (beta): {t_half_beta:.2f} h")

# Volume at steady state
vss = v1 + v2
print(f"Vss: {vss} L")
```

---

## Multiple Dosing

```python
# 500 mg every 8 hours for 3 days
doses = [{"time": i * 8.0, "amount": 500.0} for i in range(9)]

result = neopkpd.simulate_pk_twocomp_iv_bolus(
    cl=10.0, v1=20.0, q=15.0, v2=50.0,
    doses=doses,
    t0=0.0, t1=72.0,
    saveat=[i * 0.5 for i in range(145)]
)

# Steady-state metrics in last interval (64-72 h)
conc = result['observations']['conc']
t = result['t']

ss_start = next(i for i, x in enumerate(t) if x >= 64)
ss_conc = conc[ss_start:]

print(f"Cmax,ss: {max(ss_conc):.2f} mg/L")
print(f"Cmin,ss: {min(ss_conc):.2f} mg/L")
```

---

## Compartment Amounts

```python
import neopkpd

result = neopkpd.simulate_pk_twocomp_iv_bolus(
    cl=10.0, v1=20.0, q=15.0, v2=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[0, 1, 4, 8, 12, 24]
)

# Track drug distribution over time
a_central = result['states']['A_central']
a_periph = result['states']['A_peripheral']
t = result['t']

print("Time | Central | Peripheral | % Central")
print("-" * 45)
for i, time in enumerate(t):
    total = a_central[i] + a_periph[i]
    if total > 0:
        pct_central = 100 * a_central[i] / total
        print(f"{time:4.0f}h | {a_central[i]:7.1f} | {a_periph[i]:10.1f} | {pct_central:6.1f}%")
```

---

## Clinical Example: Vancomycin

```python
import neopkpd

# Typical vancomycin parameters
result = neopkpd.simulate_pk_twocomp_iv_bolus(
    cl=4.5,       # L/h
    v1=15.0,      # L (central)
    q=4.0,        # L/h
    v2=40.0,      # L (peripheral)
    doses=[{"time": 0.0, "amount": 1000.0, "duration": 1.0}],  # 1g over 1h
    t0=0.0,
    t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

conc = result['observations']['conc']
t = result['t']

# TDM sampling times
idx_2h = next(i for i, x in enumerate(t) if x >= 2.0)
idx_trough = next(i for i, x in enumerate(t) if x >= 12.0)

print(f"Peak (2h): {conc[idx_2h]:.1f} mg/L")
print(f"Trough (12h): {conc[idx_trough]:.1f} mg/L")

# Target ranges: Peak 30-40 mg/L, Trough 15-20 mg/L
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pk_profile

result = neopkpd.simulate_pk_twocomp_iv_bolus(
    cl=10.0, v1=20.0, q=15.0, v2=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

fig = plot_pk_profile(
    result['t'],
    result['observations']['conc'],
    title="Two-Compartment IV Bolus",
    xlabel="Time (h)",
    ylabel="Concentration (mg/L)",
    log_y=True  # Semi-log plot shows bi-exponential phases
)
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| k10 | CL / V1 |
| k12, k21 | Q/V1, Q/V2 |
| alpha | $(k_{10}+k_{12}+k_{21}+\sqrt{\Delta})/2$ |
| beta | $(k_{10}+k_{12}+k_{21}-\sqrt{\Delta})/2$ |
| Vss | V1 + V2 |
| Terminal t1/2 | 0.693 / beta |
| AUC | Dose / CL |

---

## See Also

- [One-Compartment IV](iv-bolus.md) - Simpler model
- [Two-Compartment Oral](twocomp-oral.md) - With absorption
- [Three-Compartment IV](threecomp.md) - More compartments
- [Population Simulation](../../population/index.md)
