# Two-Compartment Oral

Two-compartment model with first-order oral absorption, combining absorption kinetics with distribution.

---

## Function Signature

```python
neopkpd.simulate_pk_twocomp_oral(
    ka: float,
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
| `ka` | float | Absorption rate constant (1/h) |
| `cl` | float | Apparent clearance CL/F (L/h) |
| `v1` | float | Apparent central volume V1/F (L) |
| `q` | float | Apparent inter-compartmental clearance Q/F (L/h) |
| `v2` | float | Apparent peripheral volume V2/F (L) |
| `doses` | list[dict] | Dose events |
| `alag` | float | Absorption lag time (optional) |
| `bioavailability` | float | Fraction absorbed (optional) |

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "A_gut": [...],            # Amount in gut
        "A_central": [...],        # Amount in central
        "A_peripheral": [...]      # Amount in peripheral
    },
    "observations": {
        "conc": [...]              # Plasma concentration
    },
    "metadata": {...}
}
```

---

## Model Equations

$$\frac{dA_{gut}}{dt} = -K_a \cdot A_{gut}$$

$$\frac{dA_1}{dt} = K_a \cdot A_{gut} - \frac{CL}{V_1} \cdot A_1 - \frac{Q}{V_1} \cdot A_1 + \frac{Q}{V_2} \cdot A_2$$

$$\frac{dA_2}{dt} = \frac{Q}{V_1} \cdot A_1 - \frac{Q}{V_2} \cdot A_2$$

---

## Basic Examples

### Single Oral Dose

```python
import neopkpd

result = neopkpd.simulate_pk_twocomp_oral(
    ka=1.2,       # Absorption rate (1/h)
    cl=8.0,       # Clearance (L/h)
    v1=25.0,      # Central volume (L)
    q=12.0,       # Inter-compartmental clearance (L/h)
    v2=60.0,      # Peripheral volume (L)
    doses=[{"time": 0.0, "amount": 400.0}],
    t0=0.0,
    t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

# Find Cmax and Tmax
conc = result['observations']['conc']
t = result['t']
cmax_idx = max(range(len(conc)), key=lambda i: conc[i])

print(f"Cmax: {conc[cmax_idx]:.2f} mg/L")
print(f"Tmax: {t[cmax_idx]:.2f} h")
```

### With Lag Time

```python
# 30-minute lag, 80% bioavailability
result = neopkpd.simulate_pk_twocomp_oral(
    ka=1.2, cl=8.0, v1=25.0, q=12.0, v2=60.0,
    doses=[{"time": 0.0, "amount": 400.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)],
    alag=0.5,
    bioavailability=0.8
)

print(f"Effective dose: {400 * 0.8:.0f} mg")
```

---

## Absorption Rate Effects

```python
import neopkpd

# Compare fast vs slow absorption
ka_values = [0.5, 1.0, 2.0, 4.0]

print("Ka (1/h) | Cmax (mg/L) | Tmax (h)")
print("-" * 40)

for ka in ka_values:
    result = neopkpd.simulate_pk_twocomp_oral(
        ka=ka, cl=8.0, v1=25.0, q=12.0, v2=60.0,
        doses=[{"time": 0.0, "amount": 400.0}],
        t0=0.0, t1=24.0,
        saveat=[i * 0.1 for i in range(241)]
    )

    conc = result['observations']['conc']
    t = result['t']
    cmax_idx = max(range(len(conc)), key=lambda i: conc[i])

    print(f"{ka:8.1f} | {conc[cmax_idx]:11.2f} | {t[cmax_idx]:7.2f}")
```

---

## Multiple Dosing

### BID Regimen to Steady State

```python
# 400 mg every 12 hours for 7 days
doses = [{"time": i * 12.0, "amount": 400.0} for i in range(14)]

result = neopkpd.simulate_pk_twocomp_oral(
    ka=1.2, cl=8.0, v1=25.0, q=12.0, v2=60.0,
    doses=doses,
    t0=0.0, t1=168.0,
    saveat=[i * 0.5 for i in range(337)]
)

conc = result['observations']['conc']
t = result['t']

# Steady-state in last interval (156-168 h)
ss_start = next(i for i, x in enumerate(t) if x >= 156)
ss_conc = conc[ss_start:]

print(f"Cmax,ss: {max(ss_conc):.2f} mg/L")
print(f"Cmin,ss: {min(ss_conc):.2f} mg/L")
print(f"Fluctuation: {(max(ss_conc) - min(ss_conc)) / min(ss_conc) * 100:.1f}%")
```

---

## Food Effect Study Design

```python
import neopkpd

# Fasted: Fast absorption
params_fasted = {
    "ka": 2.0,
    "cl": 8.0,
    "v1": 25.0,
    "q": 12.0,
    "v2": 60.0
}

# Fed: Slower absorption, enhanced bioavailability (25%)
params_fed = {
    "ka": 0.8,
    "cl": 6.4,  # CL/F decreases as F increases
    "v1": 25.0,
    "q": 12.0,
    "v2": 60.0
}

result_fasted = neopkpd.simulate_pk_twocomp_oral(
    **params_fasted,
    doses=[{"time": 0.0, "amount": 400.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

result_fed = neopkpd.simulate_pk_twocomp_oral(
    **params_fed,
    doses=[{"time": 0.0, "amount": 400.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

# Compare exposure
import numpy as np

auc_fasted = np.trapz(result_fasted['observations']['conc'], result_fasted['t'])
auc_fed = np.trapz(result_fed['observations']['conc'], result_fed['t'])

cmax_fasted = max(result_fasted['observations']['conc'])
cmax_fed = max(result_fed['observations']['conc'])

print(f"Fed/Fasted AUC ratio: {auc_fed/auc_fasted:.2f}")
print(f"Fed/Fasted Cmax ratio: {cmax_fed/cmax_fasted:.2f}")
```

---

## Comparison: One-Comp vs Two-Comp

```python
import neopkpd
import numpy as np

# Two-compartment oral
result_2comp = neopkpd.simulate_pk_twocomp_oral(
    ka=1.2, cl=8.0, v1=25.0, q=12.0, v2=60.0,
    doses=[{"time": 0.0, "amount": 400.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

# One-compartment oral (same total volume)
result_1comp = neopkpd.simulate_pk_oral_first_order(
    ka=1.2, cl=8.0, v=85.0,  # V = V1 + V2
    doses=[{"time": 0.0, "amount": 400.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

# Key differences:
# - 2-comp has higher initial peak (smaller V1)
# - 2-comp shows distribution phase
# - Same terminal AUC
print(f"2-comp Cmax: {max(result_2comp['observations']['conc']):.2f} mg/L")
print(f"1-comp Cmax: {max(result_1comp['observations']['conc']):.2f} mg/L")
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pk_profile

result = neopkpd.simulate_pk_twocomp_oral(
    ka=1.2, cl=8.0, v1=25.0, q=12.0, v2=60.0,
    doses=[{"time": 0.0, "amount": 400.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

fig = plot_pk_profile(
    result['t'],
    result['observations']['conc'],
    title="Two-Compartment Oral",
    xlabel="Time (h)",
    ylabel="Concentration (mg/L)"
)
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| dA_gut/dt | $-K_a \cdot A_{gut}$ |
| dA1/dt | $K_a A_{gut} - (k_{10} + k_{12})A_1 + k_{21}A_2$ |
| dA2/dt | $k_{12}A_1 - k_{21}A_2$ |
| AUC | $F \cdot Dose / CL$ |
| Absorption t1/2 | $0.693 / K_a$ |
| Terminal t1/2 | $0.693 / \beta$ |

---

## See Also

- [One-Compartment Oral](oral.md) - Simpler model
- [Two-Compartment IV](twocomp-iv.md) - Without absorption
- [Transit Absorption](transit.md) - Complex absorption
- [NCA Analysis](../../nca/index.md) - Exposure calculations
