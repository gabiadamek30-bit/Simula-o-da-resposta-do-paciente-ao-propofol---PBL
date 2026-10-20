# Indirect Response Models (IRM)

Mechanism-based PD models where drug affects the production (Kin) or elimination (Kout) of a response variable.

---

## Available Functions

| Function | Model Type | Mechanism |
|----------|------------|-----------|
| `simulate_pkpd_indirect_response` | IRM-III | Inhibition of Kout |
| `simulate_pkpd_irm1` | IRM-I | Inhibition of Kin |
| `simulate_pkpd_irm2` | IRM-II | Stimulation of Kin |
| `simulate_pkpd_irm4` | IRM-IV | Stimulation of Kout |

---

## The Four IRM Types

| Model | Target | Direction | Effect on R |
|-------|--------|-----------|-------------|
| IRM-I | Kin (production) | Inhibition | Decreases R |
| IRM-II | Kin (production) | Stimulation | Increases R |
| IRM-III | Kout (elimination) | Inhibition | Increases R |
| IRM-IV | Kout (elimination) | Stimulation | Decreases R |

---

## IRM-I: Inhibition of Production

### Function Signature

```python
neopkpd.simulate_pkpd_irm1(
    cl: float, v: float, doses: list[dict],
    kin: float, kout: float, r0: float,
    imax: float, ic50: float,
    t0: float, t1: float, saveat: list[float],
    pk_kind: str = "OneCompIVBolus",
    ka: float | None = None,
    q: float | None = None,
    v2: float | None = None,
    ...
) -> dict
```

### Equation

$$\frac{dR}{dt} = K_{in} \cdot (1 - I(C)) - K_{out} \cdot R$$

Where: $I(C) = \frac{I_{max} \cdot C}{IC_{50} + C}$

### Example: Corticosteroid Effect on Cortisol

```python
import neopkpd

# Cortisol dynamics: baseline 15 mcg/dL, t1/2 ~1.5 h
kout = 0.693 / 1.5  # 0.46/h
r0 = 15.0
kin = kout * r0

result = neopkpd.simulate_pkpd_irm1(
    cl=2.0, v=50.0,
    doses=[{"time": 0.0, "amount": 10.0}],
    kin=kin, kout=kout, r0=r0,
    imax=0.9, ic50=0.01,  # Very potent
    t0=0.0, t1=48.0,
    saveat=[i * 0.5 for i in range(97)]
)

response = result['observations']['response']
print(f"Baseline cortisol: {r0} mcg/dL")
print(f"Minimum cortisol: {min(response):.1f} mcg/dL")
```

---

## IRM-II: Stimulation of Production

### Equation

$$\frac{dR}{dt} = K_{in} \cdot (1 + S(C)) - K_{out} \cdot R$$

Where: $S(C) = \frac{S_{max} \cdot C}{SC_{50} + C}$

### Example: EPO Effect on Red Blood Cells

```python
import neopkpd

# RBC dynamics: long turnover
kout = 0.693 / (30 * 24)  # 30-day half-life
r0 = 5.0  # million/mcL
kin = kout * r0

result = neopkpd.simulate_pkpd_irm2(
    cl=0.5, v=10.0,
    doses=[{"time": i * 168, "amount": 100.0} for i in range(4)],  # Weekly
    kin=kin, kout=kout, r0=r0,
    smax=3.0, sc50=0.1,
    t0=0.0, t1=672.0,  # 4 weeks
    saveat=[i * 12 for i in range(57)]
)

response = result['observations']['response']
print(f"Baseline RBC: {r0} million/mcL")
print(f"Maximum RBC: {max(response):.2f} million/mcL")
```

---

## IRM-III: Inhibition of Elimination

### Function (simpler signature)

```python
neopkpd.simulate_pkpd_indirect_response(
    cl: float, v: float, doses: list[dict],
    kin: float, kout: float, r0: float,
    imax: float, ic50: float,
    t0: float, t1: float, saveat: list[float],
    pk_kind: str = "OneCompIVBolus",
    ka: float | None = None,
    ...
) -> dict
```

### Equation

$$\frac{dR}{dt} = K_{in} - K_{out} \cdot (1 - I(C)) \cdot R$$

### Example: Warfarin Effect

```python
import neopkpd

# Clotting factor dynamics
kout = 0.693 / 36.0  # 36-hour half-life
r0 = 100.0  # % of normal
kin = kout * r0

result = neopkpd.simulate_pkpd_indirect_response(
    cl=0.1, v=10.0,
    doses=[{"time": 0.0, "amount": 5.0}],
    kin=kin, kout=kout, r0=r0,
    imax=0.95, ic50=1.5,
    t0=0.0, t1=120.0,
    saveat=[i * 1.0 for i in range(121)]
)

response = result['observations']['response']
print(f"Baseline: {r0}%")
print(f"Maximum response: {max(response):.1f}%")
print(f"Time to max: {response.index(max(response))} h")
```

---

## IRM-IV: Stimulation of Elimination

### Equation

$$\frac{dR}{dt} = K_{in} - K_{out} \cdot (1 + S(C)) \cdot R$$

### Example: Laxative Effect

```python
import neopkpd

kout = 0.693 / 12.0  # 12-hour transit
r0 = 100.0
kin = kout * r0

result = neopkpd.simulate_pkpd_irm4(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 200.0}],
    kin=kin, kout=kout, r0=r0,
    smax=5.0, sc50=0.5,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

response = result['observations']['response']
print(f"Minimum response: {min(response):.1f}")
print(f"Theoretical minimum: {r0 / (1 + 5.0):.1f}")
```

---

## Comparing IRM Types

```python
import neopkpd

# Same baseline and turnover
kout = 0.1  # 1/h
r0 = 100.0
kin = kout * r0

# Same PK
doses = [{"time": 0.0, "amount": 100.0}]

# IRM-I: Inhibit production → Decrease
result_irm1 = neopkpd.simulate_pkpd_irm1(
    cl=2.0, v=20.0, doses=doses,
    kin=kin, kout=kout, r0=r0,
    imax=0.8, ic50=2.0,
    t0=0.0, t1=48.0, saveat=[i for i in range(49)]
)

# IRM-II: Stimulate production → Increase
result_irm2 = neopkpd.simulate_pkpd_irm2(
    cl=2.0, v=20.0, doses=doses,
    kin=kin, kout=kout, r0=r0,
    smax=3.0, sc50=2.0,
    t0=0.0, t1=48.0, saveat=[i for i in range(49)]
)

# IRM-IV: Stimulate elimination → Decrease
result_irm4 = neopkpd.simulate_pkpd_irm4(
    cl=2.0, v=20.0, doses=doses,
    kin=kin, kout=kout, r0=r0,
    smax=4.0, sc50=2.0,
    t0=0.0, t1=48.0, saveat=[i for i in range(49)]
)

print("Model    | Min/Max Response")
print("-" * 30)
print(f"IRM-I    | {min(result_irm1['observations']['response']):.1f} (decrease)")
print(f"IRM-II   | {max(result_irm2['observations']['response']):.1f} (increase)")
print(f"IRM-IV   | {min(result_irm4['observations']['response']):.1f} (decrease)")
```

---

## Steady State Predictions

```python
# At high drug concentration (complete effect):

# IRM-I: R_min = R0 × (1 - Imax)
r0, imax = 100.0, 0.8
print(f"IRM-I steady state: {r0 * (1 - imax):.1f}")

# IRM-II: R_max = R0 × (1 + Smax)
r0, smax = 100.0, 3.0
print(f"IRM-II steady state: {r0 * (1 + smax):.1f}")

# IRM-III: R_max = R0 / (1 - Imax)  [limited by Imax < 1]
r0, imax = 100.0, 0.8
print(f"IRM-III steady state: {r0 / (1 - imax):.1f}")

# IRM-IV: R_min = R0 / (1 + Smax)
r0, smax = 100.0, 4.0
print(f"IRM-IV steady state: {r0 / (1 + smax):.1f}")
```

---

## Model Selection Guide

| Observation | Suggested Model |
|-------------|-----------------|
| Drug decreases biomarker | IRM-I or IRM-IV |
| Drug increases biomarker | IRM-II or IRM-III |
| Effect persists after drug washout | Any IRM |
| Known to inhibit synthesis | IRM-I |
| Known to stimulate synthesis | IRM-II |
| Known to inhibit degradation | IRM-III |
| Known to stimulate elimination | IRM-IV |

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pkpd_profile

result = neopkpd.simulate_pkpd_irm1(
    cl=2.0, v=20.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    kin=10.0, kout=0.1, r0=100.0,
    imax=0.8, ic50=2.0,
    t0=0.0, t1=48.0,
    saveat=[i * 0.5 for i in range(97)]
)

fig = plot_pkpd_profile(
    result['t'],
    result['observations']['conc'],
    result['observations']['response'],
    title="IRM-I: Inhibition of Production",
    effect_label="Response"
)
```

---

## Equations Summary

| Model | dR/dt | Steady State |
|-------|-------|--------------|
| IRM-I | $K_{in}(1-I(C)) - K_{out}R$ | $R_0(1-I_{max})$ |
| IRM-II | $K_{in}(1+S(C)) - K_{out}R$ | $R_0(1+S_{max})$ |
| IRM-III | $K_{in} - K_{out}(1-I(C))R$ | $R_0/(1-I_{max})$ |
| IRM-IV | $K_{in} - K_{out}(1+S(C))R$ | $R_0/(1+S_{max})$ |

Common: $R_0 = K_{in}/K_{out}$, Recovery $t_{1/2} = \ln(2)/K_{out}$

---

## See Also

- [Direct Emax](direct-emax.md) - Simpler model
- [Sigmoid Emax](sigmoid-emax.md) - Variable steepness
- [Effect Compartment](effect-compartment.md) - Temporal delay
