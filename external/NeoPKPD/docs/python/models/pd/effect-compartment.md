# Effect Compartment (Biophase) Model

Hypothetical effect site compartment to model temporal delays between plasma concentration and pharmacodynamic effect.

---

## Function Signature

```python
neopkpd.simulate_pkpd_biophase_equilibration(
    cl: float,
    v: float,
    doses: list[dict],
    ke0: float,
    e0: float,
    emax: float,
    ec50: float,
    t0: float,
    t1: float,
    saveat: list[float],
    pk_kind: str = "OneCompIVBolus",
    ka: float | None = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> dict
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `cl` | float | Clearance (L/h) |
| `v` | float | Volume of distribution (L) |
| `doses` | list[dict] | Dose events |
| `ke0` | float | Effect site equilibration rate (1/h) |
| `e0` | float | Baseline effect |
| `emax` | float | Maximum effect change |
| `ec50` | float | Effect site concentration at 50% Emax |
| `pk_kind` | str | PK model type |
| `ka` | float | Absorption rate (for oral) |

### Key Derived Parameter

**Equilibration half-life**: $t_{1/2,ke0} = \ln(2) / ke0$

---

## Model Equations

Effect site equilibration:
$$\frac{dC_e}{dt} = k_{e0} \cdot (C_p - C_e)$$

Effect from effect site concentration:
$$E(C_e) = E_0 + \frac{E_{max} \cdot C_e}{EC_{50} + C_e}$$

---

## Basic Example

```python
import neopkpd

result = neopkpd.simulate_pkpd_biophase_equilibration(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    ke0=0.5,       # Equilibration rate (1/h)
    e0=0.0,
    emax=100.0,
    ec50=5.0,
    t0=0.0,
    t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

conc = result['observations']['conc']
effect = result['observations']['effect']
t = result['t']

# Effect lags behind concentration
conc_max_idx = max(range(len(conc)), key=lambda i: conc[i])
effect_max_idx = max(range(len(effect)), key=lambda i: effect[i])

print(f"Cmax at t = {t[conc_max_idx]:.1f} h")
print(f"Emax at t = {t[effect_max_idx]:.1f} h")
print(f"Effect delay: {t[effect_max_idx] - t[conc_max_idx]:.1f} h")
```

---

## Effect of ke0 on Response

```python
import neopkpd

ke0_values = [0.2, 0.5, 1.0, 2.0, 5.0]

print("ke0 (1/h) | t1/2,ke0 (h) | Tmax effect (h)")
print("-" * 50)

for ke0 in ke0_values:
    result = neopkpd.simulate_pkpd_biophase_equilibration(
        cl=5.0, v=50.0,
        doses=[{"time": 0.0, "amount": 500.0}],
        ke0=ke0, e0=0.0, emax=100.0, ec50=5.0,
        t0=0.0, t1=24.0,
        saveat=[i * 0.1 for i in range(241)]
    )

    effect = result['observations']['effect']
    t = result['t']
    tmax_effect = t[max(range(len(effect)), key=lambda i: effect[i])]
    t_half_ke0 = 0.693 / ke0

    print(f"{ke0:9.1f} | {t_half_ke0:12.2f} | {tmax_effect:15.2f}")
```

**Expected Pattern:**
- Higher ke0 → Faster equilibration → Earlier Tmax,effect
- Lower ke0 → Slower equilibration → More delayed Tmax,effect

---

## Counter-Clockwise Hysteresis

```python
import neopkpd
import numpy as np

result = neopkpd.simulate_pkpd_biophase_equilibration(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    ke0=0.3, e0=0.0, emax=100.0, ec50=5.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.1 for i in range(241)]
)

conc = result['observations']['conc']
effect = result['observations']['effect']

# Plot Effect vs Concentration would show counter-clockwise loop
# - Rising limb: Effect lags behind rising concentration
# - Falling limb: Effect persists as concentration falls
print("Hysteresis visible when plotting Effect vs Concentration")
print("Counter-clockwise loop indicates effect site delay")
```

---

## Clinical Example: Propofol Anesthesia

```python
import neopkpd

# Propofol effect compartment
result = neopkpd.simulate_pkpd_biophase_equilibration(
    cl=100.0,      # Fast clearance (L/h)
    v=20.0,        # Small central volume (L)
    doses=[{"time": 0.0, "amount": 200.0}],
    ke0=2.0,       # Fast equilibration (~20 sec t1/2)
    e0=0.0,        # Awake = 0
    emax=100.0,    # Deep anesthesia = 100
    ec50=3.0,      # mcg/mL
    t0=0.0,
    t1=0.5,        # 30 minutes
    saveat=[i * 0.01 for i in range(51)]
)

effect = result['observations']['effect']
t = result['t']

# Time to loss of consciousness (effect > 50)
loc_idx = next((i for i, e in enumerate(effect) if e > 50), None)
if loc_idx:
    print(f"Time to LOC: {t[loc_idx] * 60:.1f} seconds")
```

---

## Oral Administration with Effect Delay

```python
import neopkpd

result = neopkpd.simulate_pkpd_biophase_equilibration(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    ke0=0.5, e0=0.0, emax=100.0, ec50=5.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)],
    pk_kind="OneCompOralFirstOrder",
    ka=1.5
)

conc = result['observations']['conc']
effect = result['observations']['effect']
t = result['t']

# Two sources of delay: absorption + effect compartment
cmax_t = t[max(range(len(conc)), key=lambda i: conc[i])]
emax_t = t[max(range(len(effect)), key=lambda i: effect[i])]

print(f"Tmax (concentration): {cmax_t:.2f} h")
print(f"Tmax (effect): {emax_t:.2f} h")
```

---

## Comparing Direct vs Effect Compartment

```python
import neopkpd

# Direct Emax (no delay)
result_direct = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    e0=0.0, emax=100.0, ec50=5.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Effect compartment (with delay)
result_biophase = neopkpd.simulate_pkpd_biophase_equilibration(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    ke0=0.3, e0=0.0, emax=100.0, ec50=5.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Compare Tmax
t = result_direct['t']
tmax_direct = t[max(range(len(result_direct['observations']['effect'])),
                    key=lambda i: result_direct['observations']['effect'][i])]
tmax_biophase = t[max(range(len(result_biophase['observations']['effect'])),
                      key=lambda i: result_biophase['observations']['effect'][i])]

print(f"Direct Emax Tmax: {tmax_direct:.2f} h")
print(f"Effect Compartment Tmax: {tmax_biophase:.2f} h")
print(f"Delay: {tmax_biophase - tmax_direct:.2f} h")
```

---

## Model Selection Guide

```python
# Rule of thumb:
# If t1/2,ke0 < elimination t1/2 / 10, use Direct Emax
# Otherwise, use Effect Compartment

def recommend_model(ke0, cl, v):
    t_half_ke0 = 0.693 / ke0
    kel = cl / v
    t_half_kel = 0.693 / kel

    print(f"t1/2,ke0 = {t_half_ke0:.2f} h")
    print(f"t1/2,kel = {t_half_kel:.2f} h")

    if t_half_ke0 < t_half_kel / 10:
        print("Recommend: Direct Emax (fast equilibration)")
    else:
        print("Recommend: Effect Compartment (significant delay)")

# Examples
recommend_model(ke0=5.0, cl=5.0, v=50.0)   # Fast ke0
print()
recommend_model(ke0=0.2, cl=5.0, v=50.0)   # Slow ke0
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pkpd_profile

result = neopkpd.simulate_pkpd_biophase_equilibration(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    ke0=0.5, e0=0.0, emax=100.0, ec50=5.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

fig = plot_pkpd_profile(
    result['t'],
    result['observations']['conc'],
    result['observations']['effect'],
    title="Effect Compartment PKPD",
    effect_label="Effect (%)"
)
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| dCe/dt | $k_{e0} \cdot (C_p - C_e)$ |
| Effect | $E_0 + E_{max} \cdot C_e / (EC_{50} + C_e)$ |
| t1/2,ke0 | $\ln(2) / k_{e0}$ |
| t90% equilibration | $\ln(10) / k_{e0}$ |
| Steady state | $C_e = C_p$ |

---

## See Also

- [Direct Emax](direct-emax.md) - Without delay
- [Sigmoid Emax](sigmoid-emax.md) - Variable steepness
- [Indirect Response](indirect-response.md) - Mechanism-based
