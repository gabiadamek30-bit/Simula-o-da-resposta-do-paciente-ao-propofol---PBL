# Michaelis-Menten Elimination

One-compartment model with saturable (capacity-limited) elimination kinetics.

---

## Function Signature

```python
neopkpd.simulate_pk_michaelis_menten(
    vmax: float,
    km: float,
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
| `vmax` | float | Maximum elimination rate (mg/h) |
| `km` | float | Michaelis constant (mg/L) |
| `v` | float | Volume of distribution (L) |
| `doses` | list[dict] | Dose events |

---

## Model Equation

$$\frac{dA}{dt} = -\frac{V_{max} \cdot C}{K_m + C}$$

$$C = \frac{A}{V}$$

### Kinetic Behavior

| Condition | Behavior | Apparent CL |
|-----------|----------|-------------|
| C << Km | First-order (linear) | CL ≈ Vmax/Km |
| C ≈ Km | Mixed-order | Variable |
| C >> Km | Zero-order (saturated) | Rate ≈ Vmax |

---

## Basic Example

```python
import neopkpd

result = neopkpd.simulate_pk_michaelis_menten(
    vmax=500.0,   # Maximum rate (mg/h)
    km=10.0,      # Km (mg/L)
    v=50.0,       # Volume (L)
    doses=[{"time": 0.0, "amount": 1000.0}],
    t0=0.0,
    t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

conc = result['observations']['conc']
t = result['t']

# Initial concentration
c0 = conc[0]  # 1000/50 = 20 mg/L

print(f"C0: {c0:.1f} mg/L (C >> Km: saturated)")
print(f"C at 6h: {conc[24]:.2f} mg/L")
print(f"C at 24h: {conc[96]:.3f} mg/L (C < Km: linear)")
```

---

## Dose-Dependent Half-Life

With Michaelis-Menten kinetics, half-life increases with dose:

```python
import neopkpd

doses_list = [100.0, 500.0, 1000.0, 2000.0]

print("Dose (mg) | C0 (mg/L) | t1/2 (h)")
print("-" * 40)

for dose in doses_list:
    result = neopkpd.simulate_pk_michaelis_menten(
        vmax=500.0, km=10.0, v=50.0,
        doses=[{"time": 0.0, "amount": dose}],
        t0=0.0, t1=96.0,
        saveat=[i * 0.5 for i in range(193)]
    )

    conc = result['observations']['conc']
    t = result['t']

    # Find t1/2 (time to reach half of initial)
    c0 = conc[0]
    target = c0 / 2

    t_half = None
    for i, c in enumerate(conc):
        if c <= target:
            t_half = t[i]
            break

    if t_half:
        print(f"{dose:9.0f} | {c0:9.1f} | {t_half:7.1f}")
    else:
        print(f"{dose:9.0f} | {c0:9.1f} | > 96")
```

**Expected Pattern:**
- Higher doses → Longer t1/2 (saturation)
- Low doses → Shorter t1/2 (linear kinetics)

---

## Dose Proportionality

AUC increases more than proportionally with dose:

```python
import neopkpd
import numpy as np

doses_list = [100.0, 200.0, 500.0, 1000.0]
auc_values = []

print("Dose (mg) | AUC (mg*h/L)")
print("-" * 30)

for dose in doses_list:
    result = neopkpd.simulate_pk_michaelis_menten(
        vmax=500.0, km=10.0, v=50.0,
        doses=[{"time": 0.0, "amount": dose}],
        t0=0.0, t1=120.0,
        saveat=[i * 0.1 for i in range(1201)]
    )

    conc = np.array(result['observations']['conc'])
    t = np.array(result['t'])

    auc = np.trapz(conc, t)
    auc_values.append(auc)

    print(f"{dose:9.0f} | {auc:12.1f}")

# Check proportionality
print(f"\nDose ratio (10x): {doses_list[-1] / doses_list[0]:.0f}")
print(f"AUC ratio: {auc_values[-1] / auc_values[0]:.2f}")
# AUC ratio > 10 indicates saturation
```

---

## Clinical Example: Phenytoin

```python
import neopkpd

# Phenytoin typical parameters (for 70 kg patient)
# Vmax ≈ 7 mg/kg/day, Km ≈ 4-6 mg/L, V ≈ 0.65 L/kg
vmax = 7.0 * 70 / 24  # mg/h
km = 5.0              # mg/L
v = 0.65 * 70         # L

# Loading dose
loading = 15 * 70  # 15 mg/kg = 1050 mg

result = neopkpd.simulate_pk_michaelis_menten(
    vmax=vmax, km=km, v=v,
    doses=[{"time": 0.0, "amount": loading}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.5 for i in range(49)]
)

conc = result['observations']['conc']

# Therapeutic range: 10-20 mg/L
print(f"Post-loading C: {conc[0]:.1f} mg/L")
print(f"24h C: {conc[-1]:.1f} mg/L")
print("Therapeutic range: 10-20 mg/L")
```

---

## Steady-State Considerations

At steady state with constant infusion rate R:

$$C_{ss} = \frac{R \cdot K_m}{V_{max} - R}$$

This only applies when R < Vmax.

```python
vmax = 500.0  # mg/h
km = 10.0     # mg/L

print(f"Maximum sustainable infusion rate: {vmax} mg/h")

# At different infusion rates
for fraction in [0.5, 0.8, 0.9, 0.95]:
    r = fraction * vmax
    css = (r * km) / (vmax - r)
    print(f"At R = {r:.0f} mg/h ({fraction*100:.0f}% of Vmax): Css = {css:.1f} mg/L")

# If R >= Vmax, no steady state is achievable!
print("\nWARNING: If R >= Vmax, concentrations increase indefinitely!")
```

---

## Comparison with Linear Model

```python
import neopkpd

# Michaelis-Menten at low dose (approximately linear)
result_mm = neopkpd.simulate_pk_michaelis_menten(
    vmax=500.0, km=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],  # Low dose: C0 = 2 mg/L << Km
    t0=0.0, t1=24.0,
    saveat=[i * 0.1 for i in range(241)]
)

# Equivalent linear model: CL = Vmax/Km = 50 L/h
result_linear = neopkpd.simulate_pk_iv_bolus(
    cl=50.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.1 for i in range(241)]
)

# At low dose, profiles should be nearly identical
mm_conc = result_mm['observations']['conc']
lin_conc = result_linear['observations']['conc']

print("Low dose (C << Km): MM ≈ Linear")
print(f"MM at 4h: {mm_conc[40]:.3f} mg/L")
print(f"Linear at 4h: {lin_conc[40]:.3f} mg/L")
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pk_profile

# Compare different doses
doses = [100.0, 500.0, 1000.0]

for dose in doses:
    result = neopkpd.simulate_pk_michaelis_menten(
        vmax=500.0, km=10.0, v=50.0,
        doses=[{"time": 0.0, "amount": dose}],
        t0=0.0, t1=48.0,
        saveat=[i * 0.25 for i in range(193)]
    )

    # Plot each dose level
    # fig = plot_pk_profile(result['t'], result['observations']['conc'], ...)
```

---

## Apparent Clearance

$$CL_{app} = \frac{V_{max}}{K_m + C}$$

At low concentrations: $CL_{max} = V_{max}/K_m$

```python
vmax = 500.0  # mg/h
km = 10.0     # mg/L

print("C (mg/L) | CL_app (L/h)")
print("-" * 25)

for c in [0.1, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0]:
    cl_app = vmax / (km + c)
    print(f"{c:8.1f} | {cl_app:11.1f}")

print(f"\nCL_max (C→0): {vmax/km:.1f} L/h")
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Elimination rate | $-V_{max} \cdot C / (K_m + C)$ |
| CL_apparent | $V_{max} / (K_m + C)$ |
| CL_max (C→0) | $V_{max} / K_m$ |
| t1/2 (apparent) | $0.693 \cdot V \cdot (K_m + C) / V_{max}$ |
| Css (infusion R) | $R \cdot K_m / (V_{max} - R)$ |

---

## See Also

- [One-Compartment IV](iv-bolus.md) - Linear model
- [Transit Absorption](transit.md) - Complex absorption
- [Population Simulation](../../population/index.md) - Adding variability
