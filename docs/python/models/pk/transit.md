# Transit Compartment Absorption

Transit compartment model for delayed and complex oral absorption, using a chain of compartments to model gastrointestinal transit.

---

## Function Signature

```python
neopkpd.simulate_pk_transit_absorption(
    n: int,
    ktr: float,
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
| `n` | int | Number of transit compartments (1-20) |
| `ktr` | float | Transit rate constant (1/h) |
| `ka` | float | Final absorption rate constant (1/h) |
| `cl` | float | Clearance (L/h) |
| `v` | float | Volume of distribution (L) |
| `doses` | list[dict] | Dose events |

### Key Derived Parameter

**Mean Transit Time (MTT)**:
$$MTT = \frac{N + 1}{K_{tr}}$$

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "Transit_1": [...],        # First transit
        "Transit_2": [...],        # Second transit
        ...
        "Transit_N": [...],        # Last transit
        "A_central": [...]         # Central compartment
    },
    "observations": {
        "conc": [...]              # Plasma concentration
    },
    "metadata": {...}
}
```

---

## Model Equations

$$\frac{dT_1}{dt} = -K_{tr} \cdot T_1$$

$$\frac{dT_i}{dt} = K_{tr} \cdot T_{i-1} - K_{tr} \cdot T_i \quad (i = 2..N)$$

$$\frac{dA}{dt} = K_a \cdot T_N - k \cdot A$$

The absorption input follows a gamma distribution profile.

---

## Basic Example

```python
import neopkpd

result = neopkpd.simulate_pk_transit_absorption(
    n=5,          # 5 transit compartments
    ktr=0.5,      # Transit rate (1/h)
    ka=2.0,       # Absorption rate (1/h)
    cl=10.0,      # Clearance (L/h)
    v=70.0,       # Volume (L)
    doses=[{"time": 0.0, "amount": 300.0}],
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

# Compare to Mean Transit Time
mtt = (5 + 1) / 0.5  # 12 hours
print(f"MTT: {mtt} h")
```

---

## Effect of Number of Transit Compartments

```python
import neopkpd

# Keep MTT constant at 12 hours
mtt = 12.0
n_values = [1, 3, 5, 10]

print("N | Ktr (1/h) | Cmax (mg/L) | Tmax (h)")
print("-" * 45)

for n in n_values:
    ktr = (n + 1) / mtt

    result = neopkpd.simulate_pk_transit_absorption(
        n=n, ktr=ktr, ka=2.0, cl=10.0, v=70.0,
        doses=[{"time": 0.0, "amount": 300.0}],
        t0=0.0, t1=48.0,
        saveat=[i * 0.1 for i in range(481)]
    )

    conc = result['observations']['conc']
    t = result['t']
    cmax_idx = max(range(len(conc)), key=lambda i: conc[i])

    print(f"{n:2d} | {ktr:9.3f} | {conc[cmax_idx]:11.2f} | {t[cmax_idx]:7.2f}")
```

**Expected Pattern:**
- Higher N → More delayed Tmax
- Higher N → Lower, broader Cmax
- Same AUC (same CL)

---

## Comparison: Transit vs Simple Oral

```python
import neopkpd
import numpy as np

# Transit absorption (5 compartments, MTT = 12h)
result_transit = neopkpd.simulate_pk_transit_absorption(
    n=5, ktr=0.5, ka=2.0, cl=10.0, v=70.0,
    doses=[{"time": 0.0, "amount": 300.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

# Simple first-order oral (slower Ka to approximate)
result_simple = neopkpd.simulate_pk_oral_first_order(
    ka=0.3, cl=10.0, v=70.0,
    doses=[{"time": 0.0, "amount": 300.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

# Transit model has:
# - Delayed onset (sigmoidal rise)
# - Broader peak
# - More physiological shape
print("Transit model: delayed onset, broader peak")
print("Simple model: immediate onset, sharper peak")
```

---

## Controlled-Release Formulation

```python
import neopkpd

# Extended-release tablet: Long MTT
result_er = neopkpd.simulate_pk_transit_absorption(
    n=8,          # More transit compartments
    ktr=0.5,      # MTT = 9/0.5 = 18 hours
    ka=0.5,       # Slow final absorption
    cl=10.0,
    v=70.0,
    doses=[{"time": 0.0, "amount": 300.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

# Immediate-release for comparison (once daily dose split into TID)
# ... would show higher peaks and valleys

conc = result_er['observations']['conc']
print(f"ER Cmax: {max(conc):.2f} mg/L")
print(f"ER Cmin (24h): {conc[96]:.2f} mg/L")
```

---

## Estimation Considerations

When fitting transit models:

```python
# N and Ktr are correlated (same MTT with different combinations)
# Typically fix N based on physiology (3-7 for oral drugs)
# Estimate Ktr (or MTT directly)

# From observed Tmax, estimate initial MTT
observed_tmax = 8.0  # hours

# MTT ≈ Tmax for transit models
# Choose N based on formulation type
n_guess = 5
ktr_guess = (n_guess + 1) / observed_tmax  # 0.75/h

print(f"Initial N: {n_guess}")
print(f"Initial Ktr: {ktr_guess:.3f} /h")
print(f"Initial MTT: {(n_guess + 1) / ktr_guess:.1f} h")
```

---

## Absorption Variability

The coefficient of variation of the absorption profile:

$$CV_{absorption} = \frac{1}{\sqrt{N + 1}}$$

```python
# More transit compartments = narrower, more reproducible absorption
n_values = [1, 3, 5, 10, 20]

print("N | CV_absorption")
print("-" * 20)
for n in n_values:
    cv = 1 / (n + 1)**0.5 * 100
    print(f"{n:2d} | {cv:5.1f}%")
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pk_profile

result = neopkpd.simulate_pk_transit_absorption(
    n=5, ktr=0.5, ka=2.0, cl=10.0, v=70.0,
    doses=[{"time": 0.0, "amount": 300.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

fig = plot_pk_profile(
    result['t'],
    result['observations']['conc'],
    title="Transit Compartment Absorption",
    xlabel="Time (h)",
    ylabel="Concentration (mg/L)"
)
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| MTT | $(N+1)/K_{tr}$ |
| t_max,input | $N/K_{tr}$ |
| CV_absorption | $1/\sqrt{N+1}$ |
| Input rate | $\frac{D \cdot K_{tr}^{N+1} \cdot t^N \cdot e^{-K_{tr}t}}{N!}$ |
| dT_1/dt | $-K_{tr} \cdot T_1$ |
| dT_i/dt | $K_{tr}(T_{i-1} - T_i)$ |
| dA/dt | $K_a T_N - (CL/V) A$ |

---

## See Also

- [One-Compartment Oral](oral.md) - Simple absorption
- [Two-Compartment Oral](twocomp-oral.md) - With distribution
- [Michaelis-Menten](michaelis-menten.md) - Nonlinear elimination
- [NCA Analysis](../../nca/index.md) - Exposure calculations
