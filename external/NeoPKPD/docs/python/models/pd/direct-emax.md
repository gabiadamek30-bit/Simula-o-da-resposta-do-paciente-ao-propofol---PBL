# Direct Emax Model

Simple hyperbolic concentration-effect relationship where effect is directly proportional to receptor occupancy.

---

## Function Signature

```python
neopkpd.simulate_pkpd_direct_emax(
    cl: float,
    v: float,
    doses: list[dict],
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
| `e0` | float | Baseline effect |
| `emax` | float | Maximum effect change |
| `ec50` | float | Concentration at 50% Emax (mg/L) |
| `pk_kind` | str | PK model type |
| `ka` | float | Absorption rate (for oral models) |

### Supported PK Models

- `"OneCompIVBolus"` - One-compartment IV
- `"OneCompOralFirstOrder"` - One-compartment oral (requires `ka`)

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {...},
    "observations": {
        "conc": [...],       # Plasma concentration
        "effect": [...]      # Pharmacodynamic effect
    },
    "metadata": {...}
}
```

---

## Model Equation

$$E(C) = E_0 + \frac{E_{max} \cdot C}{EC_{50} + C}$$

---

## Basic Examples

### IV Bolus with Direct Effect

```python
import neopkpd

result = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    e0=80.0,      # Baseline (e.g., heart rate)
    emax=-30.0,   # Maximum reduction
    ec50=2.0,     # mg/L
    t0=0.0,
    t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Effect tracks concentration directly
conc = result['observations']['conc']
effect = result['observations']['effect']
t = result['t']

print(f"Baseline effect (E0): {effect[0]:.1f}")  # Effect at C=Cmax
print(f"Effect at 24h: {effect[-1]:.1f}")         # Near baseline
```

### Oral Administration

```python
result = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    e0=0.0,       # Baseline
    emax=100.0,   # Maximum stimulation
    ec50=5.0,
    t0=0.0,
    t1=24.0,
    saveat=[i * 0.25 for i in range(97)],
    pk_kind="OneCompOralFirstOrder",
    ka=1.5
)

# Find maximum effect and when it occurs
effect = result['observations']['effect']
t = result['t']
max_effect_idx = max(range(len(effect)), key=lambda i: effect[i])

print(f"Maximum effect: {effect[max_effect_idx]:.1f}")
print(f"Time of max effect: {t[max_effect_idx]:.2f} h")
```

---

## Multiple Dosing

```python
import neopkpd

# 500 mg every 8 hours
doses = [{"time": i * 8.0, "amount": 500.0} for i in range(6)]

result = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    doses=doses,
    e0=100.0,      # Baseline blood pressure
    emax=-40.0,    # Maximum reduction
    ec50=3.0,
    t0=0.0, t1=48.0,
    saveat=[i * 0.5 for i in range(97)]
)

effect = result['observations']['effect']

# Effect oscillates with concentration
print(f"Max effect: {min(effect):.1f}")   # Note: negative Emax means min is max effect
print(f"Min effect: {max(effect):.1f}")   # Trough effect
```

---

## Inhibitory vs Stimulatory Effects

### Inhibitory (Emax < 0)

```python
# Drug reduces blood pressure
result = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    e0=140.0,      # Baseline BP
    emax=-40.0,    # Maximum reduction
    ec50=2.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Effect range: 140 → 100 mmHg
```

### Stimulatory (Emax > 0)

```python
# Drug increases enzyme activity
result = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    e0=100.0,      # Baseline activity
    emax=200.0,    # Maximum increase
    ec50=5.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Effect range: 100 → 300 units
```

---

## Potency vs Efficacy

```python
import neopkpd

# Drug A: High potency, moderate efficacy
result_a = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    e0=0.0, emax=80.0, ec50=0.5,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Drug B: Low potency, high efficacy
result_b = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    e0=0.0, emax=100.0, ec50=5.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

# Compare effects
print("Drug A (high potency): EC50 = 0.5, Emax = 80")
print("Drug B (low potency): EC50 = 5.0, Emax = 100")
print(f"Drug A max effect: {max(result_a['observations']['effect']):.1f}")
print(f"Drug B max effect: {max(result_b['observations']['effect']):.1f}")
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pkpd_profile

result = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    e0=80.0, emax=-30.0, ec50=2.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

fig = plot_pkpd_profile(
    result['t'],
    result['observations']['conc'],
    result['observations']['effect'],
    title="Direct Emax PKPD",
    conc_label="Concentration (mg/L)",
    effect_label="Heart Rate (bpm)"
)
```

---

## Concentration-Effect Relationship

```python
import neopkpd
import numpy as np

result = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    e0=0.0, emax=100.0, ec50=5.0,
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

conc = np.array(result['observations']['conc'])
effect = np.array(result['observations']['effect'])

# Verify Emax equation
# At C = EC50, Effect should be E0 + Emax/2 = 50
idx_ec50 = np.argmin(np.abs(conc - 5.0))
print(f"Effect at C=EC50: {effect[idx_ec50]:.1f} (expected: 50.0)")
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Effect | $E_0 + E_{max} \cdot C / (EC_{50} + C)$ |
| Fraction of Emax | $C / (EC_{50} + C)$ |
| C for target effect | $EC_{50} \cdot (E - E_0) / (E_{max} - (E - E_0))$ |
| Sensitivity | $E_{max} / EC_{50}$ |

---

## See Also

- [Sigmoid Emax](sigmoid-emax.md) - Variable steepness
- [Effect Compartment](effect-compartment.md) - With temporal delay
- [Indirect Response](indirect-response.md) - Mechanism-based
- [PKPD Visualization](../../viz/pkpd.md)
