# Sigmoid Emax (Hill) Model

Extended Emax model with Hill coefficient (gamma) controlling the steepness of the concentration-response curve.

---

## Function Signature

```python
neopkpd.simulate_pkpd_sigmoid_emax(
    cl: float,
    v: float,
    doses: list[dict],
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
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
| `ec50` | float | Concentration at 50% Emax |
| `gamma` | float | Hill coefficient (steepness) |
| `pk_kind` | str | PK model type |
| `ka` | float | Absorption rate (for oral) |

---

## Hill Equation

$$E(C) = E_0 + \frac{E_{max} \cdot C^\gamma}{EC_{50}^\gamma + C^\gamma}$$

### Effect of Gamma

| Gamma | Curve Shape |
|-------|-------------|
| < 1 | Shallow, gradual |
| = 1 | Standard Emax (hyperbolic) |
| > 1 | Steep, sigmoidal |
| >> 3 | Near-threshold |

---

## Basic Example

```python
import neopkpd

result = neopkpd.simulate_pkpd_sigmoid_emax(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    e0=0.0,
    emax=100.0,
    ec50=2.0,
    gamma=3.0,    # Steep response
    t0=0.0,
    t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

effect = result['observations']['effect']
t = result['t']

print(f"Max effect: {max(effect):.1f}")
```

---

## Comparing Gamma Values

```python
import neopkpd

gamma_values = [0.5, 1.0, 2.0, 3.0, 5.0]

print("Gamma | Max Effect | Time > 50% effect")
print("-" * 45)

for gamma in gamma_values:
    result = neopkpd.simulate_pkpd_sigmoid_emax(
        cl=5.0, v=50.0,
        doses=[{"time": 0.0, "amount": 500.0}],
        e0=0.0, emax=100.0, ec50=2.0, gamma=gamma,
        t0=0.0, t1=24.0,
        saveat=[i * 0.1 for i in range(241)]
    )

    effect = result['observations']['effect']
    t = result['t']

    max_effect = max(effect)

    # Time above 50% of Emax
    above_50 = [i for i, e in enumerate(effect) if e > 50]
    if above_50:
        duration = (t[above_50[-1]] - t[above_50[0]])
    else:
        duration = 0

    print(f"{gamma:5.1f} | {max_effect:10.1f} | {duration:6.1f} h")
```

---

## Steepness at EC50

Higher gamma = narrower therapeutic window:

```python
import neopkpd

# Calculate concentration at EC10 and EC90 for different gamma
ec50 = 2.0

print("Gamma | EC10 (mg/L) | EC90 (mg/L) | EC90/EC10")
print("-" * 50)

for gamma in [1.0, 2.0, 3.0, 5.0]:
    ec10 = ec50 * (0.1 / 0.9) ** (1/gamma)
    ec90 = ec50 * (0.9 / 0.1) ** (1/gamma)

    print(f"{gamma:5.1f} | {ec10:11.3f} | {ec90:11.2f} | {ec90/ec10:9.1f}")
```

**Key insight**: Higher gamma compresses the concentration range between EC10 and EC90.

---

## Clinical Example: Neuromuscular Blockade

```python
import neopkpd

# Rocuronium: steep response (gamma = 3-4)
result = neopkpd.simulate_pkpd_sigmoid_emax(
    cl=3.0, v=15.0,
    doses=[{"time": 0.0, "amount": 50.0}],
    e0=100.0,      # 100% baseline twitch
    emax=-100.0,   # Complete block possible
    ec50=1.0,      # mcg/mL
    gamma=3.5,     # Steep response
    t0=0.0, t1=2.0,  # 2 hours
    saveat=[i * 0.02 for i in range(101)]
)

effect = result['observations']['effect']
t = result['t']

# Time to onset (twitch < 10%)
onset_idx = next((i for i, e in enumerate(effect) if e < 10), None)
if onset_idx:
    print(f"Onset time (to <10% twitch): {t[onset_idx]*60:.1f} seconds")

# Duration of action (until twitch recovers to 25%)
recovery_idx = next((i for i, e in enumerate(effect[onset_idx:]) if e > 25), None)
if recovery_idx:
    duration = t[onset_idx + recovery_idx] - t[onset_idx]
    print(f"Duration of action: {duration*60:.1f} minutes")
```

---

## Therapeutic Index

```python
import neopkpd

# Drug with narrow therapeutic index (high gamma)
ec50_efficacy = 2.0
ec50_toxicity = 8.0  # 4× separation
gamma = 3.0

# Find therapeutic window (>80% efficacy with <20% toxicity)
def effect_at_conc(c, e0, emax, ec50, gamma):
    return e0 + emax * c**gamma / (ec50**gamma + c**gamma)

print("Concentration | Efficacy | Toxicity | Therapeutic?")
print("-" * 55)

for c in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]:
    eff = effect_at_conc(c, 0, 100, ec50_efficacy, gamma)
    tox = effect_at_conc(c, 0, 100, ec50_toxicity, gamma)

    therapeutic = "Yes" if eff > 80 and tox < 20 else "No"
    print(f"{c:13.1f} | {eff:8.1f}% | {tox:8.1f}% | {therapeutic}")
```

---

## Multiple Dosing with Steep Response

```python
import neopkpd

# With steep gamma, small concentration changes cause large effect changes
doses = [{"time": i * 8.0, "amount": 200.0} for i in range(6)]

result = neopkpd.simulate_pkpd_sigmoid_emax(
    cl=5.0, v=50.0,
    doses=doses,
    e0=0.0, emax=100.0, ec50=2.0, gamma=4.0,
    t0=0.0, t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

effect = result['observations']['effect']

# With high gamma, effect has "all-or-none" appearance
print(f"Peak effect: {max(effect):.1f}")
print(f"Trough effect: {min(effect):.1f}")
print(f"Fluctuation: {max(effect) - min(effect):.1f}")
```

---

## Visualization

```python
import neopkpd
from neopkpd.viz import plot_pkpd_profile

result = neopkpd.simulate_pkpd_sigmoid_emax(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    e0=0.0, emax=100.0, ec50=2.0, gamma=3.0,
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

fig = plot_pkpd_profile(
    result['t'],
    result['observations']['conc'],
    result['observations']['effect'],
    title="Sigmoid Emax PKPD (gamma=3)",
    effect_label="Effect (%)"
)
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Effect | $E_0 + \frac{E_{max} \cdot C^\gamma}{EC_{50}^\gamma + C^\gamma}$ |
| C at fraction f | $EC_{50} \cdot (f/(1-f))^{1/\gamma}$ |
| Slope at EC50 | $\gamma \cdot E_{max} / (4 \cdot EC_{50})$ |
| EC90/EC10 ratio | $81^{1/\gamma}$ |

---

## See Also

- [Direct Emax](direct-emax.md) - gamma = 1 case
- [Effect Compartment](effect-compartment.md) - With delay
- [Indirect Response](indirect-response.md) - Mechanism-based
