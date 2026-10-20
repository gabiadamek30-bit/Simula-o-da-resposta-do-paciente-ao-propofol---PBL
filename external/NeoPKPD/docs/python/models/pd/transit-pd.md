# Transit Compartment PD Model

PD model with a chain of transit compartments to model delayed drug effects and signal transduction cascades.

---

## Function Signature

```python
neopkpd.simulate_pkpd_transit_compartment(
    cl: float,
    v: float,
    doses: list[dict],
    n_transit: int,
    ktr: float,
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
    t0: float,
    t1: float,
    saveat: list[float],
    pk_kind: str = "OneCompIVBolus",
    ka: float | None = None,
    q: float | None = None,
    v2: float | None = None,
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
| `n_transit` | int | Number of transit compartments (1-20 typical) |
| `ktr` | float | Transit rate constant (1/h) |
| `e0` | float | Baseline effect/signal |
| `emax` | float | Maximum effect above baseline |
| `ec50` | float | Concentration at 50% of Emax |
| `gamma` | float | Hill coefficient (steepness) |

### Key Derived Parameter

**Mean Transit Time (MTT)**:
$$MTT = \frac{N + 1}{K_{tr}}$$

---

## Model Equations

Drug-induced signal:
$$Signal(C) = E_0 + \frac{E_{max} \cdot C^\gamma}{EC_{50}^\gamma + C^\gamma}$$

Transit compartment chain:
$$\frac{dA_1}{dt} = K_{tr} \cdot (Signal(C) - A_1)$$

$$\frac{dA_i}{dt} = K_{tr} \cdot (A_{i-1} - A_i) \quad \text{for } i = 2..N$$

Final effect:
$$Effect = A_N$$

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "A_central": [...],   # Drug amount
        "transit_1": [...],   # Transit compartments
        ...
    },
    "observations": {
        "conc": [...],        # Drug concentration
        "effect": [...]       # Final effect
    },
    "metadata": {...}
}
```

---

## Basic Example

```python
import neopkpd

result = neopkpd.simulate_pkpd_transit_compartment(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    n_transit=5,      # 5 transit compartments
    ktr=0.5,          # Transit rate
    e0=100.0,         # Baseline effect
    emax=-50.0,       # Maximum decrease
    ec50=5.0,
    gamma=1.0,
    t0=0.0,
    t1=168.0,         # 1 week
    saveat=[i * 1.0 for i in range(169)]
)

effect = result['observations']['effect']
conc = result['observations']['conc']
t = result['t']

# MTT = (5+1) / 0.5 = 12 hours
print(f"Mean Transit Time: {(5+1)/0.5:.1f} h")
print(f"Max effect: {min(effect):.1f}")  # Minimum because emax is negative
print(f"Time to max effect: {t[effect.index(min(effect))]} h")
```

---

## Effect of Number of Transit Compartments

```python
import neopkpd

n_transit_values = [1, 3, 5, 10, 20]
ktr = 0.5  # Fixed ktr

print("N transit | MTT (h) | T_nadir (h) | Nadir effect")
print("-" * 55)

for n in n_transit_values:
    result = neopkpd.simulate_pkpd_transit_compartment(
        cl=5.0, v=50.0,
        doses=[{"time": 0.0, "amount": 500.0}],
        n_transit=n, ktr=ktr, e0=100.0, emax=-60.0, ec50=5.0, gamma=1.0,
        t0=0.0, t1=168.0,
        saveat=[i * 0.5 for i in range(337)]
    )

    effect = result['observations']['effect']
    t = result['t']
    mtt = (n + 1) / ktr

    nadir = min(effect)
    t_nadir = t[effect.index(nadir)]

    print(f"{n:9d} | {mtt:7.1f} | {t_nadir:11.1f} | {nadir:12.1f}")
```

**Expected:** More compartments = sharper nadir, closer to MTT.

---

## Clinical Example: Myelosuppression

```python
import neopkpd

# Chemotherapy-induced neutropenia
# Typical: 5-7 transit compartments, MTT ~5 days for neutrophils

mtt_days = 5.0
n_transit = 5
ktr = (n_transit + 1) / (mtt_days * 24)  # Convert to 1/h

result = neopkpd.simulate_pkpd_transit_compartment(
    cl=3.0, v=30.0,
    doses=[{"time": 0.0, "amount": 100.0}],  # Single chemo dose
    n_transit=n_transit,
    ktr=ktr,
    e0=5.0,           # Baseline ANC (10^9/L)
    emax=-4.5,        # Max 90% reduction
    ec50=1.0,
    gamma=1.5,
    t0=0.0, t1=504.0,  # 3 weeks
    saveat=[i * 2.0 for i in range(253)]
)

effect = result['observations']['effect']
t = result['t']

# Grade neutropenia thresholds
print("Neutropenia Profile:")
print(f"  Baseline ANC: 5.0 x10^9/L")
nadir = min(effect)
t_nadir = t[effect.index(nadir)]
print(f"  Nadir ANC: {nadir:.2f} x10^9/L at day {t_nadir/24:.1f}")

# Time to recovery
for i, e in enumerate(effect):
    if i > effect.index(nadir) and e > 1.5:
        print(f"  Recovery to >1.5 x10^9/L: day {t[i]/24:.1f}")
        break

# Grade classification
if nadir < 0.5:
    print("  Grade 4 (severe)")
elif nadir < 1.0:
    print("  Grade 3 (moderate)")
elif nadir < 1.5:
    print("  Grade 2 (mild)")
```

---

## Multiple Chemotherapy Cycles

```python
import neopkpd

# Every 3 weeks (21 days)
cycle_length = 21 * 24  # hours
n_cycles = 4

doses = [{"time": i * cycle_length, "amount": 100.0} for i in range(n_cycles)]

mtt_days = 5.0
n_transit = 5
ktr = (n_transit + 1) / (mtt_days * 24)

result = neopkpd.simulate_pkpd_transit_compartment(
    cl=3.0, v=30.0,
    doses=doses,
    n_transit=n_transit, ktr=ktr,
    e0=5.0, emax=-4.0, ec50=1.0, gamma=1.5,
    t0=0.0, t1=n_cycles * cycle_length + 168,
    saveat=[i * 4.0 for i in range(500)]
)

effect = result['observations']['effect']
t = result['t']

print("Cycle-by-Cycle Nadir:")
for cycle in range(n_cycles):
    start_h = cycle * cycle_length
    end_h = (cycle + 1) * cycle_length
    start_idx = int(start_h / 4)
    end_idx = int(end_h / 4)

    cycle_effect = effect[start_idx:end_idx]
    nadir = min(cycle_effect)
    nadir_t = t[start_idx + cycle_effect.index(nadir)]

    print(f"  Cycle {cycle+1}: Nadir = {nadir:.2f} at day {nadir_t/24:.0f}")
```

---

## Thrombocytopenia Model

```python
import neopkpd

# Platelet dynamics: longer MTT than neutrophils (~8-10 days)
mtt_days = 8.0
n_transit = 5
ktr = (n_transit + 1) / (mtt_days * 24)

result = neopkpd.simulate_pkpd_transit_compartment(
    cl=3.0, v=30.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    n_transit=n_transit, ktr=ktr,
    e0=250.0,         # Baseline platelets (10^9/L)
    emax=-200.0,      # Max 80% reduction
    ec50=1.0,
    gamma=2.0,
    t0=0.0, t1=672.0,  # 4 weeks
    saveat=[i * 2.0 for i in range(337)]
)

effect = result['observations']['effect']
t = result['t']

nadir = min(effect)
t_nadir = t[effect.index(nadir)]

print("Thrombocytopenia Profile:")
print(f"  Baseline platelets: 250 x10^9/L")
print(f"  Nadir: {nadir:.0f} x10^9/L at day {t_nadir/24:.1f}")

# Recovery
for i, e in enumerate(effect):
    if i > effect.index(nadir) and e > 100:
        print(f"  Recovery to >100 x10^9/L: day {t[i]/24:.1f}")
        break
```

---

## Varying MTT with Fixed N

```python
import neopkpd

mtt_values = [12, 24, 48, 72, 120]  # hours
n_transit = 5

print("MTT (h) | Ktr (1/h) | T_nadir (h)")
print("-" * 40)

for mtt in mtt_values:
    ktr = (n_transit + 1) / mtt

    result = neopkpd.simulate_pkpd_transit_compartment(
        cl=5.0, v=50.0,
        doses=[{"time": 0.0, "amount": 500.0}],
        n_transit=n_transit, ktr=ktr,
        e0=100.0, emax=-50.0, ec50=5.0, gamma=1.0,
        t0=0.0, t1=336.0,
        saveat=[i * 1.0 for i in range(337)]
    )

    effect = result['observations']['effect']
    t = result['t']
    t_nadir = t[effect.index(min(effect))]

    print(f"{mtt:7d} | {ktr:9.4f} | {t_nadir:11.1f}")
```

---

## Stimulation Model (G-CSF Effect)

```python
import neopkpd

# G-CSF stimulates neutrophil production
# Positive Emax = stimulation

result = neopkpd.simulate_pkpd_transit_compartment(
    cl=0.5, v=5.0,     # Typical G-CSF PK
    doses=[{"time": i * 24, "amount": 5.0} for i in range(7)],  # Daily dosing
    n_transit=5,
    ktr=0.05,          # Slow effect
    e0=3.0,            # Low baseline (neutropenic patient)
    emax=10.0,         # Stimulation
    ec50=0.1,
    gamma=1.0,
    t0=0.0, t1=336.0,
    saveat=[i * 2.0 for i in range(169)]
)

effect = result['observations']['effect']
t = result['t']

print("G-CSF Effect on ANC:")
print(f"  Baseline ANC: 3.0 x10^9/L")
print(f"  Peak ANC: {max(effect):.1f} x10^9/L")
print(f"  Day 7 ANC: {effect[84]:.1f} x10^9/L")
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| MTT | $(N+1) / K_{tr}$ |
| Signal | $E_0 + E_{max} \cdot C^\gamma / (EC_{50}^\gamma + C^\gamma)$ |
| Transit rate | $K_{tr} \cdot (A_{i-1} - A_i)$ |
| CV of delay | $1 / \sqrt{N+1}$ |
| Peak effect lag | $\approx MTT$ from peak concentration |

---

## See Also

- [Indirect Response](indirect-response.md) - Simpler delay mechanism
- [Effect Compartment](effect-compartment.md) - Single delay compartment
- [Disease Progression](disease-progression.md) - Tumor growth models
