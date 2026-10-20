# Parallel First-Order Absorption

PK model for drugs with multiple absorption sites or mechanisms, each with distinct absorption rate constants.

---

## Function Signature

```python
neopkpd.simulate_pk_parallel_absorption(
    ka1: float,
    ka2: float,
    f1: float,
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
) -> dict
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `ka1` | float | First absorption rate constant (1/h) |
| `ka2` | float | Second absorption rate constant (1/h) |
| `f1` | float | Fraction to first depot (0-1, F2 = 1-F1) |
| `cl` | float | Clearance (L/h) |
| `v` | float | Volume of distribution (L) |
| `doses` | list[dict] | Dose events |

---

## Model Equations

Three-compartment system (two depot, one central):

$$\frac{dA_1}{dt} = -K_{a1} \cdot A_1$$

$$\frac{dA_2}{dt} = -K_{a2} \cdot A_2$$

$$\frac{dA_c}{dt} = K_{a1} \cdot A_1 + K_{a2} \cdot A_2 - \frac{CL}{V} \cdot A_c$$

$$C = \frac{A_c}{V}$$

Initial conditions after dose D:
- $A_1(0) = F_1 \cdot D$
- $A_2(0) = (1-F_1) \cdot D$

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "A_depot1": [...],    # Amount in first depot
        "A_depot2": [...],    # Amount in second depot
        "A_central": [...]    # Amount in central compartment
    },
    "observations": {
        "conc": [...]         # Plasma concentration
    },
    "metadata": {...}
}
```

---

## Basic Example

```python
import neopkpd

result = neopkpd.simulate_pk_parallel_absorption(
    ka1=2.0,      # Fast absorption (1/h)
    ka2=0.5,      # Slow absorption (1/h)
    f1=0.6,       # 60% to fast pathway
    cl=10.0,      # Clearance (L/h)
    v=50.0,       # Volume (L)
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

conc = result['observations']['conc']
t = result['t']

# Characteristic double-peak or shoulder
tmax_idx = max(range(len(conc)), key=lambda i: conc[i])
print(f"Tmax: {t[tmax_idx]:.2f} h")
print(f"Cmax: {conc[tmax_idx]:.2f} mg/L")
```

---

## Effect of Fraction Split

```python
import neopkpd

f1_values = [0.2, 0.4, 0.6, 0.8]

print("F1 (fast) | Cmax (mg/L) | Tmax (h)")
print("-" * 40)

for f1 in f1_values:
    result = neopkpd.simulate_pk_parallel_absorption(
        ka1=3.0, ka2=0.3, f1=f1, cl=10.0, v=50.0,
        doses=[{"time": 0.0, "amount": 500.0}],
        t0=0.0, t1=24.0,
        saveat=[i * 0.1 for i in range(241)]
    )

    conc = result['observations']['conc']
    t = result['t']

    cmax_idx = max(range(len(conc)), key=lambda i: conc[i])

    print(f"{f1:9.1f} | {conc[cmax_idx]:11.2f} | {t[cmax_idx]:8.2f}")
```

**Expected:** Higher F1 leads to higher, earlier Cmax.

---

## Biphasic Absorption Profile

```python
import neopkpd

# Strong separation between fast and slow components
result = neopkpd.simulate_pk_parallel_absorption(
    ka1=5.0,      # Very fast component
    ka2=0.2,      # Very slow component
    f1=0.3,       # 30% fast, 70% slow
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

conc = result['observations']['conc']
t = result['t']

# Look for double peak or shoulder
# Find local maxima
local_max = []
for i in range(1, len(conc)-1):
    if conc[i] > conc[i-1] and conc[i] > conc[i+1]:
        local_max.append((t[i], conc[i]))

print("Profile shape analysis:")
print(f"  Number of peaks/shoulders: {len(local_max)}")
for i, (time, c) in enumerate(local_max):
    print(f"  Peak {i+1}: t = {time:.2f} h, C = {c:.2f} mg/L")
```

---

## Clinical Example: Extended-Release Formulation

```python
import neopkpd

# ER formulation with immediate release coat + slow release core
result = neopkpd.simulate_pk_parallel_absorption(
    ka1=2.5,      # IR coat: fast
    ka2=0.15,     # ER core: slow
    f1=0.25,      # 25% IR, 75% ER
    cl=8.0,
    v=60.0,
    doses=[{"time": 0.0, "amount": 400.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.25 for i in range(97)]
)

conc = result['observations']['conc']
t = result['t']

# ER formulation characteristics
cmax_idx = max(range(len(conc)), key=lambda i: conc[i])
cmin = min(conc)

print("Extended-Release Profile:")
print(f"  Cmax: {conc[cmax_idx]:.2f} mg/L at {t[cmax_idx]:.1f} h")
print(f"  C at 12h: {conc[48]:.2f} mg/L")
print(f"  C at 24h: {conc[-1]:.2f} mg/L")
print(f"  Fluctuation: {(conc[cmax_idx] - conc[-1]) / conc[cmax_idx] * 100:.1f}%")
```

---

## Comparison with Single Absorption

```python
import neopkpd

# Parallel absorption
result_parallel = neopkpd.simulate_pk_parallel_absorption(
    ka1=3.0, ka2=0.5, f1=0.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.1 for i in range(241)]
)

# Single absorption (weighted average ka)
ka_avg = 0.5 * 3.0 + 0.5 * 0.5  # 1.75

result_single = neopkpd.simulate_pk_oral_first_order(
    ka=ka_avg, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.1 for i in range(241)]
)

parallel_conc = result_parallel['observations']['conc']
single_conc = result_single['observations']['conc']
t = result_parallel['t']

print("Parallel vs Single Absorption:")
print(f"  Parallel Cmax: {max(parallel_conc):.2f} mg/L")
print(f"  Single Cmax: {max(single_conc):.2f} mg/L")

# AUC should be similar (same total dose)
auc_parallel = sum(0.5 * (parallel_conc[i] + parallel_conc[i+1]) * 0.1
                   for i in range(len(parallel_conc)-1))
auc_single = sum(0.5 * (single_conc[i] + single_conc[i+1]) * 0.1
                 for i in range(len(single_conc)-1))
print(f"  Parallel AUC: {auc_parallel:.1f} mg*h/L")
print(f"  Single AUC: {auc_single:.1f} mg*h/L")
```

---

## Multiple Dosing

```python
import neopkpd

# Twice daily dosing
doses = [{"time": i * 12.0, "amount": 250.0} for i in range(6)]

result = neopkpd.simulate_pk_parallel_absorption(
    ka1=2.0, ka2=0.3, f1=0.4, cl=10.0, v=50.0,
    doses=doses,
    t0=0.0, t1=72.0,
    saveat=[i * 0.25 for i in range(289)]
)

conc = result['observations']['conc']
t = result['t']

# Steady-state characteristics (last dosing interval)
ss_start = 48  # Start of last full interval
ss_conc = [conc[i] for i in range(int(ss_start/0.25), len(conc))]

print("Steady-State Profile (48-72h):")
print(f"  Cmax,ss: {max(ss_conc):.2f} mg/L")
print(f"  Cmin,ss: {min(ss_conc):.2f} mg/L")
print(f"  Fluctuation: {(max(ss_conc) - min(ss_conc)) / min(ss_conc) * 100:.1f}%")
```

---

## Food Effect Simulation

```python
import neopkpd

# Fasted: faster absorption, smaller fraction
result_fasted = neopkpd.simulate_pk_parallel_absorption(
    ka1=3.0, ka2=0.8, f1=0.7, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.1 for i in range(241)]
)

# Fed: slower absorption, larger slow fraction
result_fed = neopkpd.simulate_pk_parallel_absorption(
    ka1=1.5, ka2=0.3, f1=0.4, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.1 for i in range(241)]
)

fasted_conc = result_fasted['observations']['conc']
fed_conc = result_fed['observations']['conc']
t = result_fasted['t']

print("Food Effect on Parallel Absorption:")
print(f"  Fasted Cmax: {max(fasted_conc):.2f} mg/L")
print(f"  Fed Cmax: {max(fed_conc):.2f} mg/L")
print(f"  Cmax ratio (Fed/Fasted): {max(fed_conc)/max(fasted_conc):.2f}")

# Tmax
tmax_fasted = t[max(range(len(fasted_conc)), key=lambda i: fasted_conc[i])]
tmax_fed = t[max(range(len(fed_conc)), key=lambda i: fed_conc[i])]
print(f"  Fasted Tmax: {tmax_fasted:.2f} h")
print(f"  Fed Tmax: {tmax_fed:.2f} h")
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Depot 1 initial | $A_1(0) = F_1 \cdot Dose$ |
| Depot 2 initial | $A_2(0) = (1-F_1) \cdot Dose$ |
| Total absorption rate | $K_{a1} \cdot A_1 + K_{a2} \cdot A_2$ |
| Concentration | $C = A_c / V$ |
| Elimination | $CL \cdot C = (CL/V) \cdot A_c$ |

---

## See Also

- [Oral First-Order](oral.md) - Single absorption pathway
- [Transit Absorption](transit.md) - Delayed absorption
- [Population Simulation](../../population/index.md) - Adding variability
