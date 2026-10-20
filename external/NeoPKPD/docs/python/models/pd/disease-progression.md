# Disease Progression Model

PD model for tumor growth dynamics with drug-induced cell kill, supporting multiple growth models.

---

## Function Signature

```python
neopkpd.simulate_pkpd_disease_progression(
    cl: float,
    v: float,
    doses: list[dict],
    growth_model: str,
    s0: float,
    kgrow: float,
    smax: float,
    alpha: float,
    kdrug: float,
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
| `growth_model` | str | Growth model type (see below) |
| `s0` | float | Initial tumor size |
| `kgrow` | float | Growth rate constant |
| `smax` | float | Maximum size (carrying capacity) |
| `alpha` | float | Linear growth rate (for linear model) |
| `kdrug` | float | Drug-induced cell kill rate constant |

### Growth Models

| Model | Equation |
|-------|----------|
| `"exponential"` | $dS/dt = k_{grow} \cdot S - k_{drug} \cdot C \cdot S$ |
| `"linear"` | $dS/dt = \alpha - k_{drug} \cdot C \cdot S$ |
| `"logistic"` | $dS/dt = k_{grow} \cdot S \cdot (1 - S/S_{max}) - k_{drug} \cdot C \cdot S$ |
| `"gompertz"` | $dS/dt = k_{grow} \cdot S \cdot \ln(S_{max}/S) - k_{drug} \cdot C \cdot S$ |
| `"asymptotic"` | $dS/dt = k_{grow} \cdot (S_{max} - S) - k_{drug} \cdot C \cdot S$ |

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "A_central": [...],   # Drug amount
        "S": [...]            # Tumor size
    },
    "observations": {
        "conc": [...],        # Drug concentration
        "tumor_size": [...]   # Tumor size
    },
    "metadata": {...}
}
```

---

## Basic Example: Exponential Growth

```python
import neopkpd

result = neopkpd.simulate_pkpd_disease_progression(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    growth_model="exponential",
    s0=100.0,         # Initial tumor size
    kgrow=0.02,       # Growth rate (1/h) = doubling ~35h
    smax=1000.0,      # Not used for exponential
    alpha=0.0,        # Not used for exponential
    kdrug=0.005,      # Drug kill rate
    t0=0.0,
    t1=336.0,         # 2 weeks
    saveat=[i * 2.0 for i in range(169)]
)

tumor = result['observations']['tumor_size']
t = result['t']

print("Exponential Growth Model:")
print(f"  Initial size: {s0:.0f}")
print(f"  Size at 7 days: {tumor[84]:.1f}")
print(f"  Size at 14 days: {tumor[-1]:.1f}")
```

---

## Comparing Growth Models

```python
import neopkpd

models = ["exponential", "logistic", "gompertz", "asymptotic"]
s0 = 100.0

print("Growth Model | Day 7 Size | Day 14 Size | Max Size")
print("-" * 55)

for model in models:
    result = neopkpd.simulate_pkpd_disease_progression(
        cl=5.0, v=50.0,
        doses=[],  # No treatment
        growth_model=model,
        s0=s0, kgrow=0.02, smax=1000.0, alpha=1.0, kdrug=0.0,
        t0=0.0, t1=672.0,  # 4 weeks
        saveat=[i * 4.0 for i in range(169)]
    )

    tumor = result['observations']['tumor_size']
    max_size = max(tumor)

    print(f"{model:12s} | {tumor[42]:.1f} | {tumor[84]:.1f} | {max_size:.1f}")
```

---

## Treatment Effect

```python
import neopkpd

# Weekly dosing
doses = [{"time": i * 168, "amount": 200.0} for i in range(4)]

# No treatment baseline
result_no_tx = neopkpd.simulate_pkpd_disease_progression(
    cl=5.0, v=50.0, doses=[],
    growth_model="gompertz",
    s0=100.0, kgrow=0.02, smax=1000.0, alpha=0.0, kdrug=0.0,
    t0=0.0, t1=672.0,
    saveat=[i * 4.0 for i in range(169)]
)

# With treatment
result_tx = neopkpd.simulate_pkpd_disease_progression(
    cl=5.0, v=50.0, doses=doses,
    growth_model="gompertz",
    s0=100.0, kgrow=0.02, smax=1000.0, alpha=0.0, kdrug=0.01,
    t0=0.0, t1=672.0,
    saveat=[i * 4.0 for i in range(169)]
)

no_tx = result_no_tx['observations']['tumor_size']
tx = result_tx['observations']['tumor_size']
t = result_no_tx['t']

print("Treatment Effect (Gompertz Model):")
print(f"  Day 0: {no_tx[0]:.1f}")
print(f"  Day 28 no treatment: {no_tx[-1]:.1f}")
print(f"  Day 28 with treatment: {tx[-1]:.1f}")
print(f"  Tumor growth inhibition: {(1 - tx[-1]/no_tx[-1]) * 100:.1f}%")
```

---

## Dose-Response Relationship

```python
import neopkpd

doses_list = [50, 100, 200, 400, 800]

print("Dose (mg) | Final Size | TGI (%)")
print("-" * 40)

# No treatment reference
result_ref = neopkpd.simulate_pkpd_disease_progression(
    cl=5.0, v=50.0, doses=[],
    growth_model="logistic",
    s0=100.0, kgrow=0.02, smax=1000.0, alpha=0.0, kdrug=0.0,
    t0=0.0, t1=336.0,
    saveat=[i * 4.0 for i in range(85)]
)
ref_final = result_ref['observations']['tumor_size'][-1]

for dose in doses_list:
    weekly_doses = [{"time": i * 168, "amount": float(dose)} for i in range(2)]

    result = neopkpd.simulate_pkpd_disease_progression(
        cl=5.0, v=50.0, doses=weekly_doses,
        growth_model="logistic",
        s0=100.0, kgrow=0.02, smax=1000.0, alpha=0.0, kdrug=0.01,
        t0=0.0, t1=336.0,
        saveat=[i * 4.0 for i in range(85)]
    )

    final = result['observations']['tumor_size'][-1]
    tgi = (1 - final / ref_final) * 100

    print(f"{dose:9d} | {final:10.1f} | {tgi:7.1f}")
```

---

## Complete Response Detection

```python
import neopkpd

# Intensive treatment
doses = [{"time": i * 24, "amount": 300.0} for i in range(14)]

result = neopkpd.simulate_pkpd_disease_progression(
    cl=5.0, v=50.0, doses=doses,
    growth_model="exponential",
    s0=100.0, kgrow=0.02, smax=1000.0, alpha=0.0, kdrug=0.02,
    t0=0.0, t1=672.0,
    saveat=[i * 2.0 for i in range(337)]
)

tumor = result['observations']['tumor_size']
t = result['t']

# Find if tumor shrinks below detection threshold
detection_threshold = 10.0
min_size = min(tumor)
t_min = t[tumor.index(min_size)]

print("Tumor Response:")
print(f"  Initial: {tumor[0]:.1f}")
print(f"  Minimum: {min_size:.2f} at day {t_min/24:.1f}")

if min_size < detection_threshold:
    print("  Response: Complete Response (CR)")
elif min_size < tumor[0] * 0.5:
    print("  Response: Partial Response (PR)")
elif min_size < tumor[0] * 1.25:
    print("  Response: Stable Disease (SD)")
else:
    print("  Response: Progressive Disease (PD)")
```

---

## Regrowth After Treatment

```python
import neopkpd

# Treatment for 2 weeks, then observe
doses = [{"time": i * 24, "amount": 200.0} for i in range(14)]

result = neopkpd.simulate_pkpd_disease_progression(
    cl=5.0, v=50.0, doses=doses,
    growth_model="gompertz",
    s0=100.0, kgrow=0.03, smax=1000.0, alpha=0.0, kdrug=0.015,
    t0=0.0, t1=1008.0,  # 6 weeks total
    saveat=[i * 4.0 for i in range(253)]
)

tumor = result['observations']['tumor_size']
t = result['t']

print("Tumor Dynamics:")
print(f"  Day 0: {tumor[0]:.1f}")
print(f"  Day 14 (end of treatment): {tumor[84]:.1f}")
print(f"  Day 28: {tumor[168]:.1f}")
print(f"  Day 42: {tumor[-1]:.1f}")

# Time to regrow to initial size
for i, s in enumerate(tumor):
    if i > 84 and s > tumor[0]:
        print(f"  Regrowth to initial size: day {t[i]/24:.0f}")
        break
```

---

## Combination Treatment

```python
import neopkpd

# Single agent
result_single = neopkpd.simulate_pkpd_disease_progression(
    cl=5.0, v=50.0,
    doses=[{"time": i * 168, "amount": 200.0} for i in range(4)],
    growth_model="logistic",
    s0=100.0, kgrow=0.025, smax=800.0, alpha=0.0, kdrug=0.008,
    t0=0.0, t1=672.0,
    saveat=[i * 4.0 for i in range(169)]
)

# Combination (simulate as higher effective kdrug)
result_combo = neopkpd.simulate_pkpd_disease_progression(
    cl=5.0, v=50.0,
    doses=[{"time": i * 168, "amount": 200.0} for i in range(4)],
    growth_model="logistic",
    s0=100.0, kgrow=0.025, smax=800.0, alpha=0.0, kdrug=0.015,
    t0=0.0, t1=672.0,
    saveat=[i * 4.0 for i in range(169)]
)

single = result_single['observations']['tumor_size']
combo = result_combo['observations']['tumor_size']

print("Single Agent vs Combination:")
print(f"  Single agent day 28: {single[-1]:.1f}")
print(f"  Combination day 28: {combo[-1]:.1f}")
print(f"  Additional benefit: {(single[-1] - combo[-1])/single[-1] * 100:.1f}%")
```

---

## Survival Surrogate

```python
import neopkpd

# Tumor doubling time as survival surrogate
doses = [{"time": i * 168, "amount": 200.0} for i in range(4)]

result = neopkpd.simulate_pkpd_disease_progression(
    cl=5.0, v=50.0, doses=doses,
    growth_model="exponential",
    s0=100.0, kgrow=0.025, smax=1000.0, alpha=0.0, kdrug=0.01,
    t0=0.0, t1=2016.0,  # 12 weeks
    saveat=[i * 8.0 for i in range(253)]
)

tumor = result['observations']['tumor_size']
t = result['t']

# Time to reach lethal tumor burden
lethal_burden = 1000.0
for i, s in enumerate(tumor):
    if s > lethal_burden:
        print(f"Time to lethal burden: {t[i]/24:.0f} days ({t[i]/168:.1f} weeks)")
        break
else:
    print(f"Tumor below lethal burden at end of simulation")
    print(f"Final tumor size: {tumor[-1]:.1f}")
```

---

## Equations Summary

| Model | dS/dt (without drug) | Steady State |
|-------|---------------------|--------------|
| Exponential | $k_{grow} \cdot S$ | Infinite |
| Linear | $\alpha$ | Infinite |
| Logistic | $k_{grow} \cdot S \cdot (1 - S/S_{max})$ | $S_{max}$ |
| Gompertz | $k_{grow} \cdot S \cdot \ln(S_{max}/S)$ | $S_{max}$ |
| Asymptotic | $k_{grow} \cdot (S_{max} - S)$ | $S_{max}$ |

Drug effect: $-k_{drug} \cdot C \cdot S$ (concentration-dependent cell kill)

---

## See Also

- [Transit Compartment PD](transit-pd.md) - Delayed effects
- [Indirect Response](indirect-response.md) - Turnover models
- [Population Simulation](../../population/index.md) - Adding variability
