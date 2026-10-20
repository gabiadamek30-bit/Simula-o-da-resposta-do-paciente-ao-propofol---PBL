# Enterohepatic Recirculation (EHR)

PK model for drugs that undergo biliary excretion and intestinal reabsorption, leading to secondary concentration peaks.

---

## Function Signature

```python
neopkpd.simulate_pk_enterohepatic_recirculation(
    ka: float,
    cl: float,
    v: float,
    kbile: float,
    kreab: float,
    f_reab: float,
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
| `ka` | float | Absorption rate constant (1/h) |
| `cl` | float | Clearance (L/h) |
| `v` | float | Volume of distribution (L) |
| `kbile` | float | Biliary excretion rate constant (1/h) |
| `kreab` | float | Reabsorption rate constant from bile (1/h) |
| `f_reab` | float | Fraction reabsorbed (0-1) |
| `doses` | list[dict] | Dose events |

---

## Model Equations

Three-compartment system (GI tract, central, bile):

$$\frac{dA_{gut}}{dt} = -K_a \cdot A_{gut} + F_{reab} \cdot K_{reab} \cdot A_{bile}$$

$$\frac{dA_c}{dt} = K_a \cdot A_{gut} - \frac{CL}{V} \cdot A_c - K_{bile} \cdot A_c$$

$$\frac{dA_{bile}}{dt} = K_{bile} \cdot A_c - K_{reab} \cdot A_{bile}$$

$$C = \frac{A_c}{V}$$

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "A_gut": [...],       # Amount in GI tract
        "A_central": [...],   # Amount in central compartment
        "A_bile": [...]       # Amount in biliary compartment
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

result = neopkpd.simulate_pk_enterohepatic_recirculation(
    ka=1.0,       # Absorption rate (1/h)
    cl=10.0,      # Clearance (L/h)
    v=50.0,       # Volume (L)
    kbile=0.5,    # Biliary excretion rate (1/h)
    kreab=0.3,    # Reabsorption rate (1/h)
    f_reab=0.7,   # 70% reabsorbed
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=48.0,
    saveat=[i * 0.25 for i in range(193)]
)

conc = result['observations']['conc']
t = result['t']

# Look for secondary peaks
print("Concentration profile with EHR:")
print(f"  Initial peak: {max(conc[:20]):.2f} mg/L")
```

---

## Secondary Peak Detection

```python
import neopkpd

result = neopkpd.simulate_pk_enterohepatic_recirculation(
    ka=1.5, cl=8.0, v=50.0,
    kbile=0.4, kreab=0.2, f_reab=0.8,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

conc = result['observations']['conc']
t = result['t']

# Find local maxima (peaks)
peaks = []
for i in range(2, len(conc)-2):
    if (conc[i] > conc[i-1] and conc[i] > conc[i-2] and
        conc[i] > conc[i+1] and conc[i] > conc[i+2]):
        peaks.append((t[i], conc[i]))

print("Detected peaks:")
for i, (time, c) in enumerate(peaks):
    print(f"  Peak {i+1}: t = {time:.1f} h, C = {c:.2f} mg/L")
```

---

## Effect of Reabsorption Fraction

```python
import neopkpd

f_reab_values = [0.0, 0.3, 0.6, 0.9]

print("F_reab | AUC (mg*h/L) | Secondary peaks?")
print("-" * 50)

for f_reab in f_reab_values:
    result = neopkpd.simulate_pk_enterohepatic_recirculation(
        ka=1.5, cl=10.0, v=50.0,
        kbile=0.5, kreab=0.3, f_reab=f_reab,
        doses=[{"time": 0.0, "amount": 500.0}],
        t0=0.0, t1=72.0,
        saveat=[i * 0.1 for i in range(721)]
    )

    conc = result['observations']['conc']
    t = result['t']

    # Calculate AUC
    auc = sum(0.5 * (conc[i] + conc[i+1]) * 0.1
              for i in range(len(conc)-1))

    # Count peaks after initial
    peaks = 0
    for i in range(30, len(conc)-2):  # Skip first 3 hours
        if conc[i] > conc[i-1] and conc[i] > conc[i+1]:
            if conc[i] > conc[i-1] * 1.05:  # 5% increase threshold
                peaks += 1

    print(f"{f_reab:6.1f} | {auc:12.1f} | {'Yes' if peaks > 0 else 'No'}")
```

**Expected:** Higher F_reab leads to larger AUC and more pronounced secondary peaks.

---

## Clinical Example: Digoxin

```python
import neopkpd

# Digoxin-like EHR parameters
result = neopkpd.simulate_pk_enterohepatic_recirculation(
    ka=0.8,       # Absorption rate
    cl=7.0,       # Clearance (L/h)
    v=500.0,      # Large Vd (tissue binding)
    kbile=0.15,   # Biliary excretion
    kreab=0.1,    # Reabsorption rate
    f_reab=0.6,   # 60% reabsorbed
    doses=[{"time": 0.0, "amount": 0.5}],  # 0.5 mg dose
    t0=0.0, t1=72.0,
    saveat=[i * 0.5 for i in range(145)]
)

conc = result['observations']['conc']
t = result['t']

# Convert to ng/mL (typical clinical units)
conc_ng = [c * 1000 for c in conc]

print("Digoxin-like Profile:")
print(f"  Cmax: {max(conc_ng):.2f} ng/mL")
print(f"  C at 24h: {conc_ng[48]:.2f} ng/mL")
print(f"  C at 48h: {conc_ng[96]:.2f} ng/mL")
print("  Therapeutic range: 0.8-2.0 ng/mL")
```

---

## Meal-Triggered Gallbladder Emptying

```python
import neopkpd

# Model gallbladder emptying at mealtimes
# Note: This is a simplified approach - true meal effect would need
# time-varying parameters

# Morning dose, evening meal triggers bile release
result = neopkpd.simulate_pk_enterohepatic_recirculation(
    ka=1.2, cl=10.0, v=50.0,
    kbile=0.3,    # Moderate biliary excretion
    kreab=0.5,    # Fast reabsorption when released
    f_reab=0.75,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=24.0,
    saveat=[i * 0.1 for i in range(241)]
)

conc = result['observations']['conc']
a_bile = result['states']['A_bile']
t = result['t']

print("Biliary Compartment Dynamics:")
print(f"  Max bile accumulation: {max(a_bile):.1f} mg at t={t[a_bile.index(max(a_bile))]:.1f} h")
print(f"  This bile would be released with next meal")
```

---

## Comparison: With vs Without EHR

```python
import neopkpd

# With EHR
result_ehr = neopkpd.simulate_pk_enterohepatic_recirculation(
    ka=1.5, cl=10.0, v=50.0,
    kbile=0.4, kreab=0.3, f_reab=0.7,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

# Without EHR (standard oral model)
result_no_ehr = neopkpd.simulate_pk_oral_first_order(
    ka=1.5, cl=10.0, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=48.0,
    saveat=[i * 0.1 for i in range(481)]
)

ehr_conc = result_ehr['observations']['conc']
no_ehr_conc = result_no_ehr['observations']['conc']
t = result_ehr['t']

# Calculate AUCs
auc_ehr = sum(0.5 * (ehr_conc[i] + ehr_conc[i+1]) * 0.1
              for i in range(len(ehr_conc)-1))
auc_no_ehr = sum(0.5 * (no_ehr_conc[i] + no_ehr_conc[i+1]) * 0.1
                 for i in range(len(no_ehr_conc)-1))

print("Effect of Enterohepatic Recirculation:")
print(f"  AUC without EHR: {auc_no_ehr:.1f} mg*h/L")
print(f"  AUC with EHR: {auc_ehr:.1f} mg*h/L")
print(f"  AUC increase: {(auc_ehr/auc_no_ehr - 1) * 100:.1f}%")

# Terminal half-life
# Find time to reach 10% of Cmax
cmax_ehr = max(ehr_conc)
cmax_no_ehr = max(no_ehr_conc)
print(f"  Cmax without EHR: {cmax_no_ehr:.2f} mg/L")
print(f"  Cmax with EHR: {cmax_ehr:.2f} mg/L")
```

---

## Multiple Dosing with EHR

```python
import neopkpd

# Once daily dosing
doses = [{"time": i * 24.0, "amount": 200.0} for i in range(5)]

result = neopkpd.simulate_pk_enterohepatic_recirculation(
    ka=1.0, cl=8.0, v=50.0,
    kbile=0.3, kreab=0.2, f_reab=0.65,
    doses=doses,
    t0=0.0, t1=120.0,
    saveat=[i * 0.5 for i in range(241)]
)

conc = result['observations']['conc']
t = result['t']

# Day 5 (steady state)
day5_start = int(96 / 0.5)
day5_conc = conc[day5_start:]

print("Multiple Dose Profile (Day 5):")
print(f"  Cmax,ss: {max(day5_conc):.2f} mg/L")
print(f"  Cmin,ss: {min(day5_conc):.2f} mg/L")
print(f"  Average: {sum(day5_conc)/len(day5_conc):.2f} mg/L")
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| GI absorption | $K_a \cdot A_{gut}$ |
| Biliary excretion | $K_{bile} \cdot A_c$ |
| Reabsorption | $F_{reab} \cdot K_{reab} \cdot A_{bile}$ |
| Net elimination | $CL \cdot C + (1-F_{reab}) \cdot K_{reab} \cdot A_{bile}$ |
| Effective half-life | Prolonged due to recycling |

---

## See Also

- [Oral First-Order](oral.md) - Without EHR
- [Transit Absorption](transit.md) - Delayed absorption
- [Population Simulation](../../population/index.md) - Adding variability
