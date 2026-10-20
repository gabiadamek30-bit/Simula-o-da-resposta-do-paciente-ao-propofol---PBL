# Autoinduction

PK model for drugs that induce their own metabolism, leading to time-varying clearance that increases with chronic dosing.

---

## Function Signature

```python
neopkpd.simulate_pk_autoinduction(
    cl0: float,
    v: float,
    emax: float,
    ec50: float,
    kenz: float,
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
| `cl0` | float | Baseline clearance (L/h) |
| `v` | float | Volume of distribution (L) |
| `emax` | float | Maximum enzyme induction (fold increase) |
| `ec50` | float | Concentration for 50% of max induction |
| `kenz` | float | Enzyme turnover rate constant (1/h) |
| `doses` | list[dict] | Dose events |

### Derived Parameters

- **Enzyme half-life**: $t_{1/2,enz} = \ln(2) / k_{enz}$
- **Maximum clearance**: $CL_{max} = CL_0 \cdot (1 + E_{max})$

---

## Model Equations

Two-state system (drug and enzyme):

Enzyme induction signal:
$$E_{induced} = 1 + \frac{E_{max} \cdot C}{EC_{50} + C}$$

Enzyme dynamics (indirect response):
$$\frac{dE}{dt} = k_{enz} \cdot (E_{induced} - E)$$

Drug elimination with induced clearance:
$$\frac{dA_c}{dt} = -\frac{CL_0 \cdot E}{V} \cdot A_c$$

$$C = \frac{A_c}{V}$$

Initial condition: $E(0) = 1.0$ (baseline enzyme)

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "A_central": [...],   # Amount in central compartment
        "E_enzyme": [...]     # Enzyme level (relative to baseline)
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

# Single dose - minimal induction
result = neopkpd.simulate_pk_autoinduction(
    cl0=10.0,     # Baseline clearance (L/h)
    v=50.0,       # Volume (L)
    emax=2.0,     # Max 3-fold increase in enzyme
    ec50=5.0,     # EC50 for induction
    kenz=0.1,     # Enzyme turnover (~7 hour half-life)
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=48.0,
    saveat=[i * 0.5 for i in range(97)]
)

conc = result['observations']['conc']
enzyme = result['states']['E_enzyme']
t = result['t']

print("Single Dose Profile:")
print(f"  Initial enzyme: {enzyme[0]:.2f}")
print(f"  Peak enzyme: {max(enzyme):.2f}")
print(f"  Enzyme at 48h: {enzyme[-1]:.2f}")
```

---

## Chronic Dosing - Time-Varying Clearance

```python
import neopkpd

# Daily dosing for 2 weeks
doses = [{"time": i * 24.0, "amount": 200.0} for i in range(14)]

result = neopkpd.simulate_pk_autoinduction(
    cl0=10.0, v=50.0, emax=2.0, ec50=5.0, kenz=0.05,
    doses=doses,
    t0=0.0, t1=336.0,  # 14 days
    saveat=[i * 1.0 for i in range(337)]
)

conc = result['observations']['conc']
enzyme = result['states']['E_enzyme']
t = result['t']

# Compare day 1 vs day 14
day1_cmax = max(conc[:24])
day14_cmax = max(conc[312:336])

print("Autoinduction Effect:")
print(f"  Day 1 Cmax: {day1_cmax:.2f} mg/L")
print(f"  Day 14 Cmax: {day14_cmax:.2f} mg/L")
print(f"  Cmax decrease: {(1 - day14_cmax/day1_cmax) * 100:.1f}%")
print(f"  Enzyme level day 14: {enzyme[312]:.2f}x baseline")
print(f"  Effective clearance: {10.0 * enzyme[312]:.1f} L/h")
```

---

## Clinical Example: Carbamazepine

```python
import neopkpd

# Carbamazepine-like autoinduction
# Enzyme t1/2 ~3-5 days, full induction ~2-4 weeks

kenz = 0.693 / (4 * 24)  # ~4 day enzyme half-life

# BID dosing for 4 weeks
doses = [{"time": i * 12.0, "amount": 200.0} for i in range(56)]

result = neopkpd.simulate_pk_autoinduction(
    cl0=2.0,      # Initial low clearance
    v=80.0,       # Volume
    emax=1.5,     # Can increase enzyme 2.5-fold
    ec50=4.0,     # EC50 for induction
    kenz=kenz,
    doses=doses,
    t0=0.0, t1=672.0,  # 4 weeks
    saveat=[i * 2.0 for i in range(337)]
)

conc = result['observations']['conc']
enzyme = result['states']['E_enzyme']
t = result['t']

# Weekly Cmax comparison
print("Carbamazepine-like Autoinduction:")
for week in range(4):
    start = week * 168
    end = start + 168
    start_idx = int(start / 2)
    end_idx = int(end / 2)
    week_cmax = max(conc[start_idx:end_idx])
    week_enzyme = enzyme[end_idx - 1]
    print(f"  Week {week+1}: Cmax = {week_cmax:.2f} mg/L, Enzyme = {week_enzyme:.2f}x")
```

---

## Dose Adjustment for Autoinduction

```python
import neopkpd

# Strategy: Increase dose over time to maintain target concentration
target_conc = 6.0  # Target trough

# Week 1: Start low
# Week 2: Increase
# Week 3+: Maintenance

kenz = 0.693 / (5 * 24)  # 5-day enzyme half-life

# Stepped dosing regimen
doses = (
    [{"time": i * 24.0, "amount": 200.0} for i in range(7)] +      # Week 1
    [{"time": (7 + i) * 24.0, "amount": 300.0} for i in range(7)] +  # Week 2
    [{"time": (14 + i) * 24.0, "amount": 400.0} for i in range(7)]   # Week 3
)

result = neopkpd.simulate_pk_autoinduction(
    cl0=5.0, v=50.0, emax=1.5, ec50=4.0, kenz=kenz,
    doses=doses,
    t0=0.0, t1=504.0,
    saveat=[i * 2.0 for i in range(253)]
)

conc = result['observations']['conc']
t = result['t']

# Check trough levels at end of each week
print("Dose Titration Strategy:")
for week, dose in enumerate([200, 300, 400]):
    trough_time = (week + 1) * 168 - 2  # 2h before next dose
    trough_idx = int(trough_time / 2)
    trough = conc[trough_idx]
    print(f"  Week {week+1} ({dose} mg QD): Trough = {trough:.2f} mg/L")
```

---

## Time to Steady State

```python
import neopkpd

# With autoinduction, steady state takes longer
kenz = 0.693 / (4 * 24)  # 4-day enzyme half-life

doses = [{"time": i * 24.0, "amount": 200.0} for i in range(42)]  # 6 weeks

result = neopkpd.simulate_pk_autoinduction(
    cl0=5.0, v=50.0, emax=2.0, ec50=3.0, kenz=kenz,
    doses=doses,
    t0=0.0, t1=1008.0,
    saveat=[i * 6.0 for i in range(169)]
)

enzyme = result['states']['E_enzyme']
t = result['t']

# Time to 90% of steady-state enzyme
e_ss = enzyme[-1]
e_90 = 1.0 + 0.9 * (e_ss - 1.0)

for i, e in enumerate(enzyme):
    if e >= e_90:
        print(f"Time to 90% enzyme steady state: {t[i]:.0f} h ({t[i]/24:.1f} days)")
        break

print(f"Final enzyme level: {e_ss:.2f}x baseline")
```

---

## Drug Washout - Enzyme Recovery

```python
import neopkpd

kenz = 0.693 / (5 * 24)

# 2 weeks dosing, then stop
doses = [{"time": i * 24.0, "amount": 200.0} for i in range(14)]

result = neopkpd.simulate_pk_autoinduction(
    cl0=5.0, v=50.0, emax=2.0, ec50=3.0, kenz=kenz,
    doses=doses,
    t0=0.0, t1=672.0,  # Continue simulation 2 more weeks
    saveat=[i * 4.0 for i in range(169)]
)

enzyme = result['states']['E_enzyme']
t = result['t']

print("Enzyme Recovery After Drug Discontinuation:")
print(f"  Enzyme at day 14 (stop): {enzyme[84]:.2f}x")
print(f"  Enzyme at day 21: {enzyme[126]:.2f}x")
print(f"  Enzyme at day 28: {enzyme[-1]:.2f}x")
print(f"  Recovery t1/2: ~{0.693/kenz/24:.1f} days")
```

---

## Comparison: With vs Without Autoinduction

```python
import neopkpd

# Daily dosing for 2 weeks
doses = [{"time": i * 24.0, "amount": 200.0} for i in range(14)]

# With autoinduction
result_auto = neopkpd.simulate_pk_autoinduction(
    cl0=5.0, v=50.0, emax=1.5, ec50=3.0, kenz=0.693/(3*24),
    doses=doses,
    t0=0.0, t1=336.0,
    saveat=[i * 1.0 for i in range(337)]
)

# Without autoinduction (standard IV)
result_linear = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=doses,
    t0=0.0, t1=336.0,
    saveat=[i * 1.0 for i in range(337)]
)

auto_conc = result_auto['observations']['conc']
linear_conc = result_linear['observations']['conc']

print("Day 14 Comparison:")
print(f"  Linear PK Cmax: {max(linear_conc[312:336]):.2f} mg/L")
print(f"  Autoinduction Cmax: {max(auto_conc[312:336]):.2f} mg/L")

# AUC comparison
auc_linear = sum(linear_conc[312:336])
auc_auto = sum(auto_conc[312:336])
print(f"  Linear AUC (day 14): {auc_linear:.1f} mg*h/L")
print(f"  Autoinduction AUC (day 14): {auc_auto:.1f} mg*h/L")
print(f"  AUC reduction: {(1 - auc_auto/auc_linear) * 100:.1f}%")
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Induction signal | $E_{induced} = 1 + E_{max} \cdot C / (EC_{50} + C)$ |
| Enzyme rate | $dE/dt = k_{enz} \cdot (E_{induced} - E)$ |
| Effective clearance | $CL_{eff} = CL_0 \cdot E$ |
| Max clearance | $CL_{max} = CL_0 \cdot (1 + E_{max})$ |
| Enzyme t1/2 | $\ln(2) / k_{enz}$ |

---

## See Also

- [IV Bolus](iv-bolus.md) - Linear PK model
- [Michaelis-Menten](michaelis-menten.md) - Saturable elimination
- [Population Simulation](../../population/index.md) - Adding variability
