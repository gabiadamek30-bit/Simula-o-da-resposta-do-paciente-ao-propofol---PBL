# Tolerance Models

PD models for tolerance development through counter-regulation or receptor regulation mechanisms.

---

## Counter-Regulation Tolerance

### Function Signature

```python
neopkpd.simulate_pkpd_tolerance_counter_regulation(
    cl: float,
    v: float,
    doses: list[dict],
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
    kin_mod: float,
    kout_mod: float,
    alpha_feedback: float,
    t0: float,
    t1: float,
    saveat: list[float],
    pk_kind: str = "OneCompIVBolus",
    ka: float | None = None,
    q: float | None = None,
    v2: float | None = None,
    ...
) -> dict
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `e0` | float | Baseline effect |
| `emax` | float | Maximum drug effect |
| `ec50` | float | EC50 for drug effect |
| `gamma` | float | Hill coefficient |
| `kin_mod` | float | Moderator production rate constant |
| `kout_mod` | float | Moderator elimination rate constant |
| `alpha_feedback` | float | Feedback strength coefficient |

### Model Equations

Drug effect:
$$E_{drug} = \frac{E_{max} \cdot C^\gamma}{EC_{50}^\gamma + C^\gamma}$$

Moderator dynamics (builds up with drug effect):
$$\frac{dM}{dt} = k_{in,mod} \cdot E_{drug} - k_{out,mod} \cdot M$$

Net effect with tolerance:
$$E_{net} = E_0 + E_{drug} - \alpha \cdot M$$

---

## Basic Counter-Regulation Example

```python
import neopkpd

# Chronic dosing with tolerance development
doses = [{"time": i * 8.0, "amount": 50.0} for i in range(21)]  # TID for 7 days

result = neopkpd.simulate_pkpd_tolerance_counter_regulation(
    cl=5.0,
    v=50.0,
    doses=doses,
    e0=0.0,
    emax=100.0,
    ec50=5.0,
    gamma=1.0,
    kin_mod=0.1,       # Moderator production
    kout_mod=0.05,     # Moderator elimination (t1/2 ~14h)
    alpha_feedback=1.0, # Full feedback
    t0=0.0,
    t1=168.0,          # 1 week
    saveat=[i * 0.5 for i in range(337)]
)

effect = result['observations']['effect']
t = result['t']

# Compare first dose vs last dose effect
first_dose_peak = max(effect[:16])  # First 8 hours
last_dose_peak = max(effect[-16:])  # Last 8 hours

print("Counter-Regulation Tolerance:")
print(f"  First dose peak effect: {first_dose_peak:.1f}")
print(f"  Last dose peak effect: {last_dose_peak:.1f}")
print(f"  Tolerance: {(1 - last_dose_peak/first_dose_peak) * 100:.1f}% reduction")
```

---

## Clinical Example: Opioid Tolerance

```python
import neopkpd

# Morphine-like tolerance
# Daily dosing for 2 weeks

doses = [{"time": i * 6.0, "amount": 10.0} for i in range(56)]  # Q6H for 2 weeks

result = neopkpd.simulate_pkpd_tolerance_counter_regulation(
    cl=60.0, v=200.0,    # Morphine-like PK
    doses=doses,
    e0=0.0,              # Baseline pain score
    emax=100.0,          # Max analgesia
    ec50=0.05,           # EC50 (mg/L)
    gamma=2.0,
    kin_mod=0.2,
    kout_mod=0.01,       # Slow moderator turnover (~3 day t1/2)
    alpha_feedback=0.8,
    t0=0.0, t1=336.0,
    saveat=[i * 1.0 for i in range(337)]
)

effect = result['observations']['effect']
t = result['t']

# Daily peak effect
print("Daily Peak Analgesia:")
for day in [1, 3, 7, 14]:
    start = (day - 1) * 24
    end = day * 24
    peak = max(effect[start:end])
    print(f"  Day {day:2d}: {peak:.1f}")
```

---

## Effect of Feedback Strength

```python
import neopkpd

alpha_values = [0.0, 0.25, 0.5, 1.0, 1.5]

doses = [{"time": i * 12.0, "amount": 50.0} for i in range(14)]

print("Alpha | Day 1 Peak | Day 7 Peak | Tolerance %")
print("-" * 50)

for alpha in alpha_values:
    result = neopkpd.simulate_pkpd_tolerance_counter_regulation(
        cl=5.0, v=50.0, doses=doses,
        e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
        kin_mod=0.1, kout_mod=0.05, alpha_feedback=alpha,
        t0=0.0, t1=168.0,
        saveat=[i * 0.5 for i in range(337)]
    )

    effect = result['observations']['effect']
    day1_peak = max(effect[:24])
    day7_peak = max(effect[-24:])
    tolerance = (1 - day7_peak/day1_peak) * 100 if day1_peak > 0 else 0

    print(f"{alpha:5.2f} | {day1_peak:10.1f} | {day7_peak:10.1f} | {tolerance:10.1f}")
```

---

## Receptor Regulation Model

### Function Signature

```python
neopkpd.simulate_pkpd_receptor_regulation(
    cl: float,
    v: float,
    doses: list[dict],
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
    r_baseline: float,
    kreg: float,
    rmax: float,
    kchange: float,
    direction: str,  # "down" or "up"
    t0: float,
    t1: float,
    saveat: list[float],
    ...
) -> dict
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `r_baseline` | float | Baseline receptor density (normalized, typically 1.0) |
| `kreg` | float | Receptor return-to-baseline rate constant |
| `rmax` | float | Maximum receptor density (for up-regulation) |
| `kchange` | float | Rate of receptor change |
| `direction` | str | "down" or "up" regulation |

### Model Equations

Receptor dynamics:
$$\frac{dR}{dt} = k_{reg} \cdot (R_{baseline} - R) + regulation\_effect$$

Down-regulation:
$$regulation\_effect = -k_{change} \cdot E_{drug} \cdot R$$

Up-regulation:
$$regulation\_effect = +k_{change} \cdot E_{drug} \cdot (R_{max} - R)$$

Net effect:
$$E_{net} = E_0 + R \cdot E_{drug}$$

---

## Receptor Down-Regulation Example

```python
import neopkpd

# Beta-receptor down-regulation with chronic agonist
doses = [{"time": i * 8.0, "amount": 50.0} for i in range(21)]

result = neopkpd.simulate_pkpd_receptor_regulation(
    cl=5.0, v=50.0, doses=doses,
    e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
    r_baseline=1.0,    # Normalized baseline
    kreg=0.05,         # Receptor recovery rate
    rmax=2.0,          # Not used for down-regulation
    kchange=0.02,      # Down-regulation rate
    direction="down",
    t0=0.0, t1=168.0,
    saveat=[i * 0.5 for i in range(337)]
)

effect = result['observations']['effect']
receptor = result['states']['R']
t = result['t']

print("Beta-Receptor Down-Regulation:")
print(f"  Initial receptor: {receptor[0]:.2f}")
print(f"  Day 7 receptor: {receptor[-1]:.2f}")
print(f"  First dose peak: {max(effect[:16]):.1f}")
print(f"  Day 7 dose peak: {max(effect[-16:]):.1f}")
```

---

## Receptor Up-Regulation Example

```python
import neopkpd

# Receptor up-regulation (e.g., chronic antagonist exposure)
doses = [{"time": i * 12.0, "amount": 30.0} for i in range(14)]

result = neopkpd.simulate_pkpd_receptor_regulation(
    cl=5.0, v=50.0, doses=doses,
    e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
    r_baseline=1.0,
    kreg=0.03,
    rmax=3.0,          # Max 3x receptor density
    kchange=0.01,
    direction="up",
    t0=0.0, t1=168.0,
    saveat=[i * 0.5 for i in range(337)]
)

effect = result['observations']['effect']
receptor = result['states']['R']
t = result['t']

print("Receptor Up-Regulation:")
print(f"  Initial receptor: {receptor[0]:.2f}")
print(f"  Day 7 receptor: {receptor[-1]:.2f}")
print(f"  Receptor increase: {(receptor[-1]/receptor[0] - 1) * 100:.1f}%")
```

---

## Recovery After Drug Discontinuation

```python
import neopkpd

# 1 week dosing, then 2 weeks recovery
doses = [{"time": i * 8.0, "amount": 50.0} for i in range(21)]  # 7 days only

result = neopkpd.simulate_pkpd_tolerance_counter_regulation(
    cl=5.0, v=50.0, doses=doses,
    e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
    kin_mod=0.1, kout_mod=0.02, alpha_feedback=1.0,
    t0=0.0, t1=504.0,  # 3 weeks total
    saveat=[i * 1.0 for i in range(505)]
)

effect = result['observations']['effect']
t = result['t']

print("Tolerance Recovery:")
print(f"  Day 7 (end of dosing): Moderator at maximum")
print(f"  Day 14: Partial recovery")
print(f"  Day 21: Near baseline")

# Simulate re-exposure after 2-week break
# (Would need new simulation with rechallenge dose)
```

---

## Comparing Tolerance Mechanisms

```python
import neopkpd

doses = [{"time": i * 8.0, "amount": 50.0} for i in range(21)]

# Counter-regulation
result_counter = neopkpd.simulate_pkpd_tolerance_counter_regulation(
    cl=5.0, v=50.0, doses=doses,
    e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
    kin_mod=0.1, kout_mod=0.05, alpha_feedback=1.0,
    t0=0.0, t1=168.0,
    saveat=[i * 0.5 for i in range(337)]
)

# Receptor down-regulation
result_receptor = neopkpd.simulate_pkpd_receptor_regulation(
    cl=5.0, v=50.0, doses=doses,
    e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
    r_baseline=1.0, kreg=0.05, rmax=2.0, kchange=0.02,
    direction="down",
    t0=0.0, t1=168.0,
    saveat=[i * 0.5 for i in range(337)]
)

counter_effect = result_counter['observations']['effect']
receptor_effect = result_receptor['observations']['effect']

print("Tolerance Mechanism Comparison:")
print("                  | Counter-Reg | Receptor")
print("-" * 50)
print(f"Day 1 peak effect | {max(counter_effect[:48]):11.1f} | {max(receptor_effect[:48]):.1f}")
print(f"Day 7 peak effect | {max(counter_effect[-48:]):11.1f} | {max(receptor_effect[-48:]):.1f}")
```

---

## Dose Escalation to Overcome Tolerance

```python
import neopkpd

# Escalating dose to maintain effect
base_dose = 50.0
escalation = 1.2  # 20% increase each day

doses = []
for day in range(7):
    daily_dose = base_dose * (escalation ** day)
    for dose_num in range(3):  # TID
        doses.append({"time": day * 24 + dose_num * 8.0, "amount": daily_dose})

result = neopkpd.simulate_pkpd_tolerance_counter_regulation(
    cl=5.0, v=50.0, doses=doses,
    e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
    kin_mod=0.1, kout_mod=0.05, alpha_feedback=1.0,
    t0=0.0, t1=168.0,
    saveat=[i * 0.5 for i in range(337)]
)

effect = result['observations']['effect']

print("Dose Escalation to Maintain Effect:")
for day in range(7):
    start = day * 48
    end = (day + 1) * 48
    peak = max(effect[start:end])
    dose = base_dose * (escalation ** day)
    print(f"  Day {day+1}: Dose = {dose:.0f} mg, Peak effect = {peak:.1f}")
```

---

## Equations Summary

### Counter-Regulation Model

| Quantity | Formula |
|----------|---------|
| Drug effect | $E_{max} \cdot C^\gamma / (EC_{50}^\gamma + C^\gamma)$ |
| Moderator rate | $k_{in,mod} \cdot E_{drug} - k_{out,mod} \cdot M$ |
| Net effect | $E_0 + E_{drug} - \alpha \cdot M$ |
| Moderator t1/2 | $\ln(2) / k_{out,mod}$ |

### Receptor Regulation Model

| Quantity | Formula |
|----------|---------|
| Down-regulation | $k_{reg}(R_0 - R) - k_{change} \cdot E_{drug} \cdot R$ |
| Up-regulation | $k_{reg}(R_0 - R) + k_{change} \cdot E_{drug} \cdot (R_{max} - R)$ |
| Net effect | $E_0 + R \cdot E_{drug}$ |

---

## See Also

- [Direct Emax](direct-emax.md) - Without tolerance
- [Indirect Response](indirect-response.md) - Turnover-based effects
- [Effect Compartment](effect-compartment.md) - Temporal delay
