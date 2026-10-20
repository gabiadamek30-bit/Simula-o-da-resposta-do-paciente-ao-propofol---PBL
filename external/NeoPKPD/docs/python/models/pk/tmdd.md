# Target-Mediated Drug Disposition (TMDD)

Advanced PK model for drugs that bind to their pharmacological target, forming drug-target complexes that affect both PK and PD behavior.

---

## Function Signature

```python
neopkpd.simulate_pk_tmdd_custom(
    kel: float,
    kon: float,
    koff: float,
    ksyn: float,
    kdeg: float,
    kint: float,
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
| `kel` | float | Drug elimination rate constant (1/h) |
| `kon` | float | Association rate constant (1/(concentration*h)) |
| `koff` | float | Dissociation rate constant (1/h) |
| `ksyn` | float | Receptor synthesis rate (concentration/h) |
| `kdeg` | float | Receptor degradation rate constant (1/h) |
| `kint` | float | Complex internalization rate constant (1/h) |
| `v` | float | Volume of distribution (L) |
| `doses` | list[dict] | Dose events |

### Derived Parameters

- **KD (dissociation constant)**: $K_D = k_{off} / k_{on}$
- **Receptor baseline**: $R_0 = k_{syn} / k_{deg}$

---

## Model Equations

Three-state system:
$$\frac{dL}{dt} = -k_{el} \cdot L - k_{on} \cdot L \cdot R + k_{off} \cdot RL$$

$$\frac{dR}{dt} = k_{syn} - k_{deg} \cdot R - k_{on} \cdot L \cdot R + k_{off} \cdot RL$$

$$\frac{dRL}{dt} = k_{on} \cdot L \cdot R - k_{off} \cdot RL - k_{int} \cdot RL$$

Where:
- L = Free drug (ligand) concentration
- R = Free receptor concentration
- RL = Drug-receptor complex concentration

---

## Returns

```python
{
    "t": [0.0, 1.0, ...],
    "states": {
        "L_free": [...],      # Free drug
        "R_free": [...],      # Free receptor
        "RL_complex": [...]   # Drug-receptor complex
    },
    "observations": {
        "conc": [...]         # Free drug concentration
    },
    "metadata": {...}
}
```

---

## Basic Example

```python
import neopkpd

result = neopkpd.simulate_pk_tmdd_custom(
    kel=0.1,      # Drug elimination (1/h)
    kon=0.01,     # Association rate
    koff=0.001,   # Dissociation rate
    ksyn=1.0,     # Receptor synthesis
    kdeg=0.1,     # Receptor degradation
    kint=0.05,    # Complex internalization
    v=50.0,       # Volume (L)
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=72.0,
    saveat=list(range(73))
)

conc = result['observations']['conc']
t = result['t']

print(f"Initial free drug: {conc[0]:.2f} mg/L")
print(f"Free drug at 24h: {conc[24]:.3f} mg/L")
```

---

## Non-Linear PK Behavior

TMDD causes characteristic non-linear PK:

```python
import neopkpd

doses_list = [50.0, 100.0, 200.0, 500.0, 1000.0]

print("Dose (mg) | Cmax (mg/L) | AUC ratio")
print("-" * 45)

auc_ref = None
for dose in doses_list:
    result = neopkpd.simulate_pk_tmdd_custom(
        kel=0.1, kon=0.01, koff=0.001,
        ksyn=1.0, kdeg=0.1, kint=0.05, v=50.0,
        doses=[{"time": 0.0, "amount": dose}],
        t0=0.0, t1=96.0,
        saveat=[i * 0.5 for i in range(193)]
    )

    conc = result['observations']['conc']
    t = result['t']

    cmax = max(conc)

    # Calculate AUC (trapezoidal)
    auc = sum(0.5 * (conc[i] + conc[i+1]) * (t[i+1] - t[i])
              for i in range(len(conc)-1))

    if auc_ref is None:
        auc_ref = auc

    print(f"{dose:9.0f} | {cmax:11.2f} | {auc/auc_ref:9.2f}")
```

**Expected:** More-than-proportional increase in AUC at higher doses due to receptor saturation.

---

## Target-Mediated Clearance

```python
import neopkpd

# At low concentrations: high clearance (receptor-mediated)
# At high concentrations: low clearance (receptor saturated)

result = neopkpd.simulate_pk_tmdd_custom(
    kel=0.05,      # Low linear elimination
    kon=0.1,       # Fast binding
    koff=0.01,     # Slow unbinding (high affinity)
    ksyn=0.5,      # Receptor turnover
    kdeg=0.05,
    kint=0.2,      # Fast internalization
    v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=168.0,
    saveat=[i * 0.5 for i in range(337)]
)

conc = result['observations']['conc']
t = result['t']

# Find apparent terminal half-life at different phases
# Early: rapid decline (receptor-mediated)
# Late: slower decline (linear elimination dominates)
print("Phase analysis:")
print(f"  Conc at 4h: {conc[8]:.3f} mg/L")
print(f"  Conc at 24h: {conc[48]:.4f} mg/L")
print(f"  Conc at 72h: {conc[144]:.5f} mg/L")
```

---

## Clinical Example: Monoclonal Antibody

```python
import neopkpd

# Typical mAb TMDD parameters
result = neopkpd.simulate_pk_tmdd_custom(
    kel=0.01,        # Slow linear elimination (typical mAb)
    kon=0.1,         # Fast target binding
    koff=0.001,      # Very slow unbinding (high affinity)
    ksyn=0.1,        # Target synthesis
    kdeg=0.05,       # Target degradation
    kint=0.02,       # Complex internalization
    v=3.0,           # Central volume (L, typical mAb)
    doses=[{"time": 0.0, "amount": 100.0}],  # 100 mg dose
    t0=0.0, t1=672.0,  # 4 weeks
    saveat=[i * 6.0 for i in range(113)]  # Every 6 hours
)

# Access states
l_free = result['states']['L_free']
r_free = result['states']['R_free']
rl_complex = result['states']['RL_complex']
t = result['t']

print("Monoclonal Antibody TMDD Profile:")
print(f"  Initial free drug: {l_free[0]:.2f} mg/L")
print(f"  Initial free receptor: {r_free[0]:.2f}")
print(f"  Free receptor at 24h: {r_free[4]:.3f} (suppressed)")
print(f"  Complex at 24h: {rl_complex[4]:.2f}")
```

---

## Receptor Dynamics

```python
import neopkpd

result = neopkpd.simulate_pk_tmdd_custom(
    kel=0.1, kon=0.05, koff=0.005,
    ksyn=1.0, kdeg=0.1, kint=0.1, v=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0, t1=168.0,
    saveat=[i * 1.0 for i in range(169)]
)

r_free = result['states']['R_free']
t = result['t']

# Baseline receptor level
r_baseline = 1.0 / 0.1  # ksyn / kdeg = 10

print("Receptor Dynamics:")
print(f"  Baseline R: {r_baseline:.1f}")
print(f"  R at 1h: {r_free[1]:.2f} (depleted by drug)")
print(f"  R at 24h: {r_free[24]:.2f}")
print(f"  R at 168h: {r_free[168]:.2f} (recovering)")
```

---

## Multiple Dosing

```python
import neopkpd

# Weekly dosing of mAb
doses = [{"time": i * 168.0, "amount": 200.0} for i in range(4)]

result = neopkpd.simulate_pk_tmdd_custom(
    kel=0.01, kon=0.1, koff=0.001,
    ksyn=0.1, kdeg=0.05, kint=0.02, v=3.0,
    doses=doses,
    t0=0.0, t1=672.0,
    saveat=[i * 6.0 for i in range(113)]
)

conc = result['observations']['conc']
r_free = result['states']['R_free']

print("Weekly Dosing Profile:")
print(f"  Peak after dose 1: {max(conc[:28]):.2f} mg/L")
print(f"  Trough before dose 2: {conc[27]:.3f} mg/L")
print(f"  Receptor suppression at week 4: {r_free[-1] / (0.1/0.05) * 100:.1f}% of baseline")
```

---

## Dose Selection Considerations

```python
import neopkpd

# Find dose needed for >90% receptor occupancy
r_baseline = 1.0 / 0.1  # ksyn / kdeg

for dose in [10, 50, 100, 200, 500]:
    result = neopkpd.simulate_pk_tmdd_custom(
        kel=0.1, kon=0.05, koff=0.005,
        ksyn=1.0, kdeg=0.1, kint=0.1, v=50.0,
        doses=[{"time": 0.0, "amount": float(dose)}],
        t0=0.0, t1=24.0,
        saveat=[i * 0.5 for i in range(49)]
    )

    r_free = result['states']['R_free']
    min_r = min(r_free)
    max_occupancy = (r_baseline - min_r) / r_baseline * 100

    print(f"Dose {dose:3d} mg: Max receptor occupancy = {max_occupancy:.1f}%")
```

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| KD | $k_{off} / k_{on}$ |
| Receptor baseline | $R_0 = k_{syn} / k_{deg}$ |
| Free drug rate | $-k_{el}L - k_{on}LR + k_{off}RL$ |
| Complex rate | $k_{on}LR - k_{off}RL - k_{int}RL$ |
| Total drug | $L + RL$ |

---

## See Also

- [Michaelis-Menten](michaelis-menten.md) - Saturable elimination
- [Two-Compartment IV](twocomp-iv.md) - Distribution kinetics
- [Population Simulation](../../population/index.md) - Adding variability
