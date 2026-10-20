# Models Reference

Python simulation functions for all PK and PD models supported by NeoPKPD.

---

## Model Categories

<div class="grid cards" markdown>

-   :material-cube-outline:{ .lg .middle } **Pharmacokinetic Models**

    ---

    Compartmental PK simulation functions

    [:octicons-arrow-right-24: PK Functions](#pharmacokinetic-functions)

-   :material-chart-line:{ .lg .middle } **Pharmacodynamic Models**

    ---

    Effect model simulation functions

    [:octicons-arrow-right-24: PD Functions](#pharmacodynamic-functions)

</div>

---

## Pharmacokinetic Functions

| Function | Model | Parameters |
|----------|-------|------------|
| [`simulate_pk_iv_bolus`](pk/iv-bolus.md) | One-comp IV | cl, v |
| [`simulate_pk_oral_first_order`](pk/oral.md) | One-comp oral | ka, cl, v |
| [`simulate_pk_twocomp_iv_bolus`](pk/twocomp-iv.md) | Two-comp IV | cl, v1, q, v2 |
| [`simulate_pk_twocomp_oral`](pk/twocomp-oral.md) | Two-comp oral | ka, cl, v1, q, v2 |
| [`simulate_pk_threecomp_iv_bolus`](pk/threecomp.md) | Three-comp IV | cl, v1, q2, v2, q3, v3 |
| [`simulate_pk_transit_absorption`](pk/transit.md) | Transit | n_transit, ktr, ka, cl, v |
| [`simulate_pk_michaelis_menten`](pk/michaelis-menten.md) | MM elimination | vmax, km, v |

---

## Pharmacodynamic Functions

| Function | Model | PD Parameters |
|----------|-------|---------------|
| [`simulate_pkpd_direct_emax`](pd/direct-emax.md) | Direct Emax | e0, emax, ec50 |
| [`simulate_pkpd_sigmoid_emax`](pd/sigmoid-emax.md) | Sigmoid Emax | e0, emax, ec50, gamma |
| [`simulate_pkpd_biophase_equilibration`](pd/effect-compartment.md) | Effect compartment | ke0, e0, emax, ec50 |
| [`simulate_pkpd_indirect_response`](pd/indirect-response.md) | Indirect response | kin, kout, ic50, imax |

---

## Common Parameters

All simulation functions share these parameters:

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `doses` | `list[dict]` | Yes | List of dose events |
| `t0` | `float` | Yes | Simulation start time |
| `t1` | `float` | Yes | Simulation end time |
| `saveat` | `list[float]` | Yes | Output time points |

### Dose Event Format

```python
doses = [
    {"time": 0.0, "amount": 100.0},              # Bolus
    {"time": 12.0, "amount": 100.0},             # Bolus
    {"time": 0.0, "amount": 100.0, "duration": 1.0}  # Infusion
]
```

---

## Quick Examples

### One-Compartment IV Bolus

```python
import neopkpd

neopkpd.init_julia()

result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)]
)

print("Concentrations:", result["observations"]["conc"][:5])
```

### Two-Compartment with Distribution

```python
result = neopkpd.simulate_pk_twocomp_iv_bolus(
    cl=10.0,
    v1=20.0,
    q=15.0,
    v2=50.0,
    doses=[{"time": 0.0, "amount": 500.0}],
    t0=0.0,
    t1=48.0,
    saveat=[t * 0.5 for t in range(97)]
)
```

### Direct Emax PD

```python
result = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    e0=0.0,
    emax=100.0,
    ec50=2.0,
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)]
)

print("Effect:", result["observations"]["effect"][:5])
```

---

## Return Structure

All functions return a dictionary:

```python
{
    "t": [0.0, 1.0, ...],              # Time points
    "states": {
        "A_central": [100.0, 90.5, ...]  # State variables
    },
    "observations": {
        "conc": [2.0, 1.81, ...],        # Concentrations
        "effect": [...]                   # Effects (PKPD only)
    },
    "metadata": {
        "model": "OneCompIVBolus",
        "version": "0.1.0"
    }
}
```

---

## Next Steps

- [IV Bolus Details](pk/iv-bolus.md)
- [Population Simulation](../population/index.md)
- [Visualization](../viz/index.md)
