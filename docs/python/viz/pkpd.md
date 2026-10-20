# PKPD Plots

Pharmacokinetic-pharmacodynamic relationship visualization.

---

## Overview

PKPD plots visualize the relationship between drug concentration and effect.

```python
from neopkpd import viz

fig = viz.plot_effect_conc(result, title="Concentration-Effect")
```

---

## Functions

### plot_effect_conc

Effect vs concentration relationship:

```python
def plot_effect_conc(
    result: PKPDResult | dict,
    *,
    show_emax_fit: bool = True,
    show_ec50: bool = True,
    title: str | None = None,
    xlabel: str = "Concentration",
    ylabel: str = "Effect",
    figsize: tuple = (10, 6),
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

**Usage:**

```python
result = neopkpd.simulate_pkpd_emax(
    cl=5.0, v=50.0,
    emax=100.0, ec50=5.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5
)

fig = viz.plot_effect_conc(result, title="Emax Model")
```

### plot_hysteresis

Hysteresis loop (effect compartment):

```python
def plot_hysteresis(
    result: PKPDResult | dict,
    *,
    show_arrows: bool = True,
    title: str | None = None,
    xlabel: str = "Concentration",
    ylabel: str = "Effect",
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Effect compartment model with hysteresis
result = neopkpd.simulate_pkpd_effect_compartment(
    cl=5.0, v=50.0, ke0=0.5,
    emax=100.0, ec50=5.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5
)

fig = viz.plot_hysteresis(result, show_arrows=True, title="Hysteresis Loop")
```

### plot_dose_response

Dose-response curve:

```python
def plot_dose_response(
    results: list[PKPDResult | dict],
    doses: list[float],
    metric: str = "emax",
    *,
    show_fit: bool = True,
    title: str | None = None,
    xlabel: str = "Dose",
    ylabel: str = "Effect",
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Simulate multiple doses
doses = [10, 25, 50, 100, 200, 400]
results = []
for dose in doses:
    r = neopkpd.simulate_pkpd_emax(
        cl=5.0, v=50.0, emax=100.0, ec50=5.0,
        doses=[{"time": 0.0, "amount": float(dose)}],
        t0=0.0, t1=24.0, saveat=0.5
    )
    results.append(r)

fig = viz.plot_dose_response(results, doses, title="Dose-Response")
```

---

## Complete Example

```python
import neopkpd
from neopkpd import viz

neopkpd.init_julia()
viz.set_backend("matplotlib")

# Emax model simulation
result = neopkpd.simulate_pkpd_emax(
    cl=5.0, v=50.0,
    emax=100.0, ec50=5.0,
    e0=10.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.25
)

# Concentration-effect plot
fig = viz.plot_effect_conc(
    result,
    show_ec50=True,
    title="Direct Emax Model"
)
fig.savefig("effect_conc.png", dpi=300)

# Effect compartment with hysteresis
result_ec = neopkpd.simulate_pkpd_effect_compartment(
    cl=5.0, v=50.0, ke0=0.3,
    emax=100.0, ec50=5.0, e0=10.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.25
)

fig = viz.plot_hysteresis(result_ec, title="Counter-Clockwise Hysteresis")
fig.savefig("hysteresis.png", dpi=300)
```

---

## See Also

- [PD Models](../models/pd/direct-emax.md) - PD model documentation
- [PK Plots](pk.md) - Concentration-time plots
- [Population Plots](population.md) - Population PKPD
