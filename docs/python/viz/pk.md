# PK Plots

Concentration-time profile visualization functions.

---

## Overview

PK plots display drug concentration over time, the fundamental visualization in pharmacokinetics.

```python
from neopkpd import viz

# Basic concentration-time plot
fig = viz.plot_conc_time(result, title="PK Profile")
```

---

## Functions

### plot_conc_time

Single subject concentration-time profile:

```python
def plot_conc_time(
    result: SimResult | dict,
    *,
    log_scale: bool = False,
    title: str | None = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: tuple = (10, 6),
    color: str | None = None,
    marker: str | None = None,
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

**Usage:**

```python
result = neopkpd.simulate_pk_oral(
    ka=1.5, cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5
)

# Linear scale
fig = viz.plot_conc_time(result, title="Oral PK")

# Semi-log scale
fig = viz.plot_conc_time(result, log_scale=True)

# With markers
fig = viz.plot_conc_time(result, marker="o")
```

### plot_multi_conc_time

Multiple profiles overlay:

```python
def plot_multi_conc_time(
    results: list[SimResult | dict],
    labels: list[str] | None = None,
    *,
    log_scale: bool = False,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Compare different doses
results = []
labels = []
for dose in [50, 100, 200]:
    r = neopkpd.simulate_pk_oral(
        ka=1.5, cl=5.0, v=50.0,
        doses=[{"time": 0.0, "amount": float(dose)}],
        t0=0.0, t1=24.0, saveat=0.5
    )
    results.append(r)
    labels.append(f"{dose} mg")

fig = viz.plot_multi_conc_time(results, labels=labels, title="Dose Comparison")
```

### plot_spaghetti

Population spaghetti plot:

```python
def plot_spaghetti(
    pop_result: PopulationResult | dict,
    *,
    n_subjects: int | None = None,
    alpha: float = 0.3,
    show_mean: bool = True,
    log_scale: bool = False,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
pop_result = neopkpd.simulate_population_oral(
    ka=1.5, cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100, omegas={"Ka": 0.16, "CL": 0.09, "V": 0.04},
    seed=42
)

# All subjects
fig = viz.plot_spaghetti(pop_result, title="Population Profiles")

# First 20 subjects only
fig = viz.plot_spaghetti(pop_result, n_subjects=20, alpha=0.5)
```

### plot_mean_ribbon

Mean with confidence ribbon:

```python
def plot_mean_ribbon(
    pop_result: PopulationResult | dict,
    *,
    ci_levels: list[float] = [0.05, 0.95],
    show_median: bool = True,
    show_individual: bool = False,
    log_scale: bool = False,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# 90% prediction interval
fig = viz.plot_mean_ribbon(pop_result, ci_levels=[0.05, 0.95])

# 50% and 90% intervals
fig = viz.plot_mean_ribbon(
    pop_result,
    ci_levels=[0.25, 0.75],  # Inner ribbon
    show_median=True
)
```

### plot_individual_fits

Grid of individual subject fits:

```python
def plot_individual_fits(
    pop_result: PopulationResult | dict,
    observed: dict | None = None,
    *,
    n_subjects: int = 9,
    n_cols: int = 3,
    figsize: tuple | None = None,
    log_scale: bool = False,
    title: str | None = None,
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Grid of 9 subjects (3x3)
fig = viz.plot_individual_fits(pop_result, n_subjects=9, n_cols=3)

# With observed data overlay
fig = viz.plot_individual_fits(pop_result, observed=observed_data)
```

---

## Complete Example

```python
import neopkpd
from neopkpd import viz

neopkpd.init_julia()
viz.set_backend("matplotlib")

# Single simulation
result = neopkpd.simulate_pk_oral(
    ka=1.5, cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5
)

# Basic plot
fig = viz.plot_conc_time(result, title="One-Compartment Oral PK")
fig.savefig("pk_single.png", dpi=300)

# Population simulation
pop_result = neopkpd.simulate_population_oral(
    ka=1.5, cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omegas={"Ka": 0.16, "CL": 0.09, "V": 0.04},
    seed=42
)

# Spaghetti plot
fig = viz.plot_spaghetti(pop_result, n_subjects=50, alpha=0.2)
fig.savefig("pk_spaghetti.png", dpi=300)

# Mean with ribbon
fig = viz.plot_mean_ribbon(pop_result, ci_levels=[0.05, 0.95])
fig.savefig("pk_ribbon.png", dpi=300)
```

---

## See Also

- [Backends & Themes](backends.md) - Styling options
- [Population Plots](population.md) - More population visualizations
- [VPC Plots](vpc.md) - Model validation plots
