# Population Plots

Population modeling visualization functions.

---

## Overview

Population plots visualize inter-individual variability, parameter distributions, and model diagnostics.

```python
from neopkpd import viz

fig = viz.plot_parameter_distributions(pop_result)
```

---

## Functions

### plot_vpc

Basic Visual Predictive Check:

```python
def plot_vpc(
    vpc_result: VPCResult | dict,
    *,
    log_scale: bool = False,
    show_ci: bool = True,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

See [VPC Plots](vpc.md) for detailed VPC visualization.

### plot_parameter_distributions

Parameter distribution histograms:

```python
def plot_parameter_distributions(
    pop_result: PopulationResult | dict,
    params: list[str] | None = None,
    *,
    show_typical: bool = True,
    n_cols: int = 3,
    title: str | None = None,
    figsize: tuple | None = None,
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

**Usage:**

```python
pop_result = neopkpd.simulate_population_oral(
    ka=1.5, cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omegas={"Ka": 0.16, "CL": 0.09, "V": 0.04},
    seed=42
)

fig = viz.plot_parameter_distributions(
    pop_result,
    show_typical=True,
    title="Individual Parameter Distributions"
)
```

### plot_forest

Forest plot of parameter effects:

```python
def plot_forest(
    forest_data: list[dict],
    *,
    reference_line: float = 1.0,
    show_ci: bool = True,
    ci_level: float = 0.95,
    title: str | None = None,
    xlabel: str = "Effect",
    figsize: tuple = (10, 8),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
forest_data = [
    {"label": "Age (per 10 yr)", "estimate": 0.95, "lower": 0.88, "upper": 1.02},
    {"label": "Weight (per 10 kg)", "estimate": 1.12, "lower": 1.05, "upper": 1.20},
    {"label": "Sex (Female)", "estimate": 0.85, "lower": 0.75, "upper": 0.96},
    {"label": "Renal Impairment", "estimate": 0.72, "lower": 0.58, "upper": 0.89},
]

fig = viz.plot_forest(
    forest_data,
    reference_line=1.0,
    title="Covariate Effects on Clearance"
)
```

### plot_boxplot

Box plot comparison:

```python
def plot_boxplot(
    data: dict[str, list[float]],
    *,
    show_points: bool = False,
    title: str | None = None,
    xlabel: str = "Group",
    ylabel: str = "Value",
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Compare AUC by group
data = {
    "Low Dose": [85, 92, 78, 95, 88],
    "Medium Dose": [150, 165, 142, 158, 170],
    "High Dose": [280, 310, 265, 295, 320],
}

fig = viz.plot_boxplot(data, ylabel="AUC (mg*hr/L)")
```

### plot_goodness_of_fit

4-panel GOF diagnostic:

```python
def plot_goodness_of_fit(
    est_result: EstimationResult | dict,
    *,
    log_scale: bool = False,
    title: str | None = None,
    figsize: tuple = (12, 10),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# 4-panel GOF: DV vs PRED, DV vs IPRED, CWRES vs TIME, CWRES vs PRED
fig = viz.plot_goodness_of_fit(
    est_result,
    log_scale=False,
    title="Goodness of Fit"
)
```

### plot_estimation_summary

Summary dashboard:

```python
def plot_estimation_summary(
    est_result: EstimationResult | dict,
    *,
    title: str | None = None,
    figsize: tuple = (16, 12),
    backend: str | None = None
) -> Figure:
```

### plot_sensitivity

Parameter sensitivity:

```python
def plot_sensitivity(
    pop_result: PopulationResult | dict,
    param: str,
    metric: str = "auc",
    *,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

### plot_sensitivity_tornado

Sensitivity tornado plot:

```python
def plot_sensitivity_tornado(
    sensitivity_results: list[dict],
    *,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

---

## Complete Example

```python
import neopkpd
from neopkpd import viz

neopkpd.init_julia()
viz.set_backend("matplotlib")

# Population simulation
pop_result = neopkpd.simulate_population_oral(
    ka=1.5, cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omegas={"Ka": 0.16, "CL": 0.09, "V": 0.04},
    seed=42
)

# Parameter distributions
fig = viz.plot_parameter_distributions(
    pop_result,
    show_typical=True,
    n_cols=3,
    title="Individual Parameters"
)
fig.savefig("param_dist.png", dpi=300)

# Forest plot (example data)
forest_data = [
    {"label": "Weight", "estimate": 1.15, "lower": 1.08, "upper": 1.23},
    {"label": "Age", "estimate": 0.92, "lower": 0.85, "upper": 0.99},
    {"label": "Sex", "estimate": 0.88, "lower": 0.78, "upper": 0.99},
]
fig = viz.plot_forest(forest_data, title="Covariate Effects on CL")
fig.savefig("forest.png", dpi=300)

# Box plot
import numpy as np
cmax_by_group = {
    "Young": np.random.lognormal(2.0, 0.3, 30).tolist(),
    "Elderly": np.random.lognormal(2.2, 0.35, 25).tolist(),
}
fig = viz.plot_boxplot(cmax_by_group, ylabel="Cmax (mg/L)")
fig.savefig("boxplot.png", dpi=300)
```

---

## See Also

- [VPC Plots](vpc.md) - Visual predictive checks
- [PK Plots](pk.md) - Concentration-time profiles
- [Estimation Diagnostics](estimation.md) - Fitting diagnostics
