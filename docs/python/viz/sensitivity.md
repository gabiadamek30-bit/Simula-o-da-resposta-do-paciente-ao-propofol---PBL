# Sensitivity Plots

Sensitivity analysis visualization functions.

---

## Overview

Sensitivity plots visualize how model outputs change with parameter variations.

```python
from neopkpd import viz

fig = viz.plot_tornado(sensitivity_results, title="Parameter Sensitivity")
```

---

## Functions

### plot_tornado

Tornado diagram of parameter sensitivity:

```python
def plot_tornado(
    sensitivity_results: list[dict],
    *,
    baseline_value: float = 0.0,
    sort_by: str = "range",
    title: str | None = None,
    xlabel: str = "Change from Baseline",
    figsize: tuple = (10, 6),
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

**Usage:**

```python
sensitivity_results = [
    {"parameter": "CL", "low": -0.30, "high": 0.25},
    {"parameter": "V", "low": -0.15, "high": 0.18},
    {"parameter": "Ka", "low": -0.40, "high": 0.35},
]

fig = viz.plot_tornado(
    sensitivity_results,
    baseline_value=0.0,
    title="AUC Sensitivity to Parameters"
)
```

### plot_spider

Spider/radar plot of parameter effects:

```python
def plot_spider(
    sensitivity_results: dict[str, list[float]],
    parameter_values: list[float],
    *,
    title: str | None = None,
    figsize: tuple = (10, 10),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Metric values at different parameter multipliers
sensitivity_results = {
    "CL": [1.5, 1.2, 1.0, 0.8, 0.6],
    "V": [0.9, 0.95, 1.0, 1.05, 1.1],
    "Ka": [1.3, 1.1, 1.0, 0.9, 0.7],
}
param_values = [0.5, 0.75, 1.0, 1.25, 1.5]

fig = viz.plot_spider(sensitivity_results, param_values)
```

### plot_sensitivity_heatmap

Parameter-metric sensitivity matrix:

```python
def plot_sensitivity_heatmap(
    sensitivity_matrix: np.ndarray | list[list[float]],
    parameters: list[str],
    metrics: list[str],
    *,
    annotate: bool = True,
    title: str | None = None,
    figsize: tuple = (10, 8),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
import numpy as np

# Sensitivity matrix: parameters x metrics
matrix = np.array([
    [0.8, -0.3, 0.1],   # CL effect on [AUC, Cmax, Tmax]
    [-0.2, 0.5, 0.0],   # V effect
    [0.1, 0.4, -0.6],   # Ka effect
])

fig = viz.plot_sensitivity_heatmap(
    matrix,
    parameters=["CL", "V", "Ka"],
    metrics=["AUC", "Cmax", "Tmax"],
    title="Parameter-Metric Sensitivity"
)
```

### plot_waterfall

Waterfall chart of ranked sensitivities:

```python
def plot_waterfall(
    sensitivities: dict[str, float],
    *,
    title: str | None = None,
    xlabel: str = "Sensitivity",
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
sensitivities = {
    "CL": -0.85,
    "Ka": 0.72,
    "V": -0.45,
    "F": 0.38,
    "Ke0": 0.15,
}

fig = viz.plot_waterfall(sensitivities, title="Ranked Sensitivities")
```

### plot_one_at_a_time

OFAT sensitivity curves:

```python
def plot_one_at_a_time(
    ofat_results: dict[str, tuple[list[float], list[float]]],
    *,
    baseline_value: float = 1.0,
    title: str | None = None,
    xlabel: str = "Parameter Multiplier",
    ylabel: str = "Metric Value",
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# OFAT results: param -> (multipliers, metric_values)
ofat_results = {
    "CL": ([0.5, 0.75, 1.0, 1.25, 1.5], [150, 120, 100, 85, 72]),
    "V": ([0.5, 0.75, 1.0, 1.25, 1.5], [98, 99, 100, 101, 102]),
    "Ka": ([0.5, 0.75, 1.0, 1.25, 1.5], [85, 92, 100, 108, 115]),
}

fig = viz.plot_one_at_a_time(ofat_results, baseline_value=100)
```

---

## Complete Example

```python
import neopkpd
from neopkpd import viz
import numpy as np

neopkpd.init_julia()
viz.set_backend("matplotlib")

# Run sensitivity analysis
sensitivity = neopkpd.run_sensitivity(
    model_spec=model,
    parameters=["CL", "V", "Ka"],
    ranges={"CL": (0.5, 1.5), "V": (0.5, 1.5), "Ka": (0.5, 1.5)},
    metric="auc"
)

# Tornado plot
fig = viz.plot_tornado(
    sensitivity["tornado"],
    title="AUC Sensitivity"
)
fig.savefig("tornado.png", dpi=300)

# Heatmap
fig = viz.plot_sensitivity_heatmap(
    sensitivity["matrix"],
    parameters=["CL", "V", "Ka"],
    metrics=["AUC", "Cmax", "t1/2"],
    title="Sensitivity Matrix"
)
fig.savefig("heatmap.png", dpi=300)

# Waterfall
fig = viz.plot_waterfall(
    sensitivity["ranked"],
    title="Ranked Parameter Sensitivities"
)
fig.savefig("waterfall.png", dpi=300)
```

---

## See Also

- [Population Plots](population.md) - Population sensitivity
- [Estimation Diagnostics](estimation.md) - Parameter uncertainty
- [Backends & Themes](backends.md) - Styling options
