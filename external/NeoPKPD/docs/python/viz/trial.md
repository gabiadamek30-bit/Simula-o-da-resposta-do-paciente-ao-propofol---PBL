# Trial Plots

Clinical trial visualization functions.

---

## Overview

Trial plots visualize study designs, power analyses, and endpoint distributions.

```python
from neopkpd import viz

fig = viz.plot_power_curve(power_results, title="Power Analysis")
```

---

## Functions

### plot_power_curve

Power vs sample size curve:

```python
def plot_power_curve(
    power_results: dict | list[dict],
    *,
    target_power: float = 0.8,
    show_target: bool = True,
    title: str | None = None,
    xlabel: str = "Sample Size (per arm)",
    ylabel: str = "Power",
    figsize: tuple = (10, 6),
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

**Usage:**

```python
# Power analysis results
power_results = {
    "n": [10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
    "power": [0.25, 0.45, 0.62, 0.75, 0.83, 0.89, 0.93, 0.95, 0.97, 0.98]
}

fig = viz.plot_power_curve(
    power_results,
    target_power=0.8,
    title="Sample Size Calculation"
)
```

### plot_trial_tornado

Trial design sensitivity:

```python
def plot_trial_tornado(
    sensitivity_results: list[dict],
    *,
    baseline: float = 0.8,
    title: str | None = None,
    xlabel: str = "Power",
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
sensitivity_results = [
    {"parameter": "Effect Size", "low": 0.65, "high": 0.92},
    {"parameter": "Variability", "low": 0.70, "high": 0.88},
    {"parameter": "Dropout Rate", "low": 0.75, "high": 0.82},
]

fig = viz.plot_trial_tornado(
    sensitivity_results,
    baseline=0.80,
    title="Power Sensitivity"
)
```

### plot_kaplan_meier

Survival curves with confidence intervals:

```python
def plot_kaplan_meier(
    survival_data: dict | list[dict],
    *,
    show_ci: bool = True,
    show_censored: bool = True,
    title: str | None = None,
    xlabel: str = "Time",
    ylabel: str = "Survival Probability",
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
survival_data = [
    {
        "label": "Treatment",
        "times": [1, 3, 5, 7, 10, 12, 15, 18, 20, 24],
        "survival": [0.98, 0.95, 0.90, 0.85, 0.78, 0.72, 0.65, 0.60, 0.55, 0.50],
        "ci_lower": [0.95, 0.90, 0.82, 0.75, 0.65, 0.58, 0.50, 0.45, 0.40, 0.35],
        "ci_upper": [1.00, 0.98, 0.96, 0.92, 0.88, 0.83, 0.78, 0.73, 0.68, 0.63],
    },
    {
        "label": "Control",
        "times": [1, 3, 5, 7, 10, 12, 15, 18, 20, 24],
        "survival": [0.96, 0.88, 0.78, 0.68, 0.55, 0.45, 0.38, 0.32, 0.28, 0.25],
    }
]

fig = viz.plot_kaplan_meier(survival_data, title="Overall Survival")
```

### plot_endpoint_distribution

Endpoint histogram by group:

```python
def plot_endpoint_distribution(
    endpoint_data: dict[str, list[float]],
    *,
    show_stats: bool = True,
    title: str | None = None,
    xlabel: str = "Endpoint",
    ylabel: str = "Frequency",
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
endpoint_data = {
    "Treatment": [12.5, 15.2, 11.8, 14.0, 13.5, 16.1, 12.0, 14.8],
    "Placebo": [8.2, 9.5, 7.8, 10.1, 8.9, 9.2, 8.0, 9.8],
}

fig = viz.plot_endpoint_distribution(
    endpoint_data,
    show_stats=True,
    title="Primary Endpoint Distribution"
)
```

---

## Complete Example

```python
import neopkpd
from neopkpd import viz
import numpy as np

neopkpd.init_julia()
viz.set_backend("matplotlib")

# Power analysis
power_results = neopkpd.calculate_power(
    effect_size=0.5,
    variability=0.8,
    alpha=0.05,
    n_range=range(10, 101, 10)
)

fig = viz.plot_power_curve(
    power_results,
    target_power=0.8,
    title="Sample Size for 80% Power"
)
fig.savefig("power_curve.png", dpi=300)

# Trial sensitivity
sensitivity = [
    {"parameter": "Effect Size ±20%", "low": 0.62, "high": 0.91},
    {"parameter": "CV ±25%", "low": 0.68, "high": 0.88},
    {"parameter": "Dropout 10% vs 20%", "low": 0.75, "high": 0.82},
]
fig = viz.plot_trial_tornado(sensitivity, baseline=0.80)
fig.savefig("trial_tornado.png", dpi=300)

# Simulated endpoint distribution
np.random.seed(42)
endpoints = {
    "Active": np.random.normal(15, 3, 50).tolist(),
    "Placebo": np.random.normal(12, 3, 50).tolist(),
}
fig = viz.plot_endpoint_distribution(
    endpoints,
    title="Simulated AUC Response"
)
fig.savefig("endpoint_dist.png", dpi=300)
```

---

## See Also

- [Clinical Trials](../trial/index.md) - Trial design module
- [Power Analysis](../trial/power.md) - Sample size calculation
- [Population Plots](population.md) - Population visualizations
