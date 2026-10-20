# Bootstrap Plots

Bootstrap analysis visualization for parameter uncertainty.

---

## Overview

Bootstrap plots visualize parameter distributions and confidence intervals from resampling analysis.

```python
from neopkpd import viz

fig = viz.plot_bootstrap_distributions(bootstrap_result)
```

---

## Functions

### plot_bootstrap_distributions

Histograms of bootstrap parameter estimates:

```python
def plot_bootstrap_distributions(
    bootstrap_result: BootstrapResult | dict,
    params: list[str] | None = None,
    *,
    show_original: bool = True,
    show_ci: bool = True,
    ci_level: float = 0.95,
    n_cols: int = 3,
    title: str | None = None,
    figsize: tuple | None = None,
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

**Usage:**

```python
# Bootstrap parameter distributions
fig = viz.plot_bootstrap_distributions(
    bootstrap_result,
    show_original=True,
    show_ci=True,
    ci_level=0.95,
    title="Parameter Distributions (1000 bootstraps)"
)
```

### plot_bootstrap_ci

Confidence interval comparison:

```python
def plot_bootstrap_ci(
    bootstrap_result: BootstrapResult | dict,
    params: list[str] | None = None,
    *,
    methods: list[str] = ["percentile", "bca"],
    ci_level: float = 0.95,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Compare CI methods
fig = viz.plot_bootstrap_ci(
    bootstrap_result,
    methods=["percentile", "bca"],
    ci_level=0.95,
    title="Bootstrap CI Comparison"
)
```

### plot_bootstrap_stability

Parameter stability over bootstrap runs:

```python
def plot_bootstrap_stability(
    bootstrap_result: BootstrapResult | dict,
    params: list[str] | None = None,
    *,
    show_cumulative_mean: bool = True,
    title: str | None = None,
    figsize: tuple = (12, 8),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Stability assessment
fig = viz.plot_bootstrap_stability(
    bootstrap_result,
    show_cumulative_mean=True,
    title="Bootstrap Stability"
)
```

### plot_bootstrap_correlation

Inter-parameter correlations:

```python
def plot_bootstrap_correlation(
    bootstrap_result: BootstrapResult | dict,
    *,
    annotate: bool = True,
    title: str | None = None,
    figsize: tuple = (10, 10),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Correlation matrix from bootstrap
fig = viz.plot_bootstrap_correlation(
    bootstrap_result,
    annotate=True,
    title="Parameter Correlations"
)
```

---

## Complete Example

```python
import neopkpd
from neopkpd import viz

neopkpd.init_julia()
viz.set_backend("matplotlib")

# Run bootstrap analysis
bootstrap_result = neopkpd.bootstrap_estimation(
    observed_data=data,
    model_spec=model,
    n_bootstrap=1000,
    seed=42
)

# Parameter distributions
fig = viz.plot_bootstrap_distributions(
    bootstrap_result,
    show_original=True,
    ci_level=0.95,
    n_cols=3,
    title="Bootstrap Parameter Distributions"
)
fig.savefig("bootstrap_dist.png", dpi=300)

# CI comparison
fig = viz.plot_bootstrap_ci(
    bootstrap_result,
    methods=["percentile", "bca"],
    title="Confidence Interval Methods"
)
fig.savefig("bootstrap_ci.png", dpi=300)

# Stability
fig = viz.plot_bootstrap_stability(bootstrap_result)
fig.savefig("bootstrap_stability.png", dpi=300)

# Correlations
fig = viz.plot_bootstrap_correlation(bootstrap_result)
fig.savefig("bootstrap_corr.png", dpi=300)
```

---

## See Also

- [Bootstrap Analysis](../estimation/bootstrap.md) - Bootstrap computation
- [Estimation Diagnostics](estimation.md) - More diagnostic plots
- [Population Plots](population.md) - Population visualizations
