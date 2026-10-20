# Estimation Diagnostics

Parameter estimation diagnostic visualization.

---

## Overview

Estimation diagnostics help assess model fit quality and parameter uncertainty.

```python
from neopkpd import viz

fig = viz.plot_convergence(est_result, title="FOCE Convergence")
```

---

## Functions

### plot_convergence

OFV vs iteration trace:

```python
def plot_convergence(
    est_result: EstimationResult | dict,
    *,
    show_final: bool = True,
    title: str | None = None,
    xlabel: str = "Iteration",
    ylabel: str = "Objective Function Value",
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

### plot_parameter_estimates

Forest plot of parameter estimates with CI:

```python
def plot_parameter_estimates(
    est_result: EstimationResult | dict,
    *,
    show_ci: bool = True,
    ci_level: float = 0.95,
    reference_line: float | None = None,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

### plot_omega_matrix

Omega covariance matrix heatmap:

```python
def plot_omega_matrix(
    est_result: EstimationResult | dict,
    *,
    annotate: bool = True,
    title: str | None = None,
    figsize: tuple = (8, 8),
    backend: str | None = None
) -> Figure:
```

### plot_parameter_convergence

Individual parameter traces:

```python
def plot_parameter_convergence(
    est_result: EstimationResult | dict,
    params: list[str] | None = None,
    *,
    title: str | None = None,
    figsize: tuple = (12, 8),
    backend: str | None = None
) -> Figure:
```

### plot_shrinkage

Eta shrinkage bar chart:

```python
def plot_shrinkage(
    est_result: EstimationResult | dict,
    *,
    threshold: float = 0.3,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

### plot_eta_distributions

Histograms of random effects:

```python
def plot_eta_distributions(
    est_result: EstimationResult | dict,
    *,
    show_normal: bool = True,
    n_cols: int = 3,
    title: str | None = None,
    figsize: tuple | None = None,
    backend: str | None = None
) -> Figure:
```

### plot_individual_parameters

EBE distributions vs population:

```python
def plot_individual_parameters(
    est_result: EstimationResult | dict,
    *,
    show_typical: bool = True,
    n_cols: int = 3,
    title: str | None = None,
    figsize: tuple | None = None,
    backend: str | None = None
) -> Figure:
```

### plot_correlation_matrix

Parameter correlation heatmap:

```python
def plot_correlation_matrix(
    est_result: EstimationResult | dict,
    *,
    annotate: bool = True,
    title: str | None = None,
    figsize: tuple = (10, 10),
    backend: str | None = None
) -> Figure:
```

### plot_ofv_comparison

Model comparison by OFV:

```python
def plot_ofv_comparison(
    results: list[EstimationResult | dict],
    labels: list[str],
    *,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

### plot_sigma_residuals

Residual error visualization:

```python
def plot_sigma_residuals(
    est_result: EstimationResult | dict,
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

# Run estimation (example)
est_result = neopkpd.estimate_foce(
    observed_data=data,
    model_spec=model,
    initial_params=init_params
)

# Convergence plot
fig = viz.plot_convergence(est_result, title="FOCE-I Convergence")
fig.savefig("convergence.png", dpi=300)

# Parameter estimates with CI
fig = viz.plot_parameter_estimates(est_result, title="Parameter Estimates")
fig.savefig("parameters.png", dpi=300)

# Omega matrix
fig = viz.plot_omega_matrix(est_result, title="Random Effects Covariance")
fig.savefig("omega.png", dpi=300)

# Shrinkage
fig = viz.plot_shrinkage(est_result, threshold=0.3)
fig.savefig("shrinkage.png", dpi=300)

# Eta distributions
fig = viz.plot_eta_distributions(est_result)
fig.savefig("etas.png", dpi=300)
```

---

## See Also

- [FOCE Method](../estimation/foce.md) - FOCE estimation
- [SAEM Algorithm](../estimation/saem.md) - SAEM estimation
- [Bootstrap Plots](bootstrap.md) - Uncertainty visualization
