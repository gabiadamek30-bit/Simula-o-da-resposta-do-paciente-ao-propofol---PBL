# NCA Plots

Non-compartmental analysis visualization functions.

---

## Overview

NCA plots visualize pharmacokinetic parameters derived from concentration-time data without assuming a specific model.

```python
from neopkpd import viz

fig = viz.plot_lambda_z_fit(nca_result, title="Terminal Phase")
```

---

## Functions

### plot_lambda_z_fit

Terminal phase regression visualization:

```python
def plot_lambda_z_fit(
    nca_result: NCAResult | dict,
    *,
    show_regression: bool = True,
    show_r2: bool = True,
    log_scale: bool = True,
    title: str | None = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: tuple = (10, 6),
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

**Usage:**

```python
nca_result = neopkpd.run_nca(
    times=[0, 0.5, 1, 2, 4, 8, 12, 24],
    concentrations=[0, 5.2, 8.1, 6.3, 3.8, 1.9, 0.9, 0.2],
    dose=100.0,
    route="oral"
)

# Lambda-z regression plot
fig = viz.plot_lambda_z_fit(nca_result, title="Terminal Phase Analysis")
```

### plot_auc_visualization

AUC with shaded area:

```python
def plot_auc_visualization(
    nca_result: NCAResult | dict,
    *,
    show_auc_last: bool = True,
    show_auc_extrap: bool = True,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Show AUC components
fig = viz.plot_auc_visualization(
    nca_result,
    show_auc_last=True,
    show_auc_extrap=True,
    title="AUC Breakdown"
)
```

### plot_dose_proportionality

Dose vs exposure relationship:

```python
def plot_dose_proportionality(
    nca_results: list[NCAResult | dict],
    doses: list[float],
    metric: str = "auc",
    *,
    show_regression: bool = True,
    log_scale: bool = False,
    title: str | None = None,
    figsize: tuple = (10, 6),
    backend: str | None = None
) -> Figure:
```

**Usage:**

```python
# Multiple dose NCA results
doses = [25, 50, 100, 200, 400]
nca_results = [run_nca_for_dose(d) for d in doses]

# Dose proportionality plot
fig = viz.plot_dose_proportionality(
    nca_results,
    doses=doses,
    metric="auc",
    title="Dose Proportionality"
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

# Sample PK data
times = [0, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24]
conc = [0, 3.2, 6.8, 9.1, 7.5, 4.2, 2.8, 1.9, 0.8, 0.15]

# Run NCA
nca_result = neopkpd.run_nca(
    times=times,
    concentrations=conc,
    dose=100.0,
    route="oral"
)

# Lambda-z fit
fig = viz.plot_lambda_z_fit(
    nca_result,
    title=f"Terminal Phase (t1/2 = {nca_result['half_life']:.1f} hr)"
)
fig.savefig("lambda_z.png", dpi=300)

# AUC visualization
fig = viz.plot_auc_visualization(nca_result, title="AUC Components")
fig.savefig("auc_viz.png", dpi=300)

print(f"AUC0-inf: {nca_result['auc_inf']:.1f} mg*hr/L")
print(f"Cmax: {nca_result['cmax']:.2f} mg/L")
print(f"t1/2: {nca_result['half_life']:.1f} hr")
```

---

## See Also

- [NCA Module](../nca/index.md) - NCA computation
- [Bioequivalence](../nca/bioequivalence.md) - BE analysis
- [PK Plots](pk.md) - Concentration-time plots
