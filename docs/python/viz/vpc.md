# VPC Visualization

Comprehensive guide to Visual Predictive Check visualization in Python with dual matplotlib and plotly backend support.

---

## Overview

NeoPKPD provides 5 specialized VPC visualization functions:

| Function | Description |
|----------|-------------|
| `plot_vpc` | Standard VPC with percentile ribbons |
| `plot_pcvpc` | Prediction-corrected VPC |
| `plot_stratified_vpc` | Multi-panel faceted VPC |
| `plot_vpc_with_blq` | VPC with BLQ statistics panel |
| `plot_vpc_ci` | Detailed confidence interval visualization |

```python
from neopkpd import viz

# Set backend
viz.set_backend("matplotlib")  # or "plotly"

# All VPC functions follow same pattern
fig = viz.plot_vpc(vpc_result, title="Model Validation")
```

---

## VPC Data Structures

### VPCResult

```python
from dataclasses import dataclass
from typing import List, Dict, Optional

@dataclass
class VPCResult:
    """Result of VPC computation."""
    config: VPCConfig                    # Configuration used
    bins: List[VPCBin]                   # Computed bins
    n_subjects_observed: int             # Number of subjects
    n_observations_observed: int         # Total observations
    n_simulations: int                   # Simulations performed
    strata: str                          # Strata label
    simulation_seed: int                 # Seed used

@dataclass
class VPCBin:
    """Data for a single time bin."""
    bin_id: int                          # Bin identifier
    time_min: float                      # Lower bound
    time_max: float                      # Upper bound
    time_midpoint: float                 # Midpoint for plotting
    n_observed: int                      # Observations in bin
    n_simulated: int                     # Simulated per replicate
    percentiles: List[VPCPercentileData]

@dataclass
class VPCPercentileData:
    """Percentile data for a bin."""
    percentile: float                    # Level (e.g., 0.50)
    observed: float                      # Observed percentile
    simulated_median: float              # Median of simulations
    simulated_lower: float               # Lower CI
    simulated_upper: float               # Upper CI
```

### Accessor Functions

```python
from neopkpd.vpc import (
    get_bin_midpoints,
    get_observed_percentile,
    get_simulated_median,
    get_simulated_ci,
    get_blq_observed,
    get_blq_simulated
)

# Extract data for custom plotting
times = get_bin_midpoints(vpc_result)                    # np.ndarray
obs_median = get_observed_percentile(vpc_result, 0.50)   # np.ndarray
sim_median = get_simulated_median(vpc_result, 0.50)      # np.ndarray
sim_lower, sim_upper = get_simulated_ci(vpc_result, 0.50)  # tuple

# BLQ data (if available)
blq_obs = get_blq_observed(vpc_result)
blq_sim_med, blq_sim_lo, blq_sim_hi = get_blq_simulated(vpc_result)
```

---

## plot_vpc

Standard VPC visualization with observed percentiles and simulated confidence ribbons.

### Function Signature

```python
def plot_vpc(
    vpc_result: VPCResult,
    *,
    log_scale: bool = False,
    show_observed: bool = True,
    show_simulated_median: bool = True,
    show_ci: bool = True,
    pi_levels: list[float] | None = None,
    title: str | None = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: tuple[float, float] = (10, 6),
    obs_color: str | None = None,
    sim_color: str | None = None,
    ci_alpha: float = 0.2,
    median_alpha: float = 0.5,
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

### Basic Usage

```python
from neopkpd import viz
import neopkpd

# Initialize Julia connection
neopkpd.init_julia()

# Compute VPC (example)
vpc_result = neopkpd.compute_vpc(
    observed_data=observed,
    population_spec=pop_spec,
    grid={"t0": 0.0, "t1": 24.0, "saveat": 0.5},
    config={"n_simulations": 500, "pi_levels": [0.05, 0.50, 0.95]}
)

# Create VPC plot
fig = viz.plot_vpc(
    vpc_result,
    title="Visual Predictive Check",
    xlabel="Time (hours)",
    ylabel="Concentration (mg/L)"
)

fig.savefig("vpc.png", dpi=300)
```

### Customization Options

```python
# Linear scale (default)
fig = viz.plot_vpc(vpc_result, log_scale=False)

# Semi-log scale
fig = viz.plot_vpc(vpc_result, log_scale=True)

# Show only specific elements
fig = viz.plot_vpc(
    vpc_result,
    show_observed=True,
    show_simulated_median=False,  # Hide simulated median lines
    show_ci=True
)

# Custom percentile levels display
fig = viz.plot_vpc(
    vpc_result,
    pi_levels=[0.10, 0.50, 0.90]  # Override displayed levels
)

# Custom colors
fig = viz.plot_vpc(
    vpc_result,
    obs_color="#E74C3C",    # Red for observed
    sim_color="#3498DB",    # Blue for simulated
    ci_alpha=0.3            # More opaque ribbons
)
```

### Plot Elements

The standard VPC plot shows:

1. **Observed percentiles** - Solid lines at each PI level
2. **Simulated CI ribbons** - Shaded areas showing simulation uncertainty
3. **Simulated median** - Dashed lines (optional)

```
Y-axis: Concentration
X-axis: Time
Legend:
  - Observed P5/P50/P95 (solid lines)
  - Simulated 95% CI (shaded ribbons)
  - Simulated median (dashed lines)
```

---

## plot_pcvpc

Prediction-corrected VPC for variable dosing or covariate-adjusted data.

### Function Signature

```python
def plot_pcvpc(
    vpc_result: VPCResult,
    *,
    log_scale: bool = False,
    show_observed: bool = True,
    show_ci: bool = True,
    title: str | None = None,
    xlabel: str = "Time",
    ylabel: str = "Prediction-Corrected Concentration",
    figsize: tuple[float, float] = (10, 6),
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

### Usage

```python
# Compute pcVPC
pcvpc_result = neopkpd.compute_pcvpc(
    observed_data=observed,
    population_spec=pop_spec,
    grid=grid,
    config={"n_simulations": 500}
)

# Plot pcVPC
fig = viz.plot_pcvpc(
    pcvpc_result,
    title="Prediction-Corrected VPC",
    ylabel="Prediction-Corrected Concentration"
)

fig.savefig("pcvpc.png", dpi=300)
```

### When to Use pcVPC

- Variable dosing across subjects
- Significant covariate effects on PK
- Dose escalation studies
- Different formulations/routes

```python
# Example: Dose escalation study
# Doses: 50, 100, 200, 400 mg

# Standard VPC would show dose-related variability
fig_standard = viz.plot_vpc(standard_vpc, title="Standard VPC")

# pcVPC normalizes by dose, showing only IIV
fig_pc = viz.plot_pcvpc(pcvpc_result, title="pcVPC (dose-normalized)")
```

---

## plot_stratified_vpc

Multi-panel VPC faceted by covariate strata.

### Function Signature

```python
def plot_stratified_vpc(
    stratified_result: StratifiedVPCResult,
    *,
    log_scale: bool = False,
    n_cols: int = 2,
    subplot_size: tuple[float, float] = (5, 4),
    share_y: bool = True,
    show_observed: bool = True,
    show_ci: bool = True,
    title: str | None = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

### Usage

```python
# Compute stratified VPC
stratified_result = neopkpd.compute_stratified_vpc(
    observed_data=observed,
    population_spec=pop_spec,
    grid=grid,
    strata_data=strata_data,
    config={"stratify_by": ["DOSE"], "n_simulations": 500}
)

# Plot multi-panel VPC
fig = viz.plot_stratified_vpc(
    stratified_result,
    n_cols=3,                    # 3 columns layout
    subplot_size=(4, 3),         # Size per subplot
    title="VPC by Dose Group"
)

fig.savefig("stratified_vpc.png", dpi=300, bbox_inches="tight")
```

### Layout Options

```python
# Vertical layout (1 column)
fig = viz.plot_stratified_vpc(
    stratified_result,
    n_cols=1,
    subplot_size=(8, 3)
)

# Grid layout
fig = viz.plot_stratified_vpc(
    stratified_result,
    n_cols=2,
    share_y=True  # Same Y-axis scale across panels
)

# Independent Y-axes (for very different concentrations)
fig = viz.plot_stratified_vpc(
    stratified_result,
    share_y=False
)
```

### Stratum Labels

Stratum names from the result are used as subplot titles:

```python
# Example output layout:
# ┌─────────────────┬─────────────────┐
# │   50 mg         │   150 mg        │
# │   [VPC plot]    │   [VPC plot]    │
# ├─────────────────┼─────────────────┤
# │   300 mg        │   600 mg        │
# │   [VPC plot]    │   [VPC plot]    │
# └─────────────────┴─────────────────┘
```

---

## plot_vpc_with_blq

VPC with Below Limit of Quantification statistics panel.

### Function Signature

```python
def plot_vpc_with_blq(
    vpc_result: VPCResult,
    blq_stats: list[BLQBinStats] | None = None,
    *,
    log_scale: bool = False,
    blq_panel_height: float = 0.25,
    title: str | None = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    blq_ylabel: str = "% BLQ",
    figsize: tuple[float, float] = (10, 8),
    obs_color: str | None = None,
    blq_color: str | None = None,
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

### Usage

```python
from neopkpd.vpc import BLQMethod

# Compute VPC with BLQ handling
vpc_result, blq_stats = neopkpd.compute_vpc_with_blq(
    observed_data=observed,
    population_spec=pop_spec,
    grid=grid,
    config={
        "lloq": 0.1,  # Lower limit of quantification
        "n_simulations": 500
    },
    blq_method=BLQMethod.M4  # Replace with LLOQ/2
)

# Two-panel plot: VPC + BLQ fraction
fig = viz.plot_vpc_with_blq(
    vpc_result,
    blq_stats=blq_stats,
    title="VPC with BLQ Analysis",
    blq_ylabel="% Below LOQ"
)

fig.savefig("vpc_blq.png", dpi=300)
```

### BLQ Panel

The lower panel shows:

- **Observed %BLQ** - Solid line/markers
- **Simulated %BLQ CI** - Shaded ribbon
- **Simulated %BLQ median** - Dashed line

```
Upper panel: Standard VPC (75% height)
Lower panel: %BLQ comparison (25% height)
```

### BLQ Methods

```python
from neopkpd.vpc import BLQMethod

# Available methods (Beal 2001)
BLQMethod.M1  # Discard all BLQ observations
BLQMethod.M3  # Treat as censored (keep original)
BLQMethod.M4  # Replace with LLOQ/2 (default)
BLQMethod.M5  # 0 before Tmax, LLOQ/2 after
BLQMethod.M6  # LLOQ/2 before Tmax, discard after
BLQMethod.M7  # 0 before Tmax, discard after
```

---

## plot_vpc_ci

Detailed confidence interval visualization with nested ribbons.

### Function Signature

```python
def plot_vpc_ci(
    vpc_result: VPCResult,
    *,
    log_scale: bool = False,
    ci_levels: list[float] = [0.05, 0.95],
    show_bin_boundaries: bool = False,
    title: str | None = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: tuple[float, float] = (10, 6),
    backend: str | None = None,
    save_path: str | None = None
) -> Figure:
```

### Usage

```python
# Detailed CI visualization
fig = viz.plot_vpc_ci(
    vpc_result,
    ci_levels=[0.05, 0.95],      # 5th and 95th percentile CIs
    show_bin_boundaries=True,    # Show vertical bin lines
    title="VPC Confidence Intervals"
)

fig.savefig("vpc_ci.png", dpi=300)
```

### Features

- **Nested ribbons** - Inner (median) and outer (5th/95th) CIs
- **Thicker observed lines** - Stand out from simulated
- **Bin boundaries** - Optional vertical lines at bin edges

```python
# Show bin structure
fig = viz.plot_vpc_ci(
    vpc_result,
    show_bin_boundaries=True
)
# Helps assess if binning is appropriate
```

---

## Complete Example

```python
import neopkpd
from neopkpd import viz
from neopkpd.vpc import VPCConfig, QuantileBinning
import numpy as np
import matplotlib.pyplot as plt

# ================================================
# Complete VPC Visualization Example
# ================================================

# 1. Initialize
neopkpd.init_julia()
viz.set_backend("matplotlib")

print("=== VPC Visualization Example ===\n")

# 2. Simulate population data (for demonstration)
pop_result = neopkpd.simulate_population_oral(
    ka=1.5, cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0, saveat=0.5,
    n=100,
    omegas={"Ka": 0.16, "CL": 0.09, "V": 0.04},
    seed=42
)

# 3. Create observed data from population (subset)
observed_data = {
    "subject_ids": [f"S{i}" for i in range(100) for _ in range(7)],
    "times": [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0] * 100,
    "dv": [],  # Extract from pop_result
    "dvid": ["conc"] * 700
}

# Extract concentrations at specific times
sampling_idx = [1, 2, 4, 8, 16, 24, 48]  # indices for 0.5, 1, 2, 4, 8, 12, 24 hr
for ind in pop_result["individuals"]:
    for idx in sampling_idx:
        observed_data["dv"].append(ind["concentrations"][idx] * (1 + 0.1 * np.random.randn()))

# 4. Define population model specification
pop_spec = {
    "model": "OneCompOral",
    "params": {"Ka": 1.5, "CL": 5.0, "V": 50.0},
    "omega": {"Ka": 0.16, "CL": 0.09, "V": 0.04},
    "doses": [{"time": 0.0, "amount": 100.0}],
    "n": 100
}

grid = {"t0": 0.0, "t1": 24.0, "saveat": 0.5}

# 5. Compute VPC variants
print("Computing standard VPC...")
vpc_config = {
    "pi_levels": [0.05, 0.50, 0.95],
    "n_simulations": 500,
    "binning": "quantile",
    "n_bins": 7
}

vpc_result = neopkpd.compute_vpc(
    observed_data=observed_data,
    population_spec=pop_spec,
    grid=grid,
    config=vpc_config
)

# 6. Create visualization panels
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 6a. Standard VPC
plt.sca(axes[0, 0])
viz.plot_vpc(
    vpc_result,
    title="Standard VPC",
    xlabel="Time (hr)",
    ylabel="Concentration (mg/L)"
)

# 6b. Semi-log VPC
plt.sca(axes[0, 1])
viz.plot_vpc(
    vpc_result,
    log_scale=True,
    title="VPC (Semi-log)",
    xlabel="Time (hr)",
    ylabel="Concentration (mg/L)"
)

# 6c. VPC with CI details
plt.sca(axes[1, 0])
viz.plot_vpc_ci(
    vpc_result,
    show_bin_boundaries=True,
    title="VPC with Bin Boundaries",
    xlabel="Time (hr)",
    ylabel="Concentration (mg/L)"
)

# 6d. Custom styled VPC
plt.sca(axes[1, 1])
viz.plot_vpc(
    vpc_result,
    obs_color="#E74C3C",
    sim_color="#2980B9",
    ci_alpha=0.25,
    title="Custom Styled VPC",
    xlabel="Time (hr)",
    ylabel="Concentration (mg/L)"
)

plt.tight_layout()
plt.savefig("vpc_panel.png", dpi=300, bbox_inches="tight")
print("Saved: vpc_panel.png")

# 7. Individual VPC plots
print("\nCreating individual plots...")

# Standard VPC
fig = viz.plot_vpc(
    vpc_result,
    title="Visual Predictive Check - Final Model",
    xlabel="Time (hours)",
    ylabel="Concentration (mg/L)",
    figsize=(10, 6)
)
fig.savefig("vpc_standard.png", dpi=300, bbox_inches="tight")
print("Saved: vpc_standard.png")

# Interactive plotly version
viz.set_backend("plotly")
fig = viz.plot_vpc(
    vpc_result,
    title="Interactive VPC",
    xlabel="Time (hours)",
    ylabel="Concentration (mg/L)"
)
fig.write_html("vpc_interactive.html")
print("Saved: vpc_interactive.html")

print("\n✓ VPC visualization complete")
```

---

## Plotting Parameters Reference

### Common Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `log_scale` | `bool` | `False` | Semi-log Y-axis |
| `title` | `str` | `None` | Plot title |
| `xlabel` | `str` | `"Time"` | X-axis label |
| `ylabel` | `str` | `"Concentration"` | Y-axis label |
| `figsize` | `tuple` | `(10, 6)` | Figure size (inches) |
| `backend` | `str` | Current | `"matplotlib"` or `"plotly"` |
| `save_path` | `str` | `None` | Auto-save path |

### Style Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obs_color` | `str` | Theme default | Observed line color |
| `sim_color` | `str` | Theme default | Simulated ribbon color |
| `ci_alpha` | `float` | `0.2` | Ribbon transparency |
| `median_alpha` | `float` | `0.5` | Median line transparency |

### Behavior Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `show_observed` | `bool` | `True` | Show observed percentiles |
| `show_simulated_median` | `bool` | `True` | Show simulated medians |
| `show_ci` | `bool` | `True` | Show CI ribbons |
| `pi_levels` | `list` | From result | Percentile levels to display |

---

## Backend Comparison

### Matplotlib

```python
viz.set_backend("matplotlib")

# Advantages:
# - Publication-quality static figures
# - PDF/SVG vector output
# - Fine-grained customization
# - Familiar matplotlib API

fig = viz.plot_vpc(vpc_result)
fig.savefig("vpc.pdf", bbox_inches="tight")
```

### Plotly

```python
viz.set_backend("plotly")

# Advantages:
# - Interactive zoom/pan
# - Hover information
# - Web embedding
# - Animation support

fig = viz.plot_vpc(vpc_result)
fig.write_html("vpc.html")
```

---

## Interpretation Guide

### Good Model Fit

```
✓ Observed percentiles within shaded CI regions
✓ No systematic deviation over time
✓ Similar fit at early and late times
✓ Variability captured (5th/95th within bounds)
```

### Potential Issues

| Visual Pattern | Interpretation | Action |
|----------------|----------------|--------|
| Median above upper CI | Under-prediction | Check structural model |
| 95th above CI | Under-estimated variability | Increase omega |
| Early times offset | Absorption issue | Consider transit model |
| Late times offset | Elimination issue | Check clearance |
| Widening deviation | Time-varying misfit | Add time-varying effect |

---

## See Also

- [Julia Standard VPC](../../julia/vpc/standard.md) - VPC methodology
- [Julia pcVPC](../../julia/vpc/pcvpc.md) - Prediction correction
- [Julia Stratified VPC](../../julia/vpc/stratified.md) - Stratification
- [Backends & Themes](backends.md) - Backend configuration
- [Population Plots](population.md) - Related visualizations
