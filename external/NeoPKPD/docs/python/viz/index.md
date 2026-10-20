# Visualization

NeoPKPD provides 55+ professional visualization functions with dual matplotlib and plotly backend support.

---

## Overview

The visualization module (`neopkpd.viz`) offers:

- **Dual Backends**: Static (matplotlib) and interactive (plotly)
- **Consistent API**: Same function signatures across backends
- **Publication Quality**: Professional styling and themes
- **Complete Coverage**: All analysis types supported

---

## Visualization Categories

<div class="grid cards" markdown>

-   :material-cog:{ .lg .middle } **Backends & Themes**

    ---

    Configure visualization backend and styling

    [:octicons-arrow-right-24: Backends](backends.md)

-   :material-chart-line:{ .lg .middle } **PK Plots**

    ---

    Concentration-time profiles

    [:octicons-arrow-right-24: PK Plots](pk.md)

-   :material-chart-areaspline:{ .lg .middle } **NCA Plots**

    ---

    Lambda-z fit, AUC visualization

    [:octicons-arrow-right-24: NCA Plots](nca.md)

-   :material-chart-bell-curve:{ .lg .middle } **PKPD Plots**

    ---

    Effect-concentration, hysteresis

    [:octicons-arrow-right-24: PKPD Plots](pkpd.md)

-   :material-chart-scatter-plot:{ .lg .middle } **VPC Plots**

    ---

    Visual predictive checks

    [:octicons-arrow-right-24: VPC Plots](vpc.md)

-   :material-chart-box:{ .lg .middle } **Estimation Diagnostics**

    ---

    Convergence, shrinkage, correlations

    [:octicons-arrow-right-24: Estimation](estimation.md)

-   :material-chart-histogram:{ .lg .middle } **Bootstrap Plots**

    ---

    Distribution and CI visualization

    [:octicons-arrow-right-24: Bootstrap](bootstrap.md)

-   :material-chart-waterfall:{ .lg .middle } **Sensitivity Plots**

    ---

    Tornado, spider, heatmap

    [:octicons-arrow-right-24: Sensitivity](sensitivity.md)

-   :material-account-group:{ .lg .middle } **Population Plots**

    ---

    Forest plots, distributions

    [:octicons-arrow-right-24: Population](population.md)

-   :material-flask:{ .lg .middle } **Trial Plots**

    ---

    Power curves, Kaplan-Meier

    [:octicons-arrow-right-24: Trial](trial.md)

</div>

---

## Quick Start

### Backend Selection

```python
from neopkpd import viz

# Set matplotlib backend (default)
viz.set_backend("matplotlib")

# Or use plotly for interactive plots
viz.set_backend("plotly")

# Check current backend
print(viz.get_backend())

# List available backends
print(viz.available_backends())  # ["matplotlib", "plotly"]
```

### Basic Plotting

```python
import neopkpd
from neopkpd import viz

neopkpd.init_julia()
viz.set_backend("matplotlib")

# Run simulation
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[t * 0.5 for t in range(49)]
)

# Create plot
fig = viz.plot_conc_time(result, title="PK Profile")
fig.savefig("pk_profile.png", dpi=300)
```

---

## Complete Function List

### PK Plots (5 functions)

| Function | Description |
|----------|-------------|
| `plot_conc_time` | Single concentration-time profile |
| `plot_multi_conc_time` | Multiple profiles overlay |
| `plot_spaghetti` | Population spaghetti plot |
| `plot_mean_ribbon` | Mean with confidence ribbon |
| `plot_individual_fits` | Grid of individual fits |

### NCA Plots (3 functions)

| Function | Description |
|----------|-------------|
| `plot_lambda_z_fit` | Terminal phase regression |
| `plot_auc_visualization` | AUC shaded area |
| `plot_dose_proportionality` | Dose vs exposure |

### PKPD Plots (3 functions)

| Function | Description |
|----------|-------------|
| `plot_effect_conc` | Effect vs concentration |
| `plot_hysteresis` | Hysteresis loop |
| `plot_dose_response` | Dose-response curve |

### VPC Plots (5 functions)

| Function | Description |
|----------|-------------|
| `plot_vpc_detailed` | Full VPC with percentiles and CI |
| `plot_pcvpc` | Prediction-corrected VPC |
| `plot_stratified_vpc` | VPC by strata |
| `plot_vpc_with_blq` | VPC with BLQ handling |
| `plot_vpc_ci` | VPC confidence intervals |

### Estimation Diagnostics (10 functions)

| Function | Description |
|----------|-------------|
| `plot_parameter_estimates` | Forest plot of theta |
| `plot_omega_matrix` | Omega covariance heatmap |
| `plot_convergence` | OFV vs iteration |
| `plot_parameter_convergence` | Parameter traces |
| `plot_shrinkage` | Eta shrinkage bars |
| `plot_eta_distributions` | Eta histograms |
| `plot_individual_parameters` | EBE distributions |
| `plot_ofv_comparison` | Model OFV comparison |
| `plot_correlation_matrix` | Parameter correlations |
| `plot_sigma_residuals` | Residual error |

### Bootstrap Plots (4 functions)

| Function | Description |
|----------|-------------|
| `plot_bootstrap_distributions` | Parameter histograms |
| `plot_bootstrap_ci` | CI comparison |
| `plot_bootstrap_stability` | Stability over runs |
| `plot_bootstrap_correlation` | Inter-parameter correlation |

### Sensitivity Plots (5 functions)

| Function | Description |
|----------|-------------|
| `plot_tornado` | Tornado diagram |
| `plot_spider` | Spider/radar plot |
| `plot_sensitivity_heatmap` | Parameter-metric heatmap |
| `plot_waterfall` | Ranked sensitivities |
| `plot_one_at_a_time` | OFAT curves |

### Population Plots (8 functions)

| Function | Description |
|----------|-------------|
| `plot_vpc` | Basic VPC |
| `plot_parameter_distributions` | Parameter histograms |
| `plot_forest` | Forest plot |
| `plot_boxplot` | Box plot comparison |
| `plot_goodness_of_fit` | GOF panel |
| `plot_estimation_summary` | Summary dashboard |
| `plot_sensitivity` | Parameter sensitivity |
| `plot_sensitivity_tornado` | Sensitivity tornado |

### Trial Plots (4 functions)

| Function | Description |
|----------|-------------|
| `plot_power_curve` | Power vs sample size |
| `plot_trial_tornado` | Trial sensitivity |
| `plot_kaplan_meier` | Survival curves |
| `plot_endpoint_distribution` | Endpoint histograms |

---

## Common Parameters

All visualization functions share these parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `backend` | `str` | Current | "matplotlib" or "plotly" |
| `title` | `str` | None | Plot title |
| `xlabel` | `str` | Auto | X-axis label |
| `ylabel` | `str` | Auto | Y-axis label |
| `figsize` | `tuple` | (10, 6) | Figure size (inches) |
| `theme` | `str` | "neopkpd" | Color theme |
| `save_path` | `str` | None | Path to save figure |

---

## Themes

```python
from neopkpd import viz

# Set theme
viz.set_theme("neopkpd")      # Default professional theme
viz.set_theme("publication")    # Minimal for publications
viz.set_theme("presentation")   # Bold for slides

# Available themes
print(viz.available_themes())

# Access color palette
colors = viz.NEOPKPD_COLORS
print(colors)
# {"primary": "#3498DB", "secondary": "#2ECC71", ...}
```

---

## Examples

### Population Spaghetti with Mean

```python
pop_result = neopkpd.simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=100, seed=42,
    omegas={"CL": 0.3, "V": 0.2}
)

# Spaghetti plot
fig = viz.plot_spaghetti(pop_result, alpha=0.2, n_subjects=50)
fig.savefig("spaghetti.png", dpi=300)

# Mean with 90% prediction interval
fig = viz.plot_mean_ribbon(
    pop_result,
    ci_levels=[0.05, 0.95],
    show_median=True,
    show_individual=False
)
fig.savefig("mean_ribbon.png", dpi=300)
```

### VPC Plot

```python
fig = viz.plot_vpc_detailed(
    simulated_data,
    observed_data,
    prediction_intervals=[0.05, 0.50, 0.95],
    show_ci=True,
    show_binning=True
)
fig.savefig("vpc.png", dpi=300)
```

### Tornado Plot

```python
sensitivity_results = [
    {"parameter": "CL", "low": -0.3, "high": 0.25},
    {"parameter": "V", "low": -0.15, "high": 0.18},
    {"parameter": "Ka", "low": -0.4, "high": 0.35},
]

fig = viz.plot_tornado(
    sensitivity_results,
    baseline_value=0.0,
    title="Parameter Sensitivity"
)
fig.savefig("tornado.png", dpi=300)
```

---

## Saving Plots

### Matplotlib

```python
viz.set_backend("matplotlib")
fig = viz.plot_conc_time(result)

# PNG (raster)
fig.savefig("plot.png", dpi=300, bbox_inches="tight")

# PDF (vector)
fig.savefig("plot.pdf", bbox_inches="tight")

# SVG (vector)
fig.savefig("plot.svg", bbox_inches="tight")
```

### Plotly

```python
viz.set_backend("plotly")
fig = viz.plot_conc_time(result)

# Interactive HTML
fig.write_html("plot.html")

# Static image (requires kaleido)
fig.write_image("plot.png", scale=2)
fig.write_image("plot.pdf")
fig.write_image("plot.svg")
```

---

## Next Steps

- [Backends & Themes](backends.md) - Detailed configuration
- [PK Plots](pk.md) - Concentration-time visualization
- [VPC Plots](vpc.md) - Visual predictive checks
- [Estimation Plots](estimation.md) - Diagnostic plots
