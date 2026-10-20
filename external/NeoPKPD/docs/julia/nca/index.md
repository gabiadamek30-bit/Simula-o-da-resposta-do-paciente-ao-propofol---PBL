# Non-Compartmental Analysis (NCA)

Non-compartmental analysis provides model-independent exposure metrics following FDA and EMA guidelines.

---

## Overview

NCA calculates pharmacokinetic parameters directly from concentration-time data without assuming a specific compartmental model.

### Key Metrics

| Metric | Symbol | Description |
|--------|--------|-------------|
| Maximum Concentration | Cmax | Peak observed concentration |
| Time to Maximum | Tmax | Time of Cmax |
| AUC to Last | AUC0-t | Area under curve to last observation |
| AUC to Infinity | AUC0-∞ | Extrapolated total exposure |
| Terminal Half-life | t½ | ln(2)/λz |
| Terminal Rate | λz | Slope of terminal phase |
| Clearance | CL/F | Dose/AUC (apparent) |
| Volume | Vz/F | CL/(F·λz) (apparent) |
| Mean Residence Time | MRT | AUMC/AUC |

---

## Documentation

<div class="grid cards" markdown>

-   :material-chart-areaspline:{ .lg .middle } **Exposure Metrics**

    ---

    Cmax, Tmax, AUC calculations

    [:octicons-arrow-right-24: Exposure Metrics](exposure-metrics.md)

-   :material-chart-timeline-variant:{ .lg .middle } **Terminal Phase**

    ---

    λz estimation and half-life

    [:octicons-arrow-right-24: Terminal Phase](terminal-phase.md)

-   :material-scale-balance:{ .lg .middle } **Bioequivalence**

    ---

    90% CI and TOST analysis

    [:octicons-arrow-right-24: Bioequivalence](bioequivalence.md)

-   :material-account-group:{ .lg .middle } **Population NCA**

    ---

    NCA for multiple subjects

    [:octicons-arrow-right-24: Population NCA](population-nca.md)

-   :material-pill-multiple:{ .lg .middle } **Multiple Dose**

    ---

    Steady-state metrics and accumulation

    [:octicons-arrow-right-24: Multiple Dose](multiple-dose.md)

</div>

---

## Quick Start

### Basic NCA

```julia
using NeoPKPD

# Concentration-time data
times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
conc = [0.0, 1.8, 2.0, 1.5, 1.0, 0.5, 0.25, 0.06]
dose = 100.0

# Run NCA
result = run_nca(times, conc, dose)

# Access results
println("Cmax: ", result.cmax, " mg/L")
println("Tmax: ", result.tmax, " h")
println("AUC0-t: ", result.auc_0_t, " mg·h/L")
println("AUC0-∞: ", result.auc_0_inf, " mg·h/L")
println("t½: ", result.t_half, " h")
println("CL/F: ", result.cl_f, " L/h")
```

### With Configuration

```julia
config = NCAConfig(
    method = :log_linear,           # AUC calculation method
    lambda_z_min_points = 3,        # Min points for λz
    lambda_z_r2_threshold = 0.9,    # Quality threshold
    extrapolation_max_pct = 20.0,   # Warning if >20% extrapolated
    blq_handling = :zero            # Handle BLQ as zero
)

result = run_nca(times, conc, dose; config=config)

# Check quality metrics
println("λz R²: ", result.lambda_z_r_squared)
println("AUC extrapolated: ", result.auc_extrapolated_pct, "%")
```

---

## AUC Calculation Methods

### Linear Trapezoidal

$$AUC_{t_1 \to t_2} = \frac{(C_1 + C_2)}{2} \cdot (t_2 - t_1)$$

Best for: Ascending portions of the curve

### Log-Linear Trapezoidal

$$AUC_{t_1 \to t_2} = \frac{(C_1 - C_2)}{\ln(C_1/C_2)} \cdot (t_2 - t_1)$$

Best for: Descending (elimination) portions

### Linear-Log Mixed (Recommended)

- Linear trapezoidal for ascending
- Log-linear for descending

```julia
config = NCAConfig(method = :lin_log_mixed)
```

---

## Terminal Phase Analysis

### λz Estimation

The terminal elimination rate constant is estimated by log-linear regression:

```julia
# Manual lambda_z estimation
lambda_z, r_squared, n_points = estimate_lambda_z(times, conc)

# With custom parameters
lambda_z, r_squared, n_points = estimate_lambda_z(
    times, conc;
    min_points = 3,
    r2_threshold = 0.9
)
```

### Half-life

$$t_{1/2} = \frac{\ln(2)}{\lambda_z}$$

---

## NCA Result Structure

```julia
struct NCAResult
    # Primary metrics
    cmax::Float64
    tmax::Float64
    auc_0_t::Float64
    auc_0_inf::Float64
    t_half::Float64

    # Terminal phase
    lambda_z::Float64
    lambda_z_r_squared::Float64
    lambda_z_n_points::Int

    # Clearance and volume
    cl_f::Float64
    vz_f::Float64

    # Additional metrics
    mrt::Float64                    # Mean residence time
    auc_extrapolated_pct::Float64   # % AUC extrapolated

    # Metadata
    method::Symbol
    dose::Float64
end
```

---

## FDA/EMA Compliance

NeoPKPD NCA calculations follow regulatory guidance:

- **FDA Guidance for Industry: Bioavailability and Bioequivalence Studies**
- **EMA Guideline on the Investigation of Bioequivalence**

### Key Requirements

| Requirement | NeoPKPD Implementation |
|-------------|------------------------|
| λz from ≥3 points | Configurable minimum |
| R² > 0.9 for λz | Configurable threshold |
| AUC extrapolation <20% | Warning flag |
| Log-linear interpolation | Supported method |
| BLQ handling | Multiple options |

---

## Next Steps

- [Exposure Metrics](exposure-metrics.md) - Detailed metric calculations
- [Terminal Phase](terminal-phase.md) - Lambda_z and half-life calculation
- [Bioequivalence](bioequivalence.md) - BE analysis methods
- [Population NCA](population-nca.md) - Multi-subject analysis
- [Multiple Dose](multiple-dose.md) - Steady-state analysis
