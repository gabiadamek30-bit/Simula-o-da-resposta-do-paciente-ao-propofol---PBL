# Terminal Phase Analysis

Comprehensive documentation for terminal elimination phase characterization including lambda_z estimation and half-life calculation.

---

## Overview

Terminal phase analysis estimates the elimination rate constant (lambda_z) from the terminal portion of the concentration-time profile, enabling calculation of half-life and extrapolation of AUC to infinity.

---

## Lambda_z Estimation

### Basic Usage

```julia
using NeoPKPD

times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
conc = [0.0, 1.8, 2.5, 2.0, 1.2, 0.6, 0.3, 0.075]

# From NCA result
result = run_nca(times, conc, 100.0)
println("Lambda_z: $(result.lambda_z_result.lambda_z)")
println("t1/2: $(result.t_half)")

# Direct estimation
lambda_z_result = estimate_lambda_z(times, conc, NCAConfig())
```

### Lambda_z Result Structure

```julia
struct LambdaZResult
    lambda_z::Float64           # Elimination rate constant (1/h)
    intercept::Float64          # Y-intercept of log-linear regression
    r_squared::Float64          # Coefficient of determination
    adjusted_r_squared::Float64 # Adjusted R²
    n_points::Int               # Number of points used
    start_idx::Int              # Index of first point used
    end_idx::Int                # Index of last point used
    times_used::Vector{Float64} # Time points included
    conc_used::Vector{Float64}  # Concentrations included
end
```

---

## Point Selection Methods

### MinPointsFirst (FDA/EMA Default)

Uses minimum required points starting from the latest time points:

```julia
config = NCAConfig(
    lambda_z_selection = MinPointsFirst(),
    lambda_z_min_points = 3,
    lambda_z_r2_threshold = 0.9
)

result = run_nca(times, conc, 100.0; config=config)
```

This method:
1. Starts with the last 3 points (minimum)
2. Calculates R² for log-linear regression
3. Adds earlier points if R² improves
4. Stops when R² drops below threshold

### MaxAdjR2

Selects points that maximize adjusted R²:

```julia
config = NCAConfig(
    lambda_z_selection = MaxAdjR2(),
    lambda_z_min_points = 3,
    lambda_z_max_points = 10
)

result = run_nca(times, conc, 100.0; config=config)
```

This method:
1. Tests all valid combinations of terminal points
2. Selects combination with highest adjusted R²
3. Ensures minimum number of points is met

### Manual Point Selection

Specify exact points to use:

```julia
# Use points from index 5 to 8
config = NCAConfig(
    lambda_z_start_idx = 5,
    lambda_z_end_idx = 8
)

result = run_nca(times, conc, 100.0; config=config)
```

---

## Quality Criteria

### R² Threshold

```julia
# Require R² ≥ 0.95 for lambda_z estimation
config = NCAConfig(
    lambda_z_r2_threshold = 0.95
)

result = run_nca(times, conc, 100.0; config=config)

# Check if threshold was met
if result.lambda_z_result.r_squared >= 0.95
    println("Lambda_z estimation meets quality criteria")
else
    println("WARNING: R² below threshold")
end
```

### Minimum Points

```julia
# Require at least 4 points for lambda_z
config = NCAConfig(
    lambda_z_min_points = 4
)
```

### Extrapolation Warning

```julia
# Warn if AUC extrapolation exceeds 20%
config = NCAConfig(
    extrapolation_max_pct = 20.0
)

result = run_nca(times, conc, 100.0; config=config)

if result.auc_extra_pct > 20.0
    println("WARNING: $(result.auc_extra_pct)% extrapolation exceeds threshold")
end
```

---

## Half-Life Calculation

### Terminal Half-Life

$$t_{1/2} = \frac{\ln(2)}{\lambda_z} = \frac{0.693}{\lambda_z}$$

```julia
result = run_nca(times, conc, 100.0)
println("t1/2: $(result.t_half) hours")

# Direct calculation
t_half = log(2) / result.lambda_z_result.lambda_z
```

### Effective Half-Life (Multiple Dose)

For multiple dose analysis:

```julia
result = run_nca(times, conc, 100.0; dosing_type=:multiple, tau=24.0)

# Effective half-life from accumulation
t_half_eff = nca_effective_half_life(result.accumulation_index, result.tau)
```

---

## Mean Residence Time (MRT)

### MRT for Extravascular Administration

$$MRT = \frac{AUMC_{0-\infty}}{AUC_{0-\infty}}$$

```julia
result = run_nca(times, conc, 100.0; route=:extravascular)
println("MRT: $(result.mrt) hours")
```

### MRT for IV Bolus

$$MRT = \frac{AUMC_{0-\infty}}{AUC_{0-\infty}} - \frac{T_{inf}}{2}$$

For IV bolus (instantaneous), no infusion time correction:

```julia
result = run_nca(times, conc, 100.0; route=:iv_bolus)
println("MRT (IV): $(result.mrt) hours")
```

### MRT for IV Infusion

```julia
# 1-hour infusion
result = run_nca(times, conc, 100.0; route=:iv_infusion, infusion_time=1.0)
println("MRT (corrected): $(result.mrt) hours")
```

---

## PK Parameters from Terminal Phase

### Clearance (CL/F)

$$CL/F = \frac{Dose}{AUC_{0-\infty}}$$

```julia
result = run_nca(times, conc, 100.0)
println("CL/F: $(result.cl_f) L/h")
```

### Volume of Distribution at Terminal Phase (Vz/F)

$$V_z/F = \frac{CL/F}{\lambda_z} = \frac{Dose}{AUC_{0-\infty} \cdot \lambda_z}$$

```julia
println("Vz/F: $(result.vz_f) L")
```

### Volume of Distribution at Steady State (Vss/F)

$$V_{ss}/F = MRT \cdot CL/F$$

```julia
println("Vss/F: $(result.vss_f) L")
```

---

## Terminal Phase Regression Diagnostics

### Visualizing the Fit

```julia
# Get regression data
lz = result.lambda_z_result
times_used = lz.times_used
conc_used = lz.conc_used

# Log-linear regression
log_conc = log.(conc_used)
predicted = lz.intercept .- lz.lambda_z .* times_used

println("Points used: $(lz.n_points)")
println("R²: $(round(lz.r_squared, digits=4))")
println("Adjusted R²: $(round(lz.adjusted_r_squared, digits=4))")
```

### Residual Analysis

```julia
# Calculate residuals
log_conc_obs = log.(lz.conc_used)
log_conc_pred = lz.intercept .- lz.lambda_z .* lz.times_used
residuals = log_conc_obs .- log_conc_pred

# Summary statistics
println("Mean residual: $(mean(residuals))")
println("SD residual: $(std(residuals))")
```

---

## Handling Special Cases

### Insufficient Terminal Points

```julia
# Short profile with few points
short_times = [0.0, 1.0, 2.0]
short_conc = [0.0, 2.5, 1.5]

config = NCAConfig(lambda_z_min_points = 3)
result = run_nca(short_times, short_conc, 100.0; config=config)

if isnothing(result.lambda_z_result) || isnan(result.lambda_z_result.lambda_z)
    println("WARNING: Insufficient data for lambda_z estimation")
    println("AUC0-inf cannot be calculated")
end
```

### Multi-Phasic Elimination

For drugs with multiple elimination phases:

```julia
# Focus on terminal phase only
# Exclude distribution phase points
config = NCAConfig(
    lambda_z_start_time = 4.0,  # Start after distribution
    lambda_z_min_points = 3
)

result = run_nca(times, conc, 100.0; config=config)
```

### Concentration Below LLOQ

```julia
# Handle BLQ in terminal phase
config = NCAConfig(
    lloq = 0.05,
    blq_handling = BLQMissing()  # Exclude BLQ from lambda_z
)

result = run_nca(times, conc, 100.0; config=config)
```

---

## Example: Complete Terminal Phase Analysis

```julia
using NeoPKPD

# PK data (oral administration)
times = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0]
conc = [0.0, 2.1, 4.5, 5.8, 5.2, 4.5, 3.2, 2.4, 1.3, 0.72, 0.24, 0.02]
dose = 500.0  # mg

# Configure terminal phase analysis
config = NCAConfig(
    method = LinLogMixedMethod(),
    lambda_z_selection = MinPointsFirst(),
    lambda_z_min_points = 3,
    lambda_z_r2_threshold = 0.9,
    extrapolation_max_pct = 20.0,
    lloq = 0.01
)

# Run NCA
result = run_nca(times, conc, dose; config=config, route=:extravascular)

# Terminal phase results
lz = result.lambda_z_result

println("=== Terminal Phase Analysis ===")
println("Lambda_z: $(round(lz.lambda_z, digits=4)) 1/h")
println("t1/2: $(round(result.t_half, digits=2)) h")
println("MRT: $(round(result.mrt, digits=2)) h")
println("")
println("=== Regression Quality ===")
println("R²: $(round(lz.r_squared, digits=4))")
println("Adjusted R²: $(round(lz.adjusted_r_squared, digits=4))")
println("Points used: $(lz.n_points) (indices $(lz.start_idx) to $(lz.end_idx))")
println("Times used: $(lz.times_used)")
println("")
println("=== PK Parameters ===")
println("CL/F: $(round(result.cl_f, digits=2)) L/h")
println("Vz/F: $(round(result.vz_f, digits=1)) L")
println("Vss/F: $(round(result.vss_f, digits=1)) L")
println("")
println("=== AUC Extrapolation ===")
println("AUC0-t: $(round(result.auc_0_t, digits=2)) mg·h/L")
println("AUC0-inf: $(round(result.auc_0_inf, digits=2)) mg·h/L")
println("Extrapolation: $(round(result.auc_extra_pct, digits=1))%")

# Quality assessment
println("\n=== Quality Assessment ===")
if lz.r_squared >= 0.9
    println("✓ R² meets threshold (≥0.9)")
else
    println("⚠ R² below threshold")
end

if result.auc_extra_pct <= 20.0
    println("✓ AUC extrapolation acceptable (≤20%)")
else
    println("⚠ High AUC extrapolation")
end

if lz.n_points >= 3
    println("✓ Sufficient terminal points")
else
    println("⚠ Insufficient terminal points")
end
```

---

## Formulas Summary

| Parameter | Formula |
|-----------|---------|
| Lambda_z | Slope of ln(C) vs time regression |
| t1/2 | $\ln(2) / \lambda_z$ |
| MRT (extravascular) | $AUMC_{0-\infty} / AUC_{0-\infty}$ |
| MRT (IV infusion) | $AUMC_{0-\infty} / AUC_{0-\infty} - T_{inf}/2$ |
| CL/F | $Dose / AUC_{0-\infty}$ |
| Vz/F | $Dose / (AUC_{0-\infty} \cdot \lambda_z)$ |
| Vss/F | $MRT \cdot CL/F$ |

---

## See Also

- [Exposure Metrics](exposure-metrics.md) - AUC calculations
- [Multiple Dose](multiple-dose.md) - Steady-state analysis
- [Population NCA](population-nca.md) - Multi-subject analysis
- [Bioequivalence](bioequivalence.md) - BE studies

