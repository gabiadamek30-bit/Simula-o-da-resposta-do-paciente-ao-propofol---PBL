# Exposure Metrics

Complete documentation for primary and secondary exposure metrics in NCA following FDA/EMA guidance.

---

## Overview

Exposure metrics characterize drug exposure from concentration-time profiles without assuming a specific model structure.

---

## Primary Exposure Metrics

### Cmax - Maximum Concentration

The highest observed concentration:

```julia
using NeoPKPD

times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
conc = [0.0, 1.8, 2.5, 2.0, 1.2, 0.6, 0.3, 0.075]

# From NCA result
result = run_nca(times, conc, 100.0)
println("Cmax: $(result.cmax)")

# Direct calculation
cmax = nca_cmax(conc)
```

### Tmax - Time of Maximum Concentration

Time at which Cmax occurs:

```julia
# From NCA result
println("Tmax: $(result.tmax)")

# Direct calculation
tmax = nca_tmax(times, conc)
```

When multiple points have the same maximum concentration, Tmax is the **first** occurrence.

### Cmin - Minimum Concentration

For multiple dose/steady-state analysis:

```julia
# Multiple dose NCA
result = run_nca(times, conc, 100.0; dosing_type=:multiple, tau=24.0)
println("Cmin: $(result.cmin)")

# Direct calculation
cmin = nca_cmin(conc)
```

### Clast and Tlast - Last Measurable Concentration

Last concentration above the LLOQ:

```julia
# From NCA result
println("Clast: $(result.clast)")
println("Tlast: $(result.tlast)")

# With LLOQ specification
config = NCAConfig(lloq=0.05)
result = run_nca(times, conc, 100.0; config=config)
```

### Cavg - Average Concentration

Average concentration over a dosing interval (AUC0-tau/tau):

```julia
result = run_nca(times, conc, 100.0; dosing_type=:steady_state, tau=24.0)
println("Cavg: $(result.cavg)")

# Direct calculation
cavg = nca_cavg(times, conc, 24.0, NCAConfig())
```

---

## AUC Calculations

### AUC Methods

NeoPKPD supports three AUC calculation methods:

| Method | Usage | Formula |
|--------|-------|---------|
| `LinearMethod()` | Ascending phases | $(C_1 + C_2) \cdot \Delta t / 2$ |
| `LogLinearMethod()` | Descending phases | $(C_1 - C_2) / \ln(C_1/C_2) \cdot \Delta t$ |
| `LinLogMixedMethod()` | Recommended default | Linear up, log-linear down |

```julia
# Linear trapezoidal
config_linear = NCAConfig(method=LinearMethod())

# Log-linear trapezoidal
config_log = NCAConfig(method=LogLinearMethod())

# Lin-Log Mixed (FDA/EMA recommended)
config_mixed = NCAConfig(method=LinLogMixedMethod())

result = run_nca(times, conc, 100.0; config=config_mixed)
```

### AUC0-t (AUC to Last)

Area under the curve from time 0 to the last measurable concentration:

```julia
result = run_nca(times, conc, 100.0)
println("AUC0-t: $(result.auc_0_t)")

# Direct calculation
auc_0t = auc_0_t(times, conc, NCAConfig())
```

### AUC0-inf (AUC to Infinity)

AUC extrapolated to infinity using lambda_z:

$$AUC_{0-\infty} = AUC_{0-t} + \frac{C_{last}}{\lambda_z}$$

```julia
println("AUC0-inf: $(result.auc_0_inf)")
println("Extrapolation: $(result.auc_extra_pct)%")

# Direct calculation (requires lambda_z)
lambda_z = result.lambda_z_result.lambda_z
auc_inf, extra_pct = auc_0_inf(times, conc, lambda_z, result.clast, NCAConfig())
```

**Quality Warning**: If `auc_extra_pct > 20%`, the AUC0-inf may be unreliable.

### AUC0-tau (AUC over Dosing Interval)

For multiple dose or steady-state analysis:

```julia
result = run_nca(times, conc, 100.0; dosing_type=:multiple, tau=24.0)
println("AUC0-tau: $(result.auc_0_tau)")

# Direct calculation
auc_tau = auc_0_tau(times, conc, 24.0, NCAConfig())
```

### Partial AUC

AUC between any two time points:

```julia
# AUC from t=0 to t=4
auc_0_4 = auc_partial(times, conc, 0.0, 4.0, NCAConfig())

# AUC from t=4 to t=12
auc_4_12 = auc_partial(times, conc, 4.0, 12.0, NCAConfig())
```

Useful for assessing early exposure (e.g., AUC0-4h for rapid-acting drugs).

---

## AUMC (Area Under Moment Curve)

First moment curve for MRT calculation:

$$AUMC = \int_0^{t_{last}} t \cdot C(t) \, dt$$

```julia
println("AUMC0-t: $(result.aumc_0_t)")
println("AUMC0-inf: $(result.aumc_0_inf)")

# Direct calculation
aumc_0t = aumc_0_t(times, conc, NCAConfig())
```

---

## Dose-Normalized Metrics

For dose proportionality assessment:

```julia
# Dose-normalized Cmax
cmax_dn = nca_dose_normalized_cmax(result.cmax, 100.0)

# Dose-normalized AUC
auc_dn = nca_dose_normalized_auc(result.auc_0_inf, 100.0)

# From NCA result
println("Cmax/D: $(result.cmax_dn)")
println("AUC/D: $(result.auc_dn)")
```

---

## Time Above Concentration

Time above a specified concentration threshold:

```julia
# Time above MIC of 1.0 mg/L
t_above_mic = time_above_concentration(times, conc, 1.0)
println("Time above MIC: $t_above_mic h")
```

---

## Concentration at Specific Time

Interpolated concentration at any time point:

```julia
# Concentration at t=3h (interpolated)
c_at_3h = nca_c_at_time(times, conc, 3.0)

# Concentration at t=6h
c_at_6h = nca_c_at_time(times, conc, 6.0)
```

Uses linear interpolation for ascending and log-linear for descending phases.

---

## C0 - Initial Concentration (IV)

### Back-Extrapolation

For IV bolus, estimate C0 by back-extrapolation:

```julia
# Back-extrapolation from first two points
c0 = nca_c0_backextrap(times, conc)

# From terminal phase regression
c0_reg = nca_c0_from_regression(times, conc, lambda_z, intercept)
```

### Validation

```julia
validation = validate_c0_extrapolation(times, conc, c0)
println("Extrapolation valid: $(validation.is_valid)")
println("Warnings: $(validation.warnings)")
```

---

## BLQ Handling

Below Limit of Quantification handling:

```julia
# Treat BLQ as zero
config = NCAConfig(blq_handling=BLQZero(), lloq=0.05)

# Treat BLQ as LLOQ/2
config = NCAConfig(blq_handling=BLQLLOQHalf(), lloq=0.05)

# Exclude BLQ from calculations
config = NCAConfig(blq_handling=BLQMissing(), lloq=0.05)
```

---

## Example: Complete Exposure Analysis

```julia
using NeoPKPD

# PK data
times = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0]
conc = [0.0, 2.5, 4.8, 5.2, 4.5, 3.8, 2.6, 1.9, 1.0, 0.55, 0.18, 0.02]
dose = 500.0  # mg

# Configure NCA
config = NCAConfig(
    method = LinLogMixedMethod(),
    lambda_z_min_points = 3,
    lambda_z_r2_threshold = 0.9,
    extrapolation_max_pct = 20.0,
    lloq = 0.01
)

# Run analysis
result = run_nca(times, conc, dose; config=config, route=:extravascular)

# Report primary exposure
println("=== Primary Exposure Metrics ===")
println("Cmax:     $(round(result.cmax, digits=2)) mg/L")
println("Tmax:     $(result.tmax) h")
println("AUC0-t:   $(round(result.auc_0_t, digits=2)) mg·h/L")
println("AUC0-inf: $(round(result.auc_0_inf, digits=2)) mg·h/L")
println("AUC extrapolated: $(round(result.auc_extra_pct, digits=1))%")

# Report PK parameters
println("\n=== PK Parameters ===")
println("t1/2:  $(round(result.t_half, digits=2)) h")
println("CL/F:  $(round(result.cl_f, digits=2)) L/h")
println("Vz/F:  $(round(result.vz_f, digits=1)) L")
println("MRT:   $(round(result.mrt, digits=2)) h")

# Quality check
println("\n=== Quality Assessment ===")
println("Lambda_z R²: $(round(result.lambda_z_result.r_squared, digits=4))")
println("Points used: $(result.lambda_z_result.n_points)")
println("Quality flags: $(result.quality_flags)")
if !isempty(result.warnings)
    println("Warnings: $(result.warnings)")
end
```

---

## Formulas Summary

| Metric | Formula |
|--------|---------|
| Cmax | $\max(C)$ |
| AUC (linear) | $\sum \frac{(C_i + C_{i+1})}{2} \cdot (t_{i+1} - t_i)$ |
| AUC (log-linear) | $\sum \frac{(C_i - C_{i+1})}{\ln(C_i/C_{i+1})} \cdot (t_{i+1} - t_i)$ |
| AUC0-inf | $AUC_{0-t} + C_{last}/\lambda_z$ |
| AUC%extrap | $100 \cdot (C_{last}/\lambda_z) / AUC_{0-\infty}$ |
| Cavg | $AUC_{0-\tau} / \tau$ |
| Cmax/D | $C_{max} / Dose$ |
| AUC/D | $AUC / Dose$ |

---

## See Also

- [Terminal Phase Analysis](terminal-phase.md) - Lambda_z estimation
- [Bioequivalence](bioequivalence.md) - BE analysis
- [Multiple Dose](multiple-dose.md) - Steady-state metrics
- [Population NCA](population-nca.md) - Multi-subject analysis

