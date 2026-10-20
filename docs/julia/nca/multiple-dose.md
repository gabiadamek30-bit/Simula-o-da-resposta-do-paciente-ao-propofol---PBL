# Multiple Dose NCA

Comprehensive documentation for non-compartmental analysis of multiple dose and steady-state pharmacokinetic data.

---

## Overview

Multiple dose NCA extends single-dose analysis to repeated dosing scenarios, providing additional metrics relevant to chronic therapy and steady-state characterization.

---

## Quick Start

```julia
using NeoPKPD

# Steady-state data (after multiple doses)
times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]  # Within dosing interval
conc = [2.1, 5.8, 5.2, 4.1, 2.8, 1.9, 1.3]    # Steady-state concentrations
dose = 100.0
tau = 12.0  # Dosing interval (hours)

# Run steady-state NCA
result = run_nca(times, conc, dose; dosing_type=:steady_state, tau=tau)

println("Cmax,ss: $(result.cmax)")
println("Cmin,ss: $(result.cmin)")
println("Cavg,ss: $(result.cavg)")
println("AUC0-tau: $(result.auc_0_tau)")
println("Fluctuation: $(result.fluctuation)%")
```

---

## Dosing Types

### Single Dose

```julia
# First dose analysis (default)
result = run_nca(times, conc, dose; dosing_type=:single)
```

### Multiple Dose (Non-Steady-State)

```julia
# Analysis during accumulation phase
result = run_nca(times, conc, dose; dosing_type=:multiple, tau=12.0)
```

### Steady State

```julia
# Analysis at steady state
result = run_nca(times, conc, dose; dosing_type=:steady_state, tau=12.0)
```

---

## Steady-State Metrics

### Cmax,ss and Cmin,ss

Peak and trough concentrations at steady state:

```julia
result = run_nca(times, conc, dose; dosing_type=:steady_state, tau=12.0)

# Maximum concentration at steady state
println("Cmax,ss: $(result.cmax)")

# Minimum concentration (trough) at steady state
println("Cmin,ss: $(result.cmin)")
println("Cmin at time: $(result.tmin)")
```

### Cavg (Average Concentration)

$$C_{avg} = \frac{AUC_{0-\tau}}{\tau}$$

```julia
# Average concentration over dosing interval
println("Cavg: $(result.cavg)")

# Direct calculation
cavg = result.auc_0_tau / tau
```

### AUC0-tau

Area under the curve over one dosing interval:

```julia
# AUC over dosing interval
println("AUC0-tau: $(result.auc_0_tau)")

# Relates to total exposure per dose at steady state
```

---

## Fluctuation Metrics

### Peak-Trough Fluctuation (PTF)

$$PTF = \frac{C_{max} - C_{min}}{C_{avg}} \times 100\%$$

```julia
result = run_nca(times, conc, dose; dosing_type=:steady_state, tau=12.0)

println("PTF: $(result.fluctuation)%")

# Manual calculation
ptf = (result.cmax - result.cmin) / result.cavg * 100
```

### Swing

$$Swing = \frac{C_{max} - C_{min}}{C_{min}} \times 100\%$$

```julia
println("Swing: $(result.swing)%")

# Manual calculation
swing = (result.cmax - result.cmin) / result.cmin * 100
```

### Interpretation

| PTF | Interpretation |
|-----|----------------|
| < 100% | Low fluctuation (extended release) |
| 100-200% | Moderate fluctuation (typical IR) |
| > 200% | High fluctuation (may need dosing adjustment) |

---

## Accumulation Metrics

### Accumulation Index (Racc)

Ratio of steady-state to first-dose exposure:

$$R_{acc} = \frac{AUC_{0-\tau,ss}}{AUC_{0-\tau,sd}}$$

```julia
# From steady-state result
println("Accumulation Index: $(result.accumulation_index)")

# Theoretical accumulation (from t1/2)
theoretical_racc = 1 / (1 - exp(-log(2) * tau / result.t_half))
```

### Observed vs Theoretical Accumulation

```julia
# Compare observed accumulation to theoretical
observed_racc = result.accumulation_index
theoretical_racc = nca_theoretical_accumulation(result.t_half, tau)

ratio = observed_racc / theoretical_racc
if abs(ratio - 1.0) < 0.2
    println("Accumulation consistent with linear kinetics")
else
    println("Possible time-dependent kinetics (ratio: $(ratio))")
end
```

### Time to Steady State

Approximate time to reach steady state (~5 half-lives):

```julia
t_ss = 5 * result.t_half
println("Time to steady state: ~$t_ss hours")

# More precise: time to 90% steady state
t_90_ss = log(10) / log(2) * result.t_half  # ~3.3 × t1/2
println("Time to 90% SS: ~$t_90_ss hours")
```

---

## Multiple Dose Analysis

### Non-Steady-State Multiple Dose

When steady state has not yet been achieved:

```julia
# Day 3 data (before steady state)
times_d3 = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]
conc_d3 = [1.5, 4.2, 4.0, 3.2, 2.1, 1.4, 0.9]

result = run_nca(
    times_d3, conc_d3, dose;
    dosing_type = :multiple,
    tau = 12.0,
    dose_number = 5  # 5th dose
)

# Compare to single dose
sd_result = run_nca(times_sd, conc_sd, dose)
accumulation_ratio = result.auc_0_tau / sd_result.auc_0_tau
println("Accumulation after 5 doses: $(accumulation_ratio)")
```

### Predicting Steady State from Single Dose

```julia
# Single dose analysis
sd_result = run_nca(times, conc, dose)

# Predict steady-state metrics
predicted_ss = predict_steady_state(sd_result, tau=12.0)

println("Predicted Cmax,ss: $(predicted_ss.cmax_ss)")
println("Predicted Cmin,ss: $(predicted_ss.cmin_ss)")
println("Predicted Cavg,ss: $(predicted_ss.cavg_ss)")
println("Predicted Racc: $(predicted_ss.accumulation_index)")
```

---

## Superposition Principle

For linear PK, predict multiple dose profiles from single dose:

```julia
# Single dose data
sd_times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
sd_conc = [0.0, 2.1, 2.5, 2.0, 1.2, 0.5, 0.2, 0.02]

# Predict steady-state profile using superposition
ss_times, ss_conc = superposition_predict(
    sd_times, sd_conc;
    tau = 12.0,
    n_doses = 10  # Number of doses to reach SS
)

# Run NCA on predicted SS profile
ss_result = run_nca(ss_times, ss_conc, dose; dosing_type=:steady_state, tau=12.0)
```

---

## Trough Concentration Analysis

### Ctrough Timing

```julia
# Ensure trough is at end of interval
result = run_nca(times, conc, dose; dosing_type=:steady_state, tau=12.0)

# Verify trough timing
if result.tmin ≈ tau
    println("Trough at expected time (end of interval)")
else
    println("Trough at $(result.tmin)h, expected at $(tau)h")
end
```

### Ctrough Target

For drugs with therapeutic targets:

```julia
# Check if trough meets target
target_ctrough = 1.0  # mg/L (e.g., MIC for antibiotics)

if result.cmin >= target_ctrough
    println("Trough concentration meets target")
else
    # Calculate dose adjustment needed
    dose_factor = target_ctrough / result.cmin
    new_dose = dose * dose_factor
    println("Consider dose increase to $(new_dose) mg")
end
```

---

## Time Above Threshold

### Time Above MIC (Antibiotics)

```julia
# Calculate %T>MIC over dosing interval
mic = 0.5  # mg/L
t_above_mic = time_above_concentration(times, conc, mic)
pct_above_mic = t_above_mic / tau * 100

println("Time above MIC: $t_above_mic h ($(pct_above_mic)% of interval)")

# Target interpretation (depends on antibiotic class)
# Time-dependent: Need >40-50% of interval
# Concentration-dependent: Cmax/MIC > 8-10
```

### AUC/MIC Ratio

```julia
auc_mic_ratio = result.auc_0_tau / mic
println("AUC0-tau/MIC: $auc_mic_ratio")

# Target varies by drug class
# Fluoroquinolones: typically >100-125
```

---

## Effective Half-Life

At steady state, the effective half-life describes drug elimination:

```julia
# From accumulation index
t_half_eff = nca_effective_half_life(result.accumulation_index, tau)
println("Effective t1/2: $t_half_eff h")

# Compare to terminal half-life
println("Terminal t1/2: $(result.t_half) h")

# Difference indicates multi-compartment kinetics
if t_half_eff < result.t_half * 0.8
    println("Effective t1/2 shorter than terminal - multi-compartment behavior")
end
```

---

## Example: Complete Multiple Dose Analysis

```julia
using NeoPKPD

# Drug administered Q12H at steady state
times = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
conc = [2.8, 8.5, 9.2, 7.5, 5.1, 3.8, 3.0, 2.5, 2.2]
dose = 250.0  # mg
tau = 12.0    # hours

# Configure NCA
config = NCAConfig(
    method = LinLogMixedMethod(),
    lambda_z_min_points = 3
)

# Run steady-state NCA
result = run_nca(
    times, conc, dose;
    config = config,
    dosing_type = :steady_state,
    tau = tau,
    route = :extravascular
)

# Report steady-state metrics
println("=== Steady-State PK Parameters ===")
println("Dose: $dose mg Q$(Int(tau))H")
println("")
println("Exposure Metrics:")
println("  Cmax,ss:    $(round(result.cmax, digits=2)) mg/L at $(result.tmax)h")
println("  Cmin,ss:    $(round(result.cmin, digits=2)) mg/L")
println("  Cavg,ss:    $(round(result.cavg, digits=2)) mg/L")
println("  AUC0-tau:   $(round(result.auc_0_tau, digits=2)) mg·h/L")
println("")
println("Variability Metrics:")
println("  PTF:        $(round(result.fluctuation, digits=1))%")
println("  Swing:      $(round(result.swing, digits=1))%")
println("")
println("Accumulation:")
println("  Racc:       $(round(result.accumulation_index, digits=2))")
println("  Effective t1/2: $(round(nca_effective_half_life(result.accumulation_index, tau), digits=2)) h")
println("  Terminal t1/2:  $(round(result.t_half, digits=2)) h")
println("")
println("Clearance:")
println("  CLss/F:     $(round(result.cl_f, digits=2)) L/h")

# Target assessment (e.g., for antibiotics)
mic = 0.5  # mg/L
println("\n=== Target Assessment (MIC = $mic mg/L) ===")
println("  Cmax/MIC:   $(round(result.cmax / mic, digits=1))")
println("  Cmin/MIC:   $(round(result.cmin / mic, digits=1))")
println("  AUC/MIC:    $(round(result.auc_0_tau / mic, digits=1))")

t_above = time_above_concentration(times, conc, mic)
println("  %T>MIC:     $(round(t_above / tau * 100, digits=1))%")

# Dosing recommendation
if result.cmin < mic
    println("\n⚠ Trough below MIC - consider dose adjustment")
    suggested_dose = dose * (mic / result.cmin) * 1.1  # 10% margin
    println("  Suggested dose: $(round(suggested_dose, digits=0)) mg Q$(Int(tau))H")
else
    println("\n✓ Trough concentration adequate")
end
```

---

## Formulas Summary

| Parameter | Formula |
|-----------|---------|
| Cavg | $AUC_{0-\tau} / \tau$ |
| PTF | $(C_{max} - C_{min}) / C_{avg} \times 100\%$ |
| Swing | $(C_{max} - C_{min}) / C_{min} \times 100\%$ |
| Racc | $AUC_{0-\tau,ss} / AUC_{0-\tau,sd}$ |
| Racc (theoretical) | $1 / (1 - e^{-\lambda_z \cdot \tau})$ |
| t1/2,eff | $\ln(2) \cdot \tau / \ln(R_{acc} / (R_{acc} - 1))$ |

---

## See Also

- [Exposure Metrics](exposure-metrics.md) - Single dose metrics
- [Terminal Phase](terminal-phase.md) - Lambda_z and t1/2
- [Population NCA](population-nca.md) - Multi-subject analysis
- [Bioequivalence](bioequivalence.md) - Steady-state BE studies

