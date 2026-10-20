# Prediction-Corrected VPC (pcVPC)

Comprehensive guide to prediction-corrected Visual Predictive Checks for models with variable dosing or covariates.

---

## Overview

Prediction-Corrected VPC (pcVPC) normalizes both observed and simulated data by the population prediction, removing structural model trends and enabling meaningful comparisons when:

- Doses vary between subjects
- Covariates affect predictions
- Sampling designs differ across subjects

### The Problem with Standard VPC

Standard VPC assumes homogeneous design across subjects. When doses or covariates vary:

```
Subject 1: 100 mg dose → Expected Cmax ≈ 10 mg/L
Subject 2: 200 mg dose → Expected Cmax ≈ 20 mg/L

Standard VPC bins mix these, creating artificial variability
```

### pcVPC Solution

pcVPC corrects each observation by the population prediction:

$$pcDV_{ij} = DV_{ij} \times \frac{PRED_{bin}}{PRED_{ij}}$$

Where:
- $DV_{ij}$ = Observed concentration for subject $i$ at time $j$
- $PRED_{ij}$ = Population prediction for subject $i$ at time $j$
- $PRED_{bin}$ = Median population prediction for the time bin

---

## Mathematical Foundation

### Prediction Correction Formula

For each observation:

1. **Compute individual PRED**: $PRED_{ij}$ from population model (no IIV)
2. **Compute bin PRED**: $PRED_{bin} = \text{median}(PRED_{ij})$ for all $j$ in bin
3. **Correct observation**: $pcDV_{ij} = DV_{ij} \times \frac{PRED_{bin}}{PRED_{ij}}$

### Effect of Correction

| Scenario | Before Correction | After Correction |
|----------|-------------------|------------------|
| High dose subject | High DV | Normalized to bin median |
| Low dose subject | Low DV | Normalized to bin median |
| Heavy patient (CL↑) | Low DV | Normalized |
| Light patient (CL↓) | High DV | Normalized |

---

## Computing pcVPC

### Basic Usage

```julia
using NeoPKPD

# Observed data with variable dosing
observed = ObservedData(
    subject_ids = subject_ids,
    times = obs_times,
    dv = obs_dv,
    dvid = fill(:conc, length(obs_dv))
)

# Population model
typical_params = OneCompOralParams(1.5, 5.0, 50.0)
omega = OmegaMatrix([
    0.09 0.0  0.0;
    0.0  0.09 0.0;
    0.0  0.0  0.04
])

# Variable doses per subject
doses_per_subject = Dict(
    "S1" => [DoseEvent(0.0, 100.0)],
    "S2" => [DoseEvent(0.0, 200.0)],
    "S3" => [DoseEvent(0.0, 150.0)],
    # ...
)

base_spec = ModelSpec(OneCompOral(), "pcvpc_model", typical_params, doses_per_subject)

pop_spec = PopulationSpec(
    base_spec,
    n = length(unique(observed.subject_ids)),
    omega = omega,
    seed = 12345
)

# VPC configuration with prediction correction
config = VPCConfig(
    pi_levels = [0.05, 0.50, 0.95],
    prediction_corrected = true,    # Enable pcVPC
    binning = QuantileBinning(8),
    n_simulations = 500,
    seed = 42
)

grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Compute pcVPC
result = compute_pcvpc(observed, pop_spec, grid, solver; config=config)
```

### Dedicated pcVPC Function

```julia
# Using dedicated function (equivalent to above)
result = compute_pcvpc(
    observed,
    pop_spec,
    grid,
    solver;
    config = VPCConfig(
        pi_levels = [0.05, 0.50, 0.95],
        n_simulations = 500,
        binning = QuantileBinning(8)
    )
)
```

---

## When to Use pcVPC

### Use pcVPC When

1. **Variable dosing across subjects**
   ```julia
   # Different doses
   doses = [100, 200, 300, 400, 500]  # mg
   ```

2. **Significant covariate effects**
   ```julia
   # Weight-based dosing
   covariate_model = CovariateModel([
       CovariateEffect(:CL, :WT, 70.0, :power, 0.75),
       CovariateEffect(:V, :WT, 70.0, :power, 1.0)
   ])
   ```

3. **Multiple formulations or routes**
   ```julia
   # Different bioavailability
   formulations = [:tablet, :capsule, :solution]
   ```

4. **Dose escalation studies**
   ```julia
   # Phase 1 with escalating doses
   cohorts = [10, 30, 100, 300, 1000]  # mg
   ```

### Use Standard VPC When

- All subjects receive same dose
- No significant covariate effects
- Homogeneous study design

---

## pcVPC Internals

### Prediction Computation

```julia
# Internal: Compute base model predictions
function _compute_prediction_correction_data(
    observed::ObservedData,
    pop_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    bin_defs::Vector{BinDefinition}
)
    # 1. Simulate population with IIV=0 (typical subject)
    # 2. Interpolate to exact observation times
    # 3. Compute bin median predictions
    return (obs_preds, bin_pred_medians)
end
```

### Correction Application

```julia
# Internal: Apply correction to observations
function _apply_prediction_correction(
    obs_values::Vector{Float64},
    obs_times::Vector{Float64},
    obs_preds::Vector{Float64},
    bin_defs::Vector{BinDefinition},
    bin_pred_medians::Dict{Int, Float64}
)
    corrected = similar(obs_values)
    for i in eachindex(obs_values)
        bin_id = _find_bin(obs_times[i], bin_defs)
        if bin_id !== nothing
            corrected[i] = obs_values[i] * bin_pred_medians[bin_id] / obs_preds[i]
        end
    end
    return corrected
end
```

---

## Covariate Effects in pcVPC

### With Covariate Model

```julia
# Define covariate effects
covariate_model = CovariateModel([
    CovariateEffect(:CL, :WT, 70.0, :power, 0.75),   # Allometric
    CovariateEffect(:V, :WT, 70.0, :power, 1.0),     # Allometric
    CovariateEffect(:CL, :CRCL, 100.0, :linear, 0.5) # Renal function
])

# Subject covariates
covariates = [
    Dict(:WT => 65.0, :CRCL => 90.0),
    Dict(:WT => 85.0, :CRCL => 120.0),
    Dict(:WT => 55.0, :CRCL => 60.0),
    # ...
]

# Population with covariates
pop_spec = PopulationSpec(
    base_spec,
    n = n_subjects,
    omega = omega,
    covariate_model = covariate_model,
    covariates = covariates,
    seed = 12345
)

# pcVPC automatically handles covariate-adjusted predictions
result = compute_pcvpc(observed, pop_spec, grid, solver; config=config)
```

### Covariate Impact on PRED

The prediction correction accounts for:

$$PRED_i = f(\theta, X_i)$$

Where $X_i$ are individual covariates. Heavy patients have higher PRED due to:

```julia
# For CL with power covariate
CL_i = CL_pop * (WT_i / 70)^0.75

# PRED is computed using this adjusted CL
```

---

## Result Structure

pcVPC results have the same structure as standard VPC:

```julia
struct VPCResult
    config::VPCConfig
    bins::Vector{VPCBin}
    n_subjects_observed::Int
    n_observations_observed::Int
    n_simulations::Int
    strata::String
    simulation_seed::UInt64
end
```

### Accessing pcVPC Results

```julia
# Same accessor functions work
times = bin_midpoints(result)

obs_median = observed_percentile(result, 0.50)
sim_median = simulated_median(result, 0.50)
sim_lower = simulated_lower(result, 0.50)
sim_upper = simulated_upper(result, 0.50)

# These are now prediction-corrected values
println("Prediction-corrected observed median: ", obs_median)
println("Prediction-corrected simulated CI: [$sim_lower, $sim_upper]")
```

---

## Interpretation

### pcVPC vs Standard VPC

| Aspect | Standard VPC | pcVPC |
|--------|--------------|-------|
| Y-axis | Concentration | Prediction-corrected concentration |
| Variability source | Structural + IIV | IIV only |
| Dose differences | Visible as variability | Normalized out |
| Covariate effects | Visible as variability | Normalized out |
| Model assessment | Overall fit | Random effects fit |

### What pcVPC Shows

1. **Residual variability capture** - Is omega correctly estimated?
2. **IIV distribution shape** - Log-normal assumption valid?
3. **Time-varying effects** - Any model misspecification over time?

### What pcVPC Hides

1. **Dose-response relationship** - Already normalized
2. **Covariate relationships** - Already corrected
3. **Structural model trends** - Removed by correction

---

## Complete Example

```julia
using NeoPKPD
using Random

# ================================================
# pcVPC Example: Variable Dosing Study
# ================================================

println("=== Prediction-Corrected VPC ===\n")

# 1. Simulate variable-dose study
Random.seed!(456)
n_subjects = 60

# Dose levels (mg)
dose_levels = [50.0, 100.0, 200.0]
n_per_dose = 20

subject_ids = String[]
obs_times = Float64[]
obs_dv = Float64[]
doses_dict = Dict{String, Vector{DoseEvent}}()

# True parameters
true_ka = 1.5
true_cl = 5.0
true_v = 50.0
omega_ka = 0.16
omega_cl = 0.09
omega_v = 0.04

sampling_times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]

subject_idx = 1
for dose in dose_levels
    for _ in 1:n_per_dose
        sid = "S$subject_idx"
        doses_dict[sid] = [DoseEvent(0.0, dose)]

        # Individual parameters
        ka_i = true_ka * exp(randn() * sqrt(omega_ka))
        cl_i = true_cl * exp(randn() * sqrt(omega_cl))
        v_i = true_v * exp(randn() * sqrt(omega_v))

        for t in sampling_times
            # One-compartment oral
            conc = dose * ka_i / (v_i * (ka_i - cl_i/v_i)) *
                   (exp(-cl_i/v_i * t) - exp(-ka_i * t))

            # Proportional error
            conc_obs = conc * (1 + 0.1 * randn())
            conc_obs = max(0.01, conc_obs)

            push!(subject_ids, sid)
            push!(obs_times, t)
            push!(obs_dv, conc_obs)
        end

        subject_idx += 1
    end
end

# 2. Create observed data
observed = ObservedData(
    subject_ids = subject_ids,
    times = obs_times,
    dv = obs_dv,
    dvid = fill(:conc, length(obs_dv))
)

println("Variable-dose study:")
println("  Dose levels: $dose_levels mg")
println("  Subjects per dose: $n_per_dose")
println("  Total subjects: $n_subjects")
println("  Total observations: $(length(obs_dv))")

# 3. Define population model
typical_params = OneCompOralParams(true_ka, true_cl, true_v)
omega = OmegaMatrix([
    omega_ka 0.0      0.0;
    0.0      omega_cl 0.0;
    0.0      0.0      omega_v
])

# Need subject-specific doses for pcVPC
base_spec = ModelSpec(OneCompOral(), "variable_dose_model", typical_params, doses_dict)

pop_spec = PopulationSpec(
    base_spec,
    n = n_subjects,
    omega = omega,
    seed = 12345
)

grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# 4. Compare Standard VPC vs pcVPC
println("\n--- Computing Standard VPC ---")
config_standard = VPCConfig(
    pi_levels = [0.05, 0.50, 0.95],
    prediction_corrected = false,
    binning = QuantileBinning(7),
    n_simulations = 500,
    seed = 42
)

result_standard = compute_vpc(observed, pop_spec, grid, solver; config=config_standard)

println("\n--- Computing pcVPC ---")
config_pc = VPCConfig(
    pi_levels = [0.05, 0.50, 0.95],
    prediction_corrected = true,
    binning = QuantileBinning(7),
    n_simulations = 500,
    seed = 42
)

result_pcvpc = compute_pcvpc(observed, pop_spec, grid, solver; config=config_pc)

# 5. Compare results
println("\n--- Comparison: Standard vs pcVPC ---")
println("\nStandard VPC (affected by dose variability):")
println("Bin | Obs P50 Range | Sim P50 CI Width")
for bin in result_standard.bins
    p50 = filter(p -> p.percentile == 0.50, bin.percentiles)[1]
    ci_width = p50.simulated_upper - p50.simulated_lower
    println("$(bin.bin_id)   | $(round(p50.observed, digits=2)) | $(round(ci_width, digits=2))")
end

println("\npcVPC (normalized for dose):")
println("Bin | Obs P50 Range | Sim P50 CI Width")
for bin in result_pcvpc.bins
    p50 = filter(p -> p.percentile == 0.50, bin.percentiles)[1]
    ci_width = p50.simulated_upper - p50.simulated_lower
    println("$(bin.bin_id)   | $(round(p50.observed, digits=2)) | $(round(ci_width, digits=2))")
end

# 6. Coverage comparison
println("\n--- Coverage Comparison ---")
for level in [0.05, 0.50, 0.95]
    cov_std = vpc_coverage(result_standard, level)
    cov_pc = vpc_coverage(result_pcvpc, level)
    println("P$(Int(level*100)): Standard=$(round(cov_std*100, digits=1))%, pcVPC=$(round(cov_pc*100, digits=1))%")
end

println("\n✓ pcVPC computation complete")
println("\nNote: pcVPC shows tighter CI as dose variability is normalized")
```

---

## Best Practices

### When Computing pcVPC

1. **Ensure accurate PRED** - Use final model parameters
2. **Include all covariates** - That affect PRED
3. **Match dosing exactly** - Subject-specific doses
4. **Sufficient simulations** - 500+ for stable CI

### Reporting

When presenting pcVPC:

1. **Label axes clearly** - "Prediction-Corrected Concentration"
2. **Note correction method** - "pcVPC per Bergstrand et al. (2011)"
3. **Report uncorrected alongside** - When relevant
4. **Explain interpretation** - For non-expert audience

---

## References

- Bergstrand M, Hooker AC, Wallin JE, Karlsson MO. Prediction-corrected visual predictive checks for diagnosing nonlinear mixed-effects models. AAPS J. 2011;13(2):143-151.

---

## See Also

- [Standard VPC](standard.md) - Basic VPC methodology
- [Stratified VPC](stratified.md) - Covariate stratification
- [VPC Index](index.md) - Overview
- [Python pcVPC](../../python/viz/vpc.md) - Python visualization
