# Residual Error Models

Comprehensive guide for modeling unexplained variability in observations.

---

## Overview

Residual error models describe the discrepancy between model predictions and observed data, encompassing measurement error, model misspecification, and other unexplained variability.

```julia
using NeoPKPD

# Proportional error (10% CV)
error_model = ResidualErrorSpec(
    ProportionalError(),
    ProportionalErrorParams(0.10)
)
```

---

## Mathematical Foundation

### General Form

Observed concentration relates to predicted concentration:

$$Y_{ij} = F_{ij} + g(F_{ij}, \sigma) \cdot \varepsilon_{ij}$$

Where:
- $Y_{ij}$ = Observed concentration for subject $i$ at time $j$
- $F_{ij}$ = Model-predicted concentration (IPRED)
- $g()$ = Error magnitude function
- $\varepsilon_{ij}$ = Standard normal random variable, $\varepsilon \sim N(0, 1)$
- $\sigma$ = Error parameter(s)

---

## Error Model Types

### Additive Error

Constant error magnitude regardless of concentration:

$$Y = F + \sigma_{add} \cdot \varepsilon$$

```julia
# Additive error: SD = 0.5 mg/L
additive = ResidualErrorSpec(
    AdditiveError(),
    AdditiveErrorParams(0.5)
)

# Variance is constant: Var(Y) = σ²_add = 0.25
```

**When to use:**
- Low concentration range
- Assay precision limited by detection
- When CV increases at low concentrations

### Proportional Error

Error scales with prediction (constant CV):

$$Y = F \cdot (1 + \sigma_{prop} \cdot \varepsilon)$$

```julia
# Proportional error: CV = 15%
proportional = ResidualErrorSpec(
    ProportionalError(),
    ProportionalErrorParams(0.15)
)

# SD scales with prediction: SD(Y) = F × σ_prop
# CV is constant: CV(Y) = σ_prop = 15%
```

**When to use:**
- Wide concentration range
- Consistent assay CV across range
- Most common for PK data

### Combined Error

Both additive and proportional components:

$$Y = F + \sqrt{\sigma^2_{add} + (\sigma_{prop} \cdot F)^2} \cdot \varepsilon$$

```julia
# Combined: additive SD = 0.1, proportional CV = 10%
combined = ResidualErrorSpec(
    CombinedError(),
    CombinedErrorParams(0.1, 0.10)
)

# Variance: Var(Y) = σ²_add + (σ_prop × F)²
# CV decreases at higher concentrations
```

**When to use:**
- Wide dynamic range
- Higher CV at low concentrations
- Most flexible model

### Exponential (Log-Normal) Error

Error on log scale (multiplicative):

$$Y = F \cdot e^{\sigma_{exp} \cdot \varepsilon}$$

Or equivalently:
$$\log(Y) = \log(F) + \sigma_{exp} \cdot \varepsilon$$

```julia
# Exponential error: 20% CV approximately
exponential = ResidualErrorSpec(
    ExponentialError(),
    ExponentialErrorParams(0.20)
)

# Assumes log(Y) normally distributed
# Appropriate for log-transformed data
```

**When to use:**
- Log-transformed data analysis
- Right-skewed residuals on natural scale
- Common in bioequivalence

---

## Error Parameter Structures

### AdditiveErrorParams

```julia
struct AdditiveErrorParams
    sigma::Float64    # Standard deviation (same units as observation)
end

# Example: If concentration in mg/L
params = AdditiveErrorParams(0.5)  # SD = 0.5 mg/L
```

### ProportionalErrorParams

```julia
struct ProportionalErrorParams
    sigma::Float64    # Proportional error (as fraction, not %)
end

# Example: 15% CV
params = ProportionalErrorParams(0.15)  # σ = 0.15, not 15
```

### CombinedErrorParams

```julia
struct CombinedErrorParams
    sigma_add::Float64    # Additive SD
    sigma_prop::Float64   # Proportional fraction
end

# Example: SD = 0.1 mg/L + 10% CV
params = CombinedErrorParams(0.1, 0.10)
```

### ExponentialErrorParams

```julia
struct ExponentialErrorParams
    sigma::Float64    # SD on log scale
end

# Example: ~20% CV
params = ExponentialErrorParams(0.20)
```

---

## ResidualErrorSpec Structure

### Class Definition

```julia
struct ResidualErrorSpec{K<:ErrorModelKind, P<:ErrorParams}
    kind::K           # AdditiveError(), ProportionalError(), etc.
    params::P         # Error parameters
    seed::Int         # Random seed for reproducibility
end
```

### Creating Specifications

```julia
# Proportional with seed
error_spec = ResidualErrorSpec(
    ProportionalError(),
    ProportionalErrorParams(0.12),
    seed = 54321
)

# Combined without seed (random)
error_spec = ResidualErrorSpec(
    CombinedError(),
    CombinedErrorParams(0.05, 0.08)
)
```

---

## Computing Error Components

### Variance Functions

```julia
# Compute variance at a prediction
function residual_variance(
    f::Float64,
    params::AdditiveErrorParams
)
    return params.sigma^2
end

function residual_variance(
    f::Float64,
    params::ProportionalErrorParams
)
    return (params.sigma * f)^2
end

function residual_variance(
    f::Float64,
    params::CombinedErrorParams
)
    return params.sigma_add^2 + (params.sigma_prop * f)^2
end

function residual_variance(
    f::Float64,
    params::ExponentialErrorParams
)
    # Variance on natural scale (approximation)
    return (f * params.sigma)^2
end
```

### Standard Deviation Functions

```julia
function residual_sd(f::Float64, params)
    return sqrt(residual_variance(f, params))
end

# Example usage
f = 5.0  # Predicted concentration
params = CombinedErrorParams(0.1, 0.10)

sd = residual_sd(f, params)
# SD = √(0.1² + (0.10 × 5.0)²) = √(0.01 + 0.25) = √0.26 = 0.51
```

---

## Applying Residual Error

### Adding Error to Predictions

```julia
using Random

function apply_residual_error(
    predictions::Vector{Float64},
    error_spec::ResidualErrorSpec;
    seed::Int = nothing
)
    if seed !== nothing
        Random.seed!(seed)
    end

    n = length(predictions)
    observations = similar(predictions)
    epsilons = randn(n)

    for i in 1:n
        sd = residual_sd(predictions[i], error_spec.params)
        observations[i] = predictions[i] + sd * epsilons[i]
    end

    return observations
end

# Example
preds = [10.0, 8.0, 5.0, 2.0, 1.0]
error = ResidualErrorSpec(ProportionalError(), ProportionalErrorParams(0.10))

obs = apply_residual_error(preds, error, seed=42)
```

### Type-Specific Application

```julia
# Additive
function apply_error(f::Float64, eps::Float64, params::AdditiveErrorParams)
    return f + params.sigma * eps
end

# Proportional
function apply_error(f::Float64, eps::Float64, params::ProportionalErrorParams)
    return f * (1 + params.sigma * eps)
end

# Combined
function apply_error(f::Float64, eps::Float64, params::CombinedErrorParams)
    sd = sqrt(params.sigma_add^2 + (params.sigma_prop * f)^2)
    return f + sd * eps
end

# Exponential
function apply_error(f::Float64, eps::Float64, params::ExponentialErrorParams)
    return f * exp(params.sigma * eps)
end
```

---

## Weighted Residuals

### Individual Weighted Residuals (IWRES)

Standardized residuals using individual predictions:

$$IWRES_{ij} = \frac{Y_{ij} - IPRED_{ij}}{SD(Y|IPRED_{ij})}$$

```julia
function compute_iwres(
    observed::Float64,
    ipred::Float64,
    error_params
)
    sd = residual_sd(ipred, error_params)
    return (observed - ipred) / sd
end

# Should be approximately N(0,1) for good fit
```

### Conditional Weighted Residuals (CWRES)

More sophisticated residuals accounting for random effects:

```julia
function compute_cwres(
    observed::Vector{Float64},
    pred::Vector{Float64},      # Population predictions
    ipred::Vector{Float64},     # Individual predictions
    error_params,
    omega::OmegaMatrix
)
    n = length(observed)
    cwres = similar(observed)

    for i in 1:n
        # Approximation using FOCE method
        var_y = residual_variance(pred[i], error_params)
        # Additional variance from random effects
        # (simplified - full implementation more complex)
        cwres[i] = (observed[i] - pred[i]) / sqrt(var_y)
    end

    return cwres
end
```

---

## Model Diagnostics

### Residual Plots

```julia
# Key diagnostic plots for error model assessment:

# 1. IWRES vs Time
# - Should be randomly scattered around 0
# - Patterns indicate model misspecification

# 2. IWRES vs IPRED
# - Should be randomly scattered around 0
# - Funnel shape suggests wrong error model

# 3. |IWRES| vs IPRED
# - Should be flat
# - Increasing trend: need proportional component
# - Decreasing trend: proportional over-estimated

# 4. QQ plot of IWRES
# - Should follow 45° line
# - Deviations indicate non-normal residuals
```

### Choosing Error Model

```julia
# Compare error models using likelihood
function compare_error_models(
    data::PopulationData,
    base_model::PopulationModel,
    error_specs::Vector{ResidualErrorSpec}
)
    results = []

    for error in error_specs
        model = set_error_model(base_model, error)
        fit = fit_population(data, model)

        push!(results, (
            error_type = typeof(error.kind),
            n_params = count_error_params(error),
            ofv = fit.ofv,
            aic = fit.ofv + 2 * count_error_params(error)
        ))
    end

    # Sort by AIC
    sort!(results, by = r -> r.aic)

    return results
end
```

---

## Multi-Observation Error

### Different Error for Different Outputs

```julia
# PK and PD may have different error models
struct MultiErrorSpec
    errors::Dict{Symbol, ResidualErrorSpec}
end

multi_error = MultiErrorSpec(Dict(
    :conc => ResidualErrorSpec(ProportionalError(), ProportionalErrorParams(0.10)),
    :effect => ResidualErrorSpec(AdditiveError(), AdditiveErrorParams(5.0))
))
```

### Multiple Analytes

```julia
# Parent drug and metabolite
analyte_errors = Dict(
    :parent => ResidualErrorSpec(CombinedError(), CombinedErrorParams(0.05, 0.08)),
    :metabolite => ResidualErrorSpec(ProportionalError(), ProportionalErrorParams(0.15))
)
```

---

## Complete Example

```julia
using NeoPKPD
using Statistics
using Random

# ============================================
# Residual Error Model Demonstration
# ============================================

println("=== Residual Error Models ===\n")

# 1. Generate true PK profile
typical_cl = 10.0
typical_v = 50.0
dose = 500.0
times = collect(0.0:0.5:24.0)

# True concentrations (one-compartment IV bolus)
ke = typical_cl / typical_v
true_conc = (dose / typical_v) .* exp.(-ke .* times)

println("--- True Concentration Profile ---")
println("Cmax: $(round(maximum(true_conc), digits=2)) mg/L")
println("t1/2: $(round(log(2)/ke, digits=2)) hr")

# 2. Define different error models
error_models = [
    ("Additive (0.5 mg/L)", ResidualErrorSpec(AdditiveError(), AdditiveErrorParams(0.5))),
    ("Proportional (15%)", ResidualErrorSpec(ProportionalError(), ProportionalErrorParams(0.15))),
    ("Combined (0.2 + 10%)", ResidualErrorSpec(CombinedError(), CombinedErrorParams(0.2, 0.10))),
    ("Exponential (15%)", ResidualErrorSpec(ExponentialError(), ExponentialErrorParams(0.15)))
]

# 3. Simulate observations with each error model
Random.seed!(42)

println("\n--- Error Model Comparison ---")
println("Time points: $(length(times))")
println("")

for (name, error_spec) in error_models
    # Generate observations
    obs = apply_residual_error(true_conc, error_spec, seed=42)

    # Calculate residuals
    residuals = obs .- true_conc

    # Calculate metrics
    rmse = sqrt(mean(residuals.^2))
    mae = mean(abs.(residuals))

    # CV at different concentration levels
    high_idx = findall(true_conc .> 5.0)
    low_idx = findall(0.5 .< true_conc .<= 2.0)

    cv_high = std(obs[high_idx] .- true_conc[high_idx]) / mean(true_conc[high_idx]) * 100
    cv_low = std(obs[low_idx] .- true_conc[low_idx]) / mean(true_conc[low_idx]) * 100

    println("$name:")
    println("  RMSE: $(round(rmse, digits=3)) mg/L")
    println("  MAE:  $(round(mae, digits=3)) mg/L")
    println("  CV at high conc: $(round(cv_high, digits=1))%")
    println("  CV at low conc:  $(round(cv_low, digits=1))%")
    println("")
end

# 4. Detailed analysis of combined error
println("--- Combined Error Analysis ---")
combined = ResidualErrorSpec(CombinedError(), CombinedErrorParams(0.2, 0.10))

# Show how SD changes with concentration
println("SD by concentration:")
for conc in [10.0, 5.0, 2.0, 1.0, 0.5]
    sd = residual_sd(conc, combined.params)
    cv = sd / conc * 100
    println("  Conc = $conc: SD = $(round(sd, digits=3)), CV = $(round(cv, digits=1))%")
end

# 5. Compute IWRES
println("\n--- IWRES Analysis ---")
obs = apply_residual_error(true_conc, combined, seed=123)
iwres = [(obs[i] - true_conc[i]) / residual_sd(true_conc[i], combined.params) for i in 1:length(obs)]

println("IWRES statistics:")
println("  Mean: $(round(mean(iwres), digits=4)) (should be ~0)")
println("  SD:   $(round(std(iwres), digits=4)) (should be ~1)")
println("  Min:  $(round(minimum(iwres), digits=2))")
println("  Max:  $(round(maximum(iwres), digits=2))")

# Check normality
sorted_iwres = sort(iwres)
n = length(iwres)
expected_quantiles = [quantile(Normal(), (i - 0.5)/n) for i in 1:n]

# Correlation for QQ plot (should be close to 1)
qq_corr = cor(sorted_iwres, expected_quantiles)
println("  QQ correlation: $(round(qq_corr, digits=4)) (should be ~1)")

# 6. Population simulation with error
println("\n--- Population Simulation with Error ---")

# Simulate 100 subjects
n_subjects = 100
omega = OmegaMatrix([0.09 0.0; 0.0 0.04])  # IIV on CL and V

# Sample etas
Random.seed!(42)
eta_cl = randn(n_subjects) .* sqrt(0.09)
eta_v = randn(n_subjects) .* sqrt(0.04)

# Simulate each subject
all_obs = Matrix{Float64}(undef, n_subjects, length(times))
all_pred = Matrix{Float64}(undef, n_subjects, length(times))

for i in 1:n_subjects
    ind_cl = typical_cl * exp(eta_cl[i])
    ind_v = typical_v * exp(eta_v[i])
    ind_ke = ind_cl / ind_v

    # Individual predictions
    pred = (dose / ind_v) .* exp.(-ind_ke .* times)
    all_pred[i, :] = pred

    # Add residual error
    all_obs[i, :] = apply_residual_error(pred, combined)
end

# Summary statistics
mean_obs = vec(mean(all_obs, dims=1))
sd_obs = vec(std(all_obs, dims=1))
cv_obs = sd_obs ./ mean_obs .* 100

println("Observed concentration summary:")
println("  Mean Cmax: $(round(maximum(mean_obs), digits=2)) mg/L")
println("  CV at Cmax: $(round(cv_obs[1], digits=1))%")
println("  CV at t=12h: $(round(cv_obs[findfirst(==(12.0), times)], digits=1))%")

# Partition variance
println("\n--- Variance Components ---")
# Total variance = IIV variance + Residual variance
t_idx = 10  # Example time point
total_var = var(all_obs[:, t_idx])
pred_mean = mean(all_pred[:, t_idx])
iiv_var = var(all_pred[:, t_idx])
res_var = residual_variance(pred_mean, combined.params)

println("At t=$(times[t_idx]) hr:")
println("  Total variance: $(round(total_var, digits=3))")
println("  IIV variance:   $(round(iiv_var, digits=3)) ($(round(iiv_var/total_var*100, digits=1))%)")
println("  Residual var:   $(round(res_var, digits=3)) ($(round(res_var/total_var*100, digits=1))%)")
```

---

## Error Model Selection Guide

| Scenario | Recommended Model | Rationale |
|----------|-------------------|-----------|
| Wide conc range, constant CV | Proportional | Simplest when CV uniform |
| Low conc near LLOQ | Additive or Combined | Accounts for detection limit |
| Very wide range (>100-fold) | Combined | Flexible across all levels |
| Log-transformed analysis | Exponential | Matches log-normal assumption |
| PD response data | Often Additive | Response often has fixed precision |
| Count data | Special models | Poisson or negative binomial |

---

## See Also

- [IIV](iiv.md) - Inter-individual variability
- [IOV](iov.md) - Inter-occasion variability
- [Covariates](covariates.md) - Covariate effects
- [Estimation Diagnostics](../estimation/diagnostics.md) - Model evaluation
