# Covariate Models

Comprehensive guide for incorporating patient characteristics into population PK/PD models.

---

## Overview

Covariate models explain part of the inter-individual variability by relating parameters to measurable patient characteristics.

```julia
using NeoPKPD

# Allometric scaling: CL and V based on weight
covariate_model = CovariateModel([
    CovariateEffect(:CL, :WT, 70.0, PowerCovariate(), 0.75),
    CovariateEffect(:V, :WT, 70.0, PowerCovariate(), 1.0)
])
```

---

## Mathematical Foundation

### General Covariate Model

Parameters are adjusted based on covariate values:

$$\theta_i = \theta_{pop} \cdot f(cov_i, \beta)$$

Where:
- $\theta_i$ = Individual parameter
- $\theta_{pop}$ = Population typical value
- $f(cov_i, \beta)$ = Covariate function
- $\beta$ = Covariate coefficient

### Order of Application

1. **Covariates first**: Adjust typical value for covariates
2. **IIV second**: Apply random effects to covariate-adjusted value
3. **IOV third**: Apply occasion-specific effects

$$\theta_{ij} = \theta_{pop} \cdot f(cov_i) \cdot e^{\eta_i + \kappa_{ij}}$$

---

## Covariate Effect Types

### Power Model

Most common for body size effects:

$$\theta_i = \theta_{pop} \cdot \left(\frac{cov_i}{ref}\right)^{\beta}$$

```julia
# Allometric scaling
power_effect = CovariateEffect(
    :CL,                    # Target parameter
    :WT,                    # Covariate name
    70.0,                   # Reference value
    PowerCovariate(),       # Effect type
    0.75                    # Exponent (β)
)

# Common allometric exponents:
# Clearance: 0.75 (3/4 power)
# Volume: 1.0 (proportional to body size)
# Half-life: 0.25 (1/4 power)
```

### Linear Model

For effects proportional to deviation from reference:

$$\theta_i = \theta_{pop} \cdot (1 + \beta \cdot (cov_i - ref))$$

```julia
# Age effect on clearance
linear_effect = CovariateEffect(
    :CL,
    :AGE,
    45.0,                   # Reference age
    LinearCovariate(),
    -0.01                   # -1% per year above reference
)

# For a 65-year-old:
# CL = CL_pop * (1 + (-0.01) * (65 - 45))
# CL = CL_pop * (1 - 0.20) = 0.80 * CL_pop
```

### Exponential Model

For multiplicative effects:

$$\theta_i = \theta_{pop} \cdot e^{\beta \cdot (cov_i - ref)}$$

```julia
# Renal function effect
exp_effect = CovariateEffect(
    :CL,
    :CRCL,                  # Creatinine clearance
    100.0,                  # Reference (normal)
    ExpCovariate(),
    0.005                   # Effect coefficient
)

# For CRCL = 50 mL/min:
# CL = CL_pop * exp(0.005 * (50 - 100))
# CL = CL_pop * exp(-0.25) = 0.78 * CL_pop
```

---

## CovariateEffect Structure

### Class Definition

```julia
struct CovariateEffect{K<:CovariateKind}
    param::Symbol              # Target parameter (:CL, :V, etc.)
    covariate::Symbol          # Covariate name (:WT, :AGE, etc.)
    ref::Float64               # Reference value
    kind::K                    # PowerCovariate, LinearCovariate, ExpCovariate
    beta::Float64              # Effect coefficient
end
```

### Creating Effects

```julia
# Power effect
effect1 = CovariateEffect(:CL, :WT, 70.0, PowerCovariate(), 0.75)

# Linear effect
effect2 = CovariateEffect(:CL, :AGE, 45.0, LinearCovariate(), -0.008)

# Exponential effect
effect3 = CovariateEffect(:CLR, :CRCL, 100.0, ExpCovariate(), 0.007)
```

---

## CovariateModel Structure

### Class Definition

```julia
struct CovariateModel
    name::String
    effects::Vector{CovariateEffect}
end
```

### Creating Models

```julia
# Full covariate model
covariate_model = CovariateModel(
    "standard_allometry",
    [
        # Body size effects
        CovariateEffect(:CL, :WT, 70.0, PowerCovariate(), 0.75),
        CovariateEffect(:V1, :WT, 70.0, PowerCovariate(), 1.0),
        CovariateEffect(:Q, :WT, 70.0, PowerCovariate(), 0.75),
        CovariateEffect(:V2, :WT, 70.0, PowerCovariate(), 1.0),

        # Age effect
        CovariateEffect(:CL, :AGE, 45.0, LinearCovariate(), -0.005),

        # Renal function effect on renal clearance
        CovariateEffect(:CLR, :CRCL, 100.0, LinearCovariate(), 0.006)
    ]
)
```

### Applying Covariates

```julia
function apply_covariates(
    params::Dict{Symbol, Float64},
    model::CovariateModel,
    covariates::Dict{Symbol, Float64}
)
    adjusted = copy(params)

    for effect in model.effects
        if haskey(covariates, effect.covariate)
            cov_val = covariates[effect.covariate]
            adjusted[effect.param] = apply_effect(
                adjusted[effect.param],
                cov_val,
                effect
            )
        end
    end

    return adjusted
end

function apply_effect(
    param_val::Float64,
    cov_val::Float64,
    effect::CovariateEffect{PowerCovariate}
)
    return param_val * (cov_val / effect.ref) ^ effect.beta
end

function apply_effect(
    param_val::Float64,
    cov_val::Float64,
    effect::CovariateEffect{LinearCovariate}
)
    return param_val * (1 + effect.beta * (cov_val - effect.ref))
end

function apply_effect(
    param_val::Float64,
    cov_val::Float64,
    effect::CovariateEffect{ExpCovariate}
)
    return param_val * exp(effect.beta * (cov_val - effect.ref))
end
```

---

## IndividualCovariates Structure

### Class Definition

```julia
struct IndividualCovariates
    values::Dict{Symbol, Float64}                # Static covariates
    time_varying::Union{TimeVaryingCovariates, Nothing}  # Dynamic covariates
end
```

### Creating Individual Covariates

```julia
# Static covariates
subject_covs = IndividualCovariates(
    Dict(:WT => 80.0, :AGE => 55.0, :SEX => 1.0, :CRCL => 85.0),
    nothing
)

# With time-varying covariate
time_cov = TimeVaryingCovariates(Dict(
    :CRCL => TimeCovariateSeries(
        StepTimeCovariate(),
        [0.0, 24.0, 48.0, 72.0],
        [85.0, 80.0, 75.0, 78.0]
    )
))

subject_covs = IndividualCovariates(
    Dict(:WT => 80.0, :AGE => 55.0),
    time_cov
)
```

---

## Time-Varying Covariates

### Step Interpolation

Value held constant until next time point:

```julia
# Creatinine clearance declining stepwise
crcl_series = TimeCovariateSeries(
    StepTimeCovariate(),
    [0.0, 24.0, 48.0, 72.0, 96.0],    # Times
    [100.0, 90.0, 75.0, 80.0, 85.0]   # Values
)

# At t=30: CRCL = 90.0 (from t=24 value)
# At t=50: CRCL = 75.0 (from t=48 value)
```

### Linear Interpolation

Smooth transition between values:

```julia
# Body weight changing gradually
wt_series = TimeCovariateSeries(
    LinearTimeCovariate(),
    [0.0, 168.0, 336.0],              # Weekly
    [80.0, 78.0, 76.0]                # Decreasing
)

# At t=84 (halfway between 0 and 168):
# WT = 80.0 + (78.0 - 80.0) * (84/168) = 79.0
```

### TimeVaryingCovariates Structure

```julia
struct TimeVaryingCovariates
    series::Dict{Symbol, TimeCovariateSeries}
end

# Usage
time_covs = TimeVaryingCovariates(Dict(
    :CRCL => TimeCovariateSeries(StepTimeCovariate(), t_crcl, v_crcl),
    :WT => TimeCovariateSeries(LinearTimeCovariate(), t_wt, v_wt)
))
```

### Applying at Specific Time

```julia
function apply_covariates_at_time(
    params::Dict{Symbol, Float64},
    model::CovariateModel,
    covs::IndividualCovariates,
    t::Float64
)
    # Get covariate values at time t
    current_covs = get_covariates_at_time(covs, t)

    # Apply effects
    return apply_covariates(params, model, current_covs)
end

function get_covariates_at_time(
    covs::IndividualCovariates,
    t::Float64
)
    result = copy(covs.values)

    if covs.time_varying !== nothing
        for (name, series) in covs.time_varying.series
            result[name] = interpolate_at_time(series, t)
        end
    end

    return result
end
```

---

## Common Covariate Effects

### Body Size (Allometric Scaling)

```julia
# Standard allometry based on weight
allometric = CovariateModel([
    CovariateEffect(:CL, :WT, 70.0, PowerCovariate(), 0.75),
    CovariateEffect(:V, :WT, 70.0, PowerCovariate(), 1.0),
    CovariateEffect(:Q, :WT, 70.0, PowerCovariate(), 0.75),
    CovariateEffect(:V2, :WT, 70.0, PowerCovariate(), 1.0)
])

# Alternative using lean body mass
lbm_model = CovariateModel([
    CovariateEffect(:CL, :LBM, 55.0, PowerCovariate(), 0.75),
    CovariateEffect(:V, :LBM, 55.0, PowerCovariate(), 1.0)
])

# BSA-based dosing
bsa_model = CovariateModel([
    CovariateEffect(:CL, :BSA, 1.73, PowerCovariate(), 1.0)
])
```

### Renal Function

```julia
# Creatinine clearance effect
renal_model = CovariateModel([
    # Renal clearance
    CovariateEffect(:CLR, :CRCL, 100.0, LinearCovariate(), 0.007),
    # Or as fraction of total CL
    CovariateEffect(:CL, :CRCL, 100.0, PowerCovariate(), 0.5)
])

# eGFR-based
egfr_model = CovariateModel([
    CovariateEffect(:CL, :EGFR, 90.0, LinearCovariate(), 0.005)
])
```

### Age Effects

```julia
# Linear age effect
age_linear = CovariateEffect(:CL, :AGE, 40.0, LinearCovariate(), -0.006)

# Pediatric maturation
maturation = CovariateEffect(:CL, :PMA, 40.0, PowerCovariate(), 0.75)

# Elderly decline (exponential)
elderly = CovariateEffect(:CL, :AGE, 50.0, ExpCovariate(), -0.01)
```

### Sex/Gender

```julia
# Categorical: use indicator variable
# SEX = 0 for male (reference), SEX = 1 for female
sex_effect = CovariateEffect(:CL, :SEX, 0.0, LinearCovariate(), -0.15)

# Female CL = CL_pop * (1 + (-0.15) * (1 - 0)) = 0.85 * CL_pop
```

### Genetic Polymorphisms

```julia
# CYP2D6 metabolizer status
# PM=0, IM=0.5, EM=1.0 (reference), UM=1.5
cyp2d6_effect = CovariateEffect(:CLM, :CYP2D6, 1.0, PowerCovariate(), 1.0)

# Or as categories with separate effects
# Implement using indicator variables for each category
```

---

## Population Simulation with Covariates

### Complete Example

```julia
using NeoPKPD
using Random

# 1. Define model and typical parameters
model = TwoCompOral()
typical = TwoCompOralParams(
    Ka = 1.5,
    CL = 10.0,
    V1 = 50.0,
    Q = 5.0,
    V2 = 100.0
)

# 2. Define covariate model
cov_model = CovariateModel("full_model", [
    # Allometric scaling
    CovariateEffect(:CL, :WT, 70.0, PowerCovariate(), 0.75),
    CovariateEffect(:V1, :WT, 70.0, PowerCovariate(), 1.0),
    CovariateEffect(:Q, :WT, 70.0, PowerCovariate(), 0.75),
    CovariateEffect(:V2, :WT, 70.0, PowerCovariate(), 1.0),

    # Age effect on CL
    CovariateEffect(:CL, :AGE, 45.0, LinearCovariate(), -0.005),

    # Renal function
    CovariateEffect(:CL, :CRCL, 100.0, LinearCovariate(), 0.003)
])

# 3. Generate population covariates
Random.seed!(42)
n_subjects = 100

covariates = Vector{IndividualCovariates}()
for i in 1:n_subjects
    wt = 70.0 + randn() * 15.0
    age = 45.0 + randn() * 12.0
    # CRCL depends on age
    crcl = max(30.0, 120.0 - age * 0.8 + randn() * 15.0)

    push!(covariates, IndividualCovariates(
        Dict(:WT => wt, :AGE => age, :CRCL => crcl),
        nothing
    ))
end

# 4. IIV (residual after covariate adjustment)
omega = OmegaMatrix([
    0.04 0.0;    # CL: 20% CV (reduced from 30% without covariates)
    0.04 0.0;    # V1: 20% CV
    0.0  0.16    # Ka: 40% CV
])

# 5. Create population specification
doses = [DoseEvent(0.0, 500.0)]
base_spec = ModelSpec(model, "cov_sim", typical, doses)

pop_spec = PopulationSpec(
    base_spec,
    n = n_subjects,
    omega = omega,
    covariate_model = cov_model,
    covariates = covariates,
    seed = 12345
)

# 6. Simulate
grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate_population(pop_spec, grid, solver)
```

### Analyzing Covariate Effects

```julia
using Statistics

# Extract realized parameters and covariates
cl_values = [p[:CL] for p in result.params]
wt_values = [c.values[:WT] for c in covariates]
age_values = [c.values[:AGE] for c in covariates]

# Correlation with covariates
corr_cl_wt = cor(cl_values, wt_values)
corr_cl_age = cor(cl_values, age_values)

println("Correlation CL vs WT: $(round(corr_cl_wt, digits=3))")
println("Correlation CL vs AGE: $(round(corr_cl_age, digits=3))")

# CL by weight quartiles
wt_q = quantile(wt_values, [0.25, 0.5, 0.75])
low_wt = cl_values[wt_values .<= wt_q[1]]
high_wt = cl_values[wt_values .>= wt_q[3]]

println("\nCL by weight quartile:")
println("  Low WT (≤$(round(wt_q[1], digits=1)) kg): $(round(mean(low_wt), digits=2)) L/hr")
println("  High WT (≥$(round(wt_q[3], digits=1)) kg): $(round(mean(high_wt), digits=2)) L/hr")
```

---

## Covariate Model Building

### Stepwise Selection

```julia
# Forward selection approach
function stepwise_covariate_selection(
    base_model::PopulationModel,
    candidate_effects::Vector{CovariateEffect},
    data::PopulationData;
    alpha_forward = 0.05,
    alpha_backward = 0.01
)
    selected = CovariateEffect[]
    remaining = copy(candidate_effects)

    # Forward selection
    improved = true
    while improved && !isempty(remaining)
        improved = false
        best_effect = nothing
        best_delta_ofv = 0.0

        for effect in remaining
            # Fit model with this effect
            test_model = add_covariate(base_model, effect)
            fit = fit_population(data, test_model)

            delta_ofv = base_model.ofv - fit.ofv
            p_value = 1 - cdf(Chisq(1), delta_ofv)

            if p_value < alpha_forward && delta_ofv > best_delta_ofv
                best_delta_ofv = delta_ofv
                best_effect = effect
                improved = true
            end
        end

        if improved
            push!(selected, best_effect)
            filter!(e -> e != best_effect, remaining)
            base_model = add_covariate(base_model, best_effect)
        end
    end

    return CovariateModel(selected)
end
```

### Covariate Significance

```julia
# Test significance of covariate effect
function test_covariate_significance(
    full_model::PopulationModel,
    reduced_model::PopulationModel,  # Without the covariate
    n_params_diff::Int = 1
)
    delta_ofv = reduced_model.ofv - full_model.ofv
    p_value = 1 - cdf(Chisq(n_params_diff), delta_ofv)

    return (
        delta_ofv = delta_ofv,
        p_value = p_value,
        significant = p_value < 0.05
    )
end
```

---

## Complete Example

```julia
using NeoPKPD
using Statistics
using Random

# ============================================
# Population PK with Comprehensive Covariates
# ============================================

println("=== Covariate Modeling Example ===\n")

# 1. Model setup
model = TwoCompOral()
typical = TwoCompOralParams(Ka=1.5, CL=10.0, V1=50.0, Q=5.0, V2=100.0)

# 2. Full covariate model
println("--- Covariate Model ---")
cov_model = CovariateModel("comprehensive", [
    # Allometry
    CovariateEffect(:CL, :WT, 70.0, PowerCovariate(), 0.75),
    CovariateEffect(:V1, :WT, 70.0, PowerCovariate(), 1.0),
    CovariateEffect(:Q, :WT, 70.0, PowerCovariate(), 0.75),
    CovariateEffect(:V2, :WT, 70.0, PowerCovariate(), 1.0),
    # Demographics
    CovariateEffect(:CL, :AGE, 45.0, LinearCovariate(), -0.006),
    CovariateEffect(:CL, :SEX, 0.0, LinearCovariate(), -0.12),
    # Renal function
    CovariateEffect(:CL, :CRCL, 100.0, LinearCovariate(), 0.004)
])

for eff in cov_model.effects
    println("  $(eff.param) ~ $(eff.covariate): $(eff.kind) (β=$(eff.beta), ref=$(eff.ref))")
end

# 3. Generate realistic population
Random.seed!(42)
n = 200

covariates = Vector{IndividualCovariates}()
for i in 1:n
    sex = rand() < 0.5 ? 0.0 : 1.0
    age = 20.0 + 40.0 * rand()
    wt = sex == 0 ? 75.0 + randn() * 12.0 : 65.0 + randn() * 10.0
    crcl = max(30.0, 130.0 - age * 0.9 + randn() * 15.0)

    push!(covariates, IndividualCovariates(
        Dict(:WT => wt, :AGE => age, :SEX => sex, :CRCL => crcl),
        nothing
    ))
end

# 4. Summarize covariates
println("\n--- Population Covariates ---")
wts = [c.values[:WT] for c in covariates]
ages = [c.values[:AGE] for c in covariates]
crcls = [c.values[:CRCL] for c in covariates]
sexes = [c.values[:SEX] for c in covariates]

println("  WT: $(round(mean(wts), digits=1)) ± $(round(std(wts), digits=1)) kg")
println("  AGE: $(round(mean(ages), digits=1)) ± $(round(std(ages), digits=1)) years")
println("  CRCL: $(round(mean(crcls), digits=1)) ± $(round(std(crcls), digits=1)) mL/min")
println("  Female: $(round(mean(sexes) * 100, digits=0))%")

# 5. IIV (residual variability after covariates)
omega = OmegaMatrix([
    0.0225 0.0 0.0;   # CL: 15% residual CV
    0.0225 0.0 0.0;   # V1: 15% residual CV
    0.0 0.0 0.09      # Ka: 30% CV
])

# 6. Simulate
doses = [DoseEvent(0.0, 500.0)]
base_spec = ModelSpec(model, "cov_example", typical, doses)

pop_spec = PopulationSpec(
    base_spec,
    n = n,
    omega = omega,
    covariate_model = cov_model,
    covariates = covariates,
    seed = 12345
)

grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

println("\n--- Simulation ---")
result = simulate_population(pop_spec, grid, solver)
println("Simulated $(length(result.individuals)) subjects")

# 7. Analyze covariate impact
println("\n--- Covariate Impact on CL ---")
cl_values = [p[:CL] for p in result.params]

# By sex
male_cl = cl_values[sexes .== 0]
female_cl = cl_values[sexes .== 1]
println("  Male CL: $(round(mean(male_cl), digits=2)) L/hr")
println("  Female CL: $(round(mean(female_cl), digits=2)) L/hr")
println("  Ratio (F/M): $(round(mean(female_cl)/mean(male_cl), digits=2))")

# By age tertiles
age_t = quantile(ages, [0.33, 0.67])
young_cl = cl_values[ages .< age_t[1]]
old_cl = cl_values[ages .> age_t[2]]
println("\n  Young (<$(round(age_t[1], digits=0)) y): $(round(mean(young_cl), digits=2)) L/hr")
println("  Old (>$(round(age_t[2], digits=0)) y): $(round(mean(old_cl), digits=2)) L/hr")

# By renal function
normal_crcl = cl_values[crcls .>= 90]
impaired_crcl = cl_values[crcls .< 60]
println("\n  Normal CRCL (≥90): $(round(mean(normal_crcl), digits=2)) L/hr")
println("  Impaired CRCL (<60): $(round(mean(impaired_crcl), digits=2)) L/hr")

# 8. PK outcomes
println("\n--- PK Outcomes ---")
cmax = [maximum(ind.observations[:conc]) for ind in result.individuals]
println("Cmax: $(round(mean(cmax), digits=2)) ± $(round(std(cmax), digits=2)) mg/L")
println("CV(Cmax): $(round(std(cmax)/mean(cmax)*100, digits=1))%")
```

---

## See Also

- [IIV](iiv.md) - Inter-individual variability
- [IOV](iov.md) - Inter-occasion variability
- [Residual Error](residual-error.md) - Observation error
- [Parameter Estimation](../estimation/index.md) - Estimating covariate effects
