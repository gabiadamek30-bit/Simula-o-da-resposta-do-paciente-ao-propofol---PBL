# Inter-Occasion Variability (IOV)

Comprehensive guide for modeling within-subject variability across different dosing occasions.

---

## Overview

Inter-occasion variability (IOV) captures the random fluctuations in an individual's parameters from one dosing occasion to another, representing "day-to-day" or "visit-to-visit" variability.

```julia
using NeoPKPD

# Define IOV with 15% CV on CL and 10% CV on Ka per occasion
iov = IOVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.0225, :Ka => 0.01),  # π² values
    occasion_def = OccasionDefinition(:dose_times)
)
```

---

## Mathematical Foundation

### Nested Random Effects

IOV is nested within IIV, meaning each subject has both:
- A subject-specific deviation from the population (η)
- Occasion-specific deviations from their individual value (κ)

$$P_{ij} = \theta \cdot e^{\eta_i + \kappa_{ij}}$$

Where:
- $P_{ij}$ = Parameter for subject $i$ on occasion $j$
- $\theta$ = Population typical value
- $\eta_i$ = IIV random effect, $\eta_i \sim N(0, \omega^2)$
- $\kappa_{ij}$ = IOV random effect, $\kappa_{ij} \sim N(0, \pi^2)$

### Variance Decomposition

Total variance in the log-domain:

$$Var(\log P_{ij}) = \omega^2 + \pi^2$$

Between-subject variance is $\omega^2$, within-subject variance is $\pi^2$.

---

## IOVSpec Structure

### Class Definition

```julia
struct IOVSpec{K<:RandomEffectKind}
    kind::K                              # LogNormalIIV()
    pis::Dict{Symbol, Float64}           # π² values for each parameter
    seed::Int                            # Random seed
    occasion_def::OccasionDefinition     # How to define occasions
end

struct OccasionDefinition
    mode::Symbol                         # :dose_times, :custom, :fixed_duration
    dose_compartment::Int                # CMT for dose-based occasions
    duration::Float64                    # For fixed duration mode
    custom_times::Vector{Float64}        # For custom mode
end
```

### Creating IOV Specifications

```julia
# IOV based on dose times
iov = IOVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.0225, :Ka => 0.01),
    seed = 12345,
    occasion_def = OccasionDefinition(:dose_times)
)

# IOV with fixed occasion duration (e.g., weekly)
iov = IOVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.0225),
    seed = 12345,
    occasion_def = OccasionDefinition(:fixed_duration, duration=168.0)  # 7 days
)

# IOV with custom occasion boundaries
iov = IOVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.0225),
    seed = 12345,
    occasion_def = OccasionDefinition(:custom, custom_times=[0.0, 24.0, 48.0, 72.0])
)
```

---

## Occasion Determination

### Dose-Based Occasions

Each dose defines a new occasion:

```julia
# Multiple dose regimen
doses = [
    DoseEvent(0.0, 100.0),    # Occasion 1
    DoseEvent(24.0, 100.0),   # Occasion 2
    DoseEvent(48.0, 100.0),   # Occasion 3
    DoseEvent(72.0, 100.0)    # Occasion 4
]

# Occasion boundaries at dose times
occasion_def = OccasionDefinition(:dose_times)
```

### Fixed Duration Occasions

```julia
# Weekly dosing with variable timing
# Still want occasions to align with weeks
occasion_def = OccasionDefinition(:fixed_duration, duration=168.0)

# For a simulation 0-336 hours:
# Occasion 1: 0-168 hours
# Occasion 2: 168-336 hours
```

### Finding Occasion at Time

```julia
function occasion_index_at_time(
    t::Float64,
    occasion_def::OccasionDefinition,
    doses::Vector{DoseEvent}
)
    if occasion_def.mode == :dose_times
        # Find most recent dose
        occ = 1
        for (i, dose) in enumerate(doses)
            if dose.time <= t
                occ = i
            end
        end
        return occ
    elseif occasion_def.mode == :fixed_duration
        return floor(Int, t / occasion_def.duration) + 1
    elseif occasion_def.mode == :custom
        for (i, boundary) in enumerate(occasion_def.custom_times)
            if t < boundary
                return i
            end
        end
        return length(occasion_def.custom_times)
    end
end
```

---

## Sampling IOV Kappas

### Basic Sampling

```julia
function sample_iov_kappas(
    pis::Dict{Symbol, Float64},
    n_subjects::Int,
    n_occasions::Int;
    seed::Int = nothing
)
    if seed !== nothing
        Random.seed!(seed)
    end

    kappas = Dict{Symbol, Matrix{Float64}}()
    for (param, pi_sq) in pis
        pi = sqrt(pi_sq)
        # Matrix: n_subjects × n_occasions
        kappas[param] = randn(n_subjects, n_occasions) .* pi
    end

    return kappas
end

# Usage
kappas = sample_iov_kappas(
    Dict(:CL => 0.0225, :Ka => 0.01),
    n_subjects = 50,
    n_occasions = 4,
    seed = 42
)

# kappas[:CL] is 50 × 4 matrix
# kappas[:CL][i, j] is κ for subject i, occasion j
```

### Applying IOV

```julia
function apply_iov(
    base_param::Float64,
    eta::Float64,          # Subject's IIV
    kappa::Float64         # Occasion-specific IOV
)
    return base_param * exp(eta + kappa)
end

# Full application
function get_individual_occasion_param(
    typical::Float64,
    eta::Float64,
    kappa::Float64
)
    return typical * exp(eta + kappa)
end

# Example
typical_cl = 10.0
eta_cl = 0.2       # Subject has 20% higher CL than typical
kappa_cl = -0.1    # This occasion 10% lower than subject's average

cl_this_occasion = typical_cl * exp(eta_cl + kappa_cl)
# = 10.0 * exp(0.2 - 0.1) = 10.0 * exp(0.1) = 11.05 L/hr
```

---

## Population Simulation with IOV

### Complete Example

```julia
using NeoPKPD

# 1. Define model and parameters
model = TwoCompOral()
typical = TwoCompOralParams(
    Ka = 1.5,
    CL = 10.0,
    V1 = 50.0,
    Q = 5.0,
    V2 = 100.0
)

# 2. Define IIV
omega = OmegaMatrix([
    0.09 0.0 0.0;    # CL: 30% CV
    0.0  0.04 0.0;   # V1: 20% CV
    0.0  0.0 0.16    # Ka: 40% CV
])

iiv = IIVSpec(
    LogNormalIIV(),
    omega_matrix = omega,
    seed = 12345,
    n = 50
)

# 3. Define IOV
iov = IOVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.0225, :Ka => 0.01),  # 15% on CL, 10% on Ka
    seed = 54321,
    occasion_def = OccasionDefinition(:dose_times)
)

# 4. Multiple dose regimen (4 occasions)
doses = [
    DoseEvent(0.0, 100.0),
    DoseEvent(24.0, 100.0),
    DoseEvent(48.0, 100.0),
    DoseEvent(72.0, 100.0)
]

# 5. Create specifications
base_spec = ModelSpec(model, "iov_demo", typical, doses)

pop_spec = PopulationSpec(
    base_spec,
    iiv = iiv,
    iov = iov
)

# 6. Simulate
grid = SimGrid(0.0, 96.0, collect(0.0:0.5:96.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate_population(pop_spec, grid, solver)
```

### Accessing IOV Results

```julia
# Parameters vary by occasion for each subject
for (i, subject) in enumerate(result.individuals[1:3])
    println("Subject $i:")
    for occ in 1:4
        # Get parameters at a time point in this occasion
        t = occ * 24.0 - 12.0  # Middle of each occasion
        params = result.occasion_params[i][occ]
        println("  Occasion $occ: CL=$(round(params[:CL], digits=2)), Ka=$(round(params[:Ka], digits=2))")
    end
end
```

---

## IOV vs IIV Comparison

### Typical Magnitude Relationships

| Variability | Typical CV | Represents |
|-------------|------------|------------|
| IIV (ω) | 20-50% | Between-subject differences |
| IOV (π) | 10-25% | Within-subject fluctuation |
| Residual (σ) | 10-30% | Measurement/model error |

### Rule of Thumb

IOV is typically smaller than IIV (often 30-50% of the IIV magnitude):

```julia
# If IIV CV is 30%
omega_cl = 0.09  # ω² for ~30% CV

# IOV might be 15% (half of IIV)
pi_cl = 0.0225   # π² for ~15% CV
```

---

## When to Include IOV

### Strong Candidates for IOV

1. **Absorption parameters (Ka, F)**
   - Food effects
   - Formulation dissolution variability
   - GI motility changes

2. **Clearance (CL)**
   - Enzyme activity fluctuations
   - Time-of-day effects
   - Disease state changes

3. **First-pass effect (FG, FH)**
   - Hepatic blood flow changes
   - Enzyme induction/inhibition

### Less Common

1. **Volumes (V, V2)**
   - Usually more stable
   - Body composition doesn't change quickly

2. **Distribution (Q)**
   - Relatively stable physiological parameter

---

## Model Diagnostics with IOV

### Detecting Need for IOV

```julia
# Signs that IOV may be needed:
# 1. High residual variability despite good IIV fit
# 2. Individual plots show period-to-period variation
# 3. CWRES show patterns within subjects over time

# Compare models
model_no_iov = fit_population(data, spec_no_iov)
model_with_iov = fit_population(data, spec_with_iov)

# Likelihood ratio test
delta_ofv = model_no_iov.ofv - model_with_iov.ofv
n_additional_params = 2  # Added π² parameters
p_value = 1 - cdf(Chisq(n_additional_params), delta_ofv)

if p_value < 0.05
    println("IOV significantly improves fit (p = $(round(p_value, digits=4)))")
end
```

### Kappa Shrinkage

```julia
function calculate_kappa_shrinkage(
    kappas::Matrix{Float64},  # n_subjects × n_occasions
    pi_sq::Float64
)
    # Pool all kappas
    all_kappas = vec(kappas)
    var_kappa = var(all_kappas)
    shrinkage = 1 - sqrt(var_kappa) / sqrt(pi_sq)
    return shrinkage * 100
end

# High kappa shrinkage indicates:
# - Few observations per occasion
# - Low IOV relative to residual error
# - Difficulty distinguishing IOV from noise
```

---

## Complete Example

```julia
using NeoPKPD
using Statistics

# ============================================
# Multiple Dose Simulation with IIV and IOV
# ============================================

println("=== Population PK with IIV + IOV ===\n")

# 1. Model and parameters
model = OneCompOral()
typical = OneCompOralParams(
    Ka = 1.2,
    CL = 8.0,
    V = 60.0
)

# 2. IIV specification (between-subject)
println("--- Inter-Individual Variability ---")
omega = OmegaMatrix([
    0.09 0.0;    # CL: 30% CV
    0.16 0.0;    # Ka: 40% CV
    0.0  0.04    # V:  20% CV
])

iiv = IIVSpec(LogNormalIIV(), omega_matrix=omega, seed=111, n=40)

for param in [:CL, :Ka, :V]
    cv = sqrt(exp(omega.matrix[findfirst(==(param), omega.param_names), findfirst(==(param), omega.param_names)]) - 1) * 100
    println("  $(param): $(round(cv, digits=1))% CV")
end

# 3. IOV specification (within-subject)
println("\n--- Inter-Occasion Variability ---")
pis = Dict(:CL => 0.0225, :Ka => 0.01)  # 15% on CL, 10% on Ka
iov = IOVSpec(
    LogNormalIIV(),
    pis,
    seed = 222,
    occasion_def = OccasionDefinition(:dose_times)
)

for (param, pi) in pis
    cv = sqrt(exp(pi) - 1) * 100
    println("  $(param): $(round(cv, digits=1))% CV")
end

# 4. Multiple dose regimen (7 days QD)
n_occasions = 7
doses = [DoseEvent(i * 24.0, 500.0) for i in 0:n_occasions-1]

println("\n--- Dosing Regimen ---")
println("  $(n_occasions) doses of 500 mg QD")

# 5. Simulation
base_spec = ModelSpec(model, "iov_example", typical, doses)
pop_spec = PopulationSpec(base_spec, iiv=iiv, iov=iov)

grid = SimGrid(0.0, 168.0, collect(0.0:0.5:168.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

println("\n--- Running Simulation ---")
result = simulate_population(pop_spec, grid, solver)
println("Simulated $(length(result.individuals)) subjects over $(n_occasions) occasions")

# 6. Analyze occasion-to-occasion variability for one subject
println("\n--- Example Subject (ID=1) ---")
subject = 1
println("Base CL (with IIV): $(round(result.params[subject][:CL], digits=2)) L/hr")

for occ in 1:n_occasions
    occ_params = result.occasion_params[subject][occ]
    println("  Occasion $occ: CL=$(round(occ_params[:CL], digits=2)), Ka=$(round(occ_params[:Ka], digits=2))")
end

# 7. Calculate observed variability
println("\n--- Observed Variability ---")

# Between-subject: variability in subject means
subject_mean_cmax = Float64[]
for ind in result.individuals
    # Get Cmax from last occasion (steady-state)
    t_start = (n_occasions - 1) * 24.0
    t_end = n_occasions * 24.0
    indices = findall(t -> t_start <= t <= t_end, ind.times)
    push!(subject_mean_cmax, maximum(ind.observations[:conc][indices]))
end
cv_between = std(subject_mean_cmax) / mean(subject_mean_cmax) * 100
println("Between-subject CV(Cmax at SS): $(round(cv_between, digits=1))%")

# Within-subject: variability across occasions for each subject
cv_within_values = Float64[]
for ind in result.individuals
    occasion_cmax = Float64[]
    for occ in 1:n_occasions
        t_start = (occ - 1) * 24.0
        t_end = occ * 24.0
        indices = findall(t -> t_start <= t < t_end, ind.times)
        push!(occasion_cmax, maximum(ind.observations[:conc][indices]))
    end
    cv_occ = std(occasion_cmax) / mean(occasion_cmax) * 100
    push!(cv_within_values, cv_occ)
end
cv_within = mean(cv_within_values)
println("Within-subject CV(Cmax): $(round(cv_within, digits=1))%")

# 8. Population summary at steady state
println("\n--- Steady-State Summary (Last Occasion) ---")
ss_cmax = Float64[]
for ind in result.individuals
    t_start = (n_occasions - 1) * 24.0
    indices = findall(t -> t >= t_start, ind.times)
    push!(ss_cmax, maximum(ind.observations[:conc][indices]))
end

println("Cmax,ss: $(round(mean(ss_cmax), digits=2)) ± $(round(std(ss_cmax), digits=2)) mg/L")
println("5th percentile: $(round(quantile(ss_cmax, 0.05), digits=2)) mg/L")
println("95th percentile: $(round(quantile(ss_cmax, 0.95), digits=2)) mg/L")
```

---

## See Also

- [IIV](iiv.md) - Inter-individual variability
- [Covariates](covariates.md) - Covariate effects
- [Residual Error](residual-error.md) - Observation error
- [Multiple Dose NCA](../nca/multiple-dose.md) - Steady-state analysis
