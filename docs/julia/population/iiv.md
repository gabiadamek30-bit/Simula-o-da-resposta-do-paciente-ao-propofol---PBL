# Inter-Individual Variability (IIV)

Comprehensive guide for modeling between-subject variability in population PK/PD models.

---

## Overview

Inter-individual variability (IIV) captures the differences in pharmacokinetic and pharmacodynamic parameters between subjects in a population.

```julia
using NeoPKPD

# Define IIV with 30% CV on CL and 20% CV on V
omega = OmegaMatrix([
    0.09 0.0;   # ω²_CL (30% CV)
    0.0  0.04   # ω²_V  (20% CV)
])
```

---

## Mathematical Foundation

### Log-Normal Distribution

Individual parameters are derived from population values using log-normal random effects:

$$P_i = \theta \cdot e^{\eta_i}$$

Where:
- $P_i$ = Individual parameter value for subject $i$
- $\theta$ = Typical (population) parameter value
- $\eta_i$ = Random effect, $\eta_i \sim N(0, \omega^2)$

### Omega to CV Relationship

For log-normal distributions:

$$CV = \sqrt{e^{\omega^2} - 1}$$

For small ω (< 0.5), this approximates to:

$$CV \approx \omega$$

| ω² | ω | Exact CV | Approx CV |
|----|---|----------|-----------|
| 0.01 | 0.10 | 10.0% | 10% |
| 0.04 | 0.20 | 20.1% | 20% |
| 0.09 | 0.30 | 30.9% | 30% |
| 0.16 | 0.40 | 43.3% | 40% |
| 0.25 | 0.50 | 58.4% | 50% |
| 0.36 | 0.60 | 77.3% | 60% |

---

## OmegaMatrix Structure

### Class Definition

```julia
struct OmegaMatrix
    param_names::Vector{Symbol}      # Parameter names [:CL, :V, ...]
    matrix::Matrix{Float64}          # Full variance-covariance matrix
    cholesky_L::LowerTriangular      # Lower Cholesky factor for sampling
end
```

### Creating Omega Matrices

```julia
# Diagonal (uncorrelated parameters)
omega = OmegaMatrix([
    0.09 0.0;    # CL
    0.0  0.04    # V
])

# Full covariance (correlated parameters)
omega = OmegaMatrix([
    0.09 0.03;   # CL, correlation with V
    0.03 0.04    # V
])

# With explicit parameter names
omega = OmegaMatrix(
    [:CL, :V, :Ka],
    [
        0.09 0.02 0.0;
        0.02 0.04 0.0;
        0.0  0.0  0.16
    ]
)
```

### Omega Matrix Properties

```julia
# Check if correlations exist
has_corr = has_correlations(omega)

# Get diagonal elements only
diag_omegas = get_diagonal_omegas(omega)
# Returns: Dict(:CL => 0.09, :V => 0.04)

# Get correlation matrix
corr = get_correlation_matrix(omega)
# Returns correlation coefficients (-1 to 1)

# Ensure positive definite
omega_pd = ensure_positive_definite_omega(omega)
```

---

## IIVSpec Structure

### Class Definition

```julia
struct IIVSpec{K<:RandomEffectKind}
    kind::K                          # LogNormalIIV()
    omegas::Dict{Symbol, Float64}    # Diagonal omegas (backward compatible)
    omega_matrix::Union{OmegaMatrix, Nothing}  # Full covariance (optional)
    seed::Int                        # Random seed for reproducibility
    n::Int                           # Number of subjects
end
```

### Creating IIV Specifications

```julia
# Simple diagonal IIV
iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.09, :V => 0.04),
    seed = 12345,
    n = 100
)

# Full covariance IIV
omega = OmegaMatrix([
    0.09 0.03;
    0.03 0.04
])

iiv = IIVSpec(
    LogNormalIIV(),
    omega_matrix = omega,
    seed = 12345,
    n = 100
)
```

---

## Sampling Random Effects

### Diagonal Sampling

```julia
# Sample uncorrelated etas
function sample_etas_diagonal(
    omegas::Dict{Symbol, Float64},
    n::Int;
    seed::Int = nothing
)
    if seed !== nothing
        Random.seed!(seed)
    end

    etas = Dict{Symbol, Vector{Float64}}()
    for (param, omega_sq) in omegas
        omega = sqrt(omega_sq)
        etas[param] = randn(n) .* omega
    end

    return etas
end

# Usage
etas = sample_etas_diagonal(Dict(:CL => 0.09, :V => 0.04), 100, seed=42)
# etas[:CL] is Vector{Float64} of length 100
# etas[:V] is Vector{Float64} of length 100
```

### Correlated Sampling (Cholesky)

```julia
# Sample correlated etas using Cholesky decomposition
function sample_etas_correlated(
    omega::OmegaMatrix,
    n::Int;
    seed::Int = nothing
)
    if seed !== nothing
        Random.seed!(seed)
    end

    p = size(omega.matrix, 1)
    Z = randn(n, p)                    # Standard normal samples
    etas = Z * omega.cholesky_L'       # Transform to correlated

    # Convert to Dict
    result = Dict{Symbol, Vector{Float64}}()
    for (i, param) in enumerate(omega.param_names)
        result[param] = etas[:, i]
    end

    return result
end

# Usage
omega = OmegaMatrix([0.09 0.03; 0.03 0.04])
etas = sample_etas_correlated(omega, 100, seed=42)
```

---

## Applying IIV

### Log-Normal Application

```julia
# Apply etas to typical parameters
function apply_etas_lognormal(
    typical_params::Dict{Symbol, Float64},
    etas::Dict{Symbol, Vector{Float64}},
    subject_index::Int
)
    individual_params = Dict{Symbol, Float64}()

    for (param, theta) in typical_params
        if haskey(etas, param)
            eta = etas[param][subject_index]
            individual_params[param] = theta * exp(eta)
        else
            individual_params[param] = theta
        end
    end

    return individual_params
end

# Example
typical = Dict(:CL => 10.0, :V => 50.0, :Ka => 1.5)
etas = Dict(:CL => [-0.2, 0.1, 0.3], :V => [0.1, -0.1, 0.0])

for i in 1:3
    ind = apply_etas_lognormal(typical, etas, i)
    println("Subject $i: CL=$(round(ind[:CL], digits=2)), V=$(round(ind[:V], digits=2))")
end
# Subject 1: CL=8.19, V=55.26
# Subject 2: CL=11.05, V=45.24
# Subject 3: CL=13.50, V=50.00
```

---

## Population Simulation with IIV

### Basic Simulation

```julia
using NeoPKPD

# 1. Define typical parameters
typical_params = TwoCompOralParams(
    Ka = 1.5,
    CL = 10.0,
    V1 = 50.0,
    Q = 5.0,
    V2 = 100.0
)

# 2. Define IIV
omega = OmegaMatrix([
    0.09 0.02 0.0;    # CL (correlated with V1)
    0.02 0.04 0.0;    # V1
    0.0  0.0  0.16    # Ka (uncorrelated)
])

# 3. Create base specification
doses = [DoseEvent(0.0, 100.0)]
base_spec = ModelSpec(
    TwoCompOral(),
    "pop_sim",
    typical_params,
    doses
)

# 4. Create population specification
pop_spec = PopulationSpec(
    base_spec,
    n = 100,
    omega = omega,
    seed = 12345
)

# 5. Simulate
grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate_population(pop_spec, grid, solver)
```

### Accessing Results

```julia
# Individual results
for (i, ind) in enumerate(result.individuals[1:5])
    println("Subject $i:")
    println("  Cmax = $(maximum(ind.observations[:conc]))")
    println("  Parameters: $(result.params[i])")
end

# Population summaries
summary = result.summaries[:conc]
println("\nPopulation Summary:")
println("Mean Cmax: $(maximum(summary.mean))")
println("Median Cmax: $(maximum(summary.median))")
println("5th percentile: $(maximum(summary.quantiles[0.05]))")
println("95th percentile: $(maximum(summary.quantiles[0.95]))")

# Eta values
println("\nEta distribution:")
println("η_CL: mean=$(mean(result.etas[:CL])), sd=$(std(result.etas[:CL]))")
println("η_V1: mean=$(mean(result.etas[:V1])), sd=$(std(result.etas[:V1]))")
```

---

## Correlation Between Parameters

### Interpreting Correlations

```julia
# Full omega matrix with correlation
omega = OmegaMatrix([
    0.09 0.04;    # ω²_CL = 0.09, cov(CL,V) = 0.04
    0.04 0.16     # ω²_V = 0.16
])

# Calculate correlation coefficient
omega_cl = sqrt(0.09)   # 0.30
omega_v = sqrt(0.16)    # 0.40
covariance = 0.04
correlation = covariance / (omega_cl * omega_v)  # 0.04 / (0.30 * 0.40) = 0.33

println("Correlation between CL and V: $(round(correlation, digits=2))")
# Subjects with higher CL tend to have higher V
```

### Positive vs Negative Correlations

```julia
# Positive correlation: CL and V increase together
omega_pos = OmegaMatrix([
    0.09  0.04;
    0.04  0.09
])

# Negative correlation: as CL increases, V decreases
omega_neg = OmegaMatrix([
    0.09 -0.04;
   -0.04  0.09
])

# Common physiological correlations:
# - CL and V often positively correlated (larger subjects have both)
# - Ka and F may be negatively correlated (fast absorption, lower F)
```

---

## IIV on Different Parameter Types

### Volume Parameters

```julia
# Volume typically scales with body size
# Common CV: 20-40%
omega_v = 0.04 to 0.16  # ω² values
```

### Clearance Parameters

```julia
# Clearance also scales with body size and organ function
# Common CV: 20-50%
omega_cl = 0.04 to 0.25
```

### Absorption Parameters (Ka)

```julia
# Absorption rate often highly variable
# Common CV: 30-60%
omega_ka = 0.09 to 0.36
```

### Bioavailability (F)

```julia
# Bioavailability constrained 0-1, use logit transform
# Or use proportional model with constraints
# Common CV: 20-40%
```

---

## Shrinkage

### Eta Shrinkage

Shrinkage indicates how much individual estimates are pulled toward population values:

```julia
function calculate_eta_shrinkage(etas::Vector{Float64}, omega::Float64)
    var_eta = var(etas)
    shrinkage = 1 - sqrt(var_eta) / sqrt(omega)
    return shrinkage * 100  # As percentage
end

# Interpretation:
# < 20%: Good, individual estimates reliable
# 20-40%: Moderate, some uncertainty
# > 40%: High, individual estimates unreliable
```

### Causes of High Shrinkage

- Sparse sampling (few observations per subject)
- Low IIV relative to residual error
- Parameters not well estimated from data

---

## Complete Example

```julia
using NeoPKPD
using Statistics

# ============================================
# Population PK Simulation with IIV
# ============================================

println("=== Population PK with IIV ===\n")

# 1. Define model and typical parameters
model = TwoCompOral()
typical = TwoCompOralParams(
    Ka = 1.5,     # /hr
    CL = 10.0,    # L/hr
    V1 = 50.0,    # L
    Q = 5.0,      # L/hr
    V2 = 100.0    # L
)

# 2. Define IIV with correlations
println("--- IIV Specification ---")
omega = OmegaMatrix(
    [:CL, :V1, :Ka],
    [
        0.09 0.03 0.00;   # CL: 30% CV, correlated with V1
        0.03 0.04 0.00;   # V1: 20% CV
        0.00 0.00 0.16    # Ka: 40% CV, uncorrelated
    ]
)

println("Omega matrix:")
display(omega.matrix)

corr = get_correlation_matrix(omega)
println("\nCorrelation matrix:")
display(corr)

# 3. Setup simulation
doses = [DoseEvent(0.0, 500.0)]  # 500 mg single dose
base_spec = ModelSpec(model, "iiv_demo", typical, doses)

pop_spec = PopulationSpec(
    base_spec,
    n = 200,
    omega = omega,
    seed = 42
)

grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# 4. Simulate
println("\n--- Running Simulation ---")
result = simulate_population(pop_spec, grid, solver)
println("Simulated $(length(result.individuals)) subjects")

# 5. Parameter distribution
println("\n--- Realized Parameters ---")
cl_values = [p[:CL] for p in result.params]
v1_values = [p[:V1] for p in result.params]
ka_values = [p[:Ka] for p in result.params]

println("CL: $(round(mean(cl_values), digits=2)) ± $(round(std(cl_values), digits=2)) L/hr")
println("V1: $(round(mean(v1_values), digits=2)) ± $(round(std(v1_values), digits=2)) L")
println("Ka: $(round(mean(ka_values), digits=2)) ± $(round(std(ka_values), digits=2)) /hr")

# 6. Check correlation in realized parameters
corr_cl_v1 = cor(cl_values, v1_values)
corr_cl_ka = cor(cl_values, ka_values)
println("\nParameter correlations:")
println("ρ(CL, V1) = $(round(corr_cl_v1, digits=3))")
println("ρ(CL, Ka) = $(round(corr_cl_ka, digits=3))")

# 7. PK metrics
println("\n--- PK Metrics ---")
cmax_values = [maximum(ind.observations[:conc]) for ind in result.individuals]
tmax_values = [ind.times[argmax(ind.observations[:conc])] for ind in result.individuals]

# Approximate AUC using trapezoidal rule
function approx_auc(times, conc)
    auc = 0.0
    for i in 2:length(times)
        auc += 0.5 * (conc[i] + conc[i-1]) * (times[i] - times[i-1])
    end
    return auc
end

auc_values = [approx_auc(ind.times, ind.observations[:conc]) for ind in result.individuals]

println("Cmax: $(round(mean(cmax_values), digits=2)) ± $(round(std(cmax_values), digits=2)) mg/L")
println("Tmax: $(round(mean(tmax_values), digits=2)) ± $(round(std(tmax_values), digits=2)) hr")
println("AUC:  $(round(mean(auc_values), digits=1)) ± $(round(std(auc_values), digits=1)) mg*hr/L")

# 8. Percentiles
println("\n--- Population Percentiles (Cmax) ---")
percentiles = [5, 25, 50, 75, 95]
for p in percentiles
    val = quantile(cmax_values, p/100)
    println("  $(p)th percentile: $(round(val, digits=2)) mg/L")
end

# 9. CV calculation
cv_cmax = std(cmax_values) / mean(cmax_values) * 100
cv_auc = std(auc_values) / mean(auc_values) * 100
println("\n--- Variability ---")
println("CV(Cmax): $(round(cv_cmax, digits=1))%")
println("CV(AUC): $(round(cv_auc, digits=1))%")
```

---

## See Also

- [IOV](iov.md) - Inter-occasion variability
- [Covariates](covariates.md) - Covariate effects on parameters
- [Residual Error](residual-error.md) - Observation error models
- [Parameter Estimation](../estimation/index.md) - Fitting population models
