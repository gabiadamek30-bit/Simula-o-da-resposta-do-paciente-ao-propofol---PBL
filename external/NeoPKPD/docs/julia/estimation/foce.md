# FOCE-I Method

First-Order Conditional Estimation with Interaction (FOCE-I) is the industry-standard estimation method for nonlinear mixed-effects modeling, matching NONMEM's implementation.

---

## Overview

FOCE-I linearizes the model around individual eta estimates rather than zero, providing more accurate parameter estimates for most PK/PD models.

### Key Features

- **Full Laplacian correction**: log|H_eta| computed from actual Hessian
- **Proper CWRES**: Sensitivity gradients dF/d_eta included
- **Analytic solutions**: AD-compatible for common PK models
- **Sandwich estimator**: Robust standard errors
- **Covariate effects on IIV**: Support for eta-covariate relationships
- **IOV support**: Inter-Occasion Variability estimation

---

## Usage

### Basic FOCE-I Estimation

```julia
using NeoPKPD

# Create observed data
observed = ObservedData(
    subjects = [
        Subject(
            id = "1",
            times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            observations = Dict(:conc => [8.2, 9.1, 7.8, 5.2, 2.8, 1.5, 0.4]),
            doses = [DoseEvent(0.0, 500.0)],
            covariates = Dict(:WT => 70.0, :AGE => 45.0)
        ),
        # ... more subjects
    ]
)

# Model specification
model = OneCompIVBolus()
model_spec = ModelSpec(model, "pk", nothing, nothing)

# Configure FOCE-I estimation
config = EstimationConfig(
    method = FOCEIMethod(),
    theta_init = [10.0, 50.0],           # Initial [CL, V]
    omega_init = diagm([0.09, 0.04]),    # Initial omega (30% CV CL, 20% CV V)
    sigma_init = ResidualErrorSpec(
        ProportionalError(),
        ProportionalErrorParams(0.1),     # 10% proportional error
        :conc
    ),
    max_iter = 1000,
    tol = 1e-6,
    compute_se = true,
    compute_ci = true,
    ci_level = 0.95,
    verbose = true
)

# Create grid and solver
grid = SimGrid(0.0, 24.0, 0:0.1:24)
solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^6)

# Run estimation
result = foce_estimate(observed, model_spec, config, grid, solver, StableRNG(12345))

# Access results
println("Fixed effects (theta):")
println("  CL = $(result.theta[1]) ± $(result.theta_se[1])")
println("  V  = $(result.theta[2]) ± $(result.theta_se[2])")
println("\nRandom effects (omega):")
println("  omega_CL = $(result.omega[1,1]) (CV = $(sqrt(result.omega[1,1])*100)%)")
println("  omega_V  = $(result.omega[2,2]) (CV = $(sqrt(result.omega[2,2])*100)%)")
println("\nResidual error:")
println("  sigma = $(result.sigma)")
println("\nModel fit:")
println("  OFV = $(result.ofv)")
println("  AIC = $(result.aic)")
println("  BIC = $(result.bic)")
```

---

## Configuration Options

### EstimationConfig Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `method` | FOCEIMethod | Required | FOCE-I method specification |
| `theta_init` | Vector{Float64} | Required | Initial fixed effect estimates |
| `omega_init` | Matrix{Float64} | Required | Initial omega matrix |
| `sigma_init` | ResidualErrorSpec | Required | Initial residual error specification |
| `max_iter` | Int | 1000 | Maximum iterations |
| `tol` | Float64 | 1e-6 | Convergence tolerance |
| `compute_se` | Bool | true | Compute standard errors |
| `compute_ci` | Bool | true | Compute confidence intervals |
| `ci_level` | Float64 | 0.95 | Confidence level |
| `omega_structure` | Symbol | :diagonal | Omega structure (:diagonal, :block, :full) |
| `verbose` | Bool | false | Print progress |
| `seed` | UInt64 | 12345 | Random seed for reproducibility |

### FOCEIMethod Options

```julia
FOCEIMethod(;
    centered = false,        # Use FOCE (not FOCE-I) when true
    compute_cwres = true,    # Compute conditional weighted residuals
    laplacian = true,        # Include Laplacian correction
    eta_tol = 1e-8,          # Inner optimization tolerance
    max_eta_iter = 100       # Max iterations for eta optimization
)
```

---

## Omega Structures

### Diagonal Omega

Independent random effects (most common):

```julia
config = EstimationConfig(
    # ...
    omega_init = diagm([0.09, 0.04]),  # Diagonal
    omega_structure = :diagonal
)
```

### Block Omega

Correlated random effects within blocks:

```julia
# CL and V correlated, ka independent
omega_init = [
    0.09  0.02  0.0;   # CL
    0.02  0.04  0.0;   # V
    0.0   0.0   0.16   # ka
]

config = EstimationConfig(
    omega_init = omega_init,
    omega_structure = :block
)
```

### Full Omega

All random effects correlated:

```julia
omega_init = [
    0.09  0.02  0.01;
    0.02  0.04  0.005;
    0.01  0.005 0.16
]

config = EstimationConfig(
    omega_init = omega_init,
    omega_structure = :full
)
```

---

## Residual Error Models

### Additive Error

$$\epsilon \sim N(0, \sigma^2)$$

```julia
sigma_init = ResidualErrorSpec(
    AdditiveError(),
    AdditiveErrorParams(0.5),  # sigma = 0.5
    :conc
)
```

### Proportional Error

$$\epsilon \sim N(0, (\sigma \cdot F)^2)$$

```julia
sigma_init = ResidualErrorSpec(
    ProportionalError(),
    ProportionalErrorParams(0.1),  # 10% CV
    :conc
)
```

### Combined Error

$$\epsilon \sim N(0, (\sigma_{add}^2 + (\sigma_{prop} \cdot F)^2))$$

```julia
sigma_init = ResidualErrorSpec(
    CombinedError(),
    CombinedErrorParams(0.5, 0.1),  # additive=0.5, proportional=10%
    :conc
)
```

### Exponential Error

$$\log(Y) = \log(F) + \epsilon, \quad \epsilon \sim N(0, \sigma^2)$$

```julia
sigma_init = ResidualErrorSpec(
    ExponentialError(),
    ExponentialErrorParams(0.15),  # 15% CV on log-scale
    :conc
)
```

---

## BLQ Handling

FOCE-I supports Below Limit of Quantification (BLQ) data with multiple methods:

### M1 Method (Discard)

Simply discards BLQ observations:

```julia
blq_config = BLQConfig(
    method = :M1,
    lloq = 0.1
)

config = EstimationConfig(
    # ... other options
    blq_config = blq_config
)
```

### M2 Method (Imputation)

Imputes BLQ values at LLOQ/2 or 0:

```julia
blq_config = BLQConfig(
    method = :M2,
    lloq = 0.1,
    impute_value = :half  # or :zero
)
```

### M3 Method (Censored Likelihood)

Gold standard - uses censored likelihood:

```julia
blq_config = BLQConfig(
    method = :M3,
    lloq = 0.1
)
```

M3 properly accounts for the information that observations are below LLOQ using:

$$L_{BLQ} = \Phi\left(\frac{LLOQ - F}{\sigma}\right)$$

---

## Covariate Effects on IIV

Model covariate effects on inter-individual variability:

```julia
# Define covariate effects
covariate_effects = [
    CovariateOnIIV(
        eta_name = :eta_CL,
        covariate_name = :WT,
        effect_type = :exponential,
        reference_value = 70.0
    ),
    CovariateOnIIV(
        eta_name = :eta_V,
        covariate_name = :AGE,
        effect_type = :linear,
        reference_value = 40.0
    )
]

config = EstimationConfig(
    # ... other options
    covariate_effects = covariate_effects
)
```

This implements:
- **Exponential**: $\omega_{adjusted} = \omega_{base} \cdot \exp(\theta_{cov} \cdot (COV - COV_{ref}))$
- **Linear**: $\omega_{adjusted} = \omega_{base} \cdot (1 + \theta_{cov} \cdot (COV - COV_{ref}))$

---

## Inter-Occasion Variability (IOV)

Support multiple occasions with occasion-specific random effects:

```julia
# Define IOV specification
iov_specs = [
    EstimationIOVSpec(
        eta_name = :eta_CL,
        occasion_names = [:OCC1, :OCC2, :OCC3],
        omega_iov = 0.04  # 20% IOV on CL
    )
]

config = EstimationConfig(
    # ... other options
    iov_specs = iov_specs
)
```

Total eta: $\eta_{total} = \eta_{IIV} + \eta_{IOV}[occasion]$

---

## Diagnostics

### FOCEIDiagnostics Structure

```julia
struct FOCEIDiagnostics
    eta_hessians::Vector{Matrix{Float64}}      # Per-subject Hessians
    laplacian_corrections::Vector{Float64}      # Laplacian terms
    likelihood_contributions::Vector{Float64}   # Per-subject -2LL
    prior_contributions::Vector{Float64}        # Prior penalty
    interaction_enabled::Bool                   # FOCE-I vs FOCE
    cwres_fallback_subjects::Vector{String}     # CWRES issues
end
```

### Accessing Diagnostics

```julia
result = foce_estimate(observed, model_spec, config, grid, solver, rng)

# Individual estimates
for (i, ind) in enumerate(result.individual_estimates)
    println("Subject $(ind.subject_id):")
    println("  eta = $(ind.eta)")
    println("  ipred = $(ind.ipred[1:3])...")
    println("  cwres = $(ind.cwres[1:3])...")
end

# Shrinkage
println("Eta shrinkage: $(result.eta_shrinkage)")
println("Epsilon shrinkage: $(result.epsilon_shrinkage)")
```

---

## Standard Errors

### Asymptotic Standard Errors

Computed from inverse Fisher Information Matrix:

```julia
config = EstimationConfig(
    # ...
    compute_se = true,
    se_method = :hessian  # Default
)
```

### Sandwich (Robust) Standard Errors

More robust to model misspecification:

```julia
config = EstimationConfig(
    # ...
    compute_se = true,
    se_method = :sandwich
)
```

The sandwich estimator uses: $SE = R^{-1} S R^{-1}$ where R is the Hessian and S is the score covariance.

---

## Convergence Monitoring

### Objective Function Value (OFV)

$$OFV = -2 \cdot \log(L) = \sum_i \left[ \ln|C_i| + (y_i - f_i)^T C_i^{-1} (y_i - f_i) + \ln|H_{\eta,i}| \right]$$

### Convergence Criteria

The algorithm converges when:
1. Relative parameter change < tolerance
2. Gradient norm < gradient tolerance
3. OFV change < function tolerance

```julia
result.converged      # Did it converge?
result.n_iterations   # How many iterations?
result.final_gradient # Final gradient vector
```

---

## Mathematical Details

### FOCE-I Objective Function

For subject $i$:

$$OFV_i = \ln|C_i(\hat{\eta}_i)| + (y_i - f_i(\hat{\eta}_i))^T C_i^{-1}(\hat{\eta}_i) (y_i - f_i(\hat{\eta}_i)) + \hat{\eta}_i^T \Omega^{-1} \hat{\eta}_i + \ln|H_{\eta,i}|$$

Where:
- $f_i(\eta)$ = Model predictions for subject i
- $C_i(\eta)$ = Residual variance matrix
- $\hat{\eta}_i$ = Individual eta estimates (modes)
- $H_{\eta,i}$ = Hessian of -2LL with respect to eta

### FOCE vs FOCE-I

| Method | Interaction Term | $\partial f / \partial \eta$ |
|--------|-----------------|------------------------------|
| FOCE | No | Evaluated at $\eta = 0$ |
| FOCE-I | Yes | Evaluated at $\eta = \hat{\eta}$ |

FOCE-I includes the "interaction" between random effects and residual error.

---

## Example: Full Population PK Analysis

```julia
using NeoPKPD
using StableRNGs

# Load data
data = load_estimation_data("pk_data.csv")

observed = ObservedData(
    subjects = [
        Subject(
            id = row.id,
            times = row.times,
            observations = Dict(:conc => row.dv),
            doses = row.doses,
            covariates = Dict(:WT => row.wt, :SEX => row.sex)
        )
        for row in data
    ]
)

# Two-compartment model
model = TwoCompIVBolus()
model_spec = ModelSpec(model, "pk", nothing, nothing)

# Initial estimates from NCA
config = EstimationConfig(
    method = FOCEIMethod(centered=false, compute_cwres=true),
    theta_init = [10.0, 30.0, 3.0, 100.0],  # CL, V1, Q, V2
    omega_init = diagm([0.09, 0.04, 0.16, 0.04]),
    sigma_init = ResidualErrorSpec(
        CombinedError(),
        CombinedErrorParams(0.1, 0.05),
        :conc
    ),
    omega_structure = :diagonal,
    blq_config = BLQConfig(method=:M3, lloq=0.01),
    max_iter = 2000,
    tol = 1e-8,
    compute_se = true,
    verbose = true
)

grid = SimGrid(0.0, 72.0, 0:0.1:72)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

# Run estimation
result = foce_estimate(observed, model_spec, config, grid, solver, StableRNG(42))

# Print results
println("=== FOCE-I Estimation Results ===")
println("\nFixed Effects:")
params = ["CL", "V1", "Q", "V2"]
for (i, p) in enumerate(params)
    println("  $p = $(round(result.theta[i], sigdigits=4)) (RSE = $(round(result.theta_rse[i], digits=1))%)")
end

println("\nRandom Effects (CV%):")
for (i, p) in enumerate(params)
    cv = sqrt(result.omega[i,i]) * 100
    println("  $p IIV = $(round(cv, digits=1))%")
end

println("\nEta Shrinkage:")
for (i, p) in enumerate(params)
    println("  $p = $(round(result.eta_shrinkage[i]*100, digits=1))%")
end

println("\nModel Fit:")
println("  OFV = $(round(result.ofv, digits=2))")
println("  AIC = $(round(result.aic, digits=2))")
println("  BIC = $(round(result.bic, digits=2))")
println("  Converged: $(result.converged)")
```

---

## See Also

- [SAEM Algorithm](saem.md) - Alternative for complex models
- [Laplacian Method](laplacian.md) - For sparse data
- [Diagnostics](diagnostics.md) - Model validation
- [Bootstrap](bootstrap.md) - Non-parametric uncertainty

