# Laplacian Method

The Laplacian approximation method provides fast, efficient parameter estimation particularly suited for sparse data situations.

---

## Overview

The Laplacian method approximates the marginal likelihood by:
1. Finding the mode (MAP estimate) of eta for each subject
2. Using the Laplace approximation for the integral over random effects

### Key Features

- **Fast computation**: Simpler than FOCE-I, no interaction terms
- **Sparse data friendly**: Works well with few observations per subject
- **Analytical gradients**: Efficient optimization
- **BLQ support**: All M1/M2/M3 methods supported
- **Standard errors**: From Hessian approximation

---

## When to Use Laplacian

The Laplacian method is preferred when:

| Situation | Why Laplacian Helps |
|-----------|---------------------|
| Sparse sampling | Few observations per subject (1-3) |
| Large populations | Fast computation scales well |
| Initial estimates | Quick screening before FOCE/SAEM |
| Simple models | One-compartment, few parameters |
| Exploratory analysis | Rapid model comparison |

---

## Usage

### Basic Laplacian Estimation

```julia
using NeoPKPD
using StableRNGs

# Observed data (sparse sampling)
observed = ObservedData(
    subjects = [
        Subject(
            id = "1",
            times = [1.0, 4.0],  # Only 2 observations
            observations = Dict(:conc => [8.5, 3.2]),
            doses = [DoseEvent(0.0, 500.0)]
        ),
        Subject(
            id = "2",
            times = [2.0, 8.0],
            observations = Dict(:conc => [6.8, 1.4]),
            doses = [DoseEvent(0.0, 500.0)]
        ),
        # ... more sparse subjects
    ]
)

# Model specification
model = OneCompIVBolus()
model_spec = ModelSpec(model, "pk", nothing, nothing)

# Configure Laplacian estimation
config = EstimationConfig(
    method = LaplacianMethod(
        max_inner_iter = 50,    # Max iterations for eta optimization
        inner_tol = 1e-6        # Inner optimization tolerance
    ),
    theta_init = [10.0, 50.0],        # CL, V
    omega_init = diagm([0.09, 0.04]),
    sigma_init = ResidualErrorSpec(
        ProportionalError(),
        ProportionalErrorParams(0.15),
        :conc
    ),
    max_iter = 500,
    tol = 1e-6,
    compute_se = true,
    verbose = true
)

grid = SimGrid(0.0, 24.0, 0:0.1:24)
solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^6)

result = laplacian_estimate(observed, model_spec, config, grid, solver, StableRNG(12345))

println("Theta: ", result.theta)
println("Theta SE: ", result.theta_se)
println("OFV: ", result.ofv)
```

---

## LaplacianMethod Configuration

### All Parameters

```julia
LaplacianMethod(;
    max_inner_iter = 50,     # Max iterations for finding eta mode
    inner_tol = 1e-6,        # Convergence tolerance for inner loop
    use_prior = true,        # Include prior on eta in optimization
    hessian_method = :exact  # or :finite_diff
)
```

---

## Mathematical Details

### Laplace Approximation

The marginal likelihood integrates out random effects:

$$L(\theta, \Omega, \sigma | y) = \int p(y | \eta, \theta, \sigma) p(\eta | \Omega) d\eta$$

The Laplace approximation expands around the mode $\hat{\eta}$:

$$\log L \approx \log p(y | \hat{\eta}, \theta, \sigma) + \log p(\hat{\eta} | \Omega) - \frac{1}{2}\log|H|$$

Where $H$ is the Hessian of $-\log p(\eta | y)$ at $\hat{\eta}$.

### Objective Function

$$OFV = \sum_i \left[ (y_i - f_i(\hat{\eta}_i))^T C_i^{-1} (y_i - f_i(\hat{\eta}_i)) + \ln|C_i| + \hat{\eta}_i^T \Omega^{-1} \hat{\eta}_i + \ln|\Omega| + \ln|H_i| \right]$$

### Finding Eta Mode

For each subject, minimize:

$$J_i(\eta) = (y_i - f_i(\eta))^T C_i^{-1} (y_i - f_i(\eta)) + \eta^T \Omega^{-1} \eta$$

This is solved using Newton-Raphson or L-BFGS.

---

## Comparison with FOCE-I

| Aspect | Laplacian | FOCE-I |
|--------|-----------|--------|
| Interaction term | No | Yes |
| F linearization | Around $\hat{\eta}$ | Around $\hat{\eta}$ |
| C linearization | No | Yes (at $\hat{\eta}$) |
| Computational cost | Lower | Higher |
| Accuracy | Good for sparse | Better for dense |
| Implementation | Simpler | More complex |

The main difference is that Laplacian doesn't include the "interaction" term in the residual covariance linearization.

---

## Sparse Data Example

```julia
using NeoPKPD
using StableRNGs

# Very sparse: 1-2 observations per subject
# This is common in Phase I oncology trials

subjects = [
    Subject(id="$(i)", times=[rand([1.0, 2.0, 4.0])],
            observations=Dict(:conc => [rand()*10]),
            doses=[DoseEvent(0.0, 100.0)])
    for i in 1:50
]

observed = ObservedData(subjects=subjects)

model = OneCompIVBolus()
model_spec = ModelSpec(model, "pk", nothing, nothing)

# Laplacian handles sparse data well
config = EstimationConfig(
    method = LaplacianMethod(),
    theta_init = [5.0, 30.0],
    omega_init = diagm([0.25, 0.16]),  # Wide initial IIV
    sigma_init = ResidualErrorSpec(
        CombinedError(),
        CombinedErrorParams(0.5, 0.2),
        :conc
    ),
    max_iter = 500,
    compute_se = true,
    verbose = true
)

grid = SimGrid(0.0, 24.0, 0:0.5:24)
solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^6)

result = laplacian_estimate(observed, model_spec, config, grid, solver, StableRNG(123))

println("=== Sparse Data Analysis ===")
println("CL = $(round(result.theta[1], sigdigits=3)) ± $(round(result.theta_se[1], sigdigits=2))")
println("V  = $(round(result.theta[2], sigdigits=3)) ± $(round(result.theta_se[2], sigdigits=2))")
println("OFV = $(round(result.ofv, digits=2))")
```

---

## BLQ Handling

All BLQ methods are supported:

```julia
# M3 (censored likelihood) - recommended
blq_config = BLQConfig(
    method = :M3,
    lloq = 0.1,
    report_blq_summary = true
)

config = EstimationConfig(
    method = LaplacianMethod(),
    # ... other options
    blq_config = blq_config
)

result = laplacian_estimate(observed, model_spec, config, grid, solver, rng)

# Check BLQ summary
println("BLQ observations: $(result.blq_summary.blq_observations)")
println("BLQ percentage: $(result.blq_summary.blq_percentage)%")
```

---

## Standard Errors

### Hessian-Based SEs

```julia
config = EstimationConfig(
    method = LaplacianMethod(hessian_method=:exact),
    compute_se = true,
    # ...
)

result = laplacian_estimate(...)

# Access SEs
println("Theta SE: ", result.theta_se)
println("Theta RSE%: ", result.theta_rse)

# Confidence intervals
println("95% CI for CL: [$(result.theta_ci_lower[1]), $(result.theta_ci_upper[1])]")
```

### SE Computation Details

Standard errors are computed from the inverse Hessian of the marginal log-likelihood:

$$SE(\hat{\theta}) = \sqrt{diag(H^{-1})}$$

Where H is computed using central finite differences or automatic differentiation.

---

## Parallel Execution

```julia
parallel_config = ParallelConfig(
    ThreadedBackend(4),
    seed = 12345,
    load_balance = true
)

result = laplacian_estimate(
    observed, model_spec, config, grid, solver, rng;
    parallel_config = parallel_config
)
```

---

## Model Fit Statistics

```julia
result = laplacian_estimate(...)

# Objective function value
println("OFV = ", result.ofv)

# Information criteria
println("AIC = ", result.aic)
println("BIC = ", result.bic)

# For model comparison
# Lower AIC/BIC indicates better fit (penalized for complexity)
```

---

## Algorithm Summary

1. **Initialize**: Set initial theta, omega, sigma
2. **E-step**: For each subject, find eta mode by minimizing:
   $$J(\eta) = \|y - f(\eta)\|_C^2 + \eta^T\Omega^{-1}\eta$$
3. **M-step**: Update theta, omega, sigma using the eta modes
4. **Convergence check**: If parameter change < tolerance, stop
5. **Iterate** steps 2-4 until convergence

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Eta mode | $\hat{\eta} = \arg\min_\eta [(y-f)^TC^{-1}(y-f) + \eta^T\Omega^{-1}\eta]$ |
| Laplacian correction | $-\frac{1}{2}\log\|H_\eta\|$ |
| OFV contribution | $\|y-f\|_C^2 + \ln\|C\| + \|\eta\|_\Omega^2 + \ln\|\Omega\| + \ln\|H\|$ |
| Hessian | $H = \nabla_\eta^2 J(\eta)$ |

---

## Limitations

- **Dense data**: FOCE-I may be more accurate
- **High IIV**: May underestimate variability
- **Complex models**: SAEM may be more robust
- **Interaction effects**: Not captured (no C linearization)

---

## See Also

- [FOCE-I Method](foce.md) - More accurate for dense data
- [SAEM Algorithm](saem.md) - More robust for complex models
- [Diagnostics](diagnostics.md) - Model validation

