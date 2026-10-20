# SAEM Algorithm

Stochastic Approximation Expectation Maximization (SAEM) is a robust estimation method that uses MCMC sampling to handle complex models and problematic datasets.

---

## Overview

SAEM alternates between:
1. **E-step**: MCMC sampling of random effects from their posterior distribution
2. **M-step**: Stochastic approximation update of population parameters

### Key Features

- **Robust convergence**: Less prone to local minima than FOCE
- **Handles complex models**: High IIV, multimodal likelihoods
- **Multiple MCMC chains**: Improved mixing and convergence diagnostics
- **Adaptive proposals**: Automatic tuning of Metropolis-Hastings proposals
- **Parallel execution**: Multi-threaded E-step for large datasets
- **Convergence diagnostics**: Gelman-Rubin, effective sample size

---

## When to Use SAEM

SAEM is preferred over FOCE-I when:

| Situation | Why SAEM Helps |
|-----------|----------------|
| Large IIV (>50% CV) | Better handles extreme random effects |
| Complex nonlinear models | More robust to model nonlinearity |
| Multimodal likelihood | Multiple chains explore the space |
| FOCE convergence issues | SAEM often converges where FOCE fails |
| Many random effects | Stochastic sampling scales better |

---

## Usage

### Basic SAEM Estimation

```julia
using NeoPKPD
using StableRNGs

# Create observed data
observed = ObservedData(
    subjects = [
        Subject(
            id = "1",
            times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            observations = Dict(:conc => [8.2, 9.1, 7.8, 5.2, 2.8, 1.5, 0.4]),
            doses = [DoseEvent(0.0, 500.0)],
            covariates = Dict(:WT => 70.0)
        ),
        # ... more subjects
    ]
)

# Model specification
model = TwoCompOral()
model_spec = ModelSpec(model, "pk", nothing, nothing)

# Configure SAEM estimation
config = EstimationConfig(
    method = SAEMMethod(
        n_burn = 200,          # Burn-in iterations
        n_iter = 300,          # Main iterations
        n_chains = 3,          # MCMC chains per subject
        n_mcmc_steps = 50,     # MCMC steps per E-step
        adapt_proposal = true, # Adaptive Metropolis-Hastings
        target_acceptance = 0.3
    ),
    theta_init = [15.0, 30.0, 5.0, 100.0, 1.5],  # CL, V1, Q, V2, ka
    omega_init = diagm([0.09, 0.04, 0.16, 0.04, 0.25]),
    sigma_init = ResidualErrorSpec(
        ProportionalError(),
        ProportionalErrorParams(0.15),
        :conc
    ),
    max_iter = 500,
    tol = 1e-6,
    compute_se = true,
    verbose = true,
    seed = 12345
)

# Create grid and solver
grid = SimGrid(0.0, 48.0, 0:0.1:48)
solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^6)

# Run estimation
result = saem_estimate(observed, model_spec, config, grid, solver, StableRNG(12345))

println("Theta: ", result.theta)
println("Omega diagonal: ", diag(result.omega))
println("OFV: ", result.ofv)
println("Converged: ", result.converged)
```

---

## SAEMMethod Configuration

### All Parameters

```julia
SAEMMethod(;
    # Iteration control
    n_burn = 200,              # Burn-in iterations (adaptation phase)
    n_iter = 300,              # Main iterations (averaging phase)

    # MCMC settings
    n_chains = 3,              # Number of parallel chains per subject
    n_mcmc_steps = 50,         # MCMC steps per E-step

    # Proposal adaptation
    adapt_proposal = true,     # Enable adaptive MH
    target_acceptance = 0.3,   # Target acceptance rate (0.2-0.4)
    adaptation_interval = 20,  # Adapt every N iterations

    # Step size schedule
    step_size_schedule = :standard,  # or :fast, :slow

    # Convergence monitoring
    track_diagnostics = true,  # Track MCMC diagnostics
    use_all_chains = true,     # Average across all chains
    parallel_chains = false    # Parallel chain execution (within subject)
)
```

### Step Size Schedules

The step size controls the stochastic approximation:

```julia
# Standard schedule: γ_k = 1 for k ≤ n_burn, then γ_k = 1/(k - n_burn + 1)
SAEMMethod(step_size_schedule = :standard)

# Fast convergence: γ_k = 1/(k^0.6)
SAEMMethod(step_size_schedule = :fast)

# Slow (more stable): γ_k = 1/(k^0.7)
SAEMMethod(step_size_schedule = :slow)
```

---

## MCMC Details

### Metropolis-Hastings Sampling

For each subject, SAEM samples eta from the posterior:

$$p(\eta_i | y_i, \theta, \Omega, \sigma) \propto p(y_i | \eta_i, \theta, \sigma) \cdot p(\eta_i | \Omega)$$

The algorithm uses random walk MH with proposal:

$$\eta^* \sim N(\eta_{current}, \sigma_{proposal}^2)$$

### Adaptive Proposals

During burn-in, proposals are adapted to achieve target acceptance rate:

```julia
# If acceptance rate > target + 0.05: increase proposal SD by 10%
# If acceptance rate < target - 0.05: decrease proposal SD by 10%
```

---

## Sufficient Statistics

SAEM updates parameters using sufficient statistics computed from MCMC samples:

### Fixed Effects Update (theta)

Uses gradient-based update with stochastic approximation:

$$\theta_{k+1} = \theta_k + \gamma_k \cdot \nabla_\theta Q(\theta | \eta_{sampled})$$

### Random Effects Update (omega)

Updates omega using empirical covariance:

$$S_\Omega^{(k+1)} = (1 - \gamma_k) S_\Omega^{(k)} + \gamma_k \cdot \frac{1}{N} \sum_{i=1}^N E[\eta_i \eta_i^T | y_i]$$

The expectation is approximated by averaging over MCMC samples.

### Residual Error Update (sigma)

Similar stochastic approximation for sigma parameters.

---

## Convergence Diagnostics

### SAEMDiagnostics Structure

```julia
struct SAEMDiagnostics
    acceptance_rates::Vector{Vector{Float64}}  # Per-subject, per-iteration
    mean_acceptance_rate::Float64              # Overall mean
    proposal_sds::Vector{Vector{Float64}}      # Final proposals per subject
    theta_trace::Vector{Vector{Float64}}       # Theta history
    omega_trace::Vector{Vector{Float64}}       # Omega history
    ofv_trace::Vector{Float64}                 # OFV history
    gelman_rubin::Vector{Float64}              # R-hat statistics
    effective_sample_size::Vector{Float64}     # ESS per eta
    converged::Bool                            # Convergence status
end
```

### Accessing Diagnostics

```julia
result = saem_estimate(...)

# MCMC diagnostics
diag = result.saem_diagnostics
println("Mean acceptance rate: $(diag.mean_acceptance_rate)")
println("Gelman-Rubin R-hat: $(diag.gelman_rubin)")
println("Effective sample size: $(diag.effective_sample_size)")

# Trace plots
using Plots
plot(diag.theta_trace, label=["CL" "V"], title="Parameter Trace")
plot(diag.ofv_trace, title="OFV Convergence")
```

### Gelman-Rubin Diagnostic

R-hat should be < 1.1 for convergence:

```julia
# Check all R-hat values
all(diag.gelman_rubin .< 1.1)  # true = good convergence
```

### Effective Sample Size

ESS should be > 100 for reliable estimates:

```julia
minimum(diag.effective_sample_size) > 100  # Check minimum ESS
```

---

## BLQ Handling in SAEM

SAEM supports all BLQ methods with proper likelihood integration:

```julia
blq_config = BLQConfig(
    method = :M3,  # Censored likelihood (recommended)
    lloq = 0.05
)

config = EstimationConfig(
    method = SAEMMethod(),
    # ... other options
    blq_config = blq_config
)
```

For M3, the BLQ contribution to the posterior is:

$$p(y_{BLQ} | \eta, \theta, \sigma) = \Phi\left(\frac{LLOQ - f(\eta, \theta)}{\sigma}\right)$$

---

## Parallel Execution

### Parallel E-Step

Run subjects in parallel using multiple threads:

```julia
using Base.Threads

# Set threads before Julia starts: JULIA_NUM_THREADS=8
println("Using $(nthreads()) threads")

# Enable parallel execution
parallel_config = ParallelConfig(
    ThreadedBackend(nthreads()),
    seed = 12345,
    load_balance = true,
    progress = true
)

result = saem_estimate(
    observed, model_spec, config, grid, solver, rng;
    parallel_config = parallel_config
)
```

### Parallel Chains

Run multiple chains per subject in parallel:

```julia
config = EstimationConfig(
    method = SAEMMethod(
        n_chains = 4,
        parallel_chains = true  # Chains run in parallel
    ),
    # ... other options
)
```

---

## Example: Complex Model with High IIV

```julia
using NeoPKPD
using StableRNGs

# Simulate data with high variability
# (In practice, load real data)

# Transit compartment absorption model
model = TransitAbsorption(n_transit=3)
model_spec = ModelSpec(model, "pk", nothing, nothing)

# High IIV scenario
config = EstimationConfig(
    method = SAEMMethod(
        n_burn = 300,          # More burn-in for complex models
        n_iter = 500,          # More iterations
        n_chains = 4,          # More chains for better mixing
        n_mcmc_steps = 100,    # More MCMC steps
        target_acceptance = 0.25,
        track_diagnostics = true
    ),
    theta_init = [5.0, 50.0, 0.5, 3],   # CL, V, ktr, n_transit
    omega_init = diagm([0.36, 0.16, 0.49]),  # 60%, 40%, 70% CV
    sigma_init = ResidualErrorSpec(
        CombinedError(),
        CombinedErrorParams(0.1, 0.15),
        :conc
    ),
    omega_structure = :diagonal,
    compute_se = true,
    verbose = true,
    seed = 42
)

grid = SimGrid(0.0, 72.0, 0:0.5:72)
solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^7)

result = saem_estimate(observed, model_spec, config, grid, solver, StableRNG(42))

# Check convergence
println("=== SAEM Results ===")
println("Converged: $(result.converged)")
println("Final OFV: $(result.ofv)")

println("\nMCMC Diagnostics:")
diag = result.saem_diagnostics
println("  Mean acceptance: $(round(diag.mean_acceptance_rate, digits=2))")
println("  Max R-hat: $(round(maximum(diag.gelman_rubin), digits=3))")
println("  Min ESS: $(round(minimum(diag.effective_sample_size), digits=0))")

println("\nParameter Estimates:")
println("  CL = $(round(result.theta[1], sigdigits=4)) ($(round(sqrt(result.omega[1,1])*100, digits=1))% IIV)")
println("  V  = $(round(result.theta[2], sigdigits=4)) ($(round(sqrt(result.omega[2,2])*100, digits=1))% IIV)")
println("  ktr = $(round(result.theta[3], sigdigits=4)) ($(round(sqrt(result.omega[3,3])*100, digits=1))% IIV)")
```

---

## Standard Errors from SAEM

SAEM doesn't directly produce analytical SEs. Use bootstrap for uncertainty:

```julia
# Run bootstrap after SAEM
bootstrap_spec = BootstrapSpec(
    n_bootstrap = 500,
    parallel = true,
    seed = 12345
)

bootstrap_result = run_bootstrap(
    observed, model_spec, config, grid, solver,
    result,  # Initial from SAEM
    bootstrap_spec
)

println("Bootstrap SE: ", bootstrap_result.theta_se)
println("Bootstrap 95% CI:")
for i in 1:length(result.theta)
    println("  θ$i: [$(bootstrap_result.theta_ci_lower[i]), $(bootstrap_result.theta_ci_upper[i])]")
end
```

Alternatively, use Louis' missing information principle:

```julia
# SE using Louis' method (based on MCMC samples)
se_result = compute_standard_errors_marginal(
    result.theta, result.omega, result.sigma,
    subjects, model_spec, config;
    method = LouisSE(n_importance_samples=200),
    grid = grid,
    solver = solver
)
```

---

## Algorithm Summary

### SAEM Algorithm Steps

1. **Initialize**: Set theta, omega, sigma; initialize eta chains
2. **For k = 1 to K:**
   - **E-step**: For each subject, run MCMC to sample eta from posterior
   - **M-step**: Update theta, omega, sigma using stochastic approximation
   - If k ≤ n_burn: step_size = 1 (exploration)
   - If k > n_burn: step_size = 1/(k - n_burn + 1) (averaging)
3. **Post-processing**: Compute final OFV, individual estimates, diagnostics

### Mathematical Formulation

E-step computes:
$$Q_k(\theta | \theta_k) = E_{\eta | y, \theta_k}[\log L(\theta; y, \eta)]$$

M-step updates:
$$\theta_{k+1} = \theta_k + \gamma_k (\hat{\theta}_k - \theta_k)$$

Where $\hat{\theta}_k$ maximizes the approximate Q function.

---

## Comparison with FOCE-I

| Aspect | FOCE-I | SAEM |
|--------|--------|------|
| Eta estimation | Mode (optimization) | Posterior samples (MCMC) |
| Local minima | Susceptible | More robust |
| Computation | Faster | Slower |
| Standard errors | Analytical | Bootstrap/Louis |
| High IIV | Can struggle | Handles well |
| Sparse data | Fair | Good |

---

## See Also

- [FOCE-I Method](foce.md) - Faster alternative for standard models
- [Laplacian Method](laplacian.md) - For sparse data
- [Bootstrap](bootstrap.md) - Uncertainty quantification for SAEM
- [Diagnostics](diagnostics.md) - Model validation

