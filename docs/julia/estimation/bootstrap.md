# Bootstrap Analysis

Industry-standard bootstrap methods for parameter uncertainty estimation, supporting regulatory submissions (FDA/EMA).

---

## Overview

Bootstrap analysis provides non-parametric estimates of parameter uncertainty by resampling the data and re-estimating parameters multiple times.

### Key Features

- **Case bootstrap**: Standard FDA/EMA-recommended approach
- **Parametric bootstrap**: Simulate from fitted model
- **Residual bootstrap**: Preserve covariate structure
- **Stratified resampling**: By study, dose, formulation
- **Multiple CI methods**: Percentile, BCa, Basic
- **Parallel execution**: Multi-threaded for large studies
- **Regulatory output**: Formatted tables for submissions

---

## Bootstrap Types

### Case Bootstrap (Non-Parametric)

Resample subjects with replacement - the gold standard for regulatory submissions:

```julia
using NeoPKPD
using StableRNGs

# After initial estimation
result = foce_estimate(observed, model_spec, config, grid, solver, rng)

# Configure bootstrap
bootstrap_spec = BootstrapSpec(
    n_bootstrap = 1000,           # FDA recommends ≥500
    bootstrap_type = CaseBootstrap(),
    seed = 12345,
    parallel = true,
    ci_level = 0.95,
    ci_method = PercentileCI()
)

# Run bootstrap
bootstrap_result = run_bootstrap(
    observed,
    model_spec,
    config,
    grid,
    solver,
    result,  # Original estimates as starting point
    bootstrap_spec
)

# Access results
println("Bootstrap SE: ", bootstrap_result.theta_se)
println("Bootstrap 95% CI:")
for i in 1:length(result.theta)
    println("  θ$i: [$(bootstrap_result.theta_ci_lower[i]), $(bootstrap_result.theta_ci_upper[i])]")
end
```

### Parametric Bootstrap

Simulate new data from the fitted model:

```julia
bootstrap_spec = BootstrapSpec(
    n_bootstrap = 500,
    bootstrap_type = ParametricBootstrap(
        n_simulations_per_subject = 1
    ),
    seed = 12345
)
```

Parametric bootstrap is useful when:
- Small sample size
- Want to assess model adequacy
- Case bootstrap has high failure rate

### Residual Bootstrap

Resample residuals and add to predictions:

```julia
bootstrap_spec = BootstrapSpec(
    n_bootstrap = 500,
    bootstrap_type = ResidualBootstrap(
        standardize = true  # Standardize residuals before resampling
    ),
    seed = 12345
)
```

Residual bootstrap preserves:
- Original covariate structure
- Observation times
- Dosing patterns

---

## Stratified Resampling

Maintain proportions within strata when resampling:

```julia
# Stratify by study and formulation
bootstrap_spec = BootstrapSpec(
    n_bootstrap = 1000,
    bootstrap_type = CaseBootstrap(),
    stratify_by = [:study, :formulation],  # Stratification variables
    seed = 12345
)
```

Stratification is important for:
- Pooled analyses across studies
- Bioequivalence studies
- Different dose groups

---

## Confidence Interval Methods

### Percentile CI (Default)

Simple percentile method:

$$CI = [\theta^*_{\alpha/2}, \theta^*_{1-\alpha/2}]$$

```julia
bootstrap_spec = BootstrapSpec(
    ci_method = PercentileCI(),
    ci_level = 0.95
)
```

### BCa (Bias-Corrected and Accelerated)

More accurate for skewed distributions:

```julia
bootstrap_spec = BootstrapSpec(
    ci_method = BCCI(acceleration = 0.0),  # Computed from jackknife
    ci_level = 0.95
)
```

BCa adjusts for:
- Bias in bootstrap distribution
- Skewness (acceleration)

### Basic Bootstrap CI

Also known as reverse percentile:

$$CI = [2\hat{\theta} - \theta^*_{1-\alpha/2}, 2\hat{\theta} - \theta^*_{\alpha/2}]$$

```julia
bootstrap_spec = BootstrapSpec(
    ci_method = BasicCI(),
    ci_level = 0.95
)
```

---

## BootstrapSpec Configuration

### All Parameters

```julia
BootstrapSpec(;
    # Number of replicates
    n_bootstrap = 1000,              # FDA recommends ≥500

    # Bootstrap type
    bootstrap_type = CaseBootstrap(), # CaseBootstrap, ParametricBootstrap, ResidualBootstrap

    # Stratification
    stratify_by = Symbol[],          # Variables to stratify by

    # Reproducibility
    seed = 12345,

    # Execution
    parallel = true,                 # Use multiple threads

    # Confidence intervals
    ci_level = 0.95,                 # 95% CI
    ci_method = PercentileCI(),      # PercentileCI, BCCI, BasicCI

    # What to compute
    compute_omega_ci = true,         # CI for random effects
    compute_sigma_ci = true,         # CI for residual error

    # Quality control
    min_success_rate = 0.8           # Minimum successful runs (80%)
)
```

---

## Bootstrap Results

### BootstrapResult Structure

```julia
struct BootstrapResult
    # Fixed effects
    theta_estimates::Matrix{Float64}   # All bootstrap estimates (n_boot × n_theta)
    theta_mean::Vector{Float64}        # Mean of bootstrap estimates
    theta_se::Vector{Float64}          # Bootstrap standard error
    theta_rse::Vector{Float64}         # Relative SE (%)
    theta_ci_lower::Vector{Float64}    # Lower CI bound
    theta_ci_upper::Vector{Float64}    # Upper CI bound

    # Bias correction
    original_estimate::Vector{Float64} # Original point estimate
    bias::Vector{Float64}              # Bootstrap bias
    bias_corrected::Vector{Float64}    # Bias-corrected estimate

    # Random effects
    omega_summary::OmegaBootstrapSummary
    sigma_summary::SigmaBootstrapSummary

    # Shrinkage distribution
    eta_shrinkage::Matrix{Float64}     # Shrinkage per replicate

    # Diagnostics
    diagnostics::BootstrapDiagnostics

    # CI info
    ci_level::Float64
    ci_method::String
end
```

### Accessing Results

```julia
result = run_bootstrap(...)

# Fixed effects
println("Theta estimates: ", result.theta_mean)
println("Standard errors: ", result.theta_se)
println("RSE%: ", result.theta_rse)

# CIs
for i in 1:length(result.theta_mean)
    println("θ$i: $(result.theta_mean[i]) [$(result.theta_ci_lower[i]), $(result.theta_ci_upper[i])]")
end

# Bias
println("Bias: ", result.bias)
println("Bias-corrected: ", result.bias_corrected)

# Random effects
println("Omega SE: ", result.omega_summary.se)
println("Sigma SE: ", result.sigma_summary.se)
```

---

## Bootstrap Diagnostics

### BootstrapDiagnostics Structure

```julia
struct BootstrapDiagnostics
    n_successful::Int               # Successful runs
    n_failed::Int                   # Failed runs
    convergence_rate::Float64       # Success rate (FDA requires ≥80%)
    median_iterations::Float64      # Median iterations to converge
    outlier_indices::Vector{Int}    # Potential outlier estimates
    failure_reasons::Dict{Symbol, Int}  # Why runs failed
    rse_stability::Vector{Float64}  # RSE of SE estimate
end
```

### Checking Bootstrap Quality

```julia
diag = result.diagnostics

println("Successful runs: $(diag.n_successful)/$(diag.n_successful + diag.n_failed)")
println("Success rate: $(diag.convergence_rate * 100)%")

# FDA requires ≥80% success rate
if diag.convergence_rate < 0.8
    println("WARNING: Success rate below FDA threshold")
end

# Check for outliers
if !isempty(diag.outlier_indices)
    println("Outlier estimates at indices: ", diag.outlier_indices)
end

# Failure reasons
for (reason, count) in diag.failure_reasons
    println("  $reason: $count failures")
end
```

---

## Parallel Execution

### Using Thread Pool

```julia
using Base.Threads
println("Available threads: ", nthreads())

bootstrap_spec = BootstrapSpec(
    n_bootstrap = 1000,
    parallel = true,  # Uses all available threads
    seed = 12345
)
```

### Advanced Parallel Config

```julia
parallel_config = ParallelConfig(
    ThreadedBackend(8),      # Use 8 threads
    seed = 12345,
    load_balance = true,     # Dynamic load balancing
    progress = true          # Show progress
)

bootstrap_spec = BootstrapSpec(
    n_bootstrap = 1000,
    parallel = parallel_config
)
```

---

## Regulatory Output

### Generate Formatted Table

```julia
# Generate regulatory-compliant summary table
summary_table = generate_bootstrap_summary(result)

println(summary_table)
# Output:
# Parameter  Estimate    SE      RSE%    95% CI
# ─────────────────────────────────────────────────
# CL         10.5       0.82    7.8%    [8.9, 12.1]
# V          52.3       3.41    6.5%    [45.6, 59.2]
# ...
```

### Export for Regulatory Submission

```julia
# Format as regulatory table
reg_table = format_regulatory_table(result;
    parameter_names = ["CL (L/h)", "V (L)", "ka (1/h)"],
    ci_level = 0.95,
    decimal_places = 3
)

# Save to file
write("bootstrap_table.txt", reg_table)
```

---

## Influential Subject Analysis

### Jackknife Influence

```julia
# Identify influential subjects
influence_analysis = jackknife_influence(
    result,
    observed,
    model_spec,
    config,
    grid,
    solver
)

# Subjects that strongly affect estimates
influential = influence_analysis.influential_subjects
println("Most influential subjects: ", influential)
```

### Bootstrap Coverage Check

```julia
# Validate CI coverage via simulation
coverage = compute_bootstrap_coverage(
    true_theta = [10.0, 50.0],  # True values (if known)
    bootstrap_results = [result1, result2, ...]  # Multiple runs
)

println("Empirical coverage: $(coverage * 100)%")
# Should be close to 95% for 95% CI
```

---

## Example: Full Bootstrap Analysis

```julia
using NeoPKPD
using StableRNGs

# Load and prepare data
observed = load_observed_data("pk_study.csv")

# Model specification
model = TwoCompOral()
model_spec = ModelSpec(model, "pk", nothing, nothing)

# Initial estimation
foce_config = EstimationConfig(
    method = FOCEIMethod(),
    theta_init = [10.0, 30.0, 5.0, 100.0, 1.5],
    omega_init = diagm([0.09, 0.04, 0.16, 0.04, 0.25]),
    sigma_init = ResidualErrorSpec(CombinedError(), CombinedErrorParams(0.1, 0.1), :conc),
    compute_se = true
)

grid = SimGrid(0.0, 72.0, 0:0.5:72)
solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^7)

# Run initial estimation
base_result = foce_estimate(observed, model_spec, foce_config, grid, solver, StableRNG(42))

println("=== Initial FOCE-I Estimates ===")
println("Theta: ", base_result.theta)
println("Asymptotic SE: ", base_result.theta_se)

# Configure bootstrap
bootstrap_spec = BootstrapSpec(
    n_bootstrap = 1000,
    bootstrap_type = CaseBootstrap(),
    stratify_by = [:study],  # Stratify by study
    ci_level = 0.95,
    ci_method = PercentileCI(),
    parallel = true,
    min_success_rate = 0.8,
    seed = 12345
)

# Run bootstrap (this takes time!)
println("\n=== Running Bootstrap (1000 replicates) ===")
boot_result = run_bootstrap(
    observed, model_spec, foce_config, grid, solver,
    base_result, bootstrap_spec
)

# Report results
println("\n=== Bootstrap Results ===")
println("Success rate: $(boot_result.diagnostics.convergence_rate * 100)%")

param_names = ["CL", "V1", "Q", "V2", "ka"]
println("\nParameter Estimates with Bootstrap SE and 95% CI:")
println("-" ^ 60)
println("Parameter   Estimate    SE       RSE%     95% CI")
println("-" ^ 60)

for (i, name) in enumerate(param_names)
    est = round(base_result.theta[i], sigdigits=4)
    se = round(boot_result.theta_se[i], sigdigits=3)
    rse = round(boot_result.theta_rse[i], digits=1)
    ci_lo = round(boot_result.theta_ci_lower[i], sigdigits=4)
    ci_hi = round(boot_result.theta_ci_upper[i], sigdigits=4)
    println("$name         $est      $se      $rse%    [$ci_lo, $ci_hi]")
end

# Compare asymptotic vs bootstrap SE
println("\n=== Asymptotic vs Bootstrap SE ===")
for (i, name) in enumerate(param_names)
    asymp = round(base_result.theta_se[i], sigdigits=3)
    boot = round(boot_result.theta_se[i], sigdigits=3)
    ratio = round(boot / asymp, digits=2)
    println("$name: Asymptotic=$asymp, Bootstrap=$boot, Ratio=$ratio")
end

# Check for outliers
if !isempty(boot_result.diagnostics.outlier_indices)
    println("\nWARNING: Potential outlier estimates detected")
end

# Generate regulatory table
reg_table = format_regulatory_table(boot_result,
    parameter_names = ["CL (L/h)", "V1 (L)", "Q (L/h)", "V2 (L)", "ka (1/h)"]
)
println("\n=== Regulatory Table ===")
println(reg_table)
```

---

## Best Practices

### Sample Size

- FDA recommends n ≥ 500 bootstrap replicates
- 1000 replicates is standard for regulatory submissions
- For complex models, 2000+ may improve stability

### Success Rate

- Target ≥ 80% successful runs
- If < 80%, consider:
  - Simplifying the model
  - Improving initial estimates
  - Using SAEM instead of FOCE

### Stratification

- Always stratify for pooled analyses
- Match stratification to study design
- Check strata have sufficient subjects

### CI Method Selection

| Method | When to Use |
|--------|-------------|
| Percentile | Default, symmetric distributions |
| BCa | Skewed parameters (e.g., variance components) |
| Basic | Alternative to percentile |

---

## See Also

- [FOCE-I Method](foce.md) - Primary estimation
- [SAEM Algorithm](saem.md) - Alternative estimation
- [Diagnostics](diagnostics.md) - Model validation
- [Model Comparison](comparison.md) - Comparing models

