# Bootstrap Analysis

Industry-standard bootstrap methods for parameter uncertainty estimation.

---

## Overview

```python
from neopkpd.estimation import run_bootstrap, BootstrapConfig

bootstrap_config = BootstrapConfig(
    n_bootstrap=1000,
    parallel=True,
    seed=12345
)

bootstrap_result = run_bootstrap(
    data, model, config, result, bootstrap_config
)
```

---

## Bootstrap Types

### Case Bootstrap (Default)

Resample subjects with replacement - FDA/EMA recommended:

```python
bootstrap_config = BootstrapConfig(
    n_bootstrap=1000,
    bootstrap_type="case",
    seed=12345
)
```

### Parametric Bootstrap

Simulate new data from fitted model:

```python
bootstrap_config = BootstrapConfig(
    n_bootstrap=500,
    bootstrap_type="parametric"
)
```

### Residual Bootstrap

Resample residuals:

```python
bootstrap_config = BootstrapConfig(
    n_bootstrap=500,
    bootstrap_type="residual",
    standardize_residuals=True
)
```

---

## Configuration

### BootstrapConfig Parameters

```python
from neopkpd.estimation import BootstrapConfig

bootstrap_config = BootstrapConfig(
    # Number of replicates
    n_bootstrap=1000,           # FDA recommends ≥500

    # Bootstrap type
    bootstrap_type="case",      # "case", "parametric", "residual"

    # Stratification
    stratify_by=["study"],      # Stratification variables

    # Confidence intervals
    ci_level=0.95,              # 95% CI
    ci_method="percentile",     # "percentile", "bca", "basic"

    # Execution
    parallel=True,              # Use multiple cores
    n_jobs=-1,                  # -1 = all cores

    # Quality control
    min_success_rate=0.8,       # Minimum 80% success

    # Reproducibility
    seed=12345
)
```

---

## Stratified Resampling

```python
# For pooled analyses, stratify by study
bootstrap_config = BootstrapConfig(
    n_bootstrap=1000,
    stratify_by=["study", "formulation"],
    seed=12345
)
```

---

## CI Methods

### Percentile CI (Default)

```python
bootstrap_config = BootstrapConfig(
    ci_method="percentile",
    ci_level=0.95
)
```

### BCa (Bias-Corrected and Accelerated)

Better for skewed distributions:

```python
bootstrap_config = BootstrapConfig(
    ci_method="bca",
    ci_level=0.95
)
```

### Basic Bootstrap CI

```python
bootstrap_config = BootstrapConfig(
    ci_method="basic",
    ci_level=0.95
)
```

---

## Accessing Results

```python
result = run_bootstrap(data, model, config, initial_result, bootstrap_config)

# Fixed effects
print(f"Theta mean: {result.theta_mean}")
print(f"Theta SE: {result.theta_se}")
print(f"Theta RSE%: {result.theta_rse}")
print(f"95% CI lower: {result.theta_ci_lower}")
print(f"95% CI upper: {result.theta_ci_upper}")

# Bias
print(f"Bias: {result.bias}")
print(f"Bias-corrected: {result.bias_corrected}")

# Diagnostics
print(f"Success rate: {result.success_rate:.1%}")
print(f"N successful: {result.n_successful}")
print(f"N failed: {result.n_failed}")

# All estimates (for custom analysis)
all_estimates = result.theta_estimates  # (n_bootstrap, n_params)
```

---

## Quality Checks

```python
# Check success rate (FDA requires ≥80%)
if result.success_rate < 0.8:
    print(f"WARNING: Success rate {result.success_rate:.1%} below 80%")

# Check for outliers
if result.outlier_indices:
    print(f"Outlier estimates at indices: {result.outlier_indices}")

# Check SE stability
print(f"RSE of SE estimate: {result.rse_stability}")
```

---

## Visualization

```python
import matplotlib.pyplot as plt
import numpy as np

# Plot bootstrap distributions
fig, axes = plt.subplots(1, len(result.theta_mean), figsize=(4*len(result.theta_mean), 4))

param_names = ["CL", "V"]
for i, (ax, name) in enumerate(zip(axes, param_names)):
    estimates = result.theta_estimates[:, i]

    ax.hist(estimates, bins=50, density=True, alpha=0.7)
    ax.axvline(result.theta_mean[i], color='r', linestyle='-', label='Mean')
    ax.axvline(result.theta_ci_lower[i], color='r', linestyle='--', label='95% CI')
    ax.axvline(result.theta_ci_upper[i], color='r', linestyle='--')
    ax.set_xlabel(name)
    ax.set_ylabel('Density')
    ax.set_title(f'{name} Bootstrap Distribution')
    ax.legend()

plt.tight_layout()
plt.show()
```

---

## Example: Complete Bootstrap Analysis

```python
from neopkpd.estimation import (
    estimate, run_bootstrap,
    EstimationConfig, BootstrapConfig,
    FOCEIMethod
)

# Initial estimation
config = EstimationConfig(
    method=FOCEIMethod(),
    theta_init=[10.0, 50.0],
    omega_init=[[0.09, 0], [0, 0.04]],
    sigma_init=0.1,
    compute_se=True
)

initial_result = estimate(data, "OneCompIVBolus", config)

# Bootstrap analysis
bootstrap_config = BootstrapConfig(
    n_bootstrap=1000,
    bootstrap_type="case",
    ci_level=0.95,
    ci_method="percentile",
    parallel=True,
    seed=42
)

boot_result = run_bootstrap(
    data, "OneCompIVBolus", config,
    initial_result, bootstrap_config
)

# Report
print("=== Bootstrap Results ===")
print(f"Success rate: {boot_result.success_rate:.1%}")

param_names = ["CL (L/h)", "V (L)"]
print("\nParameter  Estimate    SE       RSE%     95% CI")
print("-" * 55)

for i, name in enumerate(param_names):
    est = initial_result.theta[i]
    se = boot_result.theta_se[i]
    rse = boot_result.theta_rse[i]
    ci_lo = boot_result.theta_ci_lower[i]
    ci_hi = boot_result.theta_ci_upper[i]
    print(f"{name:<10} {est:<10.3f} {se:<8.3f} {rse:<8.1f} [{ci_lo:.3f}, {ci_hi:.3f}]")

# Compare asymptotic vs bootstrap SE
print("\n=== SE Comparison ===")
for i, name in enumerate(param_names):
    asymp = initial_result.theta_se[i]
    boot = boot_result.theta_se[i]
    ratio = boot / asymp
    print(f"{name}: Asymptotic={asymp:.3f}, Bootstrap={boot:.3f}, Ratio={ratio:.2f}")
```

---

## Best Practices

1. **Sample size**: Use n ≥ 500 (FDA) or 1000 for regulatory submissions
2. **Success rate**: Target ≥ 80%
3. **Stratification**: Always stratify for pooled analyses
4. **CI method**: Use BCa for skewed parameters

---

## See Also

- [FOCE-I Method](foce.md) - Primary estimation
- [Diagnostics](diagnostics.md) - Model validation
- [Model Comparison](comparison.md) - Comparing models

