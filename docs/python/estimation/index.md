# Parameter Estimation

NeoPKPD Python provides comprehensive parameter estimation capabilities for population PK/PD modeling through seamless integration with the Julia estimation engine.

---

## Overview

Estimate population parameters from observed data using industry-standard methods:

```python
from neopkpd.estimation import estimate, EstimationConfig, FOCEIMethod

# Configure and run estimation
config = EstimationConfig(
    method=FOCEIMethod(),
    theta_init=[10.0, 50.0],
    omega_init=[[0.09, 0], [0, 0.04]],
    sigma_init=0.1,
    compute_se=True
)

result = estimate(data, model, config)

print(f"Theta: {result.theta}")
print(f"SE: {result.theta_se}")
print(f"OFV: {result.ofv}")
```

---

## Estimation Methods

| Method | Description | Best For |
|--------|-------------|----------|
| [FOCE-I](foce.md) | First-Order Conditional Estimation with Interaction | Standard analyses |
| [SAEM](saem.md) | Stochastic Approximation EM | Complex models, high IIV |
| [Laplacian](laplacian.md) | Laplace approximation | Sparse data |

---

## Quick Start

### Installation

```bash
pip install neopkpd
```

### Basic FOCE-I Estimation

```python
from neopkpd.estimation import (
    estimate,
    EstimationConfig,
    FOCEIMethod,
    BLQConfig,
    BLQMethod,
    EstimationData
)

# Prepare data
data = EstimationData(
    subject_ids=["1", "1", "1", "2", "2", "2"],
    times=[0.5, 2.0, 8.0, 0.5, 2.0, 8.0],
    observations=[1.8, 1.2, 0.4, 2.1, 1.4, 0.5],
    doses=[
        {"time": 0.0, "amount": 100.0, "subject_id": "1"},
        {"time": 0.0, "amount": 100.0, "subject_id": "2"}
    ],
    observation_name="conc"
)

# Configure estimation
config = EstimationConfig(
    method=FOCEIMethod(),
    theta_init=[5.0, 50.0],           # Initial CL, V
    omega_init=[[0.09, 0], [0, 0.04]], # 30% CV CL, 20% CV V
    sigma_init=0.1,                    # 10% proportional error
    max_iter=1000,
    tol=1e-6,
    compute_se=True,
    compute_ci=True,
    verbose=True
)

# Run estimation
result = estimate(data, "OneCompIVBolus", config)

# Access results
print("=== Estimation Results ===")
print(f"CL = {result.theta[0]:.3f} (SE: {result.theta_se[0]:.3f})")
print(f"V  = {result.theta[1]:.3f} (SE: {result.theta_se[1]:.3f})")
print(f"OFV = {result.ofv:.2f}")
print(f"AIC = {result.aic:.2f}")
```

---

## Configuration Classes

### EstimationConfig

```python
from neopkpd.estimation import EstimationConfig, FOCEIMethod

config = EstimationConfig(
    # Method
    method=FOCEIMethod(),          # or SAEMMethod(), LaplacianMethod()

    # Initial estimates
    theta_init=[10.0, 50.0],       # Fixed effects
    omega_init=[[0.09, 0],         # Random effects variance
                [0, 0.04]],
    sigma_init=0.1,                # Residual error

    # Convergence
    max_iter=1000,                 # Maximum iterations
    tol=1e-6,                      # Convergence tolerance

    # Standard errors
    compute_se=True,               # Compute SEs
    compute_ci=True,               # Compute CIs
    ci_level=0.95,                 # 95% CI

    # Omega structure
    omega_structure="diagonal",    # "diagonal", "block", or "full"

    # BLQ handling
    blq_config=BLQConfig(
        method=BLQMethod.M3,
        lloq=0.1
    ),

    # Output
    verbose=True,
    seed=12345
)
```

### Method-Specific Configuration

```python
# FOCE-I
from neopkpd.estimation import FOCEIMethod
method = FOCEIMethod(
    centered=False,        # FOCE-I (not FOCE)
    compute_cwres=True,    # Compute CWRES
    laplacian=True         # Include Laplacian correction
)

# SAEM
from neopkpd.estimation import SAEMMethod
method = SAEMMethod(
    n_burn=200,            # Burn-in iterations
    n_iter=300,            # Main iterations
    n_chains=3,            # MCMC chains per subject
    n_mcmc_steps=50,       # MCMC steps per E-step
    target_acceptance=0.3  # Target acceptance rate
)

# Laplacian
from neopkpd.estimation import LaplacianMethod
method = LaplacianMethod(
    max_inner_iter=50,     # Inner optimization iterations
    inner_tol=1e-6         # Inner tolerance
)
```

---

## BLQ Handling

Handle Below Limit of Quantification observations:

```python
from neopkpd.estimation import BLQConfig, BLQMethod

# M1: Discard BLQ observations
blq_config = BLQConfig(method=BLQMethod.M1, lloq=0.1)

# M2: Impute at LLOQ/2
blq_config = BLQConfig(method=BLQMethod.M2, lloq=0.1, impute_value="half")

# M3: Censored likelihood (recommended)
blq_config = BLQConfig(method=BLQMethod.M3, lloq=0.1)

config = EstimationConfig(
    # ... other options
    blq_config=blq_config
)
```

---

## Estimation Result

```python
result = estimate(data, model, config)

# Fixed effects
result.theta              # Parameter estimates
result.theta_se           # Standard errors
result.theta_rse          # Relative SE (%)
result.theta_ci_lower     # Lower CI
result.theta_ci_upper     # Upper CI

# Random effects
result.omega              # Omega matrix
result.omega_se           # Omega SEs
result.omega_corr         # Correlation matrix

# Residual error
result.sigma              # Sigma estimate
result.sigma_se           # Sigma SE

# Individual estimates
result.individual_estimates  # List of IndividualEstimate objects

# Diagnostics
result.eta_shrinkage      # Eta shrinkage per parameter
result.epsilon_shrinkage  # Epsilon shrinkage
result.condition_number   # Covariance matrix condition

# Model fit
result.ofv                # Objective function value
result.aic                # Akaike Information Criterion
result.bic                # Bayesian Information Criterion

# Convergence
result.converged          # Did it converge?
result.n_iterations       # Number of iterations
result.runtime            # Execution time
```

---

## Diagnostics

### Accessing Residuals

```python
from neopkpd.estimation import compute_diagnostics

# Get diagnostics
diagnostics = compute_diagnostics(result)

# Per-subject residuals
for ind in result.individual_estimates:
    print(f"Subject {ind.subject_id}:")
    print(f"  CWRES: {ind.cwres}")
    print(f"  IWRES: {ind.iwres}")
    print(f"  IPRED: {ind.ipred}")
```

### Shrinkage

```python
# Eta shrinkage (should be < 30%)
print(f"CL shrinkage: {result.eta_shrinkage[0]*100:.1f}%")
print(f"V shrinkage: {result.eta_shrinkage[1]*100:.1f}%")

# Epsilon shrinkage
print(f"Epsilon shrinkage: {result.epsilon_shrinkage*100:.1f}%")
```

---

## Model Comparison

### Likelihood Ratio Test

```python
from neopkpd.estimation import likelihood_ratio_test

# Compare nested models
chi_sq, p_value = likelihood_ratio_test(
    ofv_full=result_full.ofv,
    ofv_reduced=result_reduced.ofv,
    df=1  # Difference in parameters
)

print(f"Chi-squared: {chi_sq:.2f}")
print(f"p-value: {p_value:.4f}")
```

### Information Criteria

```python
# Compare models using AIC/BIC
print(f"Model 1: AIC={result1.aic:.2f}, BIC={result1.bic:.2f}")
print(f"Model 2: AIC={result2.aic:.2f}, BIC={result2.bic:.2f}")

# Lower is better
delta_aic = result1.aic - result2.aic
print(f"ΔAIC = {delta_aic:.2f}")
```

---

## Bootstrap Analysis

```python
from neopkpd.estimation import run_bootstrap, BootstrapConfig

# Configure bootstrap
bootstrap_config = BootstrapConfig(
    n_bootstrap=1000,
    bootstrap_type="case",      # "case", "parametric", or "residual"
    stratify_by=["study"],      # Stratification variables
    ci_level=0.95,
    ci_method="percentile",     # "percentile", "bca", or "basic"
    parallel=True,
    seed=12345
)

# Run bootstrap
bootstrap_result = run_bootstrap(
    data=data,
    model=model,
    config=config,
    initial_result=result,
    bootstrap_config=bootstrap_config
)

# Access results
print(f"Bootstrap SE: {bootstrap_result.theta_se}")
print(f"95% CI: [{bootstrap_result.theta_ci_lower}, {bootstrap_result.theta_ci_upper}]")
print(f"Success rate: {bootstrap_result.success_rate*100:.1f}%")
```

---

## Example: Complete Analysis

```python
from neopkpd.estimation import (
    estimate, run_bootstrap,
    EstimationConfig, BootstrapConfig,
    FOCEIMethod, BLQConfig, BLQMethod
)
import pandas as pd

# Load data
df = pd.read_csv("pk_data.csv")

# Prepare estimation data
data = EstimationData.from_dataframe(
    df,
    id_col="ID",
    time_col="TIME",
    dv_col="DV",
    amt_col="AMT",
    mdv_col="MDV"
)

# Configure FOCE-I estimation
config = EstimationConfig(
    method=FOCEIMethod(compute_cwres=True),
    theta_init=[10.0, 30.0, 5.0, 100.0],  # CL, V1, Q, V2
    omega_init=[
        [0.09, 0, 0, 0],
        [0, 0.04, 0, 0],
        [0, 0, 0.16, 0],
        [0, 0, 0, 0.04]
    ],
    sigma_init={"additive": 0.1, "proportional": 0.1},
    blq_config=BLQConfig(method=BLQMethod.M3, lloq=0.01),
    compute_se=True,
    verbose=True
)

# Run estimation
result = estimate(data, "TwoCompIVBolus", config)

# Print results
print("=== FOCE-I Results ===")
param_names = ["CL", "V1", "Q", "V2"]
for i, name in enumerate(param_names):
    print(f"{name}: {result.theta[i]:.3f} (RSE: {result.theta_rse[i]:.1f}%)")

print(f"\nOFV: {result.ofv:.2f}")
print(f"Converged: {result.converged}")

# Run bootstrap for uncertainty
bootstrap_config = BootstrapConfig(
    n_bootstrap=500,
    parallel=True,
    seed=42
)

boot_result = run_bootstrap(data, "TwoCompIVBolus", config, result, bootstrap_config)

print("\n=== Bootstrap Results ===")
for i, name in enumerate(param_names):
    ci_lo = boot_result.theta_ci_lower[i]
    ci_hi = boot_result.theta_ci_upper[i]
    print(f"{name}: 95% CI [{ci_lo:.3f}, {ci_hi:.3f}]")
```

---

## API Reference

### Main Functions

| Function | Description |
|----------|-------------|
| `estimate(data, model, config)` | Run parameter estimation |
| `run_bootstrap(data, model, config, result, boot_config)` | Run bootstrap analysis |
| `likelihood_ratio_test(ofv1, ofv2, df)` | Compare nested models |
| `compare_models(results, names)` | Compare multiple models |
| `compute_diagnostics(result)` | Compute model diagnostics |

### Configuration Classes

| Class | Description |
|-------|-------------|
| `EstimationConfig` | Main estimation configuration |
| `FOCEIMethod` | FOCE-I method settings |
| `SAEMMethod` | SAEM method settings |
| `LaplacianMethod` | Laplacian method settings |
| `BLQConfig` | BLQ handling configuration |
| `BootstrapConfig` | Bootstrap configuration |

### Result Classes

| Class | Description |
|-------|-------------|
| `EstimationResult` | Full estimation results |
| `IndividualEstimate` | Per-subject estimates |
| `BootstrapResult` | Bootstrap analysis results |
| `ModelComparisonResult` | Model comparison results |

---

## See Also

- [FOCE-I Method](foce.md) - Detailed FOCE-I documentation
- [SAEM Algorithm](saem.md) - SAEM details
- [Laplacian Method](laplacian.md) - Laplacian details
- [Diagnostics](diagnostics.md) - Model diagnostics
- [Bootstrap](bootstrap.md) - Bootstrap analysis
- [Model Comparison](comparison.md) - Comparing models

