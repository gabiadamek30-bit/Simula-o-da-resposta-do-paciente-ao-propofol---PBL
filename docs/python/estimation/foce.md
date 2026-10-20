# FOCE-I Method

First-Order Conditional Estimation with Interaction (FOCE-I) is the industry-standard method for nonlinear mixed-effects modeling.

---

## Overview

FOCE-I linearizes the model around individual eta estimates, providing accurate parameter estimates for most PK/PD models.

```python
from neopkpd.estimation import estimate, EstimationConfig, FOCEIMethod

config = EstimationConfig(
    method=FOCEIMethod(),
    theta_init=[10.0, 50.0],
    omega_init=[[0.09, 0], [0, 0.04]],
    sigma_init=0.1
)

result = estimate(data, "OneCompIVBolus", config)
```

---

## Configuration

### FOCEIMethod Parameters

```python
from neopkpd.estimation import FOCEIMethod

method = FOCEIMethod(
    centered=False,        # False=FOCE-I, True=FOCE (no interaction)
    compute_cwres=True,    # Compute conditional weighted residuals
    laplacian=True,        # Include Laplacian correction
    eta_tol=1e-8,          # Inner optimization tolerance
    max_eta_iter=100       # Max iterations for eta optimization
)
```

### Full Configuration Example

```python
from neopkpd.estimation import (
    EstimationConfig, FOCEIMethod, BLQConfig, BLQMethod,
    OmegaStructure, ResidualErrorModel
)

config = EstimationConfig(
    # Method specification
    method=FOCEIMethod(
        centered=False,
        compute_cwres=True
    ),

    # Initial estimates
    theta_init=[10.0, 50.0],           # CL, V
    omega_init=[[0.09, 0], [0, 0.04]], # IIV variances

    # Residual error
    sigma_init=0.1,                    # Proportional error
    residual_error_model=ResidualErrorModel.PROPORTIONAL,

    # Omega structure
    omega_structure=OmegaStructure.DIAGONAL,

    # BLQ handling
    blq_config=BLQConfig(
        method=BLQMethod.M3,
        lloq=0.1
    ),

    # Convergence settings
    max_iter=1000,
    tol=1e-6,

    # Output options
    compute_se=True,
    compute_ci=True,
    ci_level=0.95,
    verbose=True,
    seed=12345
)
```

---

## Residual Error Models

### Additive Error

$$Y = F + \epsilon, \quad \epsilon \sim N(0, \sigma^2)$$

```python
config = EstimationConfig(
    sigma_init=0.5,
    residual_error_model=ResidualErrorModel.ADDITIVE
)
```

### Proportional Error

$$Y = F \cdot (1 + \epsilon), \quad \epsilon \sim N(0, \sigma^2)$$

```python
config = EstimationConfig(
    sigma_init=0.1,  # 10% CV
    residual_error_model=ResidualErrorModel.PROPORTIONAL
)
```

### Combined Error

$$Y = F + F \cdot \epsilon_{prop} + \epsilon_{add}$$

```python
config = EstimationConfig(
    sigma_init={"additive": 0.5, "proportional": 0.1},
    residual_error_model=ResidualErrorModel.COMBINED
)
```

### Exponential Error

$$\log(Y) = \log(F) + \epsilon$$

```python
config = EstimationConfig(
    sigma_init=0.15,
    residual_error_model=ResidualErrorModel.EXPONENTIAL
)
```

---

## Omega Structures

### Diagonal Omega

Independent random effects:

```python
config = EstimationConfig(
    omega_init=[
        [0.09, 0],
        [0, 0.04]
    ],
    omega_structure=OmegaStructure.DIAGONAL
)
```

### Block Omega

Correlated random effects:

```python
config = EstimationConfig(
    omega_init=[
        [0.09, 0.02],
        [0.02, 0.04]
    ],
    omega_structure=OmegaStructure.BLOCK
)
```

### Full Omega

All parameters correlated:

```python
config = EstimationConfig(
    omega_init=[
        [0.09, 0.02, 0.01],
        [0.02, 0.04, 0.005],
        [0.01, 0.005, 0.16]
    ],
    omega_structure=OmegaStructure.FULL
)
```

---

## BLQ Handling

```python
from neopkpd.estimation import BLQConfig, BLQMethod

# M1: Discard BLQ
blq_config = BLQConfig(method=BLQMethod.M1, lloq=0.1)

# M2: Impute at LLOQ/2
blq_config = BLQConfig(method=BLQMethod.M2, lloq=0.1, impute_value="half")

# M3: Censored likelihood (recommended)
blq_config = BLQConfig(method=BLQMethod.M3, lloq=0.1)

config = EstimationConfig(
    blq_config=blq_config,
    # ... other options
)
```

---

## Inter-Occasion Variability (IOV)

```python
from neopkpd.estimation import IOVSpec

iov_specs = [
    IOVSpec(
        eta_name="eta_CL",
        occasion_names=["OCC1", "OCC2", "OCC3"],
        omega_iov=0.04  # 20% IOV on CL
    )
]

config = EstimationConfig(
    iov_specs=iov_specs,
    # ... other options
)
```

---

## Covariate Effects on IIV

```python
from neopkpd.estimation import CovariateOnIIV

covariate_effects = [
    CovariateOnIIV(
        eta_name="eta_CL",
        covariate_name="WT",
        effect_type="exponential",
        reference_value=70.0
    ),
    CovariateOnIIV(
        eta_name="eta_V",
        covariate_name="AGE",
        effect_type="linear",
        reference_value=40.0
    )
]

config = EstimationConfig(
    covariate_effects=covariate_effects,
    # ... other options
)
```

---

## Accessing Results

```python
result = estimate(data, model, config)

# Fixed effects
print(f"Theta: {result.theta}")
print(f"SE: {result.theta_se}")
print(f"RSE%: {result.theta_rse}")
print(f"95% CI: [{result.theta_ci_lower}, {result.theta_ci_upper}]")

# Random effects
print(f"Omega:\n{result.omega}")
print(f"Omega correlation:\n{result.omega_corr}")

# Residual error
print(f"Sigma: {result.sigma}")

# Model fit
print(f"OFV: {result.ofv}")
print(f"AIC: {result.aic}")
print(f"BIC: {result.bic}")

# Convergence
print(f"Converged: {result.converged}")
print(f"Iterations: {result.n_iterations}")
print(f"Runtime: {result.runtime:.2f}s")

# Diagnostics
print(f"Eta shrinkage: {result.eta_shrinkage}")
print(f"Epsilon shrinkage: {result.epsilon_shrinkage}")
print(f"Condition number: {result.condition_number}")
```

---

## Individual Estimates

```python
for ind in result.individual_estimates:
    print(f"\nSubject {ind.subject_id}:")
    print(f"  Eta: {ind.eta}")
    print(f"  Individual CL: {result.theta[0] * np.exp(ind.eta[0]):.3f}")
    print(f"  IPRED: {ind.ipred[:3]}...")
    print(f"  CWRES: {ind.cwres[:3]}...")
```

---

## Standard Errors

### Asymptotic SEs

```python
config = EstimationConfig(
    compute_se=True,
    se_method="hessian"  # Default
)
```

### Sandwich (Robust) SEs

```python
config = EstimationConfig(
    compute_se=True,
    se_method="sandwich"  # More robust to misspecification
)
```

---

## Example: Two-Compartment Model

```python
from neopkpd.estimation import (
    estimate, EstimationConfig, FOCEIMethod,
    BLQConfig, BLQMethod, EstimationData
)
import pandas as pd

# Load data
df = pd.read_csv("two_comp_data.csv")

data = EstimationData.from_dataframe(
    df, id_col="ID", time_col="TIME",
    dv_col="DV", amt_col="AMT"
)

# Configure estimation
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
    max_iter=2000,
    compute_se=True,
    verbose=True
)

# Run estimation
result = estimate(data, "TwoCompIVBolus", config)

# Print summary
param_names = ["CL (L/h)", "V1 (L)", "Q (L/h)", "V2 (L)"]
print("=== Parameter Estimates ===")
print(f"{'Parameter':<12} {'Estimate':<10} {'SE':<10} {'RSE%':<8} {'95% CI'}")
print("-" * 60)

for i, name in enumerate(param_names):
    est = result.theta[i]
    se = result.theta_se[i]
    rse = result.theta_rse[i]
    ci = f"[{result.theta_ci_lower[i]:.2f}, {result.theta_ci_upper[i]:.2f}]"
    print(f"{name:<12} {est:<10.3f} {se:<10.3f} {rse:<8.1f} {ci}")

print(f"\nOFV: {result.ofv:.2f}")
print(f"AIC: {result.aic:.2f}")
```

---

## Mathematical Details

### FOCE-I Objective Function

$$OFV = \sum_i \left[ \ln|C_i(\hat{\eta}_i)| + (y_i - f_i(\hat{\eta}_i))^T C_i^{-1}(\hat{\eta}_i) (y_i - f_i(\hat{\eta}_i)) + \hat{\eta}_i^T \Omega^{-1} \hat{\eta}_i + \ln|H_{\eta,i}| \right]$$

### FOCE vs FOCE-I

| Aspect | FOCE | FOCE-I |
|--------|------|--------|
| Interaction | No | Yes |
| C evaluation | At η=0 | At η=η̂ |
| Accuracy | Lower | Higher |

---

## See Also

- [SAEM Algorithm](saem.md) - Alternative for complex models
- [Laplacian Method](laplacian.md) - For sparse data
- [Diagnostics](diagnostics.md) - Model validation
- [Bootstrap](bootstrap.md) - Uncertainty quantification

