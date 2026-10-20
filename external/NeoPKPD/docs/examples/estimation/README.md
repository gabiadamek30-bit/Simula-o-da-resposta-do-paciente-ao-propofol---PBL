# Parameter Estimation Examples

Examples demonstrating non-linear mixed effects (NLME) parameter estimation.

## Estimation Methods

| Method | Description | Best For |
|--------|-------------|----------|
| FOCE-I | First-Order Conditional Estimation with Interaction | Most PK/PD models, standard choice |
| SAEM | Stochastic Approximation Expectation Maximization | Complex models, multi-modal posteriors |
| Laplacian | Laplace approximation | Sparse data, categorical outcomes |

## Examples

| Example | Description | Directory |
|---------|-------------|-----------|
| FOCE-I Basic | One-compartment estimation | [01_foce_onecomp](01_foce_onecomp/README.md) |
| SAEM Two-Comp | Two-compartment with SAEM | [02_saem_twocomp](02_saem_twocomp/README.md) |
| Laplacian | Sparse pediatric data | [03_laplacian_sparse](03_laplacian_sparse/README.md) |
| Diagnostics | GOF, CWRES, ETAs | [04_estimation_diagnostics](04_estimation_diagnostics/README.md) |
| Model Comparison | AIC/BIC selection | [05_model_comparison](05_model_comparison/README.md) |

## Estimation Workflow

```
1. Data Preparation
   ├── Format observed data (CDISC or custom)
   └── Define dosing events

2. Model Specification
   ├── Choose structural model
   ├── Define IIV structure
   └── Specify residual error model

3. Initial Values
   ├── Set θ₀ (fixed effects)
   ├── Set Ω₀ (random effects variance)
   └── Set σ₀ (residual error)

4. Run Estimation
   ├── Choose method (FOCE-I, SAEM, Laplacian)
   └── Set convergence criteria

5. Evaluate Results
   ├── Convergence diagnostics
   ├── Parameter estimates ± SE
   ├── Goodness-of-fit plots
   └── Model selection (AIC/BIC)
```

## Key Parameters

### Estimation Configuration

```julia
config = EstimationConfig(
    method = FOCEI(),
    # Fixed effects
    theta_init = [5.0, 50.0, 1.5],
    theta_lower = [0.1, 1.0, 0.1],
    theta_upper = [100.0, 500.0, 10.0],
    # Random effects
    omega_init = [0.09, 0.04, 0.16],
    omega_fixed = [false, false, false],
    # Residual error
    sigma_init = [0.01],
    # Algorithm settings
    maxiter = 500,
    tol = 1e-6
)
```

### Error Models

| Model | Formula | Use Case |
|-------|---------|----------|
| Additive | Y = F + ε | Low concentrations |
| Proportional | Y = F × (1 + ε) | Wide concentration range |
| Combined | Y = F × (1 + ε₁) + ε₂ | Best of both |

## Output Interpretation

### Parameter Estimates

```
Parameter   Estimate    SE        RSE%      95% CI
CL          5.23        0.31      5.9       [4.62, 5.84]
V           48.7        3.2       6.6       [42.4, 55.0]
Ka          1.42        0.18      12.7      [1.07, 1.77]
ω_CL        0.32        0.04      12.5      [0.24, 0.40]
ω_V         0.21        0.03      14.3      [0.15, 0.27]
σ_prop      0.12        0.01      8.3       [0.10, 0.14]
```

### Objective Function

- OFV (Objective Function Value): -2 × log-likelihood
- Lower OFV = better fit
- ΔOFV > 3.84 (df=1, α=0.05) = significant improvement

## Common Issues

| Issue | Symptom | Solution |
|-------|---------|----------|
| Non-convergence | Max iterations reached | Better initial values, simpler model |
| Boundary estimates | ω → 0 or → ∞ | Fix parameter, check data |
| High RSE | SE > 50% | More data, simpler model |
| Poor GOF | Systematic bias | Check structural model |

## See Also

- [Population Examples](../population/README.md) - Population simulation
- [VPC](../vpc/README.md) - Visual Predictive Checks
- [NCA](../nca/README.md) - Non-compartmental analysis
