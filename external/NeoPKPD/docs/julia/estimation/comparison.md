# Model Comparison

Tools for comparing nested and non-nested models to select the best-fitting model for your data.

---

## Overview

Model comparison helps answer:
- Is the more complex model significantly better?
- Which model best balances fit and parsimony?
- Should I include this random effect or covariate?

---

## Objective Function Value (OFV)

The OFV is -2 × log-likelihood:

$$OFV = -2 \cdot \ln L(\theta, \Omega, \sigma | y)$$

Lower OFV = better fit (but beware overfitting).

```julia
result = foce_estimate(...)
println("OFV: ", result.ofv)
```

---

## Likelihood Ratio Test (LRT)

For **nested models** only (reduced model is special case of full model).

### Mathematical Basis

Test statistic:

$$\chi^2 = OFV_{reduced} - OFV_{full}$$

Under H₀ (reduced model is adequate), this follows χ² with df = difference in parameters.

### Usage

```julia
# Estimate both models
result_full = foce_estimate(observed, model_full, config_full, grid, solver, rng)
result_reduced = foce_estimate(observed, model_reduced, config_reduced, grid, solver, rng)

# Likelihood ratio test
chi_sq, p_value = likelihood_ratio_test(
    result_full.ofv,
    result_reduced.ofv,
    df = 1  # Number of parameters removed
)

println("Chi-squared: $chi_sq")
println("p-value: $p_value")

if p_value < 0.05
    println("Full model significantly better (p < 0.05)")
    println("Recommendation: Use full model")
else
    println("No significant improvement with full model")
    println("Recommendation: Use simpler (reduced) model")
end
```

### Significance Thresholds

| ΔOFV | df | p-value | Decision |
|------|-----|---------|----------|
| 3.84 | 1 | 0.05 | Significant |
| 6.63 | 1 | 0.01 | Highly significant |
| 10.83 | 1 | 0.001 | Very highly significant |
| 5.99 | 2 | 0.05 | Significant |
| 7.81 | 3 | 0.05 | Significant |

### Common LRT Applications

```julia
# Test IIV on CL
# Full: OMEGA on CL and V
# Reduced: OMEGA on V only (CL IIV fixed to 0)
chi_sq, p = likelihood_ratio_test(ofv_full, ofv_reduced, df=1)

# Test covariate effect
# Full: CL = θ₁ × (WT/70)^θ₂
# Reduced: CL = θ₁ (no weight effect)
chi_sq, p = likelihood_ratio_test(ofv_full, ofv_reduced, df=1)

# Test correlation between etas
# Full: Block OMEGA (CL, V correlated)
# Reduced: Diagonal OMEGA (independent)
chi_sq, p = likelihood_ratio_test(ofv_full, ofv_reduced, df=1)
```

---

## Information Criteria

For **both nested and non-nested** model comparison.

### Akaike Information Criterion (AIC)

$$AIC = OFV + 2p$$

Where p = number of estimated parameters.

```julia
aic = compute_aic(result.ofv, n_params)
println("AIC: ", result.aic)
```

### Bayesian Information Criterion (BIC)

$$BIC = OFV + p \cdot \ln(n)$$

Where n = number of observations.

```julia
bic = compute_bic(result.ofv, n_params, n_observations)
println("BIC: ", result.bic)
```

### AIC vs BIC

| Criterion | Penalty | Best For |
|-----------|---------|----------|
| AIC | 2 per parameter | Prediction |
| BIC | ln(n) per parameter | Model selection |

BIC penalizes complexity more heavily, especially for large datasets.

### Interpretation

| ΔAIC | Evidence |
|------|----------|
| 0-2 | Weak |
| 2-4 | Positive |
| 4-7 | Strong |
| >10 | Very strong |

```julia
# Compare models
delta_aic = result1.aic - result2.aic

if delta_aic > 10
    println("Strong evidence for Model 2")
elseif delta_aic > 4
    println("Good evidence for Model 2")
elseif delta_aic > 2
    println("Weak evidence for Model 2")
else
    println("Models are similar - prefer simpler")
end
```

---

## Model Comparison Table

```julia
# Compare multiple models
models = [
    ("1-comp", result_1comp),
    ("2-comp", result_2comp),
    ("2-comp + CL~WT", result_2comp_wt)
]

println("Model                n_params    OFV        AIC        BIC")
println("-" ^ 65)

for (name, result) in models
    np = count_parameters(result)
    println("$(rpad(name, 20)) $np           $(round(result.ofv, digits=2))     $(round(result.aic, digits=2))     $(round(result.bic, digits=2))")
end
```

---

## Covariate Model Selection

### Forward Selection

Add covariates one at a time:

```julia
# Start with base model
base_ofv = base_result.ofv
candidate_covariates = [:WT, :AGE, :SEX, :CRCL]

println("=== Forward Selection ===")
println("Base OFV: $base_ofv")

for cov in candidate_covariates
    # Fit model with this covariate
    result_with_cov = estimate_with_covariate(cov)
    delta_ofv = base_ofv - result_with_cov.ofv

    if delta_ofv > 3.84  # p < 0.05
        println("$cov: ΔOFV = $delta_ofv, p < 0.05 ✓")
    else
        println("$cov: ΔOFV = $delta_ofv, NS")
    end
end
```

### Backward Elimination

Remove covariates from full model:

```julia
# Start with full model (all covariates)
full_ofv = full_result.ofv

println("=== Backward Elimination ===")
println("Full OFV: $full_ofv")

for cov in included_covariates
    # Fit model without this covariate
    result_without_cov = estimate_without_covariate(cov)
    delta_ofv = result_without_cov.ofv - full_ofv

    if delta_ofv > 6.63  # p < 0.01 (stricter for removal)
        println("$cov: ΔOFV = $delta_ofv, KEEP (p < 0.01)")
    else
        println("$cov: ΔOFV = $delta_ofv, DROP")
    end
end
```

### Stepwise Covariate Modeling (SCM)

```julia
# Automated forward-backward procedure
scm_result = stepwise_covariate_modeling(
    observed,
    base_model,
    candidate_covariates = [:WT, :AGE, :SEX, :CRCL],
    forward_p = 0.05,      # p-value for inclusion
    backward_p = 0.01,     # p-value for exclusion
    max_iterations = 20
)

println("Final covariates: ", scm_result.included_covariates)
println("Final OFV: ", scm_result.final_ofv)
```

---

## Random Effects Selection

### Test IIV on Parameter

```julia
# Should we include IIV on ka?

# Model 1: IIV on CL, V only
config_no_ka_iiv = EstimationConfig(
    omega_init = diagm([0.09, 0.04]),  # 2 random effects
    # ...
)

# Model 2: IIV on CL, V, ka
config_with_ka_iiv = EstimationConfig(
    omega_init = diagm([0.09, 0.04, 0.25]),  # 3 random effects
    # ...
)

result1 = foce_estimate(..., config_no_ka_iiv, ...)
result2 = foce_estimate(..., config_with_ka_iiv, ...)

# Test
chi_sq, p = likelihood_ratio_test(result2.ofv, result1.ofv, df=1)
println("IIV on ka: p = $p")
```

### Test Correlation Between Random Effects

```julia
# Should CL and V be correlated?

# Model 1: Diagonal omega
config_diag = EstimationConfig(
    omega_init = diagm([0.09, 0.04]),
    omega_structure = :diagonal
)

# Model 2: Block omega (CL-V correlated)
config_block = EstimationConfig(
    omega_init = [0.09 0.02; 0.02 0.04],
    omega_structure = :block
)

result_diag = foce_estimate(..., config_diag, ...)
result_block = foce_estimate(..., config_block, ...)

# Test (1 df for correlation parameter)
chi_sq, p = likelihood_ratio_test(result_block.ofv, result_diag.ofv, df=1)
println("CL-V correlation: p = $p")

if p < 0.05
    corr = result_block.omega[1,2] / sqrt(result_block.omega[1,1] * result_block.omega[2,2])
    println("Estimated correlation: $corr")
end
```

---

## Residual Error Model Selection

```julia
# Compare additive vs proportional vs combined error

# Additive
config_add = EstimationConfig(
    sigma_init = ResidualErrorSpec(AdditiveError(), AdditiveErrorParams(0.5), :conc)
)

# Proportional
config_prop = EstimationConfig(
    sigma_init = ResidualErrorSpec(ProportionalError(), ProportionalErrorParams(0.1), :conc)
)

# Combined (2 parameters)
config_comb = EstimationConfig(
    sigma_init = ResidualErrorSpec(CombinedError(), CombinedErrorParams(0.5, 0.1), :conc)
)

# Fit all
result_add = foce_estimate(..., config_add, ...)
result_prop = foce_estimate(..., config_prop, ...)
result_comb = foce_estimate(..., config_comb, ...)

# Compare using AIC (non-nested models)
println("Error Model    OFV        AIC")
println("-" ^ 40)
println("Additive       $(result_add.ofv)    $(result_add.aic)")
println("Proportional   $(result_prop.ofv)    $(result_prop.aic)")
println("Combined       $(result_comb.ofv)    $(result_comb.aic)")

# Combined vs Additive (nested)
chi_sq, p = likelihood_ratio_test(result_comb.ofv, result_add.ofv, df=1)
println("\nCombined vs Additive: p = $p")
```

---

## Model Comparison Functions

### Core Functions

```julia
# Likelihood ratio test
chi_sq, p_value = likelihood_ratio_test(ofv_full, ofv_reduced, df)

# Information criteria
aic = compute_aic(ofv, n_params)
bic = compute_bic(ofv, n_params, n_obs)

# Log-likelihood
ll = loglikelihood_from_ofv(ofv)  # Returns -OFV/2
```

### Parameter Counting

```julia
# Count parameters for information criteria
n_params = count_parameters(result)
# Includes: theta + omega parameters + sigma parameters

# Detailed breakdown
n_theta = length(result.theta)
n_omega = count_omega_params(result.omega, omega_structure)
n_sigma = count_sigma_params(result.sigma)
```

---

## Example: Complete Model Selection

```julia
using NeoPKPD
using StableRNGs

# Data
observed = load_observed_data("pk_data.csv")
grid = SimGrid(0.0, 72.0, 0:0.5:72)
solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^7)
rng = StableRNG(42)

# ============================================
# Step 1: Structural Model Selection
# ============================================
println("=== Structural Model Selection ===")

# 1-compartment
model_1c = OneCompIVBolus()
result_1c = foce_estimate(observed, ModelSpec(model_1c, "pk", nothing, nothing),
    EstimationConfig(theta_init=[10.0, 50.0], omega_init=diagm([0.09, 0.04]), ...),
    grid, solver, rng)

# 2-compartment
model_2c = TwoCompIVBolus()
result_2c = foce_estimate(observed, ModelSpec(model_2c, "pk", nothing, nothing),
    EstimationConfig(theta_init=[10.0, 30.0, 5.0, 100.0], omega_init=diagm([0.09, 0.04, 0.16, 0.04]), ...),
    grid, solver, rng)

# Non-nested: use AIC
println("1-comp AIC: $(result_1c.aic)")
println("2-comp AIC: $(result_2c.aic)")
println("ΔAIC = $(result_1c.aic - result_2c.aic)")

# Select 2-comp if AIC lower by >4
best_struct = result_2c.aic < result_1c.aic - 4 ? "2-comp" : "1-comp"
println("Selected: $best_struct")

# ============================================
# Step 2: Random Effects Model
# ============================================
println("\n=== Random Effects Selection ===")

# Start with IIV on all structural parameters
base_model = best_struct == "2-comp" ? model_2c : model_1c

# Test diagonal vs block omega for CL-V
result_diag = estimate_with_diagonal_omega(...)
result_block = estimate_with_block_omega(...)

chi_sq, p = likelihood_ratio_test(result_block.ofv, result_diag.ofv, df=1)
println("CL-V correlation test: p = $p")

# ============================================
# Step 3: Covariate Selection
# ============================================
println("\n=== Covariate Selection ===")

# Forward selection
candidates = [:WT, :AGE, :SEX, :CRCL]
selected = Symbol[]
current_ofv = result_block.ofv

for cov in candidates
    result_cov = estimate_with_covariate(base_model, cov)
    delta = current_ofv - result_cov.ofv

    if delta > 3.84
        println("$cov: ΔOFV = $delta, INCLUDE")
        push!(selected, cov)
        current_ofv = result_cov.ofv
    else
        println("$cov: ΔOFV = $delta, EXCLUDE")
    end
end

# ============================================
# Step 4: Residual Error Model
# ============================================
println("\n=== Residual Error Selection ===")

result_prop = estimate_with_proportional_error(...)
result_comb = estimate_with_combined_error(...)

chi_sq, p = likelihood_ratio_test(result_comb.ofv, result_prop.ofv, df=1)
println("Combined vs Prop: p = $p")

# ============================================
# Final Model Summary
# ============================================
println("\n=== Final Model ===")
println("Structural: $best_struct")
println("Random effects: $(result_block.converged ? "Block CL-V" : "Diagonal")")
println("Covariates: $selected")
println("Residual error: $(p < 0.05 ? "Combined" : "Proportional")")
println("\nFinal OFV: $(final_result.ofv)")
println("Final AIC: $(final_result.aic)")
println("Final BIC: $(final_result.bic)")
```

---

## Best Practices

1. **Start simple**: Begin with simplest model and add complexity
2. **Nested vs non-nested**: Use LRT for nested, AIC/BIC for non-nested
3. **Multiple comparisons**: Adjust significance threshold if many tests
4. **Clinical plausibility**: Statistical significance ≠ clinical importance
5. **Validate**: Check diagnostics after selection

---

## See Also

- [FOCE-I Method](foce.md) - Estimation method
- [Diagnostics](diagnostics.md) - Model validation
- [Bootstrap](bootstrap.md) - Uncertainty quantification

