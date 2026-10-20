# Model Comparison

Tools for comparing nested and non-nested models.

---

## Overview

```python
from neopkpd.estimation import likelihood_ratio_test, compare_models

# Compare two models
chi_sq, p_value = likelihood_ratio_test(
    ofv_full=result_full.ofv,
    ofv_reduced=result_reduced.ofv,
    df=1
)
```

---

## Likelihood Ratio Test

For **nested models** (reduced is special case of full):

```python
from neopkpd.estimation import likelihood_ratio_test

# Test if adding IIV on ka improves fit
chi_sq, p_value = likelihood_ratio_test(
    ofv_full=result_with_ka_iiv.ofv,
    ofv_reduced=result_without_ka_iiv.ofv,
    df=1  # One parameter added
)

print(f"Chi-squared: {chi_sq:.2f}")
print(f"p-value: {p_value:.4f}")

if p_value < 0.05:
    print("Full model significantly better")
else:
    print("Prefer simpler model")
```

### Significance Thresholds

| ΔOFV | df | p-value |
|------|-----|---------|
| 3.84 | 1 | 0.05 |
| 6.63 | 1 | 0.01 |
| 5.99 | 2 | 0.05 |
| 7.81 | 3 | 0.05 |

---

## Information Criteria

For **both nested and non-nested** models:

```python
# Access from results
print(f"Model 1: OFV={result1.ofv:.2f}, AIC={result1.aic:.2f}, BIC={result1.bic:.2f}")
print(f"Model 2: OFV={result2.ofv:.2f}, AIC={result2.aic:.2f}, BIC={result2.bic:.2f}")

# Lower is better
delta_aic = result1.aic - result2.aic
print(f"ΔAIC = {delta_aic:.2f}")
```

### AIC vs BIC

$$AIC = OFV + 2p$$
$$BIC = OFV + p \cdot \ln(n)$$

| Criterion | Penalty | Best For |
|-----------|---------|----------|
| AIC | 2 per param | Prediction |
| BIC | ln(n) per param | Model selection |

---

## Compare Multiple Models

```python
from neopkpd.estimation import compare_models

# Compare multiple models at once
comparison = compare_models(
    results=[result_1comp, result_2comp, result_2comp_wt],
    names=["1-comp", "2-comp", "2-comp + CL~WT"]
)

print(comparison.summary_table())
# Output:
# Model          n_params   OFV      AIC      BIC      ΔAIC
# ---------------------------------------------------------
# 1-comp         5          245.3    255.3    267.1    0.0
# 2-comp         9          198.7    216.7    240.3    -38.6
# 2-comp+CL~WT   10         192.1    212.1    238.5    -43.2

print(f"Best model by AIC: {comparison.best_by_aic}")
print(f"Best model by BIC: {comparison.best_by_bic}")
```

---

## Covariate Model Selection

### Forward Selection

```python
def forward_selection(base_result, data, model, covariates, threshold=3.84):
    """Forward covariate selection."""
    selected = []
    current_ofv = base_result.ofv

    for cov in covariates:
        # Fit model with this covariate
        result_with_cov = estimate_with_covariate(data, model, cov)
        delta_ofv = current_ofv - result_with_cov.ofv

        if delta_ofv > threshold:  # p < 0.05
            print(f"{cov}: ΔOFV = {delta_ofv:.2f}, INCLUDE")
            selected.append(cov)
            current_ofv = result_with_cov.ofv
        else:
            print(f"{cov}: ΔOFV = {delta_ofv:.2f}, EXCLUDE")

    return selected

covariates = ["WT", "AGE", "SEX", "CRCL"]
selected = forward_selection(base_result, data, model, covariates)
print(f"Selected covariates: {selected}")
```

### Backward Elimination

```python
def backward_elimination(full_result, data, model, covariates, threshold=6.63):
    """Backward covariate elimination."""
    remaining = list(covariates)
    current_ofv = full_result.ofv

    for cov in covariates:
        # Fit model without this covariate
        result_without = estimate_without_covariate(data, model, cov)
        delta_ofv = result_without.ofv - current_ofv

        if delta_ofv < threshold:  # p > 0.01
            print(f"{cov}: ΔOFV = {delta_ofv:.2f}, DROP")
            remaining.remove(cov)
        else:
            print(f"{cov}: ΔOFV = {delta_ofv:.2f}, KEEP")

    return remaining
```

---

## Random Effects Selection

### Test IIV

```python
# Test if IIV on ka is needed
result_with_ka = estimate(data, model, config_with_ka_iiv)
result_without_ka = estimate(data, model, config_without_ka_iiv)

chi_sq, p = likelihood_ratio_test(
    result_with_ka.ofv,
    result_without_ka.ofv,
    df=1
)

print(f"IIV on ka: p = {p:.4f}")
if p < 0.05:
    print("Include IIV on ka")
```

### Test Correlation

```python
# Test CL-V correlation
result_block = estimate(data, model, config_block_omega)
result_diag = estimate(data, model, config_diagonal_omega)

chi_sq, p = likelihood_ratio_test(result_block.ofv, result_diag.ofv, df=1)

if p < 0.05:
    corr = result_block.omega_corr[0, 1]
    print(f"CL-V correlation significant: r = {corr:.3f}")
```

---

## Example: Complete Model Selection

```python
from neopkpd.estimation import (
    estimate, likelihood_ratio_test, compare_models,
    EstimationConfig, FOCEIMethod
)

# Step 1: Structural model selection
print("=== Structural Model Selection ===")

result_1c = estimate(data, "OneCompIVBolus", config_1c)
result_2c = estimate(data, "TwoCompIVBolus", config_2c)

print(f"1-comp: AIC = {result_1c.aic:.2f}")
print(f"2-comp: AIC = {result_2c.aic:.2f}")

best_struct = "2-comp" if result_2c.aic < result_1c.aic - 4 else "1-comp"
print(f"Selected: {best_struct}")

# Step 2: Random effects
print("\n=== Random Effects Selection ===")

result_diag = estimate(data, model, config_diag)
result_block = estimate(data, model, config_block)

chi_sq, p = likelihood_ratio_test(result_block.ofv, result_diag.ofv, df=1)
print(f"CL-V correlation: p = {p:.4f}")

# Step 3: Covariates
print("\n=== Covariate Selection ===")

covariates = ["WT", "AGE", "SEX"]
for cov in covariates:
    result_cov = estimate_with_covariate(data, model, cov)
    delta = base_ofv - result_cov.ofv
    status = "INCLUDE" if delta > 3.84 else "EXCLUDE"
    print(f"{cov}: ΔOFV = {delta:.2f}, {status}")

# Step 4: Final comparison
print("\n=== Final Model Comparison ===")
comparison = compare_models(
    [result_base, result_final],
    ["Base", "Final"]
)
print(comparison.summary_table())
```

---

## See Also

- [FOCE-I Method](foce.md) - Estimation method
- [Diagnostics](diagnostics.md) - Model validation
- [Bootstrap](bootstrap.md) - Uncertainty quantification

