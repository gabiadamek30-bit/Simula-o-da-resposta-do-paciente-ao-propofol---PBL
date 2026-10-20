# Model Diagnostics

Comprehensive diagnostic tools for assessing model fit, validating assumptions, and identifying potential issues in population PK/PD analyses.

---

## Overview

Model diagnostics help answer:
- Does the model adequately describe the data?
- Are the random effects properly specified?
- Are there systematic biases or trends?
- Which subjects/observations are influential?

---

## Residual Types

### CWRES - Conditional Weighted Residuals

The gold standard for NLME diagnostics:

$$CWRES = \frac{y - E[y|\hat{\eta}]}{\sqrt{Var(y|\hat{\eta})}}$$

CWRES accounts for:
- Individual predictions (using $\hat{\eta}$)
- Proper variance structure
- Gradient corrections (for FOCE-I)

```julia
result = foce_estimate(...)

# Access CWRES per subject
for ind in result.individual_estimates
    println("Subject $(ind.subject_id): CWRES = $(ind.cwres)")
end

# All CWRES combined
all_cwres = vcat([ind.cwres for ind in result.individual_estimates]...)
```

### IWRES - Individual Weighted Residuals

Simpler than CWRES, uses individual predictions:

$$IWRES = \frac{y - f(\hat{\eta})}{\sigma(f(\hat{\eta}))}$$

```julia
# Access IWRES
for ind in result.individual_estimates
    println("Subject $(ind.subject_id): IWRES = $(ind.iwres)")
end
```

### WRES - Population Weighted Residuals

Uses population predictions ($\eta = 0$):

$$WRES = \frac{y - f(0)}{\sigma(f(0))}$$

```julia
wres = compute_wres(observed, predictions_eta_zero, sigma)
```

### NPDE - Normalized Prediction Distribution Errors

Simulation-based diagnostic (most rigorous):

$$NPDE = \Phi^{-1}(pde_i)$$

Where $pde_i$ is the fraction of simulated values below the observation.

```julia
# Compute NPDE from simulations
npde = compute_npde(observed_values, simulated_matrix, rng)

# NPDE should follow N(0,1) if model is correct
mean(npde)  # Should be ~0
std(npde)   # Should be ~1
```

---

## Diagnostic Criteria

### Residual Expectations

| Residual | Mean | SD | Distribution |
|----------|------|-----|--------------|
| CWRES | 0 | 1 | N(0,1) |
| IWRES | 0 | 1 | N(0,1) |
| NPDE | 0 | 1 | N(0,1) |

### Warning Signs

| Issue | Symptom | Possible Cause |
|-------|---------|----------------|
| Bias | Mean ≠ 0 | Structural model misspecification |
| Over/under-dispersion | SD ≠ 1 | Residual error misspecification |
| Trends vs time | Pattern in plot | Time-varying effect |
| Trends vs pred | Pattern in plot | Nonlinearity |
| Heavy tails | Non-normal QQ | Outliers or wrong error model |

---

## Shrinkage

### Eta Shrinkage

Measures how much individual estimates regress toward population:

$$\text{Shrinkage}_\eta = 1 - \frac{SD(\hat{\eta}_i)}{\omega}$$

```julia
# From estimation result
println("Eta shrinkage: ", result.eta_shrinkage)

# Manual computation
shrinkage = shrinkage_eta(result.individual_estimates, result.omega)
```

#### Interpretation

| Shrinkage | Interpretation |
|-----------|---------------|
| < 20% | Good - individual estimates well-supported |
| 20-30% | Moderate - use caution |
| 30-50% | High - individual estimates uncertain |
| > 50% | Very high - rely on population estimates |

High shrinkage indicates:
- Sparse data per subject
- Low IIV (hard to distinguish subjects)
- Model overparameterized

### Epsilon Shrinkage

Measures information in residual error:

$$\text{Shrinkage}_\epsilon = 1 - SD(IWRES)$$

```julia
# Compute epsilon shrinkage
all_iwres = vcat([ind.iwres for ind in result.individual_estimates]...)
eps_shrinkage = shrinkage_epsilon(all_iwres)
println("Epsilon shrinkage: $(eps_shrinkage * 100)%")
```

High epsilon shrinkage suggests:
- Sparse data
- Overparameterized residual error model
- Individual predictions nearly match observations

---

## Goodness-of-Fit Diagnostics

### Standard GOF Plots

```julia
using Plots

result = foce_estimate(...)

# Collect data
obs = Float64[]
pred = Float64[]
ipred = Float64[]
cwres = Float64[]
time = Float64[]

for ind in result.individual_estimates
    append!(obs, ind.observed)
    append!(pred, ind.pred)
    append!(ipred, ind.ipred)
    append!(cwres, ind.cwres)
    append!(time, ind.times)
end

# 1. DV vs PRED
p1 = scatter(pred, obs, xlabel="Population Prediction", ylabel="Observed",
             title="DV vs PRED", legend=false)
plot!(p1, [minimum(pred), maximum(pred)], [minimum(pred), maximum(pred)],
      linewidth=2, linestyle=:dash)

# 2. DV vs IPRED
p2 = scatter(ipred, obs, xlabel="Individual Prediction", ylabel="Observed",
             title="DV vs IPRED", legend=false)
plot!(p2, [minimum(ipred), maximum(ipred)], [minimum(ipred), maximum(ipred)],
      linewidth=2, linestyle=:dash)

# 3. CWRES vs Time
p3 = scatter(time, cwres, xlabel="Time", ylabel="CWRES",
             title="CWRES vs Time", legend=false)
hline!(p3, [-2, 0, 2], linestyle=:dash)

# 4. CWRES vs PRED
p4 = scatter(pred, cwres, xlabel="Population Prediction", ylabel="CWRES",
             title="CWRES vs PRED", legend=false)
hline!(p4, [-2, 0, 2], linestyle=:dash)

# Combined 4-panel plot
plot(p1, p2, p3, p4, layout=(2,2), size=(800, 800))
```

### QQ Plot for Residuals

```julia
using StatsPlots

# QQ plot
qqnorm(cwres, title="QQ Plot of CWRES")
```

---

## Residual Summary Statistics

```julia
# Get residual summary
summary = residual_summary(all_cwres)

println("CWRES Summary:")
println("  Mean: $(summary[:mean])")
println("  SD: $(summary[:sd])")
println("  Median: $(summary[:median])")
println("  Min: $(summary[:min])")
println("  Max: $(summary[:max])")
println("  Skewness: $(summary[:skewness])")
println("  Kurtosis: $(summary[:kurtosis])")
```

---

## Individual Fit Diagnostics

### Individual Plots

```julia
# Plot individual fits
for (i, ind) in enumerate(result.individual_estimates[1:9])
    p = plot(title="Subject $(ind.subject_id)")
    scatter!(p, ind.times, ind.observed, label="Observed", markersize=6)
    plot!(p, ind.times, ind.ipred, label="IPRED", linewidth=2)
    plot!(p, ind.times, ind.pred, label="PRED", linestyle=:dash)
end
```

### Individual Parameters

```julia
# Examine individual estimates
for ind in result.individual_estimates
    println("Subject $(ind.subject_id):")
    println("  eta = $(ind.eta)")
    println("  Individual CL = $(result.theta[1] * exp(ind.eta[1]))")
    println("  Individual V  = $(result.theta[2] * exp(ind.eta[2]))")
end
```

---

## Covariance Diagnostics

### Condition Number

Large condition numbers indicate near-singularity:

```julia
println("Condition number: ", result.condition_number)

# Interpretation
# < 100: Good
# 100-1000: Monitor
# > 1000: Problematic - consider fixing parameters
```

### Eigenvalue Analysis

```julia
eigenvalues = result.eigenvalues
println("Eigenvalue ratio: ", maximum(eigenvalues) / minimum(eigenvalues))

# Check for near-zero eigenvalues
if minimum(eigenvalues) < 1e-6
    println("Warning: Near-singular covariance matrix")
end
```

### Correlation Matrix

```julia
# Parameter correlation matrix
corr_matrix = cov_to_corr(result.covariance_matrix)

# High correlations (>0.9) suggest identifiability issues
for i in 1:size(corr_matrix, 1)
    for j in i+1:size(corr_matrix, 2)
        if abs(corr_matrix[i,j]) > 0.9
            println("High correlation: param $i - param $j: $(corr_matrix[i,j])")
        end
    end
end
```

---

## Model Comparison

### Objective Function Value

```julia
# Compare nested models
ofv_full = result_full.ofv
ofv_reduced = result_reduced.ofv

delta_ofv = ofv_reduced - ofv_full
println("ΔOFV = ", delta_ofv)
```

### Likelihood Ratio Test

```julia
# For nested models: test if full model is significantly better
chi_sq, p_value = likelihood_ratio_test(
    ofv_full,
    ofv_reduced,
    df = 1  # Difference in parameters
)

println("Chi-squared: $chi_sq, p-value: $p_value")

if p_value < 0.05
    println("Full model significantly better (p < 0.05)")
else
    println("No significant difference - prefer simpler model")
end
```

### Information Criteria

```julia
println("Model 1: AIC = $(result1.aic), BIC = $(result1.bic)")
println("Model 2: AIC = $(result2.aic), BIC = $(result2.bic)")

# Lower is better
# AIC: Good for prediction
# BIC: Good for model selection (stronger penalty)
```

### AIC/BIC Formulas

$$AIC = OFV + 2p$$
$$BIC = OFV + p \cdot \ln(n)$$

Where:
- $p$ = number of parameters
- $n$ = number of observations

---

## Influential Subject Analysis

### Jackknife Influence

```julia
# Compute influence of each subject
influence = jackknife_influence(result, observed, model_spec, config, grid, solver)

# Find most influential subjects
sorted_idx = sortperm(influence, rev=true)
println("Most influential subjects:")
for i in sorted_idx[1:5]
    println("  Subject $(observed.subjects[i].id): influence = $(influence[i])")
end
```

### Cook's Distance

```julia
# Cook's distance for each subject
cooks_d = compute_cooks_distance(result, observed)

# Subjects with Cook's D > 4/n are potentially influential
threshold = 4 / n_subjects(observed)
influential = findall(cooks_d .> threshold)
println("Potentially influential subjects: ", influential)
```

---

## VPC-Based Diagnostics

### VPC Check Statistic

```julia
# Percentage of observations within prediction interval
vpc_stat = vpc_check(observed_values, pi_lower, pi_upper)
println("$(vpc_stat)% of observations within 90% PI")

# Should be close to 90% for 90% PI
```

---

## Diagnostic Functions Reference

### Residual Functions

| Function | Description |
|----------|-------------|
| `compute_cwres(obs, ipred, sigma, gradients)` | Conditional weighted residuals |
| `compute_iwres(obs, ipred, sigma)` | Individual weighted residuals |
| `compute_wres(obs, pred, sigma)` | Population weighted residuals |
| `compute_npde(obs, simulated, rng)` | Normalized prediction distribution errors |

### Shrinkage Functions

| Function | Description |
|----------|-------------|
| `shrinkage_eta(etas, omega)` | Eta shrinkage per parameter |
| `shrinkage_epsilon(iwres)` | Epsilon shrinkage |

### Model Comparison Functions

| Function | Description |
|----------|-------------|
| `likelihood_ratio_test(ofv1, ofv2, df)` | LRT between nested models |
| `compute_aic(ofv, n_params)` | Akaike Information Criterion |
| `compute_bic(ofv, n_params, n_obs)` | Bayesian Information Criterion |

---

## Diagnostic Checklist

### Pre-Estimation
- [ ] Initial estimates reasonable?
- [ ] Data formatted correctly?
- [ ] BLQ handling specified?

### Post-Estimation
- [ ] Convergence achieved?
- [ ] Condition number acceptable (<1000)?
- [ ] RSE < 50% for all parameters?
- [ ] Shrinkage < 30%?
- [ ] CWRES mean ≈ 0, SD ≈ 1?
- [ ] No trends in residual plots?
- [ ] QQ plot approximately linear?
- [ ] Individual fits adequate?

---

## See Also

- [FOCE-I Method](foce.md) - Primary estimation method
- [SAEM Algorithm](saem.md) - Alternative estimation
- [Bootstrap](bootstrap.md) - Uncertainty quantification
- [Model Comparison](comparison.md) - Comparing models

