# Model Diagnostics

Comprehensive diagnostic tools for assessing model fit and validating assumptions.

---

## Overview

```python
from neopkpd.estimation import compute_diagnostics

result = estimate(data, model, config)
diagnostics = compute_diagnostics(result)
```

---

## Residual Types

### CWRES - Conditional Weighted Residuals

```python
# Access CWRES
for ind in result.individual_estimates:
    print(f"Subject {ind.subject_id}: CWRES = {ind.cwres}")

# All CWRES combined
import numpy as np
all_cwres = np.concatenate([ind.cwres for ind in result.individual_estimates])
print(f"CWRES mean: {np.mean(all_cwres):.3f}")
print(f"CWRES std: {np.std(all_cwres):.3f}")
```

### IWRES - Individual Weighted Residuals

```python
all_iwres = np.concatenate([ind.iwres for ind in result.individual_estimates])
```

### Expected Values

| Residual | Expected Mean | Expected SD |
|----------|---------------|-------------|
| CWRES | 0 | 1 |
| IWRES | 0 | 1 |
| NPDE | 0 | 1 |

---

## Shrinkage

### Eta Shrinkage

```python
# Access eta shrinkage
print(f"Eta shrinkage: {result.eta_shrinkage}")

# Interpretation
for i, shrink in enumerate(result.eta_shrinkage):
    pct = shrink * 100
    status = "Good" if pct < 20 else ("Moderate" if pct < 30 else "High")
    print(f"  η{i+1}: {pct:.1f}% ({status})")
```

| Shrinkage | Interpretation |
|-----------|----------------|
| < 20% | Good |
| 20-30% | Moderate |
| > 30% | High - use caution |

### Epsilon Shrinkage

```python
print(f"Epsilon shrinkage: {result.epsilon_shrinkage*100:.1f}%")
```

---

## Goodness-of-Fit Plots

```python
import matplotlib.pyplot as plt
import numpy as np

# Collect data
obs = np.concatenate([ind.observed for ind in result.individual_estimates])
pred = np.concatenate([ind.pred for ind in result.individual_estimates])
ipred = np.concatenate([ind.ipred for ind in result.individual_estimates])
cwres = np.concatenate([ind.cwres for ind in result.individual_estimates])
time = np.concatenate([ind.times for ind in result.individual_estimates])

# Create 4-panel GOF plot
fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# DV vs PRED
axes[0, 0].scatter(pred, obs, alpha=0.5)
axes[0, 0].plot([pred.min(), pred.max()], [pred.min(), pred.max()], 'r--')
axes[0, 0].set_xlabel('Population Prediction')
axes[0, 0].set_ylabel('Observed')
axes[0, 0].set_title('DV vs PRED')

# DV vs IPRED
axes[0, 1].scatter(ipred, obs, alpha=0.5)
axes[0, 1].plot([ipred.min(), ipred.max()], [ipred.min(), ipred.max()], 'r--')
axes[0, 1].set_xlabel('Individual Prediction')
axes[0, 1].set_ylabel('Observed')
axes[0, 1].set_title('DV vs IPRED')

# CWRES vs Time
axes[1, 0].scatter(time, cwres, alpha=0.5)
axes[1, 0].axhline(y=0, color='r', linestyle='--')
axes[1, 0].axhline(y=-2, color='gray', linestyle=':')
axes[1, 0].axhline(y=2, color='gray', linestyle=':')
axes[1, 0].set_xlabel('Time')
axes[1, 0].set_ylabel('CWRES')
axes[1, 0].set_title('CWRES vs Time')

# CWRES vs PRED
axes[1, 1].scatter(pred, cwres, alpha=0.5)
axes[1, 1].axhline(y=0, color='r', linestyle='--')
axes[1, 1].axhline(y=-2, color='gray', linestyle=':')
axes[1, 1].axhline(y=2, color='gray', linestyle=':')
axes[1, 1].set_xlabel('Population Prediction')
axes[1, 1].set_ylabel('CWRES')
axes[1, 1].set_title('CWRES vs PRED')

plt.tight_layout()
plt.show()
```

---

## QQ Plot

```python
import scipy.stats as stats

fig, ax = plt.subplots(figsize=(6, 6))
stats.probplot(cwres, dist="norm", plot=ax)
ax.set_title("QQ Plot of CWRES")
plt.show()
```

---

## Individual Fits

```python
# Plot first 9 subjects
fig, axes = plt.subplots(3, 3, figsize=(12, 12))
axes = axes.flatten()

for i, ind in enumerate(result.individual_estimates[:9]):
    ax = axes[i]
    ax.scatter(ind.times, ind.observed, label='Observed', s=50)
    ax.plot(ind.times, ind.ipred, 'b-', label='IPRED', linewidth=2)
    ax.plot(ind.times, ind.pred, 'r--', label='PRED', alpha=0.7)
    ax.set_title(f"Subject {ind.subject_id}")
    ax.set_xlabel('Time')
    ax.set_ylabel('Concentration')
    if i == 0:
        ax.legend()

plt.tight_layout()
plt.show()
```

---

## Covariance Diagnostics

```python
# Condition number
print(f"Condition number: {result.condition_number:.1f}")

if result.condition_number > 1000:
    print("WARNING: High condition number - potential identifiability issues")

# Eigenvalue ratio
print(f"Eigenvalue ratio: {result.eigenvalue_ratio:.1f}")

# Correlation matrix
if result.covariance_matrix is not None:
    corr = result.correlation_matrix
    print("\nParameter correlation matrix:")
    print(corr)
```

---

## Diagnostic Checklist

```python
def check_diagnostics(result):
    issues = []

    # Convergence
    if not result.converged:
        issues.append("Model did not converge")

    # Condition number
    if result.condition_number > 1000:
        issues.append(f"High condition number: {result.condition_number:.0f}")

    # RSE
    for i, rse in enumerate(result.theta_rse):
        if rse > 50:
            issues.append(f"High RSE for θ{i+1}: {rse:.1f}%")

    # Shrinkage
    for i, shrink in enumerate(result.eta_shrinkage):
        if shrink > 0.3:
            issues.append(f"High η{i+1} shrinkage: {shrink*100:.1f}%")

    # CWRES
    all_cwres = np.concatenate([ind.cwres for ind in result.individual_estimates])
    cwres_mean = np.mean(all_cwres)
    cwres_std = np.std(all_cwres)

    if abs(cwres_mean) > 0.2:
        issues.append(f"CWRES mean bias: {cwres_mean:.3f}")
    if abs(cwres_std - 1) > 0.2:
        issues.append(f"CWRES SD deviation: {cwres_std:.3f}")

    return issues

# Run checks
issues = check_diagnostics(result)
if issues:
    print("Diagnostic issues found:")
    for issue in issues:
        print(f"  - {issue}")
else:
    print("All diagnostics passed!")
```

---

## See Also

- [FOCE-I Method](foce.md) - Estimation method
- [Bootstrap](bootstrap.md) - Uncertainty quantification
- [Model Comparison](comparison.md) - Comparing models

