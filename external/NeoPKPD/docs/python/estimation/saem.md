# SAEM Algorithm

Stochastic Approximation Expectation Maximization (SAEM) is a robust estimation method using MCMC sampling for complex models.

---

## Overview

SAEM alternates between MCMC sampling of random effects and stochastic parameter updates, making it more robust for challenging datasets.

```python
from neopkpd.estimation import estimate, EstimationConfig, SAEMMethod

config = EstimationConfig(
    method=SAEMMethod(
        n_burn=200,
        n_iter=300,
        n_chains=3
    ),
    theta_init=[10.0, 50.0],
    omega_init=[[0.09, 0], [0, 0.04]],
    sigma_init=0.1
)

result = estimate(data, "OneCompIVBolus", config)
```

---

## When to Use SAEM

| Situation | Why SAEM Helps |
|-----------|----------------|
| Large IIV (>50% CV) | Better handles extreme random effects |
| Complex nonlinear models | More robust to model nonlinearity |
| FOCE convergence issues | SAEM often converges where FOCE fails |
| Multimodal likelihood | Multiple chains explore the space |

---

## Configuration

### SAEMMethod Parameters

```python
from neopkpd.estimation import SAEMMethod

method = SAEMMethod(
    # Iteration control
    n_burn=200,              # Burn-in iterations
    n_iter=300,              # Main iterations

    # MCMC settings
    n_chains=3,              # Chains per subject
    n_mcmc_steps=50,         # MCMC steps per E-step

    # Proposal adaptation
    adapt_proposal=True,     # Adaptive Metropolis-Hastings
    target_acceptance=0.3,   # Target acceptance rate
    adaptation_interval=20,  # Adapt every N iterations

    # Step size
    step_size_schedule="standard",  # "standard", "fast", or "slow"

    # Diagnostics
    track_diagnostics=True,  # Track MCMC convergence
    use_all_chains=True,     # Average across all chains
    parallel_chains=False    # Parallel chain execution
)
```

### Full Configuration

```python
from neopkpd.estimation import (
    EstimationConfig, SAEMMethod, BLQConfig, BLQMethod
)

config = EstimationConfig(
    method=SAEMMethod(
        n_burn=300,
        n_iter=500,
        n_chains=4,
        n_mcmc_steps=100,
        target_acceptance=0.25,
        track_diagnostics=True
    ),
    theta_init=[15.0, 30.0, 5.0, 100.0, 1.5],
    omega_init=[
        [0.09, 0, 0, 0, 0],
        [0, 0.04, 0, 0, 0],
        [0, 0, 0.16, 0, 0],
        [0, 0, 0, 0.04, 0],
        [0, 0, 0, 0, 0.25]
    ],
    sigma_init=0.15,
    blq_config=BLQConfig(method=BLQMethod.M3, lloq=0.05),
    compute_se=True,
    verbose=True,
    seed=12345
)
```

---

## MCMC Diagnostics

### Accessing Diagnostics

```python
result = estimate(data, model, config)

# SAEM-specific diagnostics
diag = result.saem_diagnostics

print(f"Mean acceptance rate: {diag.mean_acceptance_rate:.2%}")
print(f"Gelman-Rubin R-hat: {diag.gelman_rubin}")
print(f"Effective sample size: {diag.effective_sample_size}")
print(f"Converged: {diag.converged}")
```

### Convergence Checks

```python
# Check R-hat (should be < 1.1)
max_rhat = max(diag.gelman_rubin)
if max_rhat > 1.1:
    print(f"WARNING: R-hat = {max_rhat:.3f} > 1.1")

# Check ESS (should be > 100)
min_ess = min(diag.effective_sample_size)
if min_ess < 100:
    print(f"WARNING: ESS = {min_ess:.0f} < 100")
```

### Trace Plots

```python
import matplotlib.pyplot as plt

# Plot parameter traces
fig, axes = plt.subplots(2, 1, figsize=(10, 6))

# Theta trace
theta_trace = np.array(diag.theta_trace)
for i in range(theta_trace.shape[1]):
    axes[0].plot(theta_trace[:, i], label=f'θ{i+1}')
axes[0].set_xlabel('Iteration')
axes[0].set_ylabel('Theta')
axes[0].legend()
axes[0].axvline(config.method.n_burn, color='r', linestyle='--', label='End burn-in')

# OFV trace
axes[1].plot(diag.ofv_trace)
axes[1].set_xlabel('Iteration')
axes[1].set_ylabel('OFV')

plt.tight_layout()
plt.show()
```

---

## Standard Errors for SAEM

SAEM doesn't directly produce analytical SEs. Use bootstrap:

```python
from neopkpd.estimation import run_bootstrap, BootstrapConfig

# Run bootstrap after SAEM
bootstrap_config = BootstrapConfig(
    n_bootstrap=500,
    parallel=True,
    seed=12345
)

bootstrap_result = run_bootstrap(
    data, model, config, result, bootstrap_config
)

print(f"Bootstrap SE: {bootstrap_result.theta_se}")
print(f"95% CI: [{bootstrap_result.theta_ci_lower}, {bootstrap_result.theta_ci_upper}]")
```

---

## Parallel Execution

### Parallel Subjects

```python
config = EstimationConfig(
    method=SAEMMethod(n_chains=3),
    parallel=True,     # Parallel subject processing
    n_threads=8,       # Number of threads
    # ... other options
)
```

### Parallel Chains

```python
method = SAEMMethod(
    n_chains=4,
    parallel_chains=True  # Run chains in parallel
)
```

---

## Example: Complex Model with High IIV

```python
from neopkpd.estimation import (
    estimate, EstimationConfig, SAEMMethod,
    BLQConfig, BLQMethod
)

# Complex model scenario
config = EstimationConfig(
    method=SAEMMethod(
        n_burn=400,          # More burn-in for complex models
        n_iter=600,          # More iterations
        n_chains=5,          # More chains for better mixing
        n_mcmc_steps=100,    # More MCMC steps
        target_acceptance=0.25,
        track_diagnostics=True
    ),
    theta_init=[5.0, 50.0, 0.5],    # CL, V, ktr
    omega_init=[
        [0.36, 0, 0],    # 60% CV on CL
        [0, 0.16, 0],    # 40% CV on V
        [0, 0, 0.49]     # 70% CV on ktr
    ],
    sigma_init={"additive": 0.1, "proportional": 0.15},
    blq_config=BLQConfig(method=BLQMethod.M3, lloq=0.01),
    compute_se=True,
    verbose=True,
    seed=42
)

result = estimate(data, "TransitAbsorption", config)

# Check convergence
diag = result.saem_diagnostics
print("=== SAEM Convergence ===")
print(f"Converged: {result.converged}")
print(f"Mean acceptance: {diag.mean_acceptance_rate:.2%}")
print(f"Max R-hat: {max(diag.gelman_rubin):.3f}")
print(f"Min ESS: {min(diag.effective_sample_size):.0f}")

# Parameter estimates
print("\n=== Parameter Estimates ===")
param_names = ["CL", "V", "ktr"]
for i, name in enumerate(param_names):
    cv = np.sqrt(result.omega[i, i]) * 100
    print(f"{name} = {result.theta[i]:.3f} (IIV: {cv:.1f}%)")
```

---

## Comparison with FOCE-I

| Aspect | FOCE-I | SAEM |
|--------|--------|------|
| Eta estimation | Mode (optimization) | Posterior samples (MCMC) |
| Local minima | Susceptible | More robust |
| Computation | Faster | Slower |
| Standard errors | Analytical | Bootstrap/Louis |
| High IIV | Can struggle | Handles well |

---

## See Also

- [FOCE-I Method](foce.md) - Faster alternative
- [Laplacian Method](laplacian.md) - For sparse data
- [Bootstrap](bootstrap.md) - Uncertainty for SAEM
- [Diagnostics](diagnostics.md) - Model validation

