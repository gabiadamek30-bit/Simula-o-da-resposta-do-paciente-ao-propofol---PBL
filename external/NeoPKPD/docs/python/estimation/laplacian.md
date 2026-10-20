# Laplacian Method

The Laplacian approximation provides fast, efficient parameter estimation particularly suited for sparse data.

---

## Overview

```python
from neopkpd.estimation import estimate, EstimationConfig, LaplacianMethod

config = EstimationConfig(
    method=LaplacianMethod(),
    theta_init=[10.0, 50.0],
    omega_init=[[0.09, 0], [0, 0.04]],
    sigma_init=0.1
)

result = estimate(data, "OneCompIVBolus", config)
```

---

## When to Use Laplacian

| Situation | Why Laplacian Helps |
|-----------|---------------------|
| Sparse sampling | Few observations per subject (1-3) |
| Large populations | Fast computation scales well |
| Initial estimates | Quick screening before FOCE/SAEM |
| Simple models | One-compartment, few parameters |

---

## Configuration

### LaplacianMethod Parameters

```python
from neopkpd.estimation import LaplacianMethod

method = LaplacianMethod(
    max_inner_iter=50,     # Max iterations for eta mode
    inner_tol=1e-6,        # Inner optimization tolerance
    use_prior=True,        # Include prior on eta
    hessian_method="exact" # or "finite_diff"
)
```

### Full Configuration

```python
from neopkpd.estimation import (
    EstimationConfig, LaplacianMethod, BLQConfig, BLQMethod
)

config = EstimationConfig(
    method=LaplacianMethod(
        max_inner_iter=50,
        inner_tol=1e-6
    ),
    theta_init=[10.0, 50.0],
    omega_init=[[0.25, 0], [0, 0.16]],  # Wide initial IIV
    sigma_init={"additive": 0.5, "proportional": 0.2},
    blq_config=BLQConfig(method=BLQMethod.M3, lloq=0.1),
    max_iter=500,
    compute_se=True,
    verbose=True
)
```

---

## Sparse Data Example

```python
from neopkpd.estimation import (
    estimate, EstimationConfig, LaplacianMethod, EstimationData
)

# Very sparse data: 1-2 observations per subject
data = EstimationData(
    subject_ids=["1", "1", "2", "3", "3"],
    times=[1.0, 4.0, 2.0, 0.5, 8.0],
    observations=[8.5, 3.2, 6.1, 9.2, 1.1],
    doses=[
        {"time": 0.0, "amount": 500.0, "subject_id": "1"},
        {"time": 0.0, "amount": 500.0, "subject_id": "2"},
        {"time": 0.0, "amount": 500.0, "subject_id": "3"}
    ]
)

config = EstimationConfig(
    method=LaplacianMethod(),
    theta_init=[5.0, 30.0],
    omega_init=[[0.25, 0], [0, 0.16]],
    sigma_init=0.15,
    compute_se=True
)

result = estimate(data, "OneCompIVBolus", config)

print(f"CL = {result.theta[0]:.3f} ± {result.theta_se[0]:.3f}")
print(f"V  = {result.theta[1]:.3f} ± {result.theta_se[1]:.3f}")
print(f"OFV = {result.ofv:.2f}")
```

---

## Comparison with FOCE-I

| Aspect | Laplacian | FOCE-I |
|--------|-----------|--------|
| Interaction term | No | Yes |
| Computational cost | Lower | Higher |
| Accuracy (dense) | Lower | Higher |
| Accuracy (sparse) | Good | Good |

---

## See Also

- [FOCE-I Method](foce.md) - More accurate for dense data
- [SAEM Algorithm](saem.md) - More robust for complex models
- [Diagnostics](diagnostics.md) - Model validation

