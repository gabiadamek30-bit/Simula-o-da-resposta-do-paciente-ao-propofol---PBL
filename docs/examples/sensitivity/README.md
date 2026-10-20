# Sensitivity Analysis Examples

Examples demonstrating parameter sensitivity analysis for pharmacometric models.

## Examples

| Example | Description | Directory |
|---------|-------------|-----------|
| Single Subject | Parameter perturbation | [01_single_subject](01_single_subject/README.md) |
| Population | Population-level sensitivity | [02_population](02_population/README.md) |
| Tornado Plot | Visual sensitivity display | [03_tornado_plot](03_tornado_plot/README.md) |

## What is Sensitivity Analysis?

Sensitivity analysis examines how changes in model parameters affect outputs:

- **Local sensitivity**: Small perturbations around nominal values
- **Global sensitivity**: Large parameter ranges, interactions
- **Parameter ranking**: Which parameters most affect outputs

## Key Metrics

| Metric | Description | Use |
|--------|-------------|-----|
| AUC | Area under curve | Exposure |
| Cmax | Maximum concentration | Safety |
| Cmin | Trough concentration | Efficacy |
| t_half | Half-life | Dosing interval |

## Usage

```python
from neopkpd import compute_sensitivity

result = compute_sensitivity(
    model_spec,
    parameters=["CL", "V", "Ka"],
    perturbation=0.10,  # ±10%
    metrics=["auc", "cmax", "tmax"]
)

print(result.sensitivity_matrix)
```

## See Also

- [Models](../models/README.md) - Model specifications
- [Population](../population/README.md) - Population modeling
