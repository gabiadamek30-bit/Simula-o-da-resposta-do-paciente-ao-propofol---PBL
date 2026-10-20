# Real-World Theophylline Analysis

Complete pharmacokinetic analysis workflow using the classic theophylline dataset
from the nlmixr2 package.

## Dataset

- **Source**: Boeckmann AJ, Sheiner LB, Beal SL (1994). NONMEM Users Guide: Part V
- **Subjects**: 12 healthy adults
- **Route**: Oral administration
- **Observations**: ~11 concentration measurements per subject over 24h
- **Covariates**: Body weight (WT)

## Workflow Steps

This use case demonstrates a complete PK analysis workflow:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_data_exploration.jl` | Load and explore the dataset |
| 2 | `02_nca_analysis.jl` | Non-compartmental analysis |
| 3 | `03_model_fitting.jl` | One-compartment oral model fitting |
| 4 | `04_diagnostics.jl` | Goodness-of-fit and residual diagnostics |
| 5 | `05_vpc.jl` | Visual predictive check |
| 6 | `06_report.jl` | Generate summary report |

## Quick Run

```bash
# Run complete workflow
julia --project=packages/core docs/examples/use_cases/real_world_theophylline/run.jl

# Validate outputs
julia --project=packages/core docs/examples/use_cases/real_world_theophylline/validate.jl
```

## Expected Results

### NCA Metrics (Population Mean)
| Metric | Value | Units |
|--------|-------|-------|
| Cmax | 8.2 | mg/L |
| Tmax | 1.5 | h |
| AUC_0_24 | 114 | mg*h/L |
| t_half | 8.5 | h |

### Model Parameters (Population Estimates)
| Parameter | Estimate | CV% |
|-----------|----------|-----|
| Ka | 1.5 h⁻¹ | 25 |
| CL | 0.04 L/h/kg | 20 |
| V | 0.5 L/kg | 15 |

## See Also

- [NCA Documentation](../../../julia/nca/index.md)
- [Estimation Documentation](../../../julia/estimation/index.md)
- [VPC Documentation](../../../julia/vpc/index.md)
