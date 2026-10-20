# Bioequivalence Study Use Case

Complete workflow for a 2x2 crossover bioequivalence (BE) study comparing
a test formulation against a reference formulation.

## Study Design

| Parameter | Value |
|-----------|-------|
| Design | 2-period, 2-sequence crossover |
| Subjects per sequence | 12 (24 total) |
| Washout | 7 days |
| Sampling | 0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 24 h |
| BE criteria | 90% CI for AUC and Cmax within 80-125% |

## Formulations

- **Reference (R)**: Approved product, Ka = 1.5 /h
- **Test (T)**: Generic formulation, Ka = 1.4 /h (slightly slower absorption)

## Workflow Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `run.jl` | Complete BE simulation and analysis |

## Quick Run

```bash
# Run complete BE workflow
julia --project=packages/core docs/examples/use_cases/bioequivalence_study/run.jl

# Validate outputs
julia --project=packages/core docs/examples/use_cases/bioequivalence_study/validate.jl
```

## Expected Results

### Geometric Mean Ratios (T/R)
| Metric | Point Estimate | 90% CI |
|--------|----------------|--------|
| Cmax | 0.95 | 0.88-1.02 |
| AUC_0_t | 1.00 | 0.95-1.05 |

### BE Conclusion
- Both Cmax and AUC 90% CIs within 80-125%
- **Bioequivalence demonstrated**

## Model Parameters

| Parameter | Reference | Test | IIV (CV%) |
|-----------|-----------|------|-----------|
| Ka | 1.5 /h | 1.4 /h | 30 |
| CL | 10 L/h | 10 L/h | 25 |
| V | 100 L | 100 L | 20 |

## See Also

- [Trial Examples](../../trial/README.md)
- [NCA Examples](../../nca/README.md)
