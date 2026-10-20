# Population PKPD Analysis Use Case

Complete population pharmacokinetic-pharmacodynamic (PKPD) analysis workflow
demonstrating drug exposure-response modeling with indirect response dynamics.

## Model Structure

### PK Component
- **Model**: Two-compartment oral
- **Parameters**: Ka, CL, V1, Q, V2
- **Covariates**: Weight on CL and V1

### PD Component
- **Model**: Indirect response (inhibition of production)
- **Parameters**: Kin, Kout, Imax, IC50
- **Baseline**: E0 = Kin/Kout

## Workflow Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `run.jl` | Complete PKPD simulation and analysis |

## Quick Run

```bash
# Run complete PKPD workflow
julia --project=packages/core docs/examples/use_cases/population_pkpd_analysis/run.jl

# Validate outputs
julia --project=packages/core docs/examples/use_cases/population_pkpd_analysis/validate.jl
```

## Population

| Parameter | Description |
|-----------|-------------|
| N | 50 subjects |
| Weight | 50-100 kg (uniform) |
| Dose | 50 mg oral |
| Sampling | PK: 0-24h, PD: 0-72h |

## Expected Results

### PK Parameters (Population)
| Parameter | Typical | IIV (CV%) |
|-----------|---------|-----------|
| Ka | 1.0 /h | 30 |
| CL | 5.0 L/h | 25 |
| V1 | 20 L | 20 |
| Q | 2.0 L/h | - |
| V2 | 40 L | - |

### PD Parameters (Population)
| Parameter | Typical | IIV (CV%) |
|-----------|---------|-----------|
| Kin | 1.0 /h | 20 |
| Kout | 0.1 /h | 20 |
| Imax | 0.9 | 15 |
| IC50 | 2.0 mg/L | 30 |

### Key Findings
- Cmax (PK): ~5 mg/L
- E_max inhibition: ~45% from baseline
- Time to max PD effect: ~8-12 hours post-dose
- PD recovery: ~80% by 72 hours

## See Also

- [Population Documentation](../../../python/population/index.md)
- [VPC Documentation](../../../python/vpc/index.md)
