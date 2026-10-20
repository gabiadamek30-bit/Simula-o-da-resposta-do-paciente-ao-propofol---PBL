# Population Modeling Examples

Examples demonstrating inter-individual variability (IIV), inter-occasion variability (IOV), and covariate modeling.

## Examples

| Example | Description | Directory |
|---------|-------------|-----------|
| Basic IIV | Log-normal inter-individual variability | [01_basic_iiv](01_basic_iiv/README.md) |
| IOV | Inter-occasion variability | [02_iov](02_iov/README.md) |
| Static Covariates | Weight and age effects | [03_static_covariates](03_static_covariates/README.md) |
| Time-Varying Covariates | Creatinine clearance over time | [04_time_varying_covariates](04_time_varying_covariates/README.md) |
| Complex Model | Multiple covariates, IIV, correlations | [05_complex_covariate_model](05_complex_covariate_model/README.md) |

## Concepts

### Inter-Individual Variability (IIV)

IIV describes differences between subjects:

```
θᵢ = θ_pop × exp(ηᵢ)
```

Where:
- θᵢ = individual parameter
- θ_pop = population typical value
- ηᵢ ~ N(0, ω²)

### Inter-Occasion Variability (IOV)

IOV describes within-subject variability across occasions:

```
θᵢⱼ = θ_pop × exp(ηᵢ + κⱼ)
```

Where:
- κⱼ ~ N(0, π²) is the occasion-specific random effect

### Covariate Models

Covariates explain part of the variability:

**Power model (weight on CL):**
```
CL = CL_pop × (WT/70)^0.75 × exp(η_CL)
```

**Linear model (age on V):**
```
V = V_pop × (1 + θ_age × (AGE - 40)) × exp(η_V)
```

**Categorical (sex on CL):**
```
CL = CL_pop × θ_sex^SEX × exp(η_CL)
```

## File Structure

Each example contains:
```
01_example/
├── README.md       # Detailed explanation
├── julia.jl        # Julia implementation
├── python.py       # Python implementation
└── cli.json        # CLI specification (if applicable)
```

## Running Examples

```bash
# Julia
julia --project=packages/core docs/examples/population/01_basic_iiv/julia.jl

# Python
python docs/examples/population/01_basic_iiv/python.py
```

## Key Parameters

### Omega (ω) Interpretation

| ω² | ω (SD) | CV% | Interpretation |
|----|--------|-----|----------------|
| 0.01 | 0.1 | ~10% | Low variability |
| 0.09 | 0.3 | ~30% | Moderate variability |
| 0.25 | 0.5 | ~50% | High variability |
| 0.49 | 0.7 | ~70% | Very high variability |

Note: For log-normal, CV ≈ ω for ω < 0.5

### Common Covariate Exponents

| Covariate | Parameter | Exponent | Rationale |
|-----------|-----------|----------|-----------|
| Weight | CL | 0.75 | Allometric scaling |
| Weight | V | 1.0 | Linear scaling |
| Age | CL | -0.01 to -0.02 | Linear decline |
| CrCL | CL | 0.5-1.0 | Renal function |
| Sex | CL | 0.8-1.2 | Categorical |

## See Also

- [Parameter Estimation](../estimation/README.md) - Estimate population parameters
- [VPC](../vpc/README.md) - Validate population models
- [Models](../models/README.md) - Available structural models
