# Complex Covariate Model

Full population model with multiple covariates, IIV, correlations, and residual error.

## Model

```
CL = CL_pop × (WT/70)^0.75 × (CrCL/100)^0.5 × θ_sex^SEX × exp(η_CL)
V  = V_pop × (WT/70)^1.0 × exp(η_V)
Ka = Ka_pop × exp(η_Ka)

Where:
- η ~ MVN(0, Ω)  with correlation between η_CL and η_V
- ε ~ N(0, σ²)   residual error
```

## Covariates

| Covariate | Effect On | Type | Value |
|-----------|-----------|------|-------|
| Weight | CL, V | Power | 0.75, 1.0 |
| CrCL | CL | Power | 0.5 |
| Sex | CL | Categorical | 0.85 (female) |

## Omega Matrix (with correlations)

```
        CL    V     Ka
CL   [ 0.09  0.02  0    ]
V    [ 0.02  0.04  0    ]
Ka   [ 0     0     0.16 ]
```

Correlation CL-V ≈ 0.33
