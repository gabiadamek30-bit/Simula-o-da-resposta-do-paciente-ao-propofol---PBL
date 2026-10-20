# Static Covariates

Population simulation with weight and age effects on PK parameters.

## Model

**Weight on CL and V (allometric):**
```
CL = CL_pop × (WT/70)^0.75
V  = V_pop × (WT/70)^1.0
```

**Age on CL (linear decline):**
```
CL = CL_adj × (1 - 0.01 × (AGE - 40))
```

## Covariate Distribution

| Covariate | Distribution | Range |
|-----------|--------------|-------|
| Weight | Normal(70, 15) | 40-120 kg |
| Age | Normal(50, 15) | 18-80 years |

## What This Example Shows

1. Generate realistic covariate distributions
2. Apply covariate effects to parameters
3. Residual IIV after covariate adjustment
4. Variance explained by covariates
