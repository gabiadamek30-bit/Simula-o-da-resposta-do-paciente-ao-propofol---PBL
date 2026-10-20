# FOCE-I One-Compartment Estimation

Parameter estimation using FOCE-I method on simulated oral PK data.

## Model

**Structural:**
```
dA_depot/dt = -Ka × A_depot
dA_central/dt = Ka × A_depot - CL/V × A_central
```

**Statistical:**
```
CL_i = θ_CL × exp(η_CL,i)
V_i  = θ_V × exp(η_V,i)
Ka_i = θ_Ka × exp(η_Ka,i)

Y_ij = F_ij × (1 + ε_ij)   # Proportional error
```

## Data

Simulated data for 30 subjects:
- 100 mg oral dose at t=0
- Sampling at 0.5, 1, 2, 4, 6, 8, 12, 24h
- True parameters: CL=5, V=50, Ka=1.5
- IIV: ω_CL=0.3, ω_V=0.2, ω_Ka=0.4

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia estimation |
| `python.py` | Python estimation |
| `data.csv` | Observed data |
| `expected_results.json` | Expected parameter estimates |
