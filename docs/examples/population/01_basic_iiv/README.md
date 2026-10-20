# Basic Inter-Individual Variability (IIV)

Population simulation with log-normal IIV on PK parameters.

## Model

```
θᵢ = θ_pop × exp(ηᵢ)

Where ηᵢ ~ N(0, ω²)
```

## Parameters

| Parameter | Population | ω (SD) | CV% |
|-----------|------------|--------|-----|
| CL | 5.0 L/h | 0.3 | ~30% |
| V | 50.0 L | 0.2 | ~20% |

## What This Example Shows

1. Generate 100 virtual subjects
2. Each subject has individual CL and V values
3. Simulate PK for all subjects
4. Compute population summaries (mean, median, 5th/95th percentiles)

## Expected Variability

With ω_CL = 0.3:
- 95% of subjects have CL between ~2.7 and ~9.2 L/h
- This is 0.55× to 1.8× the typical value

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
