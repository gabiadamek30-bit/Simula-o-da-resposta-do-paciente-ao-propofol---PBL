# Standard VPC

Basic Visual Predictive Check implementation.

## Method

1. Simulate n replicates (e.g., 500) using final model estimates
2. Compute quantiles (5th, 50th, 95th) of simulated data
3. Overlay observed data
4. Check if observed quantiles fall within simulated CI

## Files

| File | Description |
|------|-------------|
| julia.jl | Julia implementation |
| python.py | Python implementation |
