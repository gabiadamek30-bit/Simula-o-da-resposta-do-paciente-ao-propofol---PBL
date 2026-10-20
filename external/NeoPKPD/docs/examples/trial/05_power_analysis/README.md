# Power Analysis and Sample Size

Statistical power analysis and sample size estimation for clinical trials.

## Features

- Power calculation for given sample size
- Sample size estimation for target power
- Effect size interpretation (Cohen's d)
- Multiple comparison considerations

## Cohen's d Reference

| Effect Size | Cohen's d | Interpretation |
|-------------|-----------|----------------|
| Small | 0.2 | Subtle difference |
| Medium | 0.5 | Moderate difference |
| Large | 0.8 | Obvious difference |

## Sample Size Guide (80% power, alpha=0.05)

| Effect Size | N per arm |
|-------------|-----------|
| Small (0.2) | ~394 |
| Medium (0.5) | ~64 |
| Large (0.8) | ~26 |

## Features Demonstrated

- Analytical power calculation
- Sample size estimation
- Dropout adjustment
- Power curves

## Files

| File | Description |
|------|-------------|
| [python.py](python.py) | Python implementation |
| [julia.jl](julia.jl) | Julia implementation |

## Expected Output

```
Power Analysis
==============

Given: n=50 per arm, effect=0.5, alpha=0.05
Calculated power: 69.7%

Sample Size Estimation
======================
Target: 80% power, effect=0.5, alpha=0.05
Required N per arm: 64
Total N: 128
Achieved power: 80.1%

With 15% dropout adjustment:
Required N per arm: 76
Total N: 152
```
