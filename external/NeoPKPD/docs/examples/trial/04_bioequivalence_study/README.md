# Bioequivalence Study

Complete bioequivalence study with 2x2 crossover design and TOST analysis.

## Design

- **Type**: 2x2 crossover bioequivalence
- **Treatments**: Reference vs Test formulation
- **N**: 24 subjects (12 per sequence)
- **Washout**: 7 days
- **BE limits**: 80.00% - 125.00%
- **Parameters**: Cmax, AUC0-inf

## Regulatory Guidance

Following FDA guidance for BE studies:
- 90% CI for GMR must be within [0.80, 1.25]
- Log-transformed data
- TOST (Two One-Sided Tests) procedure

## Features Demonstrated

- BE study design
- NCA parameter calculation
- 90% CI for geometric mean ratio
- BE assessment per FDA guidance

## Files

| File | Description |
|------|-------------|
| [python.py](python.py) | Python implementation |
| [julia.jl](julia.jl) | Julia implementation |
| [cli.json](cli.json) | CLI specification |

## Expected Output

```
Bioequivalence Study
====================
Design: 2x2 crossover
N: 24 subjects
BE limits: [0.80, 1.25]

Results:
--------
Parameter   GMR      90% CI           BE?
Cmax        0.9845   [0.8912, 1.0875] Yes
AUC0-inf    1.0123   [0.9456, 1.0834] Yes

Conclusion: Bioequivalence demonstrated
```
