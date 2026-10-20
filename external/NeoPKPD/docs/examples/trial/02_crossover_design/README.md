# 2x2 Crossover Design

Two-period, two-sequence crossover study for within-subject comparison.

## Design

- **Type**: 2x2 crossover (AB, BA)
- **Treatments**: Reference vs Test formulation
- **Washout**: 7 days between periods
- **N per sequence**: 12 subjects
- **Endpoints**: PK exposure (Cmax, AUC)

## Sequences

| Sequence | Period 1 | Washout | Period 2 |
|----------|----------|---------|----------|
| AB | Reference | 7 days | Test |
| BA | Test | 7 days | Reference |

## Features Demonstrated

- 2x2 crossover design creation
- Washout period specification
- Period-specific dosing
- Within-subject analysis

## Files

| File | Description |
|------|-------------|
| [python.py](python.py) | Python implementation |
| [julia.jl](julia.jl) | Julia implementation |
| [cli.json](cli.json) | CLI specification |

## Expected Output

```
2x2 Crossover Study
===================
Design: 2 periods, 2 sequences
Washout: 7 days

Sequence AB: 12 subjects
Sequence BA: 12 subjects

Completion rate: 95.8%
```
