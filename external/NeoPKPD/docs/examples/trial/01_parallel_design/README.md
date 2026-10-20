# Parallel Group Design

Two-arm parallel group Phase II efficacy study simulation.

## Design

- **Type**: Parallel, 1:1 randomization
- **Arms**: Placebo vs Active (100 mg QD)
- **Duration**: 28 days
- **N per arm**: 50 subjects
- **Endpoints**: PK exposure

## Features Demonstrated

- Parallel design creation
- Treatment arm specification
- QD dosing regimen
- Dropout modeling
- Compliance patterns
- Trial simulation

## Files

| File | Description |
|------|-------------|
| [python.py](python.py) | Python implementation |
| [julia.jl](julia.jl) | Julia implementation |
| [cli.json](cli.json) | CLI specification |

## Expected Output

```
Phase 2 Parallel Study
======================
Design: 2-arm parallel (1:1)
Duration: 28 days

Arm Results:
  Placebo: 48/50 completed (96.0%)
  Active:  47/50 completed (94.0%)

Overall completion: 95.0%
Overall compliance: 89.5%
```
