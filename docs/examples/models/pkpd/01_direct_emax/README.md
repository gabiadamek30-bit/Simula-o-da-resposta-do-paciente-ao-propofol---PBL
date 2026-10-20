# Direct Emax PD Model

PK/PD model with direct, hyperbolic concentration-effect relationship.

## Model Schematic

```
┌─────────────┐
│     PK      │──── Conc ────┐
│ OneCompIVBolus             │
└─────────────┘              ▼
                      ┌─────────────────────┐
                      │   Direct Effect     │
                      │                     │
                      │ E = E0 + Emax × C   │
                      │       ──────────    │
                      │        EC50 + C     │
                      └─────────────────────┘
```

## PD Equation

```
Effect = E0 + (Emax × C) / (EC50 + C)
```

Where:
- E0 = baseline effect
- Emax = maximum effect above baseline
- EC50 = concentration producing 50% of Emax
- C = drug concentration

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| E0 | Baseline effect | 0-100 | varies |
| Emax | Maximum effect | 0-100 | varies |
| EC50 | Potency (C for 50% Emax) | 0.1-100 | mg/L |

### PK Parameters
| Parameter | Description | Units |
|-----------|-------------|-------|
| CL | Clearance | L/h |
| V | Volume | L |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| E0 | 0 | No baseline effect |
| Emax | 100 | Maximum response = 100 units |
| EC50 | 1.0 mg/L | Moderate potency |
| CL | 5.0 L/h | Moderate clearance |
| V | 50.0 L | Total body water |

## Key Features

- Immediate response (no hysteresis)
- Hyperbolic saturation curve
- Effect at C = EC50 is E0 + Emax/2
- Effect at C >> EC50 approaches E0 + Emax

## Use Cases

- Receptor binding studies
- Enzyme inhibition
- Blood pressure lowering
- Biomarkers with rapid equilibration

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
