# Sigmoid Emax (Hill) PD Model

PK/PD model with sigmoidal concentration-effect relationship.

## Model Schematic

```
┌─────────────┐
│     PK      │──── Conc ────┐
│ OneCompIVBolus             │
└─────────────┘              ▼
                      ┌───────────────────────────┐
                      │   Sigmoid Effect (Hill)   │
                      │                           │
                      │ E = E0 + Emax × C^γ       │
                      │       ─────────────       │
                      │       EC50^γ + C^γ        │
                      └───────────────────────────┘
```

## PD Equation

```
Effect = E0 + (Emax × C^γ) / (EC50^γ + C^γ)
```

Where:
- E0 = baseline effect
- Emax = maximum effect
- EC50 = concentration at 50% effect
- γ (gamma) = Hill coefficient (sigmoidicity factor)

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| E0 | Baseline effect | 0-100 | varies |
| Emax | Maximum effect | 0-100 | varies |
| EC50 | Potency | 0.1-100 | mg/L |
| gamma | Hill coefficient | 0.5-5 | - |

## Hill Coefficient Interpretation

| γ Value | Curve Shape | Interpretation |
|---------|-------------|----------------|
| < 1 | Shallow | Negative cooperativity |
| = 1 | Hyperbolic | No cooperativity (= Emax) |
| > 1 | Steep | Positive cooperativity |
| >> 1 | Switch-like | Threshold effect |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| E0 | 0 | No baseline |
| Emax | 100 | 100% max effect |
| EC50 | 1.0 mg/L | Moderate potency |
| gamma | 2.0 | Positive cooperativity |

## Key Features

- Steeper dose-response than simple Emax
- Common for drugs with cooperative binding
- Can model all-or-none responses (high γ)

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
