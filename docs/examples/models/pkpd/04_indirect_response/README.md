# Indirect Response (Turnover) Model

PK/PD model for drugs that affect the synthesis or degradation of a response variable.

## Model Schematic

```
┌─────────────┐
│     PK      │──── Conc ────┐
│ OneCompIVBolus             │
└─────────────┘              │
                             ▼
                      ┌─────────────────────────┐
                      │ Inhibition/Stimulation  │
                      │  of Kin or Kout         │
                      └───────────┬─────────────┘
                                  │
                                  ▼
                      ┌─────────────────────────┐
                      │   Response Variable     │
                      │   dR/dt = Kin - Kout×R  │
                      └─────────────────────────┘
```

## Four Types of Indirect Response Models

| Type | Drug Effect | Example |
|------|-------------|---------|
| IRM-I | Inhibition of Kin | Corticosteroids on cortisol |
| IRM-II | Inhibition of Kout | Warfarin on clotting factors |
| IRM-III | Stimulation of Kin | EPO on RBC production |
| IRM-IV | Stimulation of Kout | Diuretics on sodium excretion |

## Differential Equations

**Response:**
```
dR/dt = Kin × (1 - I(C)) - Kout × R    (IRM-I: inhibit production)
dR/dt = Kin - Kout × (1 - I(C)) × R    (IRM-II: inhibit loss)
dR/dt = Kin × (1 + S(C)) - Kout × R    (IRM-III: stimulate production)
dR/dt = Kin - Kout × (1 + S(C)) × R    (IRM-IV: stimulate loss)
```

Where:
- `I(C) = Imax × C / (IC50 + C)` (inhibition function)
- `S(C) = Smax × C / (SC50 + C)` (stimulation function)

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| Kin | Zero-order production rate | 0.1-100 | units/h |
| Kout | First-order loss rate | 0.01-1 | 1/h |
| Imax/Smax | Maximum inhibition/stimulation | 0.1-1 | - |
| IC50/SC50 | Potency for effect | 0.1-100 | mg/L |

## Example Values (IRM-I)

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| CL | 5.0 L/h | Drug clearance |
| V | 50.0 L | Drug volume |
| Kin | 10 units/h | Production rate |
| Kout | 0.1 1/h | Loss rate (t½ = 6.9h) |
| Imax | 0.9 | 90% max inhibition |
| IC50 | 1.0 mg/L | Potency |

## Key Features

- Baseline: R₀ = Kin/Kout
- Response half-life: t½_R = ln(2)/Kout
- Delayed response (not immediate)
- Rebound possible after drug washout

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
