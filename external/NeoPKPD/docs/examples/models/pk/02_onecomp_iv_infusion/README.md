# One-Compartment IV Infusion Model

One-compartment model with constant-rate intravenous infusion.

## Model Schematic

```
Infusion (Rate = Dose/Duration)
         │
         ▼
    ┌─────────────┐
    │   Central   │
    │      V      │──── CL ────▶ Elimination
    └─────────────┘
```

## Differential Equation

**During infusion (0 ≤ t ≤ T_inf):**
```
dA_central/dt = Rate - CL/V × A_central
```

**After infusion (t > T_inf):**
```
dA_central/dt = -CL/V × A_central
```

Where `Rate = Dose / T_inf`

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| CL | Clearance | 1-50 | L/h |
| V | Volume of distribution | 10-500 | L |
| duration | Infusion duration | 0.5-4 | h |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| CL | 5.0 L/h | Moderate clearance |
| V | 50.0 L | Total body water |
| Dose | 100 mg | Total infusion dose |
| Duration | 1.0 h | 1-hour infusion |

## Derived Parameters

| Parameter | Formula | Example Value |
|-----------|---------|---------------|
| Rate | Dose/Duration | 100 mg/h |
| Css | Rate/CL | 20 mg/L (if continued) |
| t½ | ln(2)×V/CL | 6.93 h |

## Key Features

- Cmax occurs at end of infusion (t = T_inf)
- Lower Cmax compared to bolus for same dose
- Useful for drugs with narrow therapeutic window
- Rate-limiting can reduce infusion reactions

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
