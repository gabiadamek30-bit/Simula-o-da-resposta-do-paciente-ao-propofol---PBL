# Michaelis-Menten Elimination Model

One-compartment model with saturable (nonlinear) elimination kinetics.

## Model Schematic

```
IV Bolus
    │
    ▼
┌─────────────┐
│   Central   │
│      V      │──── Vmax/(Km + C) ────▶ Saturable Elimination
└─────────────┘
```

## Differential Equation

```
dA_central/dt = -Vmax × A_central / (Km × V + A_central)
```

Or in terms of concentration:

```
dC/dt = -Vmax × C / (Km + C) / V
```

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| Vmax | Maximum elimination rate | 10-1000 | mg/h |
| Km | Michaelis constant | 1-100 | mg/L |
| V | Volume of distribution | 10-500 | L |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Vmax | 50.0 mg/h | Maximum enzyme capacity |
| Km | 10.0 mg/L | Half-saturation concentration |
| V | 50.0 L | Total body water |
| Dose | 500 mg | High dose to show saturation |

## Key Features

- **At low C (C << Km):** First-order kinetics, CL ≈ Vmax/Km
- **At high C (C >> Km):** Zero-order kinetics, rate ≈ Vmax
- Nonlinear PK: AUC not proportional to dose
- Examples: phenytoin, ethanol, aspirin (high dose)

## Derived Parameters

| Parameter | Formula | Description |
|-----------|---------|-------------|
| CLint | Vmax/Km | Intrinsic clearance (at low C) |
| t½ (low C) | 0.693 × V × Km / Vmax | Half-life at low concentrations |

## Use Cases

- Drugs with capacity-limited metabolism
- High-dose regimens
- Therapeutic drug monitoring
- Phenytoin dosing optimization

## Clinical Implications

1. **Dose escalation:** Small dose increases can cause large concentration changes
2. **Steady-state:** Takes longer to reach at higher doses
3. **Drug interactions:** Enzyme inhibitors have greater effect when near saturation

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
