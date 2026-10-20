# Two-Compartment IV Bolus Model

Two-compartment model with rapid distribution to peripheral tissues.

## Model Schematic

```
IV Bolus
    │
    ▼
┌─────────────┐     Q      ┌─────────────┐
│   Central   │ ◄────────▶ │ Peripheral  │
│     V1      │            │     V2      │
└──────┬──────┘            └─────────────┘
       │
       │ CL
       ▼
   Elimination
```

## Differential Equations

```
dA_central/dt = Q/V2 × A_peripheral - Q/V1 × A_central - CL/V1 × A_central
dA_peripheral/dt = Q/V1 × A_central - Q/V2 × A_peripheral
```

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| CL | Clearance | 1-50 | L/h |
| V1 | Central volume | 5-100 | L |
| Q | Inter-compartmental clearance | 1-50 | L/h |
| V2 | Peripheral volume | 10-500 | L |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| CL | 5.0 L/h | Moderate clearance |
| V1 | 10.0 L | Plasma/blood volume |
| Q | 10.0 L/h | Rapid distribution |
| V2 | 40.0 L | Tissue distribution |
| Dose | 100 mg | IV bolus |

## Key Features

- Bi-exponential decline (α and β phases)
- α phase: distribution + elimination
- β phase: terminal elimination
- V1 + V2 = Vss (volume at steady-state)

## Derived Parameters

| Parameter | Formula | Description |
|-----------|---------|-------------|
| Vss | V1 + V2 | Steady-state volume |
| α, β | Macro rate constants | Eigenvalues |
| t½α | ln(2)/α | Distribution half-life |
| t½β | ln(2)/β | Terminal half-life |

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
