# Three-Compartment IV Bolus Model

Three-compartment model for drugs with complex tissue distribution.

## Model Schematic

```
                  IV Bolus
                      │
                      ▼
┌─────────────┐  Q2  ┌─────────────┐  Q3  ┌─────────────┐
│ Peripheral  │◄────▶│   Central   │◄────▶│ Peripheral  │
│   Shallow   │      │     V1      │      │    Deep     │
│     V2      │      └──────┬──────┘      │     V3      │
└─────────────┘             │             └─────────────┘
                            │ CL
                            ▼
                       Elimination
```

## Differential Equations

```
dA_central/dt = Q2/V2 × A_shallow - Q2/V1 × A_central
              + Q3/V3 × A_deep - Q3/V1 × A_central
              - CL/V1 × A_central

dA_shallow/dt = Q2/V1 × A_central - Q2/V2 × A_shallow

dA_deep/dt = Q3/V1 × A_central - Q3/V3 × A_deep
```

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| CL | Clearance | 1-50 | L/h |
| V1 | Central volume | 5-50 | L |
| Q2 | Shallow distribution clearance | 5-100 | L/h |
| V2 | Shallow peripheral volume | 10-200 | L |
| Q3 | Deep distribution clearance | 0.1-10 | L/h |
| V3 | Deep peripheral volume | 50-1000 | L |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| CL | 5.0 L/h | Moderate clearance |
| V1 | 10.0 L | Plasma volume |
| Q2 | 20.0 L/h | Rapid shallow distribution |
| V2 | 20.0 L | Well-perfused tissues |
| Q3 | 2.0 L/h | Slow deep distribution |
| V3 | 100.0 L | Poorly perfused tissues/fat |
| Dose | 100 mg | IV bolus |

## Key Features

- Tri-exponential decline (α, β, γ phases)
- α phase: rapid distribution
- β phase: slow distribution
- γ phase: terminal elimination
- Vss = V1 + V2 + V3

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
