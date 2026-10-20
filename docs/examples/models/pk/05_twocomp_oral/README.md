# Two-Compartment Oral Model

Two-compartment model with first-order oral absorption.

## Model Schematic

```
Oral Dose
    │
    ▼
┌─────────┐     Ka      ┌─────────────┐     Q      ┌─────────────┐
│   Gut   │ ─────────▶  │   Central   │ ◄────────▶ │ Peripheral  │
│  (A_d)  │             │     V1      │            │     V2      │
└─────────┘             └──────┬──────┘            └─────────────┘
                               │
                               │ CL
                               ▼
                          Elimination
```

## Differential Equations

```
dA_depot/dt = -Ka × A_depot
dA_central/dt = Ka × A_depot + Q/V2 × A_peripheral - Q/V1 × A_central - CL/V1 × A_central
dA_peripheral/dt = Q/V1 × A_central - Q/V2 × A_peripheral
```

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| Ka | Absorption rate constant | 0.1-5.0 | 1/h |
| CL | Clearance | 1-50 | L/h |
| V1 | Central volume | 5-100 | L |
| Q | Inter-compartmental clearance | 1-50 | L/h |
| V2 | Peripheral volume | 10-500 | L |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Ka | 1.5 1/h | Moderate absorption |
| CL | 5.0 L/h | Moderate clearance |
| V1 | 10.0 L | Central compartment |
| Q | 10.0 L/h | Rapid distribution |
| V2 | 40.0 L | Tissue distribution |
| Dose | 100 mg | Oral dose |

## Key Features

- Tri-exponential profile (absorption + α + β phases)
- Distribution occurs during absorption
- Common for lipophilic drugs with significant tissue binding

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
