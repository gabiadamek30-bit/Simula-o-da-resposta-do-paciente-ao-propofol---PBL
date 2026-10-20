# Transit Compartment Absorption Model

Model with transit compartments for delayed/complex oral absorption.

## Model Schematic

```
Oral Dose
    │
    ▼
┌─────────┐  Ktr   ┌─────────┐  Ktr   ┌─────────┐  Ktr   ┌─────────┐  Ka   ┌─────────────┐
│Transit 1│ ────▶ │Transit 2│ ────▶ │Transit 3│ ────▶ │  ... n  │ ────▶│   Central   │
└─────────┘       └─────────┘       └─────────┘       └─────────┘      │      V      │
                                                                        └──────┬──────┘
                                                                               │ CL
                                                                               ▼
                                                                          Elimination
```

## Key Equations

The transit model uses n identical transit compartments with rate constant Ktr:

```
dA₁/dt = -Ktr × A₁
dAᵢ/dt = Ktr × Aᵢ₋₁ - Ktr × Aᵢ   (for i = 2 to n)
dA_central/dt = Ktr × Aₙ - CL/V × A_central
```

**Mean Transit Time (MTT):**
```
MTT = (n + 1) / Ktr
```

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| Ktr | Transit rate constant | 0.5-5.0 | 1/h |
| n | Number of transit compartments | 1-10 | - |
| CL | Clearance | 1-50 | L/h |
| V | Volume of distribution | 10-500 | L |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Ktr | 2.0 1/h | MTT ≈ 2h with n=3 |
| n | 3 | 3 transit compartments |
| CL | 5.0 L/h | Moderate clearance |
| V | 50.0 L | Total body water |
| Dose | 100 mg | Oral dose |

## Key Features

- Delayed Tmax compared to first-order absorption
- Sharper peak than first-order absorption
- More physiologically realistic for enteric-coated tablets
- Flexible: adjusting n and Ktr can model various profiles

## Use Cases

- Enteric-coated formulations
- Extended-release tablets
- Drugs with gastric emptying delay
- Complex absorption profiles

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
