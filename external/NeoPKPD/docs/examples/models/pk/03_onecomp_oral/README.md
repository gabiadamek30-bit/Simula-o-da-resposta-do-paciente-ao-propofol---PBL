# One-Compartment Oral First-Order Absorption Model

One-compartment model with first-order absorption from an oral dose.

## Model Schematic

```
Oral Dose
    │
    ▼
┌─────────┐     Ka      ┌─────────────┐
│   Gut   │ ─────────▶  │   Central   │
│  (A_d)  │             │      V      │──── CL ────▶ Elimination
└─────────┘             └─────────────┘
```

## Differential Equations

```
dA_depot/dt = -Ka × A_depot
dA_central/dt = Ka × A_depot - CL/V × A_central
```

With initial conditions: `A_depot(0) = Dose × F`, `A_central(0) = 0`

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| Ka | Absorption rate constant | 0.1-5.0 | 1/h |
| CL | Clearance | 1-50 | L/h |
| V | Volume of distribution | 10-500 | L |
| F | Bioavailability | 0.1-1.0 | - |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Ka | 1.5 1/h | Moderate absorption |
| CL | 5.0 L/h | Moderate clearance |
| V | 50.0 L | Total body water |
| Dose | 100 mg | Oral tablet |
| F | 1.0 | Assumed 100% |

## Key Features

- Absorption phase followed by elimination phase
- Tmax occurs when absorption rate = elimination rate
- "Flip-flop" kinetics possible when Ka < k

## Derived Parameters

| Parameter | Formula | Example Value |
|-----------|---------|---------------|
| k | CL/V | 0.1 h⁻¹ |
| Tmax | ln(Ka/k)/(Ka-k) | 1.79 h |
| Cmax | (F×Dose/V)×(Ka/(Ka-k))×[exp(-k×Tmax)-exp(-Ka×Tmax)] | 1.43 mg/L |
| t½ | ln(2)/k | 6.93 h |

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
