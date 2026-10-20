# Biophase Equilibration (Effect Compartment) Model

PK/PD model with delayed effect due to distribution to an effect site.

## Model Schematic

```
IV Bolus
    │
    ▼
┌─────────────┐     Ke0      ┌─────────────┐
│   Central   │ ───────────▶ │   Effect    │
│      V      │              │ Compartment │
└──────┬──────┘              │    (Ce)     │
       │                     └──────┬──────┘
       │ CL                         │
       ▼                            ▼
   Elimination              ┌─────────────────┐
                            │ E = E0 + Emax×Ce│
                            │     ───────────  │
                            │     EC50 + Ce   │
                            └─────────────────┘
```

## Differential Equations

**PK:**
```
dA_central/dt = -CL/V × A_central
```

**Effect Compartment:**
```
dCe/dt = Ke0 × (Cp - Ce)
```

**PD:**
```
Effect = E0 + Emax × Ce / (EC50 + Ce)
```

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| Ke0 | Effect site equilibration rate | 0.1-5 | 1/h |
| E0 | Baseline effect | 0-100 | varies |
| Emax | Maximum effect | 0-100 | varies |
| EC50 | Effect site EC50 | 0.1-100 | mg/L |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| CL | 5.0 L/h | Moderate clearance |
| V | 50.0 L | Central volume |
| Ke0 | 0.5 1/h | 1.4h equilibration half-life |
| E0 | 0 | No baseline |
| Emax | 100 | Maximum response |
| EC50 | 1.0 mg/L | Effect site potency |

## Key Features

- **Counter-clockwise hysteresis**: Effect lags behind concentration
- **t½ke0 = ln(2)/Ke0**: Time for effect compartment to equilibrate
- Peak effect occurs after peak concentration
- Effect persists longer than plasma concentration

## Use Cases

- CNS-active drugs (BBB penetration delay)
- Anesthetics (propofol, remifentanil)
- Sedatives/hypnotics
- Analgesics

## Hysteresis Interpretation

When plotting Effect vs Concentration:
- **Counter-clockwise loop**: Effect delayed (biophase model)
- **Clockwise loop**: Tolerance development

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |
