# One-Compartment IV Bolus Model

The simplest pharmacokinetic model for intravenous administration.

## Model Schematic

```
Dose (IV Bolus)
     │
     ▼
┌─────────────┐
│   Central   │
│      V      │──── CL ────▶ Elimination
└─────────────┘
```

## Differential Equation

```
dA_central/dt = -CL/V × A_central
```

With initial condition: `A_central(0) = Dose`

## Observation

```
Conc = A_central / V
```

## Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| CL | Clearance | 1-50 | L/h |
| V | Volume of distribution | 10-500 | L |

## Example Values

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| CL | 5.0 L/h | Moderate hepatic clearance |
| V | 50.0 L | Roughly total body water |
| Dose | 100 mg | Single 100mg IV bolus |

## Analytical Solution

```
C(t) = (Dose/V) × exp(-CL/V × t)
```

## Derived Parameters

| Parameter | Formula | Example Value |
|-----------|---------|---------------|
| k | CL/V | 0.1 h⁻¹ |
| t½ | ln(2)/k | 6.93 h |
| Cmax | Dose/V | 2.0 mg/L |
| AUC₀₋∞ | Dose/CL | 20 mg·h/L |

## Files

| File | Description |
|------|-------------|
| `julia.jl` | Julia implementation |
| `python.py` | Python implementation |
| `cli.json` | CLI specification |

## Running

```bash
# Julia
julia --project=packages/core julia.jl

# Python
python python.py

# CLI
./bin/neopkpd simulate --spec cli.json
```
