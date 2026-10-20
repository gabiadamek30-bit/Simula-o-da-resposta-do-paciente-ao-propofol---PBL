# Quickstart Examples

Get started with NeoPKPD in 5 minutes. These examples demonstrate a basic one-compartment IV bolus simulation.

## Files

| File | Description |
|------|-------------|
| `julia_first_simulation.jl` | Julia quickstart |
| `python_first_simulation.py` | Python quickstart |
| `cli_first_simulation.sh` | CLI quickstart |
| `spec.json` | CLI specification file |

## What You'll Learn

1. Create a basic PK model specification
2. Configure simulation parameters
3. Run the simulation
4. Extract key metrics (Cmax, Tmax, AUC)

## The Model

**One-Compartment IV Bolus**

```
Dose (100 mg IV)
     │
     ▼
┌─────────────┐
│   Central   │
│   V = 50 L  │──── CL = 5 L/h ────▶ Elimination
└─────────────┘
```

**Parameters:**
- CL = 5 L/h (clearance)
- V = 50 L (volume of distribution)
- Dose = 100 mg at t=0

**Expected Results:**
- Cmax = 2.0 mg/L (at t=0)
- AUC₀₋₂₄ ≈ 18.1 mg·h/L
- t½ ≈ 6.9 h

## Running the Examples

### Julia

```bash
julia --project=packages/core docs/examples/quickstart/julia_first_simulation.jl
```

### Python

```bash
source packages/python/.venv/bin/activate
python docs/examples/quickstart/python_first_simulation.py
```

### CLI

```bash
./bin/neopkpd simulate --spec docs/examples/quickstart/spec.json --out result.json
./bin/neopkpd metrics --artifact result.json --metrics cmax,tmax,auc
```

## Next Steps

After completing the quickstart:

1. **Explore Models**: See [models](../models/README.md) for all available PK/PD models
2. **Add Variability**: See [population](../population/README.md) for IIV and covariates
3. **Fit to Data**: See [estimation](../estimation/README.md) for parameter estimation
4. **Real Analysis**: See [use_cases](../use_cases/fih_dose_exploration/README.md) for complete workflows
