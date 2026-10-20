# Model Examples

Complete examples for all NeoPKPD PK and PK/PD models.

## PK Models

| Model | Description | Directory |
|-------|-------------|-----------|
| OneCompIVBolus | One-compartment IV bolus | [01_onecomp_iv_bolus](pk/01_onecomp_iv_bolus/README.md) |
| OneCompIVInfusion | One-compartment IV infusion | [02_onecomp_iv_infusion](pk/02_onecomp_iv_infusion/README.md) |
| OneCompOralFirstOrder | One-compartment oral (first-order absorption) | [03_onecomp_oral](pk/03_onecomp_oral/README.md) |
| TwoCompIVBolus | Two-compartment IV bolus | [04_twocomp_iv](pk/04_twocomp_iv/README.md) |
| TwoCompOral | Two-compartment oral | [05_twocomp_oral](pk/05_twocomp_oral/README.md) |
| ThreeCompIVBolus | Three-compartment IV bolus | [06_threecomp_iv](pk/06_threecomp_iv/README.md) |
| TransitAbsorption | Transit compartment absorption | [07_transit_absorption](pk/07_transit_absorption/README.md) |
| MichaelisMentenElimination | Saturable (nonlinear) elimination | [08_michaelis_menten](pk/08_michaelis_menten/README.md) |

## PK/PD Models

| Model | Description | Directory |
|-------|-------------|-----------|
| DirectEmax | Direct Emax response | [01_direct_emax](pkpd/01_direct_emax/README.md) |
| SigmoidEmax | Sigmoidal Emax response | [02_sigmoid_emax](pkpd/02_sigmoid_emax/README.md) |
| BiophaseEquilibration | Effect compartment model | [03_biophase_equilibration](pkpd/03_biophase_equilibration/README.md) |
| IndirectResponse | Indirect response (turnover) | [04_indirect_response](pkpd/04_indirect_response/README.md) |

## File Structure

Each model directory contains:

```
01_model_name/
├── README.md       # Model equations, parameters, typical values
├── julia.jl        # Julia implementation
├── python.py       # Python implementation
└── cli.json        # CLI specification
```

## Running Examples

### All Models (Julia)

```bash
julia --project=packages/core docs/examples/models/run_all.jl
```

### Single Model

```bash
# Julia
julia --project=packages/core docs/examples/models/pk/01_onecomp_iv_bolus/julia.jl

# Python
python docs/examples/models/pk/01_onecomp_iv_bolus/python.py

# CLI
./bin/neopkpd simulate --spec docs/examples/models/pk/01_onecomp_iv_bolus/cli.json
```

## Model Selection Guide

### By Administration Route

| Route | Single Compartment | Multi-Compartment |
|-------|-------------------|-------------------|
| IV Bolus | OneCompIVBolus | TwoCompIVBolus, ThreeCompIVBolus |
| IV Infusion | OneCompIVInfusion | TwoCompIVInfusion |
| Oral | OneCompOralFirstOrder | TwoCompOral |
| Oral (complex) | TransitAbsorption | - |

### By Elimination Type

| Elimination | Model |
|-------------|-------|
| Linear (first-order) | OneCompIVBolus, TwoCompIVBolus, etc. |
| Nonlinear (saturable) | MichaelisMentenElimination |

### By PD Response

| Response Type | Model | Use Case |
|---------------|-------|----------|
| Direct, hyperbolic | DirectEmax | Immediate response, receptor binding |
| Direct, sigmoidal | SigmoidEmax | Steep dose-response, threshold effects |
| Delayed | BiophaseEquilibration | Effect compartment, hysteresis |
| Indirect | IndirectResponse | Enzyme/receptor turnover, biomarkers |

## Parameter Conventions

### PK Parameters

| Parameter | Description | Typical Units |
|-----------|-------------|---------------|
| CL | Clearance | L/h |
| V, V1 | Central volume | L |
| Q | Inter-compartmental clearance | L/h |
| V2, V3 | Peripheral volumes | L |
| Ka | Absorption rate constant | 1/h |
| Ktr | Transit rate constant | 1/h |
| Vmax | Maximum elimination rate | mg/h |
| Km | Michaelis constant | mg/L |

### PD Parameters

| Parameter | Description | Typical Units |
|-----------|-------------|---------------|
| E0 | Baseline effect | varies |
| Emax | Maximum effect | varies |
| EC50 | Concentration at 50% Emax | mg/L |
| gamma | Hill coefficient | dimensionless |
| Ke0 | Effect compartment rate constant | 1/h |
| Kin | Zero-order production rate | units/h |
| Kout | First-order elimination rate | 1/h |
| IC50/SC50 | Concentration for 50% inhibition/stimulation | mg/L |
