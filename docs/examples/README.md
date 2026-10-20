# NeoPKPD Examples

Comprehensive examples covering all NeoPKPD features across Julia, Python, and CLI interfaces.

## Quick Navigation

| Category | Description | Languages |
|----------|-------------|-----------|
| [Quickstart](quickstart/README.md) | Get started in 5 minutes | Julia, Python, CLI |
| [Models](models/README.md) | All 12 PK/PD models | Julia, Python, CLI |
| [Population](population/README.md) | IIV, IOV, covariates | Julia, Python, CLI |
| [Estimation](estimation/README.md) | FOCE-I, SAEM, diagnostics | Julia, Python |
| [NCA](nca/README.md) | Non-compartmental analysis | Julia, Python |
| [VPC](vpc/README.md) | Visual Predictive Checks | Julia, Python, CLI |
| [Trial](trial/README.md) | Clinical trial simulation | Julia, Python, CLI |
| [Import](import/README.md) | NONMEM/Monolix import | Julia, Python, CLI |
| [Data](data/README.md) | CDISC data import | Julia, Python |
| [Visualization](visualization/README.md) | Plotting and figures | Python |
| [Sensitivity](sensitivity/README.md) | Parameter sensitivity | Julia, Python, CLI |
| [Reproducibility](reproducibility/README.md) | Artifacts and replay | Julia, Python, CLI |

## End-to-End Use Cases

Complete workflows from data to analysis:

| Use Case | Description |
|----------|-------------|
| [FIH Dose Exploration](use_cases/fih_dose_exploration/README.md) | First-in-human dose selection |
| [PKPD Biomarker](use_cases/pkpd_biomarker_turnover/README.md) | Biomarker turnover modeling |
| [NONMEM Migration](use_cases/nonmem_migration/README.md) | Migrate NONMEM models to NeoPKPD |
| [Theophylline Analysis](use_cases/real_world_theophylline/README.md) | Real-world PK analysis |
| [Bioequivalence Study](use_cases/bioequivalence_study/README.md) | Complete BE workflow |
| [Population PKPD](use_cases/population_pkpd_analysis/README.md) | Full population analysis |

## Real-World Validation

Validation against published datasets:

| Study | Dataset | Features |
|-------|---------|----------|
| [Theophylline SD](real_world_validation/studies/theophylline_theo_sd/README.md) | Single dose PK | NCA, estimation |
| [Theophylline MD](real_world_validation/studies/theophylline_theo_md/README.md) | Multiple dose PK | Population, VPC |
| [Warfarin PKPD](real_world_validation/studies/warfarin_pkpd/README.md) | PK/PD model | Indirect response |

---

## Quickstart

### Julia (5 minutes)

```julia
using NeoPKPD

# Create a one-compartment IV bolus model
model = create_model_spec("OneCompIVBolus",
    params = Dict("CL" => 5.0, "V" => 50.0),
    doses = [DoseEvent(time=0.0, amount=100.0)]
)

# Simulate
grid = SimulationGrid(t0=0.0, t1=24.0, saveat=0:0.5:24)
result = simulate(model, grid)

# Extract metrics
println("Cmax: ", maximum(result.observations["conc"]))
```

### Python (5 minutes)

```python
import neopkpd

# Initialize Julia backend (once per session)
neopkpd.init_julia()

# Simulate one-compartment IV bolus
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0,
    v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
)

# Extract metrics
print(f"Cmax: {max(result['observations']['conc'])}")
```

### CLI (5 minutes)

```bash
# Run simulation
./bin/neopkpd simulate --spec quickstart/spec.json --out result.json

# Compute metrics
./bin/neopkpd metrics --artifact result.json --metrics cmax,tmax,auc
```

---

## Models Reference

### PK Models

| Model | Compartments | Administration | Key Parameters |
|-------|--------------|----------------|----------------|
| OneCompIVBolus | 1 | IV Bolus | CL, V |
| OneCompIVInfusion | 1 | IV Infusion | CL, V, duration |
| OneCompOralFirstOrder | 1 | Oral | Ka, CL, V |
| TwoCompIVBolus | 2 | IV Bolus | CL, V1, Q, V2 |
| TwoCompOral | 2 | Oral | Ka, CL, V1, Q, V2 |
| ThreeCompIVBolus | 3 | IV Bolus | CL, V1, Q2, V2, Q3, V3 |
| TransitAbsorption | 1 + transit | Oral | Ktr, n, CL, V |
| MichaelisMentenElimination | 1 | IV | Vmax, Km, V |

### PKPD Models

| Model | Response Type | Key Parameters |
|-------|---------------|----------------|
| DirectEmax | Direct effect | Emax, EC50, E0 |
| SigmoidEmax | Sigmoidal | Emax, EC50, gamma, E0 |
| BiophaseEquilibration | Effect compartment | Emax, EC50, Ke0 |
| IndirectResponse | Turnover | Kin, Kout, Imax/Smax, IC50/SC50 |

---

## Population Modeling

### IIV (Inter-Individual Variability)

```julia
# Log-normal IIV on CL and V
pop_spec = create_population_spec(
    base_model = model,
    iiv = LogNormalIIV(
        omegas = Dict("CL" => 0.3, "V" => 0.2),
        seed = 12345
    ),
    n = 100
)
result = simulate_population(pop_spec, grid)
```

### Covariates

```julia
# Weight-based allometric scaling
covariate_model = create_covariate_model(
    effects = [
        CovariateEffect(:WT, :CL, :power, exponent=0.75, reference=70.0),
        CovariateEffect(:WT, :V, :power, exponent=1.0, reference=70.0)
    ]
)
```

---

## Feature Examples

### Parameter Estimation (FOCE-I)

```julia
# Estimate parameters from observed data
config = EstimationConfig(
    method = FOCEI(),
    theta_init = [5.0, 50.0, 1.5],
    theta_lower = [0.1, 1.0, 0.1],
    theta_upper = [100.0, 500.0, 10.0],
    omega_init = [0.09, 0.04, 0.16]
)
result = estimate(observed_data, model_spec, config)
```

### Visual Predictive Check (VPC)

```julia
# Generate VPC
vpc_result = compute_vpc(
    observed_data,
    pop_spec,
    n_simulations = 500,
    quantiles = [0.05, 0.5, 0.95]
)
```

### Non-Compartmental Analysis (NCA)

```julia
# Compute NCA metrics
nca_result = compute_nca(
    times = observed_times,
    concentrations = observed_conc,
    dose = 100.0,
    route = :oral
)
# Returns: Cmax, Tmax, AUC_0_t, AUC_0_inf, t_half, CL_F, Vz_F
```

### Trial Simulation

```julia
# Simulate bioequivalence study
trial_spec = TrialSpec(
    design = CrossoverDesign(periods=2, sequences=2, washout=7),
    arm_specs = [
        ArmSpec(treatment="Reference", n=12),
        ArmSpec(treatment="Test", n=12)
    ],
    endpoints = [BioequivalenceEndpoint(metric=:AUC, acceptance=[0.8, 1.25])]
)
result = simulate_trial(trial_spec, model_spec)
```

---

## Directory Structure

```
docs/examples/
├── README.md                    # This file
├── quickstart/                  # Getting started
├── models/                      # All PK/PD models
│   ├── pk/                      # PK models (8)
│   └── pkpd/                    # PKPD models (4)
├── population/                  # Population modeling
├── estimation/                  # Parameter estimation
├── nca/                         # Non-compartmental analysis
├── vpc/                         # Visual Predictive Checks
├── trial/                       # Trial simulation
├── import/                      # Model import
│   ├── nonmem/                  # NONMEM import
│   └── monolix/                 # Monolix import
├── data/                        # Data import
│   └── cdisc/                   # CDISC format
├── visualization/               # Plotting
├── sensitivity/                 # Sensitivity analysis
├── reproducibility/             # Artifacts and replay
├── use_cases/                   # End-to-end workflows
└── real_world_validation/       # Published datasets
```

---

## Running Examples

### Run All Examples

```bash
# Run all example validation
./docs/examples/run_all.sh

# Run specific category
julia --project=packages/core docs/examples/models/run_all.jl
python docs/examples/estimation/run_all.py
```

### Validate Outputs

```bash
# Validate against expected outputs
julia --project=packages/core docs/examples/validate_outputs.jl
```

---

## Contributing Examples

1. Each example should be in its own directory
2. Include Julia, Python, and CLI versions where applicable
3. Add expected output files for CI validation
4. Include a README.md explaining the example
5. Follow the naming convention: `01_descriptive_name/`

See the [CONTRIBUTING guide](https://github.com/shramish2057/NeoPKPD/blob/main/CONTRIBUTING.md) for guidelines.
