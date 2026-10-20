# NeoPKPD

**High-performance Julia engine for pharmacokinetic and pharmacodynamic modeling**

NeoPKPD is the computational core of the NeoPKPD platform, providing industry-standard ODE-based simulation, population modeling, parameter estimation, and comprehensive model diagnostics for pharmaceutical research and drug development.

## Features

### Pharmacokinetic Models
| Model | Parameters | Description |
|-------|------------|-------------|
| `OneCompIVBolus` | CL, V | One-compartment IV bolus/infusion |
| `OneCompOralFirstOrder` | Ka, CL, V | One-compartment oral absorption |
| `TwoCompIVBolus` | CL, V1, Q, V2 | Two-compartment IV |
| `TwoCompOral` | Ka, CL, V1, Q, V2 | Two-compartment oral |
| `ThreeCompIVBolus` | CL, V1, Q2, V2, Q3, V3 | Three-compartment IV |
| `TransitAbsorption` | N, Ktr, Ka, CL, V | Transit compartment absorption |
| `MichaelisMentenElimination` | Vmax, Km, V | Saturable (nonlinear) elimination |

### Pharmacodynamic Models
| Model | Parameters | Description |
|-------|------------|-------------|
| `DirectEmax` | E0, Emax, EC50 | Direct effect model (hyperbolic) |
| `SigmoidEmax` | E0, Emax, EC50, gamma | Sigmoid Emax (Hill equation) |
| `BiophaseEquilibration` | ke0, E0, Emax, EC50 | Effect compartment model |
| `IndirectResponseTurnover` | Kin, Kout, R0, Imax, IC50 | Indirect response turnover |
| `IndirectResponseIRM1-4` | Various | IRM1, IRM2, IRM3, IRM4 variants |
| `TransitCompartmentPD` | N, Ktr, ... | Transit compartment PD |
| `DiseaseProgression` | Various | Disease progression models |
| `ToleranceCounterRegulation` | Various | Tolerance and adaptation |

### Target-Mediated Drug Disposition (TMDD)
| Model | Description |
|-------|-------------|
| `FullTMDD` | Complete TMDD model for biologics |
| `TMDD_QSS` | Quasi-steady-state approximation |
| `TMDD_MM` | Michaelis-Menten approximation |

### Drug Interaction Models
- Bliss independence
- Competitive inhibition
- Receptor regulation

### Population Modeling
- **Inter-individual variability (IIV)** with log-normal distribution
- **Inter-occasion variability (IOV)** with occasion definitions
- **Static covariates**: Linear, power, exponential effects
- **Time-varying covariates**: Step and linear interpolation
- **BLQ handling**: M3/M4 methods for below limit of quantification data

### Parameter Estimation (NLME)
| Method | Description |
|--------|-------------|
| `FOCEIMethod` | First-Order Conditional Estimation with Interaction |
| `SAEMMethod` | Stochastic Approximation Expectation Maximization |
| `LaplacianMethod` | Laplacian approximation |
| `BayesianEstimation` | Full Bayesian with AdvancedHMC |

**Estimation Features:**
- Standard error computation (sandwich estimator, bootstrap)
- Covariance matrix estimation (diagonal, block-diagonal, full)
- Objective function value (OFV), AIC, BIC
- Eta shrinkage calculation
- Bootstrap analysis with confidence intervals
- Mixture models for subpopulation identification
- Model averaging (AIC, BIC, stacked weighting with NNLS)
- Stepwise covariate modeling (SCM) - forward/backward selection

### Non-Compartmental Analysis (NCA)
FDA/EMA-compliant NCA with:
- AUC methods: linear, log-linear, linear-up/log-down trapezoidal
- Cmax, Tmax, half-life, clearance, volume of distribution
- Bioavailability and accumulation ratios
- Lambda-z estimation with R-squared criteria
- Partial AUC calculations
- Steady-state metrics (AUC_tau, Cavg)

### Bioequivalence Analysis
- Average bioequivalence (ABE) with 90% CI
- Reference-scaled average bioequivalence (RSABE)
- Narrow therapeutic index drugs (NTID)
- Crossover study designs (2x2, 3x3, Williams)

### Visual Predictive Check (VPC)
- Standard VPC with prediction intervals
- Prediction-corrected VPC (pcVPC)
- Stratification by covariates
- Binning strategies: quantile, equal width, K-means, Jenks natural breaks
- Bootstrap confidence intervals
- LLOQ handling

### Residual Diagnostics
- Conditional weighted residuals (CWRES)
- Individual weighted residuals (IWRES)
- Normalized prediction distribution errors (NPDE)
- Population predictions (PRED) and individual predictions (IPRED)

### Sensitivity Analysis
| Type | Methods |
|------|---------|
| Local | Single parameter perturbation |
| Global | Sobol' indices (first, second, total order) |
| Screening | Morris method (elementary effects) |

### Clinical Trial Simulation
| Design | Description |
|--------|-------------|
| `ParallelDesign` | Parallel group studies |
| `CrossoverDesign` | Crossover studies (2x2, 3x3, Williams) |
| `DoseEscalationDesign` | Phase I (3+3, CRM, mTPI, BOIN) |
| `BioequivalenceDesign` | BE studies |
| `AdaptiveDesign` | Response-adaptive, sample size re-estimation |

**Trial Features:**
- Virtual population generation
- Power analysis with replicates
- Dosing regimen builders (QD, BID, TID, loading dose)
- Interim analysis with alpha spending (O'Brien-Fleming, Pocock)
- Futility boundaries

### Model Import
- NONMEM control stream (.ctl) parsing
- Monolix model (.mlxtran) parsing
- Automatic model structure detection
- Parameter extraction and mapping

### Data Import
- CDISC/SDTM format support
  - PC (Pharmacokinetic Concentrations) domain
  - EX (Exposure) domain
  - DM (Demographics) domain
- XPT (SAS Transport) file reader

### Residual Error Models
| Model | Formula |
|-------|---------|
| `AdditiveError` | Y = F + eps |
| `ProportionalError` | Y = F * (1 + eps) |
| `CombinedError` | Y = F + F*eps1 + eps2 |
| `ExponentialError` | Y = F * exp(eps) |

### Compliance (FDA 21 CFR Part 11)
- Digital signatures for artifacts
- Audit trail generation
- Environment capture (Julia version, package versions, OS)
- Validation reports
- Schema validation with integrity checks
- Content hashing (SHA-256)

### Serialization & Reproducibility
- JSON artifact format with semantic versioning
- Deterministic replay system
- Golden artifact validation
- Schema migration support
- Compliance metadata in artifacts

## Installation

```julia
using Pkg
Pkg.add("NeoPKPD")

using NeoPKPD
```

For development:
```julia
using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD
```

## Quick Start

### Single Subject Simulation

```julia
using NeoPKPD

# Define model specification
spec = ModelSpec(
    OneCompIVBolus(),
    "example",
    OneCompIVBolusParams(5.0, 50.0),  # CL=5 L/h, V=50 L
    [DoseEvent(0.0, 100.0)]            # 100 mg bolus at t=0
)

# Define simulation grid and solver
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run simulation
result = simulate(spec, grid, solver)

# Access results
println(result.t)                    # Time points
println(result.observations[:conc])  # Concentrations
println(result.states[:A_central])   # Amount in central compartment
```

### IV Infusion

```julia
# 100 mg infused over 1 hour (duration > 0)
doses = [DoseEvent(0.0, 100.0, 1.0)]  # time, amount, duration

spec = ModelSpec(
    OneCompIVBolus(),
    "infusion",
    OneCompIVBolusParams(5.0, 50.0),
    doses
)

result = simulate(spec, grid, solver)
```

### Population Simulation with IIV

```julia
# Define IIV specification
iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.3, :V => 0.2),  # Omega values (CV on log scale)
    UInt64(12345),                 # Seed for reproducibility
    100                            # Number of subjects
)

# Create population spec
pop_spec = PopulationSpec(spec, iiv, nothing, nothing, IndividualCovariates[])

# Run population simulation
pop_result = simulate_population(pop_spec, grid, solver)

# Access summaries
println(pop_result.summaries[:conc].mean)
println(pop_result.summaries[:conc].quantiles["0.05"])
println(pop_result.summaries[:conc].quantiles["0.95"])
```

### PKPD Simulation

```julia
# PK model
pk_spec = ModelSpec(
    OneCompIVBolus(),
    "pk",
    OneCompIVBolusParams(5.0, 50.0),
    [DoseEvent(0.0, 100.0)]
)

# PD model (Sigmoid Emax)
pd_spec = PDSpec(
    SigmoidEmax(),
    "pd",
    SigmoidEmaxParams(0.0, 100.0, 10.0, 2.0),  # E0, Emax, EC50, gamma
    :conc,       # Input observation (from PK)
    :effect      # Output observation name
)

# Coupled simulation
result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
println(result.observations[:effect])
```

### Parameter Estimation (FOCE-I)

```julia
using LinearAlgebra

config = EstimationConfig(
    FOCEIMethod(max_inner_iter=100, inner_tol=1e-6, centered=false),
    theta_init=[5.0, 50.0],      # Initial CL, V
    theta_lower=[0.1, 1.0],
    theta_upper=[100.0, 500.0],
    omega_init=diagm([0.09, 0.04]),
    omega_structure=:diagonal,
    sigma_init=ResidualErrorSpec(ProportionalError(), (sigma=0.1,), :conc, UInt64(0)),
    max_iter=500,
    compute_se=true
)

result = estimate(observed_data, model_spec, config; grid=grid, solver=solver)

println("Theta: ", result.theta)
println("Theta SE: ", result.theta_se)
println("OFV: ", result.ofv)
println("AIC: ", result.aic)
```

### Non-Compartmental Analysis

```julia
times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
concs = [0.0, 8.5, 12.3, 9.8, 6.2, 3.1, 1.5, 0.4]
dose = 100.0

nca_result = run_nca(times, concs, dose; method=LinearUpLogDownMethod())

println("Cmax: ", nca_result.cmax)
println("Tmax: ", nca_result.tmax)
println("AUC(0-inf): ", nca_result.auc_0_inf)
println("Half-life: ", nca_result.half_life)
println("CL/F: ", nca_result.cl_f)
```

### Visual Predictive Check

```julia
vpc_config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95],
    ci_level=0.95,
    binning=QuantileBinning(n_bins=10),
    prediction_corrected=false,
    stratify_by=Symbol[],
    lloq=nothing,
    n_bootstrap=500,
    seed=UInt64(12345)
)

vpc_result = compute_vpc(observed_data, pop_spec, grid, solver; config=vpc_config)
```

### Global Sensitivity Analysis (Sobol')

```julia
# Define parameter ranges
param_ranges = Dict(
    :CL => (1.0, 20.0),
    :V => (10.0, 100.0)
)

# Run Sobol' analysis
sobol_result = run_sobol_sensitivity(
    model_spec, grid, solver,
    param_ranges, :conc;
    n_samples=1024,
    seed=UInt64(42)
)

println("First-order indices: ", sobol_result.first_order)
println("Total-order indices: ", sobol_result.total_order)
```

### CDISC Data Import

```julia
# Read CDISC domains
dataset = read_cdisc_csv("pc.csv", "ex.csv", "dm.csv")

# Or read SAS transport files
dataset = read_cdisc_xpt("pc.xpt", "ex.xpt", "dm.xpt")

# Validate dataset
warnings = validate_cdisc_dataset(dataset)

# Convert to NeoPKPD format
pop_spec, observed = cdisc_to_population(dataset, model_spec)
```

### NONMEM/Monolix Import

```julia
# Parse NONMEM control file
nmctl = parse_nonmem_control("run001.ctl")
model_spec, pop_spec, mapping = convert_nonmem_to_neopkpd(nmctl)

# Parse Monolix project
mlx = parse_monolix_project("project.mlxtran")
model_spec, pop_spec, mapping = convert_monolix_to_neopkpd(mlx)
```

### Covariate Models

```julia
# Static covariates
cov_model = CovariateModel([
    CovariateEffect(:CL, :WT, PowerCovariate(0.75, 70.0)),    # Allometric
    CovariateEffect(:V, :WT, PowerCovariate(1.0, 70.0)),
    CovariateEffect(:CL, :CRCL, LinearCovariate(0.5, 100.0)), # Renal function
])

covariates = [
    IndividualCovariates(Dict(:WT => 80.0, :CRCL => 90.0), nothing),
    IndividualCovariates(Dict(:WT => 65.0, :CRCL => 120.0), nothing),
]

pop_spec = PopulationSpec(base_spec, iiv, nothing, cov_model, covariates)

# Time-varying covariates
tv_cov = TimeVaryingCovariate(:WT, [0.0, 24.0, 48.0], [70.0, 71.0, 72.0], :linear)
```

### Clinical Trial Simulation

```julia
# Define treatment arms
arm1 = TreatmentArm("Placebo", model_spec_placebo, DosingRegimen(QD(), 0.0, 28); n_subjects=50)
arm2 = TreatmentArm("Active", model_spec_active, DosingRegimen(QD(), 100.0, 28); n_subjects=50)

# Create trial specification
trial_spec = TrialSpec(
    "Phase 2 Study",
    ParallelDesign(2),
    [arm1, arm2];
    duration_days=28,
    pk_sampling_times=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
    endpoints=[PKEndpoint(:auc_0_inf)],
    n_replicates=100,
    seed=UInt64(12345)
)

# Run power simulation
result = run_power_simulation(trial_spec)
println("Power: ", result.power_estimates[:auc_0_inf])
```

### Serialization

```julia
# Serialize to artifact
artifact = serialize_single_artifact(spec, grid, solver, result)
write_artifact("simulation.json", artifact)

# Deserialize and replay
loaded = read_artifact("simulation.json")
replayed = replay_artifact(loaded)

# Population artifacts
pop_artifact = serialize_population_artifact(pop_spec, grid, solver, pop_result)
```

## Key Types

| Type | Description |
|------|-------------|
| `ModelSpec` | Complete model specification (kind, params, doses) |
| `SimGrid` | Time domain (t0, t1, saveat points) |
| `SolverSpec` | ODE solver configuration |
| `SimResult` | Single simulation output |
| `PopulationSpec` | Population model (IIV, IOV, covariates) |
| `PopulationResult` | Population simulation output with summaries |
| `DoseEvent` | Dose administration (time, amount, duration) |
| `EstimationConfig` | NLME estimation settings |
| `EstimationResult` | Fitted parameters with diagnostics |
| `VPCConfig` | VPC computation settings |
| `VPCResult` | VPC bins and statistics |
| `NCAResult` | Non-compartmental analysis results |
| `TrialSpec` | Clinical trial specification |
| `TrialResult` | Trial simulation results |

## Module Structure

```
NeoPKPD/
├── src/
│   ├── pk/              # PK model implementations
│   ├── pd/              # PD model implementations
│   ├── tmdd/            # TMDD models for biologics
│   ├── engine/          # Core simulation engine
│   ├── specs/           # Type definitions
│   ├── estimation/      # FOCE, SAEM, Laplacian, Bayesian
│   ├── nca/             # Non-compartmental analysis
│   ├── trial/           # Clinical trial simulation
│   ├── analysis/        # VPC, sensitivity, residuals
│   ├── import/          # NONMEM/Monolix parsers
│   ├── data/            # CDISC handling
│   ├── compliance/      # FDA 21 CFR Part 11
│   └── serialization/   # JSON I/O
└── test/                # Test suite (5400+ tests)
```

## Semantic Versions

| Component | Version | Description |
|-----------|---------|-------------|
| Event Semantics | 1.0.0 | Dose event handling behavior |
| Solver Semantics | 1.0.0 | ODE solver configuration |
| Artifact Schema | 1.0.0 | JSON serialization format |

## Testing

```bash
julia --project=packages/core -e 'using Pkg; Pkg.test()'
```

## License

MIT License - see repository root for details.
