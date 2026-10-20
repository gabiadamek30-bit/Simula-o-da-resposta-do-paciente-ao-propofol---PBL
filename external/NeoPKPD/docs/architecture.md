# Architecture

NeoPKPD is designed with transparency, reproducibility, and validation as core principles.

## Design Philosophy

### 1. Pure Data Specifications

All simulation inputs are **immutable data structures** with no hidden state:

```julia
# Everything needed is explicit in the spec
spec = ModelSpec(
    OneCompIVBolus(),           # Model type
    "my_model",                 # Name
    OneCompIVBolusParams(5, 50), # Parameters
    [DoseEvent(0, 100)]         # Doses
)
```

**Benefits**:

- Complete traceability in artifacts
- No side effects or global state
- Reproducible across versions

### 2. Validated Models

Every model implementation includes:

- Mathematical equations documented in code
- Parameter validation (constraints checked at simulation time)
- Unit tests against analytical solutions
- Golden artifacts for regression testing

### 3. Semantic Versioning

Three independent version numbers track numerical behavior:

| Version | Scope | Example Change |
|---------|-------|----------------|
| Event Semantics | Dose handling | How doses at t=0 are applied |
| Solver Semantics | ODE solving | Tolerance interpretation |
| Artifact Schema | JSON format | Field names, structure |

**Version Bump Policy**: Any change to numerical output requires a semantic version bump.

### 4. Artifact-Driven Validation

Every simulation can produce a JSON artifact that:

- Contains complete inputs and outputs
- Includes semantic version fingerprint
- Can be replayed to verify reproducibility

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        User Interfaces                          │
├─────────────────┬─────────────────┬─────────────────────────────┤
│   Julia API     │   Python API    │         CLI                 │
│   (direct)      │   (juliacall)   │    (NeoPKPDCLI)           │
└────────┬────────┴────────┬────────┴────────────┬────────────────┘
         │                 │                     │
         ▼                 ▼                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                     NeoPKPD                                │
├─────────────────────────────────────────────────────────────────┤
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────────┐ │
│  │   Models    │  │   Engine    │  │    Serialization        │ │
│  │  (PK/PD)    │  │  (Solver)   │  │    (JSON Artifacts)     │ │
│  └─────────────┘  └─────────────┘  └─────────────────────────┘ │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────────┐ │
│  │ Population  │  │ Sensitivity │  │    Analysis             │ │
│  │ (IIV/IOV)   │  │ (Perturb)   │  │    (Metrics)            │ │
│  └─────────────┘  └─────────────┘  └─────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                  DifferentialEquations.jl                       │
│                      (ODE Solver)                               │
└─────────────────────────────────────────────────────────────────┘
```

---

## Module Structure

### Core (`packages/core/`)

```
src/
├── NeoPKPD.jl          # Module definition, exports
├── specs/
│   ├── specs.jl             # Core type definitions
│   ├── sensitivity.jl       # Perturbation types
│   └── time_covariates.jl   # Time-varying covariate types
├── models/
│   ├── onecomp_iv_bolus.jl      # IV bolus model
│   ├── onecomp_oral_first_order.jl  # Oral model
│   └── pk_interface.jl          # Common PK interface
├── pd/
│   ├── direct_emax.jl           # Direct Emax PD
│   └── indirect_response_turnover.jl  # Indirect response PD
├── engine/
│   ├── events.jl            # Dose normalization
│   ├── callbacks.jl         # ODE callbacks
│   ├── solve.jl             # Main simulation
│   ├── pkpd.jl              # Sequential PKPD
│   ├── pkpd_coupled.jl      # Coupled PKPD
│   ├── population.jl        # Population simulation
│   ├── iov.jl               # Inter-occasion variability
│   ├── covariates.jl        # Covariate application
│   ├── time_covariates.jl   # Time-varying evaluation
│   ├── segment_sim.jl       # Segmented PK simulation
│   ├── segment_sim_pkpd.jl  # Segmented PKPD simulation
│   ├── sensitivity.jl       # Single sensitivity
│   ├── sensitivity_population.jl  # Population sensitivity
│   ├── sensitivity_metrics.jl     # Metric computation
│   ├── perturb.jl           # Perturbation application
│   ├── semantics.jl         # Event semantics
│   ├── solver_semantics.jl  # Solver semantics
│   └── semantics_fingerprint.jl  # Version fingerprinting
├── serialization/
│   ├── schema.jl            # Schema version
│   ├── serialize.jl         # Single serialization
│   ├── deserialize.jl       # Single deserialization
│   ├── serialize_population.jl    # Population serialization
│   ├── deserialize_population.jl  # Population deserialization
│   ├── serialize_sensitivity.jl   # Sensitivity serialization
│   └── deserialize_sensitivity.jl # Sensitivity deserialization
└── analysis/
    └── exposure.jl          # Cmax, AUC
```

### CLI (`cli/NeoPKPDCLI/`)

```
src/
└── NeoPKPDCLI.jl           # CLI commands
```

### Python (`python/neopkpd/`)

```
neopkpd/
├── __init__.py              # Package exports
└── bridge.py                # Julia-Python bridge
```

---

## Data Flow

### Single Simulation

```
ModelSpec + SimGrid + SolverSpec
            │
            ▼
     ┌──────────────┐
     │   validate   │  Parameter constraints
     └──────┬───────┘
            │
            ▼
     ┌──────────────┐
     │ normalize_   │  Dose preprocessing
     │ doses_for_   │  (t0 doses → u0)
     │ sim          │
     └──────┬───────┘
            │
            ▼
     ┌──────────────┐
     │  callbacks   │  Dose callbacks for (t0,t1]
     └──────┬───────┘
            │
            ▼
     ┌──────────────┐
     │   solve      │  DifferentialEquations.jl
     └──────┬───────┘
            │
            ▼
     ┌──────────────┐
     │  extract     │  States → Observations
     │  results     │  Output grid alignment
     └──────┬───────┘
            │
            ▼
        SimResult
```

### Population Simulation

```
PopulationSpec + SimGrid + SolverSpec
            │
            ▼
     ┌──────────────┐
     │  sample_iiv  │  Generate η per individual
     └──────┬───────┘
            │
            ▼
     ┌──────────────┐
     │ apply_covs   │  Static covariate effects
     └──────┬───────┘
            │
            ▼
     ┌──────────────┐
     │ derive_occ   │  Occasion boundaries
     │ (if IOV)     │
     └──────┬───────┘
            │
            ▼
     ┌──────────────────────────────────┐
     │  For each individual:            │
     │  ┌────────────────────────────┐  │
     │  │ sample_iov (if IOV)        │  │
     │  │ segment_boundaries         │  │
     │  │ simulate_segmented_pk/pkpd │  │
     │  └────────────────────────────┘  │
     └──────────────┬───────────────────┘
                    │
                    ▼
     ┌──────────────┐
     │  compute_    │  Mean, median, quantiles
     │  summary     │
     └──────┬───────┘
            │
            ▼
     PopulationResult
```

---

## Type Hierarchy

### Model Kinds

```
ModelKind (abstract)
├── OneCompIVBolus
└── OneCompOralFirstOrder

PDModelKind (abstract)
├── DirectEmax
└── IndirectResponseTurnover
```

### Random Effect Kinds

```
RandomEffectKind (abstract)
└── LogNormalIIV
```

### Covariate Effect Kinds

```
CovariateEffectKind (abstract)
├── LinearCovariate
├── PowerCovariate
└── ExpCovariate
```

### Time Covariate Kinds

```
TimeCovariateKind (abstract)
├── StepTimeCovariate
└── LinearTimeCovariate
```

### Perturbation Kinds

```
PerturbationKind (abstract)
├── RelativePerturbation
└── AbsolutePerturbation
```

---

## Extension Points

### Adding a New PK Model

1. Create `src/models/newmodel.jl`
2. Define `NewModel <: ModelKind`
3. Define `NewModelParams` struct
4. Implement required interface:
   - `pk_validate(spec::ModelSpec{NewModel})`
   - `pk_param_tuple(spec::ModelSpec{NewModel})`
   - `pk_state_symbols(::NewModel)`
   - `pk_u0(spec::ModelSpec{NewModel}, grid)`
   - `pk_ode!(du, u, p, t, ::NewModel)`
   - `pk_conc(u, p, ::NewModel)`
   - `pk_dose_target_index(::NewModel)`
5. Add exports to `NeoPKPD.jl`
6. Add serialization support
7. Add tests and golden artifacts

### Adding a New PD Model

1. Create `src/pd/newpd.jl`
2. Define `NewPD <: PDModelKind`
3. Define `NewPDParams` struct
4. Implement `validate(spec::PDSpec{NewPD})`
5. For direct models: implement `evaluate(spec, input_series)`
6. For dynamic models: implement ODE interface
7. Add serialization support

### Adding a New Covariate Effect

1. Define `NewEffect <: CovariateEffectKind`
2. Add dispatch in `apply_covariates`
3. Add serialization support

---

## Validation Strategy

### Unit Tests

```julia
# test/models/test_onecomp_iv_bolus.jl
@test simulate(...) produces expected analytical solution
```

### Golden Artifacts

```
validation/golden/
├── pk_iv_bolus.json
├── pk_oral_first_order.json
├── population_iv_bolus.json
├── pkpd_direct_emax.json
└── sensitivity_single.json
```

### Replay Validation

```julia
artifact = read_execution_json("golden.json")
result = replay_execution(artifact)
@test result matches stored_result
```

### CI Pipeline

1. Run unit tests
2. Generate fresh artifacts
3. Replay all golden artifacts
4. Compare against stored results
5. Fail if any mismatch

---

## Performance Considerations

### ODE Solver

- Uses DifferentialEquations.jl for production-quality solvers
- Default Tsit5 (5th order Runge-Kutta) for non-stiff problems
- Rosenbrock23 available for stiff problems

### Population Simulation

- Sequential individual simulation (not parallel by default)
- Segmented simulation adds overhead for IOV/time-varying covariates
- Consider reducing population size for exploratory work

### Memory

- Each individual result stored separately
- Large populations with fine grids can consume significant memory
- Use sparse `saveat` grids when possible

---

## Thread Safety

- All simulation functions are thread-safe
- No global mutable state
- Safe to run parallel simulations from different threads

---

## Dependencies

### Julia

| Package | Purpose |
|---------|---------|
| DifferentialEquations.jl | ODE solving |
| SciMLBase | Scientific ML interface |
| JSON | Serialization |
| ArgParse | CLI parsing |

### Python

| Package | Purpose |
|---------|---------|
| juliacall | Julia-Python bridge |
| juliapkg | Julia environment management |
