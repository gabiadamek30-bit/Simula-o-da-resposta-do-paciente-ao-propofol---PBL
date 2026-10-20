# Julia Documentation

Welcome to the Julia documentation for NeoPKPD. The Julia core (`NeoPKPD.jl`) provides the foundation for all PK/PD modeling capabilities.

---

## Why Julia?

NeoPKPD uses Julia as its core language for several reasons:

- **Performance**: Near-native speed through JIT compilation
- **Mathematical Expressiveness**: Natural syntax for differential equations
- **Ecosystem**: World-class ODE solvers via DifferentialEquations.jl
- **Multiple Dispatch**: Flexible, extensible type system

---

## Quick Start

```julia
using NeoPKPD

# One-compartment IV bolus simulation
params = OneCompIVBolusParams(5.0, 50.0)  # CL=5 L/h, V=50 L
doses = [DoseEvent(0.0, 100.0)]            # 100 mg at t=0
spec = ModelSpec(OneCompIVBolus(), "example", params, doses)

grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

result = simulate(spec, grid, solver)
println(result.observations[:conc])
```

---

## Documentation Sections

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } **Tutorial**

    ---

    Step-by-step introduction to NeoPKPD Julia API

    [:octicons-arrow-right-24: Start Tutorial](tutorial.md)

-   :material-cube-outline:{ .lg .middle } **Models**

    ---

    Complete reference for PK and PD models

    [:octicons-arrow-right-24: Models Reference](models/index.md)

-   :material-account-group:{ .lg .middle } **Population Modeling**

    ---

    IIV, IOV, covariates, and residual error

    [:octicons-arrow-right-24: Population](population/index.md)

-   :material-calculator:{ .lg .middle } **NCA**

    ---

    Non-compartmental analysis

    [:octicons-arrow-right-24: NCA Reference](nca/index.md)

-   :material-chart-box:{ .lg .middle } **Parameter Estimation**

    ---

    FOCE-I, SAEM, and Laplacian methods

    [:octicons-arrow-right-24: Estimation](estimation/index.md)

-   :material-chart-scatter-plot:{ .lg .middle } **Visual Predictive Check**

    ---

    VPC, pcVPC, and stratification

    [:octicons-arrow-right-24: VPC Reference](vpc/index.md)

-   :material-file-import:{ .lg .middle } **Model Import**

    ---

    NONMEM and Monolix file parsing

    [:octicons-arrow-right-24: Import Reference](import/index.md)

-   :material-flask:{ .lg .middle } **Clinical Trials**

    ---

    Trial simulation and power analysis

    [:octicons-arrow-right-24: Trial Reference](trial/index.md)

</div>

---

## Core Types

### ModelSpec

The central type for defining a simulation:

```julia
struct ModelSpec{M<:ModelKind, P<:AbstractParams}
    model::M           # Model type (e.g., OneCompIVBolus)
    name::String       # Simulation identifier
    params::P          # Model parameters
    doses::Vector{DoseEvent}
end
```

### SimGrid

Defines the time domain:

```julia
struct SimGrid
    t0::Float64        # Start time
    t1::Float64        # End time
    saveat::Vector{Float64}  # Output time points
end
```

### SolverSpec

Configures the ODE solver:

```julia
struct SolverSpec
    alg::Symbol        # Algorithm (:Tsit5, :Rosenbrock23, etc.)
    reltol::Float64    # Relative tolerance
    abstol::Float64    # Absolute tolerance
    maxiters::Int      # Maximum iterations
end
```

### SimResult

Simulation output:

```julia
struct SimResult
    t::Vector{Float64}
    states::Dict{Symbol, Vector{Float64}}
    observations::Dict{Symbol, Vector{Float64}}
    metadata::Dict{Symbol, Any}
end
```

---

## Available Models

### Pharmacokinetic Models

| Model Type | Parameters | Description |
|------------|------------|-------------|
| `OneCompIVBolus` | CL, V | IV bolus, first-order elimination |
| `OneCompOralFirstOrder` | Ka, CL, V | Oral with first-order absorption |
| `TwoCompIVBolus` | CL, V1, Q, V2 | Two-compartment IV |
| `TwoCompOral` | Ka, CL, V1, Q, V2 | Two-compartment oral |
| `ThreeCompIVBolus` | CL, V1, Q2, V2, Q3, V3 | Three-compartment IV |
| `TransitAbsorption` | N, Ktr, Ka, CL, V | Transit compartment absorption |
| `MichaelisMentenElimination` | Vmax, Km, V | Saturable elimination |

### Pharmacodynamic Models

| Model Type | Parameters | Description |
|------------|------------|-------------|
| `DirectEmax` | E0, Emax, EC50 | Direct effect model |
| `SigmoidEmax` | E0, Emax, EC50, gamma | Hill equation |
| `BiophaseEquilibration` | ke0, E0, Emax, EC50 | Effect compartment |
| `IndirectResponseTurnover` | Kin, Kout, R0, Imax, IC50 | Indirect response |

---

## Key Functions

### Simulation

```julia
# Single subject simulation
result = simulate(spec::ModelSpec, grid::SimGrid, solver::SolverSpec)

# Population simulation
result = simulate_population(pop_spec::PopulationSpec, grid, solver)
```

### Parameter Estimation

```julia
# FOCE-I estimation
result = estimate(data, spec, FOCEConfig())

# SAEM estimation
result = estimate(data, spec, SAEMConfig())
```

### NCA

```julia
# Non-compartmental analysis
result = run_nca(times, conc, dose; config=NCAConfig())
```

### VPC

```julia
# Visual predictive check
vpc_result = compute_vpc(observed, simulated; config=VPCConfig())
```

---

## Next Steps

- [Start the Tutorial →](tutorial.md)
- [Explore PK Models →](models/pk/onecomp-iv-bolus.md)
- [Learn Population Modeling →](population/index.md)
