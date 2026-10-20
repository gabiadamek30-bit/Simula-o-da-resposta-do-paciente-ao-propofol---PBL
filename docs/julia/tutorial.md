# Tutorial: Getting Started with Julia

This tutorial provides a comprehensive introduction to using NeoPKPD with Julia. By the end, you'll be able to run simulations, analyze populations, and perform parameter estimation.

---

## Prerequisites

Ensure you have:

- Julia 1.10 or later installed
- NeoPKPD repository cloned
- Julia dependencies installed

```bash
julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'
```

---

## Part 1: Your First Simulation

### Loading the Package

```julia
# Start Julia with the project
# julia --project=packages/core

using NeoPKPD
```

### One-Compartment IV Bolus

The simplest PK model: a single compartment with first-order elimination.

```julia
# Define parameters
# CL = 5 L/h (clearance)
# V = 50 L (volume of distribution)
params = OneCompIVBolusParams(5.0, 50.0)

# Define a single 100 mg dose at time 0
doses = [DoseEvent(0.0, 100.0)]

# Create model specification
spec = ModelSpec(
    OneCompIVBolus(),      # Model type
    "tutorial_sim",        # Simulation name
    params,                # Parameters
    doses                  # Dose events
)

# Define time grid: 0 to 24 hours, output every hour
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))

# Configure ODE solver
# Tsit5 is the default 5th order Runge-Kutta
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run simulation
result = simulate(spec, grid, solver)

# Access results
println("Time points: ", length(result.t))
println("First 5 concentrations: ", result.observations[:conc][1:5])
```

**Expected output:**
```
Time points: 25
First 5 concentrations: [2.0, 1.8097..., 1.6375..., 1.4816..., 1.3406...]
```

### Understanding the Result

```julia
# The SimResult contains:
# - t: Time points
# - states: Internal ODE state variables
# - observations: Derived quantities (concentrations, effects)
# - metadata: Additional information

# Time points
result.t  # [0.0, 1.0, 2.0, ..., 24.0]

# Concentrations (primary observation)
result.observations[:conc]  # [2.0, 1.81, 1.64, ...]

# State variables (amount in central compartment)
result.states[:A_central]  # [100.0, 90.5, 81.9, ...]
```

### PK Metrics

```julia
# Calculate half-life
t_half = log(2) * params.V / params.CL
println("Half-life: ", t_half, " hours")

# Find Cmax (for IV bolus, it's at t=0)
cmax = maximum(result.observations[:conc])
println("Cmax: ", cmax, " mg/L")

# Calculate AUC using trapezoidal rule
function auc_trap(t, c)
    auc = 0.0
    for i in 2:length(t)
        auc += (t[i] - t[i-1]) * (c[i] + c[i-1]) / 2
    end
    return auc
end

auc = auc_trap(result.t, result.observations[:conc])
println("AUC0-24: ", auc, " mg·h/L")
```

---

## Part 2: Multiple Doses

### Repeated Dosing

```julia
# Multiple doses: 100 mg every 12 hours for 3 days
doses = [
    DoseEvent(0.0, 100.0),
    DoseEvent(12.0, 100.0),
    DoseEvent(24.0, 100.0),
    DoseEvent(36.0, 100.0),
    DoseEvent(48.0, 100.0),
    DoseEvent(60.0, 100.0),
]

spec = ModelSpec(OneCompIVBolus(), "multiple_dose", params, doses)

# Extend simulation to 72 hours
grid = SimGrid(0.0, 72.0, collect(0.0:0.5:72.0))

result = simulate(spec, grid, solver)

# Find steady-state trough
# Trough is just before the next dose
trough_times = [11.5, 23.5, 35.5, 47.5, 59.5, 71.5]
for t in trough_times
    idx = findfirst(x -> x ≈ t, result.t)
    if idx !== nothing
        println("Trough at t=$(t): $(result.observations[:conc][idx]) mg/L")
    end
end
```

### IV Infusion

```julia
# 100 mg infused over 1 hour
doses = [DoseEvent(0.0, 100.0, 1.0)]  # duration = 1.0 hour

spec = ModelSpec(OneCompIVBolus(), "infusion", params, doses)
grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))

result = simulate(spec, grid, solver)

# Cmax occurs at end of infusion
idx_end_infusion = findfirst(x -> x ≈ 1.0, result.t)
println("Cmax at end of infusion: ", result.observations[:conc][idx_end_infusion])
```

---

## Part 3: Two-Compartment Model

### Bi-Exponential Kinetics

```julia
# Two-compartment model
# CL = 10 L/h, V1 = 20 L (central)
# Q = 15 L/h (inter-compartmental clearance)
# V2 = 50 L (peripheral)
params = TwoCompIVBolusParams(10.0, 20.0, 15.0, 50.0)

doses = [DoseEvent(0.0, 500.0)]
spec = ModelSpec(TwoCompIVBolus(), "twocomp", params, doses)

# Dense output to see distribution phase
grid = SimGrid(0.0, 48.0, collect(0.0:0.1:48.0))

result = simulate(spec, grid, solver)

# You'll see bi-exponential decline:
# - Fast initial decline (distribution phase, α)
# - Slower terminal decline (elimination phase, β)
println("C at t=0.5h: ", result.observations[:conc][6])   # During distribution
println("C at t=24h: ", result.observations[:conc][241])  # Terminal phase
```

### Accessing Both Compartments

```julia
# Central compartment amount
A_central = result.states[:A_central]

# Peripheral compartment amount
A_peripheral = result.states[:A_peripheral]

# Verify mass balance
total_amount_t0 = A_central[1] + A_peripheral[1]
println("Total amount at t=0: ", total_amount_t0, " mg")
```

---

## Part 4: Oral Absorption

### First-Order Absorption

```julia
# One-compartment oral
# Ka = 1.5 /h, CL = 5 L/h, V = 50 L
params = OneCompOralFirstOrderParams(1.5, 5.0, 50.0)

doses = [DoseEvent(0.0, 200.0)]  # Oral dose of 200 mg
spec = ModelSpec(OneCompOralFirstOrder(), "oral", params, doses)

grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))

result = simulate(spec, grid, solver)

# Find Tmax and Cmax
cmax, idx_cmax = findmax(result.observations[:conc])
tmax = result.t[idx_cmax]

println("Tmax: ", tmax, " hours")
println("Cmax: ", cmax, " mg/L")
```

### Transit Compartment Absorption

For drugs with complex absorption (delayed peak, GI transit):

```julia
# Transit absorption model
# 5 transit compartments, Ktr = 0.5/h, Ka = 2.0/h
# CL = 10 L/h, V = 70 L
params = TransitAbsorptionParams(5, 0.5, 2.0, 10.0, 70.0)

doses = [DoseEvent(0.0, 300.0)]
spec = ModelSpec(TransitAbsorption(), "transit", params, doses)

grid = SimGrid(0.0, 24.0, collect(0.0:0.1:24.0))

result = simulate(spec, grid, solver)

# Delayed Tmax due to transit compartments
cmax, idx_cmax = findmax(result.observations[:conc])
println("Tmax: ", result.t[idx_cmax], " hours")  # Will be later than simple oral
```

---

## Part 5: PK-PD Modeling

### Direct Emax Model

```julia
# Direct Emax: effect is immediate function of concentration
# PK: CL = 5, V = 50
# PD: E0 = 0 (baseline), Emax = 100, EC50 = 2 mg/L

pk_params = OneCompIVBolusParams(5.0, 50.0)
pd_params = DirectEmaxParams(0.0, 100.0, 2.0)

# Combined PKPD parameters
params = PKPDDirectEmaxParams(pk_params, pd_params)

doses = [DoseEvent(0.0, 100.0)]
spec = ModelSpec(DirectEmax(), "pkpd", params, doses)

grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))

result = simulate(spec, grid, solver)

# Both concentration and effect are available
println("Concentration: ", result.observations[:conc][1:5])
println("Effect: ", result.observations[:effect][1:5])
```

### Indirect Response Model

```julia
# Indirect response (inhibition of Kout)
# Kin = 10, Kout = 0.5, baseline R0 = 20
# Imax = 0.8, IC50 = 1.0 mg/L

pk_params = OneCompIVBolusParams(5.0, 50.0)
pd_params = IndirectResponseTurnoverParams(10.0, 0.5, 20.0, 0.8, 1.0)

params = PKPDIndirectResponseParams(pk_params, pd_params)

doses = [DoseEvent(0.0, 100.0)]
spec = ModelSpec(IndirectResponseTurnover(), "indirect", params, doses)

grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))

result = simulate(spec, grid, solver)

# Response changes slowly (indirect mechanism)
println("Baseline response: ", result.observations[:response][1])
println("Response at 24h: ", result.observations[:response][25])
println("Response at 72h: ", result.observations[:response][73])
```

---

## Part 6: Population Simulation

### Inter-Individual Variability (IIV)

```julia
# Define typical parameters
typical_params = OneCompIVBolusParams(5.0, 50.0)

# Define omega matrix (variance of log-normal random effects)
# 30% CV on CL, 20% CV on V
# CV = sqrt(exp(omega^2) - 1) ≈ omega for small omega
omega = OmegaMatrix([
    0.09 0.0;   # omega_CL^2 = 0.09 -> ~30% CV
    0.0  0.04   # omega_V^2 = 0.04 -> ~20% CV
])

# Create population specification
doses = [DoseEvent(0.0, 100.0)]
base_spec = ModelSpec(OneCompIVBolus(), "pop", typical_params, doses)

pop_spec = PopulationSpec(
    base_spec,
    100,           # 100 individuals
    omega,
    12345          # seed for reproducibility
)

grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))

# Simulate population
pop_result = simulate_population(pop_spec, grid, solver)

# Access individual results
for i in 1:5
    cmax = maximum(pop_result.individuals[i].observations[:conc])
    println("Individual $i Cmax: $(round(cmax, digits=2)) mg/L")
end

# Population summary statistics
summary = pop_result.summaries[:conc]
println("\nPopulation median: ", summary.median)
println("5th percentile: ", summary.quantiles[0.05])
println("95th percentile: ", summary.quantiles[0.95])
```

### With Covariates

```julia
# Weight-based dosing with allometric scaling
# CL = 5 * (WT/70)^0.75
# V = 50 * (WT/70)^1.0

# Define covariate model
covariate_model = CovariateModel([
    CovariateEffect(:CL, :WT, 70.0, :power, 0.75),
    CovariateEffect(:V, :WT, 70.0, :power, 1.0)
])

# Generate subjects with varying weights
using Random
Random.seed!(42)
n_subjects = 50
weights = 50.0 .+ 40.0 .* rand(n_subjects)  # 50-90 kg

covariates = [Dict(:WT => w) for w in weights]

# Population with covariates
pop_spec = PopulationSpec(
    base_spec,
    n_subjects,
    omega,
    12345,
    covariate_model,
    covariates
)

pop_result = simulate_population(pop_spec, grid, solver)

# Show effect of weight on CL
for i in 1:5
    wt = covariates[i][:WT]
    cl = pop_result.realized_params[i][:CL]
    println("Subject $i: WT=$(round(wt, digits=1)) kg, CL=$(round(cl, digits=2)) L/h")
end
```

---

## Part 7: Non-Compartmental Analysis

### Basic NCA

```julia
# Run NCA on simulation result
times = result.t
conc = result.observations[:conc]
dose = 100.0

nca_result = run_nca(times, conc, dose)

println("Cmax: ", nca_result.cmax)
println("Tmax: ", nca_result.tmax)
println("AUC0-t: ", nca_result.auc_0_t)
println("AUC0-inf: ", nca_result.auc_0_inf)
println("Half-life: ", nca_result.t_half)
println("CL/F: ", nca_result.cl_f)
println("Vz/F: ", nca_result.vz_f)
```

### NCA Configuration

```julia
config = NCAConfig(
    method = :log_linear,           # AUC calculation method
    lambda_z_min_points = 3,        # Minimum points for lambda_z
    lambda_z_r2_threshold = 0.9,    # R² threshold
    extrapolation_max_pct = 20.0,   # Warning threshold
    blq_handling = :zero            # BLQ handling
)

nca_result = run_nca(times, conc, dose; config=config)
```

---

## Part 8: Artifacts and Reproducibility

### Writing Artifacts

```julia
# Save simulation as artifact
artifact = write_artifact("simulation.json", spec, grid, solver, result)

println("Artifact saved with version: ", artifact.schema_version)
```

### Replaying Artifacts

```julia
# Load and replay artifact
replayed = replay_artifact("simulation.json")

# Verify results match
@assert replayed.observations[:conc] == result.observations[:conc]
println("Replay successful - results match!")
```

---

## Summary

In this tutorial, you learned:

1. **Basic Simulation**: `ModelSpec`, `SimGrid`, `SolverSpec`, `simulate()`
2. **Multiple Doses**: Using `DoseEvent` arrays, IV infusions
3. **Multi-Compartment Models**: Two-compartment, three-compartment
4. **Oral Absorption**: First-order and transit compartment
5. **PK-PD**: Direct Emax and indirect response models
6. **Population Simulation**: IIV with omega matrix, covariates
7. **NCA**: Exposure metrics and configuration
8. **Reproducibility**: Artifacts and replay

---

## Next Steps

- [PK Models Reference](models/pk/onecomp-iv-bolus.md) - Detailed model documentation
- [PD Models Reference](models/pd/direct-emax.md) - Pharmacodynamic models
- [Population Modeling](population/index.md) - Advanced IIV/IOV/covariates
- [Parameter Estimation](estimation/index.md) - FOCE-I, SAEM, Laplacian
