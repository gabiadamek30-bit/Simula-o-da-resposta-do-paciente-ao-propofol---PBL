# Sensitivity Analysis

NeoPKPD provides tools for analyzing how parameter changes affect model outputs, essential for understanding model behavior and identifying critical parameters.

## Overview

Sensitivity analysis answers questions like:

- How much does concentration change if clearance increases by 10%?
- Which parameters have the largest impact on drug exposure?
- How does population variability affect therapeutic outcomes?

NeoPKPD supports both **single-subject** and **population-level** sensitivity analysis.

---

## Perturbation Types

### Relative Perturbation

Changes parameter by a fraction of its base value:

$$\theta_{new} = \theta_{base} \cdot (1 + \delta)$$

```julia
RelativePerturbation()
```

**Example**: `delta = 0.1` increases parameter by 10%

### Absolute Perturbation

Changes parameter by a fixed amount:

$$\theta_{new} = \theta_{base} + \delta$$

```julia
AbsolutePerturbation()
```

**Example**: `delta = 1.0` increases parameter by 1 unit

---

## Perturbation Specification

### Single Perturbation

```julia
struct Perturbation{K<:PerturbationKind}
    kind::K           # RelativePerturbation or AbsolutePerturbation
    param::Symbol     # Parameter to perturb
    delta::Float64    # Change magnitude
end
```

### Perturbation Plan

A plan groups multiple perturbations to be applied together:

```julia
struct PerturbationPlan
    name::String
    perturbations::Vector{Perturbation}  # Applied in sequence
end
```

---

## Sensitivity Metrics

NeoPKPD computes three metrics comparing base and perturbed output series:

```julia
struct SensitivityMetric
    max_abs_delta::Float64   # Maximum |base - pert|
    max_rel_delta::Float64   # Maximum |base - pert| / |base|
    l2_norm_delta::Float64   # sqrt(sum((base - pert)^2))
end
```

| Metric | Formula | Use Case |
|--------|---------|----------|
| `max_abs_delta` | $\max\|y_{base} - y_{pert}\|$ | Maximum deviation |
| `max_rel_delta` | $\max\frac{\|y_{base} - y_{pert}\|}{\|y_{base}\|}$ | Relative impact |
| `l2_norm_delta` | $\sqrt{\sum(y_{base} - y_{pert})^2}$ | Overall deviation |

---

## Single-Subject Sensitivity

### Function

```julia
run_sensitivity(spec, grid, solver; plan, observation=:conc)
```

| Parameter | Description |
|-----------|-------------|
| `spec` | Model specification |
| `grid` | Simulation grid |
| `solver` | Solver specification |
| `plan` | Perturbation plan |
| `observation` | Output to analyze (default `:conc`) |

### Result Structure

```julia
struct SensitivityResult
    plan::PerturbationPlan
    observation::Symbol
    base_metric_series::Vector{Float64}
    pert_metric_series::Vector{Float64}
    metrics::SensitivityMetric
    metadata::Dict{String, Any}
end
```

### Example: CL Sensitivity

```julia
using NeoPKPD

# Define model
params = OneCompIVBolusParams(5.0, 50.0)
spec = ModelSpec(
    OneCompIVBolus(),
    "sensitivity_example",
    params,
    [DoseEvent(0.0, 100.0)]
)

grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Perturbation plan: 10% increase in CL
plan = PerturbationPlan(
    "CL_10pct_increase",
    [Perturbation(RelativePerturbation(), :CL, 0.1)]
)

# Run sensitivity analysis
result = run_sensitivity(spec, grid, solver; plan=plan, observation=:conc)

println("Base concentrations: ", result.base_metric_series)
println("Perturbed concentrations: ", result.pert_metric_series)
println("Max absolute delta: ", result.metrics.max_abs_delta)
println("Max relative delta: ", result.metrics.max_rel_delta)
println("L2 norm delta: ", result.metrics.l2_norm_delta)
```

### Example: Multiple Parameter Perturbation

```julia
# Perturbation plan: increase CL by 10% and decrease V by 5%
plan = PerturbationPlan(
    "CL_up_V_down",
    [
        Perturbation(RelativePerturbation(), :CL, 0.1),   # +10% CL
        Perturbation(RelativePerturbation(), :V, -0.05)  # -5% V
    ]
)

result = run_sensitivity(spec, grid, solver; plan=plan)
```

### Example: Absolute Perturbation

```julia
# Perturbation plan: increase CL by 1 L/h
plan = PerturbationPlan(
    "CL_plus_1",
    [Perturbation(AbsolutePerturbation(), :CL, 1.0)]
)

result = run_sensitivity(spec, grid, solver; plan=plan)
```

---

## Population Sensitivity

Population sensitivity analysis evaluates how parameter perturbations affect population-level summaries (mean, quantiles).

### Function

```julia
run_population_sensitivity(pop, grid, solver; plan, observation=:conc, probs=[0.05, 0.95])
```

| Parameter | Description |
|-----------|-------------|
| `pop` | Population specification |
| `grid` | Simulation grid |
| `solver` | Solver specification |
| `plan` | Perturbation plan |
| `observation` | Output to analyze (default `:conc`) |
| `probs` | Quantile probabilities (default `[0.05, 0.95]`) |

### Result Structure

```julia
struct PopulationSensitivityResult
    plan::PerturbationPlan
    observation::Symbol
    probs::Vector{Float64}
    base_summary_mean::Vector{Float64}
    pert_summary_mean::Vector{Float64}
    base_quantiles::Dict{Float64, Vector{Float64}}
    pert_quantiles::Dict{Float64, Vector{Float64}}
    metrics_mean::SensitivityMetric
    metadata::Dict{String, Any}
end
```

### Example: Population CL Sensitivity

```julia
using NeoPKPD

# Define population
params = OneCompIVBolusParams(5.0, 50.0)
base_spec = ModelSpec(
    OneCompIVBolus(),
    "pop_sens",
    params,
    [DoseEvent(0.0, 100.0)]
)

iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.3, :V => 0.2),
    UInt64(12345),
    100
)

pop = PopulationSpec(base_spec, iiv, nothing, nothing, [])

grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Perturbation plan: 20% increase in typical CL
plan = PerturbationPlan(
    "CL_20pct",
    [Perturbation(RelativePerturbation(), :CL, 0.2)]
)

# Run population sensitivity
result = run_population_sensitivity(
    pop, grid, solver;
    plan=plan,
    observation=:conc,
    probs=[0.05, 0.50, 0.95]
)

println("Base mean: ", result.base_summary_mean)
println("Perturbed mean: ", result.pert_summary_mean)
println("Mean metrics: ", result.metrics_mean)

# Compare quantiles
println("Base 5th percentile: ", result.base_quantiles[0.05])
println("Pert 5th percentile: ", result.pert_quantiles[0.05])
```

---

## Helper Functions

### Apply Single Perturbation

```julia
new_params = apply_perturbation(params, perturbation)
```

### Apply Perturbation Plan

```julia
new_params = apply_plan(params, plan)
```

### Compute Metrics

```julia
metrics = compute_metrics(base_series, pert_series)
```

---

## Sensitivity Analysis Workflow

### 1. Define Baseline Model

```julia
params = OneCompIVBolusParams(5.0, 50.0)
spec = ModelSpec(OneCompIVBolus(), "baseline", params, [DoseEvent(0.0, 100.0)])
```

### 2. Create Perturbation Plans

```julia
plans = [
    PerturbationPlan("CL+10%", [Perturbation(RelativePerturbation(), :CL, 0.1)]),
    PerturbationPlan("CL-10%", [Perturbation(RelativePerturbation(), :CL, -0.1)]),
    PerturbationPlan("V+10%", [Perturbation(RelativePerturbation(), :V, 0.1)]),
    PerturbationPlan("V-10%", [Perturbation(RelativePerturbation(), :V, -0.1)]),
]
```

### 3. Run Sensitivity for Each Plan

```julia
results = Dict{String, SensitivityResult}()
for plan in plans
    results[plan.name] = run_sensitivity(spec, grid, solver; plan=plan)
end
```

### 4. Analyze Results

```julia
println("Parameter Sensitivity Summary:")
println("=" ^ 50)
for (name, res) in results
    println("$name:")
    println("  Max absolute delta: $(res.metrics.max_abs_delta)")
    println("  Max relative delta: $(res.metrics.max_rel_delta)")
    println()
end
```

---

## Serialization

Sensitivity results can be saved to JSON artifacts.

### Single Sensitivity Artifact

```julia
write_sensitivity_json(
    "sensitivity_result.json";
    model_spec=spec,
    grid=grid,
    solver=solver,
    result=result
)
```

### Population Sensitivity Artifact

```julia
write_population_sensitivity_json(
    "pop_sensitivity_result.json";
    population_spec=pop,
    grid=grid,
    solver=solver,
    result=result
)
```

### Replay Sensitivity Artifact

```julia
artifact = read_execution_json("sensitivity_result.json")
replayed = replay_sensitivity_execution(artifact)
```

---

## Complete Example: Parameter Ranking

```julia
using NeoPKPD

# Model setup
params = OneCompOralFirstOrderParams(1.5, 5.0, 50.0)
spec = ModelSpec(
    OneCompOralFirstOrder(),
    "oral_sensitivity",
    params,
    [DoseEvent(0.0, 200.0)]
)

grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Parameter perturbations (all +10%)
param_names = [:Ka, :CL, :V]
sensitivity_results = Dict{Symbol, Float64}()

for param in param_names
    plan = PerturbationPlan(
        "$(param)_sensitivity",
        [Perturbation(RelativePerturbation(), param, 0.1)]
    )
    result = run_sensitivity(spec, grid, solver; plan=plan)
    sensitivity_results[param] = result.metrics.max_rel_delta
end

# Rank parameters by sensitivity
sorted_params = sort(collect(sensitivity_results), by=x->x[2], rev=true)

println("Parameter Sensitivity Ranking (10% perturbation):")
println("=" ^ 50)
for (i, (param, sens)) in enumerate(sorted_params)
    println("$i. $param: max relative delta = $(round(sens, digits=4))")
end
```

---

## Best Practices

1. **Use Relative Perturbations** for comparing parameters with different scales
2. **Match Perturbation Size** to expected variability (e.g., IIV omega values)
3. **Consider Multiple Metrics** - max_abs for safety margins, max_rel for proportional impact
4. **Population Sensitivity** for understanding variability-parameter interactions
5. **Document Plans** with meaningful names for reproducibility
