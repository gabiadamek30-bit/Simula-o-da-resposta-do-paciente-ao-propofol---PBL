# Monolix Import

Comprehensive guide for importing Monolix project files (.mlxtran) into NeoPKPD.

---

## Overview

NeoPKPD can parse Monolix project files and convert them to native NeoPKPD model specifications, enabling migration from Monolix workflows.

```julia
using NeoPKPD

result = import_monolix("project.mlxtran")
println("Model: $(result.model_type)")
println("Parameters: $(result.parameters)")
```

---

## Quick Start

### Basic Import

```julia
using NeoPKPD

# Import Monolix project file
result = import_monolix("project.mlxtran")

# Access the converted model
model_spec = result.spec

# Simulate with the imported model
times = 0.0:0.5:24.0
dose = DoseEvent(time=0.0, amount=100.0)
sim = simulate_single(model_spec, collect(times), [dose])
```

### With Dose Events

```julia
# Specify doses if needed
doses = [
    DoseEvent(time=0.0, amount=100.0, route=:oral)
]

result = import_monolix("project.mlxtran"; doses=doses)
```

---

## Supported Model Types

### Pharmacokinetic Models

| Monolix Model | NeoPKPD Model |
|---------------|----------------|
| `pk_bolus1cpt_Vk_PLASMA` | `:OneCompIVBolus` |
| `pk_bolus1cpt_VCl_PLASMA` | `:OneCompIVBolus` |
| `pk_oral1cpt_kaVk_PLASMA` | `:OneCompOralFirstOrder` |
| `pk_oral1cpt_kaVCl_PLASMA` | `:OneCompOralFirstOrder` |
| `pk_bolus2cpt_V1k12k21k_PLASMA` | `:TwoCompIVBolus` |
| `pk_bolus2cpt_V1ClQ2V2_PLASMA` | `:TwoCompIVBolus` |
| `pk_oral2cpt_kaV1k12k21k_PLASMA` | `:TwoCompOral` |
| `pk_oral2cpt_kaV1ClQ2V2_PLASMA` | `:TwoCompOral` |
| `pk_bolus3cpt_V1ClQ2V2Q3V3_PLASMA` | `:ThreeCompIVBolus` |
| `pk_bolus1cpt_VVmKm_PLASMA` | `:MichaelisMentenElimination` |

### Model Name Pattern Recognition

```julia
# Monolix model names follow patterns:
# pk_{route}{n}cpt_{parameters}_{observation}

# Route: bolus, oral, infusion
# n: 1, 2, 3 (compartments)
# Parameters: V, Cl, k, ka, Q, Vm, Km
# Observation: PLASMA, EFFECT, etc.
```

---

## Monolix Project Structure

### Project File Format (.mlxtran)

```xml
<DATAFILE>
  <FILENAME path="../data/pk_data.csv"/>
  <HEADER>
    ID, TIME, DV, AMT, EVID, MDV, WT
  </HEADER>
</DATAFILE>

<STRUCTURAL_MODEL>
  <FILE path="lib:pk_oral1cpt_kaVCl_PLASMA"/>
</STRUCTURAL_MODEL>

<PARAMETER>
  <POPULATION name="ka_pop" value="1.5"/>
  <POPULATION name="V_pop" value="50.0"/>
  <POPULATION name="Cl_pop" value="5.0"/>
</PARAMETER>

<INDIVIDUAL>
  <PARAMETER name="ka" variability="none"/>
  <PARAMETER name="V" variability="lognormal" value="0.25"/>
  <PARAMETER name="Cl" variability="lognormal" value="0.30"/>
</INDIVIDUAL>

<OBSERVATION>
  <ERROR type="proportional" value="0.1"/>
</OBSERVATION>
```

---

## Parsed Structures

### MonolixProject

```julia
struct MonolixProject
    data::MonolixDataset              # Data file specification
    structural_model::MonolixStructuralModel  # Model definition
    parameters::Vector{MonolixParameter}      # Population parameters
    individual::Vector{MonolixIndividual}     # IIV specifications
    observation::MonolixObservation   # Error model
    estimation::Dict{String,Any}      # Estimation settings
    raw_text::String                  # Original file content
end
```

### MonolixParameter

```julia
struct MonolixParameter
    name::String           # Parameter name (ka_pop, V_pop, Cl_pop)
    value::Float64         # Population value
    fixed::Bool            # Fixed or estimated
    lower_bound::Float64   # Lower constraint
    upper_bound::Float64   # Upper constraint
end
```

### MonolixIndividual

```julia
struct MonolixIndividual
    parameter::String      # Base parameter name
    variability::Symbol    # :none, :lognormal, :normal, :logitnormal
    omega::Float64         # Standard deviation
    correlation::Dict{String,Float64}  # Correlations with other params
end
```

---

## Parameter Extraction

### Population Parameters

```julia
result = import_monolix("project.mlxtran")

# Access population parameters
for (param, value) in result.parameters
    println("$param = $value")
end

# Example output:
# Ka = 1.5
# V = 50.0
# CL = 5.0
```

### Random Effects

```julia
if !isnothing(result.iiv)
    println("IIV Structure:")
    for (i, param) in enumerate(result.iiv.parameters)
        omega = result.iiv.omega[i, i]
        cv = sqrt(exp(omega^2) - 1) * 100  # For lognormal
        println("  $param: ω=$(omega), CV≈$(round(cv, digits=1))%")
    end
end
```

### Parameter Transformations

| Monolix Variability | NeoPKPD Transformation | Formula |
|---------------------|-------------------------|---------|
| `lognormal` | `:exponential` | $\theta_i = \theta_{pop} \cdot e^{\eta_i}$ |
| `normal` | `:additive` | $\theta_i = \theta_{pop} + \eta_i$ |
| `logitnormal` | `:logit` | Logit-normal transform |
| `none` | No IIV | Fixed to population value |

---

## Error Model Import

### Supported Error Models

```julia
# Proportional error
# Y = F * (1 + b*ε)
struct ProportionalError
    b::Float64  # Proportional coefficient
end

# Additive error
# Y = F + a*ε
struct AdditiveError
    a::Float64  # Additive coefficient
end

# Combined error
# Y = F * (1 + b*ε₁) + a*ε₂
struct CombinedError
    a::Float64  # Additive
    b::Float64  # Proportional
end
```

### Accessing Error Model

```julia
result = import_monolix("project.mlxtran")

if !isnothing(result.error)
    println("Error type: $(result.error.type)")
    println("Error SD: $(result.error.sigma)")
end
```

---

## Covariate Models

### Supported Covariate Patterns

```xml
<COVARIATE>
  <CONTINUOUS name="WT" transformation="none"/>
  <CONTINUOUS name="AGE" transformation="log"/>
  <CATEGORICAL name="SEX" categories="M,F"/>
</COVARIATE>

<INDIVIDUAL>
  <PARAMETER name="Cl">
    <COVARIATE name="WT" coefficient="0.75" type="power" reference="70"/>
  </PARAMETER>
</INDIVIDUAL>
```

### Extracting Covariate Effects

```julia
result = import_monolix("project.mlxtran")

for effect in result.covariate_effects
    println("$(effect.covariate) on $(effect.parameter):")
    println("  Type: $(effect.effect_type)")
    println("  Coefficient: $(effect.coefficient)")
    println("  Reference: $(effect.reference_value)")
end
```

---

## Unsupported Features

### Models Not Supported

| Model Type | Reason | Workaround |
|------------|--------|------------|
| PD models | Complex dynamics | Manual definition |
| Turnover models | Indirect response | Use NeoPKPD IRM |
| Transit compartment | Absorption | Use NeoPKPD transit |
| Mixture models | Subpopulations | Not supported |
| Markov models | State transitions | Not supported |
| Time-to-event | Survival | Not supported |
| Count data | Poisson/NegBin | Not supported |
| Categorical | Ordered/unordered | Not supported |

### Features Imported with Warnings

| Feature | Handling | Warning |
|---------|----------|---------|
| Lag time (Tlag) | Ignored | "Lag time not imported, assuming Tlag=0" |
| Bioavailability (F) | Assumes F=1 | "Bioavailability not imported, assuming F=1" |
| Complex covariates | Partial | "Complex covariate effect simplified" |

---

## Validation

### Checking Import Quality

```julia
result = import_monolix("project.mlxtran")

# Check for warnings
if !isempty(result.warnings)
    println("⚠️ Import warnings:")
    for w in result.warnings
        println("  - $w")
    end
end

# Validate model type was recognized
if result.model_type == :Unknown
    println("❌ Model type not recognized")
else
    println("✓ Model type: $(result.model_type)")
end

# Check parameter completeness
expected_params = [:Ka, :CL, :V]  # For 1-comp oral
for param in expected_params
    if haskey(result.parameters, param)
        println("✓ $param = $(result.parameters[param])")
    else
        println("❌ Missing parameter: $param")
    end
end
```

---

## Complete Example

### Monolix Project File (project.mlxtran)

```xml
<?xml version="1.0" encoding="UTF-8"?>
<monolix>
  <project name="pk_analysis" version="2023R1">

    <DATAFILE>
      <FILENAME path="../data/pk_data.csv"/>
      <HEADER>ID, TIME, DV, AMT, EVID, WT</HEADER>
      <COLUMNMAPPING>
        <COLUMN name="ID" type="ID"/>
        <COLUMN name="TIME" type="TIME"/>
        <COLUMN name="DV" type="OBSERVATION"/>
        <COLUMN name="AMT" type="AMOUNT"/>
        <COLUMN name="EVID" type="EVID"/>
        <COLUMN name="WT" type="COVARIATE"/>
      </COLUMNMAPPING>
    </DATAFILE>

    <STRUCTURAL_MODEL>
      <FILE path="lib:pk_oral2cpt_kaV1ClQ2V2_PLASMA"/>
    </STRUCTURAL_MODEL>

    <PARAMETER>
      <POPULATION name="ka_pop" value="1.5" method="MLE"/>
      <POPULATION name="V1_pop" value="50.0" method="MLE"/>
      <POPULATION name="Cl_pop" value="5.0" method="MLE"/>
      <POPULATION name="Q_pop" value="3.0" method="MLE"/>
      <POPULATION name="V2_pop" value="80.0" method="MLE"/>
    </PARAMETER>

    <INDIVIDUAL>
      <PARAMETER name="ka" variability="lognormal">
        <OMEGA value="0.3"/>
      </PARAMETER>
      <PARAMETER name="V1" variability="lognormal">
        <OMEGA value="0.25"/>
      </PARAMETER>
      <PARAMETER name="Cl" variability="lognormal">
        <OMEGA value="0.30"/>
        <COVARIATE name="WT" coefficient="0.75" type="power" reference="70"/>
      </PARAMETER>
      <PARAMETER name="Q" variability="none"/>
      <PARAMETER name="V2" variability="none"/>
    </INDIVIDUAL>

    <OBSERVATION>
      <PREDICTION name="Cc"/>
      <ERROR type="proportional">
        <PARAMETER name="b" value="0.1"/>
      </ERROR>
    </OBSERVATION>

    <ESTIMATION>
      <METHOD name="SAEM"/>
      <NBCHAINS value="5"/>
      <NBITERATIONS value="500"/>
    </ESTIMATION>

  </project>
</monolix>
```

### Julia Import Code

```julia
using NeoPKPD

# Import Monolix project
result = import_monolix("project.mlxtran")

# Display results
println("=" ^ 50)
println("Monolix Import Results")
println("=" ^ 50)

println("\n--- Model Information ---")
println("Source: Monolix")
println("Model type: $(result.model_type)")

println("\n--- Population Parameters ---")
for (param, value) in result.parameters
    println("  $param = $value")
end

println("\n--- Inter-Individual Variability ---")
if !isnothing(result.iiv)
    println("  Parameters with IIV: $(result.iiv.parameters)")
    println("  Transformations: $(result.iiv.transformations)")
    println("\n  Omega matrix (SD scale):")
    for (i, p) in enumerate(result.iiv.parameters)
        omega = result.iiv.omega[i, i]
        println("    ω_$p = $omega")
    end
end

println("\n--- Residual Error ---")
if !isnothing(result.error)
    println("  Type: $(result.error.type)")
    println("  Coefficient: $(result.error.sigma)")
end

println("\n--- Covariate Effects ---")
if !isempty(result.covariate_effects)
    for eff in result.covariate_effects
        println("  $(eff.covariate) on $(eff.parameter):")
        println("    Effect: $(eff.effect_type)")
        println("    Coefficient: $(eff.coefficient)")
    end
else
    println("  None imported")
end

println("\n--- Warnings ---")
if isempty(result.warnings)
    println("  None")
else
    for w in result.warnings
        println("  ⚠️ $w")
    end
end

# Validate with simulation
println("\n--- Validation Simulation ---")
times = collect(0.0:0.5:72.0)
doses = [DoseEvent(time=0.0, amount=500.0, route=:oral)]

sim = simulate_single(result.spec, times, doses)
println("Simulation completed successfully")
println("Cmax: $(round(maximum(sim.observations[:conc]), digits=2))")
println("Tmax: $(times[argmax(sim.observations[:conc])])")
```

### Expected Output

```
==================================================
Monolix Import Results
==================================================

--- Model Information ---
Source: Monolix
Model type: TwoCompOral

--- Population Parameters ---
  Ka = 1.5
  V1 = 50.0
  CL = 5.0
  Q = 3.0
  V2 = 80.0

--- Inter-Individual Variability ---
  Parameters with IIV: [:Ka, :V1, :CL]
  Transformations: [:exponential, :exponential, :exponential]

  Omega matrix (SD scale):
    ω_Ka = 0.3
    ω_V1 = 0.25
    ω_CL = 0.3

--- Residual Error ---
  Type: proportional
  Coefficient: 0.1

--- Covariate Effects ---
  WT on CL:
    Effect: power
    Coefficient: 0.75

--- Warnings ---
  None

--- Validation Simulation ---
Simulation completed successfully
Cmax: 4.87
Tmax: 1.5
```

---

## Model Library Reference

### One-Compartment Models

| Library Model | Parameters | Route |
|---------------|------------|-------|
| `pk_bolus1cpt_Vk_PLASMA` | V, k | IV bolus |
| `pk_bolus1cpt_VCl_PLASMA` | V, Cl | IV bolus |
| `pk_oral1cpt_kaVk_PLASMA` | ka, V, k | Oral |
| `pk_oral1cpt_kaVCl_PLASMA` | ka, V, Cl | Oral |
| `pk_infusion1cpt_VCl_PLASMA` | V, Cl | IV infusion |

### Two-Compartment Models

| Library Model | Parameters | Route |
|---------------|------------|-------|
| `pk_bolus2cpt_V1k12k21k_PLASMA` | V1, k, k12, k21 | IV bolus |
| `pk_bolus2cpt_V1ClQ2V2_PLASMA` | V1, Cl, Q, V2 | IV bolus |
| `pk_oral2cpt_kaV1k12k21k_PLASMA` | ka, V1, k, k12, k21 | Oral |
| `pk_oral2cpt_kaV1ClQ2V2_PLASMA` | ka, V1, Cl, Q, V2 | Oral |

### Three-Compartment Models

| Library Model | Parameters | Route |
|---------------|------------|-------|
| `pk_bolus3cpt_V1ClQ2V2Q3V3_PLASMA` | V1, Cl, Q2, V2, Q3, V3 | IV bolus |

### Special Models

| Library Model | Parameters | Description |
|---------------|------------|-------------|
| `pk_bolus1cpt_VVmKm_PLASMA` | V, Vm, Km | Michaelis-Menten |

---

## See Also

- [NONMEM Import](nonmem.md) - Import NONMEM control files
- [CDISC Import](cdisc.md) - Import CDISC data
- [Model Specification](../models/index.md) - NeoPKPD model types

