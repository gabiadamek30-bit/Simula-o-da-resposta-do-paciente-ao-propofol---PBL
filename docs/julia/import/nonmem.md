# NONMEM Import

Comprehensive guide for importing NONMEM control stream files (.ctl, .mod) into NeoPKPD.

---

## Overview

NeoPKPD can parse NONMEM control files and convert them to native NeoPKPD model specifications, enabling seamless migration from NONMEM workflows.

```julia
using NeoPKPD

result = import_nonmem("run001.ctl")
println("Model: $(result.model_type)")
println("Parameters: $(result.parameters)")
```

---

## Quick Start

### Basic Import

```julia
using NeoPKPD

# Import NONMEM control file
result = import_nonmem("run001.ctl")

# Access the converted model
model_spec = result.spec

# Simulate with the imported model
times = 0.0:0.5:24.0
dose = DoseEvent(time=0.0, amount=100.0)
sim = simulate_single(model_spec, collect(times), [dose])
```

### With Dose Events

```julia
# Specify doses if not in control file
doses = [
    DoseEvent(time=0.0, amount=100.0),
    DoseEvent(time=12.0, amount=100.0)
]

result = import_nonmem("run001.ctl"; doses=doses)
```

---

## Import Result Structure

```julia
struct ImportResult
    model_type::Symbol              # NeoPKPD model kind
    parameters::Dict{Symbol,Float64} # Parameter values
    spec::ModelSpec                 # Complete model specification
    iiv::Union{Nothing,IIVSpec}     # Inter-individual variability
    error::Union{Nothing,ResidualErrorSpec}  # Residual error model
    source::String                  # "NONMEM"
    warnings::Vector{String}        # Import warnings
end
```

### Accessing Results

```julia
result = import_nonmem("run001.ctl")

# Model type (e.g., :OneCompOralFirstOrder)
println("Model: $(result.model_type)")

# Fixed effect parameters
for (param, value) in result.parameters
    println("  $param = $value")
end

# IIV specification
if !isnothing(result.iiv)
    println("IIV on: $(result.iiv.parameters)")
    println("Omega: $(result.iiv.omega)")
end

# Residual error
if !isnothing(result.error)
    println("Error type: $(result.error.type)")
    println("Error param: $(result.error.sigma)")
end

# Check for warnings
for warning in result.warnings
    println("⚠️ $warning")
end
```

---

## Supported ADVAN/TRANS Combinations

### ADVAN1 - One-Compartment IV Bolus

```
$SUBROUTINES ADVAN1 TRANS2
```

| TRANS | Parameters | NeoPKPD Model |
|-------|------------|----------------|
| TRANS1 | K, V | `:OneCompIVBolus` |
| TRANS2 | CL, V | `:OneCompIVBolus` |

### ADVAN2 - One-Compartment Oral

```
$SUBROUTINES ADVAN2 TRANS2
```

| TRANS | Parameters | NeoPKPD Model |
|-------|------------|----------------|
| TRANS1 | KA, K, V | `:OneCompOralFirstOrder` |
| TRANS2 | KA, CL, V | `:OneCompOralFirstOrder` |

### ADVAN3 - Two-Compartment IV Bolus

```
$SUBROUTINES ADVAN3 TRANS4
```

| TRANS | Parameters | NeoPKPD Model |
|-------|------------|----------------|
| TRANS1 | K, K12, K21, V | `:TwoCompIVBolus` |
| TRANS3 | CL, V, Q, VSS | `:TwoCompIVBolus` |
| TRANS4 | CL, V1, Q, V2 | `:TwoCompIVBolus` |

### ADVAN4 - Two-Compartment Oral

```
$SUBROUTINES ADVAN4 TRANS4
```

| TRANS | Parameters | NeoPKPD Model |
|-------|------------|----------------|
| TRANS1 | KA, K, K23, K32, V | `:TwoCompOral` |
| TRANS3 | KA, CL, V, Q, VSS | `:TwoCompOral` |
| TRANS4 | KA, CL, V1, Q, V2 | `:TwoCompOral` |

### ADVAN10 - Michaelis-Menten Elimination

```
$SUBROUTINES ADVAN10
```

| Parameters | NeoPKPD Model |
|------------|----------------|
| VM, KM, V | `:MichaelisMentenElimination` |

### ADVAN11 - Three-Compartment IV Bolus

```
$SUBROUTINES ADVAN11 TRANS4
```

| TRANS | Parameters | NeoPKPD Model |
|-------|------------|----------------|
| TRANS1 | K, K12, K21, K13, K31, V | `:ThreeCompIVBolus` |
| TRANS4 | CL, V1, Q2, V2, Q3, V3 | `:ThreeCompIVBolus` |

---

## $THETA Parsing

### Supported Formats

```fortran
; Simple initial estimate
$THETA 5.0

; With bounds
$THETA (0, 5.0, 100)

; Lower bound only
$THETA (0, 5.0)

; Fixed parameter
$THETA 5.0 FIX

; Multiple on one line
$THETA 5.0 50.0 0.1

; With comments
$THETA 5.0 ; CL (L/h)
```

### Accessing THETA Values

```julia
result = import_nonmem("run001.ctl")

# Parsed control file contains THETA specs
parsed = result.parsed_control_file
for (i, theta) in enumerate(parsed.thetas)
    println("THETA($i):")
    println("  Initial: $(theta.init)")
    println("  Lower: $(theta.lower)")
    println("  Upper: $(theta.upper)")
    println("  Fixed: $(theta.fixed)")
end
```

---

## $OMEGA Parsing

### DIAGONAL Structure

```fortran
$OMEGA
0.09       ; IIV on CL (30% CV)
0.0625     ; IIV on V (25% CV)
```

### BLOCK Structure

```fortran
$OMEGA BLOCK(2)
0.09                    ; Variance CL
0.03 0.0625             ; Covariance, Variance V
```

### Mixed Structures

```fortran
$OMEGA BLOCK(2)
0.09
0.03 0.0625

$OMEGA
0.04      ; Separate IIV on Ka
```

### Conversion to IIV

```julia
result = import_nonmem("run001.ctl")

if !isnothing(result.iiv)
    # Variances converted to standard deviations
    println("IIV parameters: $(result.iiv.parameters)")
    println("Omega matrix (SD): $(result.iiv.omega)")
    println("Transformations: $(result.iiv.transformations)")
end
```

---

## $SIGMA Parsing

### Proportional Error

```fortran
$SIGMA
0.01      ; Proportional error (10% CV)
```

### Additive Error

```fortran
$SIGMA
0.1       ; Additive error (SD = 0.316)
```

### Combined Error

```fortran
$SIGMA BLOCK(2)
0.01      ; Proportional
0.0 0.1   ; Additive
```

---

## $PK Block Parsing

### Supported Patterns

#### Typical Value Definitions

```fortran
$PK
TVCL = THETA(1)
TVV = THETA(2)
TVKA = THETA(3)
```

#### Parameter Assignments with ETA

```fortran
; Exponential IIV (default)
CL = TVCL * EXP(ETA(1))
V = TVV * EXP(ETA(2))

; Additive IIV
CL = TVCL + ETA(1)

; Proportional IIV
CL = TVCL * (1 + ETA(1))
```

#### Covariate Effects

```fortran
; Power model (weight effect)
TVCL = THETA(1) * (WT/70)**THETA(4)

; Linear model (age effect)
TVCL = THETA(1) * (1 + THETA(5)*(AGE-40))

; Exponential model
TVCL = THETA(1) * EXP(THETA(6)*(CRCL-100))
```

#### Scaling Factors

```fortran
; Volume scaling for central compartment
S1 = V
S2 = V1
```

### Parsed PK Block Structure

```julia
struct PKBlock
    tv_definitions::Dict{Symbol,Int}    # TV to THETA mapping
    assignments::Vector{PKAssignment}   # Parameter assignments
    scaling::Vector{ScalingFactor}      # Scaling factors
    unsupported_lines::Vector{String}   # Lines not parsed
end

struct PKAssignment
    parameter::Symbol         # CL, V, Ka, etc.
    tv_symbol::Symbol         # TVCL, TVV, etc.
    eta_index::Union{Nothing,Int}  # ETA index or nothing
    transformation::Symbol    # :exponential, :additive, :proportional
end
```

### Accessing Parsed $PK

```julia
result = import_nonmem("run001.ctl")
pk = result.parsed_control_file.pk_block

# TV definitions
for (tv, theta_idx) in pk.tv_definitions
    println("$tv = THETA($theta_idx)")
end

# Parameter assignments
for assign in pk.assignments
    println("$(assign.parameter): TV=$(assign.tv_symbol), ETA=$(assign.eta_index)")
end

# Check for unsupported lines
if !isempty(pk.unsupported_lines)
    println("⚠️ Unsupported $PK lines:")
    for line in pk.unsupported_lines
        println("  $line")
    end
end
```

---

## $ERROR Block Parsing

### Proportional Error

```fortran
$ERROR
IPRED = F
W = IPRED * THETA(4)
Y = IPRED + W * ERR(1)
```

### Additive Error

```fortran
$ERROR
IPRED = F
W = THETA(4)
Y = IPRED + W * ERR(1)
```

### Combined Error

```fortran
$ERROR
IPRED = F
W = SQRT(THETA(4)**2 + (THETA(5)*IPRED)**2)
Y = IPRED + W * ERR(1)
```

### Exponential Error

```fortran
$ERROR
IPRED = F
Y = IPRED * EXP(ERR(1))
```

### Detection Logic

```julia
result = import_nonmem("run001.ctl")

if !isnothing(result.error)
    println("Error type: $(result.error.type)")  # :proportional, :additive, :combined, :exponential
    println("Sigma value: $(result.error.sigma)")
end
```

---

## Covariate Effects Extraction

### Supported Covariate Patterns

```julia
struct PKCovariateEffect
    parameter::Symbol           # Affected parameter (CL, V, etc.)
    covariate::Symbol           # Covariate name (WT, AGE, CRCL)
    effect_type::Symbol         # :power, :linear, :exponential
    theta_index::Int            # THETA index for effect
    reference_value::Float64    # Centering value (70 for WT, etc.)
end
```

### Example: Weight on CL

```fortran
; Power model
TVCL = THETA(1) * (WT/70)**THETA(4)
```

Extracted as:
```julia
PKCovariateEffect(
    parameter = :CL,
    covariate = :WT,
    effect_type = :power,
    theta_index = 4,
    reference_value = 70.0
)
```

### Accessing Covariate Effects

```julia
result = import_nonmem("run001.ctl")

for effect in result.covariate_effects
    println("$(effect.covariate) on $(effect.parameter):")
    println("  Type: $(effect.effect_type)")
    println("  THETA: $(effect.theta_index)")
    println("  Reference: $(effect.reference_value)")
end
```

---

## Validation and Quality Checks

### Automatic Validation

The import process performs these checks:

1. **ADVAN/TRANS compatibility** - Unsupported combinations flagged
2. **THETA index bounds** - References beyond defined THETAs
3. **ETA index bounds** - References beyond OMEGA dimension
4. **Scaling factor validity** - S1, S2 must reference valid compartments
5. **Unused parameters** - THETAs defined but not used

### Conversion Result

```julia
struct NONMEMConversionResult
    model_spec::Union{Nothing,ModelSpec}
    iiv_spec::Union{Nothing,IIVSpec}
    error_spec::Union{Nothing,ResidualErrorSpec}
    warnings::Vector{String}
    errors::Vector{String}
    parameter_mapping::Dict{Int,Symbol}
    covariate_effects::Vector{PKCovariateEffect}
end
```

### Handling Warnings

```julia
result = import_nonmem("run001.ctl")

if !isempty(result.warnings)
    println("Import completed with warnings:")
    for w in result.warnings
        println("  ⚠️ $w")
    end
end

# Check for critical errors
if isnothing(result.spec)
    println("Import failed - check errors")
end
```

---

## Unsupported Features

### Not Currently Supported

| Feature | Description | Workaround |
|---------|-------------|------------|
| Custom $DES | User-defined ODEs | Manual model definition |
| Complex IF | Conditional logic | Simplify before import |
| ALAG | Absorption lag | Manual specification |
| F1, F2 | Bioavailability | Assume F=1 or manual |
| R1, D1 | Infusion params | Specify in dose events |
| MTIME | Model event time | Not supported |
| $MIX | Mixture models | Manual definition |
| $MODEL | Custom compartments | Use predefined models |

### Detection of Unsupported Constructs

```julia
result = import_nonmem("run001.ctl")

# Warnings include unsupported feature notices
for warning in result.warnings
    if contains(warning, "unsupported")
        println("Feature not imported: $warning")
    end
end
```

---

## Complete Example

### NONMEM Control File (run001.ctl)

```fortran
$PROBLEM Two-compartment oral PK model
$DATA ../data/pk_data.csv IGNORE=@
$INPUT ID TIME DV AMT EVID MDV WT AGE

$SUBROUTINES ADVAN4 TRANS4

$PK
; Typical values with weight covariate on CL
TVCL = THETA(1) * (WT/70)**0.75
TVV1 = THETA(2)
TVQ = THETA(3)
TVV2 = THETA(4)
TVKA = THETA(5)

; Individual parameters
CL = TVCL * EXP(ETA(1))
V1 = TVV1 * EXP(ETA(2))
Q = TVQ
V2 = TVV2
KA = TVKA * EXP(ETA(3))

S2 = V1

$ERROR
IPRED = F
W = IPRED * THETA(6)
Y = IPRED + W * ERR(1)

$THETA
(0, 10.0)     ; CL (L/h)
(0, 50.0)     ; V1 (L)
(0, 5.0)      ; Q (L/h)
(0, 100.0)    ; V2 (L)
(0, 1.5)      ; Ka (1/h)
(0, 0.1)      ; Proportional error

$OMEGA BLOCK(2)
0.09          ; IIV CL
0.03 0.0625   ; IIV V1

$OMEGA
0.04          ; IIV Ka

$SIGMA
1 FIX         ; SIGMA fixed to 1 (error in THETA(6))

$ESTIMATION METHOD=1 INTER MAXEVAL=9999
$COV PRINT=E
```

### Julia Import Code

```julia
using NeoPKPD

# Import the control file
result = import_nonmem("run001.ctl")

# Display import results
println("=" ^ 50)
println("NONMEM Import Results")
println("=" ^ 50)

println("\n--- Model Type ---")
println("NeoPKPD model: $(result.model_type)")

println("\n--- Fixed Effects (THETA) ---")
for (param, value) in result.parameters
    println("  $param = $value")
end

println("\n--- Random Effects (OMEGA) ---")
if !isnothing(result.iiv)
    println("  Parameters: $(result.iiv.parameters)")
    println("  Transformations: $(result.iiv.transformations)")
    println("  Omega matrix:")
    display(result.iiv.omega)
end

println("\n--- Residual Error (SIGMA) ---")
if !isnothing(result.error)
    println("  Type: $(result.error.type)")
    println("  Value: $(result.error.sigma)")
end

println("\n--- Covariate Effects ---")
for effect in result.covariate_effects
    println("  $(effect.covariate) on $(effect.parameter): $(effect.effect_type)")
end

println("\n--- Warnings ---")
if isempty(result.warnings)
    println("  None")
else
    for w in result.warnings
        println("  ⚠️ $w")
    end
end

# Simulate with imported model
println("\n--- Validation Simulation ---")
times = collect(0.0:0.5:48.0)
doses = [DoseEvent(time=0.0, amount=500.0)]

sim = simulate_single(result.spec, times, doses)
println("Cmax: $(maximum(sim.observations[:conc]))")
println("Tmax: $(times[argmax(sim.observations[:conc])])")
```

### Expected Output

```
==================================================
NONMEM Import Results
==================================================

--- Model Type ---
NeoPKPD model: TwoCompOral

--- Fixed Effects (THETA) ---
  CL = 10.0
  V1 = 50.0
  Q = 5.0
  V2 = 100.0
  Ka = 1.5

--- Random Effects (OMEGA) ---
  Parameters: [:CL, :V1, :Ka]
  Transformations: [:exponential, :exponential, :exponential]
  Omega matrix:
3×3 Matrix{Float64}:
 0.3   0.173  0.0
 0.173 0.25   0.0
 0.0   0.0    0.2

--- Residual Error (SIGMA) ---
  Type: proportional
  Value: 0.1

--- Covariate Effects ---
  WT on CL: power

--- Warnings ---
  None

--- Validation Simulation ---
Cmax: 4.23
Tmax: 1.5
```

---

## CLI Usage

### Basic Import

```bash
./bin/neopkpd import --input run001.ctl --format nonmem
```

### Export to JSON

```bash
./bin/neopkpd import --input run001.ctl --format nonmem --out model.json
```

### Validate Import

```bash
./bin/neopkpd import --input run001.ctl --format nonmem --validate
```

---

## See Also

- [Monolix Import](monolix.md) - Import Monolix projects
- [CDISC Import](cdisc.md) - Import CDISC data
- [Model Specification](../models/index.md) - NeoPKPD model types

