# NONMEM Model Import

Import NONMEM control stream files (.ctl, .mod) into NeoPKPD Python.

---

## Overview

```python
from neopkpd.import_ import import_nonmem

model = import_nonmem("run001.ctl")
print(f"Model: {model.model_kind}")
print(f"Parameters: {model.params}")
```

---

## Quick Start

### Basic Import

```python
from neopkpd.import_ import import_nonmem

# Import NONMEM control file
model = import_nonmem("run001.ctl")

# Access model information
print(f"Model type: {model.model_kind}")
print(f"Source: {model.source_format}")

# Access parameters
for name, value in model.params.items():
    print(f"  {name} = {value}")
```

### With Dose Events

```python
# Specify doses if not in control file
doses = [
    {"time": 0.0, "amount": 100.0},
    {"time": 12.0, "amount": 100.0}
]

model = import_nonmem("run001.ctl", doses=doses)
```

---

## ImportedModel Class

```python
@dataclass
class ImportedModel:
    """Result of importing a model from external format."""

    # Source information
    source_format: str          # "nonmem" or "monolix"
    source_file: str            # Path to source file

    # Model type
    model_kind: str             # NeoPKPD model name

    # Fixed effects
    params: dict[str, float]    # Parameter name → value
    theta_init: list[float]     # Initial THETA values
    theta_names: list[str]      # Parameter names

    # Random effects
    omega_init: list[list[float]]  # Omega matrix
    omega_names: list[str]         # Random effect names

    # Residual error
    sigma_type: str             # "proportional", "additive", etc.
    sigma_init: float           # Error parameter value

    # Quality
    warnings: list[str]         # Import warnings

    # Metadata
    metadata: dict[str, Any]    # Additional information
```

### Accessing Model Details

```python
model = import_nonmem("run001.ctl")

# Fixed effects
print("Fixed Effects:")
for name, value in model.params.items():
    print(f"  {name} = {value}")

# Random effects
print("\nRandom Effects:")
for i, name in enumerate(model.omega_names):
    omega = model.omega_init[i][i]
    cv = (math.exp(omega) - 1) ** 0.5 * 100  # CV for lognormal
    print(f"  {name}: ω² = {omega:.4f} (CV ≈ {cv:.1f}%)")

# Residual error
print(f"\nError Model: {model.sigma_type}")
print(f"Sigma: {model.sigma_init}")

# Warnings
if model.warnings:
    print("\nWarnings:")
    for w in model.warnings:
        print(f"  ⚠️ {w}")
```

---

## Supported ADVAN/TRANS

### One-Compartment Models

| ADVAN | TRANS | Parameters | Model |
|-------|-------|------------|-------|
| ADVAN1 | TRANS1 | K, V | `OneCompIVBolus` |
| ADVAN1 | TRANS2 | CL, V | `OneCompIVBolus` |
| ADVAN2 | TRANS1 | KA, K, V | `OneCompOralFirstOrder` |
| ADVAN2 | TRANS2 | KA, CL, V | `OneCompOralFirstOrder` |

### Two-Compartment Models

| ADVAN | TRANS | Parameters | Model |
|-------|-------|------------|-------|
| ADVAN3 | TRANS1 | K, K12, K21, V | `TwoCompIVBolus` |
| ADVAN3 | TRANS4 | CL, V1, Q, V2 | `TwoCompIVBolus` |
| ADVAN4 | TRANS1 | KA, K, K23, K32, V | `TwoCompOral` |
| ADVAN4 | TRANS4 | KA, CL, V1, Q, V2 | `TwoCompOral` |

### Three-Compartment Models

| ADVAN | TRANS | Parameters | Model |
|-------|-------|------------|-------|
| ADVAN11 | TRANS4 | CL, V1, Q2, V2, Q3, V3 | `ThreeCompIVBolus` |

### Special Models

| ADVAN | Parameters | Model |
|-------|------------|-------|
| ADVAN10 | VM, KM, V | `MichaelisMentenElimination` |

---

## Parsed Control File

### Accessing Raw Parsed Data

```python
model = import_nonmem("run001.ctl")

# Access parsed control file structure
parsed = model.metadata.get("parsed_control_file")

if parsed:
    # Problem statement
    print(f"Problem: {parsed.problem}")

    # Subroutines
    print(f"ADVAN: {parsed.subroutines.advan}")
    print(f"TRANS: {parsed.subroutines.trans}")

    # THETA specifications
    for i, theta in enumerate(parsed.thetas):
        print(f"THETA({i+1}): init={theta.init}, bounds=({theta.lower}, {theta.upper})")
```

### THETASpec

```python
@dataclass
class THETASpec:
    init: float           # Initial estimate
    lower: float          # Lower bound
    upper: float          # Upper bound
    fixed: bool           # Fixed parameter
    name: str             # Optional name
```

### OMEGABlock

```python
@dataclass
class OMEGABlock:
    values: list[float]          # Variance/covariance values
    structure: OmegaStructure    # DIAGONAL or BLOCK
    dimension: int               # Block dimension
    fixed: bool                  # Fixed block
```

---

## $PK Block Parsing

### Supported Patterns

```python
# Typical value definitions
# TVCL = THETA(1)

# Parameter assignments with ETA
# CL = TVCL * EXP(ETA(1))     → exponential
# CL = TVCL + ETA(1)          → additive
# CL = TVCL * (1 + ETA(1))    → proportional

# Covariate effects
# TVCL = THETA(1) * (WT/70)**THETA(4)      → power
# TVCL = THETA(1) * (1 + THETA(5)*(AGE-40)) → linear
```

### Accessing PK Block

```python
model = import_nonmem("run001.ctl")
parsed = model.metadata.get("parsed_control_file")

if parsed and parsed.pk_block:
    pk = parsed.pk_block

    # TV definitions
    for tv, theta_idx in pk.tv_definitions.items():
        print(f"{tv} = THETA({theta_idx})")

    # Parameter assignments
    for assign in pk.assignments:
        print(f"{assign.parameter}: TV={assign.tv_symbol}, ETA={assign.eta_index}")

    # Unsupported lines
    if pk.unsupported_lines:
        print("Unsupported $PK lines:")
        for line in pk.unsupported_lines:
            print(f"  {line}")
```

---

## $ERROR Block Parsing

### Supported Error Models

| Type | Pattern | Formula |
|------|---------|---------|
| `proportional` | `W = IPRED * THETA(n)` | Y = F · (1 + ε) |
| `additive` | `W = THETA(n)` | Y = F + ε |
| `combined` | `W = SQRT(...)` | Y = F · (1 + ε₁) + ε₂ |
| `exponential` | `Y = F * EXP(ERR(1))` | Y = F · e^ε |

### Accessing Error Model

```python
model = import_nonmem("run001.ctl")

print(f"Error type: {model.sigma_type}")
print(f"Sigma value: {model.sigma_init}")

# Detailed error block
parsed = model.metadata.get("parsed_control_file")
if parsed and parsed.error_block:
    error = parsed.error_block
    print(f"THETA indices in error: {error.theta_indices}")
    print(f"SIGMA fixed to 1: {error.sigma_fixed_to_1}")
```

---

## Covariate Effects

### Extraction

```python
model = import_nonmem("run001.ctl")

# Covariate effects extracted from $PK
effects = model.metadata.get("covariate_effects", [])

for effect in effects:
    print(f"{effect['covariate']} on {effect['parameter']}:")
    print(f"  Type: {effect['type']}")
    print(f"  THETA: {effect['theta_index']}")
    print(f"  Reference: {effect['reference_value']}")
```

### Supported Effect Types

| Type | NONMEM Pattern | Example |
|------|----------------|---------|
| `power` | `(COV/REF)**THETA(n)` | `(WT/70)**0.75` |
| `linear` | `(1 + THETA(n)*(COV-REF))` | `(1 + 0.01*(AGE-40))` |
| `exponential` | `EXP(THETA(n)*(COV-REF))` | `EXP(0.005*(CRCL-100))` |

---

## Using Imported Models

### Simulation

```python
from neopkpd.import_ import import_nonmem
from neopkpd import simulate

# Import model
model = import_nonmem("run001.ctl")

# Create model spec for simulation
# (Automatic conversion)
times = list(range(0, 49))  # 0-48 hours
doses = [{"time": 0.0, "amount": 500.0, "route": "oral"}]

result = simulate(
    model_kind=model.model_kind,
    params=model.params,
    times=times,
    doses=doses
)

print(f"Cmax: {max(result.concentrations):.2f}")
```

### Population Simulation

```python
from neopkpd.import_ import import_nonmem
from neopkpd import simulate_population

model = import_nonmem("run001.ctl")

# Simulate population with imported IIV
pop_result = simulate_population(
    model_kind=model.model_kind,
    params=model.params,
    omega=model.omega_init,
    sigma=model.sigma_init,
    sigma_type=model.sigma_type,
    n_subjects=100,
    times=list(range(0, 49)),
    doses=[{"time": 0.0, "amount": 500.0}]
)
```

---

## Validation

### Check Import Quality

```python
model = import_nonmem("run001.ctl")

# Check for warnings
if model.warnings:
    print("Import warnings:")
    for w in model.warnings:
        print(f"  ⚠️ {w}")
else:
    print("✓ No warnings")

# Verify model type
if model.model_kind == "Unknown":
    print("❌ Model type not recognized")
else:
    print(f"✓ Model type: {model.model_kind}")

# Check parameters
print(f"✓ {len(model.params)} parameters imported")
print(f"✓ {len(model.omega_names)} random effects imported")
```

---

## Unsupported Features

### Not Supported

| Feature | Workaround |
|---------|------------|
| Custom $DES | Manual model definition |
| Complex IF statements | Simplify before import |
| ALAG (lag time) | Manual specification |
| F1, F2 (bioavailability) | Assumes F=1 |
| R1, D1 (infusion params) | Specify in dose events |
| MTIME | Not supported |
| $MIX (mixture models) | Not supported |
| $MODEL (custom compartments) | Use predefined models |

### Detection

```python
model = import_nonmem("run001.ctl")

# Warnings include unsupported features
for w in model.warnings:
    if "unsupported" in w.lower():
        print(f"Unsupported: {w}")
```

---

## Complete Example

```python
from neopkpd.import_ import import_nonmem
from neopkpd import simulate
import math

# Import NONMEM control file
print("Importing NONMEM control file...")
model = import_nonmem("run001.ctl")

# Display results
print("=" * 50)
print("NONMEM Import Results")
print("=" * 50)

print(f"\nSource: {model.source_file}")
print(f"Model type: {model.model_kind}")

print("\n--- Fixed Effects ---")
for name, value in model.params.items():
    print(f"  {name} = {value}")

print("\n--- Random Effects ---")
for i, name in enumerate(model.omega_names):
    omega_sq = model.omega_init[i][i]
    cv = (math.exp(omega_sq) - 1) ** 0.5 * 100
    print(f"  {name}: ω² = {omega_sq:.4f} (CV ≈ {cv:.1f}%)")

print(f"\n--- Residual Error ---")
print(f"  Type: {model.sigma_type}")
print(f"  Value: {model.sigma_init}")

print("\n--- Warnings ---")
if model.warnings:
    for w in model.warnings:
        print(f"  ⚠️ {w}")
else:
    print("  None")

# Validate with simulation
print("\n--- Validation Simulation ---")
times = list(range(0, 49, 1))
doses = [{"time": 0.0, "amount": 500.0, "route": "oral"}]

result = simulate(
    model_kind=model.model_kind,
    params=model.params,
    times=times,
    doses=doses
)

cmax = max(result.concentrations)
tmax = times[result.concentrations.index(cmax)]
print(f"Cmax: {cmax:.2f}")
print(f"Tmax: {tmax} h")
print("✓ Simulation successful")
```

---

## See Also

- [Monolix Import](monolix.md) - Import Monolix projects
- [CDISC Import](cdisc.md) - Import CDISC data
- [Simulation](../models/index.md) - Model simulation

