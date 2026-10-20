# Monolix Model Import

Import Monolix project files (.mlxtran) into NeoPKPD Python.

---

## Overview

```python
from neopkpd.import_ import import_monolix

model = import_monolix("project.mlxtran")
print(f"Model: {model.model_kind}")
print(f"Parameters: {model.params}")
```

---

## Quick Start

### Basic Import

```python
from neopkpd.import_ import import_monolix

# Import Monolix project
model = import_monolix("project.mlxtran")

# Access model information
print(f"Model type: {model.model_kind}")
print(f"Source: {model.source_format}")

# Access parameters
for name, value in model.params.items():
    print(f"  {name} = {value}")
```

### With Dose Events

```python
# Specify doses if needed
doses = [
    {"time": 0.0, "amount": 100.0, "route": "oral"}
]

model = import_monolix("project.mlxtran", doses=doses)
```

---

## Supported Model Types

### One-Compartment Models

| Monolix Library Model | NeoPKPD Model |
|----------------------|----------------|
| `pk_bolus1cpt_Vk_PLASMA` | `OneCompIVBolus` |
| `pk_bolus1cpt_VCl_PLASMA` | `OneCompIVBolus` |
| `pk_oral1cpt_kaVk_PLASMA` | `OneCompOralFirstOrder` |
| `pk_oral1cpt_kaVCl_PLASMA` | `OneCompOralFirstOrder` |
| `pk_infusion1cpt_VCl_PLASMA` | `OneCompIVInfusion` |

### Two-Compartment Models

| Monolix Library Model | NeoPKPD Model |
|----------------------|----------------|
| `pk_bolus2cpt_V1k12k21k_PLASMA` | `TwoCompIVBolus` |
| `pk_bolus2cpt_V1ClQ2V2_PLASMA` | `TwoCompIVBolus` |
| `pk_oral2cpt_kaV1k12k21k_PLASMA` | `TwoCompOral` |
| `pk_oral2cpt_kaV1ClQ2V2_PLASMA` | `TwoCompOral` |

### Three-Compartment Models

| Monolix Library Model | NeoPKPD Model |
|----------------------|----------------|
| `pk_bolus3cpt_V1ClQ2V2Q3V3_PLASMA` | `ThreeCompIVBolus` |

### Special Models

| Monolix Library Model | NeoPKPD Model |
|----------------------|----------------|
| `pk_bolus1cpt_VVmKm_PLASMA` | `MichaelisMentenElimination` |

---

## ImportedModel Class

```python
@dataclass
class ImportedModel:
    """Result of importing a Monolix model."""

    source_format: str          # "monolix"
    source_file: str            # Path to mlxtran file
    model_kind: str             # NeoPKPD model name

    # Fixed effects
    params: dict[str, float]    # Parameter values
    theta_init: list[float]     # Initial values
    theta_names: list[str]      # Parameter names

    # Random effects
    omega_init: list[list[float]]  # Omega matrix
    omega_names: list[str]         # IIV parameter names

    # Residual error
    sigma_type: str             # Error model type
    sigma_init: float           # Error parameter

    # Quality
    warnings: list[str]         # Import warnings

    # Metadata
    metadata: dict[str, Any]    # Additional info
```

---

## Parameter Extraction

### Population Parameters

```python
model = import_monolix("project.mlxtran")

# Access population parameters
print("Population Parameters:")
for name, value in model.params.items():
    print(f"  {name} = {value}")

# Example output:
# Ka = 1.5
# V = 50.0
# CL = 5.0
```

### Random Effects

```python
import math

model = import_monolix("project.mlxtran")

print("Random Effects:")
for i, name in enumerate(model.omega_names):
    omega_sq = model.omega_init[i][i]
    # For lognormal IIV, convert to CV%
    cv = (math.exp(omega_sq) - 1) ** 0.5 * 100
    print(f"  {name}: ω² = {omega_sq:.4f} (CV ≈ {cv:.1f}%)")
```

### Variability Transformations

| Monolix Variability | NeoPKPD Transform | Formula |
|---------------------|-------------------|---------|
| `lognormal` | `exponential` | θᵢ = θ_pop · e^ηᵢ |
| `normal` | `additive` | θᵢ = θ_pop + ηᵢ |
| `logitnormal` | `logit` | Logit transform |
| `none` | No IIV | Fixed to pop value |

---

## Error Model Import

### Supported Error Models

```python
model = import_monolix("project.mlxtran")

print(f"Error type: {model.sigma_type}")
print(f"Sigma value: {model.sigma_init}")

# Supported types:
# - "proportional": Y = F * (1 + b*ε)
# - "additive": Y = F + a*ε
# - "combined": Y = F * (1 + b*ε₁) + a*ε₂
```

### Error Model Mapping

| Monolix Error | NeoPKPD Type |
|---------------|---------------|
| `proportional` | `proportional` |
| `constant` | `additive` |
| `combined1` | `combined` |
| `combined2` | `combined` |

---

## Covariate Models

### Extraction

```python
model = import_monolix("project.mlxtran")

# Covariate effects in metadata
effects = model.metadata.get("covariate_effects", [])

for effect in effects:
    print(f"{effect['covariate']} on {effect['parameter']}:")
    print(f"  Type: {effect['type']}")
    print(f"  Coefficient: {effect['coefficient']}")
    print(f"  Reference: {effect['reference_value']}")
```

### Supported Covariate Types

| Type | Description | Formula |
|------|-------------|---------|
| `power` | Power model | θ · (COV/REF)^β |
| `linear` | Linear model | θ · (1 + β·(COV-REF)) |
| `exponential` | Exponential | θ · exp(β·(COV-REF)) |
| `categorical` | Category effect | θ · exp(β·I) |

---

## Using Imported Models

### Simulation

```python
from neopkpd.import_ import import_monolix
from neopkpd import simulate

# Import model
model = import_monolix("project.mlxtran")

# Simulate
times = list(range(0, 49))
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
from neopkpd import simulate_population

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

## Unsupported Features

### Models Not Supported

| Model Type | Description |
|------------|-------------|
| PD models | Pharmacodynamic |
| Turnover models | Indirect response |
| Transit compartment | Absorption chain |
| Mixture models | Subpopulations |
| Markov models | State transitions |
| Time-to-event | Survival |
| Count data | Poisson/NegBin |
| Categorical | Ordered response |

### Features with Warnings

| Feature | Handling |
|---------|----------|
| Lag time (Tlag) | Ignored, assumes Tlag=0 |
| Bioavailability (F) | Ignored, assumes F=1 |
| Complex covariates | May be simplified |

---

## Validation

### Check Import Quality

```python
model = import_monolix("project.mlxtran")

# Check warnings
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
expected = ["Ka", "CL", "V"]  # For 1-comp oral
for param in expected:
    if param in model.params:
        print(f"✓ {param} = {model.params[param]}")
    else:
        print(f"❌ Missing: {param}")
```

---

## Complete Example

```python
from neopkpd.import_ import import_monolix
from neopkpd import simulate
import math

# Import Monolix project
print("Importing Monolix project...")
model = import_monolix("project.mlxtran")

# Display results
print("=" * 50)
print("Monolix Import Results")
print("=" * 50)

print(f"\nSource: {model.source_file}")
print(f"Model type: {model.model_kind}")

print("\n--- Population Parameters ---")
for name, value in model.params.items():
    print(f"  {name} = {value}")

print("\n--- Inter-Individual Variability ---")
if model.omega_names:
    for i, name in enumerate(model.omega_names):
        omega_sq = model.omega_init[i][i]
        cv = (math.exp(omega_sq) - 1) ** 0.5 * 100
        print(f"  {name}: ω² = {omega_sq:.4f} (CV ≈ {cv:.1f}%)")
else:
    print("  None")

print(f"\n--- Residual Error ---")
print(f"  Type: {model.sigma_type}")
print(f"  Value: {model.sigma_init}")

# Covariate effects
effects = model.metadata.get("covariate_effects", [])
if effects:
    print("\n--- Covariate Effects ---")
    for eff in effects:
        print(f"  {eff['covariate']} on {eff['parameter']}: {eff['type']}")

print("\n--- Warnings ---")
if model.warnings:
    for w in model.warnings:
        print(f"  ⚠️ {w}")
else:
    print("  None")

# Validation simulation
print("\n--- Validation Simulation ---")
times = list(range(0, 73, 1))
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

## Model Library Reference

### Standard Naming Convention

```
pk_{route}{n}cpt_{parameters}_{observation}
```

- **Route**: `bolus`, `oral`, `infusion`
- **n**: Number of compartments (1, 2, 3)
- **Parameters**: `V`, `Cl`, `k`, `ka`, `Q`, `Vm`, `Km`
- **Observation**: `PLASMA`, `EFFECT`

### Examples

```
pk_oral1cpt_kaVCl_PLASMA     → 1-comp oral, ka/V/Cl parameterization
pk_bolus2cpt_V1ClQ2V2_PLASMA → 2-comp IV, CL/V1/Q/V2 parameterization
pk_bolus1cpt_VVmKm_PLASMA    → Michaelis-Menten elimination
```

---

## Auto-Detection

### Using import_model()

```python
from neopkpd.import_ import import_model

# Auto-detect format from extension
model = import_model("run001.ctl")      # → NONMEM
model = import_model("project.mlxtran") # → Monolix

# Or specify format explicitly
model = import_model("model.txt", format="nonmem")
model = import_model("project.xml", format="monolix")
```

---

## See Also

- [NONMEM Import](nonmem.md) - Import NONMEM models
- [CDISC Import](cdisc.md) - Import CDISC data
- [Simulation](../models/index.md) - Model simulation

