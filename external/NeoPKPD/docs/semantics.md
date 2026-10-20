# Semantics

NeoPKPD uses **semantic versioning** for numerical behavior to ensure reproducibility across versions.

## Version Constants

Three independent version numbers track different aspects of system behavior:

```julia
NEOPKPD_VERSION = "0.1.0"           # Software version
EVENT_SEMANTICS_VERSION = "1.0.0"    # Dose event handling
SOLVER_SEMANTICS_VERSION = "1.0.0"   # ODE solver behavior
ARTIFACT_SCHEMA_VERSION = "1.0.0"    # JSON artifact format
```

---

## Event Semantics (v1.0.0)

Event semantics define how dose events are normalized and applied.

### Dose Normalization Rules

Given doses and simulation window `[t0, t1]`:

1. **Doses at t0**: Added to initial state `u0`
2. **Doses in (t0, t1]**: Applied via ODE callbacks
3. **Doses outside [t0, t1]**: Ignored
4. **Duplicate times**: Summed together

```julia
# Example: doses at [0, 0, 12] with t0=0, t1=24
# - Doses at t=0 (100 + 50 = 150) → added to u0
# - Dose at t=12 → callback
doses = [
    DoseEvent(0.0, 100.0),
    DoseEvent(0.0, 50.0),
    DoseEvent(12.0, 75.0)
]
```

### Dose Target

Each model defines which compartment receives doses:

| Model | Dose Target |
|-------|-------------|
| `OneCompIVBolus` | Central (index 1) |
| `OneCompOralFirstOrder` | Gut (index 1) |

### Callback Timing

Doses in `(t0, t1]` are applied as **preset time callbacks**:

- Exact time precision (no interpolation)
- State discontinuity at dose time
- Solver restarts after callback

---

## Solver Semantics (v1.0.0)

Solver semantics define ODE solving behavior and defaults.

### Supported Algorithms

| Algorithm | Type | Use Case |
|-----------|------|----------|
| `Tsit5` | Explicit Runge-Kutta (5th order) | Non-stiff problems (default) |
| `Rosenbrock23` | Rosenbrock method | Stiff problems |

### Tolerance Interpretation

- `reltol`: Relative error tolerance per step
- `abstol`: Absolute error tolerance per step
- Error estimate: `err ≤ abstol + reltol * |u|`

### Golden Standard Settings

For reproducibility, use:

```julia
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)
```

### Output Grid Alignment

Results are aligned exactly to `saveat` times:

- No interpolation artifacts
- Deterministic output regardless of internal step size
- Callback discontinuities preserved

### State Continuity

For segmented simulations (IOV, time-varying covariates):

- End state of segment N → initial state of segment N+1
- No discontinuity except at dose times

---

## Artifact Schema (v1.0.0)

Artifact schema defines the JSON structure for serialized simulations.

### Common Fields

Every artifact contains:

```json
{
  "artifact_schema_version": "1.0.0",
  "semantics_fingerprint": {
    "artifact_schema_version": "1.0.0",
    "event_semantics_version": "1.0.0",
    "solver_semantics_version": "1.0.0"
  }
}
```

### Artifact Types

| Type | `artifact_type` Value |
|------|----------------------|
| Single execution | `"single"` or absent |
| Population | `"population"` |
| Sensitivity (single) | `"sensitivity_single"` |
| Sensitivity (population) | `"sensitivity_population"` |

### Numeric Precision

- All floating-point values stored with full precision
- No rounding or truncation
- JSON uses IEEE 754 double representation

---

## Semantics Fingerprint

Every artifact includes a fingerprint for version tracking:

```julia
semantics_fingerprint() = Dict(
    "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
    "event_semantics_version" => EVENT_SEMANTICS_VERSION,
    "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
)
```

### Fingerprint Verification

When replaying artifacts:

1. Read stored fingerprint
2. Compare to current versions
3. Warn if mismatch (results may differ)

---

## Version Bump Policy

### When to Bump Event Semantics

- Change to dose normalization rules
- Change to dose application order
- Change to callback timing
- Change to dose target assignment

### When to Bump Solver Semantics

- Change to default tolerances
- Change to algorithm selection logic
- Change to output grid alignment
- Change to segmentation behavior

### When to Bump Artifact Schema

- Change to JSON field names
- Change to JSON structure
- Addition of required fields
- Removal of fields

### Process for Bumping

1. **Identify Change**: Document what behavior changes
2. **Bump Version**: Increment appropriate version number
3. **Update Golden**: Re-generate all affected golden artifacts
4. **Update Tests**: Ensure tests pass with new behavior
5. **Document**: Add entry to changelog

---

## Compatibility

### Forward Compatibility

Newer NeoPKPD versions can replay older artifacts:

- Missing fields use defaults
- Schema upgrades are automatic

### Backward Compatibility

Older NeoPKPD versions may not replay newer artifacts:

- Unknown fields are ignored
- Missing required fields cause errors

### Cross-Platform Compatibility

Same artifact should produce identical results on:

- Different operating systems
- Different Julia versions (≥1.9)
- Different hardware architectures

---

## IIV/IOV Semantics

### Random Number Generation

- Uses `Random.Xoshiro` PRNG
- Seed fully determines sequence
- Separate seeds for IIV and IOV

### Log-Normal Transform

```
theta_i = theta_pop * exp(eta)
eta ~ Normal(0, omega^2)
```

### Occasion Semantics (v1)

Occasion boundaries derived from dose times:

1. `t0` always starts occasion 1
2. Each unique dose time > t0 starts new occasion
3. Times sorted and deduplicated
4. Only times in [t0, t1] considered

---

## Covariate Semantics

### Application Order

Covariate effects applied in declaration order:

```julia
effects = [
    CovariateEffect(PowerCovariate(), :CL, :WT, 0.75, 70.0),
    CovariateEffect(LinearCovariate(), :CL, :AGE, -0.01, 40.0)
]
# Applied: CL → (WT effect) → (AGE effect)
```

### Time-Varying Interpolation

**Step**: Value from rightmost knot ≤ t
```
series.times = [0, 24, 48]
series.values = [1, 2, 3]
value_at(15) = 1  # Uses knot at 0
value_at(24) = 2  # Uses knot at 24
```

**Linear**: Interpolate between knots
```
value_at(12) = 1.5  # Interpolate 0→24
```

---

## Sensitivity Semantics

### Perturbation Application

Perturbations applied in declaration order:

```julia
plan.perturbations = [
    Perturbation(RelativePerturbation(), :CL, 0.1),
    Perturbation(AbsolutePerturbation(), :V, 5.0)
]
# Applied: CL → CL * 1.1, then V → V + 5
```

### Metric Computation

All metrics computed element-wise on output series:

```
max_abs_delta = max(|base[i] - pert[i]|)
max_rel_delta = max(|base[i] - pert[i]| / |base[i]|)
l2_norm_delta = sqrt(sum((base[i] - pert[i])^2))
```

---

## Examples

### Verifying Semantics Match

```julia
using NeoPKPD

# Check current versions
println("Event semantics: ", EVENT_SEMANTICS_VERSION)
println("Solver semantics: ", SOLVER_SEMANTICS_VERSION)
println("Artifact schema: ", ARTIFACT_SCHEMA_VERSION)

# Check artifact fingerprint
artifact = read_execution_json("my_artifact.json")
fingerprint = artifact["semantics_fingerprint"]

if fingerprint["event_semantics_version"] != EVENT_SEMANTICS_VERSION
    @warn "Event semantics version mismatch"
end
```

### Ensuring Reproducibility

```julia
# Use explicit, high-precision settings
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Use deterministic seeds
iiv = IIVSpec(LogNormalIIV(), omegas, UInt64(12345), n)

# Save artifact for future validation
write_execution_json("reproducible.json"; ...)
```
