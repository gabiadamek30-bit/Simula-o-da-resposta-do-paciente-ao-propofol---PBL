# Reproducibility

NeoPKPD is designed for **complete reproducibility** of simulation results across time, platforms, and versions.

## Core Principles

1. **Deterministic Execution**: Same inputs always produce identical outputs
2. **Complete Artifacts**: All information needed to reproduce results is stored
3. **Semantic Versioning**: Changes to numerical behavior are tracked and versioned
4. **Replay Validation**: Artifacts can be replayed to verify reproducibility

---

## Execution Artifacts

Every simulation can be serialized to a JSON artifact containing everything needed for reproduction.

### What Artifacts Contain

| Component | Description |
|-----------|-------------|
| Model Specification | Kind, parameters, doses |
| Simulation Grid | t0, t1, saveat times |
| Solver Settings | Algorithm, tolerances, max iterations |
| Results | Time series, states, observations |
| Semantics Fingerprint | Version information |
| Metadata | Engine version, settings snapshot |

### Creating Artifacts

**Julia**:

```julia
using NeoPKPD

spec = ModelSpec(...)
grid = SimGrid(...)
solver = SolverSpec(...)
result = simulate(spec, grid, solver)

write_execution_json(
    "my_simulation.json";
    model_spec=spec,
    grid=grid,
    solver=solver,
    result=result
)
```

**Python**:

```python
import neopkpd

neopkpd.write_single_artifact(
    "my_simulation.json",
    model={...},
    grid={...},
    solver={...}
)
```

### Artifact Types

| Type | Function | Contents |
|------|----------|----------|
| Single | `write_execution_json` | Single simulation |
| Population | `write_population_json` | Population simulation |
| Sensitivity (Single) | `write_sensitivity_json` | Single sensitivity analysis |
| Sensitivity (Population) | `write_population_sensitivity_json` | Population sensitivity |

---

## Replay System

Artifacts can be replayed to reproduce simulations.

### Replay Functions

```julia
# Read artifact
artifact = read_execution_json("simulation.json")

# Replay single execution
result = replay_execution(artifact)

# Replay population
result = replay_population_execution(artifact)

# Replay sensitivity
result = replay_sensitivity_execution(artifact)

# Replay population sensitivity
result = replay_population_sensitivity_execution(artifact)
```

### Python Replay

```python
import neopkpd

result = neopkpd.replay_artifact("simulation.json")
```

### CLI Replay

```bash
./bin/neopkpd replay --artifact simulation.json
./bin/neopkpd replay --artifact simulation.json --out replayed.json
```

---

## Golden Artifacts

Golden artifacts are the **regression contract** for NeoPKPD.

### Location

```
validation/golden/
├── pk_iv_bolus.json
├── pk_oral_first_order.json
├── pkpd_direct_emax.json
├── pkpd_indirect_coupled.json
├── population_iv_bolus.json
├── population_iov.json
├── population_covariates.json
├── population_time_varying.json
├── sensitivity_single.json
└── sensitivity_population.json
```

### Purpose

- **Regression Testing**: Verify numerical behavior hasn't changed
- **Contract Documentation**: Define expected outputs for given inputs
- **Cross-Version Validation**: Ensure compatibility across releases

### Generating Golden Artifacts

```bash
# Run generation script
julia validation/scripts/generate_golden_artifacts.jl
```

### Validating Golden Artifacts

```bash
# CLI validation
./bin/neopkpd validate-golden

# Script validation
julia validation/scripts/run_golden_validation.jl
```

---

## Validation Process

### What Validation Checks

1. **Artifact Readability**: JSON parses correctly
2. **Schema Compliance**: Required fields present
3. **Replay Success**: Simulation runs without error
4. **Result Match**: Replayed results match stored results

### Comparison Criteria

Results are compared with exact floating-point equality:

```julia
@test replayed.t == stored.t
@test replayed.observations[:conc] == stored.observations[:conc]
```

### Handling Failures

When validation fails:

1. **Bug**: Fix the code, validation should pass
2. **Intentional Change**:
   - Bump semantic version
   - Re-generate golden artifacts
   - Document change

---

## Determinism Guarantees

### Random Number Generation

IIV and IOV use deterministic seeding:

```julia
iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.3),
    UInt64(12345),  # Deterministic seed
    100
)
```

**Guarantees**:

- Same seed → same random sequence
- Independent of simulation order
- Reproducible across platforms

### Output Grid

Results are aligned exactly to `saveat` times:

```julia
grid = SimGrid(0.0, 24.0, [0.0, 1.0, 2.0])
result = simulate(spec, grid, solver)

# result.t is exactly [0.0, 1.0, 2.0]
```

### Floating-Point Precision

- IEEE 754 double precision throughout
- No intermediate rounding
- Full precision in artifacts

---

## CI/CD Integration

### GitHub Actions Example

```yaml
name: Golden Validation
on: [push, pull_request]

jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - uses: julia-actions/setup-julia@v1
        with:
          version: '1.9'

      - name: Install dependencies
        run: |
          julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'

      - name: Run tests
        run: |
          julia --project=packages/core -e 'using Pkg; Pkg.test()'

      - name: Validate golden artifacts
        run: |
          ./bin/neopkpd validate-golden
```

### Pre-Commit Hook

```bash
#!/bin/bash
# .git/hooks/pre-commit

# Validate golden artifacts before commit
./bin/neopkpd validate-golden
if [ $? -ne 0 ]; then
    echo "Golden validation failed. Commit aborted."
    exit 1
fi
```

---

## Troubleshooting

### "Results don't match stored values"

**Causes**:

1. Code change affected numerical output
2. Tolerance settings changed
3. Random seed handling changed

**Solutions**:

- Review recent changes
- Check semantic versions
- Re-generate golden artifacts if intentional

### "Artifact can't be read"

**Causes**:

1. Corrupted JSON
2. Schema version mismatch
3. Missing required fields

**Solutions**:

- Validate JSON syntax
- Check artifact schema version
- Regenerate artifact

### "Replay produces different results"

**Causes**:

1. Semantic version mismatch
2. Solver configuration differences
3. Platform-specific floating-point behavior

**Solutions**:

- Use same NeoPKPD version
- Match solver settings exactly
- Report platform-specific issues

---

## Best Practices

### 1. Always Save Artifacts

```julia
# Save every important simulation
write_execution_json("research/exp_001.json"; ...)
```

### 2. Use Explicit Seeds

```julia
# Never rely on implicit randomness
iiv = IIVSpec(..., UInt64(42), n)  # Explicit seed
```

### 3. Use High Precision

```julia
# Golden-standard tolerances
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)
```

### 4. Version Control Artifacts

```bash
# Track golden artifacts in git
git add validation/golden/*.json
git commit -m "Update golden artifacts for v1.1.0"
```

### 5. Document Changes

When numerical behavior changes:

1. Update `CHANGELOG.md`
2. Bump appropriate semantic version
3. Regenerate affected artifacts
4. Update documentation

---

## Artifact Lifecycle

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│   Create    │────▶│    Store    │────▶│   Replay    │
│ Simulation  │     │  Artifact   │     │   Verify    │
└─────────────┘     └─────────────┘     └─────────────┘
       │                   │                   │
       ▼                   ▼                   ▼
  SimResult           JSON File          SimResult
  (in memory)       (persistent)        (reproduced)
                          │
                          ▼
                    ┌─────────────┐
                    │   Version   │
                    │   Control   │
                    └─────────────┘
```

---

## Example: Complete Reproducibility Workflow

```julia
using NeoPKPD

# 1. Define simulation (all parameters explicit)
spec = ModelSpec(
    OneCompIVBolus(),
    "reproducibility_example",
    OneCompIVBolusParams(5.0, 50.0),
    [DoseEvent(0.0, 100.0)]
)

grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# 2. Run simulation
result = simulate(spec, grid, solver)

# 3. Save artifact
write_execution_json(
    "reproducible_sim.json";
    model_spec=spec,
    grid=grid,
    solver=solver,
    result=result
)

# 4. Later: replay and verify
artifact = read_execution_json("reproducible_sim.json")
replayed = replay_execution(artifact)

# 5. Confirm exact match
@assert result.t == replayed.t
@assert result.observations[:conc] == replayed.observations[:conc]
println("Reproducibility verified!")
```
