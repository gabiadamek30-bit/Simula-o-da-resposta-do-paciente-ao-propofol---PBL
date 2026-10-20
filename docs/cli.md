# CLI Reference

NeoPKPD provides a comprehensive command-line interface for PK/PD simulations, estimation, VPC, NCA, trial simulation, and model import.

## Installation

The CLI is located in `packages/cli/bin/neopkpd` and requires Julia with the NeoPKPD and NeoPKPDCLI packages.

```bash
# Install dependencies
julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'
julia --project=packages/cli -e 'using Pkg; Pkg.instantiate()'

# Make executable
chmod +x packages/cli/bin/neopkpd

# Run CLI
./packages/cli/bin/neopkpd <command> [options]
```

---

## Commands

### version

Display version information for all NeoPKPD components.

```bash
./packages/cli/bin/neopkpd version
```

**Output**:

```
NeoPKPD 0.1.0
Event semantics: 1.0.0
Solver semantics: 1.0.0
Artifact schema: 1.0.0
```

---

### simulate

Run a PK/PD simulation from a JSON specification.

```bash
./packages/cli/bin/neopkpd simulate --spec <path> [--out <output_path>]
```

**Options**:

| Option | Required | Description |
|--------|----------|-------------|
| `--spec` | Yes | Path to simulation specification JSON |
| `--out` | No | Output path for result artifact |

**Example**:

```bash
./packages/cli/bin/neopkpd simulate --spec simulation.json --out result.json
```

---

### population

Run a population simulation with IIV/IOV.

```bash
./packages/cli/bin/neopkpd population --spec <path> [--out <output_path>]
```

---

### estimate

Run NLME parameter estimation (FOCE-I, SAEM, or Laplacian).

```bash
./packages/cli/bin/neopkpd estimate --spec <path> [--out <output_path>]
```

**Specification includes**: observed data, model, estimation method, initial values, bounds.

See [Parameter Estimation](julia/estimation/index.md) for specification format.

---

### nca

Run non-compartmental analysis.

```bash
./packages/cli/bin/neopkpd nca --spec <path> [--out <output_path>]
```

See [NCA](julia/nca/index.md) for specification format.

---

### vpc

Compute Visual Predictive Check.

```bash
./packages/cli/bin/neopkpd vpc --spec <path> [--out <output_path>]
```

See [VPC](julia/vpc/index.md) for specification format.

---

### trial

Run clinical trial simulation.

```bash
./packages/cli/bin/neopkpd trial --spec <path> [--out <output_path>]
```

Supports parallel, crossover, dose-escalation, and bioequivalence designs.

See [Trial Simulation](julia/trial/index.md) for specification format.

---

### import

Import models from NONMEM or Monolix.

```bash
./packages/cli/bin/neopkpd import --input <path> --format <format> [--out <output_path>]
```

**Options**:

| Option | Required | Description |
|--------|----------|-------------|
| `--input` | Yes | Path to model file (.ctl or .mlxtran) |
| `--format` | Yes | Format: `nonmem` or `monolix` |
| `--out` | No | Output path for converted model |

**Examples**:

```bash
# Import NONMEM control file
./packages/cli/bin/neopkpd import --input run001.ctl --format nonmem --out model.json

# Import Monolix project
./packages/cli/bin/neopkpd import --input project.mlxtran --format monolix --out model.json
```

See [Model Import](julia/import/index.md) for details.

---

### sensitivity

Run parameter sensitivity analysis.

```bash
./packages/cli/bin/neopkpd sensitivity --spec <path> [--out <output_path>]
```

---

### metrics

Compute PK/PD metrics from simulation results.

```bash
./packages/cli/bin/neopkpd metrics --artifact <path> --metrics <list>
```

**Examples**:

```bash
./packages/cli/bin/neopkpd metrics --artifact result.json --metrics cmax,tmax,auc
```

---

### replay

Replay an execution artifact to verify reproducibility.

```bash
./packages/cli/bin/neopkpd replay --artifact <path> [--out <output_path>]
```

**Options**:

| Option | Required | Description |
|--------|----------|-------------|
| `--artifact` | Yes | Path to artifact JSON file |
| `--out` | No | Path to write replayed artifact |

**Supported Artifact Types**:

- Single execution (`artifact_type: "single"` or missing)
- Population execution (`artifact_type: "population"`)
- Single sensitivity (`artifact_type: "sensitivity_single"`)
- Population sensitivity (`artifact_type: "sensitivity_population"`)
- Estimation results (`artifact_type: "estimation"`)

**Examples**:

```bash
# Replay and verify
./packages/cli/bin/neopkpd replay --artifact validation/golden/pk_iv_bolus.json

# Replay and save output
./packages/cli/bin/neopkpd replay --artifact my_simulation.json --out replayed.json

# Replay population artifact
./packages/cli/bin/neopkpd replay --artifact validation/golden/population_iv_bolus.json
```

---

### validate-golden

Run the full golden artifact validation suite.

```bash
./packages/cli/bin/neopkpd validate-golden
```

**What It Does**:

1. Finds all golden artifacts in `validation/golden/`
2. Replays each artifact
3. Compares replayed results to stored results
4. Reports any discrepancies

**Output**:

```
Validating golden artifacts...
  pk_iv_bolus.json: PASS
  pk_oral_first_order.json: PASS
  population_iv_bolus.json: PASS
  pkpd_direct_emax.json: PASS
  sensitivity_single.json: PASS
  ...
All golden artifacts validated successfully.
```

**Exit Codes**:

| Code | Meaning |
|------|---------|
| 0 | All validations passed |
| 1 | One or more validations failed |

---

### help

Display help for any command.

```bash
./packages/cli/bin/neopkpd help
./packages/cli/bin/neopkpd help simulate
./packages/cli/bin/neopkpd help estimate
```

---

## Artifact Format

All NeoPKPD artifacts are JSON files with a consistent structure.

### Common Fields

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

### Single Execution Artifact

```json
{
  "artifact_schema_version": "1.0.0",
  "execution_mode": "pk",
  "model_spec": {
    "kind": "OneCompIVBolus",
    "name": "example",
    "params": {"CL": 5.0, "V": 50.0},
    "doses": [{"time": 0.0, "amount": 100.0}]
  },
  "grid": {
    "t0": 0.0,
    "t1": 24.0,
    "saveat": [0.0, 1.0, 2.0, ...]
  },
  "solver": {
    "alg": "Tsit5",
    "reltol": 1e-10,
    "abstol": 1e-12,
    "maxiters": 10000000
  },
  "result": {
    "t": [0.0, 1.0, 2.0, ...],
    "states": {"A_central": [...]},
    "observations": {"conc": [...]},
    "metadata": {...}
  }
}
```

### Population Artifact

```json
{
  "artifact_type": "population",
  "population_spec": {
    "base_model_spec": {...},
    "iiv": {
      "kind": "LogNormalIIV",
      "omegas": {"CL": 0.3, "V": 0.2},
      "seed": 12345,
      "n": 100
    },
    "iov": null,
    "covariate_model": null,
    "covariates": []
  },
  "result": {
    "individuals": [...],
    "params": [...],
    "summaries": {
      "conc": {
        "observation": "conc",
        "probs": [0.05, 0.95],
        "mean": [...],
        "median": [...],
        "quantiles": {...}
      }
    },
    "metadata": {...}
  }
}
```

---

## Integration Examples

### CI/CD Pipeline

```yaml
# .github/workflows/validate.yml
name: Golden Validation
on: [push, pull_request]
jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: julia-actions/setup-julia@v1
      - run: julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'
      - run: ./bin/neopkpd validate-golden
```

### Batch Replay

```bash
#!/bin/bash
# replay_all.sh - Replay all artifacts in a directory

for artifact in artifacts/*.json; do
    echo "Replaying: $artifact"
    ./bin/neopkpd replay --artifact "$artifact" --out "replayed/$(basename $artifact)"
done
```

### Compare Artifacts

```bash
# Generate new artifact
julia --project=packages/core my_simulation.jl

# Compare with golden
diff <(jq -S . validation/golden/pk_iv_bolus.json) <(jq -S . my_artifact.json)
```

---

## Troubleshooting

### Julia Not Found

```
Error: julia: command not found
```

**Solution**: Ensure Julia is installed and in your PATH.

### Package Not Installed

```
Error: ArgumentError: Package NeoPKPD not found
```

**Solution**: Install dependencies:

```bash
julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'
julia --project=packages/cli -e 'using Pkg; Pkg.instantiate()'
```

### Artifact Schema Mismatch

```
Warning: Artifact schema version mismatch
```

**Cause**: Artifact was created with a different version of NeoPKPD.

**Solution**: Re-generate the artifact or update your NeoPKPD installation.

---

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `JULIA_PROJECT` | Julia project path | Auto-detected |
| `NEOPKPD_ROOT` | Repository root | Auto-detected |

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error |
| 2 | Invalid arguments |
| 3 | File not found |
| 4 | Validation failure |
