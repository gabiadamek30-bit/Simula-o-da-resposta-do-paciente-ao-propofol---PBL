# NeoPKPD CLI

**Command-line interface for NeoPKPD pharmacometric simulations**

The NeoPKPD CLI provides command-line access to all core functionality including simulation, estimation, VPC, NCA, trial simulation, and model import.

## Installation

```bash
# Ensure Julia is installed and in PATH
julia --version

# Install dependencies
julia --project=packages/cli -e 'using Pkg; Pkg.instantiate()'

# Make CLI executable
chmod +x packages/cli/bin/neopkpd

# Verify installation
./packages/cli/bin/neopkpd version
```

## Commands

### version

Display version information for all NeoPKPD components.

```bash
./packages/cli/bin/neopkpd version
```

Output:
```
NeoPKPD 0.1.0
Event semantics: 1.0.0
Solver semantics: 1.0.0
Artifact schema: 1.0.0
```

### simulate

Run a PK/PD simulation from a JSON specification.

```bash
./packages/cli/bin/neopkpd simulate --spec simulation.json --out result.json
```

**Options:**
| Option | Required | Description |
|--------|----------|-------------|
| `--spec` | Yes | Path to simulation specification JSON |
| `--out` | No | Output path for result artifact |

**Specification Format:**
```json
{
  "model": {
    "kind": "OneCompIVBolus",
    "params": {"CL": 5.0, "V": 50.0},
    "doses": [{"time": 0.0, "amount": 100.0, "duration": 0.0}]
  },
  "grid": {
    "t0": 0.0,
    "t1": 24.0,
    "saveat": [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
  },
  "solver": {
    "alg": "Tsit5",
    "reltol": 1e-10,
    "abstol": 1e-12,
    "maxiters": 10000000
  }
}
```

**Supported Model Kinds:**
- `OneCompIVBolus`, `OneCompOralFirstOrder`
- `TwoCompIVBolus`, `TwoCompOral`
- `ThreeCompIVBolus`
- `TransitAbsorption`, `MichaelisMentenElimination`
- `DirectEmax`, `SigmoidEmax`
- `BiophaseEquilibration`, `IndirectResponseTurnover`

### population

Run a population simulation with IIV/IOV.

```bash
./packages/cli/bin/neopkpd population --spec pop_spec.json --out pop_result.json
```

**Specification Format:**
```json
{
  "base_model": {
    "kind": "OneCompIVBolus",
    "params": {"CL": 5.0, "V": 50.0},
    "doses": [{"time": 0.0, "amount": 100.0}]
  },
  "iiv": {
    "kind": "LogNormalIIV",
    "omegas": {"CL": 0.3, "V": 0.2},
    "seed": 12345,
    "n": 100
  },
  "grid": {"t0": 0.0, "t1": 24.0, "saveat": [0, 1, 2, 4, 8, 12, 24]},
  "solver": {"alg": "Tsit5", "reltol": 1e-10, "abstol": 1e-12}
}
```

### estimate

Run NLME parameter estimation.

```bash
./packages/cli/bin/neopkpd estimate --spec estimate_spec.json --out estimate_result.json
```

**Specification Format:**
```json
{
  "data": "observed_data.csv",
  "model": {
    "kind": "OneCompOralFirstOrder",
    "params": {"Ka": 1.0, "CL": 5.0, "V": 50.0}
  },
  "estimation": {
    "method": "FOCE-I",
    "theta_init": [1.0, 5.0, 50.0],
    "theta_lower": [0.1, 0.1, 1.0],
    "theta_upper": [10.0, 100.0, 500.0],
    "omega_init": [[0.09, 0, 0], [0, 0.09, 0], [0, 0, 0.04]],
    "omega_structure": "diagonal",
    "sigma": {"kind": "proportional", "value": 0.1},
    "max_iter": 500,
    "compute_se": true
  },
  "grid": {"t0": 0.0, "t1": 24.0},
  "solver": {"alg": "Tsit5", "reltol": 1e-8, "abstol": 1e-10}
}
```

**Estimation Methods:**
- `FOCE-I` - First-Order Conditional Estimation with Interaction
- `SAEM` - Stochastic Approximation EM
- `Laplacian` - Laplacian approximation

### nca

Run non-compartmental analysis.

```bash
./packages/cli/bin/neopkpd nca --spec nca_spec.json --out nca_result.json
```

**Specification Format:**
```json
{
  "data": {
    "times": [0, 0.5, 1, 2, 4, 8, 12, 24],
    "conc": [0, 1.8, 2.0, 1.5, 1.0, 0.5, 0.25, 0.06],
    "dose": 100.0
  },
  "config": {
    "method": "log_linear",
    "lambda_z_min_points": 3,
    "extrapolation_max_pct": 20.0
  }
}
```

### vpc

Compute Visual Predictive Check.

```bash
./packages/cli/bin/neopkpd vpc --spec vpc_spec.json --out vpc_result.json
```

**Specification Format:**
```json
{
  "observed_data": "observed.csv",
  "population_spec": {
    "base_model": {"kind": "OneCompIVBolus", "params": {"CL": 5.0, "V": 50.0}},
    "iiv": {"kind": "LogNormalIIV", "omegas": {"CL": 0.3, "V": 0.2}, "n": 1000}
  },
  "config": {
    "pi_levels": [0.05, 0.50, 0.95],
    "ci_level": 0.95,
    "binning": {"strategy": "quantile", "n_bins": 10},
    "prediction_corrected": false,
    "n_bootstrap": 500,
    "seed": 12345
  },
  "grid": {"t0": 0.0, "t1": 24.0, "saveat": [0, 1, 2, 4, 8, 12, 24]},
  "solver": {"alg": "Tsit5", "reltol": 1e-8, "abstol": 1e-10}
}
```

### trial

Run clinical trial simulation.

```bash
./packages/cli/bin/neopkpd trial --spec trial_spec.json --out trial_result.json
```

**Specification Format:**
```json
{
  "design": {
    "type": "parallel",
    "n_arms": 2,
    "randomization_ratio": [1, 1]
  },
  "arms": [
    {
      "name": "Placebo",
      "model": {"kind": "OneCompIVBolus", "params": {"CL": 5.0, "V": 50.0}},
      "regimen": {"type": "qd", "dose": 0.0, "duration_days": 28}
    },
    {
      "name": "Treatment",
      "model": {"kind": "OneCompIVBolus", "params": {"CL": 5.0, "V": 50.0}},
      "regimen": {"type": "qd", "dose": 100.0, "duration_days": 28}
    }
  ],
  "population": {
    "n": 100,
    "seed": 42
  },
  "endpoints": [
    {"type": "pk", "observation": "conc", "metric": "auc"}
  ],
  "n_replicates": 1
}
```

**Design Types:**
- `parallel` - Parallel group design
- `crossover` - Crossover design
- `dose_escalation` - 3+3 dose escalation
- `bioequivalence` - Bioequivalence study

### import

Import models from NONMEM or Monolix.

```bash
# Import NONMEM control file
./packages/cli/bin/neopkpd import --input run001.ctl --format nonmem --out model.json

# Import Monolix project
./packages/cli/bin/neopkpd import --input project.mlxtran --format monolix --out model.json
```

**Options:**
| Option | Required | Description |
|--------|----------|-------------|
| `--input` | Yes | Path to input model file |
| `--format` | Yes | Format: `nonmem` or `monolix` |
| `--out` | No | Output path for converted model |

### replay

Replay an execution artifact to verify reproducibility.

```bash
./packages/cli/bin/neopkpd replay --artifact simulation.json --out replayed.json
```

**Supported Artifact Types:**
- Single execution (`artifact_type: "single"`)
- Population execution (`artifact_type: "population"`)
- Sensitivity analysis (`artifact_type: "sensitivity_single"` or `"sensitivity_population"`)
- Estimation results (`artifact_type: "estimation"`)

### validate-golden

Run the full golden artifact validation suite.

```bash
./packages/cli/bin/neopkpd validate-golden
```

Validates all artifacts in `validation/golden/` directory against stored results.

### metrics

Compute PK/PD metrics from simulation results.

```bash
./packages/cli/bin/neopkpd metrics --artifact result.json --metrics cmax,tmax,auc
```

### sensitivity

Run parameter sensitivity analysis.

```bash
./packages/cli/bin/neopkpd sensitivity --spec sensitivity_spec.json --out sensitivity_result.json
```

**Specification Format:**
```json
{
  "model": {
    "kind": "OneCompIVBolus",
    "params": {"CL": 5.0, "V": 50.0},
    "doses": [{"time": 0.0, "amount": 100.0}]
  },
  "parameter": "CL",
  "perturbation": 0.1,
  "observation": "conc",
  "grid": {"t0": 0.0, "t1": 24.0, "saveat": [0, 1, 2, 4, 8, 12, 24]},
  "solver": {"alg": "Tsit5", "reltol": 1e-10, "abstol": 1e-12}
}
```

### help

Display help for any command.

```bash
./packages/cli/bin/neopkpd help
./packages/cli/bin/neopkpd help simulate
./packages/cli/bin/neopkpd help estimate
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error |
| 2 | Invalid arguments |
| 3 | File not found |
| 4 | Validation failure |

## Examples

### Run Population Simulation and Compute VPC

```bash
# Step 1: Run population simulation
./packages/cli/bin/neopkpd population --spec pop.json --out pop_result.json

# Step 2: Compute VPC using population results
./packages/cli/bin/neopkpd vpc --spec vpc.json --out vpc_result.json
```

### Import NONMEM Model and Estimate

```bash
# Step 1: Import NONMEM control file
./packages/cli/bin/neopkpd import --input run001.ctl --format nonmem --out model.json

# Step 2: Run estimation with imported model
./packages/cli/bin/neopkpd estimate --spec estimate.json --out fit.json
```

### CI/CD Integration

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
      - run: julia --project=packages/cli -e 'using Pkg; Pkg.instantiate()'
      - run: ./packages/cli/bin/neopkpd validate-golden
```

## Troubleshooting

### Julia Not Found

```
Error: julia: command not found
```

Ensure Julia is installed and in your PATH:
```bash
export PATH="$PATH:/path/to/julia/bin"
```

### Package Not Installed

```
Error: ArgumentError: Package NeoPKPD not found
```

Install dependencies:
```bash
julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'
julia --project=packages/cli -e 'using Pkg; Pkg.instantiate()'
```

## License

MIT License - see repository root for details.
