<p align="center">
  <h1 align="center">NeoPKPD</h1>
  <p align="center">
    <strong>Transparent, validated pharmacokinetics and pharmacodynamics modeling infrastructure</strong>
  </p>
  <p align="center">
    <a href="https://doi.org/10.5281/zenodo.18215969"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.18215969.svg" alt="DOI"></a>
    <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT"></a>
    <a href="https://julialang.org/"><img src="https://img.shields.io/badge/Julia-1.10+-purple.svg" alt="Julia"></a>
    <a href="https://python.org/"><img src="https://img.shields.io/badge/Python-3.10+-blue.svg" alt="Python"></a>
  </p>
  <p align="center">
    <a href="#features">Features</a> •
    <a href="#installation">Installation</a> •
    <a href="#quick-start">Quick Start</a> •
    <a href="#documentation">Documentation</a> •
    <a href="#contributing">Contributing</a>
  </p>
</p>

---

## Overview

**NeoPKPD** is a reference-grade PK/PD simulation platform designed for research, method development, and reproducible scientific computation. It emphasizes deterministic execution, transparent numerical semantics, and complete artifact serialization.

## Features

| Category | Features |
|----------|----------|
| **PK Models** | One/Two/Three-compartment IV & oral, transit absorption, Michaelis-Menten elimination |
| **PD Models** | Direct Emax, sigmoid Emax, biophase equilibration, indirect response (IRM1-4), transit compartment PD, disease progression, tolerance/counter-regulation |
| **TMDD Models** | Full TMDD, QSS approximation, Michaelis-Menten approximation for biologics |
| **Drug Interactions** | Bliss independence, competitive inhibition, receptor regulation |
| **IV Infusion** | Zero-order infusion with duration support, overlapping infusions |
| **Population** | IIV, IOV, static & time-varying covariates, BLQ handling |
| **Parameter Estimation** | NLME with FOCE-I, SAEM, Laplacian, and Bayesian methods |
| **Advanced Estimation** | Bootstrap, mixture models, model averaging, stepwise covariate modeling (SCM) |
| **NCA** | FDA/EMA-compliant non-compartmental analysis |
| **Trial Simulation** | Parallel, crossover, dose-escalation, bioequivalence, adaptive designs |
| **Adaptive Trials** | Response-adaptive randomization, sample size re-estimation, treatment selection, enrichment |
| **Sensitivity Analysis** | Local sensitivity, Sobol' indices (first/second/total order), Morris screening |
| **VPC & Diagnostics** | Visual Predictive Checks (pcVPC), CWRES, IWRES, NPDE |
| **Model Import** | NONMEM (.ctl) and Monolix (.mlxtran) model parsing |
| **Data Import** | CDISC/SDTM format support (PC, EX, DM domains), XPT reader |
| **Residual Error** | Additive, proportional, combined, exponential models |
| **Compliance** | FDA 21 CFR Part 11 support, digital signatures, audit trails |
| **Visualization** | Matplotlib/Plotly backends with estimation diagnostics |
| **Interfaces** | Julia API, Python bindings, CLI |
| **Reproducibility** | Versioned artifacts with deterministic replay, schema validation |

## Installation

### Julia (Core)

```bash
git clone https://github.com/shramish2057/neopkpd.git
cd neopkpd

# Install dependencies
julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'

# Verify installation
julia --project=packages/core -e 'using NeoPKPD; println("v", NEOPKPD_VERSION)'
```

### Python (Optional)

```bash
cd packages/python
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

### CLI

```bash
./packages/cli/bin/neopkpd version
```

## Quick Start

### Julia

```julia
using NeoPKPD

# Define model
spec = ModelSpec(
    OneCompIVBolus(),
    "example",
    OneCompIVBolusParams(5.0, 50.0),  # CL=5 L/h, V=50 L
    [DoseEvent(0.0, 100.0)]            # 100 mg at t=0
)

# Configure simulation
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run and access results
result = simulate(spec, grid, solver)
println(result.observations[:conc])
```

### Python

```python
import neopkpd

neopkpd.init_julia()
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)]
)
print(result["observations"]["conc"])
```

### CLI

```bash
# Check version
./packages/cli/bin/neopkpd version

# Run simulation from JSON spec
./packages/cli/bin/neopkpd simulate --spec simulation.json --out result.json

# Run NCA analysis
./packages/cli/bin/neopkpd nca --spec nca_spec.json --out nca_result.json

# Run clinical trial simulation
./packages/cli/bin/neopkpd trial --spec trial_spec.json --out trial_result.json

# Compute VPC
./packages/cli/bin/neopkpd vpc --spec vpc_spec.json --out vpc_result.json

# Import NONMEM model
./packages/cli/bin/neopkpd import --input model.ctl --format nonmem --out model.json

# Replay an artifact
./packages/cli/bin/neopkpd replay --artifact validation/golden/pk_iv_bolus.json

# Validate all golden artifacts
./packages/cli/bin/neopkpd validate-golden
```

## Repository Structure

```
neopkpd/
├── packages/
│   ├── core/                 # Julia simulation engine (NeoPKPD)
│   │   ├── src/
│   │   │   ├── pk/           # PK model definitions (1/2/3-comp, transit, MM)
│   │   │   ├── pd/           # PD models (Emax, indirect response, etc.)
│   │   │   ├── tmdd/         # TMDD models for biologics
│   │   │   ├── engine/       # ODE solving, population, sensitivity
│   │   │   ├── specs/        # Type specifications (ModelSpec, etc.)
│   │   │   ├── estimation/   # NLME: FOCE-I, SAEM, Laplacian, bootstrap
│   │   │   ├── trial/        # Clinical trial simulation & adaptive designs
│   │   │   ├── nca/          # Non-compartmental analysis
│   │   │   ├── import/       # NONMEM/Monolix parsers
│   │   │   ├── data/         # CDISC data handling
│   │   │   ├── analysis/     # VPC, NPDE, residuals
│   │   │   ├── compliance/   # FDA 21 CFR Part 11 compliance
│   │   │   └── serialization/# JSON artifact I/O with schema validation
│   │   └── test/             # Comprehensive test suite (5400+ tests)
│   ├── python/               # Python bindings (neopkpd)
│   │   └── neopkpd/
│   │       ├── simulations/  # PK/PD simulation wrappers
│   │       ├── nca/          # Non-compartmental analysis
│   │       ├── trial/        # Clinical trial simulation
│   │       ├── estimation/   # NLME Python interface
│   │       ├── vpc/          # Visual Predictive Checks
│   │       ├── import_/      # Model import utilities
│   │       ├── data/         # CDISC data utilities
│   │       └── viz/          # Visualization (matplotlib/plotly)
│   └── cli/                  # Command-line interface (NeoPKPDCLI)
│       ├── src/              # CLI commands
│       └── bin/              # Entry point script
├── validation/               # Golden artifacts & validation
│   ├── golden/               # Reference outputs (deterministic)
│   └── scripts/              # Validation runners
├── docs/                     # Documentation (MkDocs)
│   └── *.md                  # Architecture, semantics, CLI reference
└── scripts/                  # Development tools
```

## Testing

```bash
# Julia unit tests
julia --project=packages/core -e 'using Pkg; Pkg.test()'

# Golden artifact validation
./packages/cli/bin/neopkpd validate-golden

# Python tests
cd packages/python && source .venv/bin/activate && pytest tests/

# Documentation build
mkdocs build --strict
```

## Semantic Versioning

NeoPKPD uses independent version numbers for numerical behavior:

| Version | Current | Scope |
|---------|---------|-------|
| Event Semantics | 1.0.0 | Dose handling |
| Solver Semantics | 1.0.0 | ODE solver behavior |
| Artifact Schema | 1.0.0 | JSON format |

Any change to numerical output requires a version bump.

## Documentation

Full documentation: [shramish2057.github.io/NeoPKPD](https://shramish2057.github.io/NeoPKPD/)

| Guide | Description |
|-------|-------------|
| [Getting Started](docs/index.md) | Overview and quick start guide |
| [Architecture](docs/architecture.md) | System design and module structure |
| [Sensitivity Analysis](docs/sensitivity.md) | Local and global sensitivity methods |
| [Reproducibility](docs/reproducibility.md) | Artifact versioning and deterministic replay |
| [Semantics](docs/semantics.md) | Numerical semantics and versioning |
| [CLI Reference](docs/cli.md) | Command-line interface |

Build locally:
```bash
pip install -r docs/requirements.txt
mkdocs serve
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

```bibtex
@software{neopkpd,
  title = {NeoPKPD: Transparent PK/PD Modeling Infrastructure},
  url = {https://github.com/shramish2057/neopkpd},
  version = {0.1.0}
}
```

---

<p align="center">
  <sub>Built for reproducible pharmacometric research</sub>
</p>
