# NeoPKPD Benchmark Suite

Comprehensive benchmarks for NeoPKPD performance evaluation against alternative pharmacometric platforms.

## Overview

This benchmark suite provides:
- **Reproducible** timing measurements with statistical analysis
- **Cross-platform** comparison (Julia, Python, R packages)
- **Publication-ready** figures and tables for SoftwareX paper

## Directory Structure

```
benchmarks/
├── README.md                 # This file
├── scripts/
│   ├── benchmark_neopkpd.jl        # Julia core benchmarks
│   ├── benchmark_neopkpd_python.py # Python bindings benchmarks
│   ├── benchmark_mrgsolve.R        # mrgsolve comparison
│   ├── benchmark_nlmixr2.R         # nlmixr2 comparison
│   └── plot_benchmarks.py          # Generate figures
├── results/                  # CSV output files (generated)
│   ├── neopkpd_julia_benchmarks.csv
│   ├── neopkpd_python_benchmarks.csv
│   ├── mrgsolve_benchmarks.csv
│   └── nlmixr2_benchmarks.csv
└── figures/                  # Publication figures (generated)
    ├── figure_simulation_comparison.png
    ├── figure_population_scaling.png
    ├── figure_speedup_summary.png
    ├── benchmark_summary_table.csv
    └── benchmark_summary_table.tex
```

## Quick Start

### 1. Run All Benchmarks

```bash
# From repository root
cd packages/core

# Run NeoPKPD Julia benchmarks
julia --project=. benchmarks/scripts/benchmark_neopkpd.jl

# Run NeoPKPD Python benchmarks
cd ../python
pip install -e .
python ../core/benchmarks/scripts/benchmark_neopkpd_python.py

# Run R comparison benchmarks (requires R packages)
cd ../core
Rscript benchmarks/scripts/benchmark_mrgsolve.R
Rscript benchmarks/scripts/benchmark_nlmixr2.R

# Generate figures
python benchmarks/scripts/plot_benchmarks.py
```

### 2. Quick Julia-Only Benchmark

```bash
julia --project=packages/core -e '
    include("packages/core/benchmarks/scripts/benchmark_neopkpd.jl")
    run_all_benchmarks()
'
```

## Prerequisites

### Julia (Required)
- Julia 1.10+
- NeoPKPD package (this repository)

### Python (Optional, for Python benchmarks)
```bash
pip install neopkpd matplotlib pandas numpy
```

### R (Optional, for comparison benchmarks)
```R
install.packages(c("mrgsolve", "nlmixr2", "rxode2", "dplyr", "microbenchmark"))
```

## Benchmark Categories

### 1. Single Simulation Speed
Measures time to simulate a single subject for various model types:
- One-compartment IV bolus
- Two-compartment IV bolus
- Two-compartment oral
- Transit absorption
- Michaelis-Menten elimination

### 2. Population Simulation Scaling
Measures how simulation time scales with population size:
- N = 50, 100, 500, 1000 subjects
- Two-compartment model with IIV (CL: 30% CV, V: 20% CV)

### 3. PKPD Coupled Simulation
Measures PK-PD model simulation performance:
- Direct Emax
- Sigmoid Emax
- Indirect response (IRM3)

### 4. NCA Computation
Measures non-compartmental analysis speed:
- Cmax, Tmax calculation
- AUC computation
- Lambda-z estimation
- Half-life calculation

### 5. Sensitivity Analysis
Measures sensitivity analysis performance:
- Local (single parameter) perturbation
- Multi-parameter perturbation

### 6. Serialization
Measures artifact I/O performance:
- JSON serialization
- JSON deserialization

## Configuration

Each benchmark script has configurable parameters:

```julia
# Julia (benchmark_neopkpd.jl)
const CONFIG = (
    n_runs = 100,           # Number of timed runs per benchmark
    n_warmup = 5,           # Warmup runs (JIT compilation)
    random_seed = 12345,    # For reproducibility
    output_dir = "benchmarks/results",
)
```

```python
# Python (benchmark_neopkpd_python.py)
CONFIG = {
    "n_runs": 100,
    "n_warmup": 5,
    "random_seed": 12345,
    "output_dir": "benchmarks/results",
}
```

## Output Format

All benchmark scripts produce CSV files with the following columns:

| Column | Description |
|--------|-------------|
| `category` | Benchmark category (simulation, population, pkpd, nca, etc.) |
| `name` | Specific benchmark name |
| `model` | Model type used |
| `n_subjects` | Number of subjects (for population benchmarks) |
| `n_runs` | Number of timing runs |
| `mean_ms` | Mean execution time in milliseconds |
| `std_ms` | Standard deviation in milliseconds |
| `median_ms` | Median execution time |
| `min_ms` | Minimum execution time |
| `max_ms` | Maximum execution time |
| `ci_lower_ms` | 2.5th percentile (95% CI lower bound) |
| `ci_upper_ms` | 97.5th percentile (95% CI upper bound) |
| `timestamp` | When benchmark was run |
| `*_version` | Software version information |

## Generating Publication Figures

After running benchmarks, generate figures:

```bash
python benchmarks/scripts/plot_benchmarks.py
```

This produces:

1. **figure_simulation_comparison.png**: Bar chart comparing single simulation times
2. **figure_population_scaling.png**: Line plot showing scaling behavior
3. **figure_speedup_summary.png**: Summary of speedup factors
4. **benchmark_summary_table.tex**: LaTeX table for paper

## Hardware Specifications

For reproducible benchmarks, document your hardware:

```
Processor: [e.g., AMD Ryzen 9 5900X, Intel Core i9-12900K]
Cores/Threads: [e.g., 12/24]
RAM: [e.g., 64 GB DDR4-3200]
OS: [e.g., Ubuntu 22.04, macOS 14.0, Windows 11]
Julia Version: [e.g., 1.10.0]
Python Version: [e.g., 3.11.0]
R Version: [e.g., 4.3.2]
```

## Troubleshooting

### Julia JIT Warmup
Julia's first run includes JIT compilation overhead. The benchmark automatically
runs warmup iterations before timing. For accurate results, ensure at least
`n_warmup = 5`.

### Python-Julia Bridge
Python benchmarks include juliacall bridge overhead. This is intentional to
measure real-world Python API performance. For pure Julia performance, use
the Julia benchmarks.

### R Package Installation
If nlmixr2 installation fails:
```R
# Try installing dependencies first
install.packages("rxode2")
install.packages("nlmixr2", dependencies = TRUE)
```

### Memory Issues
For large population benchmarks (N > 1000), ensure sufficient RAM:
```bash
# Increase Julia memory limit
export JULIA_NUM_THREADS=4
julia --heap-size-hint=8G benchmarks/scripts/benchmark_neopkpd.jl
```

## Contributing

To add new benchmarks:

1. Add benchmark function in appropriate script
2. Follow naming convention: `benchmark_<category>()`
3. Return `BenchmarkResult` objects
4. Update this README

## Citation

If you use these benchmarks in your research, please cite:

```bibtex
@article{kafle2025neopkpd,
  title={NeoPKPD: A Transparent, Reproducible Infrastructure for
         Pharmacometric Modeling and Clinical Trial Simulation},
  author={Kafle, Shramish and Kafle, Sakrit},
  journal={SoftwareX},
  year={2025}
}
```

## License

MIT License - see repository root for details.
