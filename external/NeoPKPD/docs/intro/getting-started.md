# Getting Started

This guide will help you install NeoPKPD and run your first simulation in under 10 minutes.

---

## Prerequisites

### System Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| **Julia** | 1.10+ | 1.11+ |
| **Python** | 3.10+ | 3.12+ |
| **Memory** | 4 GB | 8 GB+ |
| **OS** | macOS, Linux, Windows | macOS, Linux |

### Installing Julia

=== "macOS"

    ```bash
    # Using Homebrew
    brew install julia

    # Or download from julialang.org
    curl -fsSL https://install.julialang.org | sh
    ```

=== "Linux"

    ```bash
    # Using juliaup (recommended)
    curl -fsSL https://install.julialang.org | sh

    # Or using apt (Ubuntu/Debian)
    sudo apt install julia
    ```

=== "Windows"

    ```powershell
    # Using winget
    winget install julia

    # Or download installer from julialang.org
    ```

Verify installation:

```bash
julia --version
# Julia Version 1.11.0
```

### Installing Python

=== "macOS"

    ```bash
    # Using Homebrew
    brew install python@3.12

    # Or using pyenv
    pyenv install 3.12.0
    ```

=== "Linux"

    ```bash
    # Ubuntu/Debian
    sudo apt install python3.12 python3.12-venv

    # Or using pyenv
    pyenv install 3.12.0
    ```

=== "Windows"

    ```powershell
    # Using winget
    winget install Python.Python.3.12
    ```

---

## Installation

### Clone the Repository

```bash
git clone https://github.com/neopkpd/neopkpd.git
cd neopkpd
```

### Install Julia Core

```bash
# Activate and instantiate the Julia project
julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'
```

This installs all Julia dependencies including:
- DifferentialEquations.jl (ODE solvers)
- Distributions.jl (Statistical distributions)
- JSON3.jl (Serialization)

### Install Python Bindings (Optional)

```bash
cd packages/python

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install package with all dependencies
pip install -e ".[all]"
```

This installs:
- `neopkpd` - Core Python bindings
- `matplotlib`, `plotly` - Visualization backends
- `numpy`, `pandas` - Data manipulation
- `juliacall` - Julia-Python bridge

### Verify Installation

=== "Julia"

    ```bash
    julia --project=packages/core

    julia> using NeoPKPD
    julia> println(version())
    # 0.1.0
    ```

=== "Python"

    ```python
    import neopkpd
    neopkpd.init_julia()
    print(neopkpd.version())
    # 0.1.0
    ```

=== "CLI"

    ```bash
    ./packages/cli/bin/neopkpd version
    # NeoPKPD version 0.1.0
    ```

---

## Your First Simulation

### Julia: One-Compartment IV Bolus

```julia
using NeoPKPD

# Define model parameters
# CL = 5 L/h, V = 50 L
params = OneCompIVBolusParams(5.0, 50.0)

# Define dose: 100 mg at time 0
doses = [DoseEvent(0.0, 100.0)]

# Create model specification
spec = ModelSpec(OneCompIVBolus(), "my_first_sim", params, doses)

# Define simulation time grid (0 to 24 hours, hourly output)
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))

# Define solver settings
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run simulation
result = simulate(spec, grid, solver)

# Access results
println("Time points: ", result.t)
println("Concentrations: ", result.observations[:conc])
```

**Expected Output:**
```
Time points: [0.0, 1.0, 2.0, ..., 24.0]
Concentrations: [2.0, 1.81, 1.64, ..., 0.18]
```

### Python: One-Compartment IV Bolus

```python
import neopkpd

# Initialize Julia (required once per session)
neopkpd.init_julia()

# Run IV bolus simulation
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0,              # Clearance (L/h)
    v=50.0,              # Volume (L)
    doses=[{"time": 0.0, "amount": 100.0}],  # 100 mg at t=0
    t0=0.0,              # Start time
    t1=24.0,             # End time
    saveat=[float(t) for t in range(25)]  # Hourly output
)

# Access results
print("Time points:", result["t"])
print("Concentrations:", result["observations"]["conc"])
```

### CLI: Simulate from Spec File

Create a spec file `simulation.json`:

```json
{
  "model": {
    "kind": "OneCompIVBolus",
    "params": {"CL": 5.0, "V": 50.0},
    "doses": [{"time": 0.0, "amount": 100.0}]
  },
  "grid": {
    "t0": 0.0,
    "t1": 24.0,
    "saveat": [0, 1, 2, 3, 4, 6, 8, 12, 24]
  },
  "solver": {
    "alg": "Tsit5",
    "reltol": 1e-10,
    "abstol": 1e-12
  }
}
```

Run simulation:

```bash
./packages/cli/bin/neopkpd simulate --spec simulation.json --out result.json
```

---

## Your First Population Simulation

### Python: Population with IIV

```python
import neopkpd

neopkpd.init_julia()

# Simulate 100 subjects with inter-individual variability
result = neopkpd.simulate_population_iv_bolus(
    cl=5.0,              # Typical CL
    v=50.0,              # Typical V
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=100,               # Number of subjects
    seed=12345,          # For reproducibility
    omegas={             # IIV (CV)
        "CL": 0.3,       # 30% CV on CL
        "V": 0.2         # 20% CV on V
    }
)

# Access individual results
for i in range(3):
    cmax = max(result["individuals"][i]["observations"]["conc"])
    print(f"Subject {i+1} Cmax: {cmax:.2f} mg/L")

# Access population summary
summary = result["summaries"]["conc"]
print(f"\nPopulation median Cmax: {max(summary['median']):.2f} mg/L")
print(f"90% prediction interval: {max(summary['quantiles']['0.05']):.2f} - {max(summary['quantiles']['0.95']):.2f}")
```

---

## Your First Visualization

### Python: Concentration-Time Plot

```python
import neopkpd
from neopkpd import viz

neopkpd.init_julia()

# Run simulation
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)]
)

# Set visualization backend
viz.set_backend("matplotlib")

# Create concentration-time plot
fig = viz.plot_conc_time(result, title="One-Compartment IV Bolus")
fig.savefig("concentration_time.png", dpi=300)
```

### Python: Population Spaghetti Plot

```python
import neopkpd
from neopkpd import viz

neopkpd.init_julia()

# Run population simulation
pop_result = neopkpd.simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=50, seed=42,
    omegas={"CL": 0.3, "V": 0.2}
)

# Create spaghetti plot
viz.set_backend("matplotlib")
fig = viz.plot_spaghetti(pop_result, alpha=0.3)
fig.savefig("population_spaghetti.png", dpi=300)

# Create mean with confidence ribbon
fig = viz.plot_mean_ribbon(pop_result, ci_levels=[0.05, 0.95])
fig.savefig("population_ribbon.png", dpi=300)
```

---

## Next Steps

Now that you have NeoPKPD running, explore these topics:

### Learn the Basics

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } **Julia Tutorial**

    Complete walkthrough of Julia API

    [:octicons-arrow-right-24: Julia Tutorial](../julia/tutorial.md)

-   :material-language-python:{ .lg .middle } **Python Tutorial**

    Complete walkthrough of Python bindings

    [:octicons-arrow-right-24: Python Tutorial](../python/tutorial.md)

</div>

### Explore Models

<div class="grid cards" markdown>

-   :material-cube-outline:{ .lg .middle } **PK Models**

    One, two, three-compartment, transit, Michaelis-Menten

    [:octicons-arrow-right-24: PK Models](../julia/models/pk/onecomp-iv-bolus.md)

-   :material-chart-line:{ .lg .middle } **PD Models**

    Emax, sigmoid, effect compartment, indirect response

    [:octicons-arrow-right-24: PD Models](../julia/models/pd/direct-emax.md)

</div>

### Advanced Features

<div class="grid cards" markdown>

-   :material-account-group:{ .lg .middle } **Population Modeling**

    IIV, IOV, and covariate effects

    [:octicons-arrow-right-24: Population](../julia/population/index.md)

-   :material-chart-box:{ .lg .middle } **Parameter Estimation**

    FOCE-I, SAEM, Laplacian

    [:octicons-arrow-right-24: Estimation](../julia/estimation/index.md)

-   :material-flask:{ .lg .middle } **Clinical Trials**

    Trial simulation and power analysis

    [:octicons-arrow-right-24: Trials](../julia/trial/index.md)

-   :material-palette:{ .lg .middle } **Visualization**

    55+ plotting functions

    [:octicons-arrow-right-24: Visualization](../python/viz/index.md)

</div>

---

## Troubleshooting

### Julia Not Found

```
Error: Julia executable not found
```

**Solution:** Ensure Julia is in your PATH:

```bash
# Check Julia location
which julia

# Add to PATH if needed (add to ~/.bashrc or ~/.zshrc)
export PATH="$PATH:/path/to/julia/bin"
```

### Python Package Import Error

```
ModuleNotFoundError: No module named 'neopkpd'
```

**Solution:** Ensure virtual environment is activated and package is installed:

```bash
source packages/python/.venv/bin/activate
pip install -e ".[all]"
```

### Julia Initialization Slow

First Julia call takes 2-5 seconds due to JIT compilation. This is normal. Subsequent calls are fast (~1-10ms).

### Memory Issues with Large Populations

For populations >1000 subjects, consider:

- Using sparser `saveat` time points
- Running in batches
- Increasing system memory

---

## Getting Help

- **Documentation**: You're reading it!
- **Examples**: See `docs/examples/` for runnable code
- **Issues**: [GitHub Issues](https://github.com/neopkpd/neopkpd/issues)
- **Discussions**: [GitHub Discussions](https://github.com/neopkpd/neopkpd/discussions)
