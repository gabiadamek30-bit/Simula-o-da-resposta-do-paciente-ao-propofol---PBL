# Python Documentation

Welcome to the Python documentation for NeoPKPD. The `neopkpd` package provides Pythonic access to all core functionality with seamless NumPy/Pandas integration.

---

## Why Python?

While NeoPKPD's core is written in Julia for performance, Python bindings enable:

- **Data Science Integration** - NumPy, Pandas, SciPy ecosystem
- **Visualization** - Matplotlib and Plotly dual backends
- **Accessibility** - Familiar syntax for most scientists
- **Reproducibility** - Jupyter notebook workflows

---

## Quick Start

```python
import neopkpd

# Initialize Julia (required once per session)
neopkpd.init_julia()

# One-compartment IV bolus simulation
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0,              # Clearance (L/h)
    v=50.0,              # Volume (L)
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[float(t) for t in range(25)]
)

print("Concentrations:", result["observations"]["conc"])
```

---

## Documentation Sections

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } **Tutorial**

    ---

    Step-by-step introduction to the Python API

    [:octicons-arrow-right-24: Start Tutorial](tutorial.md)

-   :material-cube-outline:{ .lg .middle } **Models**

    ---

    PK and PD simulation functions

    [:octicons-arrow-right-24: Models Reference](models/index.md)

-   :material-account-group:{ .lg .middle } **Population Modeling**

    ---

    Population simulation with IIV

    [:octicons-arrow-right-24: Population](population/index.md)

-   :material-calculator:{ .lg .middle } **NCA**

    ---

    Non-compartmental analysis

    [:octicons-arrow-right-24: NCA Reference](nca/index.md)

-   :material-chart-box:{ .lg .middle } **Parameter Estimation**

    ---

    NLME estimation methods

    [:octicons-arrow-right-24: Estimation](estimation/index.md)

-   :material-flask:{ .lg .middle } **Clinical Trials**

    ---

    Trial simulation and power analysis

    [:octicons-arrow-right-24: Trial Reference](trial/index.md)

-   :material-palette:{ .lg .middle } **Visualization**

    ---

    55+ plotting functions with dual backends

    [:octicons-arrow-right-24: Visualization](viz/index.md)

-   :material-database-import:{ .lg .middle } **Data Import**

    ---

    CDISC and CSV data handling

    [:octicons-arrow-right-24: Data Import](data/index.md)

</div>

---

## Installation

```bash
cd packages/python

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install with all dependencies
pip install -e ".[all]"
```

### Optional Dependencies

| Extra | Packages | Install Command |
|-------|----------|-----------------|
| `viz` | matplotlib, plotly | `pip install -e ".[viz]"` |
| `data` | pandas, pyarrow | `pip install -e ".[data]"` |
| `all` | All optional deps | `pip install -e ".[all]"` |

---

## Module Overview

### Core Module (`neopkpd`)

```python
import neopkpd

neopkpd.init_julia()              # Initialize Julia runtime
neopkpd.version()                 # Get version string

# Simulation functions
neopkpd.simulate_pk_iv_bolus(...)
neopkpd.simulate_pk_oral_first_order(...)
neopkpd.simulate_population_iv_bolus(...)
```

### NCA Module (`neopkpd.nca`)

```python
from neopkpd import nca

result = nca.run_nca(times, conc, dose)
summary = nca.summarize_population_nca(pop_results)
be_result = nca.bioequivalence_90ci(test, reference)
```

### Trial Module (`neopkpd.trial`)

```python
from neopkpd import trial

design = trial.parallel_design(n_arms=2)
regimen = trial.dosing_qd(dose=100.0, duration_days=28)
pop = trial.generate_virtual_population(n=100)
```

### Visualization Module (`neopkpd.viz`)

```python
from neopkpd import viz

viz.set_backend("matplotlib")
fig = viz.plot_conc_time(result)
fig = viz.plot_vpc_detailed(vpc_result)
```

---

## Result Structure

All simulation functions return dictionaries with consistent structure:

```python
result = {
    "t": [0.0, 1.0, 2.0, ...],           # Time points
    "states": {
        "A_central": [100.0, 90.5, ...]   # State variables
    },
    "observations": {
        "conc": [2.0, 1.81, ...]          # Concentrations
    },
    "metadata": {
        "model": "OneCompIVBolus",
        "version": "0.1.0"
    }
}
```

### Population Results

```python
pop_result = {
    "individuals": [                       # List of individual results
        {"t": [...], "observations": {...}},
        ...
    ],
    "params": [                            # Realized parameters
        {"CL": 4.8, "V": 52.1},
        ...
    ],
    "summaries": {                         # Population statistics
        "conc": {
            "mean": [...],
            "median": [...],
            "quantiles": {
                "0.05": [...],
                "0.95": [...]
            }
        }
    }
}
```

---

## Integration Examples

### With NumPy

```python
import numpy as np
import neopkpd

neopkpd.init_julia()

result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=list(np.linspace(0, 24, 49))
)

conc = np.array(result["observations"]["conc"])
print(f"Cmax: {np.max(conc):.2f}")
print(f"AUC: {np.trapz(conc, result['t']):.2f}")
```

### With Pandas

```python
import pandas as pd
import neopkpd

neopkpd.init_julia()

result = neopkpd.simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=100, seed=12345,
    omegas={"CL": 0.3, "V": 0.2}
)

# Parameters DataFrame
params_df = pd.DataFrame(result["params"])
print(params_df.describe())

# Concentration DataFrame
data = []
for i, ind in enumerate(result["individuals"]):
    for t, c in zip(ind["t"], ind["observations"]["conc"]):
        data.append({"id": i, "time": t, "conc": c})

conc_df = pd.DataFrame(data)
print(conc_df.groupby("time")["conc"].describe())
```

---

## Performance Tips

1. **Initialize Once**: Call `init_julia()` once at session start
2. **Batch Simulations**: Julia JIT makes subsequent calls faster
3. **Sparse Output**: Use fewer `saveat` points for large populations
4. **Parallel Populations**: Population simulations are internally parallelized

```python
# Good: Initialize once
import neopkpd
neopkpd.init_julia()

for params in parameter_sets:
    result = neopkpd.simulate_pk_iv_bolus(...)
```

---

## Next Steps

- [Start the Tutorial →](tutorial.md)
- [Explore Simulation Functions →](models/index.md)
- [Learn Visualization →](viz/index.md)
