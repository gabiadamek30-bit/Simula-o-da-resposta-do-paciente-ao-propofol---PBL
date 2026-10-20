# neopkpd

[![PyPI version](https://badge.fury.io/py/neopkpd.svg)](https://badge.fury.io/py/neopkpd)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**Transparent, validated pharmacokinetics and pharmacodynamics modeling infrastructure**

The `neopkpd` Python package provides a comprehensive interface to the NeoPKPD Julia simulation engine, enabling seamless integration with Python data science workflows.

## Installation

### From PyPI (Recommended)

```bash
pip install neopkpd

# With visualization support
pip install neopkpd[viz]

# With all optional dependencies
pip install neopkpd[all]
```

### From Source

```bash
git clone https://github.com/shramish2057/openpkpd.git
cd openpkpd/packages/python
pip install -e ".[all]"
```

**Requirements:**
- Python 3.10+
- Julia 1.10+ (automatically managed via juliacall)
- numpy, scipy (installed automatically)

## Quick Start

```python
import neopkpd

# Initialize Julia (required once per session)
neopkpd.init_julia()

# Run a simple simulation
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)]
)

print("Concentrations:", result["observations"]["conc"])
```

## Features

### PK/PD Simulation

```python
# One-compartment IV bolus
result = neopkpd.simulate_pk_iv_bolus(cl=5.0, v=50.0, ...)

# IV infusion (specify duration)
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0, "duration": 1.0}],  # 1-hour infusion
    ...
)

# One-compartment oral
result = neopkpd.simulate_pk_oral_first_order(ka=1.5, cl=5.0, v=50.0, ...)

# Two-compartment models
result = neopkpd.simulate_pk_twocomp_iv_bolus(cl=5.0, v1=50.0, q=2.0, v2=100.0, ...)
result = neopkpd.simulate_pk_twocomp_oral(ka=1.5, cl=5.0, v1=50.0, q=2.0, v2=100.0, ...)

# Three-compartment
result = neopkpd.simulate_pk_threecomp_iv_bolus(cl=5.0, v1=10.0, q2=20.0, v2=50.0, q3=2.0, v3=200.0, ...)

# Advanced models
result = neopkpd.simulate_pk_transit_absorption(ktr=2.0, n_transit=3, cl=5.0, v=50.0, ...)
result = neopkpd.simulate_pk_michaelis_menten(vmax=10.0, km=2.0, v=50.0, ...)

# PK-PD models
result = neopkpd.simulate_pkpd_direct_emax(cl=5.0, v=50.0, e0=0.0, emax=100.0, ec50=2.0, ...)
result = neopkpd.simulate_pkpd_sigmoid_emax(cl=5.0, v=50.0, e0=10.0, emax=40.0, ec50=0.8, gamma=2.0, ...)
result = neopkpd.simulate_pkpd_biophase_equilibration(cl=5.0, v=50.0, ke0=0.5, e0=0.0, emax=100.0, ec50=1.0, ...)
result = neopkpd.simulate_pkpd_indirect_response(cl=5.0, v=50.0, kin=10.0, kout=0.5, ic50=2.0, imax=0.9, ...)
```

### Population Simulation

```python
# Population with IIV
result = neopkpd.simulate_population_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)],
    n=100,
    seed=12345,
    omegas={"CL": 0.3, "V": 0.2}
)

# Access summaries
print("Mean:", result["summaries"]["conc"]["mean"])
print("5th percentile:", result["summaries"]["conc"]["quantiles"]["0.05"])
```

### Parameter Estimation (NLME)

```python
from neopkpd.estimation import estimate, EstimationConfig, FOCEIMethod

config = EstimationConfig(
    method=FOCEIMethod(max_inner_iter=100, inner_tol=1e-6),
    theta_init=[5.0, 50.0],
    theta_lower=[0.1, 1.0],
    theta_upper=[100.0, 500.0],
    omega_init=[[0.09, 0], [0, 0.04]],
    sigma_init={"kind": "proportional", "value": 0.1},
    max_iter=500,
    compute_se=True
)

result = estimate(observed_data, model_spec, config)

print("Theta:", result.theta)
print("Theta SE:", result.theta_se)
print("OFV:", result.ofv)
print("AIC:", result.aic)
```

### Non-Compartmental Analysis (NCA)

```python
from neopkpd.nca import run_nca, NCAConfig

result = run_nca(
    times=[0, 0.5, 1, 2, 4, 8, 12, 24],
    conc=[0, 1.8, 2.0, 1.5, 1.0, 0.5, 0.25, 0.06],
    dose=100.0,
    config=NCAConfig(method="log_linear")
)

print(f"Cmax: {result.cmax:.2f}")
print(f"Tmax: {result.tmax:.2f}")
print(f"AUC0-inf: {result.auc_0_inf:.2f}")
print(f"t1/2: {result.t_half:.2f}")
```

### Visual Predictive Check (VPC)

```python
from neopkpd.analysis import compute_vpc, VPCConfig

config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95],
    ci_level=0.95,
    n_bins=10,
    prediction_corrected=False,
    n_bootstrap=500,
    seed=12345
)

vpc_result = compute_vpc(observed_data, pop_spec, grid, solver, config=config)

# Visualize
from neopkpd import viz
viz.plot_vpc(vpc_result)
```

### Clinical Trial Simulation

```python
from neopkpd import trial

# Define design
design = trial.parallel_design(n_arms=2)

# Define dosing
regimen = trial.dosing_qd(dose=100.0, duration_days=28)

# Generate population
pop = trial.generate_virtual_population(n=100, seed=42)

# Create trial spec
spec = trial.TrialSpec(
    name="Phase 2",
    design=design,
    arms=[
        trial.TreatmentArm("Placebo", regimen=trial.dosing_qd(0.0, 28)),
        trial.TreatmentArm("Treatment", regimen=regimen),
    ],
    population=pop
)

# Run simulation
result = trial.simulate_trial(spec, seed=12345)

# Power analysis (multiple replicates)
power_result = trial.run_power_simulation(spec, n_replicates=100, seed=42)
print(f"Power: {power_result.power:.1%}")
```

### Model Import

```python
from neopkpd.import_ import parse_nonmem, parse_monolix

# Import NONMEM control file
model_spec, pop_spec, mapping = parse_nonmem("run001.ctl")

# Import Monolix project
model_spec, pop_spec, mapping = parse_monolix("project.mlxtran")
```

### CDISC Data Import

```python
from neopkpd.data import read_cdisc_csv, cdisc_to_population

# Read CDISC domains
dataset = read_cdisc_csv("pc.csv", "ex.csv", "dm.csv")

# Convert to NeoPKPD format
pop_spec, observed = cdisc_to_population(dataset, model_spec)
```

### Sensitivity Analysis

```python
result = neopkpd.run_sensitivity(
    model_kind="OneCompIVBolus",
    params={"CL": 5.0, "V": 50.0},
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)],
    sensitivity_param="CL",
    perturbation=0.1,
    observation="conc"
)

print("Max relative delta:", result.metrics.max_rel_delta)
```

### Visualization

```python
from neopkpd import viz

# Set backend
viz.set_backend("matplotlib")  # or "plotly"

# PK plots
viz.plot_conc_time(result)
viz.plot_spaghetti(pop_result, n_subjects=50)
viz.plot_mean_ribbon(pop_result, ci_levels=[0.05, 0.95])

# VPC
viz.plot_vpc(vpc_result)

# Estimation diagnostics
viz.plot_goodness_of_fit(estimation_result, plot_type="obs_vs_pred")
viz.plot_estimation_summary(estimation_result)  # 4-panel summary

# Sensitivity
viz.plot_sensitivity(sensitivity_result)
viz.plot_sensitivity_tornado(sensitivity_results)

# NCA
viz.plot_lambda_z_fit(nca_result, times, conc)
viz.plot_auc_visualization(times, conc, nca_result)

# Trial
viz.plot_power_curve(power_results, target_power=0.80)
viz.plot_tornado(sensitivity_analysis)
```

## Module Structure

```
neopkpd/
├── __init__.py           # Core simulation functions
├── _core.py              # Julia bridge utilities
├── simulations/          # PK/PD simulation wrappers
│   ├── pk_onecomp.py
│   ├── pk_twocomp.py
│   ├── pk_threecomp.py
│   ├── pk_advanced.py
│   └── pkpd.py
├── nca/                  # Non-compartmental analysis
│   ├── core.py
│   ├── config.py
│   └── bioequivalence.py
├── trial/                # Clinical trial simulation
│   ├── designs.py
│   ├── regimens.py
│   ├── population.py
│   └── power.py
├── estimation/           # NLME parameter estimation
│   ├── config.py
│   ├── nlme.py
│   └── diagnostics.py
├── import_/              # Model import
│   ├── nonmem.py
│   └── monolix.py
├── data/                 # Data handling
│   └── cdisc.py
├── analysis/             # VPC, sensitivity
│   ├── vpc.py
│   └── sensitivity.py
└── viz/                  # Visualization
    ├── pk.py
    ├── nca.py
    ├── pkpd.py
    ├── population.py
    ├── trial.py
    └── backends.py
```

## Integration with Data Science Tools

### NumPy

```python
import numpy as np
import neopkpd

result = neopkpd.simulate_pk_iv_bolus(...)
t = np.array(result["t"])
conc = np.array(result["observations"]["conc"])

print(f"Cmax: {np.max(conc):.2f}")
print(f"AUC: {np.trapz(conc, t):.2f}")
```

### Pandas

```python
import pandas as pd

# Population results to DataFrame
params_df = pd.DataFrame(pop_result["params"])
print(params_df.describe())

# Concentration data
data = []
for i, ind in enumerate(pop_result["individuals"]):
    for t, c in zip(ind["t"], ind["observations"]["conc"]):
        data.append({"id": i, "time": t, "conc": c})
conc_df = pd.DataFrame(data)
```

## Performance Tips

1. **Initialize once**: Call `init_julia()` once at startup
2. **Batch simulations**: Julia JIT makes subsequent calls faster
3. **Use appropriate tolerances**: Higher tolerance = faster, less accurate
4. **Sparse output**: Use fewer `saveat` points for large populations

## Testing

```bash
cd packages/python
source .venv/bin/activate
pytest tests/
```

## Documentation

Full documentation: [shramish2057.github.io/NeoPKPD](https://shramish2057.github.io/NeoPKPD/)

## License

MIT License - see [LICENSE](https://github.com/shramish2057/openpkpd/blob/main/LICENSE) for details.

## Citation

```bibtex
@software{neopkpd,
  title = {NeoPKPD: Transparent PK/PD Modeling Infrastructure},
  url = {https://github.com/shramish2057/openpkpd},
  version = {0.1.0}
}
```
