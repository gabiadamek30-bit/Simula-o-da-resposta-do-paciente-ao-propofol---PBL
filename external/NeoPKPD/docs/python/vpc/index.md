# Visual Predictive Check (VPC)

Comprehensive guide to VPC computation and analysis in Python.

---

## Overview

Visual Predictive Checks (VPCs) compare observed data distributions to model simulations, providing a visual diagnostic for population model adequacy.

```python
from neopkpd import compute_vpc, compute_pcvpc, compute_stratified_vpc
from neopkpd.vpc import VPCConfig, QuantileBinning

# Basic VPC computation
vpc_result = compute_vpc(
    observed_data=observed,
    population_spec=pop_spec,
    grid={"t0": 0.0, "t1": 24.0, "saveat": 0.5},
    config=VPCConfig(n_simulations=500)
)
```

---

## VPC Types

<div class="grid cards" markdown>

-   :material-chart-scatter-plot:{ .lg .middle } **Standard VPC**

    ---

    Compare observed vs simulated percentiles

    [:octicons-arrow-right-24: Standard VPC](#standard-vpc)

-   :material-chart-timeline:{ .lg .middle } **Prediction-Corrected VPC**

    ---

    Normalize for variable dosing and covariates

    [:octicons-arrow-right-24: pcVPC](#prediction-corrected-vpc)

-   :material-view-split-vertical:{ .lg .middle } **Stratified VPC**

    ---

    Separate VPC by covariate subgroups

    [:octicons-arrow-right-24: Stratified](#stratified-vpc)

-   :material-chart-box:{ .lg .middle } **VPC with BLQ**

    ---

    Handle below-quantification data

    [:octicons-arrow-right-24: BLQ Handling](blq.md)

</div>

---

## VPCConfig

Configuration dataclass for VPC computation:

```python
from neopkpd.vpc import VPCConfig, QuantileBinning

@dataclass
class VPCConfig:
    """VPC configuration parameters."""
    pi_levels: list[float] = field(default_factory=lambda: [0.05, 0.50, 0.95])
    ci_level: float = 0.95
    binning: BinningStrategy = field(default_factory=lambda: QuantileBinning(10))
    prediction_corrected: bool = False
    stratify_by: list[str] = field(default_factory=list)
    lloq: float | None = None
    n_simulations: int = 200
    n_bootstrap: int = 500
    seed: int = 12345
```

### Configuration Examples

```python
# Default configuration
config = VPCConfig()

# Custom configuration
config = VPCConfig(
    pi_levels=[0.10, 0.50, 0.90],     # 10th, 50th, 90th percentiles
    ci_level=0.90,                     # 90% confidence interval
    binning=QuantileBinning(8),        # 8 bins with equal obs
    n_simulations=500,                 # 500 simulation replicates
    n_bootstrap=1000,                  # 1000 bootstrap samples
    seed=42                            # Reproducibility
)

# pcVPC configuration
config = VPCConfig(
    prediction_corrected=True,
    n_simulations=500
)

# Stratified VPC configuration
config = VPCConfig(
    stratify_by=["DOSE_GROUP"],
    n_simulations=500
)

# VPC with BLQ
config = VPCConfig(
    lloq=0.1,                         # Lower limit of quantification
    n_simulations=500
)
```

### Parameter Details

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pi_levels` | `list[float]` | `[0.05, 0.50, 0.95]` | Prediction interval percentile levels |
| `ci_level` | `float` | `0.95` | Confidence interval level |
| `binning` | `BinningStrategy` | `QuantileBinning(10)` | Binning strategy |
| `prediction_corrected` | `bool` | `False` | Enable prediction correction |
| `stratify_by` | `list[str]` | `[]` | Stratification variables |
| `lloq` | `float | None` | `None` | Lower limit of quantification |
| `n_simulations` | `int` | `200` | Number of simulation replicates |
| `n_bootstrap` | `int` | `500` | Bootstrap samples for CI |
| `seed` | `int` | `12345` | Random seed |

---

## Standard VPC

### compute_vpc

```python
def compute_vpc(
    observed_data: dict[str, Any],
    population_spec: dict[str, Any],
    grid: dict[str, Any],
    config: VPCConfig | None = None,
    solver: dict[str, Any] | None = None,
    error_spec: dict[str, Any] | None = None
) -> VPCResult:
    """
    Compute standard Visual Predictive Check.

    Parameters
    ----------
    observed_data : dict
        Observed data with keys: subject_ids, times, dv, dvid
    population_spec : dict
        Population model specification
    grid : dict
        Simulation grid {t0, t1, saveat}
    config : VPCConfig, optional
        VPC configuration
    solver : dict, optional
        ODE solver specification
    error_spec : dict, optional
        Residual error model

    Returns
    -------
    VPCResult
        VPC computation results
    """
```

### Usage

```python
import neopkpd
from neopkpd.vpc import VPCConfig, QuantileBinning

neopkpd.init_julia()

# Observed data structure
observed_data = {
    "subject_ids": ["S1", "S1", "S1", "S2", "S2", "S2", ...],
    "times": [0.5, 2.0, 8.0, 0.5, 2.0, 8.0, ...],
    "dv": [10.5, 5.2, 1.1, 12.3, 6.1, 1.3, ...],
    "dvid": ["conc"] * len(times)
}

# Population model
population_spec = {
    "model": "OneCompOral",
    "params": {"Ka": 1.5, "CL": 5.0, "V": 50.0},
    "omega": {"Ka": 0.16, "CL": 0.09, "V": 0.04},
    "doses": [{"time": 0.0, "amount": 100.0}],
    "n": 50  # Number of subjects
}

# Simulation grid
grid = {"t0": 0.0, "t1": 24.0, "saveat": 0.5}

# VPC configuration
config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95],
    binning=QuantileBinning(7),
    n_simulations=500,
    seed=42
)

# Compute VPC
vpc_result = neopkpd.compute_vpc(
    observed_data=observed_data,
    population_spec=population_spec,
    grid=grid,
    config=config
)

# Access results
print(f"Bins: {len(vpc_result.bins)}")
print(f"Simulations: {vpc_result.n_simulations}")
```

### With Residual Error

```python
# Add residual error model
error_spec = {
    "kind": "combined",
    "sigma_add": 0.1,
    "sigma_prop": 0.1
}

vpc_result = neopkpd.compute_vpc(
    observed_data=observed_data,
    population_spec=population_spec,
    grid=grid,
    config=config,
    error_spec=error_spec
)
```

---

## Prediction-Corrected VPC

### compute_pcvpc

```python
def compute_pcvpc(
    observed_data: dict[str, Any],
    population_spec: dict[str, Any],
    grid: dict[str, Any],
    config: VPCConfig | None = None,
    solver: dict[str, Any] | None = None,
    error_spec: dict[str, Any] | None = None
) -> VPCResult:
    """
    Compute prediction-corrected Visual Predictive Check.

    Normalizes observations by population prediction:
        pcDV = DV × (PRED_bin / PRED_individual)

    Use when:
    - Variable dosing across subjects
    - Significant covariate effects
    - Dose escalation studies
    """
```

### Usage

```python
from neopkpd.vpc import VPCConfig

# pcVPC for variable-dose study
config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95],
    n_simulations=500
)

pcvpc_result = neopkpd.compute_pcvpc(
    observed_data=observed_data,
    population_spec=population_spec,
    grid=grid,
    config=config
)

# Results are prediction-corrected
# Y-axis represents deviation from population prediction
```

### When to Use pcVPC

| Scenario | Use pcVPC? | Reason |
|----------|------------|--------|
| Fixed dose study | No | Standard VPC sufficient |
| Variable dosing | Yes | Normalize dose differences |
| Weight-based dosing | Yes | Normalize covariate effects |
| Dose escalation | Yes | Compare across dose levels |
| Multiple formulations | Consider | If bioavailability differs |

---

## Stratified VPC

### compute_stratified_vpc

```python
def compute_stratified_vpc(
    observed_data: dict[str, Any],
    population_spec: dict[str, Any],
    grid: dict[str, Any],
    strata_data: dict[str, dict[str, Any]],
    config: VPCConfig | None = None,
    solver: dict[str, Any] | None = None,
    error_spec: dict[str, Any] | None = None
) -> StratifiedVPCResult:
    """
    Compute stratified Visual Predictive Check.

    Computes separate VPC for each stratum defined by
    stratification variables.

    Parameters
    ----------
    strata_data : dict
        Maps subject ID to stratum values
        {"S1": {"DOSE": "100mg"}, "S2": {"DOSE": "200mg"}, ...}
    """
```

### Usage

```python
from neopkpd.vpc import VPCConfig

# Define strata for each subject
strata_data = {
    "S1": {"DOSE_GROUP": "Low"},
    "S2": {"DOSE_GROUP": "Low"},
    "S3": {"DOSE_GROUP": "Medium"},
    "S4": {"DOSE_GROUP": "Medium"},
    "S5": {"DOSE_GROUP": "High"},
    "S6": {"DOSE_GROUP": "High"},
    # ...
}

# Configuration with stratification
config = VPCConfig(
    stratify_by=["DOSE_GROUP"],
    pi_levels=[0.05, 0.50, 0.95],
    n_simulations=500
)

# Compute stratified VPC
stratified_result = neopkpd.compute_stratified_vpc(
    observed_data=observed_data,
    population_spec=population_spec,
    grid=grid,
    strata_data=strata_data,
    config=config
)

# Access results by stratum
for i, vpc in enumerate(stratified_result.results):
    stratum = stratified_result.strata_names[i]
    print(f"{stratum}: {vpc.n_subjects_observed} subjects")
```

### Multiple Stratification Variables

```python
# Stratify by dose AND formulation
strata_data = {
    "S1": {"DOSE": "100mg", "FORM": "tablet"},
    "S2": {"DOSE": "100mg", "FORM": "capsule"},
    "S3": {"DOSE": "200mg", "FORM": "tablet"},
    # ...
}

config = VPCConfig(
    stratify_by=["DOSE", "FORM"],  # Creates combined strata
    n_simulations=500
)
```

---

## VPCResult Structure

```python
@dataclass
class VPCResult:
    """Complete VPC computation result."""
    config: VPCConfig                    # Configuration used
    bins: list[VPCBin]                   # Computed bins
    n_subjects_observed: int             # Number of subjects
    n_observations_observed: int         # Total observations
    n_simulations: int                   # Simulations performed
    strata: str                          # Strata label
    simulation_seed: int                 # Seed used

@dataclass
class VPCBin:
    """Data for a single time bin."""
    bin_id: int
    time_min: float
    time_max: float
    time_midpoint: float
    n_observed: int
    n_simulated: int
    percentiles: list[VPCPercentileData]

@dataclass
class VPCPercentileData:
    """Percentile data within a bin."""
    percentile: float                    # Level (0.05, 0.50, 0.95)
    observed: float                      # Observed percentile
    simulated_median: float              # Median of simulated
    simulated_lower: float               # Lower CI
    simulated_upper: float               # Upper CI

@dataclass
class StratifiedVPCResult:
    """Stratified VPC results."""
    results: list[VPCResult]             # VPC per stratum
    stratify_by: list[str]               # Stratification variables
    strata_names: list[str]              # Stratum labels
```

---

## Accessor Functions

Convenient functions to extract data from VPCResult:

```python
from neopkpd.vpc import (
    get_bin_midpoints,
    get_observed_percentile,
    get_simulated_median,
    get_simulated_ci
)

# Extract arrays for plotting
times = get_bin_midpoints(vpc_result)           # np.ndarray of midpoints
obs_p50 = get_observed_percentile(vpc_result, 0.50)  # Observed median
sim_p50 = get_simulated_median(vpc_result, 0.50)     # Simulated median
sim_lo, sim_hi = get_simulated_ci(vpc_result, 0.50)  # CI bounds
```

---

## Pure Python VPC

For standalone use without Julia:

```python
from neopkpd.vpc import compute_vpc_python

# Compute VPC from pre-simulated data
vpc_result = compute_vpc_python(
    observed_times=obs_times,           # np.ndarray
    observed_values=obs_values,         # np.ndarray
    simulated_data=[                    # List of (times, values) tuples
        (sim_times_1, sim_values_1),
        (sim_times_2, sim_values_2),
        # ... n_simulations tuples
    ],
    config=VPCConfig(n_bins=8)
)
```

---

## Complete Example

```python
import neopkpd
from neopkpd import viz
from neopkpd.vpc import VPCConfig, QuantileBinning
import numpy as np

# ================================================
# Complete VPC Workflow
# ================================================

# 1. Initialize
neopkpd.init_julia()

print("=== VPC Computation Example ===\n")

# 2. Generate observed data (from clinical study)
np.random.seed(42)
n_subjects = 50
sampling_times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]

# Simulate "observed" data
subject_ids = []
times = []
dv = []

true_ka, true_cl, true_v = 1.5, 5.0, 50.0
omega_ka, omega_cl, omega_v = 0.16, 0.09, 0.04
dose = 100.0

for i in range(n_subjects):
    # Individual parameters
    ka_i = true_ka * np.exp(np.random.randn() * np.sqrt(omega_ka))
    cl_i = true_cl * np.exp(np.random.randn() * np.sqrt(omega_cl))
    v_i = true_v * np.exp(np.random.randn() * np.sqrt(omega_v))

    for t in sampling_times:
        # One-compartment oral
        conc = dose * ka_i / (v_i * (ka_i - cl_i/v_i)) * \
               (np.exp(-cl_i/v_i * t) - np.exp(-ka_i * t))

        # Add residual error
        conc_obs = conc * (1 + 0.1 * np.random.randn())
        conc_obs = max(0.01, conc_obs)

        subject_ids.append(f"S{i+1}")
        times.append(t)
        dv.append(conc_obs)

observed_data = {
    "subject_ids": subject_ids,
    "times": times,
    "dv": dv,
    "dvid": ["conc"] * len(dv)
}

print(f"Observed data: {n_subjects} subjects, {len(dv)} observations")

# 3. Define population model (final estimates)
population_spec = {
    "model": "OneCompOral",
    "params": {"Ka": true_ka, "CL": true_cl, "V": true_v},
    "omega": {"Ka": omega_ka, "CL": omega_cl, "V": omega_v},
    "doses": [{"time": 0.0, "amount": dose}],
    "n": n_subjects
}

grid = {"t0": 0.0, "t1": 24.0, "saveat": 0.5}

# 4. Compute standard VPC
print("\nComputing standard VPC...")
config = VPCConfig(
    pi_levels=[0.05, 0.50, 0.95],
    binning=QuantileBinning(7),
    n_simulations=500,
    seed=42
)

vpc_result = neopkpd.compute_vpc(
    observed_data=observed_data,
    population_spec=population_spec,
    grid=grid,
    config=config
)

print(f"  Bins: {len(vpc_result.bins)}")
print(f"  Simulations: {vpc_result.n_simulations}")

# 5. Analyze results
print("\n--- Bin Analysis ---")
for bin in vpc_result.bins:
    p50 = next(p for p in bin.percentiles if p.percentile == 0.50)
    in_ci = p50.simulated_lower <= p50.observed <= p50.simulated_upper
    status = "✓" if in_ci else "✗"
    print(f"  t={bin.time_midpoint:5.1f}: obs={p50.observed:6.2f}, "
          f"CI=[{p50.simulated_lower:6.2f}, {p50.simulated_upper:6.2f}] {status}")

# 6. Compute coverage
def compute_coverage(result, level):
    n_in = sum(
        1 for bin in result.bins
        for p in bin.percentiles
        if p.percentile == level and
           p.simulated_lower <= p.observed <= p.simulated_upper
    )
    return n_in / len(result.bins)

print("\n--- Coverage ---")
for level in [0.05, 0.50, 0.95]:
    cov = compute_coverage(vpc_result, level)
    print(f"  P{int(level*100)}: {cov*100:.1f}%")

# 7. Visualize
viz.set_backend("matplotlib")
fig = viz.plot_vpc(
    vpc_result,
    title="Visual Predictive Check",
    xlabel="Time (hours)",
    ylabel="Concentration (mg/L)"
)
fig.savefig("vpc_result.png", dpi=300, bbox_inches="tight")
print("\nSaved: vpc_result.png")

print("\n✓ VPC computation complete")
```

---

## Next Steps

- [Binning Strategies](binning.md) - Binning methods and selection
- [BLQ Handling](blq.md) - Below quantification data
- [VPC Visualization](../viz/vpc.md) - Plotting VPC results
- [Julia VPC](../../julia/vpc/index.md) - Julia implementation details
