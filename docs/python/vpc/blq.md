# BLQ Handling

Comprehensive guide to handling Below Limit of Quantification (BLQ) data in VPC.

---

## Overview

BLQ observations occur when drug concentrations fall below the assay's Lower Limit of Quantification (LLOQ). Proper handling is essential for accurate VPC interpretation.

```python
from neopkpd.vpc import BLQMethod, handle_blq
from neopkpd import compute_vpc_with_blq

# Compute VPC with BLQ handling
vpc_result, blq_stats = compute_vpc_with_blq(
    observed_data=observed,
    population_spec=pop_spec,
    grid=grid,
    config={"lloq": 0.1}  # LLOQ = 0.1 mg/L
)
```

---

## BLQ Methods

Seven methods for handling BLQ data (Beal 2001):

```python
from neopkpd.vpc import BLQMethod

class BLQMethod(Enum):
    M1 = "M1"  # Discard all BLQ
    M3 = "M3"  # Treat as censored (keep original)
    M4 = "M4"  # Replace with LLOQ/2
    M5 = "M5"  # 0 before Tmax, LLOQ/2 after
    M6 = "M6"  # LLOQ/2 before Tmax, discard after
    M7 = "M7"  # 0 before Tmax, discard after
```

### Method Comparison

| Method | Before Tmax | After Tmax | Use Case |
|--------|-------------|------------|----------|
| M1 | Discard | Discard | Simple exclusion |
| M3 | Keep original | Keep original | Censored likelihood |
| M4 | LLOQ/2 | LLOQ/2 | Default, general use |
| M5 | 0 | LLOQ/2 | Pre-dose zeros expected |
| M6 | LLOQ/2 | Discard | Late BLQ uninformative |
| M7 | 0 | Discard | Strict late exclusion |

### Method Selection Guide

```python
# General recommendation
blq_method = BLQMethod.M4  # LLOQ/2 substitution

# When pre-dose should be zero
blq_method = BLQMethod.M5  # or M7

# When late-phase BLQ is uninformative
blq_method = BLQMethod.M6  # or M7

# For censored likelihood estimation
blq_method = BLQMethod.M3  # Keep for special handling
```

---

## handle_blq Function

Apply BLQ handling to concentration data:

```python
from neopkpd.vpc import handle_blq, BLQMethod
import numpy as np

def handle_blq(
    values: np.ndarray,
    times: np.ndarray,
    lloq: float,
    method: BLQMethod = BLQMethod.M4,
    tmax: float | None = None
) -> np.ndarray:
    """
    Handle BLQ observations.

    Parameters
    ----------
    values : np.ndarray
        Concentration values
    times : np.ndarray
        Observation times
    lloq : float
        Lower limit of quantification
    method : BLQMethod
        BLQ handling method
    tmax : float, optional
        Time of maximum concentration (auto-detected if None)

    Returns
    -------
    np.ndarray
        Handled values (NaN for discarded, modified for others)
    """
```

### Usage

```python
from neopkpd.vpc import handle_blq, BLQMethod
import numpy as np

# Sample data with BLQ
times = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0])
values = np.array([0.05, 2.5, 5.0, 3.2, 1.5, 0.8, 0.3, 0.08])
lloq = 0.1

# Original BLQ values
blq_mask = values < lloq
print(f"BLQ observations: {times[blq_mask]}")  # [0.0, 24.0]

# Method M4: Replace with LLOQ/2
values_m4 = handle_blq(values, times, lloq, BLQMethod.M4)
print(f"M4 result: {values_m4}")
# [0.05, 2.5, 5.0, 3.2, 1.5, 0.8, 0.3, 0.05]

# Method M5: 0 before Tmax, LLOQ/2 after
values_m5 = handle_blq(values, times, lloq, BLQMethod.M5)
print(f"M5 result: {values_m5}")
# [0.0, 2.5, 5.0, 3.2, 1.5, 0.8, 0.3, 0.05]

# Method M7: 0 before Tmax, discard after
values_m7 = handle_blq(values, times, lloq, BLQMethod.M7)
print(f"M7 result: {values_m7}")
# [0.0, 2.5, 5.0, 3.2, 1.5, 0.8, 0.3, nan]
```

---

## compute_vpc_with_blq

Compute VPC with BLQ statistics:

```python
from neopkpd import compute_vpc_with_blq
from neopkpd.vpc import VPCConfig, BLQMethod

def compute_vpc_with_blq(
    observed_data: dict,
    population_spec: dict,
    grid: dict,
    config: VPCConfig | None = None,
    solver: dict | None = None,
    error_spec: dict | None = None,
    blq_method: BLQMethod = BLQMethod.M4
) -> tuple[VPCResult, list[BLQBinStats]]:
    """
    Compute VPC with BLQ handling.

    Returns
    -------
    vpc_result : VPCResult
        Standard VPC result
    blq_stats : list[BLQBinStats]
        BLQ statistics per bin
    """
```

### Usage

```python
import neopkpd
from neopkpd.vpc import VPCConfig, BLQMethod

neopkpd.init_julia()

# Configuration with LLOQ
config = VPCConfig(
    lloq=0.1,                    # 0.1 mg/L
    pi_levels=[0.05, 0.50, 0.95],
    n_simulations=500
)

# Compute VPC with BLQ
vpc_result, blq_stats = neopkpd.compute_vpc_with_blq(
    observed_data=observed_data,
    population_spec=population_spec,
    grid=grid,
    config=config,
    blq_method=BLQMethod.M4
)

# Access BLQ statistics
for stat in blq_stats:
    print(f"Bin {stat.bin_id}:")
    print(f"  Total: {stat.n_total}")
    print(f"  BLQ: {stat.n_blq} ({stat.pct_blq_observed:.1f}%)")
    print(f"  Sim BLQ: {stat.pct_blq_simulated_median:.1f}% "
          f"[{stat.pct_blq_simulated_lower:.1f}, {stat.pct_blq_simulated_upper:.1f}]")
```

---

## BLQBinStats Structure

```python
from dataclasses import dataclass

@dataclass
class BLQBinStats:
    """BLQ statistics for a time bin."""
    bin_id: int                          # Bin identifier
    n_total: int                         # Total observations in bin
    n_blq: int                           # Number of BLQ observations
    pct_blq_observed: float              # Observed %BLQ
    pct_blq_simulated_median: float      # Simulated median %BLQ
    pct_blq_simulated_lower: float       # Lower CI bound
    pct_blq_simulated_upper: float       # Upper CI bound
```

### Accessing BLQ Data

```python
from neopkpd.vpc import get_blq_observed, get_blq_simulated

# Extract BLQ data for plotting
blq_obs = get_blq_observed(vpc_result)  # np.ndarray of observed %BLQ
blq_sim_med, blq_sim_lo, blq_sim_hi = get_blq_simulated(vpc_result)
```

---

## Visualization

### VPC with BLQ Panel

```python
from neopkpd import viz

# Two-panel plot: VPC + BLQ fraction
fig = viz.plot_vpc_with_blq(
    vpc_result,
    blq_stats=blq_stats,
    title="VPC with BLQ Analysis",
    blq_ylabel="% Below LOQ",
    figsize=(10, 8)
)

fig.savefig("vpc_blq.png", dpi=300)
```

### Custom BLQ Plot

```python
import matplotlib.pyplot as plt
import numpy as np

# Create BLQ comparison plot
fig, ax = plt.subplots(figsize=(10, 4))

times = [stat.bin_id for stat in blq_stats]  # Or use midpoints
obs_blq = [stat.pct_blq_observed for stat in blq_stats]
sim_med = [stat.pct_blq_simulated_median for stat in blq_stats]
sim_lo = [stat.pct_blq_simulated_lower for stat in blq_stats]
sim_hi = [stat.pct_blq_simulated_upper for stat in blq_stats]

# Simulated CI ribbon
ax.fill_between(times, sim_lo, sim_hi, alpha=0.3, color='blue', label='Sim 95% CI')
ax.plot(times, sim_med, 'b--', label='Sim Median')
ax.plot(times, obs_blq, 'ko-', label='Observed')

ax.set_xlabel('Bin')
ax.set_ylabel('% BLQ')
ax.set_title('BLQ Fraction Comparison')
ax.legend()
ax.set_ylim(0, max(max(sim_hi), max(obs_blq)) * 1.1)

plt.tight_layout()
plt.savefig("blq_comparison.png", dpi=300)
```

---

## Interpretation

### Good Model Fit

```
✓ Observed %BLQ within simulated CI
✓ Similar %BLQ patterns over time
✓ No systematic over/under-prediction of BLQ
```

### Potential Issues

| Pattern | Interpretation | Action |
|---------|----------------|--------|
| Obs %BLQ > Sim CI | Model over-predicts concentrations | Check elimination model |
| Obs %BLQ < Sim CI | Model under-predicts concentrations | Check absorption/bioavailability |
| Early BLQ mismatch | Pre-systemic or lag issues | Consider lag time model |
| Late BLQ mismatch | Terminal phase misspecification | Check terminal half-life |

---

## Complete Example

```python
import neopkpd
from neopkpd import viz
from neopkpd.vpc import VPCConfig, BLQMethod, handle_blq
import numpy as np

# ================================================
# VPC with BLQ Handling Example
# ================================================

neopkpd.init_julia()

print("=== VPC with BLQ Analysis ===\n")

# 1. Generate data with BLQ
np.random.seed(42)
n_subjects = 60
lloq = 0.5  # mg/L

# Simulate concentrations
times_per_subj = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 36.0, 48.0]
subject_ids = []
times = []
dv = []

true_ka, true_cl, true_v = 1.2, 8.0, 60.0
dose = 100.0

for i in range(n_subjects):
    ka_i = true_ka * np.exp(0.4 * np.random.randn())
    cl_i = true_cl * np.exp(0.3 * np.random.randn())
    v_i = true_v * np.exp(0.2 * np.random.randn())

    for t in times_per_subj:
        if t == 0:
            conc = 0.0
        else:
            conc = dose * ka_i / (v_i * (ka_i - cl_i/v_i)) * \
                   (np.exp(-cl_i/v_i * t) - np.exp(-ka_i * t))

        conc_obs = conc * (1 + 0.15 * np.random.randn())
        conc_obs = max(0.0, conc_obs)

        subject_ids.append(f"S{i+1}")
        times.append(t)
        dv.append(conc_obs)

# Count BLQ
n_blq = sum(1 for c in dv if c < lloq)
print(f"Total observations: {len(dv)}")
print(f"BLQ observations: {n_blq} ({100*n_blq/len(dv):.1f}%)")
print(f"LLOQ: {lloq} mg/L")

# 2. Prepare observed data
observed_data = {
    "subject_ids": subject_ids,
    "times": times,
    "dv": dv,
    "dvid": ["conc"] * len(dv)
}

# 3. Population model
population_spec = {
    "model": "OneCompOral",
    "params": {"Ka": true_ka, "CL": true_cl, "V": true_v},
    "omega": {"Ka": 0.16, "CL": 0.09, "V": 0.04},
    "doses": [{"time": 0.0, "amount": dose}],
    "n": n_subjects
}

grid = {"t0": 0.0, "t1": 48.0, "saveat": 0.5}

# 4. VPC configuration with LLOQ
config = VPCConfig(
    lloq=lloq,
    pi_levels=[0.05, 0.50, 0.95],
    n_simulations=500,
    seed=42
)

# 5. Compute VPC with BLQ handling
print("\nComputing VPC with BLQ handling (M4)...")
vpc_result, blq_stats = neopkpd.compute_vpc_with_blq(
    observed_data=observed_data,
    population_spec=population_spec,
    grid=grid,
    config=config,
    blq_method=BLQMethod.M4
)

# 6. Report BLQ statistics
print("\n--- BLQ Statistics by Bin ---")
print("Bin | Time Range  | N Tot | N BLQ | Obs%BLQ | Sim%BLQ CI")
print("-" * 70)

for stat in blq_stats:
    bin_vpc = vpc_result.bins[stat.bin_id - 1]
    time_range = f"[{bin_vpc.time_min:5.1f}, {bin_vpc.time_max:5.1f}]"
    ci_str = f"[{stat.pct_blq_simulated_lower:4.1f}, {stat.pct_blq_simulated_upper:4.1f}]"

    # Check if observed within CI
    in_ci = stat.pct_blq_simulated_lower <= stat.pct_blq_observed <= stat.pct_blq_simulated_upper
    status = "✓" if in_ci else "✗"

    print(f"{stat.bin_id:3} | {time_range} | {stat.n_total:5} | {stat.n_blq:5} | "
          f"{stat.pct_blq_observed:6.1f}% | {ci_str} {status}")

# 7. Assess overall BLQ handling
print("\n--- BLQ Assessment ---")
n_bins_ok = sum(
    1 for stat in blq_stats
    if stat.pct_blq_simulated_lower <= stat.pct_blq_observed <= stat.pct_blq_simulated_upper
)
print(f"Bins with BLQ within CI: {n_bins_ok}/{len(blq_stats)}")

# 8. Visualize
viz.set_backend("matplotlib")
fig = viz.plot_vpc_with_blq(
    vpc_result,
    blq_stats=blq_stats,
    title=f"VPC with BLQ (LLOQ={lloq} mg/L, M4)",
    blq_ylabel="% Below LOQ"
)
fig.savefig("vpc_with_blq.png", dpi=300, bbox_inches="tight")
print("\nSaved: vpc_with_blq.png")

print("\n✓ BLQ analysis complete")
```

---

## References

- Beal SL. Ways to fit a PK model with some data below the quantification limit. J Pharmacokinet Pharmacodyn. 2001;28(5):481-504.
- Bergstrand M, Karlsson MO. Handling data below the limit of quantification in mixed effect models. AAPS J. 2009;11(2):371-380.

---

## See Also

- [VPC Index](index.md) - VPC overview
- [Binning Strategies](binning.md) - Binning methods
- [VPC Visualization](../viz/vpc.md) - Plot options
- [Julia VPC](../../julia/vpc/index.md) - Julia implementation
