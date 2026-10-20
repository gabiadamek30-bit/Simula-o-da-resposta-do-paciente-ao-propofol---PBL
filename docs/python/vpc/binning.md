# Binning Strategies

Comprehensive guide to VPC binning methods for optimal time stratification.

---

## Overview

Binning groups observations by time for percentile calculation. The choice of binning strategy significantly affects VPC interpretability.

```python
from neopkpd.vpc import (
    QuantileBinning,
    EqualWidthBinning,
    KMeansBinning,
    BinDefinition
)

# Default binning
binning = QuantileBinning(10)  # 10 bins with equal obs count
```

---

## Binning Classes

### QuantileBinning

Creates bins with approximately equal numbers of observations:

```python
from neopkpd.vpc import QuantileBinning

# 10 bins with equal observation count
binning = QuantileBinning(n_bins=10)

# Fewer bins for sparse data
binning = QuantileBinning(n_bins=5)
```

**Algorithm:**
1. Sort observations by time
2. Divide into n_bins groups with equal count
3. Bin boundaries at group transitions

**Best for:**
- Sparse or uneven sampling designs
- Rich sampling at specific times
- Most population PK studies

**Properties:**
- Stable percentile estimates (equal n per bin)
- May have uneven time coverage
- Robust to outliers

### EqualWidthBinning

Creates bins with equal time ranges:

```python
from neopkpd.vpc import EqualWidthBinning

# 10 bins with equal time width
binning = EqualWidthBinning(n_bins=10)

# For 24-hour study: each bin covers 2.4 hours
# Bin 1: 0.0 - 2.4 hr
# Bin 2: 2.4 - 4.8 hr
# ...
```

**Algorithm:**
1. Determine time range [t_min, t_max]
2. Divide into n_bins equal intervals
3. Width = (t_max - t_min) / n_bins

**Best for:**
- Dense, uniformly sampled data
- Steady-state trough sampling
- When time coverage matters

**Properties:**
- Even time coverage
- Variable observations per bin
- May have empty bins if sampling sparse

### KMeansBinning

Uses k-means clustering to find natural time groupings:

```python
from neopkpd.vpc import KMeansBinning

# K-means with 8 clusters
binning = KMeansBinning(n_bins=8)

# With custom max iterations
binning = KMeansBinning(n_bins=8, max_iter=200)
```

**Algorithm:**
1. Initialize k cluster centers
2. Assign times to nearest center
3. Update centers as cluster means
4. Repeat until convergence or max_iter

**Best for:**
- Data with natural time clusters
- Sparse PK sampling designs
- When observation times cluster naturally

**Properties:**
- Adaptive to data structure
- May find unintuitive boundaries
- Results can vary with initialization

---

## BinDefinition

Structure representing a computed bin:

```python
from dataclasses import dataclass

@dataclass
class BinDefinition:
    """Definition of a time bin."""
    bin_id: int              # Bin identifier (1-indexed)
    time_min: float          # Lower bound (inclusive)
    time_max: float          # Upper bound (exclusive)
    time_midpoint: float     # Midpoint for plotting
```

### Accessing Bin Definitions

```python
from neopkpd.vpc import compute_bins

# Compute bins from times
times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
binning = QuantileBinning(5)
bins = compute_bins(times, binning)

for bin_def in bins:
    print(f"Bin {bin_def.bin_id}: [{bin_def.time_min:.1f}, {bin_def.time_max:.1f})")
    print(f"  Midpoint: {bin_def.time_midpoint:.1f}")
```

---

## Choosing Binning Strategy

### Decision Guide

| Data Characteristic | Recommended Binning | Reason |
|---------------------|---------------------|--------|
| Sparse sampling | QuantileBinning(5-8) | Ensures sufficient obs per bin |
| Dense sampling | EqualWidthBinning(8-12) | Even time coverage |
| Natural clusters | KMeansBinning | Adaptive to structure |
| Mixed (sparse + dense) | QuantileBinning | More robust |
| Steady-state only | EqualWidthBinning | Regular intervals |

### Sample Size Considerations

| Subjects | Obs/Subject | Recommended n_bins |
|----------|-------------|-------------------|
| < 20 | 5-7 | 3-5 bins |
| 20-50 | 5-7 | 5-7 bins |
| 50-100 | 5-10 | 7-10 bins |
| > 100 | 7-15 | 8-12 bins |

**Rule of thumb:** Minimum 10-20 observations per bin for stable percentiles.

---

## Comparison Example

```python
from neopkpd.vpc import QuantileBinning, EqualWidthBinning, KMeansBinning
import numpy as np

# Sparse PK sampling times (typical clinical design)
times = np.array([
    # Pre-dose and absorption phase (dense)
    0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0,
    # Distribution/elimination phase (sparse)
    6.0, 8.0, 12.0, 24.0
] * 50)  # 50 subjects

print("=== Binning Strategy Comparison ===\n")

# QuantileBinning
print("QuantileBinning(7):")
from neopkpd.vpc import compute_bins
bins_q = compute_bins(times, QuantileBinning(7))
for b in bins_q:
    n_in_bin = sum(1 for t in times if b.time_min <= t < b.time_max)
    print(f"  [{b.time_min:5.1f}, {b.time_max:5.1f}): n={n_in_bin}")

# EqualWidthBinning
print("\nEqualWidthBinning(7):")
bins_e = compute_bins(times, EqualWidthBinning(7))
for b in bins_e:
    n_in_bin = sum(1 for t in times if b.time_min <= t < b.time_max)
    print(f"  [{b.time_min:5.1f}, {b.time_max:5.1f}): n={n_in_bin}")

# KMeansBinning
print("\nKMeansBinning(7):")
bins_k = compute_bins(times, KMeansBinning(7))
for b in bins_k:
    n_in_bin = sum(1 for t in times if b.time_min <= t < b.time_max)
    print(f"  [{b.time_min:5.1f}, {b.time_max:5.1f}): n={n_in_bin}")
```

Expected output patterns:
- **QuantileBinning**: Even n per bin, uneven time widths
- **EqualWidthBinning**: Even time widths, highly variable n
- **KMeansBinning**: Adapts to dense early sampling

---

## Advanced: Custom Binning

For special requirements, you can create custom bin definitions:

```python
from neopkpd.vpc import BinDefinition

# Manual bin boundaries matching sampling design
manual_bins = [
    BinDefinition(1, 0.0, 0.75, 0.375),      # Pre-Cmax
    BinDefinition(2, 0.75, 2.5, 1.625),      # Around Cmax
    BinDefinition(3, 2.5, 5.0, 3.75),        # Early elimination
    BinDefinition(4, 5.0, 10.0, 7.5),        # Mid elimination
    BinDefinition(5, 10.0, 25.0, 17.5),      # Terminal phase
]

# Use in VPC computation
# (Implementation depends on VPC function accepting custom bins)
```

### Bin Boundary Recommendations

| Sampling Design | Suggested Boundaries |
|-----------------|---------------------|
| Rich PK (0-24h) | 0, 0.5, 1, 2, 4, 8, 12, 24 |
| Sparse PK | 0, 2, 8, 24 |
| Steady-state QD | 0, 4, 8, 12, 16, 20, 24 |
| Steady-state BID | 0, 2, 4, 6, 8, 10, 12 |

---

## Binning Effect on VPC

### Too Few Bins

```
Problem: Wide time ranges, mixed populations
Visual: Smooth but may miss local misfit
Risk: Masks time-varying model errors
```

### Too Many Bins

```
Problem: Few observations per bin
Visual: Noisy percentiles, wide CIs
Risk: False positive model misfit
```

### Optimal Binning

```
Goal: Balance precision and resolution
Target: 15-30 observations per bin (minimum 10)
Visual: Clear percentile trends, reasonable CI width
```

---

## Visualization

```python
from neopkpd import viz
from neopkpd.vpc import VPCConfig, QuantileBinning, EqualWidthBinning
import matplotlib.pyplot as plt

# Compare VPCs with different binning
configs = [
    ("Quantile (5 bins)", VPCConfig(binning=QuantileBinning(5))),
    ("Quantile (10 bins)", VPCConfig(binning=QuantileBinning(10))),
    ("Equal Width (10 bins)", VPCConfig(binning=EqualWidthBinning(10))),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, (name, config) in zip(axes, configs):
    vpc_result = neopkpd.compute_vpc(
        observed_data=observed_data,
        population_spec=pop_spec,
        grid=grid,
        config=config
    )
    plt.sca(ax)
    viz.plot_vpc(vpc_result, title=name)

plt.tight_layout()
plt.savefig("binning_comparison.png", dpi=300)
```

---

## Best Practices

1. **Start with QuantileBinning** - Most robust choice
2. **Match bins to sampling** - Consider nominal times
3. **Check observations per bin** - Minimum 10-20
4. **Avoid empty bins** - Reduce n_bins if needed
5. **Use VPC CI plot** - Visualize bin boundaries

```python
# Recommended workflow
from neopkpd.vpc import VPCConfig, QuantileBinning

# 1. Estimate optimal bins
n_obs = len(observed_data["times"])
n_bins = max(5, min(12, n_obs // 50))  # ~50 obs per bin

# 2. Configure VPC
config = VPCConfig(
    binning=QuantileBinning(n_bins),
    n_simulations=500
)

# 3. Compute and visualize
vpc_result = neopkpd.compute_vpc(...)
fig = viz.plot_vpc_ci(vpc_result, show_bin_boundaries=True)
```

---

## See Also

- [VPC Index](index.md) - VPC overview
- [BLQ Handling](blq.md) - Below quantification data
- [VPC Visualization](../viz/vpc.md) - Plotting options
