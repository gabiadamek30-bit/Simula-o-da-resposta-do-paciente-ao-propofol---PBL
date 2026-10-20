# Non-Compartmental Analysis

The `neopkpd.nca` module provides FDA/EMA-compliant non-compartmental analysis with comprehensive exposure metrics.

---

## Overview

<div class="grid cards" markdown>

-   :material-function:{ .lg .middle } **run_nca**

    ---

    Complete NCA analysis function

    [:octicons-arrow-right-24: run_nca](run-nca.md)

-   :material-cog:{ .lg .middle } **Configuration**

    ---

    NCAConfig options

    [:octicons-arrow-right-24: Config](config.md)

-   :material-account-group:{ .lg .middle } **Population NCA**

    ---

    Multi-subject analysis

    [:octicons-arrow-right-24: Population](population.md)

-   :material-scale-balance:{ .lg .middle } **Bioequivalence**

    ---

    90% CI and TOST

    [:octicons-arrow-right-24: BE Analysis](bioequivalence.md)

</div>

---

## Quick Start

### Basic NCA

```python
from neopkpd.nca import run_nca

# Concentration-time data
times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
conc = [0.0, 1.8, 2.0, 1.5, 1.0, 0.5, 0.25, 0.06]
dose = 100.0

# Run NCA
result = run_nca(times, conc, dose)

# Access results
print(f"Cmax: {result.cmax:.2f} mg/L")
print(f"Tmax: {result.tmax:.2f} h")
print(f"AUC0-t: {result.auc_0_t:.2f} mg·h/L")
print(f"AUC0-inf: {result.auc_0_inf:.2f} mg·h/L")
print(f"t½: {result.t_half:.2f} h")
print(f"CL/F: {result.cl_f:.2f} L/h")
```

### With Configuration

```python
from neopkpd.nca import run_nca, NCAConfig

config = NCAConfig(
    method="log_linear",           # AUC calculation method
    lambda_z_min_points=3,         # Min points for λz
    lambda_z_r2_threshold=0.9,     # R² threshold
    extrapolation_max_pct=20.0,    # Warning threshold
    blq_handling="zero"            # BLQ handling
)

result = run_nca(times, conc, dose, config=config)
print(f"λz R²: {result.lambda_z_r_squared:.4f}")
```

---

## NCA Metrics

| Metric | Attribute | Units | Description |
|--------|-----------|-------|-------------|
| Cmax | `.cmax` | mg/L | Maximum observed concentration |
| Tmax | `.tmax` | h | Time of Cmax |
| AUC0-t | `.auc_0_t` | mg·h/L | AUC to last observation |
| AUC0-inf | `.auc_0_inf` | mg·h/L | AUC extrapolated to infinity |
| t½ | `.t_half` | h | Terminal half-life |
| λz | `.lambda_z` | 1/h | Terminal elimination rate |
| CL/F | `.cl_f` | L/h | Apparent clearance |
| Vz/F | `.vz_f` | L | Apparent volume |
| MRT | `.mrt` | h | Mean residence time |

---

## AUC Methods

| Method | Description | Use Case |
|--------|-------------|----------|
| `"linear"` | Linear trapezoidal | Ascending phase |
| `"log_linear"` | Log-linear trapezoidal | Descending phase |
| `"lin_log_mixed"` | Mixed (recommended) | General purpose |

```python
# Mixed method (default)
config = NCAConfig(method="lin_log_mixed")
```

---

## Individual Functions

```python
from neopkpd import nca

# Peak metrics
cmax = nca.nca_cmax(times, conc)
tmax = nca.nca_tmax(times, conc)

# AUC calculations
auc_0_t = nca.auc_0_t(times, conc, method="log_linear")
auc_0_inf, extra_pct = nca.auc_0_inf(times, conc, lambda_z)

# Terminal phase
lambda_z, t_half, r_squared = nca.estimate_lambda_z(times, conc)
```

---

## Population NCA

```python
from neopkpd.nca import run_population_nca, summarize_population_nca

# Run NCA for all subjects
pop_results = run_population_nca(population_result, dose=100.0)

# Summarize
summary = summarize_population_nca(pop_results)

print(f"Cmax: {summary['cmax']['mean']:.2f} (CV: {summary['cmax']['cv']:.1f}%)")
print(f"AUC: {summary['auc_0_inf']['mean']:.2f} (CV: {summary['auc_0_inf']['cv']:.1f}%)")
```

---

## Bioequivalence

```python
from neopkpd.nca import bioequivalence_90ci, tost_analysis, be_conclusion

# 90% CI for geometric mean ratio
lower, upper = bioequivalence_90ci(test_auc, reference_auc)
print(f"90% CI: ({lower:.3f}, {upper:.3f})")

# TOST analysis
result = tost_analysis(
    test_auc,
    reference_auc,
    theta_lower=0.80,
    theta_upper=1.25
)

# Conclusion
is_be = be_conclusion(lower, upper, theta_lower=0.80, theta_upper=1.25)
print(f"Bioequivalent: {is_be}")
```

---

## NCA Result Object

```python
class NCAResult:
    # Primary metrics
    cmax: float
    tmax: float
    auc_0_t: float
    auc_0_inf: float
    t_half: float

    # Terminal phase
    lambda_z: float
    lambda_z_r_squared: float
    lambda_z_n_points: int

    # PK parameters
    cl_f: float
    vz_f: float
    mrt: float

    # Quality metrics
    auc_extrapolated_pct: float

    # Metadata
    dose: float
    method: str
```

---

## Next Steps

- [run_nca Details](run-nca.md)
- [NCA Configuration](config.md)
- [Bioequivalence Analysis](bioequivalence.md)
- [NCA Visualization](../viz/nca.md)
