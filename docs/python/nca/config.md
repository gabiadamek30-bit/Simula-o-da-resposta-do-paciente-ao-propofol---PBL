# NCA Configuration

Complete documentation for `NCAConfig` options controlling NCA calculations.

---

## Overview

```python
from neopkpd.nca import NCAConfig, run_nca

config = NCAConfig(
    method="lin_log_mixed",
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9
)

result = run_nca(times, conc, dose, config=config)
```

---

## NCAConfig Parameters

### Complete Parameter List

```python
from neopkpd.nca import NCAConfig

config = NCAConfig(
    # AUC Calculation Method
    method="lin_log_mixed",          # "linear", "log_linear", "lin_log_mixed"

    # Lambda_z Estimation
    lambda_z_min_points=3,           # Minimum points for terminal regression
    lambda_z_max_points=None,        # Maximum points (None = all valid)
    lambda_z_r2_threshold=0.9,       # R² quality threshold
    lambda_z_selection="min_points_first",  # Point selection method
    lambda_z_start_time=None,        # Earliest time to include
    lambda_z_start_idx=None,         # First index to include
    lambda_z_end_idx=None,           # Last index to include

    # BLQ Handling
    lloq=None,                       # Lower limit of quantification
    blq_handling="missing",          # "zero", "lloq_half", "missing"

    # Quality Thresholds
    extrapolation_max_pct=20.0,      # Warning if AUC extrap > this %
)
```

---

## AUC Calculation Methods

### Method Options

| Method | Description | Formula |
|--------|-------------|---------|
| `"linear"` | Linear trapezoidal | $(C_1 + C_2) \cdot \Delta t / 2$ |
| `"log_linear"` | Log-linear trapezoidal | $(C_1 - C_2) / \ln(C_1/C_2) \cdot \Delta t$ |
| `"lin_log_mixed"` | Linear up, log down (recommended) | Mixed |

### Linear Trapezoidal

```python
config = NCAConfig(method="linear")
```

Best for ascending concentration phases.

### Log-Linear Trapezoidal

```python
config = NCAConfig(method="log_linear")
```

Best for descending (elimination) phases when concentrations are decreasing.

### Linear-Log Mixed (Recommended)

```python
config = NCAConfig(method="lin_log_mixed")  # Default
```

FDA/EMA recommended method:
- Uses **linear** trapezoidal when concentration is increasing
- Uses **log-linear** trapezoidal when concentration is decreasing

---

## Lambda_z Configuration

### Minimum Points

```python
# Require at least 4 points for λz estimation
config = NCAConfig(lambda_z_min_points=4)
```

FDA/EMA typically require minimum 3 points.

### R² Threshold

```python
# Require R² ≥ 0.95
config = NCAConfig(lambda_z_r2_threshold=0.95)
```

### Point Selection Methods

```python
# Method 1: MinPointsFirst (FDA/EMA default)
# Starts with minimum points from end, adds more if R² improves
config = NCAConfig(lambda_z_selection="min_points_first")

# Method 2: MaxAdjR2
# Tests all combinations, selects best adjusted R²
config = NCAConfig(lambda_z_selection="max_adj_r2")
```

### Manual Point Selection

```python
# Specify exact points to use
config = NCAConfig(
    lambda_z_start_idx=5,   # Start from index 5
    lambda_z_end_idx=10     # End at index 10
)

# Or by time
config = NCAConfig(
    lambda_z_start_time=4.0  # Only use times ≥ 4h
)
```

---

## BLQ Handling

### LLOQ Setting

```python
# Set lower limit of quantification
config = NCAConfig(lloq=0.05)  # 0.05 mg/L
```

### BLQ Handling Methods

| Method | Description | Use Case |
|--------|-------------|----------|
| `"zero"` | Replace BLQ with 0 | Pre-dose samples |
| `"lloq_half"` | Replace BLQ with LLOQ/2 | Mid-profile BLQ |
| `"missing"` | Exclude from calculations | General use |

```python
# Treat BLQ as zero
config = NCAConfig(lloq=0.05, blq_handling="zero")

# Treat BLQ as LLOQ/2
config = NCAConfig(lloq=0.05, blq_handling="lloq_half")

# Exclude BLQ values
config = NCAConfig(lloq=0.05, blq_handling="missing")
```

---

## Quality Thresholds

### AUC Extrapolation Warning

```python
# Warn if >15% of AUC is extrapolated
config = NCAConfig(extrapolation_max_pct=15.0)
```

FDA/EMA typically flag studies where extrapolation exceeds 20%.

---

## Preset Configurations

### FDA-Compliant

```python
def fda_config():
    return NCAConfig(
        method="lin_log_mixed",
        lambda_z_min_points=3,
        lambda_z_r2_threshold=0.9,
        lambda_z_selection="min_points_first",
        extrapolation_max_pct=20.0,
        blq_handling="missing"
    )

config = fda_config()
```

### EMA-Compliant

```python
def ema_config():
    return NCAConfig(
        method="lin_log_mixed",
        lambda_z_min_points=3,
        lambda_z_r2_threshold=0.9,
        lambda_z_selection="min_points_first",
        extrapolation_max_pct=20.0,
        blq_handling="missing"
    )

config = ema_config()
```

### Conservative (Strict QC)

```python
def strict_config():
    return NCAConfig(
        method="lin_log_mixed",
        lambda_z_min_points=4,
        lambda_z_r2_threshold=0.95,
        extrapolation_max_pct=15.0
    )

config = strict_config()
```

---

## Example: Custom Configuration

```python
from neopkpd.nca import run_nca, NCAConfig

# Bioanalytical assay has LLOQ of 0.1 mg/L
# Study has sparse terminal sampling

config = NCAConfig(
    # Use mixed method per FDA guidance
    method="lin_log_mixed",

    # Lambda_z settings
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9,
    lambda_z_selection="min_points_first",

    # BLQ handling
    lloq=0.1,
    blq_handling="missing",  # Exclude BLQ from calculations

    # Quality thresholds
    extrapolation_max_pct=20.0
)

# Run NCA with configuration
result = run_nca(times, conc, dose, config=config)

# Check quality
print(f"λz R²: {result.lambda_z_result.r_squared:.4f}")
print(f"AUC extrapolated: {result.auc_extra_pct:.1f}%")

if result.lambda_z_result.r_squared < config.lambda_z_r2_threshold:
    print("WARNING: λz R² below threshold")

if result.auc_extra_pct > config.extrapolation_max_pct:
    print("WARNING: High AUC extrapolation")
```

---

## Configuration for Study Types

### Single Dose PK Study

```python
config = NCAConfig(
    method="lin_log_mixed",
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9,
    extrapolation_max_pct=20.0
)

result = run_nca(times, conc, dose, config=config, dosing_type="single")
```

### Steady-State Study

```python
config = NCAConfig(
    method="lin_log_mixed",
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9
)

result = run_nca(
    times, conc, dose,
    config=config,
    dosing_type="steady_state",
    tau=12.0
)
```

### Bioequivalence Study

```python
config = NCAConfig(
    method="lin_log_mixed",
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9,
    extrapolation_max_pct=20.0,  # Critical for BE
    lloq=0.05,
    blq_handling="missing"
)
```

---

## Validating Configuration

```python
from neopkpd.nca import NCAConfig

config = NCAConfig(
    lambda_z_min_points=3,
    lambda_z_r2_threshold=0.9
)

# Check configuration is valid
print(f"Method: {config.method}")
print(f"Min λz points: {config.lambda_z_min_points}")
print(f"R² threshold: {config.lambda_z_r2_threshold}")
print(f"LLOQ: {config.lloq}")
print(f"BLQ handling: {config.blq_handling}")
```

---

## See Also

- [run_nca Function](run-nca.md) - Main NCA function
- [Population NCA](population.md) - Multi-subject analysis
- [Bioequivalence](bioequivalence.md) - BE studies

