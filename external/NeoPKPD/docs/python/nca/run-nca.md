# run_nca Function

The primary function for performing non-compartmental analysis in Python.

---

## Overview

```python
from neopkpd.nca import run_nca

result = run_nca(times, conc, dose)
```

---

## Function Signature

```python
def run_nca(
    times: list[float] | np.ndarray,
    conc: list[float] | np.ndarray,
    dose: float,
    *,
    config: NCAConfig | None = None,
    route: str = "extravascular",
    dosing_type: str = "single",
    tau: float | None = None,
    infusion_time: float | None = None,
) -> NCAResult:
    """
    Perform non-compartmental analysis on concentration-time data.

    Parameters
    ----------
    times : array-like
        Time points (hours)
    conc : array-like
        Concentrations (mass/volume, e.g., mg/L)
    dose : float
        Administered dose (mass, e.g., mg)
    config : NCAConfig, optional
        NCA configuration options
    route : str
        "extravascular", "iv_bolus", or "iv_infusion"
    dosing_type : str
        "single", "multiple", or "steady_state"
    tau : float, optional
        Dosing interval (hours) for multiple dose
    infusion_time : float, optional
        IV infusion duration (hours)

    Returns
    -------
    NCAResult
        Object containing all NCA parameters
    """
```

---

## Basic Usage

### Single Dose Oral

```python
from neopkpd.nca import run_nca

times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
conc = [0.0, 1.8, 2.5, 2.0, 1.2, 0.6, 0.3, 0.075]
dose = 100.0  # mg

result = run_nca(times, conc, dose)

print(f"Cmax: {result.cmax:.2f} mg/L")
print(f"Tmax: {result.tmax:.2f} h")
print(f"AUC0-t: {result.auc_0_t:.2f} mg·h/L")
print(f"AUC0-inf: {result.auc_0_inf:.2f} mg·h/L")
print(f"t½: {result.t_half:.2f} h")
print(f"CL/F: {result.cl_f:.2f} L/h")
print(f"Vz/F: {result.vz_f:.1f} L")
```

### IV Bolus

```python
result = run_nca(times, conc, dose, route="iv_bolus")

print(f"CL: {result.cl:.2f} L/h")  # Absolute clearance
print(f"Vz: {result.vz:.1f} L")    # Absolute volume
print(f"Vss: {result.vss:.1f} L")  # Volume at steady state
```

### IV Infusion

```python
# 1-hour infusion
result = run_nca(
    times, conc, dose,
    route="iv_infusion",
    infusion_time=1.0
)

print(f"MRT (corrected): {result.mrt:.2f} h")
```

---

## Administration Routes

### Route Parameter Options

| Route | Description | MRT Correction |
|-------|-------------|----------------|
| `"extravascular"` | Oral, SC, IM, etc. | None |
| `"iv_bolus"` | IV bolus injection | None |
| `"iv_infusion"` | IV infusion | Subtract Tinf/2 |

```python
# Extravascular (default)
result = run_nca(times, conc, dose, route="extravascular")

# IV bolus
result = run_nca(times, conc, dose, route="iv_bolus")

# IV infusion (requires infusion_time)
result = run_nca(times, conc, dose, route="iv_infusion", infusion_time=2.0)
```

---

## Dosing Types

### Single Dose (Default)

```python
result = run_nca(times, conc, dose, dosing_type="single")
```

### Multiple Dose (Non-Steady-State)

```python
# Analysis during accumulation
result = run_nca(
    times, conc, dose,
    dosing_type="multiple",
    tau=12.0  # Dosing interval
)

print(f"AUC0-tau: {result.auc_0_tau:.2f}")
```

### Steady State

```python
# Steady-state analysis
result = run_nca(
    times, conc, dose,
    dosing_type="steady_state",
    tau=12.0
)

print(f"Cmax,ss: {result.cmax:.2f}")
print(f"Cmin,ss: {result.cmin:.2f}")
print(f"Cavg,ss: {result.cavg:.2f}")
print(f"AUC0-tau: {result.auc_0_tau:.2f}")
print(f"Fluctuation: {result.fluctuation:.1f}%")
print(f"Swing: {result.swing:.1f}%")
print(f"Accumulation Index: {result.accumulation_index:.2f}")
```

---

## Configuration Options

### With NCAConfig

```python
from neopkpd.nca import run_nca, NCAConfig

config = NCAConfig(
    method="lin_log_mixed",        # AUC calculation method
    lambda_z_min_points=3,         # Minimum points for λz
    lambda_z_r2_threshold=0.9,     # R² threshold
    extrapolation_max_pct=20.0,    # Warning threshold
    lloq=0.05,                     # Lower limit of quantification
    blq_handling="missing"         # BLQ handling method
)

result = run_nca(times, conc, dose, config=config)
```

See [NCA Configuration](config.md) for complete options.

---

## NCAResult Attributes

### Primary Exposure Metrics

```python
result.cmax         # Maximum concentration
result.tmax         # Time of Cmax
result.cmin         # Minimum concentration (multiple dose)
result.clast        # Last measurable concentration
result.tlast        # Time of last measurable concentration
result.cavg         # Average concentration (steady state)
```

### AUC Metrics

```python
result.auc_0_t      # AUC from 0 to tlast
result.auc_0_inf    # AUC from 0 to infinity
result.auc_0_tau    # AUC over dosing interval
result.aumc_0_t     # AUMC from 0 to tlast
result.aumc_0_inf   # AUMC from 0 to infinity
result.auc_extra_pct  # % AUC extrapolated
```

### Terminal Phase

```python
result.t_half       # Terminal half-life
result.lambda_z_result.lambda_z       # Elimination rate constant
result.lambda_z_result.r_squared      # R² of terminal regression
result.lambda_z_result.n_points       # Points used
result.lambda_z_result.intercept      # Y-intercept
```

### PK Parameters

```python
result.cl_f         # Apparent clearance (CL/F)
result.cl           # Clearance (IV route)
result.vz_f         # Apparent volume at terminal phase (Vz/F)
result.vz           # Volume at terminal phase (IV)
result.vss_f        # Apparent Vss
result.vss          # Volume at steady state (IV)
result.mrt          # Mean residence time
```

### Multiple Dose Metrics

```python
result.fluctuation         # Peak-trough fluctuation %
result.swing               # Swing %
result.accumulation_index  # Racc
```

### Dose-Normalized Metrics

```python
result.cmax_dn      # Cmax/Dose
result.auc_dn       # AUC/Dose
```

### Quality Indicators

```python
result.quality_flags   # List of quality issues
result.warnings        # Warning messages
```

---

## Example: Complete Analysis

```python
from neopkpd.nca import run_nca, NCAConfig

# PK data from oral administration
times = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0]
conc = [0.0, 2.5, 4.8, 5.2, 4.5, 3.8, 2.6, 1.9, 1.0, 0.55, 0.18, 0.02]
dose = 500.0  # mg

# Configure with FDA/EMA-recommended settings
config = NCAConfig(
    method="lin_log_mixed",      # Linear-up, log-linear down
    lambda_z_min_points=3,       # At least 3 points for λz
    lambda_z_r2_threshold=0.9,   # R² ≥ 0.9
    extrapolation_max_pct=20.0,  # Warn if >20% extrapolated
    lloq=0.01                    # LLOQ = 0.01 mg/L
)

# Run NCA
result = run_nca(times, conc, dose, config=config, route="extravascular")

# Report results
print("=" * 50)
print("Non-Compartmental Analysis Report")
print("=" * 50)

print("\n--- Primary Exposure Metrics ---")
print(f"Cmax:     {result.cmax:.2f} mg/L")
print(f"Tmax:     {result.tmax:.2f} h")
print(f"AUC0-t:   {result.auc_0_t:.2f} mg·h/L")
print(f"AUC0-inf: {result.auc_0_inf:.2f} mg·h/L")
print(f"AUC extrapolated: {result.auc_extra_pct:.1f}%")

print("\n--- Terminal Phase ---")
print(f"t½:       {result.t_half:.2f} h")
print(f"λz:       {result.lambda_z_result.lambda_z:.4f} 1/h")
print(f"λz R²:    {result.lambda_z_result.r_squared:.4f}")
print(f"Points:   {result.lambda_z_result.n_points}")

print("\n--- PK Parameters ---")
print(f"CL/F:     {result.cl_f:.2f} L/h")
print(f"Vz/F:     {result.vz_f:.1f} L")
print(f"Vss/F:    {result.vss_f:.1f} L")
print(f"MRT:      {result.mrt:.2f} h")

print("\n--- Dose-Normalized ---")
print(f"Cmax/D:   {result.cmax_dn:.4f} mg/L/mg")
print(f"AUC/D:    {result.auc_dn:.4f} h/L")

# Quality assessment
print("\n--- Quality Assessment ---")
if result.lambda_z_result.r_squared >= 0.9:
    print("✓ λz R² meets threshold")
else:
    print("⚠ λz R² below threshold")

if result.auc_extra_pct <= 20:
    print("✓ AUC extrapolation acceptable")
else:
    print("⚠ High AUC extrapolation")

if result.warnings:
    print(f"Warnings: {result.warnings}")
```

---

## Working with NumPy Arrays

```python
import numpy as np
from neopkpd.nca import run_nca

# NumPy arrays work directly
times = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0])
conc = np.array([0.0, 1.8, 2.5, 2.0, 1.2, 0.6, 0.3, 0.075])

result = run_nca(times, conc, 100.0)
```

---

## Working with Pandas

```python
import pandas as pd
from neopkpd.nca import run_nca

# From DataFrame
df = pd.DataFrame({
    'time': [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
    'conc': [0.0, 1.8, 2.5, 2.0, 1.2, 0.6, 0.3, 0.075]
})

result = run_nca(
    df['time'].values,
    df['conc'].values,
    dose=100.0
)
```

---

## Error Handling

```python
from neopkpd.nca import run_nca, NCAConfig
import math

# Handle potential issues
result = run_nca(times, conc, dose)

# Check for valid λz estimation
if math.isnan(result.lambda_z_result.lambda_z):
    print("WARNING: λz could not be estimated")
    print("AUC0-inf and t½ are not reliable")

# Check extrapolation
if result.auc_extra_pct > 20:
    print(f"WARNING: {result.auc_extra_pct:.1f}% of AUC extrapolated")

# Access quality flags
for flag in result.quality_flags:
    print(f"Quality issue: {flag}")
```

---

## See Also

- [NCA Configuration](config.md) - Configuration options
- [Population NCA](population.md) - Multi-subject analysis
- [Bioequivalence](bioequivalence.md) - BE analysis

