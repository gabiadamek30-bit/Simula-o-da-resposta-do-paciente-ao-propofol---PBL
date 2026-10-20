# Non-Compartmental Analysis (NCA)

Examples demonstrating NCA metrics computation from concentration-time data.

## Computed Metrics

| Metric | Description | Units |
|--------|-------------|-------|
| Cmax | Maximum concentration | mg/L |
| Tmax | Time of maximum concentration | h |
| AUC_0_t | AUC from 0 to last observation | mg·h/L |
| AUC_0_inf | AUC extrapolated to infinity | mg·h/L |
| t_half | Terminal half-life | h |
| CL_F | Apparent clearance | L/h |
| Vz_F | Apparent volume | L |
| MRT | Mean residence time | h |

## Examples

| Example | Description | Directory |
|---------|-------------|-----------|
| Basic NCA | Single subject NCA | [01_basic_nca](01_basic_nca/README.md) |
| Population NCA | Summary statistics | [02_population_nca](02_population_nca/README.md) |

## Usage

```julia
using NeoPKPD

# Compute NCA metrics
result = compute_nca(
    times = [0, 0.5, 1, 2, 4, 8, 12, 24],
    concentrations = [2.0, 1.8, 1.5, 1.2, 0.8, 0.4, 0.2, 0.05],
    dose = 100.0,
    route = :iv
)

println("Cmax: $(result.cmax)")
println("AUC_0_inf: $(result.auc_0_inf)")
println("t_half: $(result.t_half)")
```

```python
from neopkpd.nca import compute_nca

result = compute_nca(
    times=[0, 0.5, 1, 2, 4, 8, 12, 24],
    concentrations=[2.0, 1.8, 1.5, 1.2, 0.8, 0.4, 0.2, 0.05],
    dose=100.0,
    route="iv"
)
```
