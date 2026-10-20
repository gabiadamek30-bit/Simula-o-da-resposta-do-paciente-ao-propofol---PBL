# Visualization Examples

Examples demonstrating plotting and visualization for pharmacometric analyses.

## Plot Types

| Plot | Use Case | Directory |
|------|----------|-----------|
| Concentration-Time | Single subject PK | [01_concentration_time](01_concentration_time/README.md) |
| Spaghetti Plot | Population PK | [02_population_spaghetti](02_population_spaghetti/README.md) |
| VPC Plot | Model validation | [03_vpc_plots](03_vpc_plots/README.md) |
| Diagnostics | Estimation GOF | [04_estimation_diagnostics](04_estimation_diagnostics/README.md) |

## Quick Start

```python
from neopkpd.viz import plot_concentration_time, plot_vpc, plot_gof
import matplotlib.pyplot as plt

# Concentration-time plot
fig, ax = plot_concentration_time(result, title="PK Profile")
plt.savefig("pk_profile.png")

# VPC plot
fig, ax = plot_vpc(vpc_result, title="Visual Predictive Check")
plt.savefig("vpc.png")

# Goodness-of-fit
fig = plot_gof(estimation_result)
plt.savefig("gof.png")
```

## Customization

All plotting functions return matplotlib figures that can be customized:

```python
fig, ax = plot_concentration_time(result)
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Concentration (ng/mL)")
ax.set_xlim(0, 24)
ax.grid(True, alpha=0.3)
fig.savefig("custom_plot.png", dpi=300, bbox_inches="tight")
```

## See Also

- [NCA Examples](../nca/README.md) - NCA metrics
- [VPC Examples](../vpc/README.md) - VPC computation
- [Estimation Examples](../estimation/README.md) - Parameter estimation
