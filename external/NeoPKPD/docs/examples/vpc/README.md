# Visual Predictive Check (VPC)

Examples demonstrating VPC generation for model validation.

## What is VPC?

VPC compares simulated data from the model to observed data:
1. Simulate n replicates (e.g., 500) using final model estimates
2. Compute quantiles (5th, 50th, 95th) of simulated data at each time
3. Overlay observed data
4. Check if observed quantiles fall within simulated confidence intervals

## VPC Types

| Type | Description | Use Case |
|------|-------------|----------|
| Standard VPC | Basic VPC | Initial validation |
| pcVPC | Prediction-corrected | Variable dosing |
| Stratified | By covariate | Covariate effects |
| VPC with BLQ | Below LOQ handling | Assay limits |

## Examples

| Example | Description | Directory |
|---------|-------------|-----------|
| Standard VPC | Basic implementation | [01_standard_vpc](01_standard_vpc/README.md) |
| pcVPC | Prediction-corrected | [02_prediction_corrected_vpc](02_prediction_corrected_vpc/README.md) |

## Usage

```julia
using NeoPKPD

# Generate VPC
vpc_result = compute_vpc(
    observed_data,
    pop_spec,
    n_simulations = 500,
    quantiles = [0.05, 0.5, 0.95],
    seed = 12345
)

# Access results
sim_quantiles = vpc_result.simulated_quantiles  # Simulated CI
obs_quantiles = vpc_result.observed_quantiles   # Observed data
```

## Interpretation

**Good VPC:**
- Observed median close to simulated median
- 90% of observed data within 90% simulated PI
- No systematic deviations

**Poor VPC signals:**
- Median bias → Structural model misspecification
- Over-prediction → IIV too large
- Under-prediction → IIV too small
