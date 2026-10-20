# Key Features

NeoPKPD provides a comprehensive suite of pharmacometrics capabilities designed for research, clinical development, and regulatory submissions.

---

## Pharmacokinetic Models

### Compartmental Models

NeoPKPD supports all standard compartmental PK models with rigorous mathematical formulations:

=== "One-Compartment"

    ```julia
    # IV Bolus
    params = OneCompIVBolusParams(CL=5.0, V=50.0)

    # Oral First-Order
    params = OneCompOralFirstOrderParams(Ka=1.5, CL=5.0, V=50.0)
    ```

    $$\frac{dA}{dt} = -\frac{CL}{V} \cdot A$$

=== "Two-Compartment"

    ```julia
    params = TwoCompIVBolusParams(CL=5.0, V1=20.0, Q=15.0, V2=50.0)
    ```

    $$\frac{dA_1}{dt} = -k_{10} A_1 - k_{12} A_1 + k_{21} A_2$$

=== "Three-Compartment"

    ```julia
    params = ThreeCompIVBolusParams(CL=5.0, V1=10.0, Q2=20.0, V2=30.0, Q3=5.0, V3=100.0)
    ```

    Tri-exponential decline with shallow and deep peripheral compartments.

### Advanced Absorption

| Model | Description | Use Case |
|-------|-------------|----------|
| Transit Compartments | Chain of N compartments | Delayed/complex absorption |
| Lag Time | Fixed delay before absorption | Enteric-coated tablets |
| Zero-Order | Constant-rate absorption | Controlled-release |

### Nonlinear Kinetics

- **Michaelis-Menten Elimination**: Saturable metabolism
- **TMDD Models**: Target-mediated drug disposition for biologics

---

## Pharmacodynamic Models

### Direct Response Models

```python
# Direct Emax
result = neopkpd.simulate_pkpd_direct_emax(
    cl=5.0, v=50.0,
    e0=0.0, emax=100.0, ec50=2.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0
)

# Sigmoid Emax (Hill Equation)
result = neopkpd.simulate_pkpd_sigmoid_emax(
    cl=5.0, v=50.0,
    e0=0.0, emax=100.0, ec50=2.0, gamma=2.5,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0
)
```

### Indirect Response Models

All four indirect response types are supported:

| Type | Mechanism | Drug Effect |
|------|-----------|-------------|
| **Type I** | Inhibition of Kin | Decreases response |
| **Type II** | Inhibition of Kout | Increases response |
| **Type III** | Stimulation of Kin | Increases response |
| **Type IV** | Stimulation of Kout | Decreases response |

### Effect Compartment (Biophase)

Models temporal delay between plasma and effect site concentrations:

$$\frac{dC_e}{dt} = k_{e0} \cdot (C_p - C_e)$$

---

## Population Modeling

### Inter-Individual Variability (IIV)

Log-normal distribution of parameters across subjects:

```julia
# Define omega matrix (variance of random effects)
omega = OmegaMatrix([
    0.09 0.0;   # 30% CV on CL
    0.0  0.04   # 20% CV on V
])

pop_spec = PopulationSpec(
    model_spec,
    100,           # n subjects
    omega,
    12345          # seed
)
```

### Inter-Occasion Variability (IOV)

Parameter variation between dosing occasions within the same subject.

### Covariate Models

Multiple covariate functional forms:

| Type | Formula | Example |
|------|---------|---------|
| Linear | $\theta \cdot (1 + \theta_{cov} \cdot (COV - COV_{ref}))$ | Age effect |
| Power | $\theta \cdot (COV/COV_{ref})^{\theta_{cov}}$ | Weight on CL |
| Exponential | $\theta \cdot e^{\theta_{cov} \cdot COV}$ | Enzyme induction |
| Categorical | $\theta \cdot \theta_{cov}$ if category | Sex effect |

### Residual Error Models

- Additive: $Y = F + \epsilon$
- Proportional: $Y = F \cdot (1 + \epsilon)$
- Combined: $Y = F \cdot (1 + \epsilon_1) + \epsilon_2$
- Exponential: $Y = F \cdot e^\epsilon$

---

## Parameter Estimation

### FOCE-I (First-Order Conditional Estimation with Interaction)

The gold standard for NLME estimation:

```julia
result = estimate(
    data,
    model_spec,
    FOCEConfig(
        max_iterations=1000,
        tolerance=1e-6
    )
)
```

### SAEM (Stochastic Approximation EM)

Robust estimation for complex models:

```julia
result = estimate(
    data,
    model_spec,
    SAEMConfig(
        n_iterations=500,
        n_burn_in=100,
        n_chains=3
    )
)
```

### Laplacian Estimation

For sparse data or complex likelihoods.

### Estimation Diagnostics

- Convergence assessment
- Parameter uncertainty (SE, RSE, CI)
- Shrinkage calculation
- Correlation matrices
- Objective function comparison

---

## Non-Compartmental Analysis

### FDA/EMA-Compliant Metrics

| Metric | Description |
|--------|-------------|
| Cmax | Maximum concentration |
| Tmax | Time to maximum concentration |
| AUC0-t | Area under curve to last observation |
| AUC0-inf | Area extrapolated to infinity |
| t½ | Terminal half-life |
| CL/F | Apparent clearance |
| Vz/F | Apparent volume |
| MRT | Mean residence time |

### Bioequivalence Analysis

```python
from neopkpd.nca import bioequivalence_90ci, tost_analysis

# 90% CI for geometric mean ratio
lower, upper = bioequivalence_90ci(test_auc, reference_auc)

# TOST analysis
result = tost_analysis(
    test_auc, reference_auc,
    theta_lower=0.80, theta_upper=1.25
)
```

---

## Visual Predictive Check (VPC)

### Standard VPC

Comparison of observed vs simulated percentiles:

```python
from neopkpd import viz

fig = viz.plot_vpc_detailed(
    simulated_data,
    observed_data,
    prediction_intervals=[0.05, 0.50, 0.95]
)
```

### Prediction-Corrected VPC (pcVPC)

Corrects for design differences across bins.

### Stratified VPC

Separate VPC by covariate strata (e.g., dose group, renal function).

### BLQ Handling

Proper visualization of below-limit-of-quantification observations.

---

## Clinical Trial Simulation

### Study Designs

| Design | Description |
|--------|-------------|
| Parallel | Independent treatment groups |
| Crossover (2×2) | Standard two-period crossover |
| Williams | Three-period, three-treatment |
| 3+3 | Traditional dose escalation |
| mTPI/CRM | Model-based escalation |
| Adaptive | Interim analysis with adaptation |

### Power Analysis

```python
from neopkpd import trial

power = trial.estimate_power_analytical(
    n_per_arm=50,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05
)

sample_size = trial.estimate_sample_size(
    target_power=0.80,
    effect_size=0.5,
    sd=1.0,
    alpha=0.05
)
```

### Virtual Population Generation

Generate realistic virtual subjects with specified demographics.

---

## Visualization

### Dual Backend Support

=== "Matplotlib"

    ```python
    from neopkpd import viz
    viz.set_backend("matplotlib")

    fig = viz.plot_conc_time(result)
    fig.savefig("plot.png", dpi=300)
    ```

=== "Plotly"

    ```python
    from neopkpd import viz
    viz.set_backend("plotly")

    fig = viz.plot_conc_time(result)
    fig.write_html("plot.html")
    ```

### 55+ Visualization Functions

| Category | Functions |
|----------|-----------|
| PK Plots | `plot_conc_time`, `plot_spaghetti`, `plot_mean_ribbon` |
| NCA Plots | `plot_lambda_z_fit`, `plot_auc_visualization` |
| PKPD Plots | `plot_effect_conc`, `plot_hysteresis` |
| VPC Plots | `plot_vpc_detailed`, `plot_pcvpc`, `plot_stratified_vpc` |
| Estimation | `plot_convergence`, `plot_shrinkage`, `plot_eta_distributions` |
| Bootstrap | `plot_bootstrap_distributions`, `plot_bootstrap_ci` |
| Sensitivity | `plot_tornado`, `plot_spider`, `plot_waterfall` |
| Trial | `plot_power_curve`, `plot_kaplan_meier` |

---

## Model Import

### NONMEM

Parse NONMEM control stream files:

```bash
./bin/neopkpd import --input run001.ctl --format nonmem --out model.json
```

Supported ADVAN subroutines: ADVAN1, ADVAN2, ADVAN3, ADVAN4, ADVAN11, ADVAN12

### Monolix

Parse Monolix project files (.mlxtran).

### CDISC Data

Import PC, EX, and DM domains in CSV or XPT format.

---

## Reproducibility

### Artifact System

Every simulation produces a JSON artifact containing:

- Complete model specification
- Solver settings
- Input data hash
- Full results with numerical precision
- Semantic version fingerprint

### Golden Artifact Validation

```bash
./bin/neopkpd validate-golden
```

Ensures bit-exact reproducibility across code changes.

### Replay Capability

```bash
./bin/neopkpd replay --artifact simulation.json
```

Reproduce any previous simulation exactly.

---

## Performance

### Julia Core

- Just-in-time compilation for near-native speed
- Efficient ODE solvers (DifferentialEquations.jl)
- Parallelized population simulations

### Python Integration

- Seamless Julia-Python bridge via juliacall
- NumPy/Pandas integration
- Lazy initialization for fast startup

---

## Next Steps

- [Getting Started Guide](getting-started.md) - Installation and first simulation
- [Architecture Overview](architecture.md) - System design details
- [Julia Tutorial](../julia/tutorial.md) - Complete Julia walkthrough
- [Python Tutorial](../python/tutorial.md) - Complete Python walkthrough
