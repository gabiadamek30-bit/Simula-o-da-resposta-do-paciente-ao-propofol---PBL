# Tolerance Models

PD models for tolerance development through counter-regulation or receptor regulation mechanisms.

---

## Counter-Regulation Tolerance

### Usage

```julia
using NeoPKPD

# PK setup
pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)
doses = [DoseEvent(i * 8.0, 50.0) for i in 0:20]  # TID for 7 days
pk_spec = ModelSpec(pk_model, "pk", pk_params, doses)

# Tolerance PD
pd_model = ToleranceCounterRegulation()
pd_params = ToleranceCounterRegulationParams(
    0.0,      # e0
    100.0,    # emax
    5.0,      # ec50
    1.0,      # gamma
    0.1,      # kin_mod
    0.05,     # kout_mod
    1.0       # alpha_feedback
)

pd_spec = PDSpec(pd_model, "tol", pd_params, :conc, :effect)

grid = SimGrid(0.0, 168.0, 0:0.5:168)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
```

---

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `e0` | Float64 | Baseline effect |
| `emax` | Float64 | Maximum drug effect |
| `ec50` | Float64 | EC50 for drug effect |
| `gamma` | Float64 | Hill coefficient |
| `kin_mod` | Float64 | Moderator production rate |
| `kout_mod` | Float64 | Moderator elimination rate |
| `alpha_feedback` | Float64 | Feedback strength |

---

### Model Equations

Drug effect:
$$E_{drug} = \frac{E_{max} \cdot C^\gamma}{EC_{50}^\gamma + C^\gamma}$$

Moderator dynamics:
$$\frac{dM}{dt} = k_{in,mod} \cdot E_{drug} - k_{out,mod} \cdot M$$

Net effect:
$$E_{net} = E_0 + E_{drug} - \alpha \cdot M$$

---

## Basic Example

```julia
using NeoPKPD

pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)
doses = [DoseEvent(i * 8.0, 50.0) for i in 0:20]
pk_spec = ModelSpec(pk_model, "pk", pk_params, doses)

pd_model = ToleranceCounterRegulation()
pd_params = ToleranceCounterRegulationParams(0.0, 100.0, 5.0, 1.0, 0.1, 0.05, 1.0)
pd_spec = PDSpec(pd_model, "tol", pd_params, :conc, :effect)

grid = SimGrid(0.0, 168.0, 0:0.5:168)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)

effect = result.observations[:effect]

# Compare first vs last dose
first_peak = maximum(effect[1:16])
last_peak = maximum(effect[end-16:end])

println("First dose peak: $first_peak")
println("Last dose peak: $last_peak")
println("Tolerance: $((1 - last_peak/first_peak) * 100)%")
```

---

## Receptor Regulation Model

### Usage

```julia
using NeoPKPD

# PK setup
pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)
doses = [DoseEvent(i * 8.0, 50.0) for i in 0:20]
pk_spec = ModelSpec(pk_model, "pk", pk_params, doses)

# Receptor regulation PD
pd_model = ReceptorRegulation()
pd_params = ReceptorRegulationParams(
    0.0,      # e0
    100.0,    # emax
    5.0,      # ec50
    1.0,      # gamma
    1.0,      # r_baseline
    0.05,     # kreg
    2.0,      # rmax
    0.02,     # kchange
    :down     # direction
)

pd_spec = PDSpec(pd_model, "receptor", pd_params, :conc, :effect)

grid = SimGrid(0.0, 168.0, 0:0.5:168)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
```

---

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `r_baseline` | Float64 | Baseline receptor density |
| `kreg` | Float64 | Receptor recovery rate |
| `rmax` | Float64 | Maximum receptor density |
| `kchange` | Float64 | Rate of receptor change |
| `direction` | Symbol | `:down` or `:up` |

---

### Receptor Equations

Down-regulation:
$$\frac{dR}{dt} = k_{reg} \cdot (R_0 - R) - k_{change} \cdot E_{drug} \cdot R$$

Up-regulation:
$$\frac{dR}{dt} = k_{reg} \cdot (R_0 - R) + k_{change} \cdot E_{drug} \cdot (R_{max} - R)$$

Net effect:
$$E_{net} = E_0 + R \cdot E_{drug}$$

---

## Down-Regulation Example

```julia
using NeoPKPD

pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)
doses = [DoseEvent(i * 8.0, 50.0) for i in 0:20]
pk_spec = ModelSpec(pk_model, "pk", pk_params, doses)

pd_model = ReceptorRegulation()
pd_params = ReceptorRegulationParams(0.0, 100.0, 5.0, 1.0, 1.0, 0.05, 2.0, 0.02, :down)
pd_spec = PDSpec(pd_model, "beta_receptor", pd_params, :conc, :effect)

grid = SimGrid(0.0, 168.0, 0:0.5:168)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)

effect = result.observations[:effect]
receptor = result.states[:R]

println("Initial receptor: $(receptor[1])")
println("Day 7 receptor: $(receptor[end])")
println("First peak effect: $(maximum(effect[1:16]))")
println("Day 7 peak effect: $(maximum(effect[end-16:end]))")
```

---

## Clinical Applications

### Counter-Regulation
- Opioid tolerance
- Benzodiazepine tolerance
- Beta-blocker tolerance

### Receptor Regulation
- Beta-receptor down-regulation
- Opioid receptor adaptation
- Hormone receptor changes

---

## Equations Summary

### Counter-Regulation

| Quantity | Formula |
|----------|---------|
| Drug effect | $E_{max} \cdot C^\gamma / (EC_{50}^\gamma + C^\gamma)$ |
| Moderator | $dM/dt = k_{in} \cdot E_{drug} - k_{out} \cdot M$ |
| Net effect | $E_0 + E_{drug} - \alpha \cdot M$ |

### Receptor Regulation

| Direction | Rate Equation |
|-----------|---------------|
| Down | $k_{reg}(R_0 - R) - k_{change} \cdot E \cdot R$ |
| Up | $k_{reg}(R_0 - R) + k_{change} \cdot E \cdot (R_{max} - R)$ |

---

## See Also

- [Direct Emax](direct-emax.md) - Without tolerance
- [Indirect Response](indirect-response.md) - Turnover models
- [Effect Compartment](effect-compartment.md) - Temporal delay
