# Disease Progression Model

PD model for tumor growth dynamics with drug-induced cell kill, supporting multiple growth models.

---

## Usage

```julia
using NeoPKPD

# PK setup
pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)
doses = [DoseEvent(0.0, 500.0)]
pk_spec = ModelSpec(pk_model, "pk", pk_params, doses)

# Disease Progression PD
pd_model = DiseaseProgressionPD(GompertzGrowth)
pd_params = DiseaseProgressionPDParams(
    100.0,    # s0: initial tumor size
    0.02,     # kgrow: growth rate
    1000.0,   # smax: carrying capacity
    0.0,      # alpha: linear rate (not used)
    0.01      # kdrug: drug kill rate
)

pd_spec = PDSpec(pd_model, "tumor", pd_params, :conc, :tumor_size)

grid = SimGrid(0.0, 336.0, 0:2:336)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `s0` | Float64 | Initial tumor size |
| `kgrow` | Float64 | Growth rate constant |
| `smax` | Float64 | Maximum size (carrying capacity) |
| `alpha` | Float64 | Linear growth rate |
| `kdrug` | Float64 | Drug-induced cell kill rate |

---

## Growth Models

| Model | Constructor | Equation |
|-------|-------------|----------|
| Exponential | `ExponentialGrowth` | $dS/dt = k_{grow} \cdot S$ |
| Linear | `LinearGrowth` | $dS/dt = \alpha$ |
| Logistic | `LogisticGrowth` | $dS/dt = k_{grow} \cdot S \cdot (1 - S/S_{max})$ |
| Gompertz | `GompertzGrowth` | $dS/dt = k_{grow} \cdot S \cdot \ln(S_{max}/S)$ |
| Asymptotic | `AsymptoticGrowth` | $dS/dt = k_{grow} \cdot (S_{max} - S)$ |

Drug effect: $-k_{drug} \cdot C \cdot S$

---

## Basic Example

```julia
using NeoPKPD

# PK setup
pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)
doses = [DoseEvent(i * 168.0, 200.0) for i in 0:3]  # Weekly
pk_spec = ModelSpec(pk_model, "pk", pk_params, doses)

# Gompertz tumor growth
pd_model = DiseaseProgressionPD(GompertzGrowth)
pd_params = DiseaseProgressionPDParams(100.0, 0.02, 1000.0, 0.0, 0.01)
pd_spec = PDSpec(pd_model, "tumor", pd_params, :conc, :tumor_size)

grid = SimGrid(0.0, 672.0, 0:4:672)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)

tumor = result.observations[:tumor_size]
t = result.t

println("Initial: $(tumor[1])")
println("Day 28: $(tumor[end])")
```

---

## Comparing Growth Models

```julia
using NeoPKPD

models = [
    ("Exponential", ExponentialGrowth),
    ("Logistic", LogisticGrowth),
    ("Gompertz", GompertzGrowth),
    ("Asymptotic", AsymptoticGrowth)
]

pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)
pk_spec = ModelSpec(pk_model, "pk", pk_params, DoseEvent[])

println("Model | Day 14 | Day 28")
println("-" ^ 35)

for (name, growth_type) in models
    pd_model = DiseaseProgressionPD(growth_type)
    pd_params = DiseaseProgressionPDParams(100.0, 0.02, 1000.0, 1.0, 0.0)
    pd_spec = PDSpec(pd_model, "tumor", pd_params, :conc, :tumor_size)

    grid = SimGrid(0.0, 672.0, 0:4:672)
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    tumor = result.observations[:tumor_size]

    println("$name | $(tumor[85]) | $(tumor[end])")
end
```

---

## Tumor Growth Inhibition

```julia
using NeoPKPD

# Reference (no treatment)
pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)
pk_spec_ref = ModelSpec(pk_model, "pk", pk_params, DoseEvent[])

pd_model = DiseaseProgressionPD(LogisticGrowth)
pd_params = DiseaseProgressionPDParams(100.0, 0.02, 1000.0, 0.0, 0.0)
pd_spec = PDSpec(pd_model, "tumor", pd_params, :conc, :tumor_size)

grid = SimGrid(0.0, 336.0, 0:4:336)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result_ref = simulate_pkpd_coupled(pk_spec_ref, pd_spec, grid, solver)
ref_final = result_ref.observations[:tumor_size][end]

# With treatment
doses = [DoseEvent(i * 168.0, 200.0) for i in 0:1]
pk_spec_tx = ModelSpec(pk_model, "pk", pk_params, doses)

pd_params_tx = DiseaseProgressionPDParams(100.0, 0.02, 1000.0, 0.0, 0.01)
pd_spec_tx = PDSpec(pd_model, "tumor", pd_params_tx, :conc, :tumor_size)

result_tx = simulate_pkpd_coupled(pk_spec_tx, pd_spec_tx, grid, solver)
tx_final = result_tx.observations[:tumor_size][end]

tgi = (1 - tx_final / ref_final) * 100
println("Control final: $ref_final")
println("Treatment final: $tx_final")
println("TGI: $tgi%")
```

---

## Clinical Applications

- **Oncology dose-response modeling**
- **Tumor growth inhibition studies**
- **Combination therapy evaluation**
- **Survival surrogate endpoints**

---

## Equations Summary

| Model | dS/dt (without drug) | Steady State |
|-------|---------------------|--------------|
| Exponential | $k_{grow} \cdot S$ | ∞ |
| Linear | $\alpha$ | ∞ |
| Logistic | $k_{grow} \cdot S \cdot (1 - S/S_{max})$ | $S_{max}$ |
| Gompertz | $k_{grow} \cdot S \cdot \ln(S_{max}/S)$ | $S_{max}$ |
| Asymptotic | $k_{grow} \cdot (S_{max} - S)$ | $S_{max}$ |

---

## See Also

- [Transit Compartment PD](transit-pd.md) - Delayed effects
- [Indirect Response](indirect-response.md) - Turnover models
