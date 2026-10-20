# Autoinduction

PK model for drugs that induce their own metabolism, leading to time-varying clearance that increases with chronic dosing.

---

## Usage

```julia
using NeoPKPD

# Create autoinduction model
model = autoinduction()
params = CustomODEParams(
    CL0 = 10.0,   # Baseline clearance (L/h)
    V = 50.0,     # Volume (L)
    Emax = 2.0,   # Max enzyme induction (fold)
    EC50 = 5.0,   # EC50 for induction
    kenz = 0.1    # Enzyme turnover rate (1/h)
)

# Daily dosing for 1 week
doses = [DoseEvent(i * 24.0, 200.0) for i in 0:6]
spec = ModelSpec(model, "autoinduction", params, doses)
grid = SimGrid(0.0, 168.0, 0:1:168)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `CL0` | Float64 | Baseline clearance (L/h) |
| `V` | Float64 | Volume of distribution (L) |
| `Emax` | Float64 | Maximum enzyme induction (fold increase) |
| `EC50` | Float64 | Concentration for 50% of max induction |
| `kenz` | Float64 | Enzyme turnover rate constant (1/h) |

### Derived Parameters

- **Enzyme half-life**: $t_{1/2,enz} = \ln(2) / k_{enz}$
- **Maximum clearance**: $CL_{max} = CL_0 \cdot (1 + E_{max})$

---

## Model Equations

Enzyme induction signal:
$$E_{induced} = 1 + \frac{E_{max} \cdot C}{EC_{50} + C}$$

Enzyme dynamics:
$$\frac{dE}{dt} = k_{enz} \cdot (E_{induced} - E)$$

Drug elimination:
$$\frac{dA_c}{dt} = -\frac{CL_0 \cdot E}{V} \cdot A_c$$

Initial: $E(0) = 1.0$

---

## Basic Example

```julia
using NeoPKPD

model = autoinduction()
params = CustomODEParams(
    CL0 = 10.0, V = 50.0, Emax = 2.0, EC50 = 5.0, kenz = 0.1
)

# Daily dosing
doses = [DoseEvent(i * 24.0, 200.0) for i in 0:13]
spec = ModelSpec(model, "auto", params, doses)
grid = SimGrid(0.0, 336.0, 0:1:336)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)

conc = result.observations[:conc]
enzyme = result.states[:E_enzyme]
t = result.t

# Compare day 1 vs day 14
day1_cmax = maximum(conc[1:24])
day14_cmax = maximum(conc[313:336])

println("Day 1 Cmax: $day1_cmax mg/L")
println("Day 14 Cmax: $day14_cmax mg/L")
println("Cmax decrease: $((1 - day14_cmax/day1_cmax) * 100)%")
println("Enzyme level day 14: $(enzyme[313])x baseline")
```

---

## Clinical Example: Carbamazepine

```julia
using NeoPKPD

# Carbamazepine-like autoinduction
kenz = 0.693 / (4 * 24)  # 4-day enzyme half-life

model = autoinduction()
params = CustomODEParams(
    CL0 = 2.0, V = 80.0, Emax = 1.5, EC50 = 4.0, kenz = kenz
)

# BID dosing for 4 weeks
doses = [DoseEvent(i * 12.0, 200.0) for i in 0:55]
spec = ModelSpec(model, "carbamazepine", params, doses)
grid = SimGrid(0.0, 672.0, 0:2:672)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)

conc = result.observations[:conc]
enzyme = result.states[:E_enzyme]
t = result.t

# Weekly comparison
println("Week | Cmax (mg/L) | Enzyme")
println("-" ^ 35)
for week in 1:4
    start_h = (week - 1) * 168
    end_h = week * 168
    start_idx = findfirst(x -> x >= start_h, t)
    end_idx = findfirst(x -> x >= end_h, t)

    week_cmax = maximum(conc[start_idx:end_idx])
    week_enzyme = enzyme[end_idx]

    println("$week | $week_cmax | $(week_enzyme)x")
end
```

---

## Dose Adjustment Strategy

```julia
using NeoPKPD

kenz = 0.693 / (5 * 24)

model = autoinduction()

# Stepped dosing: increase weekly
doses = vcat(
    [DoseEvent(i * 24.0, 200.0) for i in 0:6],      # Week 1
    [DoseEvent((7 + i) * 24.0, 300.0) for i in 0:6],  # Week 2
    [DoseEvent((14 + i) * 24.0, 400.0) for i in 0:6]  # Week 3
)

params = CustomODEParams(
    CL0 = 5.0, V = 50.0, Emax = 1.5, EC50 = 4.0, kenz = kenz
)

spec = ModelSpec(model, "titration", params, doses)
grid = SimGrid(0.0, 504.0, 0:2:504)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)

conc = result.observations[:conc]
t = result.t

println("Dose Titration:")
for (week, dose) in enumerate([200, 300, 400])
    trough_idx = findfirst(x -> x >= week * 168 - 2, t)
    println("Week $week ($dose mg QD): Trough = $(conc[trough_idx]) mg/L")
end
```

---

## Clinical Applications

- **Carbamazepine** - classic autoinducer
- **Phenytoin** - at high doses
- **Rifampicin** - potent CYP3A4 inducer
- **Phenobarbital** - enzyme inducer

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Induction signal | $E_{induced} = 1 + E_{max} \cdot C / (EC_{50} + C)$ |
| Enzyme rate | $dE/dt = k_{enz} \cdot (E_{induced} - E)$ |
| Effective clearance | $CL_{eff} = CL_0 \cdot E$ |
| Max clearance | $CL_{max} = CL_0 \cdot (1 + E_{max})$ |
| Enzyme t1/2 | $\ln(2) / k_{enz}$ |

---

## See Also

- [One-Compartment IV](onecomp-iv-bolus.md) - Linear PK
- [Michaelis-Menten](michaelis-menten.md) - Saturable elimination
