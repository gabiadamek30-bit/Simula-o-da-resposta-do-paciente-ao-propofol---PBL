# Parallel First-Order Absorption

PK model for drugs with multiple absorption sites or mechanisms, each with distinct absorption rate constants.

---

## Usage

```julia
using NeoPKPD

# Create parallel absorption model
model = parallel_first_order_absorption()
params = CustomODEParams(
    Ka1 = 2.0,    # Fast absorption (1/h)
    Ka2 = 0.5,    # Slow absorption (1/h)
    F1 = 0.6,     # Fraction to fast pathway
    CL = 10.0,    # Clearance (L/h)
    V = 50.0      # Volume (L)
)

doses = [DoseEvent(0.0, 500.0)]
spec = ModelSpec(model, "parallel_abs", params, doses)
grid = SimGrid(0.0, 24.0, 0:0.25:24)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `Ka1` | Float64 | First absorption rate constant (1/h) |
| `Ka2` | Float64 | Second absorption rate constant (1/h) |
| `F1` | Float64 | Fraction to first depot (0-1) |
| `CL` | Float64 | Clearance (L/h) |
| `V` | Float64 | Volume of distribution (L) |

---

## Model Equations

Three-compartment system (two depots, one central):

$$\frac{dA_1}{dt} = -K_{a1} \cdot A_1$$

$$\frac{dA_2}{dt} = -K_{a2} \cdot A_2$$

$$\frac{dA_c}{dt} = K_{a1} \cdot A_1 + K_{a2} \cdot A_2 - \frac{CL}{V} \cdot A_c$$

Initial conditions after dose D:
- $A_1(0) = F_1 \cdot D$
- $A_2(0) = (1-F_1) \cdot D$

---

## Basic Example

```julia
using NeoPKPD

model = parallel_first_order_absorption()
params = CustomODEParams(
    Ka1 = 2.0,
    Ka2 = 0.5,
    F1 = 0.6,
    CL = 10.0,
    V = 50.0
)

doses = [DoseEvent(0.0, 500.0)]
spec = ModelSpec(model, "parallel", params, doses)
grid = SimGrid(0.0, 24.0, 0:0.25:24)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)

conc = result.observations[:conc]
t = result.t

# Find Tmax and Cmax
cmax, idx = findmax(conc)
tmax = t[idx]
println("Tmax: $tmax h")
println("Cmax: $cmax mg/L")
```

---

## Effect of Fraction Split

```julia
using NeoPKPD

model = parallel_first_order_absorption()
f1_values = [0.2, 0.4, 0.6, 0.8]

println("F1 (fast) | Cmax (mg/L) | Tmax (h)")
println("-" ^ 40)

for f1 in f1_values
    params = CustomODEParams(
        Ka1 = 3.0, Ka2 = 0.3, F1 = f1, CL = 10.0, V = 50.0
    )

    doses = [DoseEvent(0.0, 500.0)]
    spec = ModelSpec(model, "parallel", params, doses)
    grid = SimGrid(0.0, 24.0, 0:0.1:24)
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    result = simulate(spec, grid, solver)
    conc = result.observations[:conc]
    t = result.t

    cmax, idx = findmax(conc)
    tmax = t[idx]

    println("$f1 | $cmax | $tmax")
end
```

---

## Clinical Applications

- **Extended-release formulations** with immediate release coating
- **Drugs with multiple absorption sites** in GI tract
- **Food effects** changing absorption pathways
- **Modified-release dosage forms**

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| Depot 1 initial | $A_1(0) = F_1 \cdot Dose$ |
| Depot 2 initial | $A_2(0) = (1-F_1) \cdot Dose$ |
| Total absorption | $K_{a1} \cdot A_1 + K_{a2} \cdot A_2$ |
| Concentration | $C = A_c / V$ |

---

## See Also

- [One-Compartment Oral](onecomp-oral.md) - Single absorption
- [Transit Absorption](transit-absorption.md) - Delayed absorption
