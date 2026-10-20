# Enterohepatic Recirculation (EHR)

PK model for drugs that undergo biliary excretion and intestinal reabsorption, leading to secondary concentration peaks.

---

## Usage

```julia
using NeoPKPD

# Create EHR model
model = enterohepatic_recirculation()
params = CustomODEParams(
    Ka = 1.0,       # Absorption rate (1/h)
    CL = 10.0,      # Clearance (L/h)
    V = 50.0,       # Volume (L)
    Kbile = 0.5,    # Biliary excretion rate (1/h)
    Kreab = 0.3,    # Reabsorption rate (1/h)
    F_reab = 0.7    # Fraction reabsorbed
)

doses = [DoseEvent(0.0, 500.0)]
spec = ModelSpec(model, "ehr_sim", params, doses)
grid = SimGrid(0.0, 48.0, 0:0.25:48)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `Ka` | Float64 | Absorption rate constant (1/h) |
| `CL` | Float64 | Clearance (L/h) |
| `V` | Float64 | Volume of distribution (L) |
| `Kbile` | Float64 | Biliary excretion rate constant (1/h) |
| `Kreab` | Float64 | Reabsorption rate constant (1/h) |
| `F_reab` | Float64 | Fraction reabsorbed (0-1) |

---

## Model Equations

Three-compartment system:

$$\frac{dA_{gut}}{dt} = -K_a \cdot A_{gut} + F_{reab} \cdot K_{reab} \cdot A_{bile}$$

$$\frac{dA_c}{dt} = K_a \cdot A_{gut} - \frac{CL}{V} \cdot A_c - K_{bile} \cdot A_c$$

$$\frac{dA_{bile}}{dt} = K_{bile} \cdot A_c - K_{reab} \cdot A_{bile}$$

---

## Basic Example

```julia
using NeoPKPD

model = enterohepatic_recirculation()
params = CustomODEParams(
    Ka = 1.5, CL = 8.0, V = 50.0,
    Kbile = 0.4, Kreab = 0.2, F_reab = 0.8
)

doses = [DoseEvent(0.0, 500.0)]
spec = ModelSpec(model, "ehr", params, doses)
grid = SimGrid(0.0, 48.0, 0:0.1:48)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)

conc = result.observations[:conc]
t = result.t

# Find secondary peaks
peaks = Int[]
for i in 3:length(conc)-2
    if conc[i] > conc[i-1] && conc[i] > conc[i+1]
        push!(peaks, i)
    end
end

println("Number of peaks: $(length(peaks))")
for (n, idx) in enumerate(peaks)
    println("Peak $n: t = $(t[idx]) h, C = $(conc[idx]) mg/L")
end
```

---

## Effect of Reabsorption Fraction

```julia
using NeoPKPD

model = enterohepatic_recirculation()
f_reab_values = [0.0, 0.3, 0.6, 0.9]

println("F_reab | AUC (mg*h/L)")
println("-" ^ 30)

for f_reab in f_reab_values
    params = CustomODEParams(
        Ka = 1.5, CL = 10.0, V = 50.0,
        Kbile = 0.5, Kreab = 0.3, F_reab = f_reab
    )

    doses = [DoseEvent(0.0, 500.0)]
    spec = ModelSpec(model, "ehr", params, doses)
    grid = SimGrid(0.0, 72.0, 0:0.1:72)
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    result = simulate(spec, grid, solver)
    conc = result.observations[:conc]

    # Trapezoidal AUC
    auc = sum((conc[i] + conc[i+1]) * 0.05 for i in 1:length(conc)-1)

    println("$f_reab | $auc")
end
```

---

## Clinical Applications

- **Digoxin** - undergoes EHR
- **NSAIDs** - some undergo biliary cycling
- **Hormones** - estrogens, bile acids
- **Antibiotics** - some macrolides

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| GI absorption | $K_a \cdot A_{gut}$ |
| Biliary excretion | $K_{bile} \cdot A_c$ |
| Reabsorption | $F_{reab} \cdot K_{reab} \cdot A_{bile}$ |
| Net fecal elimination | $(1-F_{reab}) \cdot K_{reab} \cdot A_{bile}$ |

---

## See Also

- [One-Compartment Oral](onecomp-oral.md) - Without EHR
- [Transit Absorption](transit-absorption.md) - Delayed absorption
