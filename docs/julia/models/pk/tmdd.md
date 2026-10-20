# Target-Mediated Drug Disposition (TMDD)

Advanced PK model for drugs that bind to their pharmacological target, forming drug-target complexes that affect both PK and PD behavior.

---

## Usage

```julia
using NeoPKPD

# Create TMDD model specification
model = target_mediated_drug_disposition()
params = CustomODEParams(
    kel = 0.1,    # Drug elimination rate (1/h)
    kon = 0.01,   # Association rate
    koff = 0.001, # Dissociation rate
    ksyn = 1.0,   # Receptor synthesis
    kdeg = 0.1,   # Receptor degradation
    kint = 0.05,  # Complex internalization
    V = 50.0      # Volume (L)
)

doses = [DoseEvent(0.0, 500.0)]
spec = ModelSpec(model, "tmdd_sim", params, doses)
grid = SimGrid(0.0, 72.0, 0:1:72)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `kel` | Float64 | Drug elimination rate constant (1/h) |
| `kon` | Float64 | Association rate constant (1/(conc*h)) |
| `koff` | Float64 | Dissociation rate constant (1/h) |
| `ksyn` | Float64 | Receptor synthesis rate (conc/h) |
| `kdeg` | Float64 | Receptor degradation rate constant (1/h) |
| `kint` | Float64 | Complex internalization rate constant (1/h) |
| `V` | Float64 | Volume of distribution (L) |

### Derived Parameters

- **KD (dissociation constant)**: $K_D = k_{off} / k_{on}$
- **Receptor baseline**: $R_0 = k_{syn} / k_{deg}$

---

## Model Equations

Three-state ODE system:

$$\frac{dL}{dt} = -k_{el} \cdot L - k_{on} \cdot L \cdot R + k_{off} \cdot RL$$

$$\frac{dR}{dt} = k_{syn} - k_{deg} \cdot R - k_{on} \cdot L \cdot R + k_{off} \cdot RL$$

$$\frac{dRL}{dt} = k_{on} \cdot L \cdot R - k_{off} \cdot RL - k_{int} \cdot RL$$

Where:
- L = Free drug (ligand) concentration
- R = Free receptor concentration
- RL = Drug-receptor complex concentration

---

## Basic Example

```julia
using NeoPKPD

model = target_mediated_drug_disposition()
params = CustomODEParams(
    kel = 0.1,
    kon = 0.01,
    koff = 0.001,
    ksyn = 1.0,
    kdeg = 0.1,
    kint = 0.05,
    V = 50.0
)

doses = [DoseEvent(0.0, 500.0)]
spec = ModelSpec(model, "tmdd", params, doses)
grid = SimGrid(0.0, 72.0, 0:0.5:72)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate(spec, grid, solver)

# Access states
conc = result.observations[:conc]
println("Initial free drug: $(conc[1]) mg/L")
println("Free drug at 24h: $(conc[49]) mg/L")
```

---

## Non-Linear PK Behavior

```julia
using NeoPKPD

model = target_mediated_drug_disposition()

doses_list = [50.0, 100.0, 200.0, 500.0, 1000.0]

println("Dose (mg) | Cmax (mg/L) | Apparent t1/2")
println("-" ^ 45)

for dose in doses_list
    params = CustomODEParams(
        kel = 0.1, kon = 0.01, koff = 0.001,
        ksyn = 1.0, kdeg = 0.1, kint = 0.05, V = 50.0
    )

    doses = [DoseEvent(0.0, dose)]
    spec = ModelSpec(model, "tmdd", params, doses)
    grid = SimGrid(0.0, 96.0, 0:0.5:96)
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    result = simulate(spec, grid, solver)
    conc = result.observations[:conc]

    cmax = maximum(conc)
    println("$dose | $cmax | ...")
end
```

---

## Clinical Applications

TMDD is relevant for:

- **Monoclonal antibodies** binding to soluble targets
- **Therapeutic proteins** with receptor-mediated clearance
- **Small molecules** with high-affinity target binding
- **Biologics** with target-mediated disposition

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| KD | $k_{off} / k_{on}$ |
| Receptor baseline | $R_0 = k_{syn} / k_{deg}$ |
| Free drug rate | $-k_{el}L - k_{on}LR + k_{off}RL$ |
| Complex rate | $k_{on}LR - k_{off}RL - k_{int}RL$ |
| Total drug | $L + RL$ |

---

## See Also

- [Michaelis-Menten](michaelis-menten.md) - Saturable elimination
- [Two-Compartment IV](twocomp-iv.md) - Distribution kinetics
