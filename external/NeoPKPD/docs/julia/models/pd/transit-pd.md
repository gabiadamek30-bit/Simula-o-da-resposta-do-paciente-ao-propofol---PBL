# Transit Compartment PD Model

PD model with a chain of transit compartments to model delayed drug effects and signal transduction cascades.

---

## Usage

```julia
using NeoPKPD

# Create PK model
pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)  # CL, V
doses = [DoseEvent(0.0, 500.0)]
pk_spec = ModelSpec(pk_model, "pk", pk_params, doses)

# Create Transit Compartment PD
pd_model = TransitCompartmentPD()
pd_params = TransitCompartmentPDParams(
    5,        # n_transit
    0.5,      # ktr
    100.0,    # e0
    -50.0,    # emax
    5.0,      # ec50
    1.0       # gamma
)

pd_spec = PDSpec(pd_model, "pd", pd_params, :conc, :effect)

grid = SimGrid(0.0, 168.0, 0:1:168)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `n_transit` | Int | Number of transit compartments (1-20) |
| `ktr` | Float64 | Transit rate constant (1/h) |
| `e0` | Float64 | Baseline effect |
| `emax` | Float64 | Maximum effect change |
| `ec50` | Float64 | EC50 for effect |
| `gamma` | Float64 | Hill coefficient |

### Key Derived Parameter

**Mean Transit Time (MTT)**:
$$MTT = \frac{N + 1}{K_{tr}}$$

---

## Model Equations

Drug-induced signal:
$$Signal(C) = E_0 + \frac{E_{max} \cdot C^\gamma}{EC_{50}^\gamma + C^\gamma}$$

Transit compartment chain:
$$\frac{dA_1}{dt} = K_{tr} \cdot (Signal(C) - A_1)$$

$$\frac{dA_i}{dt} = K_{tr} \cdot (A_{i-1} - A_i) \quad \text{for } i = 2..N$$

Final effect:
$$Effect = A_N$$

---

## Basic Example

```julia
using NeoPKPD

# PK setup
pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(5.0, 50.0)
doses = [DoseEvent(0.0, 500.0)]
pk_spec = ModelSpec(pk_model, "pk", pk_params, doses)

# PD setup: 5 transit compartments
pd_model = TransitCompartmentPD()
pd_params = TransitCompartmentPDParams(5, 0.5, 100.0, -50.0, 5.0, 1.0)
pd_spec = PDSpec(pd_model, "pd", pd_params, :conc, :effect)

grid = SimGrid(0.0, 168.0, 0:1:168)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)

effect = result.observations[:effect]
t = result.t

# MTT = (5+1) / 0.5 = 12 hours
mtt = (5 + 1) / 0.5
println("MTT: $mtt h")

nadir, idx = findmin(effect)
println("Nadir: $nadir at t = $(t[idx]) h")
```

---

## Clinical Example: Myelosuppression

```julia
using NeoPKPD

# Neutropenia model (MTT ~5 days)
mtt_days = 5.0
n_transit = 5
ktr = (n_transit + 1) / (mtt_days * 24)

pk_model = OneCompIVBolus()
pk_params = OneCompIVBolusParams(3.0, 30.0)
doses = [DoseEvent(0.0, 100.0)]
pk_spec = ModelSpec(pk_model, "chemo", pk_params, doses)

pd_model = TransitCompartmentPD()
pd_params = TransitCompartmentPDParams(n_transit, ktr, 5.0, -4.5, 1.0, 1.5)
pd_spec = PDSpec(pd_model, "neutro", pd_params, :conc, :effect)

grid = SimGrid(0.0, 504.0, 0:2:504)
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)

effect = result.observations[:effect]
t = result.t

nadir, idx = findmin(effect)
println("Baseline ANC: 5.0 x10^9/L")
println("Nadir ANC: $nadir x10^9/L at day $(t[idx]/24)")
```

---

## Clinical Applications

- **Chemotherapy-induced neutropenia**
- **Thrombocytopenia**
- **Delayed biomarker responses**
- **Signal transduction cascades**

---

## Equations Summary

| Quantity | Formula |
|----------|---------|
| MTT | $(N+1) / K_{tr}$ |
| Signal | $E_0 + E_{max} \cdot C^\gamma / (EC_{50}^\gamma + C^\gamma)$ |
| Transit rate | $K_{tr} \cdot (A_{i-1} - A_i)$ |
| CV of delay | $1 / \sqrt{N+1}$ |

---

## See Also

- [Indirect Response](indirect-response.md) - Simpler delay
- [Effect Compartment](effect-compartment.md) - Single delay
- [Disease Progression](disease-progression.md) - Tumor models
