# Specifications are pure data. No solver logic and no hidden defaults.

using LinearAlgebra

# -------------------------
# Shared validation helpers
# -------------------------

function _require_positive(name::String, x::Float64)
    if !(x > 0.0)
        error("Expected positive value for $(name), got $(x)")
    end
    return nothing
end

# -------------------------
# PK specifications
# -------------------------
export ModelKind,
    OneCompIVBolus, OneCompOralFirstOrder, OneCompIVBolusParams, OneCompOralFirstOrderParams
export TwoCompIVBolus, TwoCompIVBolusParams, TwoCompOral, TwoCompOralParams
export ThreeCompIVBolus, ThreeCompIVBolusParams
export TransitAbsorption, TransitAbsorptionParams
export MichaelisMentenElimination, MichaelisMentenEliminationParams
export DoseEvent, ModelSpec, ModelSpecWithModifiers
export SolverSpec, SimGrid, SimResult
export AbsorptionModifiers, validate_absorption_modifiers

abstract type ModelKind end

struct OneCompIVBolus <: ModelKind end
struct OneCompOralFirstOrder <: ModelKind end

"""
Dose event specification supporting both bolus and infusion administration.

For IV bolus: duration = 0.0 (default), dose is instantaneously added to compartment
For IV infusion: duration > 0.0, dose is delivered at constant rate over duration

Fields:
- time: Start time of dose administration
- amount: Total drug amount to be administered
- duration: Infusion duration in time units (0.0 = instantaneous bolus)

The infusion rate is computed as: rate = amount / duration

Examples:
- DoseEvent(0.0, 100.0) - 100 mg bolus at t=0
- DoseEvent(0.0, 100.0, 0.0) - same as above, explicit bolus
- DoseEvent(0.0, 100.0, 1.0) - 100 mg infused over 1 hour (rate = 100 mg/h)
- DoseEvent(0.0, 100.0, 0.5) - 100 mg infused over 30 min (rate = 200 mg/h)
"""
struct DoseEvent
    time::Float64
    amount::Float64
    duration::Float64

    function DoseEvent(time::Float64, amount::Float64, duration::Float64=0.0)
        duration >= 0.0 || error("DoseEvent duration must be >= 0, got $(duration)")
        amount >= 0.0 || error("DoseEvent amount must be >= 0, got $(amount)")
        new(time, amount, duration)
    end
end

# Note: Two-argument constructor DoseEvent(time, amount) works via the default
# duration=0.0 in the inner constructor above.

"""
Check if dose event is a bolus (instantaneous) administration.
"""
is_bolus(dose::DoseEvent) = dose.duration == 0.0

"""
Check if dose event is an infusion (zero-order input).
"""
is_infusion(dose::DoseEvent) = dose.duration > 0.0

"""
Compute the infusion rate for a dose event.
Returns the rate in amount/time units. For bolus doses, returns Inf.
"""
function infusion_rate(dose::DoseEvent)
    if is_bolus(dose)
        return Inf
    else
        return dose.amount / dose.duration
    end
end

"""
Get the end time of a dose event (when infusion stops).
For bolus doses, this equals the start time.
"""
dose_end_time(dose::DoseEvent) = dose.time + dose.duration

export is_bolus, is_infusion, infusion_rate, dose_end_time

struct OneCompIVBolusParams
    CL::Float64
    V::Float64
end

struct OneCompOralFirstOrderParams
    Ka::Float64
    CL::Float64
    V::Float64
end

# -------------------------
# Two-compartment PK models
# -------------------------

"""
Two-compartment IV bolus PK model.

States:
- A_central: Amount in central compartment
- A_peripheral: Amount in peripheral compartment

Parameters:
- CL: Clearance from central compartment (volume/time)
- V1: Volume of central compartment
- Q: Inter-compartmental clearance (volume/time)
- V2: Volume of peripheral compartment

Micro-constants:
- k10 = CL/V1 (elimination rate constant)
- k12 = Q/V1 (central to peripheral rate constant)
- k21 = Q/V2 (peripheral to central rate constant)

Dynamics:
dA_central/dt = -k10*A_central - k12*A_central + k21*A_peripheral
dA_peripheral/dt = k12*A_central - k21*A_peripheral
"""
struct TwoCompIVBolus <: ModelKind end

struct TwoCompIVBolusParams
    CL::Float64   # Clearance
    V1::Float64   # Central volume
    Q::Float64    # Inter-compartmental clearance
    V2::Float64   # Peripheral volume
end

"""
Two-compartment oral first-order absorption PK model.

States:
- A_gut: Amount in gut compartment
- A_central: Amount in central compartment
- A_peripheral: Amount in peripheral compartment

Parameters:
- Ka: Absorption rate constant (1/time)
- CL: Clearance from central compartment (volume/time)
- V1: Volume of central compartment
- Q: Inter-compartmental clearance (volume/time)
- V2: Volume of peripheral compartment

Dynamics:
dA_gut/dt = -Ka*A_gut
dA_central/dt = Ka*A_gut - (CL/V1)*A_central - (Q/V1)*A_central + (Q/V2)*A_peripheral
dA_peripheral/dt = (Q/V1)*A_central - (Q/V2)*A_peripheral
"""
struct TwoCompOral <: ModelKind end

struct TwoCompOralParams
    Ka::Float64   # Absorption rate constant
    CL::Float64   # Clearance
    V1::Float64   # Central volume
    Q::Float64    # Inter-compartmental clearance
    V2::Float64   # Peripheral volume
end

# -------------------------
# Three-compartment PK model
# -------------------------

"""
Three-compartment IV bolus PK model (mammillary).

States:
- A_central: Amount in central compartment
- A_periph1: Amount in first peripheral (shallow) compartment
- A_periph2: Amount in second peripheral (deep) compartment

Parameters:
- CL: Clearance from central compartment (volume/time)
- V1: Volume of central compartment
- Q2: Inter-compartmental clearance to shallow peripheral (volume/time)
- V2: Volume of shallow peripheral compartment
- Q3: Inter-compartmental clearance to deep peripheral (volume/time)
- V3: Volume of deep peripheral compartment

Dynamics:
dA_central/dt = -(CL/V1)*A_central - (Q2/V1)*A_central + (Q2/V2)*A_periph1
                - (Q3/V1)*A_central + (Q3/V3)*A_periph2
dA_periph1/dt = (Q2/V1)*A_central - (Q2/V2)*A_periph1
dA_periph2/dt = (Q3/V1)*A_central - (Q3/V3)*A_periph2
"""
struct ThreeCompIVBolus <: ModelKind end

struct ThreeCompIVBolusParams
    CL::Float64   # Clearance
    V1::Float64   # Central volume
    Q2::Float64   # Inter-compartmental clearance (shallow)
    V2::Float64   # Shallow peripheral volume
    Q3::Float64   # Inter-compartmental clearance (deep)
    V3::Float64   # Deep peripheral volume
end

# -------------------------
# Transit absorption model
# -------------------------

"""
Transit compartment absorption model.

This model implements a chain of transit compartments before the absorption
compartment, providing a delayed and more physiological absorption profile.
Based on Savic et al. (2007) transit compartment model.

States:
- Transit[1:N]: Amount in each transit compartment
- A_central: Amount in central compartment

Parameters:
- N: Number of transit compartments (integer >= 1)
- Ktr: Transit rate constant (1/time) - same for all transit compartments
- Ka: Absorption rate constant from last transit to central (1/time)
- CL: Clearance (volume/time)
- V: Volume of distribution

Dynamics:
For i = 1: dTransit[1]/dt = -Ktr * Transit[1]  (receives dose)
For i > 1: dTransit[i]/dt = Ktr * Transit[i-1] - Ktr * Transit[i]
dA_central/dt = Ka * Transit[N] - (CL/V) * A_central

Note: Mean transit time (MTT) ≈ (N+1) / Ktr
"""
struct TransitAbsorption <: ModelKind end

struct TransitAbsorptionParams
    N::Int        # Number of transit compartments
    Ktr::Float64  # Transit rate constant
    Ka::Float64   # Absorption rate constant
    CL::Float64   # Clearance
    V::Float64    # Volume of distribution
end

# -------------------------
# Michaelis-Menten elimination model
# -------------------------

"""
One-compartment PK model with Michaelis-Menten (saturable) elimination.

This model describes nonlinear pharmacokinetics where the elimination
pathway becomes saturated at higher concentrations.

States:
- A_central: Amount in central compartment

Parameters:
- Vmax: Maximum elimination rate (mass/time)
- Km: Michaelis constant - concentration at half Vmax (mass/volume)
- V: Volume of distribution

Dynamics:
C = A_central / V
dA_central/dt = -Vmax * C / (Km + C)
             = -Vmax * A_central / (Km * V + A_central)

At low concentrations (C << Km): Approximates first-order with CL ≈ Vmax/Km
At high concentrations (C >> Km): Approximates zero-order with rate ≈ Vmax
"""
struct MichaelisMentenElimination <: ModelKind end

struct MichaelisMentenEliminationParams
    Vmax::Float64  # Maximum elimination rate
    Km::Float64    # Michaelis constant
    V::Float64     # Volume of distribution
end

struct ModelSpec{K<:ModelKind,P}
    kind::K
    name::String
    params::P
    doses::Vector{DoseEvent}
end

# -------------------------
# Absorption Modifiers (ALAG and F)
# -------------------------

"""
    AbsorptionModifiers

Industry-standard absorption modifiers following NONMEM conventions.

# Fields
- `alag::Float64`: Absorption lag time (ALAG) - delays dose onset (time units, >= 0)
- `bioavailability::Float64`: Bioavailability fraction (F) - fraction reaching systemic circulation

# NONMEM Equivalence
In NONMEM \$PK block:
```
ALAG1 = THETA(5)
F1 = THETA(6)
```

In NeoPKPD:
```julia
AbsorptionModifiers(alag=THETA5, bioavailability=THETA6)
```

# Clinical Applications

## ALAG (Absorption Lag Time)
- Gastric emptying delay for oral drugs
- Delayed-release formulations
- Enteric-coated tablets
- Typical range: 0.1 - 2.0 hours for oral drugs

## Bioavailability (F)
- Incomplete absorption (F < 1)
- First-pass metabolism (F < 1)
- Enhanced absorption with food (F may vary)
- Typical range: 0.1 - 1.0 for most oral drugs

# Example
```julia
# Oral drug with 30 min lag and 70% bioavailability
modifiers = AbsorptionModifiers(alag=0.5, bioavailability=0.7)

# Modified release formulation
modifiers = AbsorptionModifiers(alag=2.0, bioavailability=0.85)

# IV administration (default - no modifications)
modifiers = AbsorptionModifiers()
```
"""
struct AbsorptionModifiers
    alag::Float64           # Absorption lag time (hours or time units)
    bioavailability::Float64  # Bioavailability fraction (0-1 typical, can exceed 1)

    function AbsorptionModifiers(; alag::Float64=0.0, bioavailability::Float64=1.0)
        new(alag, bioavailability)
    end
end

# Convenience constructor with positional arguments
AbsorptionModifiers(alag::Float64, bioavailability::Float64) = AbsorptionModifiers(alag=alag, bioavailability=bioavailability)

# Default modifiers (no lag, 100% bioavailability)
const DEFAULT_ABSORPTION_MODIFIERS = AbsorptionModifiers()

"""
    validate_absorption_modifiers(modifiers::AbsorptionModifiers)

Validate absorption modifiers for physiological plausibility.

# Validation Rules
- ALAG must be non-negative (cannot have negative delay)
- Bioavailability must be positive
- Warning issued if F > 1.0 (unusual but valid for supersaturation/nonlinear absorption)
- Warning issued if ALAG > 24 hours (unusually long for most drugs)
"""
function validate_absorption_modifiers(modifiers::AbsorptionModifiers)
    if modifiers.alag < 0.0
        error("ALAG (absorption lag time) must be non-negative, got $(modifiers.alag)")
    end

    if modifiers.bioavailability <= 0.0
        error("Bioavailability (F) must be positive, got $(modifiers.bioavailability)")
    end

    # Warnings for unusual values
    if modifiers.bioavailability > 1.0
        @warn "Bioavailability (F) = $(modifiers.bioavailability) exceeds 100%. " *
              "This is unusual - verify this is intended (e.g., supersaturation, reference formulation)."
    end

    if modifiers.alag > 24.0
        @warn "ALAG = $(modifiers.alag) hours is unusually long. " *
              "Verify time units are correct."
    end

    return nothing
end

"""
    ModelSpecWithModifiers{K<:ModelKind,P}

Extended model specification that includes absorption modifiers (ALAG and F).

This is the recommended way to specify models with oral absorption
or any route requiring lag time or incomplete bioavailability.

# Fields
- `kind::K`: Model type (e.g., OneCompOralFirstOrder)
- `name::String`: Model name for identification
- `params::P`: Model parameters struct
- `doses::Vector{DoseEvent}`: Dosing schedule
- `absorption_modifiers::AbsorptionModifiers`: ALAG and F values

# Example
```julia
# One-compartment oral model with lag and bioavailability
spec = ModelSpecWithModifiers(
    OneCompOralFirstOrder(),
    "Theophylline PK",
    OneCompOralFirstOrderParams(Ka=1.5, CL=2.0, V=30.0),
    [DoseEvent(0.0, 300.0)],
    AbsorptionModifiers(alag=0.5, bioavailability=0.8)
)
```

# IIV on ALAG and F
For population simulation with IIV on absorption modifiers, use PopulationSpec
with the extended IIV specification that includes :ALAG and :F parameters.

# Notes
- ModelSpec (without modifiers) can still be used and is equivalent to
  ModelSpecWithModifiers with default modifiers (ALAG=0, F=1)
- The simulation engine automatically handles both spec types
"""
struct ModelSpecWithModifiers{K<:ModelKind,P}
    kind::K
    name::String
    params::P
    doses::Vector{DoseEvent}
    absorption_modifiers::AbsorptionModifiers

    function ModelSpecWithModifiers(
        kind::K,
        name::String,
        params::P,
        doses::Vector{DoseEvent},
        absorption_modifiers::AbsorptionModifiers
    ) where {K<:ModelKind,P}
        validate_absorption_modifiers(absorption_modifiers)
        new{K,P}(kind, name, params, doses, absorption_modifiers)
    end
end

# Convenience constructor with default modifiers
function ModelSpecWithModifiers(
    kind::K,
    name::String,
    params::P,
    doses::Vector{DoseEvent}
) where {K<:ModelKind,P}
    return ModelSpecWithModifiers(kind, name, params, doses, DEFAULT_ABSORPTION_MODIFIERS)
end

"""
Convert ModelSpec to ModelSpecWithModifiers with default absorption modifiers.
"""
function with_modifiers(spec::ModelSpec{K,P}) where {K,P}
    return ModelSpecWithModifiers(
        spec.kind,
        spec.name,
        spec.params,
        spec.doses,
        DEFAULT_ABSORPTION_MODIFIERS
    )
end

"""
Convert ModelSpec to ModelSpecWithModifiers with specified absorption modifiers.
"""
function with_modifiers(spec::ModelSpec{K,P}, modifiers::AbsorptionModifiers) where {K,P}
    return ModelSpecWithModifiers(
        spec.kind,
        spec.name,
        spec.params,
        spec.doses,
        modifiers
    )
end

export with_modifiers, DEFAULT_ABSORPTION_MODIFIERS

"""
Check if a model spec has non-default absorption modifiers.
"""
has_absorption_modifiers(spec::ModelSpec) = false
has_absorption_modifiers(spec::ModelSpecWithModifiers) =
    spec.absorption_modifiers.alag != 0.0 || spec.absorption_modifiers.bioavailability != 1.0

"""
Get absorption modifiers from a spec (returns defaults for ModelSpec).
"""
get_absorption_modifiers(spec::ModelSpec) = DEFAULT_ABSORPTION_MODIFIERS
get_absorption_modifiers(spec::ModelSpecWithModifiers) = spec.absorption_modifiers

export has_absorption_modifiers, get_absorption_modifiers

struct SolverSpec
    alg::Symbol
    reltol::Float64
    abstol::Float64
    maxiters::Int
end

struct SimGrid
    t0::Float64
    t1::Float64
    saveat::Vector{Float64}
end

struct SimResult
    t::Vector{Float64}
    states::Dict{Symbol,Vector{Float64}}
    observations::Dict{Symbol,Vector{Float64}}
    metadata::Dict{String,Any}
end

# -------------------------
# PD specifications
# -------------------------

export PDModelKind, DirectEmax, DirectEmaxParams, PDSpec
export SigmoidEmax, SigmoidEmaxParams
export BiophaseEquilibration, BiophaseEquilibrationParams

abstract type PDModelKind end

"""
Direct Emax PD model.

Effect(C) = E0 + (Emax * C) / (EC50 + C)
"""
struct DirectEmax <: PDModelKind end

struct DirectEmaxParams
    E0::Float64
    Emax::Float64
    EC50::Float64
end

"""
Sigmoid Emax PD model (Hill equation).

This model extends the direct Emax model with a Hill coefficient (gamma)
that controls the steepness of the concentration-effect relationship.

Effect(C) = E0 + (Emax * C^gamma) / (EC50^gamma + C^gamma)

Parameters:
- E0: Baseline effect (no drug)
- Emax: Maximum effect above baseline
- EC50: Concentration at 50% of maximum effect
- gamma: Hill coefficient (steepness parameter)
  - gamma = 1: Standard Emax model (hyperbolic)
  - gamma > 1: Steeper (more switch-like) response
  - gamma < 1: More gradual response

Note: gamma is also known as the Hill coefficient or slope factor.
Typical range is 0.5 to 5.
"""
struct SigmoidEmax <: PDModelKind end

struct SigmoidEmaxParams
    E0::Float64      # Baseline effect
    Emax::Float64    # Maximum effect
    EC50::Float64    # Concentration at 50% Emax
    gamma::Float64   # Hill coefficient
end

"""
Biophase equilibration (effect compartment) PD model.

This model introduces a hypothetical effect compartment to account for
temporal delays between plasma concentration and observed effect.
The effect compartment has no volume (doesn't affect PK) and equilibrates
with the plasma concentration via first-order kinetics.

States:
- Ce: Effect site concentration (hypothetical)

Parameters:
- ke0: Effect site equilibration rate constant (1/time)
- E0: Baseline effect
- Emax: Maximum effect
- EC50: Effect site concentration at 50% Emax

Dynamics:
dCe/dt = ke0 * (Cp - Ce)

where Cp is plasma concentration (from PK model)

Effect:
E(Ce) = E0 + (Emax * Ce) / (EC50 + Ce)

Note: t1/2,ke0 = ln(2)/ke0 is the equilibration half-life.
When ke0 is large, effect follows plasma concentration closely (direct effect).
When ke0 is small, there is significant hysteresis between PK and PD.
"""
struct BiophaseEquilibration <: PDModelKind end

struct BiophaseEquilibrationParams
    ke0::Float64    # Effect site equilibration rate constant
    E0::Float64     # Baseline effect
    Emax::Float64   # Maximum effect
    EC50::Float64   # Effect site EC50
end

export IndirectResponseTurnover, IndirectResponseTurnoverParams

"""
Indirect response turnover PD model with inhibition of Kout.

States:
- R(t): response

Effect:
I(C) = (Imax * C) / (IC50 + C)

Dynamics:
dR/dt = Kin - Kout * (1 - I(C)) * R
"""
struct IndirectResponseTurnover <: PDModelKind end

struct IndirectResponseTurnoverParams
    Kin::Float64
    Kout::Float64
    R0::Float64
    Imax::Float64
    IC50::Float64
end

# Alias for backward compatibility (IRM-III is inhibition of Kout)
const IndirectResponseIRM3 = IndirectResponseTurnover
const IndirectResponseIRM3Params = IndirectResponseTurnoverParams

export IndirectResponseIRM1, IndirectResponseIRM1Params
export IndirectResponseIRM2, IndirectResponseIRM2Params
export IndirectResponseIRM3, IndirectResponseIRM3Params
export IndirectResponseIRM4, IndirectResponseIRM4Params

"""
Indirect response model type I: Inhibition of Kin (production).

States:
- R(t): response

Dynamics:
dR/dt = Kin × (1 - I(C)) - Kout × R

where I(C) = (Imax × C) / (IC50 + C)

Effect: Drug inhibits production → Response decreases below baseline.

Clinical Applications:
- Corticosteroids on cortisol production
- Statins on cholesterol synthesis
- Immunosuppressants on cytokine production
"""
struct IndirectResponseIRM1 <: PDModelKind end

struct IndirectResponseIRM1Params
    Kin::Float64      # Zero-order production rate
    Kout::Float64     # First-order elimination rate constant
    R0::Float64       # Baseline response (= Kin/Kout at steady state)
    Imax::Float64     # Maximum inhibition [0, 1]
    IC50::Float64     # Concentration at 50% inhibition
end

"""
Indirect response model type II: Stimulation of Kin (production).

States:
- R(t): response

Dynamics:
dR/dt = Kin × (1 + S(C)) - Kout × R

where S(C) = (Smax × C) / (SC50 + C)

Effect: Drug stimulates production → Response increases above baseline.

Clinical Applications:
- EPO on red blood cell production
- G-CSF on neutrophil production
- Growth factors on tissue growth
"""
struct IndirectResponseIRM2 <: PDModelKind end

struct IndirectResponseIRM2Params
    Kin::Float64      # Zero-order production rate
    Kout::Float64     # First-order elimination rate constant
    R0::Float64       # Baseline response (= Kin/Kout at steady state)
    Smax::Float64     # Maximum stimulation (can exceed 1)
    SC50::Float64     # Concentration at 50% stimulation
end

"""
Indirect response model type IV: Stimulation of Kout (elimination).

States:
- R(t): response

Dynamics:
dR/dt = Kin - Kout × (1 + S(C)) × R

where S(C) = (Smax × C) / (SC50 + C)

Effect: Drug stimulates elimination → Response decreases below baseline.

Clinical Applications:
- Diuretics on sodium excretion
- Thyroid hormone on metabolic rate
- Laxatives on bowel motility
"""
struct IndirectResponseIRM4 <: PDModelKind end

struct IndirectResponseIRM4Params
    Kin::Float64      # Zero-order production rate
    Kout::Float64     # First-order elimination rate constant
    R0::Float64       # Baseline response (= Kin/Kout at steady state)
    Smax::Float64     # Maximum stimulation (can exceed 1)
    SC50::Float64     # Concentration at 50% stimulation
end

# -------------------------
# Transit Compartment PD
# -------------------------

export TransitCompartmentPD, TransitCompartmentPDParams

"""
Transit compartment PD model for signal transduction delays.

Implements a chain of transit compartments to model delayed drug effects.

States:
- A1, A2, ..., AN: Transit compartment amounts
- Effect = AN (final transit compartment)

Signal generation:
- Signal(C) = E0 + Emax × C^γ / (EC50^γ + C^γ)

Transit dynamics (N compartments):
- dA1/dt = ktr × (Signal(C) - A1)
- dAi/dt = ktr × (A(i-1) - Ai)  for i = 2..N
- Effect = AN

Mean transit time (MTT):
- MTT = (N + 1) / ktr

Clinical Applications:
- Delayed myelosuppression (neutropenia, thrombocytopenia)
- Delayed biomarker responses
- Signal transduction cascades
"""
struct TransitCompartmentPD <: PDModelKind end

struct TransitCompartmentPDParams
    N::Int              # Number of transit compartments (1-20 typical)
    ktr::Float64        # Transit rate constant (1/time)
    E0::Float64         # Baseline effect/signal
    Emax::Float64       # Maximum effect above baseline
    EC50::Float64       # Concentration at 50% of maximum effect
    gamma::Float64      # Hill coefficient (steepness)
end

# -------------------------
# Disease Progression PD
# -------------------------

export DiseaseProgressionPD, DiseaseProgressionPDParams
export GrowthModelType, LinearGrowth, AsymptoticGrowth, GompertzGrowth, LogisticGrowth, ExponentialGrowth

"""
Growth model types for disease progression.
"""
@enum GrowthModelType begin
    LinearGrowth        # dS/dt = alpha
    AsymptoticGrowth    # dS/dt = kgrow * (Smax - S)
    GompertzGrowth      # dS/dt = kgrow * S * log(Smax / S)
    LogisticGrowth      # dS/dt = kgrow * S * (1 - S/Smax)
    ExponentialGrowth   # dS/dt = kgrow * S
end

"""
Disease progression PD model with tumor growth dynamics.

Implements various growth models commonly used in oncology:
- Linear: dS/dt = alpha - Drug_effect
- Asymptotic: dS/dt = kgrow × (Smax - S) - Drug_effect
- Gompertz: dS/dt = kgrow × S × log(Smax/S) - Drug_effect
- Logistic: dS/dt = kgrow × S × (1 - S/Smax) - Drug_effect
- Exponential: dS/dt = kgrow × S - Drug_effect

Drug effect:
- Drug_effect = kdrug × C × S  (cytotoxic, cell-kill)

where:
- S: Tumor size/volume
- C: Drug concentration
- kgrow: Growth rate constant
- Smax: Maximum tumor size (carrying capacity)
- alpha: Linear growth rate
- kdrug: Drug-induced cell kill rate constant

Clinical Applications:
- Tumor growth modeling
- Oncology dose-response
- Survival analysis
"""
struct DiseaseProgressionPD <: PDModelKind
    growth_model::GrowthModelType
end

struct DiseaseProgressionPDParams
    S0::Float64         # Initial tumor size
    kgrow::Float64      # Growth rate constant
    Smax::Float64       # Maximum size (carrying capacity, for bounded models)
    alpha::Float64      # Linear growth rate (for LinearGrowth only)
    kdrug::Float64      # Drug-induced cell kill rate constant
end

# -------------------------
# Combination Effect Models
# -------------------------

export BlissIndependence, BlissIndependenceParams
export CompetitiveInhibition, CompetitiveInhibitionParams
export DrugInteraction, DrugInteractionParams

"""
Bliss Independence combination effect model.

Models the combined effect of two drugs assuming independent action.
Used when drugs act on independent pathways.

Combined Effect:
E_combined = E_A + E_B - E_A × E_B

For fractional effects (0 to 1):
- If E_A = 0.5 and E_B = 0.5: E_combined = 0.75
- Represents probability of effect if mechanisms are independent

Individual drug effects use sigmoid Emax:
E_A = Emax_A × C_A^γ_A / (EC50_A^γ_A + C_A^γ_A)
E_B = Emax_B × C_B^γ_B / (EC50_B^γ_B + C_B^γ_B)

Parameters:
- E0: Baseline effect (no drug)
- Emax_A, EC50_A, gamma_A: Drug A dose-response parameters
- Emax_B, EC50_B, gamma_B: Drug B dose-response parameters
- input_A, input_B: Observation keys for drug concentrations

Clinical Applications:
- Combination chemotherapy
- Antibacterial combinations
- Multi-target therapies
"""
struct BlissIndependence <: PDModelKind end

struct BlissIndependenceParams
    E0::Float64           # Baseline effect
    Emax_A::Float64       # Maximum effect for drug A
    EC50_A::Float64       # EC50 for drug A
    gamma_A::Float64      # Hill coefficient for drug A
    Emax_B::Float64       # Maximum effect for drug B
    EC50_B::Float64       # EC50 for drug B
    gamma_B::Float64      # Hill coefficient for drug B
    input_A::Symbol       # Observation key for drug A concentration
    input_B::Symbol       # Observation key for drug B concentration
end

"""
Competitive inhibition PD model.

Models the effect of a competitive inhibitor on drug action.
The inhibitor competes for the same binding site, shifting the EC50.

Effect equation:
E = E0 + Emax × C^γ / ((EC50 × (1 + I/Ki))^γ + C^γ)

where:
- C: Drug (agonist) concentration
- I: Inhibitor concentration
- Ki: Inhibitor binding constant
- EC50_apparent = EC50 × (1 + I/Ki)

Parameters:
- E0: Baseline effect
- Emax: Maximum effect
- EC50: Intrinsic EC50 (without inhibitor)
- gamma: Hill coefficient
- Ki: Inhibitor binding constant
- input_drug: Observation key for drug concentration
- input_inhibitor: Observation key for inhibitor concentration

Clinical Applications:
- Drug-drug interactions
- Receptor antagonism
- Enzyme inhibition
"""
struct CompetitiveInhibition <: PDModelKind end

struct CompetitiveInhibitionParams
    E0::Float64              # Baseline effect
    Emax::Float64            # Maximum effect
    EC50::Float64            # Intrinsic EC50 (without inhibitor)
    gamma::Float64           # Hill coefficient
    Ki::Float64              # Inhibitor binding constant
    input_drug::Symbol       # Observation key for drug concentration
    input_inhibitor::Symbol  # Observation key for inhibitor concentration
end

"""
Drug interaction model (Greco model) for synergy/antagonism.

Models combined drug effects with an interaction parameter (ψ) that
quantifies synergy or antagonism.

Greco Equation:
(C_A/EC50_A)^γ + (C_B/EC50_B)^γ + ψ×(C_A/EC50_A)×(C_B/EC50_B) = E/(Emax-E)

Simplified for symmetric case (γ=1):
E = Emax × (C_A/EC50_A + C_B/EC50_B + ψ×C_A×C_B/(EC50_A×EC50_B)) /
    (1 + C_A/EC50_A + C_B/EC50_B + ψ×C_A×C_B/(EC50_A×EC50_B))

Interaction parameter ψ:
- ψ > 0: Synergy (greater than additive effect)
- ψ = 0: Additivity (Loewe additivity)
- ψ < 0: Antagonism (less than additive effect)

Parameters:
- E0: Baseline effect
- Emax: Maximum combined effect
- EC50_A: EC50 for drug A alone
- EC50_B: EC50 for drug B alone
- psi: Interaction parameter (synergy/antagonism)
- input_A, input_B: Observation keys for drug concentrations

Clinical Applications:
- Synergy detection in combination therapy
- Drug interaction studies
- Isobologram analysis
"""
struct DrugInteraction <: PDModelKind end

struct DrugInteractionParams
    E0::Float64           # Baseline effect
    Emax::Float64         # Maximum combined effect
    EC50_A::Float64       # EC50 for drug A
    EC50_B::Float64       # EC50 for drug B
    psi::Float64          # Interaction parameter (synergy: >0, antagonism: <0)
    input_A::Symbol       # Observation key for drug A concentration
    input_B::Symbol       # Observation key for drug B concentration
end

# -------------------------
# Tolerance/Sensitization Models
# -------------------------

export ToleranceCounterRegulation, ToleranceCounterRegulationParams
export ReceptorRegulation, ReceptorRegulationParams

"""
Tolerance model with counter-regulatory feedback.

Models the development of tolerance through a feedback mechanism that
opposes the drug effect over time.

States:
- M: Moderator/feedback variable (0 at baseline, increases with drug effect)

Dynamics:
- dM/dt = kin_mod × E_drug - kout_mod × M

Drug effect (before moderation):
- E_drug = Emax × C^γ / (EC50^γ + C^γ)

Net effect (with tolerance):
- E_net = E0 + E_drug - alpha × M

Parameters:
- E0: Baseline effect (no drug)
- Emax: Maximum drug effect
- EC50: Concentration at 50% of maximum effect
- gamma: Hill coefficient
- kin_mod: Rate of moderator production (driven by drug effect)
- kout_mod: Rate of moderator elimination
- alpha: Moderator effect on response (feedback strength)

Clinical Applications:
- Opioid tolerance development
- Beta-blocker tolerance
- Benzodiazepine adaptation
"""
struct ToleranceCounterRegulation <: PDModelKind end

struct ToleranceCounterRegulationParams
    E0::Float64       # Baseline effect
    Emax::Float64     # Maximum drug effect
    EC50::Float64     # EC50 for drug effect
    gamma::Float64    # Hill coefficient
    kin_mod::Float64  # Moderator production rate constant
    kout_mod::Float64 # Moderator elimination rate constant
    alpha::Float64    # Feedback strength (moderator effect coefficient)
end

"""
Receptor regulation (up/down-regulation) tolerance model.

Models tolerance through changes in receptor density in response to
sustained drug exposure.

States:
- R: Receptor density (normalized, baseline = 1.0)

Dynamics:
- dR/dt = kreg × (R_baseline - R) + regulation_effect

Regulation effect (drug-induced change):
- Down-regulation: regulation_effect = -kdown × E_drug × R
- Up-regulation: regulation_effect = +kup × E_drug × (Rmax - R)

Net effect:
- E_net = E0 + R × E_drug  (receptor amplifies/attenuates drug effect)

Parameters:
- E0: Baseline effect
- Emax: Maximum drug effect on receptors
- EC50: EC50 for drug effect
- gamma: Hill coefficient
- R_baseline: Baseline receptor density (normalized, typically 1.0)
- kreg: Receptor return-to-baseline rate constant
- Rmax: Maximum receptor density (for up-regulation, typically > R_baseline)
- kchange: Rate of receptor change (kdown for down-regulation, kup for up-regulation)
- direction: :down for down-regulation, :up for up-regulation

Clinical Applications:
- Beta-receptor down-regulation with chronic agonist exposure
- Opioid receptor adaptation
- Hormone receptor regulation
"""
struct ReceptorRegulation <: PDModelKind end

struct ReceptorRegulationParams
    E0::Float64           # Baseline effect
    Emax::Float64         # Maximum drug effect
    EC50::Float64         # EC50 for drug effect
    gamma::Float64        # Hill coefficient
    R_baseline::Float64   # Baseline receptor density (normalized, typically 1.0)
    kreg::Float64         # Rate of return to baseline
    Rmax::Float64         # Maximum receptor density (for up-regulation)
    kchange::Float64      # Rate of receptor change (down or up)
    direction::Symbol     # :down or :up
end

"""
PD specification container.

input_observation:
- which observation key from the upstream system is used as input, usually :conc

output_observation:
- name of the produced PD observable, default :effect is typical
"""
struct PDSpec{K<:PDModelKind,P}
    kind::K
    name::String
    params::P
    input_observation::Symbol
    output_observation::Symbol
end

# -------------------------
# Population specifications
# -------------------------

export RandomEffectKind, LogNormalIIV, IIVSpec, PopulationSpec, IndividualCovariates
export OmegaMatrix, ensure_positive_definite_omega
export get_diagonal_omegas, get_correlation_matrix, has_correlations, get_omega_matrix

abstract type RandomEffectKind end

"""
Log-normal inter-individual variability (IIV).

Parameter transform:
theta_i = theta_pop * exp(eta_i)

eta_i ~ Normal(0, omega^2)
"""
struct LogNormalIIV <: RandomEffectKind end

"""
Full covariance matrix for random effects with parameter ordering.

Supports:
- Full covariance matrix with correlations between all random effects
- Positive definiteness enforcement
- Cholesky decomposition for efficient sampling
"""
struct OmegaMatrix
    # Parameter names in order (corresponds to matrix rows/columns)
    param_names::Vector{Symbol}

    # Full covariance matrix (n_params x n_params)
    # Diagonal elements are variances, off-diagonal are covariances
    matrix::Matrix{Float64}

    # Cholesky factor (lower triangular) for efficient sampling
    # η = L * z where z ~ N(0, I)
    cholesky_L::Matrix{Float64}

    function OmegaMatrix(param_names::Vector{Symbol}, matrix::Matrix{Float64})
        n = length(param_names)
        if size(matrix) != (n, n)
            error("Matrix dimensions $(size(matrix)) don't match parameter count $n")
        end

        # Ensure symmetric
        matrix_sym = (matrix + matrix') / 2

        # Ensure positive definite
        matrix_pd = _ensure_pd(matrix_sym)

        # Compute Cholesky decomposition
        L = try
            cholesky(Symmetric(matrix_pd)).L |> Matrix
        catch e
            # Fallback: use eigendecomposition for numerical stability
            eig = eigen(Symmetric(matrix_pd))
            vals = max.(eig.values, 1e-10)
            sqrt_cov = eig.vectors * Diagonal(sqrt.(vals))
            sqrt_cov
        end

        new(param_names, matrix_pd, L)
    end
end

"""
Ensure a matrix is positive definite by adding a small ridge if needed.
"""
function _ensure_pd(matrix::Matrix{Float64}; min_eigenvalue::Float64=1e-8)::Matrix{Float64}
    eigenvalues = eigvals(Symmetric(matrix))
    min_eig = minimum(real.(eigenvalues))

    if min_eig <= 0
        # Add ridge to make positive definite
        ridge = abs(min_eig) + min_eigenvalue
        return matrix + ridge * I
    end

    return matrix
end

"""
Public function to ensure positive definiteness of omega matrix.
"""
function ensure_positive_definite_omega(matrix::Matrix{Float64}; min_eigenvalue::Float64=1e-8)::Matrix{Float64}
    return _ensure_pd(matrix; min_eigenvalue=min_eigenvalue)
end

"""
Create OmegaMatrix from diagonal variances (backward compatibility).
"""
function OmegaMatrix(omegas::Dict{Symbol,Float64})
    param_names = collect(keys(omegas))
    sort!(param_names)  # Ensure consistent ordering
    n = length(param_names)

    # Create diagonal covariance matrix (variances on diagonal)
    matrix = zeros(n, n)
    for (i, p) in enumerate(param_names)
        # omegas stores standard deviations, so square for variance
        matrix[i, i] = omegas[p]^2
    end

    return OmegaMatrix(param_names, matrix)
end

"""
Get diagonal elements (standard deviations) as a Dict for backward compatibility.
"""
function get_diagonal_omegas(omega::OmegaMatrix)::Dict{Symbol,Float64}
    result = Dict{Symbol,Float64}()
    for (i, p) in enumerate(omega.param_names)
        result[p] = sqrt(omega.matrix[i, i])
    end
    return result
end

"""
Get correlation matrix from OmegaMatrix.
"""
function get_correlation_matrix(omega::OmegaMatrix)::Matrix{Float64}
    n = length(omega.param_names)
    corr = zeros(n, n)
    sds = sqrt.(diag(omega.matrix))

    for i in 1:n
        for j in 1:n
            if sds[i] > 1e-12 && sds[j] > 1e-12
                corr[i, j] = omega.matrix[i, j] / (sds[i] * sds[j])
            elseif i == j
                corr[i, j] = 1.0
            end
        end
    end

    return corr
end

"""
IIV specification for a set of parameters.

Supports both:
1. Diagonal-only IIV (backward compatible): omegas Dict{Symbol,Float64}
2. Full covariance IIV: omega_matrix OmegaMatrix

omegas:
- Dict mapping parameter symbol to omega (standard deviation of eta)
- Used for diagonal-only IIV (backward compatible)

omega_matrix:
- Full covariance matrix with correlations
- If provided, omegas is derived from diagonal for compatibility

seed:
- deterministic seed for RNG

n:
- number of individuals
"""
struct IIVSpec{K<:RandomEffectKind}
    kind::K
    omegas::Dict{Symbol,Float64}  # Backward compatible: diagonal elements (SDs)
    omega_matrix::Union{Nothing,OmegaMatrix}  # Full covariance (optional)
    seed::UInt64
    n::Int

    # Constructor with full covariance matrix
    function IIVSpec(kind::K, omega_matrix::OmegaMatrix, seed::UInt64, n::Int) where K<:RandomEffectKind
        omegas = get_diagonal_omegas(omega_matrix)
        new{K}(kind, omegas, omega_matrix, seed, n)
    end

    # Constructor with diagonal-only (backward compatible)
    function IIVSpec(kind::K, omegas::Dict{Symbol,Float64}, seed::UInt64, n::Int) where K<:RandomEffectKind
        new{K}(kind, omegas, nothing, seed, n)
    end
end

"""
Check if IIV has correlations (full covariance matrix).
"""
has_correlations(iiv::IIVSpec)::Bool = iiv.omega_matrix !== nothing

"""
Get the effective OmegaMatrix (creates one from diagonal if needed).
"""
function get_omega_matrix(iiv::IIVSpec)::OmegaMatrix
    if iiv.omega_matrix !== nothing
        return iiv.omega_matrix
    else
        return OmegaMatrix(iiv.omegas)
    end
end

export IOVSpec, OccasionDefinition

"""
OccasionDefinition defines how dosing occasions are determined.

Supported v1 mode:
- :dose_times -> each unique dose time strictly greater than t0 starts a new occasion
- t0 is occasion 1
"""
struct OccasionDefinition
    mode::Symbol
end

"""
IOV specification.

pis:
- Dict mapping parameter symbol to pi (std dev of kappa)

seed:
- deterministic seed for IOV RNG stream (separate from IIV)

occasion_def:
- how occasions are determined
"""
struct IOVSpec{K<:RandomEffectKind}
    kind::K
    pis::Dict{Symbol,Float64}
    seed::UInt64
    occasion_def::OccasionDefinition
end

export CovariateEffectKind, LinearCovariate, PowerCovariate, ExpCovariate
export CovariateEffect, CovariateModel

abstract type CovariateEffectKind end

"""
Linear covariate model:
theta_i = theta_pop * (1 + beta * (cov - ref))
"""
struct LinearCovariate <: CovariateEffectKind end

"""
Power covariate model:
theta_i = theta_pop * (cov / ref) ^ beta
"""
struct PowerCovariate <: CovariateEffectKind end

"""
Exponential covariate model:
theta_i = theta_pop * exp(beta * (cov - ref))
"""
struct ExpCovariate <: CovariateEffectKind end

struct CovariateEffect{K<:CovariateEffectKind}
    kind::K
    param::Symbol
    covariate::Symbol
    beta::Float64
    ref::Float64
end

"""
A set of covariate effects applied to parameters.
"""
struct CovariateModel
    name::String
    effects::Vector{CovariateEffect}
end

"""
Optional covariates per individual.

values:
- static covariates as Dict{Symbol, Float64}

time_varying:
- optional TimeVaryingCovariates for time-dependent covariates
"""
struct IndividualCovariates
    values::Dict{Symbol,Float64}
    time_varying::Union{Nothing,TimeVaryingCovariates}
end

"""
Population simulation specification.

base_model_spec:
- the typical value model (population thetas) with doses, grid, solver handled outside
iiv:
- optional IIV spec (can be nothing)
covariates:
- optional vector aligned with n individuals (can be empty)
"""
struct PopulationSpec{MS}
    base_model_spec::MS
    iiv::Union{Nothing,IIVSpec}
    iov::Union{Nothing,IOVSpec}
    covariate_model::Union{Nothing,CovariateModel}
    covariates::Vector{IndividualCovariates}
end
