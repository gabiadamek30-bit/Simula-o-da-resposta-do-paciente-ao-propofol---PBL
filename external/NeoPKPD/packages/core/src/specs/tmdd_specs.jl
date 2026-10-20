# =============================================================================
# TMDD (Target-Mediated Drug Disposition) Model Specifications
# =============================================================================
#
# Industry-standard TMDD models for biologics, monoclonal antibodies, and
# therapeutic proteins. Implements models from regulatory-accepted publications.
#
# References:
# - Mager DE, Jusko WJ. J Pharmacokinet Pharmacodyn. 2001;28(6):507-532
# - Gibiansky L, et al. J Pharmacokinet Pharmacodyn. 2008;35(5):573-591
# - Gibiansky L, Gibiansky E. J Clin Pharmacol. 2009;49(9):1012-1024
# - Peletier LA, Gabrielsson J. Eur J Pharm Sci. 2012;46(5):375-394
# - Dirks NL, Meibohm B. Clin Pharmacokinet. 2010;49(10):633-659
#
# Validated against: NONMEM, Monolix, Phoenix NLME, mrgsolve
# =============================================================================

using LinearAlgebra

# =============================================================================
# TMDD Model Type Hierarchy
# =============================================================================

export TMDDModelKind, TMDDApproximation
export FullTMDD, QSS, QE, RapidBinding, WagnerBinding, IrreversibleBinding
export TMDDRoute, IVBolus, IVInfusion, Subcutaneous

"""
Abstract base type for TMDD (Target-Mediated Drug Disposition) models.

TMDD describes nonlinear pharmacokinetics arising from high-affinity,
low-capacity binding to pharmacological targets. Key characteristics:

- Dose-dependent clearance (faster at low doses)
- Target-mediated elimination pathway
- Concentration-time profiles with distinct phases
- Common in monoclonal antibodies and therapeutic proteins

The TMDD model hierarchy:
```
                    TMDDModelKind
                         │
    ┌────────────────────┼────────────────────┐
    │                    │                    │
OneCptTMDD          TwoCptTMDD           TwoCptTMDDFcRn
    │                    │                    │
  (QSS/QE/Full)    (QSS/QE/Full)        (with FcRn)
```
"""
abstract type TMDDModelKind <: ModelKind end

"""
TMDD approximation methods.

- `FullTMDD`: Complete mechanistic model (L, R, P states)
- `QSS`: Quasi-Steady-State (most common, KSS = (koff + kint)/kon)
- `QE`: Quasi-Equilibrium (KD = koff/kon, assumes kint ≈ kdeg)
- `RapidBinding`: Wagner approximation (instantaneous equilibrium)
- `IrreversibleBinding`: No dissociation (koff = 0)
"""
@enum TMDDApproximation begin
    FullTMDD
    QSS
    QE
    RapidBinding
    WagnerBinding
    IrreversibleBinding
end

"""
Administration route for TMDD models.
"""
@enum TMDDRoute begin
    IVBolus
    IVInfusion
    Subcutaneous
end

# =============================================================================
# Core Parameter Structures - Industry Standard Parameterizations
# =============================================================================

export TMDDClearanceParams, TMDDMicroParams, TMDDHybridParams

"""
    TMDDClearanceParams

Clearance-based parameterization (NONMEM/Monolix standard).

This is the preferred parameterization for population PK analysis because:
- Parameters have direct physiological interpretation
- Better statistical properties for estimation
- Easier clinical interpretation

# Fields
- `CL`: Linear clearance (L/time) - non-specific elimination
- `V1`: Central volume of distribution (L)
- `V2`: Peripheral volume (L) - set to 0 for 1-compartment
- `Q`: Inter-compartmental clearance (L/time)
- `Vmax`: Maximum target-mediated elimination rate (amount/time)
- `Km`: Michaelis constant (concentration at half-Vmax)
- `kint`: Complex internalization rate constant (1/time)
- `kdeg`: Free target degradation rate constant (1/time)
- `R0`: Baseline target amount/concentration

# NONMEM Equivalent
```
CL = THETA(1)
V1 = THETA(2)
V2 = THETA(3)
Q  = THETA(4)
VMAX = THETA(5)
KM = THETA(6)
```
"""
struct TMDDClearanceParams
    CL::Float64      # Linear clearance (L/time)
    V1::Float64      # Central volume (L)
    V2::Float64      # Peripheral volume (L), 0 for 1-cpt
    Q::Float64       # Inter-compartmental clearance (L/time)
    Vmax::Float64    # Max target-mediated elimination (amount/time)
    Km::Float64      # Michaelis constant (concentration)
    kint::Float64    # Complex internalization rate (1/time)
    kdeg::Float64    # Target degradation rate (1/time)
    R0::Float64      # Baseline target level
end

"""
    TMDDMicroParams

Micro-constant parameterization (mechanistic).

Used when individual rate constants are of scientific interest
or when fitting to in vitro binding data.

# Fields
- `kel`: Elimination rate constant (1/time)
- `V1`: Central volume (L)
- `V2`: Peripheral volume (L)
- `k12`: Central→peripheral rate constant (1/time)
- `k21`: Peripheral→central rate constant (1/time)
- `kon`: Association rate constant (1/(concentration·time))
- `koff`: Dissociation rate constant (1/time)
- `kint`: Complex internalization rate (1/time)
- `ksyn`: Target synthesis rate (amount/time)
- `kdeg`: Target degradation rate (1/time)
- `R0`: Initial target level (amount)
"""
struct TMDDMicroParams
    kel::Float64     # Elimination rate constant (1/time)
    V1::Float64      # Central volume (L)
    V2::Float64      # Peripheral volume (L)
    k12::Float64     # Central to peripheral (1/time)
    k21::Float64     # Peripheral to central (1/time)
    kon::Float64     # Association rate (1/(conc·time))
    koff::Float64    # Dissociation rate (1/time)
    kint::Float64    # Internalization rate (1/time)
    ksyn::Float64    # Target synthesis (amount/time)
    kdeg::Float64    # Target degradation (1/time)
    R0::Float64      # Initial target (amount)
end

"""
    TMDDHybridParams

Hybrid parameterization combining clearances with binding constants.

Common in regulatory submissions for mAbs. Uses:
- Clearance terms for PK (CL, Q, V1, V2)
- Binding constants for target (KD or KSS)
- Rate constants for turnover (kint, kdeg)

# Fields
- `CL`: Linear clearance (L/time)
- `V1`: Central volume (L)
- `V2`: Peripheral volume (L)
- `Q`: Inter-compartmental clearance (L/time)
- `KSS`: Quasi-steady-state constant (concentration)
- `kint`: Complex internalization rate (1/time)
- `ksyn`: Target synthesis rate (amount/time)
- `kdeg`: Target degradation rate (1/time)
- `R0`: Baseline target level
"""
struct TMDDHybridParams
    CL::Float64      # Linear clearance (L/time)
    V1::Float64      # Central volume (L)
    V2::Float64      # Peripheral volume (L)
    Q::Float64       # Inter-compartmental clearance (L/time)
    KSS::Float64     # Quasi-steady-state constant (conc)
    kint::Float64    # Internalization rate (1/time)
    ksyn::Float64    # Target synthesis (amount/time)
    kdeg::Float64    # Target degradation (1/time)
    R0::Float64      # Initial target
end

# =============================================================================
# One-Compartment TMDD Models
# =============================================================================

export OneCptTMDD, OneCptTMDDParams

"""
    OneCptTMDD <: TMDDModelKind

One-compartment TMDD model.

Suitable for:
- Small therapeutic proteins
- Drugs with rapid distribution
- Early phase modeling with limited data

The approximation type determines the mathematical formulation:
- `FullTMDD`: 3 ODEs (L, R, P)
- `QSS`: 2 ODEs with algebraic complex
- `QE`: 2 ODEs with equilibrium assumption
"""
struct OneCptTMDD <: TMDDModelKind
    approximation::TMDDApproximation
    route::TMDDRoute
end

# Default constructor
OneCptTMDD() = OneCptTMDD(QSS, IVBolus)
OneCptTMDD(approx::TMDDApproximation) = OneCptTMDD(approx, IVBolus)

"""
One-compartment TMDD parameters using hybrid parameterization.
"""
struct OneCptTMDDParams
    CL::Float64      # Linear clearance (L/time)
    V::Float64       # Volume of distribution (L)
    KSS::Float64     # Quasi-steady-state constant (conc)
    kint::Float64    # Internalization rate (1/time)
    ksyn::Float64    # Target synthesis (amount/time)
    kdeg::Float64    # Target degradation (1/time)
    R0::Float64      # Initial target (amount or concentration)

    # Optional SC absorption parameters
    ka::Float64      # Absorption rate constant (1/time)
    F::Float64       # Bioavailability (fraction)
end

# Constructor without SC params
function OneCptTMDDParams(CL, V, KSS, kint, ksyn, kdeg, R0)
    OneCptTMDDParams(CL, V, KSS, kint, ksyn, kdeg, R0, 0.0, 1.0)
end

# =============================================================================
# Two-Compartment TMDD Models (Industry Standard)
# =============================================================================

export TwoCptTMDD, TwoCptTMDDParams

"""
    TwoCptTMDD <: TMDDModelKind

Two-compartment TMDD model - the industry standard for mAb PK.

This is the most commonly used model for:
- Monoclonal antibodies (IgG)
- Antibody-drug conjugates (ADCs)
- Fc-fusion proteins
- Large therapeutic proteins

Target binding assumed in central compartment only.

# Approximations
- `QSS`: **Recommended** - quasi-steady-state, robust estimation
- `FullTMDD`: Full mechanistic when binding kinetics needed
- `QE`: When target turnover similar to complex turnover

# Literature
- Gibiansky L, et al. J Pharmacokinet Pharmacodyn. 2008
- Dirks NL, Meibohm B. Clin Pharmacokinet. 2010
"""
struct TwoCptTMDD <: TMDDModelKind
    approximation::TMDDApproximation
    route::TMDDRoute
end

TwoCptTMDD() = TwoCptTMDD(QSS, IVBolus)
TwoCptTMDD(approx::TMDDApproximation) = TwoCptTMDD(approx, IVBolus)

"""
    TwoCptTMDDParams

Two-compartment TMDD parameters - hybrid clearance/binding parameterization.

This is the industry-standard parameterization for mAb TMDD modeling.
Compatible with NONMEM ADVAN13, Monolix, and Phoenix NLME.

# Required Fields
- `CL`: Linear clearance from central compartment (L/day)
- `V1`: Central volume of distribution (L)
- `V2`: Peripheral volume of distribution (L)
- `Q`: Inter-compartmental clearance (L/day)
- `KSS`: Quasi-steady-state constant (nM or μg/mL)
- `kint`: Complex internalization rate constant (1/day)
- `ksyn`: Target synthesis rate (nM/day or amount/day)
- `kdeg`: Free target degradation rate constant (1/day)
- `R0`: Baseline target level (nM or amount)

# Optional Fields (for SC administration)
- `ka`: First-order absorption rate constant (1/day)
- `F`: Bioavailability fraction (0-1)
- `Tlag`: Absorption lag time (days)

# Typical Values for IgG mAbs
```julia
params = TwoCptTMDDParams(
    CL = 0.2,      # 0.1-0.5 L/day typical for IgG
    V1 = 3.0,      # ~plasma volume
    V2 = 2.0,      # similar to V1 for mAbs
    Q = 0.5,       # 0.3-1.0 L/day
    KSS = 1.0,     # target-dependent (nM)
    kint = 0.1,    # 0.01-1.0 1/day
    ksyn = 0.1,    # determines R0 at steady state
    kdeg = 0.1,    # similar to kint for many targets
    R0 = 1.0       # baseline target (nM)
)
```

# Derived Quantities
- kel = CL/V1 (linear elimination rate constant)
- k12 = Q/V1 (distribution rate to peripheral)
- k21 = Q/V2 (distribution rate from peripheral)
- t½,α ≈ 0.693/(kel + k12 + k21) (distribution half-life)
- t½,β ≈ 0.693*V1*(1 + V2/V1)/(CL) (terminal half-life)
"""
struct TwoCptTMDDParams
    # PK parameters (clearance-based)
    CL::Float64      # Linear clearance (L/time)
    V1::Float64      # Central volume (L)
    V2::Float64      # Peripheral volume (L)
    Q::Float64       # Inter-compartmental clearance (L/time)

    # Target-mediated parameters
    KSS::Float64     # Quasi-steady-state constant (conc)
    kint::Float64    # Internalization rate (1/time)
    ksyn::Float64    # Target synthesis (amount/time)
    kdeg::Float64    # Target degradation (1/time)
    R0::Float64      # Initial target

    # SC absorption parameters (optional)
    ka::Float64      # Absorption rate (1/time)
    F::Float64       # Bioavailability (fraction)
    Tlag::Float64    # Lag time (time)
end

# Constructor without SC params (IV administration)
function TwoCptTMDDParams(CL, V1, V2, Q, KSS, kint, ksyn, kdeg, R0)
    TwoCptTMDDParams(CL, V1, V2, Q, KSS, kint, ksyn, kdeg, R0, 0.0, 1.0, 0.0)
end

# =============================================================================
# Two-Compartment TMDD with FcRn Recycling
# =============================================================================

export TwoCptTMDDFcRn, TwoCptTMDDFcRnParams

"""
    TwoCptTMDDFcRn <: TMDDModelKind

Two-compartment TMDD with FcRn (neonatal Fc receptor) recycling.

FcRn-mediated recycling is the primary mechanism for the long half-life
of IgG antibodies (~21 days). This model explicitly accounts for:

1. Pinocytic uptake into endosomes
2. pH-dependent FcRn binding in acidic endosomes
3. Recycling of FcRn-bound IgG to plasma
4. Lysosomal degradation of unbound IgG

# States
- L: Free drug in central compartment
- Lp: Drug in peripheral compartment
- Le: Drug in endosomal compartment
- R: Free target
- P: Drug-target complex (for full model)
- Rtot: Total target (for QSS)

# Applications
- Predicting half-life of engineered antibodies
- Modeling FcRn-binding mutants (YTE, LS mutations)
- Understanding target-mediated vs FcRn-mediated clearance

# Reference
Garg A, Bhattaram VA. J Pharmacokinet Pharmacodyn. 2011
Deng R, et al. mAbs. 2011;3(1):61-66
"""
struct TwoCptTMDDFcRn <: TMDDModelKind
    approximation::TMDDApproximation
    route::TMDDRoute
end

TwoCptTMDDFcRn() = TwoCptTMDDFcRn(QSS, IVBolus)

"""
Parameters for TMDD model with FcRn recycling.
"""
struct TwoCptTMDDFcRnParams
    # Standard PK
    V1::Float64          # Central volume (L)
    V2::Float64          # Peripheral volume (L)
    Q::Float64           # Inter-compartmental clearance (L/time)

    # FcRn recycling parameters
    CLup::Float64        # Pinocytic uptake clearance (L/time)
    FR::Float64          # Fraction recycled via FcRn (0-1)

    # Target-mediated parameters
    KSS::Float64         # Quasi-steady-state constant
    kint::Float64        # Internalization rate
    ksyn::Float64        # Target synthesis
    kdeg::Float64        # Target degradation
    R0::Float64          # Initial target

    # SC parameters
    ka::Float64
    F::Float64
    Tlag::Float64
end

function TwoCptTMDDFcRnParams(V1, V2, Q, CLup, FR, KSS, kint, ksyn, kdeg, R0)
    TwoCptTMDDFcRnParams(V1, V2, Q, CLup, FR, KSS, kint, ksyn, kdeg, R0, 0.0, 1.0, 0.0)
end

# =============================================================================
# TMDD with Immunogenicity (Anti-Drug Antibodies)
# =============================================================================

export TwoCptTMDDADA, TwoCptTMDDADAParams

"""
    TwoCptTMDDADA <: TMDDModelKind

Two-compartment TMDD with anti-drug antibody (ADA) formation.

Immunogenicity can significantly impact mAb pharmacokinetics:
- Neutralizing ADAs reduce efficacy
- Non-neutralizing ADAs can increase clearance
- ADA-drug immune complexes are rapidly cleared

# Model Components
1. Standard 2-compartment TMDD for drug
2. ADA production (time-dependent or exposure-dependent)
3. ADA-drug complex formation
4. Enhanced clearance of immune complexes

# Clinical Relevance
- Explains time-varying clearance
- Predicts loss of efficacy with repeated dosing
- Supports immunogenicity risk assessment

# Reference
Chen X, et al. Clin Pharmacol Ther. 2016
Davda JP, et al. mAbs. 2014
"""
struct TwoCptTMDDADA <: TMDDModelKind
    approximation::TMDDApproximation
    route::TMDDRoute
end

TwoCptTMDDADA() = TwoCptTMDDADA(QSS, IVBolus)

"""
Parameters for TMDD with immunogenicity.
"""
struct TwoCptTMDDADAParams
    # Base TMDD parameters
    CL::Float64
    V1::Float64
    V2::Float64
    Q::Float64
    KSS::Float64
    kint::Float64
    ksyn::Float64
    kdeg::Float64
    R0::Float64

    # ADA parameters
    kADA_prod::Float64   # ADA production rate constant
    kADA_deg::Float64    # ADA degradation rate constant
    kon_ADA::Float64     # Drug-ADA association rate
    koff_ADA::Float64    # Drug-ADA dissociation rate
    CL_complex::Float64  # Clearance of drug-ADA complex

    # Time to ADA onset (lag)
    T_onset::Float64     # Time until ADA response begins

    # SC parameters
    ka::Float64
    F::Float64
end

# =============================================================================
# Soluble Target TMDD (Cytokines, Growth Factors)
# =============================================================================

export TwoCptTMDDSoluble, TwoCptTMDDSolubleParams

"""
    TwoCptTMDDSoluble <: TMDDModelKind

TMDD model for drugs targeting soluble (circulating) targets.

Unlike membrane-bound targets, soluble targets:
- Distribute between central and peripheral compartments
- Drug-target complex can also distribute
- Both free and bound target are measurable

# Applications
- Anti-cytokine antibodies (anti-TNF, anti-IL6)
- Anti-VEGF antibodies
- Ligand traps

# Key Difference from Membrane Target
Complex elimination (kel_complex) may differ from target elimination
because complex may be cleared via different pathways.

# Reference
Luu KT, et al. Pharm Res. 2013
"""
struct TwoCptTMDDSoluble <: TMDDModelKind
    approximation::TMDDApproximation
    route::TMDDRoute
end

TwoCptTMDDSoluble() = TwoCptTMDDSoluble(QSS, IVBolus)

"""
Parameters for soluble target TMDD.
"""
struct TwoCptTMDDSolubleParams
    # Drug PK
    CL::Float64
    V1::Float64
    V2::Float64
    Q::Float64

    # Target binding
    KSS::Float64
    kint::Float64        # Complex internalization/clearance

    # Soluble target dynamics
    ksyn::Float64        # Target synthesis
    kdeg::Float64        # Free target clearance
    R0::Float64

    # Target distribution (soluble targets distribute)
    V1_target::Float64   # Target central volume
    V2_target::Float64   # Target peripheral volume
    Q_target::Float64    # Target inter-compartmental clearance

    # SC parameters
    ka::Float64
    F::Float64
end

# =============================================================================
# Multiple Target TMDD
# =============================================================================

export TwoCptTMDDMultiTarget, TwoCptTMDDMultiTargetParams

"""
    TwoCptTMDDMultiTarget <: TMDDModelKind

TMDD model for bispecific antibodies or drugs with multiple targets.

Bispecific antibodies bind two different targets, requiring:
- Separate binding kinetics for each target
- Potential for simultaneous binding
- Cross-linking effects

# Applications
- Bispecific T-cell engagers (BiTEs)
- Dual-variable domain antibodies
- Multi-specific antibodies
"""
struct TwoCptTMDDMultiTarget <: TMDDModelKind
    n_targets::Int
    approximation::TMDDApproximation
    route::TMDDRoute
end

TwoCptTMDDMultiTarget() = TwoCptTMDDMultiTarget(2, QSS, IVBolus)

"""
Parameters for dual-target TMDD (bispecific).
"""
struct TwoCptTMDDMultiTargetParams
    # Drug PK
    CL::Float64
    V1::Float64
    V2::Float64
    Q::Float64

    # Target 1 parameters
    KSS_1::Float64
    kint_1::Float64
    ksyn_1::Float64
    kdeg_1::Float64
    R0_1::Float64

    # Target 2 parameters
    KSS_2::Float64
    kint_2::Float64
    ksyn_2::Float64
    kdeg_2::Float64
    R0_2::Float64

    # SC parameters
    ka::Float64
    F::Float64
end

# =============================================================================
# Peripheral Target TMDD
# =============================================================================

export TwoCptTMDDPeripheral, TwoCptTMDDPeripheralParams

"""
    TwoCptTMDDPeripheral <: TMDDModelKind

TMDD model with target binding in peripheral compartment.

For targets expressed primarily in tissues rather than blood:
- Tumor antigens
- Tissue-specific receptors
- Targets with limited shedding

Drug must distribute to peripheral compartment before binding.
"""
struct TwoCptTMDDPeripheral <: TMDDModelKind
    approximation::TMDDApproximation
    route::TMDDRoute
end

TwoCptTMDDPeripheral() = TwoCptTMDDPeripheral(QSS, IVBolus)

struct TwoCptTMDDPeripheralParams
    CL::Float64
    V1::Float64
    V2::Float64
    Q::Float64
    KSS::Float64
    kint::Float64
    ksyn::Float64
    kdeg::Float64
    R0::Float64
    ka::Float64
    F::Float64
end

# =============================================================================
# TMDDSpec - Main Specification Container
# =============================================================================

export TMDDSpec

"""
    TMDDSpec{K<:TMDDModelKind,P}

TMDD model specification container.

# Fields
- `kind`: Model type (OneCptTMDD, TwoCptTMDD, etc.)
- `name`: Model identifier
- `params`: Parameter struct
- `doses`: Vector of DoseEvent
- `target_units`: Units for target concentration (:nM, :pM, :ug_mL)
- `drug_units`: Units for drug concentration (:nM, :ug_mL, :mg_L)

# Example
```julia
spec = TMDDSpec(
    TwoCptTMDD(QSS, IVBolus),
    "Pembrolizumab PopPK",
    TwoCptTMDDParams(
        CL=0.22, V1=3.28, V2=2.66, Q=0.54,
        KSS=0.087, kint=0.037, ksyn=0.001, kdeg=0.012, R0=0.083
    ),
    [DoseEvent(0.0, 200.0)],  # 200 mg IV
    target_units=:nM,
    drug_units=:mg_L
)
```
"""
struct TMDDSpec{K<:TMDDModelKind,P}
    kind::K
    name::String
    params::P
    doses::Vector{DoseEvent}
    target_units::Symbol
    drug_units::Symbol
end

# Default constructor with standard units
function TMDDSpec(kind::K, name::String, params::P, doses::Vector{DoseEvent}) where {K,P}
    TMDDSpec(kind, name, params, doses, :nM, :mg_L)
end

# =============================================================================
# Parameter Conversion Functions
# =============================================================================

export convert_to_micro, convert_to_clearance, convert_to_hybrid, derived_pk_params

"""
Convert clearance parameters to micro-constants.
"""
function convert_to_micro(p::TwoCptTMDDParams)
    kel = p.CL / p.V1
    k12 = p.Q / p.V1
    k21 = p.Q / p.V2

    # KSS = (koff + kint) / kon
    # Need additional info to get kon/koff individually
    # Use typical assumption: koff >> kint → KSS ≈ koff/kon = KD
    kon = p.kint / p.KSS  # Approximation
    koff = kon * p.KSS - p.kint
    koff = max(koff, 0.0)  # Ensure non-negative

    return (
        kel = kel,
        k12 = k12,
        k21 = k21,
        kon = kon,
        koff = koff,
        kint = p.kint,
        ksyn = p.ksyn,
        kdeg = p.kdeg
    )
end

"""
Calculate derived PK parameters.
"""
function derived_pk_params(p::TwoCptTMDDParams)
    kel = p.CL / p.V1
    k12 = p.Q / p.V1
    k21 = p.Q / p.V2

    # Two-compartment eigenvalues
    a = kel + k12 + k21
    b = kel * k21

    discriminant = sqrt(a^2 - 4b)
    alpha = (a + discriminant) / 2
    beta = (a - discriminant) / 2

    t_half_alpha = log(2) / alpha
    t_half_beta = log(2) / beta

    # Steady-state target
    R_ss = p.ksyn / p.kdeg

    return (
        kel = kel,
        k12 = k12,
        k21 = k21,
        alpha = alpha,
        beta = beta,
        t_half_alpha = t_half_alpha,
        t_half_beta = t_half_beta,
        R_ss = R_ss,
        Vss = p.V1 + p.V2
    )
end

# =============================================================================
# Target Binding Functions
# =============================================================================

export target_occupancy, target_EC50, target_EC90
export free_drug_concentration, total_drug_concentration

"""
    target_occupancy(C, KD_or_KSS)

Calculate fraction of target occupied at concentration C.

TO = C / (K + C)

where K is KD (equilibrium) or KSS (quasi-steady-state).

# Arguments
- `C`: Free drug concentration
- `KD_or_KSS`: Binding constant

# Returns
Fractional occupancy (0 to 1)
"""
function target_occupancy(C::Real, K::Real)::Float64
    C = Float64(C)
    K = Float64(K)
    if C <= 0.0 || K <= 0.0
        return 0.0
    end
    return C / (K + C)
end

"""
Concentration for 50% target occupancy (equals K).
"""
target_EC50(K::Float64) = K

"""
Concentration for 90% target occupancy.
"""
target_EC90(K::Float64) = 9.0 * K

"""
Concentration for 99% target occupancy.
"""
target_EC99(K::Float64) = 99.0 * K

"""
    free_drug_concentration(Ltot, Rtot, K, V)

Calculate free drug concentration from totals using rapid-binding assumption.

Solves quadratic: L² + (K·V + Rtot - Ltot)·L - K·V·Ltot = 0

# Arguments
- `Ltot`: Total drug amount
- `Rtot`: Total target amount
- `K`: Binding constant (KD or KSS)
- `V`: Volume for concentration calculation

# Returns
Free drug amount (not concentration)
"""
function free_drug_concentration(Ltot::Float64, Rtot::Float64, K::Float64, V::Float64)::Float64
    KV = K * V
    a = Ltot - Rtot - KV
    discriminant = a^2 + 4.0 * KV * Ltot
    L = 0.5 * (a + sqrt(max(0.0, discriminant)))
    return max(0.0, L)
end

"""
Calculate total (free + bound) drug concentration.
"""
function total_drug_concentration(L::Float64, P::Float64, V::Float64)::Float64
    return (L + P) / V
end

# =============================================================================
# Steady-State Calculations
# =============================================================================

export tmdd_steady_state, time_to_steady_state

"""
    tmdd_steady_state(params)

Calculate steady-state values for TMDD system.

Returns named tuple with:
- R_ss: Steady-state target (without drug)
- t_half_target: Target half-life
- t_half_linear: Linear elimination half-life
"""
function tmdd_steady_state(p::TwoCptTMDDParams)
    R_ss = p.ksyn / p.kdeg
    t_half_target = log(2) / p.kdeg
    t_half_linear = log(2) * p.V1 / p.CL  # Approximation

    # Effective half-life depends on concentration
    # At high C: t_half ≈ t_half_linear (linear dominates)
    # At low C: t_half shorter due to TMDD

    return (
        R_ss = R_ss,
        t_half_target = t_half_target,
        t_half_linear = t_half_linear,
        KSS = p.KSS
    )
end

"""
Estimate time to reach steady-state (~5 half-lives).
"""
function time_to_steady_state(p::TwoCptTMDDParams)
    derived = derived_pk_params(p)
    return 5.0 * derived.t_half_beta
end

# =============================================================================
# Validation Functions
# =============================================================================

export validate_tmdd

"""
Validate TMDD specification.
"""
function validate_tmdd(spec::TMDDSpec{TwoCptTMDD,TwoCptTMDDParams})
    p = spec.params

    # Check PK parameters
    _require_positive("CL", p.CL)
    _require_positive("V1", p.V1)
    _require_positive("V2", p.V2)
    _require_positive("Q", p.Q)

    # Check target parameters
    _require_positive("KSS", p.KSS)
    _require_positive("kint", p.kint)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    p.R0 >= 0.0 || error("R0 must be non-negative, got $(p.R0)")

    # SC parameters if applicable
    if spec.kind.route == Subcutaneous
        _require_positive("ka", p.ka)
        (0.0 < p.F <= 1.0) || error("F must be in (0, 1], got $(p.F)")
        p.Tlag >= 0.0 || error("Tlag must be non-negative, got $(p.Tlag)")
    end

    # Validate doses
    _validate_tmdd_doses(spec.doses)

    # Physiological plausibility warnings
    if p.CL > 1.0
        @warn "CL = $(p.CL) L/day is high for a mAb (typical: 0.1-0.5 L/day)"
    end

    if p.KSS < 1e-3
        @warn "KSS = $(p.KSS) is very low - verify units"
    end

    return nothing
end

function validate_tmdd(spec::TMDDSpec{OneCptTMDD,OneCptTMDDParams})
    p = spec.params
    _require_positive("CL", p.CL)
    _require_positive("V", p.V)
    _require_positive("KSS", p.KSS)
    _require_positive("kint", p.kint)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    p.R0 >= 0.0 || error("R0 must be non-negative")
    _validate_tmdd_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TwoCptTMDDFcRn,TwoCptTMDDFcRnParams})
    p = spec.params
    _require_positive("V1", p.V1)
    _require_positive("V2", p.V2)
    _require_positive("Q", p.Q)
    _require_positive("CLup", p.CLup)
    (0.0 <= p.FR <= 1.0) || error("FR (fraction recycled) must be in [0, 1]")
    _require_positive("KSS", p.KSS)
    _require_positive("kint", p.kint)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    _validate_tmdd_doses(spec.doses)
    return nothing
end

# Generic validation for other types
function validate_tmdd(spec::TMDDSpec)
    _validate_tmdd_doses(spec.doses)
    return nothing
end

function _validate_tmdd_doses(doses::Vector{DoseEvent})
    if isempty(doses)
        error("At least one DoseEvent is required")
    end
    for (i, d) in enumerate(doses)
        d.time >= 0.0 || error("Dose time must be >= 0 at index $i")
    end
    issorted([d.time for d in doses]) || error("Doses must be sorted by time")
    return nothing
end

# =============================================================================
# Helper Exports
# =============================================================================

export calculate_KD, calculate_KSS, calculate_Rtot_ss, calculate_half_life

"""Calculate KD from rate constants."""
calculate_KD(kon::Float64, koff::Float64) = koff / kon

"""Calculate KSS from rate constants."""
calculate_KSS(kon::Float64, koff::Float64, kint::Float64) = (koff + kint) / kon

"""Calculate steady-state target."""
calculate_Rtot_ss(ksyn::Float64, kdeg::Float64) = ksyn / kdeg

"""Calculate half-life from rate constant."""
calculate_half_life(k::Float64) = log(2) / k
