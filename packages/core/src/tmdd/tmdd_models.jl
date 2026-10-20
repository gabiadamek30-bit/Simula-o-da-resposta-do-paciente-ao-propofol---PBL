# =============================================================================
# TMDD Model ODE Systems
# =============================================================================
#
# Industry-standard ODE implementations for TMDD models.
# Optimized for numerical stability with stiff solvers.
#
# All ODEs use amount-based formulations to avoid division by zero.
# Concentrations are computed in the observation functions.
# =============================================================================

using DifferentialEquations
using SciMLBase

# =============================================================================
# One-Compartment TMDD ODEs
# =============================================================================

"""
One-compartment QSS TMDD ODE system.

States:
- u[1] = L: Free drug amount in central
- u[2] = Rtot: Total target amount (R + complex)

Parameters accessed via params struct (OneCptTMDDParams)
"""
function _ode_1cpt_tmdd_qss!(du, u, p, t)
    L = max(u[1], 0.0)      # Free drug
    Rtot = max(u[2], 0.0)   # Total target

    # Parameters
    CL = p.CL
    V = p.V
    KSS = p.KSS
    kint = p.kint
    ksyn = p.ksyn
    kdeg = p.kdeg

    # Derived
    kel = CL / V
    KSS_V = KSS * V

    # QSS approximation: complex at quasi-steady-state
    # RL = L * Rtot / (KSS*V + L)
    denom = KSS_V + L
    RL = (denom > 0.0) ? (L * Rtot / denom) : 0.0
    R = max(Rtot - RL, 0.0)

    # Free drug: linear elimination + target-mediated elimination
    du[1] = -kel * L - kint * RL

    # Total target: synthesis - free degradation - complex internalization
    du[2] = ksyn - kdeg * R - kint * RL

    return nothing
end

"""
One-compartment Full TMDD ODE system (3 states).

States:
- u[1] = L: Free drug amount
- u[2] = R: Free target amount
- u[3] = P: Complex amount
"""
function _ode_1cpt_tmdd_full!(du, u, p, t)
    L = max(u[1], 0.0)
    R = max(u[2], 0.0)
    P = max(u[3], 0.0)

    CL = p.CL
    V = p.V
    KSS = p.KSS
    kint = p.kint
    ksyn = p.ksyn
    kdeg = p.kdeg

    kel = CL / V

    # For full model, need kon and koff
    # Derive from KSS assuming rapid binding dominates: kon ≈ kint/KSS
    kon = kint / KSS
    koff = kon * KSS - kint
    koff = max(koff, 0.0)

    # Binding (concentration-based)
    C = L / V
    binding_rate = kon * C * R

    # Free drug
    du[1] = -kel * L - binding_rate + koff * P

    # Free target
    du[2] = ksyn - kdeg * R - binding_rate + koff * P

    # Complex
    du[3] = binding_rate - koff * P - kint * P

    return nothing
end

# =============================================================================
# Two-Compartment TMDD ODEs (Industry Standard)
# =============================================================================

"""
Two-compartment QSS TMDD ODE system - Industry Standard.

This is the most commonly used TMDD model for mAbs.

States:
- u[1] = L: Free drug in central compartment
- u[2] = Lp: Drug in peripheral compartment
- u[3] = Rtot: Total target (R + complex) in central

Parameters: TwoCptTMDDParams
"""
function _ode_2cpt_tmdd_qss!(du, u, p, t)
    L = max(u[1], 0.0)      # Free drug (central)
    Lp = max(u[2], 0.0)     # Drug (peripheral)
    Rtot = max(u[3], 0.0)   # Total target

    # PK parameters
    CL = p.CL
    V1 = p.V1
    V2 = p.V2
    Q = p.Q

    # Target parameters
    KSS = p.KSS
    kint = p.kint
    ksyn = p.ksyn
    kdeg = p.kdeg

    # Derived micro-constants
    kel = CL / V1
    k12 = Q / V1
    k21 = Q / V2
    KSS_V1 = KSS * V1

    # QSS: Complex at quasi-steady-state
    denom = KSS_V1 + L
    RL = (denom > 0.0) ? (L * Rtot / denom) : 0.0
    R = max(Rtot - RL, 0.0)

    # Free drug in central: elimination + TMDD + distribution
    du[1] = -kel * L - kint * RL - k12 * L + k21 * Lp

    # Drug in peripheral: distribution only (no target binding)
    du[2] = k12 * L - k21 * Lp

    # Total target: synthesis - degradation - internalization
    du[3] = ksyn - kdeg * R - kint * RL

    return nothing
end

"""
Two-compartment Full TMDD ODE system.

States:
- u[1] = L: Free drug in central
- u[2] = Lp: Drug in peripheral
- u[3] = R: Free target
- u[4] = P: Drug-target complex
"""
function _ode_2cpt_tmdd_full!(du, u, p, t)
    L = max(u[1], 0.0)
    Lp = max(u[2], 0.0)
    R = max(u[3], 0.0)
    P = max(u[4], 0.0)

    CL = p.CL
    V1 = p.V1
    V2 = p.V2
    Q = p.Q
    KSS = p.KSS
    kint = p.kint
    ksyn = p.ksyn
    kdeg = p.kdeg

    # Micro-constants
    kel = CL / V1
    k12 = Q / V1
    k21 = Q / V2

    # Binding rates (derive from KSS)
    kon = kint / KSS
    koff = max(kon * KSS - kint, 0.0)

    C = L / V1
    binding_rate = kon * C * R

    # Central drug
    du[1] = -kel * L - binding_rate + koff * P - k12 * L + k21 * Lp

    # Peripheral drug
    du[2] = k12 * L - k21 * Lp

    # Free target
    du[3] = ksyn - kdeg * R - binding_rate + koff * P

    # Complex
    du[4] = binding_rate - koff * P - kint * P

    return nothing
end

"""
Two-compartment QE (Quasi-Equilibrium) TMDD ODE system.

Uses KD (equilibrium constant) instead of KSS.
Assumes kint ≈ kdeg (similar turnover of complex and free target).
"""
function _ode_2cpt_tmdd_qe!(du, u, p, t)
    L = max(u[1], 0.0)
    Lp = max(u[2], 0.0)
    Rtot = max(u[3], 0.0)

    CL = p.CL
    V1 = p.V1
    V2 = p.V2
    Q = p.Q
    KSS = p.KSS  # Used as KD in QE
    kint = p.kint
    ksyn = p.ksyn
    kdeg = p.kdeg

    kel = CL / V1
    k12 = Q / V1
    k21 = Q / V2
    KD_V1 = KSS * V1  # KD approximation

    # QE: Complex at equilibrium
    denom = KD_V1 + L
    RL = (denom > 0.0) ? (L * Rtot / denom) : 0.0

    # Free drug
    du[1] = -kel * L - kint * RL - k12 * L + k21 * Lp

    # Peripheral drug
    du[2] = k12 * L - k21 * Lp

    # Total target (QE: same degradation rate for R and RL)
    du[3] = ksyn - kdeg * Rtot

    return nothing
end

"""
Two-compartment Rapid Binding (Wagner) TMDD ODE system.

Assumes instantaneous equilibrium. Tracks total drug and total target.

States:
- u[1] = Ltot: Total drug (L + P) in central
- u[2] = Lp: Drug in peripheral
- u[3] = Rtot: Total target (R + P)
"""
function _ode_2cpt_tmdd_rapid!(du, u, p, t)
    Ltot = max(u[1], 0.0)
    Lp = max(u[2], 0.0)
    Rtot = max(u[3], 0.0)

    CL = p.CL
    V1 = p.V1
    V2 = p.V2
    Q = p.Q
    KSS = p.KSS
    kint = p.kint
    ksyn = p.ksyn
    kdeg = p.kdeg

    kel = CL / V1
    k12 = Q / V1
    k21 = Q / V2
    KSS_V1 = KSS * V1

    # Calculate free drug from quadratic solution
    a = Ltot - Rtot - KSS_V1
    discriminant = a^2 + 4.0 * KSS_V1 * Ltot
    L = 0.5 * (a + sqrt(max(0.0, discriminant)))
    L = max(0.0, min(L, Ltot))  # Bound to physical limits

    P = Ltot - L
    R = max(Rtot - P, 0.0)

    # Total drug in central
    du[1] = -kel * L - kint * P - k12 * L + k21 * Lp

    # Drug in peripheral
    du[2] = k12 * L - k21 * Lp

    # Total target
    du[3] = ksyn - kdeg * R - kint * P

    return nothing
end

# =============================================================================
# Two-Compartment TMDD with FcRn Recycling
# =============================================================================

"""
Two-compartment TMDD with FcRn recycling.

States:
- u[1] = L: Free drug in central
- u[2] = Lp: Drug in peripheral
- u[3] = Le: Drug in endosomal compartment
- u[4] = Rtot: Total target

FcRn recycling mechanism:
- CLup: Pinocytic uptake into endosomes
- FR: Fraction recycled back to plasma via FcRn
- (1-FR): Fraction degraded in lysosomes
"""
function _ode_2cpt_tmdd_fcrn!(du, u, p, t)
    L = max(u[1], 0.0)
    Lp = max(u[2], 0.0)
    Le = max(u[3], 0.0)
    Rtot = max(u[4], 0.0)

    V1 = p.V1
    V2 = p.V2
    Q = p.Q
    CLup = p.CLup
    FR = p.FR
    KSS = p.KSS
    kint = p.kint
    ksyn = p.ksyn
    kdeg = p.kdeg

    k12 = Q / V1
    k21 = Q / V2

    # Endosomal rate constants
    kup = CLup / V1        # Uptake rate
    krec = 1.0             # Recycling rate (can be parameterized)
    kdeg_endo = krec * (1.0 - FR) / FR  # Degradation to balance recycling

    KSS_V1 = KSS * V1
    denom = KSS_V1 + L
    RL = (denom > 0.0) ? (L * Rtot / denom) : 0.0
    R = max(Rtot - RL, 0.0)

    # Central drug: -uptake + recycling - TMDD - distribution
    du[1] = -kup * L + krec * FR * Le - kint * RL - k12 * L + k21 * Lp

    # Peripheral drug
    du[2] = k12 * L - k21 * Lp

    # Endosomal drug: uptake - recycling - degradation
    du[3] = kup * L - krec * Le

    # Total target
    du[4] = ksyn - kdeg * R - kint * RL

    return nothing
end

# =============================================================================
# Two-Compartment TMDD with Subcutaneous Absorption
# =============================================================================

"""
Two-compartment QSS TMDD with SC absorption.

States:
- u[1] = A_depot: Amount in SC depot
- u[2] = L: Free drug in central
- u[3] = Lp: Drug in peripheral
- u[4] = Rtot: Total target
"""
function _ode_2cpt_tmdd_qss_sc!(du, u, p, t)
    A_depot = max(u[1], 0.0)
    L = max(u[2], 0.0)
    Lp = max(u[3], 0.0)
    Rtot = max(u[4], 0.0)

    CL = p.CL
    V1 = p.V1
    V2 = p.V2
    Q = p.Q
    KSS = p.KSS
    kint = p.kint
    ksyn = p.ksyn
    kdeg = p.kdeg
    ka = p.ka
    F = p.F

    kel = CL / V1
    k12 = Q / V1
    k21 = Q / V2
    KSS_V1 = KSS * V1

    # QSS complex
    denom = KSS_V1 + L
    RL = (denom > 0.0) ? (L * Rtot / denom) : 0.0
    R = max(Rtot - RL, 0.0)

    # SC depot (bioavailability applied at dosing)
    du[1] = -ka * A_depot

    # Central drug: absorption + elimination + TMDD + distribution
    du[2] = ka * A_depot - kel * L - kint * RL - k12 * L + k21 * Lp

    # Peripheral drug
    du[3] = k12 * L - k21 * Lp

    # Total target
    du[4] = ksyn - kdeg * R - kint * RL

    return nothing
end

# =============================================================================
# TMDD with Immunogenicity (ADA)
# =============================================================================

"""
Two-compartment TMDD with anti-drug antibody (ADA) formation.

States:
- u[1] = L: Free drug in central
- u[2] = Lp: Drug in peripheral
- u[3] = Rtot: Total target
- u[4] = ADA: Anti-drug antibody amount
- u[5] = LADA: Drug-ADA complex

ADA develops after T_onset and affects drug clearance.
"""
function _ode_2cpt_tmdd_ada!(du, u, p, t)
    L = max(u[1], 0.0)
    Lp = max(u[2], 0.0)
    Rtot = max(u[3], 0.0)
    ADA = max(u[4], 0.0)
    LADA = max(u[5], 0.0)

    CL = p.CL
    V1 = p.V1
    V2 = p.V2
    Q = p.Q
    KSS = p.KSS
    kint = p.kint
    ksyn = p.ksyn
    kdeg = p.kdeg

    kADA_prod = p.kADA_prod
    kADA_deg = p.kADA_deg
    kon_ADA = p.kon_ADA
    koff_ADA = p.koff_ADA
    CL_complex = p.CL_complex
    T_onset = p.T_onset

    kel = CL / V1
    k12 = Q / V1
    k21 = Q / V2
    KSS_V1 = KSS * V1

    # Target-mediated binding
    denom = KSS_V1 + L
    RL = (denom > 0.0) ? (L * Rtot / denom) : 0.0
    R = max(Rtot - RL, 0.0)

    # ADA binding
    C = L / V1
    ADA_conc = ADA / V1
    ADA_binding_rate = kon_ADA * C * ADA

    # ADA production (starts after T_onset)
    ADA_production = (t >= T_onset) ? kADA_prod : 0.0

    # Drug in central
    du[1] = -kel * L - kint * RL - ADA_binding_rate + koff_ADA * LADA - k12 * L + k21 * Lp

    # Drug in peripheral
    du[2] = k12 * L - k21 * Lp

    # Total target
    du[3] = ksyn - kdeg * R - kint * RL

    # ADA
    du[4] = ADA_production - kADA_deg * ADA - ADA_binding_rate + koff_ADA * LADA

    # Drug-ADA complex (cleared rapidly)
    kel_complex = CL_complex / V1
    du[5] = ADA_binding_rate - koff_ADA * LADA - kel_complex * LADA

    return nothing
end

# =============================================================================
# Model Interface Functions
# =============================================================================

export n_states, state_names, get_ode_function, get_initial_state
export dosing_compartment, get_central_volume

"""
Get number of ODE states for a TMDD model.
"""
function n_states(kind::OneCptTMDD)
    if kind.approximation == FullTMDD
        return 3  # L, R, P
    else
        return 2  # L, Rtot
    end
end

function n_states(kind::TwoCptTMDD)
    if kind.approximation == FullTMDD
        return 4  # L, Lp, R, P
    elseif kind.approximation == RapidBinding
        return 3  # Ltot, Lp, Rtot
    else
        return 3  # L, Lp, Rtot (QSS, QE)
    end
end

function n_states(kind::TwoCptTMDDFcRn)
    return 4  # L, Lp, Le, Rtot
end

function n_states(kind::TwoCptTMDDADA)
    return 5  # L, Lp, Rtot, ADA, LADA
end

"""
Get state variable names.
"""
function state_names(kind::OneCptTMDD)
    if kind.approximation == FullTMDD
        return [:L, :R, :P]
    else
        return [:L, :Rtot]
    end
end

function state_names(kind::TwoCptTMDD)
    if kind.approximation == FullTMDD
        return [:L, :Lp, :R, :P]
    elseif kind.approximation == RapidBinding
        return [:Ltot, :Lp, :Rtot]
    else
        return [:L, :Lp, :Rtot]
    end
end

function state_names(kind::TwoCptTMDDFcRn)
    return [:L, :Lp, :Le, :Rtot]
end

function state_names(kind::TwoCptTMDDADA)
    return [:L, :Lp, :Rtot, :ADA, :LADA]
end

# SC versions add depot
function n_states_with_sc(kind::TwoCptTMDD)
    return n_states(kind) + 1
end

function state_names_with_sc(kind::TwoCptTMDD)
    return vcat([:A_depot], state_names(kind))
end

"""
Get ODE function for a TMDD model.
"""
function get_ode_function(kind::OneCptTMDD)
    if kind.approximation == FullTMDD
        return _ode_1cpt_tmdd_full!
    else
        return _ode_1cpt_tmdd_qss!
    end
end

function get_ode_function(kind::TwoCptTMDD)
    if kind.route == Subcutaneous
        return _ode_2cpt_tmdd_qss_sc!  # SC version
    end

    if kind.approximation == FullTMDD
        return _ode_2cpt_tmdd_full!
    elseif kind.approximation == QE
        return _ode_2cpt_tmdd_qe!
    elseif kind.approximation == RapidBinding
        return _ode_2cpt_tmdd_rapid!
    else  # QSS (default)
        return _ode_2cpt_tmdd_qss!
    end
end

function get_ode_function(kind::TwoCptTMDDFcRn)
    return _ode_2cpt_tmdd_fcrn!
end

function get_ode_function(kind::TwoCptTMDDADA)
    return _ode_2cpt_tmdd_ada!
end

"""
Get initial state vector for a TMDD model.
"""
function get_initial_state(spec::TMDDSpec{OneCptTMDD,OneCptTMDDParams})
    p = spec.params
    if spec.kind.approximation == FullTMDD
        return [0.0, p.R0, 0.0]  # L=0, R=R0, P=0
    else
        return [0.0, p.R0]  # L=0, Rtot=R0
    end
end

function get_initial_state(spec::TMDDSpec{TwoCptTMDD,TwoCptTMDDParams})
    p = spec.params
    if spec.kind.route == Subcutaneous
        if spec.kind.approximation == FullTMDD
            return [0.0, 0.0, 0.0, p.R0, 0.0]  # Depot, L, Lp, R, P
        else
            return [0.0, 0.0, 0.0, p.R0]  # Depot, L, Lp, Rtot
        end
    else
        if spec.kind.approximation == FullTMDD
            return [0.0, 0.0, p.R0, 0.0]  # L, Lp, R, P
        elseif spec.kind.approximation == RapidBinding
            return [0.0, 0.0, p.R0]  # Ltot, Lp, Rtot
        else
            return [0.0, 0.0, p.R0]  # L, Lp, Rtot
        end
    end
end

function get_initial_state(spec::TMDDSpec{TwoCptTMDDFcRn,TwoCptTMDDFcRnParams})
    p = spec.params
    return [0.0, 0.0, 0.0, p.R0]  # L, Lp, Le, Rtot
end

function get_initial_state(spec::TMDDSpec{TwoCptTMDDADA,TwoCptTMDDADAParams})
    p = spec.params
    return [0.0, 0.0, p.R0, 0.0, 0.0]  # L, Lp, Rtot, ADA, LADA
end

"""
Get dosing compartment index.
"""
function dosing_compartment(kind::OneCptTMDD)
    return 1  # Central
end

function dosing_compartment(kind::TwoCptTMDD)
    if kind.route == Subcutaneous
        return 1  # SC depot
    else
        return 1  # Central (for IV)
    end
end

function dosing_compartment(kind::TwoCptTMDDFcRn)
    return 1  # Central
end

function dosing_compartment(kind::TwoCptTMDDADA)
    return 1  # Central
end

"""
Get central compartment index (for concentration calculation).
"""
function central_compartment(kind::TwoCptTMDD)
    if kind.route == Subcutaneous
        return 2  # After depot
    else
        return 1  # Direct IV
    end
end

"""
Get central volume from parameters.
"""
get_central_volume(p::OneCptTMDDParams) = p.V
get_central_volume(p::TwoCptTMDDParams) = p.V1
get_central_volume(p::TwoCptTMDDFcRnParams) = p.V1
get_central_volume(p::TwoCptTMDDADAParams) = p.V1
