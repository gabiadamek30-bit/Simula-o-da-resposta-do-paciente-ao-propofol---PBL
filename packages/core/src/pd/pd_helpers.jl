# Shared PD Helper Functions
# Common functions used across indirect response and other PD models

export inhibition, stimulation, emax_effect, sigmoid_emax_effect

"""
    inhibition(C, Imax, IC50) -> Float64

Calculate inhibition term I(C) ∈ [0, Imax] using Emax-type function.

Mathematical form:
    I(C) = (Imax × C) / (IC50 + C)

# Arguments
- `C::Float64`: Drug concentration
- `Imax::Float64`: Maximum inhibition (typically 0 ≤ Imax ≤ 1)
- `IC50::Float64`: Concentration at 50% of maximum inhibition

# Returns
- Inhibition value in [0, Imax]

# Clinical Applications
- IRM-I: Inhibition of production (Kin)
- IRM-III: Inhibition of loss (Kout)
"""
function inhibition(C::Float64, Imax::Float64, IC50::Float64)::Float64
    C <= 0.0 && return 0.0
    return (Imax * C) / (IC50 + C)
end

"""
    stimulation(C, Smax, SC50) -> Float64

Calculate stimulation term S(C) ∈ [0, Smax] using Emax-type function.

Mathematical form:
    S(C) = (Smax × C) / (SC50 + C)

# Arguments
- `C::Float64`: Drug concentration
- `Smax::Float64`: Maximum stimulation (Smax > 0, can exceed 1)
- `SC50::Float64`: Concentration at 50% of maximum stimulation

# Returns
- Stimulation value in [0, Smax]

# Clinical Applications
- IRM-II: Stimulation of production (Kin)
- IRM-IV: Stimulation of loss (Kout)
"""
function stimulation(C::Float64, Smax::Float64, SC50::Float64)::Float64
    C <= 0.0 && return 0.0
    return (Smax * C) / (SC50 + C)
end

"""
    emax_effect(C, E0, Emax, EC50) -> Float64

Calculate direct Emax effect.

Mathematical form:
    E(C) = E0 + (Emax × C) / (EC50 + C)

# Arguments
- `C::Float64`: Drug concentration
- `E0::Float64`: Baseline effect
- `Emax::Float64`: Maximum effect above baseline
- `EC50::Float64`: Concentration at 50% of maximum effect

# Returns
- Effect value
"""
function emax_effect(C::Float64, E0::Float64, Emax::Float64, EC50::Float64)::Float64
    C <= 0.0 && return E0
    return E0 + (Emax * C) / (EC50 + C)
end

"""
    sigmoid_emax_effect(C, E0, Emax, EC50, gamma) -> Float64

Calculate sigmoid Emax (Hill) effect.

Mathematical form:
    E(C) = E0 + (Emax × C^γ) / (EC50^γ + C^γ)

Uses numerically stable form to avoid overflow with large gamma values.

# Arguments
- `C::Float64`: Drug concentration
- `E0::Float64`: Baseline effect
- `Emax::Float64`: Maximum effect above baseline
- `EC50::Float64`: Concentration at 50% of maximum effect
- `gamma::Float64`: Hill coefficient (steepness parameter)

# Returns
- Effect value
"""
function sigmoid_emax_effect(C::Float64, E0::Float64, Emax::Float64, EC50::Float64, gamma::Float64)::Float64
    C <= 0.0 && return E0
    # Numerically stable form
    ratio = (EC50 / C)^gamma
    return E0 + Emax / (1.0 + ratio)
end

"""
    validate_irm_base_params(Kin, Kout, R0, prefix::String="")

Validate base parameters common to all IRM models.
"""
function validate_irm_base_params(Kin::Float64, Kout::Float64, R0::Float64; prefix::String="")
    _require_positive("$(prefix)Kin", Kin)
    _require_positive("$(prefix)Kout", Kout)
    _require_positive("$(prefix)R0", R0)

    # Check baseline consistency (warning only)
    expected_R0 = Kin / Kout
    if abs(R0 - expected_R0) / expected_R0 > 0.01
        @warn "R0 ($(R0)) differs from Kin/Kout ($(expected_R0)) by more than 1%. " *
              "Baseline may not be at steady state."
    end

    return nothing
end

"""
    validate_inhibition_params(Imax, IC50, prefix::String="")

Validate inhibition parameters.
"""
function validate_inhibition_params(Imax::Float64, IC50::Float64; prefix::String="")
    if Imax < 0.0 || Imax > 1.0
        error("$(prefix)Imax must be in [0, 1], got $(Imax)")
    end
    _require_positive("$(prefix)IC50", IC50)
    return nothing
end

"""
    validate_stimulation_params(Smax, SC50, prefix::String="")

Validate stimulation parameters.
Smax can be 0 (no stimulation) or positive, SC50 must be positive.
"""
function validate_stimulation_params(Smax::Float64, SC50::Float64; prefix::String="")
    if Smax < 0.0
        error("$(prefix)Smax must be non-negative, got $(Smax)")
    end
    _require_positive("$(prefix)SC50", SC50)
    return nothing
end
