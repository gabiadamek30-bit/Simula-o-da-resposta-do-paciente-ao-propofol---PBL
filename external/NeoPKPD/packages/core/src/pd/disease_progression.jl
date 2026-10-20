# Disease Progression PD Model
# Implements tumor growth dynamics with drug effect
#
# Growth models:
# - Linear: dS/dt = alpha
# - Asymptotic: dS/dt = kgrow * (Smax - S)
# - Gompertz: dS/dt = kgrow * S * log(Smax / S)
# - Logistic: dS/dt = kgrow * S * (1 - S/Smax)
# - Exponential: dS/dt = kgrow * S
#
# Drug effect: Drug_effect = kdrug * C * S (cell-kill model)

export validate

"""
    validate(spec::PDSpec{DiseaseProgressionPD, DiseaseProgressionPDParams})

Validate disease progression PD parameters.

# Validation Rules
- S0: must be positive (initial tumor size)
- kgrow: must be non-negative (growth rate)
- Smax: must be positive and > S0 for bounded models
- alpha: must be non-negative (linear growth rate)
- kdrug: must be non-negative (drug effect)
"""
function validate(spec::PDSpec{DiseaseProgressionPD,DiseaseProgressionPDParams})
    p = spec.params
    model = spec.kind.growth_model

    _require_positive("S0", p.S0)

    if p.kgrow < 0.0
        error("kgrow must be non-negative, got $(p.kgrow)")
    end

    if p.kdrug < 0.0
        error("kdrug must be non-negative, got $(p.kdrug)")
    end

    # Model-specific validation
    if model == LinearGrowth
        if p.alpha < 0.0
            error("alpha must be non-negative for LinearGrowth, got $(p.alpha)")
        end
    elseif model in (AsymptoticGrowth, GompertzGrowth, LogisticGrowth)
        _require_positive("Smax", p.Smax)
        if p.Smax <= p.S0
            @warn "Smax ($(p.Smax)) <= S0 ($(p.S0)). Tumor may shrink even without drug."
        end
    end

    return nothing
end

"""
    growth_rate(S, model, kgrow, Smax, alpha) -> Float64

Calculate the natural growth rate for the specified growth model.
"""
function growth_rate(S::Float64, model::GrowthModelType, kgrow::Float64, Smax::Float64, alpha::Float64)::Float64
    if model == LinearGrowth
        return alpha
    elseif model == AsymptoticGrowth
        return kgrow * (Smax - S)
    elseif model == GompertzGrowth
        # Gompertz: kgrow * S * ln(Smax/S)
        # Numerically stable for S close to Smax
        if S <= 0.0 || S >= Smax
            return 0.0
        end
        return kgrow * S * log(Smax / S)
    elseif model == LogisticGrowth
        return kgrow * S * (1.0 - S / Smax)
    elseif model == ExponentialGrowth
        return kgrow * S
    else
        error("Unknown growth model: $model")
    end
end

"""
    drug_effect(C, S, kdrug) -> Float64

Calculate the drug-induced cell kill effect.
Uses linear cell-kill model: Drug_effect = kdrug * C * S
"""
function drug_effect(C::Float64, S::Float64, kdrug::Float64)::Float64
    C <= 0.0 && return 0.0
    S <= 0.0 && return 0.0
    return kdrug * C * S
end
