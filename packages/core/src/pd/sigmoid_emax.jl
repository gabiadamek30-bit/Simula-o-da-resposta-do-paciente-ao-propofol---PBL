"""
Sigmoid Emax (Hill equation) PD model implementation.

This model extends the hyperbolic Emax model with a Hill coefficient
(gamma) that controls the steepness of the concentration-response curve.

The Hill equation is one of the most widely used pharmacodynamic models
and can describe a wide range of dose-response relationships.

Key characteristics:
- gamma = 1: Standard Emax model (hyperbolic)
- gamma > 1: Steeper, more "switch-like" response (threshold effect)
- gamma < 1: More gradual response

The model can describe:
- Receptor binding with cooperativity (gamma ≠ 1)
- Enzyme inhibition/activation
- Ion channel effects
- Most clinical pharmacodynamic endpoints

Typical gamma values: 0.5 to 5 for most drugs
Steep drugs (e.g., neuromuscular blockers): gamma = 3-6

Reference: Hill AV, J Physiol 1910;40:iv-vii
"""

export validate, evaluate

function validate(spec::PDSpec{SigmoidEmax,SigmoidEmaxParams})
    p = spec.params

    # E0 can be any value (including 0 or negative for suppression endpoints)
    # Emax can be negative for inhibitory effects
    _require_positive("EC50", p.EC50)
    _require_positive("gamma", p.gamma)

    if p.gamma > 10.0
        @warn "Hill coefficient gamma > 10 may cause numerical instability"
    end

    return nothing
end

"""
Calculate Sigmoid Emax effect for a given concentration.

Effect(C) = E0 + (Emax * C^gamma) / (EC50^gamma + C^gamma)

For numerical stability with large gamma, we use:
Effect(C) = E0 + Emax / (1 + (EC50/C)^gamma)  when C > 0
Effect(0) = E0
"""
function sigmoid_emax_effect(C::Float64, p::SigmoidEmaxParams)
    if C <= 0.0
        return p.E0
    end

    # Use numerically stable form
    ratio = (p.EC50 / C)^p.gamma
    return p.E0 + p.Emax / (1.0 + ratio)
end

"""
Alternative calculation using the standard Hill equation form.
"""
function sigmoid_emax_effect_standard(C::Float64, E0::Float64, Emax::Float64, EC50::Float64, gamma::Float64)
    if C <= 0.0
        return E0
    end

    C_gamma = C^gamma
    EC50_gamma = EC50^gamma

    return E0 + (Emax * C_gamma) / (EC50_gamma + C_gamma)
end

"""
Calculate the slope of the concentration-response curve at EC50.

The slope at EC50 (n × Emax / (4 × EC50)) is a useful metric
for comparing the sensitivity of different drugs.
"""
function sigmoid_emax_slope_at_ec50(Emax::Float64, EC50::Float64, gamma::Float64)
    return gamma * Emax / (4.0 * EC50)
end

"""
Calculate the concentration at which a given fraction of Emax is achieved.

For fraction f (0 < f < 1):
C_f = EC50 * (f / (1 - f))^(1/gamma)

Common use: EC90 (f=0.9), EC10 (f=0.1)
"""
function concentration_at_fraction(EC50::Float64, gamma::Float64, fraction::Float64)
    if fraction <= 0.0 || fraction >= 1.0
        error("Fraction must be between 0 and 1 (exclusive)")
    end
    return EC50 * (fraction / (1.0 - fraction))^(1.0 / gamma)
end

"""
Evaluate the Sigmoid Emax model for a series of concentrations.

Returns effect values for each concentration in the input series.
"""
function evaluate(spec::PDSpec{SigmoidEmax,SigmoidEmaxParams}, input_series::Vector{Float64})
    validate(spec)

    p = spec.params

    out = Vector{Float64}(undef, length(input_series))
    for i in eachindex(input_series)
        C = input_series[i]
        out[i] = sigmoid_emax_effect(C, p)
    end
    return out
end
