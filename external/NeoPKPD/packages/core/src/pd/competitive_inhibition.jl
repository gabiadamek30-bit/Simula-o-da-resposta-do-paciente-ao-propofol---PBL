# Competitive Inhibition PD Model
# Models the effect of a competitive inhibitor on drug action
#
# Effect: E = E0 + Emax × C^γ / ((EC50 × (1 + I/Ki))^γ + C^γ)
#
# Clinical Use:
# - Drug-drug interactions
# - Receptor antagonism
# - Enzyme inhibition

export validate, evaluate

"""
    validate(spec::PDSpec{CompetitiveInhibition, CompetitiveInhibitionParams})

Validate Competitive Inhibition PD parameters.

# Validation Rules
- EC50: must be positive
- gamma: must be positive
- Ki: must be positive (inhibitor binding constant)
"""
function validate(spec::PDSpec{CompetitiveInhibition,CompetitiveInhibitionParams})
    p = spec.params

    _require_positive("EC50", p.EC50)
    _require_positive("Ki", p.Ki)

    if p.gamma <= 0.0
        error("gamma must be positive, got $(p.gamma)")
    end

    return nothing
end

"""
    evaluate(spec::PDSpec{CompetitiveInhibition,CompetitiveInhibitionParams}, observations::Dict{Symbol,Vector{Float64}}) -> Vector{Float64}

Evaluate Competitive Inhibition model at given drug and inhibitor concentrations.

# Arguments
- `spec`: PD specification with Competitive Inhibition parameters
- `observations`: Dict with concentration time series for drug and inhibitor
  - Must contain keys matching `input_drug` and `input_inhibitor` from params

# Returns
- Vector of effect values

# Effect Equation
EC50_apparent = EC50 × (1 + I/Ki)
E = E0 + Emax × C^γ / (EC50_apparent^γ + C^γ)
"""
function evaluate(spec::PDSpec{CompetitiveInhibition,CompetitiveInhibitionParams}, observations::Dict{Symbol,Vector{Float64}})::Vector{Float64}
    p = spec.params

    # Get concentration vectors for drug and inhibitor
    conc_drug = get(observations, p.input_drug, Float64[])
    conc_inhib = get(observations, p.input_inhibitor, Float64[])

    if isempty(conc_drug) || isempty(conc_inhib)
        error("Missing concentration data: input_drug=$(p.input_drug) or input_inhibitor=$(p.input_inhibitor)")
    end

    if length(conc_drug) != length(conc_inhib)
        error("Concentration vectors must have same length: $(length(conc_drug)) vs $(length(conc_inhib))")
    end

    n = length(conc_drug)
    effect = Vector{Float64}(undef, n)

    for i in 1:n
        C = max(conc_drug[i], 0.0)
        I = max(conc_inhib[i], 0.0)

        # Calculate apparent EC50 with competitive inhibition
        EC50_apparent = p.EC50 * (1.0 + I / p.Ki)

        # Calculate effect using sigmoid Emax with apparent EC50
        if C <= 0.0
            effect[i] = p.E0
        else
            C_gamma = C^p.gamma
            EC50_app_gamma = EC50_apparent^p.gamma
            effect[i] = p.E0 + p.Emax * C_gamma / (EC50_app_gamma + C_gamma)
        end
    end

    return effect
end
