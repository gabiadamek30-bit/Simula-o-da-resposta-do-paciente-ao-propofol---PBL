# Drug Interaction Model (Greco Model)
# Models combined drug effects with synergy/antagonism parameter
#
# Greco Equation for symmetric case (γ=1):
# E = Emax × (C_A/EC50_A + C_B/EC50_B + ψ×C_A×C_B/(EC50_A×EC50_B)) /
#     (1 + C_A/EC50_A + C_B/EC50_B + ψ×C_A×C_B/(EC50_A×EC50_B))
#
# Interaction parameter ψ:
# - ψ > 0: Synergy
# - ψ = 0: Additivity (Loewe additivity)
# - ψ < 0: Antagonism
#
# Clinical Use:
# - Synergy detection in combination therapy
# - Drug interaction studies
# - Isobologram analysis

export validate, evaluate

"""
    validate(spec::PDSpec{DrugInteraction, DrugInteractionParams})

Validate Drug Interaction (Greco) PD parameters.

# Validation Rules
- EC50_A, EC50_B: must be positive
- psi: can be any real number (positive=synergy, negative=antagonism)
"""
function validate(spec::PDSpec{DrugInteraction,DrugInteractionParams})
    p = spec.params

    _require_positive("EC50_A", p.EC50_A)
    _require_positive("EC50_B", p.EC50_B)

    # psi can be any real number, no validation needed

    return nothing
end

"""
    evaluate(spec::PDSpec{DrugInteraction,DrugInteractionParams}, observations::Dict{Symbol,Vector{Float64}}) -> Vector{Float64}

Evaluate Drug Interaction (Greco) model at given drug concentrations.

# Arguments
- `spec`: PD specification with Drug Interaction parameters
- `observations`: Dict with concentration time series for both drugs
  - Must contain keys matching `input_A` and `input_B` from params

# Returns
- Vector of combined effect values

# Effect Equation (Greco model, γ=1)
Let:
- a = C_A / EC50_A
- b = C_B / EC50_B
- interaction = ψ × a × b

E = E0 + Emax × (a + b + interaction) / (1 + a + b + interaction)
"""
function evaluate(spec::PDSpec{DrugInteraction,DrugInteractionParams}, observations::Dict{Symbol,Vector{Float64}})::Vector{Float64}
    p = spec.params

    # Get concentration vectors for both drugs
    conc_A = get(observations, p.input_A, Float64[])
    conc_B = get(observations, p.input_B, Float64[])

    if isempty(conc_A) || isempty(conc_B)
        error("Missing concentration data: input_A=$(p.input_A) or input_B=$(p.input_B)")
    end

    if length(conc_A) != length(conc_B)
        error("Concentration vectors must have same length: $(length(conc_A)) vs $(length(conc_B))")
    end

    n = length(conc_A)
    effect = Vector{Float64}(undef, n)

    for i in 1:n
        C_A = max(conc_A[i], 0.0)
        C_B = max(conc_B[i], 0.0)

        # Normalized concentrations
        a = C_A / p.EC50_A
        b = C_B / p.EC50_B

        # Interaction term (Greco model)
        interaction = p.psi * a * b

        # Combined effect
        numerator = a + b + interaction
        denominator = 1.0 + a + b + interaction

        if denominator <= 0.0
            # Handle edge case (strong antagonism can cause negative denominator)
            effect[i] = p.E0
        else
            effect[i] = p.E0 + p.Emax * numerator / denominator
        end
    end

    return effect
end
