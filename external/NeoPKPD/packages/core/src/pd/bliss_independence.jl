# Bliss Independence Combination Effect Model
# Models combined drug effects assuming independent action
#
# Combined Effect: E = E_A + E_B - E_A × E_B
#
# Clinical Use:
# - Combination chemotherapy
# - Antibacterial combinations
# - Multi-target therapies

export validate, evaluate

"""
    validate(spec::PDSpec{BlissIndependence, BlissIndependenceParams})

Validate Bliss Independence PD parameters.

# Validation Rules
- EC50_A, EC50_B: must be positive
- gamma_A, gamma_B: must be positive
- Emax_A, Emax_B: typically [0, 1] for fractional effect, but can be any value
"""
function validate(spec::PDSpec{BlissIndependence,BlissIndependenceParams})
    p = spec.params

    # Validate drug A parameters
    _require_positive("EC50_A", p.EC50_A)
    if p.gamma_A <= 0.0
        error("gamma_A must be positive, got $(p.gamma_A)")
    end

    # Validate drug B parameters
    _require_positive("EC50_B", p.EC50_B)
    if p.gamma_B <= 0.0
        error("gamma_B must be positive, got $(p.gamma_B)")
    end

    return nothing
end

"""
    evaluate(spec::PDSpec{BlissIndependence,BlissIndependenceParams}, observations::Dict{Symbol,Vector{Float64}}) -> Vector{Float64}

Evaluate Bliss Independence model at given drug concentrations.

# Arguments
- `spec`: PD specification with Bliss Independence parameters
- `observations`: Dict with concentration time series for both drugs
  - Must contain keys matching `input_A` and `input_B` from params

# Returns
- Vector of combined effect values

# Effect Equation
E_A = Emax_A × C_A^γ_A / (EC50_A^γ_A + C_A^γ_A)
E_B = Emax_B × C_B^γ_B / (EC50_B^γ_B + C_B^γ_B)
E_combined = E0 + E_A + E_B - E_A × E_B
"""
function evaluate(spec::PDSpec{BlissIndependence,BlissIndependenceParams}, observations::Dict{Symbol,Vector{Float64}})::Vector{Float64}
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

        # Calculate fractional effects (normalized to [0, 1])
        # f = C^gamma / (EC50^gamma + C^gamma)
        if C_A <= 0.0
            f_A = 0.0
        else
            ratio_A = (p.EC50_A / C_A)^p.gamma_A
            f_A = 1.0 / (1.0 + ratio_A)
        end

        if C_B <= 0.0
            f_B = 0.0
        else
            ratio_B = (p.EC50_B / C_B)^p.gamma_B
            f_B = 1.0 / (1.0 + ratio_B)
        end

        # Bliss Independence for fractional effects: f_combined = f_A + f_B - f_A * f_B
        f_combined = f_A + f_B - f_A * f_B

        # Scale by maximum possible combined effect
        Emax_combined = max(p.Emax_A, p.Emax_B)
        effect[i] = p.E0 + f_combined * Emax_combined
    end

    return effect
end
