# Transit Compartment PD Model
# Implements signal transduction delay through a chain of transit compartments
#
# Signal: S(C) = E0 + Emax × C^γ / (EC50^γ + C^γ)
# Transit chain:
#   dA1/dt = ktr × (Signal - A1)
#   dAi/dt = ktr × (A(i-1) - Ai)  for i = 2..N
# Effect = AN

export validate

"""
    validate(spec::PDSpec{TransitCompartmentPD, TransitCompartmentPDParams})

Validate transit compartment PD parameters.

# Validation Rules
- N: must be positive integer (1-20 typical)
- ktr: must be positive
- E0: can be any value (baseline)
- Emax: can be any value (can be negative for inhibition)
- EC50: must be positive
- gamma: must be positive
"""
function validate(spec::PDSpec{TransitCompartmentPD,TransitCompartmentPDParams})
    p = spec.params

    if p.N < 1
        error("N (number of transit compartments) must be >= 1, got $(p.N)")
    end
    if p.N > 20
        @warn "N = $(p.N) is unusually large for transit compartments. " *
              "Typical values are 1-10."
    end

    _require_positive("ktr", p.ktr)
    _require_positive("EC50", p.EC50)
    _require_positive("gamma", p.gamma)

    # Inform about MTT
    mtt = (p.N + 1) / p.ktr
    @info "Transit compartment model: N=$(p.N), ktr=$(p.ktr), MTT=$(round(mtt, digits=2)) time units"

    return nothing
end

"""
    transit_signal(C, E0, Emax, EC50, gamma) -> Float64

Calculate the input signal to the transit chain using sigmoid Emax model.
"""
function transit_signal(C::Float64, E0::Float64, Emax::Float64, EC50::Float64, gamma::Float64)::Float64
    C <= 0.0 && return E0
    return sigmoid_emax_effect(C, E0, Emax, EC50, gamma)
end
