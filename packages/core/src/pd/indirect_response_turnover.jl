# Indirect Response Model Type III (IRM-III): Inhibition of Kout
# dR/dt = Kin - Kout × (1 - I(C)) × R
# Note: inhibition() and stimulation() functions are in pd_helpers.jl

export validate

"""
    validate(spec::PDSpec{IndirectResponseTurnover, IndirectResponseTurnoverParams})

Validate IRM-III (Inhibition of Kout) parameters.
IndirectResponseTurnover is the original name for IRM-III model.

# Validation Rules
- Kin, Kout, R0, IC50: must be positive
- Imax: must be in [0, 1]
"""
function validate(spec::PDSpec{IndirectResponseTurnover,IndirectResponseTurnoverParams})
    p = spec.params

    validate_irm_base_params(p.Kin, p.Kout, p.R0)
    validate_inhibition_params(p.Imax, p.IC50)

    return nothing
end
