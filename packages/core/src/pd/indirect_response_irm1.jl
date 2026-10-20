# Indirect Response Model Type I: Inhibition of Kin (Production)
# dR/dt = Kin × (1 - I(C)) - Kout × R

export validate

"""
    validate(spec::PDSpec{IndirectResponseIRM1, IndirectResponseIRM1Params})

Validate IRM-I (Inhibition of Kin) parameters.

# Validation Rules
- Kin, Kout, R0, IC50: must be positive
- Imax: must be in [0, 1]
- Warns if R0 differs from Kin/Kout by >1%
"""
function validate(spec::PDSpec{IndirectResponseIRM1,IndirectResponseIRM1Params})
    p = spec.params

    validate_irm_base_params(p.Kin, p.Kout, p.R0)
    validate_inhibition_params(p.Imax, p.IC50)

    return nothing
end
