# Indirect Response Model Type II: Stimulation of Kin (Production)
# dR/dt = Kin × (1 + S(C)) - Kout × R

export validate

"""
    validate(spec::PDSpec{IndirectResponseIRM2, IndirectResponseIRM2Params})

Validate IRM-II (Stimulation of Kin) parameters.

# Validation Rules
- Kin, Kout, R0, Smax, SC50: must be positive
- Warns if R0 differs from Kin/Kout by >1%
"""
function validate(spec::PDSpec{IndirectResponseIRM2,IndirectResponseIRM2Params})
    p = spec.params

    validate_irm_base_params(p.Kin, p.Kout, p.R0)
    validate_stimulation_params(p.Smax, p.SC50)

    return nothing
end
