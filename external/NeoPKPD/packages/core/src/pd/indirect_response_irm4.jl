# Indirect Response Model Type IV: Stimulation of Kout (Elimination)
# dR/dt = Kin - Kout × (1 + S(C)) × R

export validate

"""
    validate(spec::PDSpec{IndirectResponseIRM4, IndirectResponseIRM4Params})

Validate IRM-IV (Stimulation of Kout) parameters.

# Validation Rules
- Kin, Kout, R0, Smax, SC50: must be positive
- Warns if R0 differs from Kin/Kout by >1%
"""
function validate(spec::PDSpec{IndirectResponseIRM4,IndirectResponseIRM4Params})
    p = spec.params

    validate_irm_base_params(p.Kin, p.Kout, p.R0)
    validate_stimulation_params(p.Smax, p.SC50)

    return nothing
end
