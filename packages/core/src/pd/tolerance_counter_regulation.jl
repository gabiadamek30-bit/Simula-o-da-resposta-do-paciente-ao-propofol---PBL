# Tolerance Counter-Regulation PD Model
# Models tolerance development through feedback mechanism
#
# States:
# - M: Moderator/feedback variable
#
# Dynamics:
# - dM/dt = kin_mod × E_drug - kout_mod × M
# - E_net = E0 + E_drug - alpha × M
#
# Clinical Use:
# - Opioid tolerance
# - Beta-blocker tolerance
# - Benzodiazepine adaptation

export validate

"""
    validate(spec::PDSpec{ToleranceCounterRegulation, ToleranceCounterRegulationParams})

Validate Tolerance Counter-Regulation PD parameters.

# Validation Rules
- EC50: must be positive
- gamma: must be positive
- kin_mod, kout_mod: must be positive (rate constants)
- alpha: must be non-negative (feedback strength)
"""
function validate(spec::PDSpec{ToleranceCounterRegulation,ToleranceCounterRegulationParams})
    p = spec.params

    _require_positive("EC50", p.EC50)

    if p.gamma <= 0.0
        error("gamma must be positive, got $(p.gamma)")
    end

    _require_positive("kin_mod", p.kin_mod)
    _require_positive("kout_mod", p.kout_mod)

    if p.alpha < 0.0
        error("alpha must be non-negative, got $(p.alpha)")
    end

    return nothing
end
