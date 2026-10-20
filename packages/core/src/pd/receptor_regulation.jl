# Receptor Regulation PD Model
# Models tolerance through receptor up/down-regulation
#
# States:
# - R: Receptor density (normalized, baseline = 1.0)
#
# Dynamics (down-regulation):
# - dR/dt = kreg × (R_baseline - R) - kchange × E_drug × R
#
# Dynamics (up-regulation):
# - dR/dt = kreg × (R_baseline - R) + kchange × E_drug × (Rmax - R)
#
# Net effect:
# - E_net = E0 + R × E_drug
#
# Clinical Use:
# - Beta-receptor down-regulation
# - Opioid receptor adaptation
# - Hormone receptor regulation

export validate

"""
    validate(spec::PDSpec{ReceptorRegulation, ReceptorRegulationParams})

Validate Receptor Regulation PD parameters.

# Validation Rules
- EC50: must be positive
- gamma: must be positive
- R_baseline: must be positive (typically 1.0)
- kreg: must be positive (rate constant)
- Rmax: must be > 0 (for up-regulation, should be >= R_baseline)
- kchange: must be non-negative
- direction: must be :down or :up
"""
function validate(spec::PDSpec{ReceptorRegulation,ReceptorRegulationParams})
    p = spec.params

    _require_positive("EC50", p.EC50)

    if p.gamma <= 0.0
        error("gamma must be positive, got $(p.gamma)")
    end

    _require_positive("R_baseline", p.R_baseline)
    _require_positive("kreg", p.kreg)
    _require_positive("Rmax", p.Rmax)

    if p.kchange < 0.0
        error("kchange must be non-negative, got $(p.kchange)")
    end

    if p.direction != :down && p.direction != :up
        error("direction must be :down or :up, got $(p.direction)")
    end

    # Validate Rmax for up-regulation
    if p.direction == :up && p.Rmax < p.R_baseline
        @warn "Rmax ($(p.Rmax)) < R_baseline ($(p.R_baseline)) for up-regulation. " *
              "Receptor density cannot increase above Rmax."
    end

    return nothing
end
