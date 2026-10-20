using DifferentialEquations
using SciMLBase

export validate, _ode_onecomp_oral_first_order!, _observe_conc_onecomp_oral_first_order

function validate(spec::ModelSpec{OneCompOralFirstOrder,OneCompOralFirstOrderParams})
    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V

    _require_positive("Ka", Ka)
    _require_positive("CL", CL)
    _require_positive("V", V)

    if isempty(spec.doses)
        error("At least one DoseEvent is required")
    end

    for (i, d) in enumerate(spec.doses)
        if d.time < 0.0
            error("DoseEvent time must be >= 0 at index $(i), got $(d.time)")
        end
        _require_positive("DoseEvent amount at index $(i)", d.amount)
    end

    if !issorted([d.time for d in spec.doses])
        error("Dose events must be sorted by time ascending")
    end

    return nothing
end

"""
Oral one-compartment first-order absorption model.

States:
- u[1] = Agut
- u[2] = Acent

Dynamics:
dAgut/dt = -Ka * Agut
dAcent/dt =  Ka * Agut - (CL/V) * Acent

Bolus doses add to Agut.
"""
function _ode_onecomp_oral_first_order!(du, u, p, t)
    Ka = p.Ka
    CL = p.CL
    V = p.V

    Agut = u[1]
    Acent = u[2]

    du[1] = -Ka * Agut
    du[2] = Ka * Agut - (CL / V) * Acent
    return nothing
end

function _observe_conc_onecomp_oral_first_order(u, p)
    V = p.V
    return u[2] / V
end
