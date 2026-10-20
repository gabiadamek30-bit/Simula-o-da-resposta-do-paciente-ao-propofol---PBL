using DifferentialEquations
using SciMLBase

export validate, _ode_onecomp_ivbolus!

function validate(spec::ModelSpec{OneCompIVBolus,OneCompIVBolusParams})
    CL = spec.params.CL
    V = spec.params.V

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

function _ode_onecomp_ivbolus!(dA, A, p, t)
    CL = p.CL
    V = p.V
    dA[1] = -(CL / V) * A[1]
    return nothing
end
