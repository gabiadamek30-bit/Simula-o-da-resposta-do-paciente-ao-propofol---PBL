"""
Two-compartment IV bolus PK model implementation.

This model describes drug distribution between a central (plasma) compartment
and a peripheral (tissue) compartment, with elimination from the central compartment.

The bi-exponential concentration-time profile is characterized by:
- Distribution phase (alpha): Rapid initial decline as drug distributes to tissues
- Elimination phase (beta): Slower terminal decline as drug is eliminated

Reference parameterization uses clearance terms (CL, Q) rather than micro-constants
for better physiological interpretation and identifiability.
"""

using DifferentialEquations
using SciMLBase

export validate

function validate(spec::ModelSpec{TwoCompIVBolus,TwoCompIVBolusParams})
    p = spec.params

    _require_positive("CL", p.CL)
    _require_positive("V1", p.V1)
    _require_positive("Q", p.Q)
    _require_positive("V2", p.V2)

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
Two-compartment IV bolus ODE system.

States:
- u[1] = A_central (amount in central compartment)
- u[2] = A_peripheral (amount in peripheral compartment)

Micro-rate constants derived from clearance parameterization:
- k10 = CL/V1 (elimination)
- k12 = Q/V1 (central to peripheral)
- k21 = Q/V2 (peripheral to central)
"""
function _ode_twocomp_ivbolus!(du, u, p, t)
    CL = p.CL
    V1 = p.V1
    Q = p.Q
    V2 = p.V2

    A_central = u[1]
    A_peripheral = u[2]

    # Micro-rate constants
    k10 = CL / V1
    k12 = Q / V1
    k21 = Q / V2

    # Central compartment: elimination + distribution
    du[1] = -k10 * A_central - k12 * A_central + k21 * A_peripheral

    # Peripheral compartment: distribution only
    du[2] = k12 * A_central - k21 * A_peripheral

    return nothing
end
