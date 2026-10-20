"""
Three-compartment IV bolus PK model implementation.

This mammillary model describes drug distribution among:
- Central compartment (plasma/blood)
- Shallow peripheral compartment (rapidly equilibrating tissues)
- Deep peripheral compartment (slowly equilibrating tissues)

The tri-exponential concentration-time profile shows three phases:
- Alpha phase: Rapid initial decline (distribution to shallow peripheral)
- Beta phase: Intermediate decline (distribution to deep peripheral)
- Gamma phase: Terminal elimination phase

This model is commonly used for drugs with extensive tissue distribution,
such as many lipophilic compounds and some biologics.

Reference: Gibaldi & Perrier, Pharmacokinetics, 2nd ed.
"""

using DifferentialEquations
using SciMLBase

export validate

function validate(spec::ModelSpec{ThreeCompIVBolus,ThreeCompIVBolusParams})
    p = spec.params

    _require_positive("CL", p.CL)
    _require_positive("V1", p.V1)
    _require_positive("Q2", p.Q2)
    _require_positive("V2", p.V2)
    _require_positive("Q3", p.Q3)
    _require_positive("V3", p.V3)

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
Three-compartment IV bolus ODE system.

States:
- u[1] = A_central (amount in central compartment)
- u[2] = A_periph1 (amount in shallow peripheral compartment)
- u[3] = A_periph2 (amount in deep peripheral compartment)

Micro-rate constants:
- k10 = CL/V1 (elimination from central)
- k12 = Q2/V1 (central to shallow peripheral)
- k21 = Q2/V2 (shallow peripheral to central)
- k13 = Q3/V1 (central to deep peripheral)
- k31 = Q3/V3 (deep peripheral to central)
"""
function _ode_threecomp_ivbolus!(du, u, p, t)
    CL = p.CL
    V1 = p.V1
    Q2 = p.Q2
    V2 = p.V2
    Q3 = p.Q3
    V3 = p.V3

    A_central = u[1]
    A_periph1 = u[2]
    A_periph2 = u[3]

    # Micro-rate constants
    k10 = CL / V1
    k12 = Q2 / V1
    k21 = Q2 / V2
    k13 = Q3 / V1
    k31 = Q3 / V3

    # Central compartment: elimination + distribution to both peripheral compartments
    du[1] = -k10 * A_central - k12 * A_central + k21 * A_periph1 - k13 * A_central + k31 * A_periph2

    # Shallow peripheral compartment
    du[2] = k12 * A_central - k21 * A_periph1

    # Deep peripheral compartment
    du[3] = k13 * A_central - k31 * A_periph2

    return nothing
end
