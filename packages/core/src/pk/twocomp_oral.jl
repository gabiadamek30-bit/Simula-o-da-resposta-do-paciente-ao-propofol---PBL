"""
Two-compartment oral first-order absorption PK model implementation.

This model extends the two-compartment model with first-order absorption
from a gut (depot) compartment. It captures both the absorption phase
and the bi-exponential disposition commonly seen with oral dosing.

The concentration-time profile shows:
- Absorption phase: Rising concentrations as drug absorbs from gut
- Distribution phase: Peak and initial decline as drug distributes
- Elimination phase: Terminal decline dominated by beta phase

Reference parameterization uses clearance terms for physiological interpretation.
"""

using DifferentialEquations
using SciMLBase

export validate

function validate(spec::ModelSpec{TwoCompOral,TwoCompOralParams})
    p = spec.params

    _require_positive("Ka", p.Ka)
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
Two-compartment oral first-order absorption ODE system.

States:
- u[1] = A_gut (amount in gut/depot compartment)
- u[2] = A_central (amount in central compartment)
- u[3] = A_peripheral (amount in peripheral compartment)

Doses are added to A_gut (first-order absorption into central).
"""
function _ode_twocomp_oral!(du, u, p, t)
    Ka = p.Ka
    CL = p.CL
    V1 = p.V1
    Q = p.Q
    V2 = p.V2

    A_gut = u[1]
    A_central = u[2]
    A_peripheral = u[3]

    # Micro-rate constants
    k10 = CL / V1
    k12 = Q / V1
    k21 = Q / V2

    # Gut compartment: first-order absorption
    du[1] = -Ka * A_gut

    # Central compartment: absorption input, elimination, distribution
    du[2] = Ka * A_gut - k10 * A_central - k12 * A_central + k21 * A_peripheral

    # Peripheral compartment: distribution only
    du[3] = k12 * A_central - k21 * A_peripheral

    return nothing
end
