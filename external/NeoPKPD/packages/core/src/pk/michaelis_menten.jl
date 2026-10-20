"""
One-compartment PK model with Michaelis-Menten (saturable) elimination.

This model describes nonlinear pharmacokinetics where the elimination
pathway becomes saturated at higher drug concentrations. This is common
for drugs eliminated by capacity-limited metabolic enzymes or transporters.

Key characteristics:
- At low concentrations (C << Km): Linear (first-order) elimination
  Apparent CL ≈ Vmax/Km
- At high concentrations (C >> Km): Zero-order (constant rate) elimination
  Rate ≈ Vmax
- At intermediate concentrations: Mixed-order kinetics

Clinical implications:
- Dose-dependent half-life (increases with dose)
- Disproportionate increase in AUC with dose
- Time to steady-state varies with dose

Common examples: Phenytoin, Ethanol, Aspirin (high doses)

Reference: Gibaldi & Perrier, Pharmacokinetics, 2nd ed., Chapter 10
"""

using DifferentialEquations
using SciMLBase

export validate

function validate(spec::ModelSpec{MichaelisMentenElimination,MichaelisMentenEliminationParams})
    p = spec.params

    _require_positive("Vmax", p.Vmax)
    _require_positive("Km", p.Km)
    _require_positive("V", p.V)

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
Michaelis-Menten elimination ODE system.

States:
- u[1] = A_central (amount in central compartment)

Dynamics:
C = A_central / V
dA_central/dt = -Vmax * C / (Km + C)
             = -Vmax * A_central / (Km * V + A_central)

Note: The rate of elimination is in amount/time (not concentration/time).
At C = Km, the elimination rate is Vmax/2.
"""
function _ode_michaelis_menten!(du, u, p, t)
    Vmax = p.Vmax
    Km = p.Km
    V = p.V

    A_central = u[1]

    # Michaelis-Menten elimination
    # Rate = Vmax * C / (Km + C) = Vmax * A / (Km * V + A)
    du[1] = -Vmax * A_central / (Km * V + A_central)

    return nothing
end

"""
Calculate apparent clearance at a given concentration.
At low C, this approaches Vmax/Km (intrinsic clearance).
"""
function apparent_clearance(C::Float64, Vmax::Float64, Km::Float64)
    return Vmax / (Km + C)
end

"""
Calculate apparent half-life at a given concentration.
This demonstrates the dose-dependent kinetics.
"""
function apparent_half_life(C::Float64, Vmax::Float64, Km::Float64, V::Float64)
    CL_app = apparent_clearance(C, Vmax, Km)
    return log(2) * V / CL_app
end
