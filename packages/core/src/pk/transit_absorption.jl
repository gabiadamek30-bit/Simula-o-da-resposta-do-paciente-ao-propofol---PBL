"""
Transit compartment absorption model implementation.

This model implements a chain of transit compartments to describe delayed
and/or complex absorption profiles. It provides a more physiological
representation of drug transit through the GI tract.

Key features:
- Flexible number of transit compartments (N)
- Uniform transit rate constant (Ktr) through the chain
- Final absorption via Ka into the central compartment
- Mean transit time (MTT) ≈ (N+1)/Ktr

This approach can capture:
- Delayed absorption (lag time)
- Complex absorption profiles
- Enterohepatic recirculation (with modifications)
- Gastric emptying effects

Reference: Savic RM et al. J Pharmacokinet Pharmacodyn 2007;34:711-726
"""

using DifferentialEquations
using SciMLBase

export validate

function validate(spec::ModelSpec{TransitAbsorption,TransitAbsorptionParams})
    p = spec.params

    if p.N < 1
        error("Number of transit compartments (N) must be >= 1, got $(p.N)")
    end
    if p.N > 20
        error("Number of transit compartments (N) must be <= 20 for numerical stability, got $(p.N)")
    end

    _require_positive("Ktr", p.Ktr)
    _require_positive("Ka", p.Ka)
    _require_positive("CL", p.CL)
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
Transit compartment absorption ODE system.

States:
- u[1:N] = Transit compartments 1 through N
- u[N+1] = A_central (amount in central compartment)

The dose is added to the first transit compartment (u[1]).
Drug flows through the transit chain with rate Ktr,
then absorbs into central with rate Ka from the last transit compartment.

This creates a gamma-like absorption profile that can model
delayed peak times and variable absorption shapes.
"""
function _ode_transit_absorption!(du, u, p, t)
    N = p.N
    Ktr = p.Ktr
    Ka = p.Ka
    CL = p.CL
    V = p.V

    # First transit compartment (receives dose, loses to second transit)
    du[1] = -Ktr * u[1]

    # Middle transit compartments (N-1 of them, if N > 1)
    for i in 2:N
        du[i] = Ktr * u[i-1] - Ktr * u[i]
    end

    # Central compartment: absorption from last transit, elimination
    # Note: For N=1, Ka absorbs from u[1] (which is Transit[1])
    #       For N>1, Ka absorbs from u[N] (which is Transit[N])
    A_central = u[N+1]
    du[N+1] = Ka * u[N] - (CL / V) * A_central

    return nothing
end

"""
Generate state symbols for transit absorption model.
Returns [:Transit_1, :Transit_2, ..., :Transit_N, :A_central]
"""
function transit_state_symbols(N::Int)
    transit_syms = [Symbol("Transit_$i") for i in 1:N]
    return vcat(transit_syms, [:A_central])
end
