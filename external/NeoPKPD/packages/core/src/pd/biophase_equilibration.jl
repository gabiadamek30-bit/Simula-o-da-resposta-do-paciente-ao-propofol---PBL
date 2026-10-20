"""
Biophase equilibration (effect compartment) PD model implementation.

This model introduces a hypothetical effect site compartment to account for
the temporal delay between plasma concentration changes and pharmacodynamic
effects. This is one of the most important PD models for understanding
PK-PD relationships.

Key concepts:
- The effect compartment is a hypothetical compartment with negligible volume
- It equilibrates with plasma concentration via first-order kinetics
- The rate constant ke0 determines the speed of equilibration
- Effect is calculated from effect site concentration, not plasma concentration

Clinical applications:
- Anesthetic agents (e.g., propofol, remifentanil)
- Neuromuscular blocking agents
- CNS-active drugs with delayed onset
- Any drug showing hysteresis in concentration-effect relationship

The equilibration half-life t1/2,ke0 = ln(2)/ke0 indicates:
- Small ke0 (long t1/2): Slow equilibration, significant hysteresis
- Large ke0 (short t1/2): Fast equilibration, near-direct effect

Reference: Sheiner LB et al. J Pharmacokinet Biopharm 1979;7:115-134
"""

using DifferentialEquations
using SciMLBase

export validate, evaluate, simulate_biophase

function validate(spec::PDSpec{BiophaseEquilibration,BiophaseEquilibrationParams})
    p = spec.params

    _require_positive("ke0", p.ke0)
    _require_positive("EC50", p.EC50)
    # E0 can be any real value (depending on baseline definition)
    # Emax can be negative for inhibitory effects

    if p.ke0 > 100.0
        @warn "ke0 > 100 may indicate very rapid equilibration; consider DirectEmax model"
    end

    return nothing
end

"""
Calculate effect site concentration at steady state for a given plasma concentration.
At steady state, Ce = Cp.
"""
function effect_site_concentration_ss(Cp::Float64)
    return Cp
end

"""
Calculate effect from effect site concentration using Emax model.

E(Ce) = E0 + (Emax * Ce) / (EC50 + Ce)
"""
function biophase_effect(Ce::Float64, p::BiophaseEquilibrationParams)
    if Ce <= 0.0
        return p.E0
    end
    return p.E0 + (p.Emax * Ce) / (p.EC50 + Ce)
end

"""
ODE for effect site compartment.

The effect compartment equilibrates with plasma concentration:
dCe/dt = ke0 * (Cp - Ce)

This is driven by the plasma concentration time course from the PK model.
"""
function _ode_biophase!(du, u, p, t)
    ke0 = p.ke0
    Cp_at_t = p.Cp_func(t)
    Ce = u[1]

    du[1] = ke0 * (Cp_at_t - Ce)

    return nothing
end

"""
Simulate biophase equilibration given a PK concentration time series.

Arguments:
- times: Time points from PK simulation
- concentrations: Plasma concentrations at those times
- params: BiophaseEquilibrationParams
- solver: SolverSpec

Returns:
- Ce: Effect site concentrations
- Effect: Pharmacodynamic effect

This function creates an interpolating function for Cp(t) and solves
the effect compartment ODE.
"""
function simulate_biophase(
    times::Vector{Float64},
    concentrations::Vector{Float64},
    params::BiophaseEquilibrationParams,
    solver::SolverSpec
)
    # Create interpolating function for plasma concentration
    # Use linear interpolation for efficiency and stability
    function Cp_func(t)
        if t <= times[1]
            return concentrations[1]
        elseif t >= times[end]
            return concentrations[end]
        else
            # Find bracketing indices
            idx = searchsortedlast(times, t)
            if idx == 0
                return concentrations[1]
            elseif idx >= length(times)
                return concentrations[end]
            else
                # Linear interpolation
                t0, t1 = times[idx], times[idx+1]
                C0, C1 = concentrations[idx], concentrations[idx+1]
                frac = (t - t0) / (t1 - t0)
                return C0 + frac * (C1 - C0)
            end
        end
    end

    # Initial condition: Ce = 0 (or could use Cp at t=0 for quasi-steady state)
    Ce0 = 0.0
    u0 = [Ce0]

    tspan = (times[1], times[end])

    # Parameters including the concentration function
    p = (ke0=params.ke0, Cp_func=Cp_func)

    prob = ODEProblem(_ode_biophase!, u0, tspan, p)

    sol = solve(
        prob,
        Tsit5();
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=times,
    )

    # Extract effect site concentrations
    Ce = [u[1] for u in sol.u]

    # Calculate effects
    effects = [biophase_effect(ce, params) for ce in Ce]

    return Ce, effects
end

"""
Evaluate the biophase equilibration model for a concentration series.

This is the quasi-steady state approximation, valid when ke0 >> elimination rate.
For rapid equilibration, Ce ≈ Cp, and this reduces to DirectEmax.

For true effect compartment dynamics with hysteresis, use simulate_biophase.
"""
function evaluate(
    spec::PDSpec{BiophaseEquilibration,BiophaseEquilibrationParams},
    input_series::Vector{Float64}
)
    validate(spec)

    p = spec.params

    # Assume quasi-steady state: Ce = Cp
    out = Vector{Float64}(undef, length(input_series))
    for i in eachindex(input_series)
        C = input_series[i]
        out[i] = biophase_effect(C, p)
    end
    return out
end

"""
Calculate the equilibration half-life.

t1/2,ke0 = ln(2) / ke0

This represents the time for the effect site concentration to reach
50% of the change toward a new plasma concentration steady state.
"""
function equilibration_half_life(ke0::Float64)
    return log(2) / ke0
end

"""
Calculate the time to 90% equilibration.

t90 = ln(10) / ke0 ≈ 2.3 / ke0

This is useful for understanding when effect will catch up to concentration.
"""
function time_to_90_percent_equilibration(ke0::Float64)
    return log(10) / ke0
end
