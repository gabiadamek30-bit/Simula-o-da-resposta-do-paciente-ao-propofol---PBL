export pk_validate, pk_param_tuple, pk_state_symbols, pk_u0, pk_ode!, pk_conc
export pk_dose_target_index, pk_ode_with_infusion!

"""
Internal PK interface.

Each supported PK model kind must implement these methods:
- pk_validate(spec)
- pk_param_tuple(spec)
- pk_state_symbols(kind)
- pk_u0(spec, grid)
- pk_ode!(du, u, p, t, kind)
- pk_ode_with_infusion!(du, u, p, t, kind, infusion_rate) - ODE with zero-order input
- pk_conc(u, p, kind)
- pk_dose_target_index(kind)

The infusion_rate parameter represents zero-order drug input (amount/time)
that is added to the dose target compartment during IV infusion.
"""

# -------------------------
# Generic fallbacks
# -------------------------

function pk_validate(spec::ModelSpec)
    validate(spec)
end

# -------------------------
# OneCompIVBolus
# -------------------------

function pk_param_tuple(spec::ModelSpec{OneCompIVBolus,OneCompIVBolusParams})
    return (CL=spec.params.CL, V=spec.params.V)
end

pk_state_symbols(::OneCompIVBolus) = [:A_central]

pk_dose_target_index(::OneCompIVBolus) = 1

function pk_u0(spec::ModelSpec{OneCompIVBolus,OneCompIVBolusParams}, grid::SimGrid)
    return [0.0]
end

function pk_ode!(du, u, p, t, ::OneCompIVBolus)
    A = u[1]
    du[1] = -(p.CL / p.V) * A
    return nothing
end

function pk_conc(u, p, ::OneCompIVBolus)
    return u[1] / p.V
end

function pk_ode_with_infusion!(du, u, p, t, kind::OneCompIVBolus, infusion_rate::Float64)
    pk_ode!(du, u, p, t, kind)
    du[pk_dose_target_index(kind)] += infusion_rate
    return nothing
end

# -------------------------
# OneCompOralFirstOrder
# -------------------------

function pk_param_tuple(spec::ModelSpec{OneCompOralFirstOrder,OneCompOralFirstOrderParams})
    return (Ka=spec.params.Ka, CL=spec.params.CL, V=spec.params.V)
end

pk_state_symbols(::OneCompOralFirstOrder) = [:A_gut, :A_central]

pk_dose_target_index(::OneCompOralFirstOrder) = 1

function pk_u0(
    spec::ModelSpec{OneCompOralFirstOrder,OneCompOralFirstOrderParams}, grid::SimGrid
)
    return [0.0, 0.0]
end

function pk_ode!(du, u, p, t, ::OneCompOralFirstOrder)
    Agut = u[1]
    Acent = u[2]

    du[1] = -p.Ka * Agut
    du[2] = p.Ka * Agut - (p.CL / p.V) * Acent

    return nothing
end

function pk_conc(u, p, ::OneCompOralFirstOrder)
    return u[2] / p.V
end

function pk_ode_with_infusion!(du, u, p, t, kind::OneCompOralFirstOrder, infusion_rate::Float64)
    pk_ode!(du, u, p, t, kind)
    du[pk_dose_target_index(kind)] += infusion_rate
    return nothing
end

# -------------------------
# TwoCompIVBolus
# -------------------------

function pk_param_tuple(spec::ModelSpec{TwoCompIVBolus,TwoCompIVBolusParams})
    return (CL=spec.params.CL, V1=spec.params.V1, Q=spec.params.Q, V2=spec.params.V2)
end

pk_state_symbols(::TwoCompIVBolus) = [:A_central, :A_peripheral]

pk_dose_target_index(::TwoCompIVBolus) = 1

function pk_u0(spec::ModelSpec{TwoCompIVBolus,TwoCompIVBolusParams}, grid::SimGrid)
    return [0.0, 0.0]
end

function pk_ode!(du, u, p, t, ::TwoCompIVBolus)
    A_central = u[1]
    A_peripheral = u[2]

    k10 = p.CL / p.V1
    k12 = p.Q / p.V1
    k21 = p.Q / p.V2

    du[1] = -k10 * A_central - k12 * A_central + k21 * A_peripheral
    du[2] = k12 * A_central - k21 * A_peripheral

    return nothing
end

function pk_conc(u, p, ::TwoCompIVBolus)
    return u[1] / p.V1
end

function pk_ode_with_infusion!(du, u, p, t, kind::TwoCompIVBolus, infusion_rate::Float64)
    pk_ode!(du, u, p, t, kind)
    du[pk_dose_target_index(kind)] += infusion_rate
    return nothing
end

# -------------------------
# TwoCompOral
# -------------------------

function pk_param_tuple(spec::ModelSpec{TwoCompOral,TwoCompOralParams})
    return (Ka=spec.params.Ka, CL=spec.params.CL, V1=spec.params.V1, Q=spec.params.Q, V2=spec.params.V2)
end

pk_state_symbols(::TwoCompOral) = [:A_gut, :A_central, :A_peripheral]

pk_dose_target_index(::TwoCompOral) = 1

function pk_u0(spec::ModelSpec{TwoCompOral,TwoCompOralParams}, grid::SimGrid)
    return [0.0, 0.0, 0.0]
end

function pk_ode!(du, u, p, t, ::TwoCompOral)
    A_gut = u[1]
    A_central = u[2]
    A_peripheral = u[3]

    k10 = p.CL / p.V1
    k12 = p.Q / p.V1
    k21 = p.Q / p.V2

    du[1] = -p.Ka * A_gut
    du[2] = p.Ka * A_gut - k10 * A_central - k12 * A_central + k21 * A_peripheral
    du[3] = k12 * A_central - k21 * A_peripheral

    return nothing
end

function pk_conc(u, p, ::TwoCompOral)
    return u[2] / p.V1
end

function pk_ode_with_infusion!(du, u, p, t, kind::TwoCompOral, infusion_rate::Float64)
    pk_ode!(du, u, p, t, kind)
    du[pk_dose_target_index(kind)] += infusion_rate
    return nothing
end

# -------------------------
# ThreeCompIVBolus
# -------------------------

function pk_param_tuple(spec::ModelSpec{ThreeCompIVBolus,ThreeCompIVBolusParams})
    return (CL=spec.params.CL, V1=spec.params.V1, Q2=spec.params.Q2, V2=spec.params.V2, Q3=spec.params.Q3, V3=spec.params.V3)
end

pk_state_symbols(::ThreeCompIVBolus) = [:A_central, :A_periph1, :A_periph2]

pk_dose_target_index(::ThreeCompIVBolus) = 1

function pk_u0(spec::ModelSpec{ThreeCompIVBolus,ThreeCompIVBolusParams}, grid::SimGrid)
    return [0.0, 0.0, 0.0]
end

function pk_ode!(du, u, p, t, ::ThreeCompIVBolus)
    A_central = u[1]
    A_periph1 = u[2]
    A_periph2 = u[3]

    k10 = p.CL / p.V1
    k12 = p.Q2 / p.V1
    k21 = p.Q2 / p.V2
    k13 = p.Q3 / p.V1
    k31 = p.Q3 / p.V3

    du[1] = -k10 * A_central - k12 * A_central + k21 * A_periph1 - k13 * A_central + k31 * A_periph2
    du[2] = k12 * A_central - k21 * A_periph1
    du[3] = k13 * A_central - k31 * A_periph2

    return nothing
end

function pk_conc(u, p, ::ThreeCompIVBolus)
    return u[1] / p.V1
end

function pk_ode_with_infusion!(du, u, p, t, kind::ThreeCompIVBolus, infusion_rate::Float64)
    pk_ode!(du, u, p, t, kind)
    du[pk_dose_target_index(kind)] += infusion_rate
    return nothing
end

# -------------------------
# TransitAbsorption
# -------------------------

function pk_param_tuple(spec::ModelSpec{TransitAbsorption,TransitAbsorptionParams})
    return (N=spec.params.N, Ktr=spec.params.Ktr, Ka=spec.params.Ka, CL=spec.params.CL, V=spec.params.V)
end

function pk_state_symbols(::TransitAbsorption)
    # Generic symbols - actual number depends on N parameter
    return [:Transit_chain, :A_central]
end

pk_dose_target_index(::TransitAbsorption) = 1

function pk_u0(spec::ModelSpec{TransitAbsorption,TransitAbsorptionParams}, grid::SimGrid)
    N = spec.params.N
    # N transit compartments + 1 central compartment
    return zeros(N + 1)
end

function pk_ode!(du, u, p, t, ::TransitAbsorption)
    N = p.N
    Ktr = p.Ktr
    Ka = p.Ka
    CL = p.CL
    V = p.V

    # First transit compartment
    du[1] = -Ktr * u[1]

    # Middle transit compartments
    for i in 2:N
        du[i] = Ktr * u[i-1] - Ktr * u[i]
    end

    # Central compartment: absorption from last transit, elimination
    A_central = u[N+1]
    du[N+1] = Ka * u[N] - (CL / V) * A_central

    return nothing
end

function pk_conc(u, p, ::TransitAbsorption)
    N = p.N
    return u[N+1] / p.V
end

function pk_ode_with_infusion!(du, u, p, t, kind::TransitAbsorption, infusion_rate::Float64)
    pk_ode!(du, u, p, t, kind)
    du[pk_dose_target_index(kind)] += infusion_rate
    return nothing
end

# -------------------------
# MichaelisMentenElimination
# -------------------------

function pk_param_tuple(spec::ModelSpec{MichaelisMentenElimination,MichaelisMentenEliminationParams})
    return (Vmax=spec.params.Vmax, Km=spec.params.Km, V=spec.params.V)
end

pk_state_symbols(::MichaelisMentenElimination) = [:A_central]

pk_dose_target_index(::MichaelisMentenElimination) = 1

function pk_u0(spec::ModelSpec{MichaelisMentenElimination,MichaelisMentenEliminationParams}, grid::SimGrid)
    return [0.0]
end

function pk_ode!(du, u, p, t, ::MichaelisMentenElimination)
    Vmax = p.Vmax
    Km = p.Km
    V = p.V

    A_central = u[1]

    # Michaelis-Menten elimination
    # Rate = Vmax * C / (Km + C) = Vmax * A / (Km * V + A)
    du[1] = -Vmax * A_central / (Km * V + A_central)

    return nothing
end

function pk_conc(u, p, ::MichaelisMentenElimination)
    return u[1] / p.V
end

function pk_ode_with_infusion!(du, u, p, t, kind::MichaelisMentenElimination, infusion_rate::Float64)
    pk_ode!(du, u, p, t, kind)
    du[pk_dose_target_index(kind)] += infusion_rate
    return nothing
end
