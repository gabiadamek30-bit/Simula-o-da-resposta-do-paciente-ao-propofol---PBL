# =============================================================================
# TMDD Simulation/Solve Functions
# =============================================================================
#
# Industry-standard TMDD simulation interface.
# Integrates with NeoPKPD's solver and callback infrastructure.
#
# Supports all TMDD model types:
# - OneCptTMDD, TwoCptTMDD (QSS, QE, Full, Rapid)
# - TwoCptTMDDFcRn (FcRn recycling)
# - TwoCptTMDDADA (immunogenicity)
# - TwoCptTMDDSoluble (soluble targets)
# - TwoCptTMDDMultiTarget (bispecifics)
# - TwoCptTMDDPeripheral (tissue targets)
# =============================================================================

using DifferentialEquations
using SciMLBase

export solve_tmdd, TMDDSimResult

# =============================================================================
# TMDD Simulation Result
# =============================================================================

"""
    TMDDSimResult

Result container for TMDD simulation.

# Fields
- `t`: Time points
- `states`: All state variables by name
- `observations`: Derived quantities (concentrations, occupancy, etc.)
- `metadata`: Additional information about the simulation
"""
struct TMDDSimResult
    t::Vector{Float64}
    states::Dict{Symbol,Vector{Float64}}
    observations::Dict{Symbol,Vector{Float64}}
    metadata::Dict{String,Any}
end

# =============================================================================
# Dosing Callbacks
# =============================================================================

"""
Create PresetTimeCallback for bolus dosing.
"""
function _create_bolus_callbacks(doses::Vector{DoseEvent}, dose_cpt::Int, F::Float64=1.0)
    bolus_doses = filter(is_bolus, doses)
    if isempty(bolus_doses)
        return nothing
    end

    dose_times = [d.time for d in bolus_doses]
    dose_amounts = [d.amount * F for d in bolus_doses]

    affect! = function(integrator)
        idx = findfirst(t -> isapprox(t, integrator.t, atol=1e-10), dose_times)
        if idx !== nothing
            integrator.u[dose_cpt] += dose_amounts[idx]
        end
    end

    return PresetTimeCallback(dose_times, affect!)
end

"""
Create callbacks for SC dosing with Tlag support.
"""
function _create_sc_callbacks(doses::Vector{DoseEvent}, dose_cpt::Int, Tlag::Float64, F::Float64)
    bolus_doses = filter(is_bolus, doses)
    if isempty(bolus_doses)
        return nothing
    end

    # Apply lag time to dose times
    dose_times = [d.time + Tlag for d in bolus_doses]
    dose_amounts = [d.amount * F for d in bolus_doses]

    affect! = function(integrator)
        idx = findfirst(t -> isapprox(t, integrator.t, atol=1e-10), dose_times)
        if idx !== nothing
            integrator.u[dose_cpt] += dose_amounts[idx]
        end
    end

    return PresetTimeCallback(dose_times, affect!)
end

# =============================================================================
# Main Solve Function
# =============================================================================

"""
    solve_tmdd(spec::TMDDSpec, grid::SimGrid, solver::SolverSpec) -> TMDDSimResult

Simulate a TMDD model.

# Arguments
- `spec::TMDDSpec`: TMDD model specification with parameters and dosing
- `grid::SimGrid`: Simulation time grid (t0, t1, saveat)
- `solver::SolverSpec`: ODE solver configuration

# Returns
- `TMDDSimResult`: Simulation results with states and derived observations

# Example
```julia
spec = TMDDSpec(
    TwoCptTMDD(QSS, IVBolus),
    "Pembrolizumab",
    TwoCptTMDDParams(
        CL=0.22, V1=3.28, V2=2.66, Q=0.54,
        KSS=0.087, kint=0.037, ksyn=0.001, kdeg=0.012, R0=0.083
    ),
    [DoseEvent(0.0, 200.0)]
)

grid = SimGrid(0.0, 672.0, collect(0.0:1.0:672.0))  # 28 days
solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

result = solve_tmdd(spec, grid, solver)
```
"""
function solve_tmdd(spec::TMDDSpec{K,P}, grid::SimGrid, solver::SolverSpec) where {K<:TMDDModelKind,P}
    # Validate
    validate_tmdd(spec)

    # Get model components
    ode_fn! = get_ode_function(spec.kind)
    u0 = get_initial_state(spec)
    dose_cpt = dosing_compartment(spec.kind)
    names = _get_state_names(spec)
    n_states_model = length(u0)

    # Handle SC with lag time
    is_sc = _is_subcutaneous(spec.kind)
    Tlag = _get_tlag(spec)
    F = _get_bioavailability(spec)

    # Apply first dose at t=0 if applicable (IV bolus only)
    first_dose_at_zero = false
    if !isempty(spec.doses) && spec.doses[1].time == 0.0 && is_bolus(spec.doses[1])
        if !is_sc && Tlag == 0.0
            # IV bolus - apply immediately
            first_dose_at_zero = true
            u0[dose_cpt] += spec.doses[1].amount
        end
    end

    # Create callbacks for remaining doses
    remaining_doses = first_dose_at_zero ? spec.doses[2:end] : spec.doses

    if is_sc
        callbacks = _create_sc_callbacks(remaining_doses, dose_cpt, Tlag, F)
    else
        callbacks = _create_bolus_callbacks(remaining_doses, dose_cpt)
    end

    # Get solver algorithm (prefer stiff solvers for TMDD)
    alg = _get_tmdd_algorithm(solver.alg)

    # Setup ODE problem
    tspan = (grid.t0, grid.t1)
    prob = ODEProblem(ode_fn!, u0, tspan, spec.params)

    # Solve with appropriate settings for stiff TMDD systems
    sol = solve(
        prob,
        alg;
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=callbacks,
        dtmax=_get_dtmax(spec),  # Limit step size for TMDD
    )

    # Check for successful solve
    if sol.retcode != :Success && sol.retcode != ReturnCode.Success
        @warn "TMDD solve may have issues: $(sol.retcode)"
    end

    # Extract results
    t = sol.t
    n_times = length(t)

    # Build states dictionary
    states = Dict{Symbol,Vector{Float64}}()
    for (i, name) in enumerate(names)
        if i <= n_states_model
            states[name] = [max(0.0, sol.u[j][i]) for j in 1:n_times]
        end
    end

    # Compute derived observations
    observations = _compute_observations(spec, states, t)

    # Metadata
    metadata = Dict{String,Any}(
        "model_type" => string(K),
        "model_name" => spec.name,
        "approximation" => string(_get_approximation(spec.kind)),
        "route" => string(_get_route(spec.kind)),
        "n_doses" => length(spec.doses),
        "solver" => string(solver.alg),
        "retcode" => string(sol.retcode),
        "target_units" => string(spec.target_units),
        "drug_units" => string(spec.drug_units)
    )

    return TMDDSimResult(t, states, observations, metadata)
end

# =============================================================================
# Solver Algorithm Selection
# =============================================================================

"""
Get DifferentialEquations algorithm optimized for TMDD.
TMDD systems are often stiff due to rapid binding kinetics.
"""
function _get_tmdd_algorithm(alg::Symbol)
    if alg == :Tsit5
        return Tsit5()
    elseif alg == :Rodas4
        return Rodas4()
    elseif alg == :Rodas5
        return Rodas5()
    elseif alg == :Rosenbrock23
        return Rosenbrock23()
    elseif alg == :TRBDF2
        return TRBDF2()
    elseif alg == :KenCarp4
        return KenCarp4()
    elseif alg == :Rodas4P
        return Rodas4P()
    elseif alg == :RadauIIA5
        return RadauIIA5()
    else
        # Default to robust stiff solver
        return Rodas5()
    end
end

"""
Get maximum time step for TMDD simulation.
"""
function _get_dtmax(spec::TMDDSpec)
    # Limit step size based on fastest process
    if hasproperty(spec.params, :kdeg)
        t_half_min = log(2) / spec.params.kdeg
        return min(t_half_min / 10.0, 1.0)  # At least 10 steps per half-life
    end
    return 1.0  # Default 1 day max step
end

# =============================================================================
# Helper Functions for Type Dispatch
# =============================================================================

_is_subcutaneous(kind::OneCptTMDD) = kind.route == Subcutaneous
_is_subcutaneous(kind::TwoCptTMDD) = kind.route == Subcutaneous
_is_subcutaneous(kind::TwoCptTMDDFcRn) = kind.route == Subcutaneous
_is_subcutaneous(kind::TwoCptTMDDADA) = kind.route == Subcutaneous
_is_subcutaneous(kind::TwoCptTMDDSoluble) = kind.route == Subcutaneous
_is_subcutaneous(kind::TwoCptTMDDMultiTarget) = kind.route == Subcutaneous
_is_subcutaneous(kind::TwoCptTMDDPeripheral) = kind.route == Subcutaneous

_get_approximation(kind::OneCptTMDD) = kind.approximation
_get_approximation(kind::TwoCptTMDD) = kind.approximation
_get_approximation(kind::TwoCptTMDDFcRn) = kind.approximation
_get_approximation(kind::TwoCptTMDDADA) = kind.approximation
_get_approximation(kind::TwoCptTMDDSoluble) = kind.approximation
_get_approximation(kind::TwoCptTMDDMultiTarget) = kind.approximation
_get_approximation(kind::TwoCptTMDDPeripheral) = kind.approximation

_get_route(kind::OneCptTMDD) = kind.route
_get_route(kind::TwoCptTMDD) = kind.route
_get_route(kind::TwoCptTMDDFcRn) = kind.route
_get_route(kind::TwoCptTMDDADA) = kind.route
_get_route(kind::TwoCptTMDDSoluble) = kind.route
_get_route(kind::TwoCptTMDDMultiTarget) = kind.route
_get_route(kind::TwoCptTMDDPeripheral) = kind.route

function _get_tlag(spec::TMDDSpec)
    if hasproperty(spec.params, :Tlag)
        return spec.params.Tlag
    end
    return 0.0
end

function _get_bioavailability(spec::TMDDSpec)
    if hasproperty(spec.params, :F)
        return spec.params.F
    end
    return 1.0
end

function _get_state_names(spec::TMDDSpec{OneCptTMDD,OneCptTMDDParams})
    return state_names(spec.kind)
end

function _get_state_names(spec::TMDDSpec{TwoCptTMDD,TwoCptTMDDParams})
    if spec.kind.route == Subcutaneous
        return vcat([:A_depot], state_names(spec.kind))
    end
    return state_names(spec.kind)
end

function _get_state_names(spec::TMDDSpec{TwoCptTMDDFcRn,TwoCptTMDDFcRnParams})
    return state_names(spec.kind)
end

function _get_state_names(spec::TMDDSpec{TwoCptTMDDADA,TwoCptTMDDADAParams})
    return state_names(spec.kind)
end

# Generic fallback
function _get_state_names(spec::TMDDSpec)
    if hasmethod(state_names, Tuple{typeof(spec.kind)})
        return state_names(spec.kind)
    end
    return [:u1, :u2, :u3]  # Fallback
end

# =============================================================================
# Observation Computation
# =============================================================================

"""
Compute derived observations from TMDD states.
Dispatches based on model type for correct calculations.
"""
function _compute_observations(spec::TMDDSpec{OneCptTMDD,OneCptTMDDParams}, states, t)
    p = spec.params
    V = p.V
    KSS = p.KSS

    observations = Dict{Symbol,Vector{Float64}}()

    if spec.kind.approximation == FullTMDD
        L = states[:L]
        R = states[:R]
        P = states[:P]

        observations[:conc] = L ./ V
        observations[:conc_free] = L ./ V
        observations[:conc_bound] = P ./ V
        observations[:conc_total] = (L .+ P) ./ V

        observations[:R_free] = R ./ V
        Rtot = R .+ P
        observations[:R_total] = Rtot ./ V
        observations[:target_occupancy] = [Rtot[i] > 0 ? P[i] / Rtot[i] : 0.0 for i in 1:length(t)]
    else
        # QSS model
        L = states[:L]
        Rtot = states[:Rtot]

        KSS_V = KSS * V
        RL = [L[i] * Rtot[i] / (KSS_V + L[i]) for i in 1:length(t)]

        observations[:conc] = L ./ V
        observations[:conc_free] = L ./ V
        observations[:conc_bound] = RL ./ V
        observations[:conc_total] = (L .+ RL) ./ V

        R_free = Rtot .- RL
        observations[:R_free] = R_free ./ V
        observations[:R_total] = Rtot ./ V
        observations[:target_occupancy] = [Rtot[i] > 0 ? RL[i] / Rtot[i] : 0.0 for i in 1:length(t)]
    end

    return observations
end

function _compute_observations(spec::TMDDSpec{TwoCptTMDD,TwoCptTMDDParams}, states, t)
    p = spec.params
    V1 = p.V1
    V2 = p.V2
    KSS = p.KSS

    observations = Dict{Symbol,Vector{Float64}}()

    # Handle SC depot
    is_sc = spec.kind.route == Subcutaneous
    L_key = is_sc ? :L : :L
    Lp_key = :Lp

    if spec.kind.approximation == FullTMDD
        L = haskey(states, :L) ? states[:L] : get(states, :L, zeros(length(t)))
        Lp = haskey(states, :Lp) ? states[:Lp] : zeros(length(t))
        R = haskey(states, :R) ? states[:R] : zeros(length(t))
        P = haskey(states, :P) ? states[:P] : zeros(length(t))

        observations[:conc] = L ./ V1
        observations[:conc_free] = L ./ V1
        observations[:conc_bound] = P ./ V1
        observations[:conc_total] = (L .+ P) ./ V1
        observations[:conc_peripheral] = Lp ./ V2

        observations[:R_free] = R ./ V1
        Rtot = R .+ P
        observations[:R_total] = Rtot ./ V1
        observations[:target_occupancy] = [Rtot[i] > 0 ? P[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    elseif spec.kind.approximation == RapidBinding
        Ltot = haskey(states, :Ltot) ? states[:Ltot] : zeros(length(t))
        Lp = haskey(states, :Lp) ? states[:Lp] : zeros(length(t))
        Rtot = haskey(states, :Rtot) ? states[:Rtot] : zeros(length(t))

        KSS_V1 = KSS * V1
        L_free = Vector{Float64}(undef, length(t))
        for i in 1:length(t)
            a = Ltot[i] - Rtot[i] - KSS_V1
            discriminant = a^2 + 4.0 * KSS_V1 * Ltot[i]
            L_free[i] = max(0.0, 0.5 * (a + sqrt(max(0.0, discriminant))))
        end
        P = Ltot .- L_free

        observations[:conc] = L_free ./ V1
        observations[:conc_free] = L_free ./ V1
        observations[:conc_bound] = P ./ V1
        observations[:conc_total] = Ltot ./ V1
        observations[:conc_peripheral] = Lp ./ V2

        R_free = Rtot .- P
        observations[:R_free] = R_free ./ V1
        observations[:R_total] = Rtot ./ V1
        observations[:target_occupancy] = [Rtot[i] > 0 ? P[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    else  # QSS or QE
        L = haskey(states, :L) ? states[:L] : zeros(length(t))
        Lp = haskey(states, :Lp) ? states[:Lp] : zeros(length(t))
        Rtot = haskey(states, :Rtot) ? states[:Rtot] : zeros(length(t))

        KSS_V1 = KSS * V1
        RL = [(KSS_V1 + L[i]) > 0 ? L[i] * Rtot[i] / (KSS_V1 + L[i]) : 0.0 for i in 1:length(t)]

        observations[:conc] = L ./ V1
        observations[:conc_free] = L ./ V1
        observations[:conc_bound] = RL ./ V1
        observations[:conc_total] = (L .+ RL) ./ V1
        observations[:conc_peripheral] = Lp ./ V2

        R_free = Rtot .- RL
        observations[:R_free] = R_free ./ V1
        observations[:R_total] = Rtot ./ V1
        observations[:target_occupancy] = [Rtot[i] > 0 ? RL[i] / Rtot[i] : 0.0 for i in 1:length(t)]
    end

    # Add SC depot amount if applicable
    if is_sc && haskey(states, :A_depot)
        observations[:depot_amount] = states[:A_depot]
    end

    return observations
end

function _compute_observations(spec::TMDDSpec{TwoCptTMDDFcRn,TwoCptTMDDFcRnParams}, states, t)
    p = spec.params
    V1 = p.V1
    V2 = p.V2
    KSS = p.KSS

    L = states[:L]
    Lp = states[:Lp]
    Le = states[:Le]
    Rtot = states[:Rtot]

    observations = Dict{Symbol,Vector{Float64}}()

    KSS_V1 = KSS * V1
    RL = [(KSS_V1 + L[i]) > 0 ? L[i] * Rtot[i] / (KSS_V1 + L[i]) : 0.0 for i in 1:length(t)]

    observations[:conc] = L ./ V1
    observations[:conc_free] = L ./ V1
    observations[:conc_bound] = RL ./ V1
    observations[:conc_total] = (L .+ RL) ./ V1
    observations[:conc_peripheral] = Lp ./ V2
    observations[:conc_endosomal] = Le ./ V1  # Relative to central volume

    R_free = Rtot .- RL
    observations[:R_free] = R_free ./ V1
    observations[:R_total] = Rtot ./ V1
    observations[:target_occupancy] = [Rtot[i] > 0 ? RL[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    # FcRn-specific observations
    observations[:endosomal_fraction] = Le ./ (L .+ Le .+ eps())

    return observations
end

function _compute_observations(spec::TMDDSpec{TwoCptTMDDADA,TwoCptTMDDADAParams}, states, t)
    p = spec.params
    V1 = p.V1
    V2 = p.V2
    KSS = p.KSS

    L = states[:L]
    Lp = states[:Lp]
    Rtot = states[:Rtot]
    ADA = states[:ADA]
    LADA = states[:LADA]

    observations = Dict{Symbol,Vector{Float64}}()

    KSS_V1 = KSS * V1
    RL = [(KSS_V1 + L[i]) > 0 ? L[i] * Rtot[i] / (KSS_V1 + L[i]) : 0.0 for i in 1:length(t)]

    observations[:conc] = L ./ V1
    observations[:conc_free] = L ./ V1
    observations[:conc_bound] = RL ./ V1
    observations[:conc_total] = (L .+ RL) ./ V1
    observations[:conc_peripheral] = Lp ./ V2

    R_free = Rtot .- RL
    observations[:R_free] = R_free ./ V1
    observations[:R_total] = Rtot ./ V1
    observations[:target_occupancy] = [Rtot[i] > 0 ? RL[i] / Rtot[i] : 0.0 for i in 1:length(t)]

    # ADA-specific observations
    observations[:ADA] = ADA ./ V1
    observations[:LADA] = LADA ./ V1  # Drug-ADA complex
    observations[:total_drug_with_ADA] = (L .+ RL .+ LADA) ./ V1

    # Neutralization ratio
    total_drug = L .+ RL .+ LADA
    observations[:neutralized_fraction] = [total_drug[i] > 0 ? LADA[i] / total_drug[i] : 0.0 for i in 1:length(t)]

    return observations
end

# Generic fallback for other model types
function _compute_observations(spec::TMDDSpec, states, t)
    observations = Dict{Symbol,Vector{Float64}}()

    # Try to compute basic observations
    if haskey(states, :L)
        V = _get_central_volume_from_spec(spec)
        L = states[:L]
        observations[:conc] = L ./ V
        observations[:conc_free] = L ./ V
    end

    if haskey(states, :Rtot)
        observations[:R_total] = states[:Rtot] ./ _get_central_volume_from_spec(spec)
    end

    return observations
end

function _get_central_volume_from_spec(spec::TMDDSpec)
    p = spec.params
    if hasproperty(p, :V1)
        return p.V1
    elseif hasproperty(p, :V)
        return p.V
    else
        return 1.0
    end
end
