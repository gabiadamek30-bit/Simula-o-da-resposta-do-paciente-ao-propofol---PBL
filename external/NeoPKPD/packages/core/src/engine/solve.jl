using DifferentialEquations
using SciMLBase

export simulate

function validate(grid::SimGrid)
    if grid.t0 < 0.0
        error("t0 must be >= 0, got $(grid.t0)")
    end
    if !(grid.t1 > grid.t0)
        error("t1 must be > t0")
    end
    if isempty(grid.saveat)
        error("saveat must not be empty")
    end
    if any(t -> t < grid.t0 || t > grid.t1, grid.saveat)
        error("All saveat values must be within [t0, t1]")
    end
    if !issorted(grid.saveat)
        error("saveat must be sorted ascending")
    end
    return nothing
end

function validate(solver::SolverSpec)
    if !(solver.reltol > 0.0)
        error("Expected positive value for reltol, got $(solver.reltol)")
    end
    if !(solver.abstol > 0.0)
        error("Expected positive value for abstol, got $(solver.abstol)")
    end
    if solver.maxiters < 1
        error("maxiters must be >= 1")
    end
    return nothing
end

const _SOLVER_MAP = Dict{Symbol,Any}(:Tsit5 => Tsit5, :Rosenbrock23 => Rosenbrock23)

function _solver_alg(alg::Symbol)
    if !haskey(_SOLVER_MAP, alg)
        error("Unsupported solver alg: $(alg). Supported: $(collect(keys(_SOLVER_MAP)))")
    end
    return _SOLVER_MAP[alg]()
end

"""
Check if any doses in the vector are infusions (duration > 0).
"""
function _has_any_infusion(doses::Vector{DoseEvent})::Bool
    for d in doses
        if is_infusion(d)
            return true
        end
    end
    return false
end

function _dose_callback(doses::Vector{DoseEvent}, t0::Float64, t1::Float64)
    _, dose_times, dose_amounts = normalize_doses_for_sim(doses, t0, t1)

    if isempty(dose_times)
        return nothing
    end

    function affect!(integrator)
        idx = findfirst(==(integrator.t), dose_times)
        if idx === nothing
            error("Internal error: dose time not found for t=$(integrator.t)")
        end
        integrator.u[1] += dose_amounts[idx]
    end

    return PresetTimeCallback(dose_times, affect!)
end

function simulate(
    spec::ModelSpec{OneCompIVBolus,OneCompIVBolusParams}, grid::SimGrid, solver::SolverSpec
)
    validate(spec)
    validate(grid)
    validate(solver)

    CL = spec.params.CL
    V = spec.params.V
    p = (CL=CL, V=V)
    tspan = (grid.t0, grid.t1)

    # Check if we have any infusion doses
    if _has_any_infusion(spec.doses)
        # Use infusion-aware simulation
        schedule = build_infusion_schedule(spec.doses, grid.t0, grid.t1)
        u0 = [schedule.initial_amount]

        # Create ODE that includes infusion rate
        rate_fn = make_infusion_rate_function(schedule)
        function ode_with_infusion!(du, u, params, t)
            pk_ode_with_infusion!(du, u, params, t, OneCompIVBolus(), rate_fn(t))
        end

        prob = ODEProblem(ode_with_infusion!, u0, tspan, p)
        cb = build_bolus_callback(schedule, pk_dose_target_index(OneCompIVBolus()))
        tstops = get_infusion_tstops(schedule, grid.t0, grid.t1)

        sol = solve(
            prob,
            _solver_alg(solver.alg);
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            callback=cb,
            tstops=tstops,
        )
    else
        # Original bolus-only path
        a0_add, _, _ = normalize_doses_for_sim(spec.doses, grid.t0, grid.t1)
        A0 = a0_add
        u0 = [A0]

        prob = ODEProblem(_ode_onecomp_ivbolus!, u0, tspan, p)
        cb = _dose_callback(spec.doses, grid.t0, grid.t1)

        sol = solve(
            prob,
            _solver_alg(solver.alg);
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            callback=cb,
        )
    end

    A = [u[1] for u in sol.u]
    C = [a / V for a in A]

    states = Dict(:A_central => A)
    observations = Dict(:conc => C)

    # Include infusion info in metadata if present
    dose_schedule = [(d.time, d.amount, d.duration) for d in spec.doses]

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "OneCompIVBolus",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => dose_schedule,
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

function simulate(
    spec::ModelSpec{OneCompOralFirstOrder,OneCompOralFirstOrderParams},
    grid::SimGrid,
    solver::SolverSpec,
)
    validate(spec)
    validate(grid)
    validate(solver)

    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V

    # Oral bolus doses add to Agut
    a0_add, _, _ = normalize_doses_for_sim(spec.doses, grid.t0, grid.t1)
    Agut0 = a0_add

    p = (Ka=Ka, CL=CL, V=V)
    u0 = [Agut0, 0.0]
    tspan = (grid.t0, grid.t1)

    prob = ODEProblem(_ode_onecomp_oral_first_order!, u0, tspan, p)

    # Dose callback adds to gut compartment
    dose_times = Float64[]
    dose_amounts = Float64[]
    for d in spec.doses
        if d.time > grid.t0 && d.time <= grid.t1
            push!(dose_times, d.time)
            push!(dose_amounts, d.amount)
        end
    end

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx === nothing
                error("Internal error: dose time not found for t=$(integrator.t)")
            end
            integrator.u[1] += dose_amounts[idx]
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    Agut = [u[1] for u in sol.u]
    Acent = [u[2] for u in sol.u]
    C = [a / V for a in Acent]

    states = Dict(:A_gut => Agut, :A_central => Acent)

    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "OneCompOralFirstOrder",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount) for d in spec.doses],
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

# -------------------------
# TwoCompIVBolus
# -------------------------

function simulate(
    spec::ModelSpec{TwoCompIVBolus,TwoCompIVBolusParams}, grid::SimGrid, solver::SolverSpec
)
    validate(spec)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)
    tspan = (grid.t0, grid.t1)

    # Check if we have any infusion doses
    if _has_any_infusion(spec.doses)
        # Use infusion-aware simulation
        schedule = build_infusion_schedule(spec.doses, grid.t0, grid.t1)
        u0 = [schedule.initial_amount, 0.0]

        # Create ODE that includes infusion rate
        rate_fn = make_infusion_rate_function(schedule)
        function ode_with_infusion!(du, u, params, t)
            pk_ode_with_infusion!(du, u, params, t, TwoCompIVBolus(), rate_fn(t))
        end

        prob = ODEProblem(ode_with_infusion!, u0, tspan, p)
        cb = build_bolus_callback(schedule, pk_dose_target_index(TwoCompIVBolus()))
        tstops = get_infusion_tstops(schedule, grid.t0, grid.t1)

        sol = solve(
            prob,
            _solver_alg(solver.alg);
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            callback=cb,
            tstops=tstops,
        )
    else
        # Original bolus-only path
        a0_add, dose_times, dose_amounts = normalize_doses_for_sim(spec.doses, grid.t0, grid.t1)
        u0 = [a0_add, 0.0]

        function ode!(du, u, params, t)
            pk_ode!(du, u, params, t, TwoCompIVBolus())
        end

        prob = ODEProblem(ode!, u0, tspan, p)

        cb = nothing
        if !isempty(dose_times)
            function affect!(integrator)
                idx = findfirst(==(integrator.t), dose_times)
                if idx !== nothing
                    integrator.u[1] += dose_amounts[idx]
                end
            end
            cb = PresetTimeCallback(dose_times, affect!)
        end

        sol = solve(
            prob,
            _solver_alg(solver.alg);
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            callback=cb,
        )
    end

    A_central = [u[1] for u in sol.u]
    A_peripheral = [u[2] for u in sol.u]
    C = [a / p.V1 for a in A_central]

    states = Dict(:A_central => A_central, :A_peripheral => A_peripheral)
    observations = Dict(:conc => C)

    dose_schedule = [(d.time, d.amount, d.duration) for d in spec.doses]

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "TwoCompIVBolus",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => dose_schedule,
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

# -------------------------
# TwoCompOral
# -------------------------

function simulate(
    spec::ModelSpec{TwoCompOral,TwoCompOralParams}, grid::SimGrid, solver::SolverSpec
)
    validate(spec)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)

    a0_add, dose_times, dose_amounts = normalize_doses_for_sim(spec.doses, grid.t0, grid.t1)

    u0 = [a0_add, 0.0, 0.0]  # [A_gut, A_central, A_peripheral]
    tspan = (grid.t0, grid.t1)

    function ode!(du, u, params, t)
        pk_ode!(du, u, params, t, TwoCompOral())
    end

    prob = ODEProblem(ode!, u0, tspan, p)

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx !== nothing
                integrator.u[1] += dose_amounts[idx]  # Dose goes to gut compartment
            end
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    A_gut = [u[1] for u in sol.u]
    A_central = [u[2] for u in sol.u]
    A_peripheral = [u[3] for u in sol.u]
    C = [a / p.V1 for a in A_central]

    states = Dict(:A_gut => A_gut, :A_central => A_central, :A_peripheral => A_peripheral)
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "TwoCompOral",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount) for d in spec.doses],
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

# -------------------------
# ThreeCompIVBolus
# -------------------------

function simulate(
    spec::ModelSpec{ThreeCompIVBolus,ThreeCompIVBolusParams}, grid::SimGrid, solver::SolverSpec
)
    validate(spec)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)
    tspan = (grid.t0, grid.t1)

    # Check if we have any infusion doses
    if _has_any_infusion(spec.doses)
        # Use infusion-aware simulation
        schedule = build_infusion_schedule(spec.doses, grid.t0, grid.t1)
        u0 = [schedule.initial_amount, 0.0, 0.0]

        # Create ODE that includes infusion rate
        rate_fn = make_infusion_rate_function(schedule)
        function ode_with_infusion!(du, u, params, t)
            pk_ode_with_infusion!(du, u, params, t, ThreeCompIVBolus(), rate_fn(t))
        end

        prob = ODEProblem(ode_with_infusion!, u0, tspan, p)
        cb = build_bolus_callback(schedule, pk_dose_target_index(ThreeCompIVBolus()))
        tstops = get_infusion_tstops(schedule, grid.t0, grid.t1)

        sol = solve(
            prob,
            _solver_alg(solver.alg);
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            callback=cb,
            tstops=tstops,
        )
    else
        # Original bolus-only path
        a0_add, dose_times, dose_amounts = normalize_doses_for_sim(spec.doses, grid.t0, grid.t1)
        u0 = [a0_add, 0.0, 0.0]  # [A_central, A_periph1, A_periph2]

        function ode!(du, u, params, t)
            pk_ode!(du, u, params, t, ThreeCompIVBolus())
        end

        prob = ODEProblem(ode!, u0, tspan, p)

        cb = nothing
        if !isempty(dose_times)
            function affect!(integrator)
                idx = findfirst(==(integrator.t), dose_times)
                if idx !== nothing
                    integrator.u[1] += dose_amounts[idx]
                end
            end
            cb = PresetTimeCallback(dose_times, affect!)
        end

        sol = solve(
            prob,
            _solver_alg(solver.alg);
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            callback=cb,
        )
    end

    A_central = [u[1] for u in sol.u]
    A_periph1 = [u[2] for u in sol.u]
    A_periph2 = [u[3] for u in sol.u]
    C = [a / p.V1 for a in A_central]

    states = Dict(:A_central => A_central, :A_periph1 => A_periph1, :A_periph2 => A_periph2)
    observations = Dict(:conc => C)

    dose_schedule = [(d.time, d.amount, d.duration) for d in spec.doses]

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "ThreeCompIVBolus",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => dose_schedule,
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

# -------------------------
# TransitAbsorption
# -------------------------

function simulate(
    spec::ModelSpec{TransitAbsorption,TransitAbsorptionParams}, grid::SimGrid, solver::SolverSpec
)
    validate(spec)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)
    N = spec.params.N

    a0_add, dose_times, dose_amounts = normalize_doses_for_sim(spec.doses, grid.t0, grid.t1)

    # N transit compartments + 1 central compartment
    u0 = zeros(N + 1)
    u0[1] = a0_add  # Initial dose goes to first transit compartment

    tspan = (grid.t0, grid.t1)

    function ode!(du, u, params, t)
        pk_ode!(du, u, params, t, TransitAbsorption())
    end

    prob = ODEProblem(ode!, u0, tspan, p)

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx !== nothing
                integrator.u[1] += dose_amounts[idx]  # Dose goes to first transit
            end
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    # Extract all transit compartments and central
    transit_states = Dict{Symbol,Vector{Float64}}()
    for i in 1:N
        transit_states[Symbol("Transit_$i")] = [u[i] for u in sol.u]
    end

    A_central = [u[N+1] for u in sol.u]
    C = [a / p.V for a in A_central]

    states = merge(transit_states, Dict(:A_central => A_central))
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "TransitAbsorption",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount) for d in spec.doses],
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
        "N_transit" => N,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

# -------------------------
# MichaelisMentenElimination
# -------------------------

function simulate(
    spec::ModelSpec{MichaelisMentenElimination,MichaelisMentenEliminationParams}, grid::SimGrid, solver::SolverSpec
)
    validate(spec)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)
    tspan = (grid.t0, grid.t1)

    # Use Rosenbrock23 for stiff nonlinear elimination by default
    alg = solver.alg == :Tsit5 ? Rosenbrock23() : _solver_alg(solver.alg)

    # Check if we have any infusion doses
    if _has_any_infusion(spec.doses)
        # Use infusion-aware simulation
        schedule = build_infusion_schedule(spec.doses, grid.t0, grid.t1)
        u0 = [schedule.initial_amount]

        # Create ODE that includes infusion rate
        rate_fn = make_infusion_rate_function(schedule)
        function ode_with_infusion!(du, u, params, t)
            pk_ode_with_infusion!(du, u, params, t, MichaelisMentenElimination(), rate_fn(t))
        end

        prob = ODEProblem(ode_with_infusion!, u0, tspan, p)
        cb = build_bolus_callback(schedule, pk_dose_target_index(MichaelisMentenElimination()))
        tstops = get_infusion_tstops(schedule, grid.t0, grid.t1)

        sol = solve(
            prob,
            alg;
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            callback=cb,
            tstops=tstops,
        )
    else
        # Original bolus-only path
        a0_add, dose_times, dose_amounts = normalize_doses_for_sim(spec.doses, grid.t0, grid.t1)
        u0 = [a0_add]

        function ode!(du, u, params, t)
            pk_ode!(du, u, params, t, MichaelisMentenElimination())
        end

        prob = ODEProblem(ode!, u0, tspan, p)

        cb = nothing
        if !isempty(dose_times)
            function affect!(integrator)
                idx = findfirst(==(integrator.t), dose_times)
                if idx !== nothing
                    integrator.u[1] += dose_amounts[idx]
                end
            end
            cb = PresetTimeCallback(dose_times, affect!)
        end

        sol = solve(
            prob,
            alg;
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            callback=cb,
        )
    end

    A = [u[1] for u in sol.u]
    C = [a / p.V for a in A]

    states = Dict(:A_central => A)
    observations = Dict(:conc => C)

    dose_schedule = [(d.time, d.amount, d.duration) for d in spec.doses]

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "MichaelisMentenElimination",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => dose_schedule,
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

# ============================================================================
# ModelSpecWithModifiers Support (ALAG and Bioavailability)
# ============================================================================

"""
Apply absorption modifiers (ALAG and F) to doses for simulation.

Returns:
- initial_amount: Amount for initial condition (doses at or before effective t0)
- dose_times: Times of callback doses (after ALAG applied)
- dose_amounts: Amounts of callback doses (after F applied)
"""
function _normalize_doses_with_absorption_modifiers(
    doses::Vector{DoseEvent},
    modifiers::AbsorptionModifiers,
    t0::Float64,
    t1::Float64
)
    initial_amount = 0.0
    dose_times_dict = Dict{Float64, Float64}()

    alag = modifiers.alag
    F = modifiers.bioavailability

    for dose in doses
        # Apply modifiers
        effective_time = dose.time + alag
        effective_amount = dose.amount * F

        if is_bolus(dose)
            if effective_time < t0
                # Dose completed before t0, add to initial
                initial_amount += effective_amount
            elseif effective_time == t0
                # Dose at exactly t0
                initial_amount += effective_amount
            elseif effective_time <= t1
                # Dose within (t0, t1]
                dose_times_dict[effective_time] = get(dose_times_dict, effective_time, 0.0) + effective_amount
            end
        else
            # Infusion handling with ALAG/F
            effective_end = effective_time + dose.duration

            if effective_end <= t0
                # Infusion completed before t0
                initial_amount += effective_amount
            elseif effective_time >= t1
                # Infusion starts after simulation ends - ignore
            else
                # Infusion overlaps with simulation
                dose_times_dict[effective_time] = get(dose_times_dict, effective_time, 0.0) + effective_amount
            end
        end
    end

    dose_times = sort(collect(keys(dose_times_dict)))
    dose_amounts = [dose_times_dict[t] for t in dose_times]

    return initial_amount, dose_times, dose_amounts
end

# -------------------------
# OneCompOralFirstOrder with Modifiers
# -------------------------

function simulate(
    spec::ModelSpecWithModifiers{OneCompOralFirstOrder,OneCompOralFirstOrderParams},
    grid::SimGrid,
    solver::SolverSpec,
)
    # Validate inputs
    validate(spec.params)
    validate(grid)
    validate(solver)

    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V
    modifiers = spec.absorption_modifiers

    # Apply ALAG and bioavailability to doses
    a0_add, dose_times, dose_amounts = _normalize_doses_with_absorption_modifiers(
        spec.doses, modifiers, grid.t0, grid.t1
    )

    p = (Ka=Ka, CL=CL, V=V)
    u0 = [a0_add, 0.0]  # [A_gut, A_central]
    tspan = (grid.t0, grid.t1)

    prob = ODEProblem(_ode_onecomp_oral_first_order!, u0, tspan, p)

    # Build callback for ALAG-shifted doses
    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx !== nothing
                integrator.u[1] += dose_amounts[idx]  # Add to gut
            end
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    # Include ALAG-shifted times in tstops for accuracy
    tstops = filter(t -> t > grid.t0 && t < grid.t1, dose_times)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
        tstops=tstops,
    )

    Agut = [u[1] for u in sol.u]
    Acent = [u[2] for u in sol.u]
    C = [a / V for a in Acent]

    states = Dict(:A_gut => Agut, :A_central => Acent)
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "OneCompOralFirstOrder",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount) for d in spec.doses],
        "absorption_modifiers" => Dict(
            "ALAG" => modifiers.alag,
            "F" => modifiers.bioavailability
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

# Validate method for params (reuse existing)
function validate(params::OneCompOralFirstOrderParams)
    _require_positive("Ka", params.Ka)
    _require_positive("CL", params.CL)
    _require_positive("V", params.V)
    return nothing
end

# -------------------------
# TwoCompOral with Modifiers
# -------------------------

function simulate(
    spec::ModelSpecWithModifiers{TwoCompOral,TwoCompOralParams},
    grid::SimGrid,
    solver::SolverSpec
)
    validate(spec.params)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)
    modifiers = spec.absorption_modifiers

    # Apply ALAG and bioavailability to doses
    a0_add, dose_times, dose_amounts = _normalize_doses_with_absorption_modifiers(
        spec.doses, modifiers, grid.t0, grid.t1
    )

    u0 = [a0_add, 0.0, 0.0]  # [A_gut, A_central, A_peripheral]
    tspan = (grid.t0, grid.t1)

    function ode!(du, u, params, t)
        pk_ode!(du, u, params, t, TwoCompOral())
    end

    prob = ODEProblem(ode!, u0, tspan, p)

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx !== nothing
                integrator.u[1] += dose_amounts[idx]  # Add to gut
            end
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    tstops = filter(t -> t > grid.t0 && t < grid.t1, dose_times)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
        tstops=tstops,
    )

    A_gut = [u[1] for u in sol.u]
    A_central = [u[2] for u in sol.u]
    A_peripheral = [u[3] for u in sol.u]
    C = [a / p.V1 for a in A_central]

    states = Dict(:A_gut => A_gut, :A_central => A_central, :A_peripheral => A_peripheral)
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "TwoCompOral",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount) for d in spec.doses],
        "absorption_modifiers" => Dict(
            "ALAG" => modifiers.alag,
            "F" => modifiers.bioavailability
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

function validate(params::TwoCompOralParams)
    _require_positive("Ka", params.Ka)
    _require_positive("CL", params.CL)
    _require_positive("V1", params.V1)
    _require_positive("Q", params.Q)
    _require_positive("V2", params.V2)
    return nothing
end

# pk_param_tuple for ModelSpecWithModifiers
function pk_param_tuple(spec::ModelSpecWithModifiers{TwoCompOral,TwoCompOralParams})
    return (Ka=spec.params.Ka, CL=spec.params.CL, V1=spec.params.V1, Q=spec.params.Q, V2=spec.params.V2)
end

# -------------------------
# TransitAbsorption with Modifiers
# -------------------------

function simulate(
    spec::ModelSpecWithModifiers{TransitAbsorption,TransitAbsorptionParams},
    grid::SimGrid,
    solver::SolverSpec
)
    validate(spec.params)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)
    N = spec.params.N
    modifiers = spec.absorption_modifiers

    # Apply ALAG and bioavailability to doses
    a0_add, dose_times, dose_amounts = _normalize_doses_with_absorption_modifiers(
        spec.doses, modifiers, grid.t0, grid.t1
    )

    # N transit compartments + 1 central compartment
    u0 = zeros(N + 1)
    u0[1] = a0_add  # Initial dose goes to first transit compartment

    tspan = (grid.t0, grid.t1)

    function ode!(du, u, params, t)
        pk_ode!(du, u, params, t, TransitAbsorption())
    end

    prob = ODEProblem(ode!, u0, tspan, p)

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx !== nothing
                integrator.u[1] += dose_amounts[idx]  # Dose goes to first transit
            end
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    tstops = filter(t -> t > grid.t0 && t < grid.t1, dose_times)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
        tstops=tstops,
    )

    # Extract all transit compartments and central
    transit_states = Dict{Symbol,Vector{Float64}}()
    for i in 1:N
        transit_states[Symbol("Transit_$i")] = [u[i] for u in sol.u]
    end

    A_central = [u[N+1] for u in sol.u]
    C = [a / p.V for a in A_central]

    states = merge(transit_states, Dict(:A_central => A_central))
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "TransitAbsorption",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount) for d in spec.doses],
        "absorption_modifiers" => Dict(
            "ALAG" => modifiers.alag,
            "F" => modifiers.bioavailability
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
        "N_transit" => N,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

function validate(params::TransitAbsorptionParams)
    if params.N < 1
        error("N (number of transit compartments) must be >= 1, got $(params.N)")
    end
    _require_positive("Ktr", params.Ktr)
    _require_positive("Ka", params.Ka)
    _require_positive("CL", params.CL)
    _require_positive("V", params.V)
    return nothing
end

function pk_param_tuple(spec::ModelSpecWithModifiers{TransitAbsorption,TransitAbsorptionParams})
    return (N=spec.params.N, Ktr=spec.params.Ktr, Ka=spec.params.Ka, CL=spec.params.CL, V=spec.params.V)
end

# -------------------------
# IV Models with Modifiers (ALAG only, F=1 typical)
# For completeness, allow modifiers on IV models too
# -------------------------

function simulate(
    spec::ModelSpecWithModifiers{OneCompIVBolus,OneCompIVBolusParams},
    grid::SimGrid,
    solver::SolverSpec
)
    validate(spec.params)
    validate(grid)
    validate(solver)

    CL = spec.params.CL
    V = spec.params.V
    p = (CL=CL, V=V)
    modifiers = spec.absorption_modifiers

    # Apply modifiers to doses
    a0_add, dose_times, dose_amounts = _normalize_doses_with_absorption_modifiers(
        spec.doses, modifiers, grid.t0, grid.t1
    )

    tspan = (grid.t0, grid.t1)
    u0 = [a0_add]

    prob = ODEProblem(_ode_onecomp_ivbolus!, u0, tspan, p)

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx !== nothing
                integrator.u[1] += dose_amounts[idx]
            end
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    tstops = filter(t -> t > grid.t0 && t < grid.t1, dose_times)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
        tstops=tstops,
    )

    A = [u[1] for u in sol.u]
    C = [a / V for a in A]

    states = Dict(:A_central => A)
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "OneCompIVBolus",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount, d.duration) for d in spec.doses],
        "absorption_modifiers" => Dict(
            "ALAG" => modifiers.alag,
            "F" => modifiers.bioavailability
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

function validate(params::OneCompIVBolusParams)
    _require_positive("CL", params.CL)
    _require_positive("V", params.V)
    return nothing
end

function simulate(
    spec::ModelSpecWithModifiers{TwoCompIVBolus,TwoCompIVBolusParams},
    grid::SimGrid,
    solver::SolverSpec
)
    validate(spec.params)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)
    modifiers = spec.absorption_modifiers

    a0_add, dose_times, dose_amounts = _normalize_doses_with_absorption_modifiers(
        spec.doses, modifiers, grid.t0, grid.t1
    )

    tspan = (grid.t0, grid.t1)
    u0 = [a0_add, 0.0]

    function ode!(du, u, params, t)
        pk_ode!(du, u, params, t, TwoCompIVBolus())
    end

    prob = ODEProblem(ode!, u0, tspan, p)

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx !== nothing
                integrator.u[1] += dose_amounts[idx]
            end
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    tstops = filter(t -> t > grid.t0 && t < grid.t1, dose_times)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
        tstops=tstops,
    )

    A_central = [u[1] for u in sol.u]
    A_peripheral = [u[2] for u in sol.u]
    C = [a / p.V1 for a in A_central]

    states = Dict(:A_central => A_central, :A_peripheral => A_peripheral)
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "TwoCompIVBolus",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount, d.duration) for d in spec.doses],
        "absorption_modifiers" => Dict(
            "ALAG" => modifiers.alag,
            "F" => modifiers.bioavailability
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

function validate(params::TwoCompIVBolusParams)
    _require_positive("CL", params.CL)
    _require_positive("V1", params.V1)
    _require_positive("Q", params.Q)
    _require_positive("V2", params.V2)
    return nothing
end

function pk_param_tuple(spec::ModelSpecWithModifiers{TwoCompIVBolus,TwoCompIVBolusParams})
    return (CL=spec.params.CL, V1=spec.params.V1, Q=spec.params.Q, V2=spec.params.V2)
end

function simulate(
    spec::ModelSpecWithModifiers{ThreeCompIVBolus,ThreeCompIVBolusParams},
    grid::SimGrid,
    solver::SolverSpec
)
    validate(spec.params)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)
    modifiers = spec.absorption_modifiers

    a0_add, dose_times, dose_amounts = _normalize_doses_with_absorption_modifiers(
        spec.doses, modifiers, grid.t0, grid.t1
    )

    tspan = (grid.t0, grid.t1)
    u0 = [a0_add, 0.0, 0.0]

    function ode!(du, u, params, t)
        pk_ode!(du, u, params, t, ThreeCompIVBolus())
    end

    prob = ODEProblem(ode!, u0, tspan, p)

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx !== nothing
                integrator.u[1] += dose_amounts[idx]
            end
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    tstops = filter(t -> t > grid.t0 && t < grid.t1, dose_times)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
        tstops=tstops,
    )

    A_central = [u[1] for u in sol.u]
    A_periph1 = [u[2] for u in sol.u]
    A_periph2 = [u[3] for u in sol.u]
    C = [a / p.V1 for a in A_central]

    states = Dict(:A_central => A_central, :A_periph1 => A_periph1, :A_periph2 => A_periph2)
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "ThreeCompIVBolus",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount, d.duration) for d in spec.doses],
        "absorption_modifiers" => Dict(
            "ALAG" => modifiers.alag,
            "F" => modifiers.bioavailability
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

function validate(params::ThreeCompIVBolusParams)
    _require_positive("CL", params.CL)
    _require_positive("V1", params.V1)
    _require_positive("Q2", params.Q2)
    _require_positive("V2", params.V2)
    _require_positive("Q3", params.Q3)
    _require_positive("V3", params.V3)
    return nothing
end

function pk_param_tuple(spec::ModelSpecWithModifiers{ThreeCompIVBolus,ThreeCompIVBolusParams})
    return (CL=spec.params.CL, V1=spec.params.V1, Q2=spec.params.Q2, V2=spec.params.V2, Q3=spec.params.Q3, V3=spec.params.V3)
end

function simulate(
    spec::ModelSpecWithModifiers{MichaelisMentenElimination,MichaelisMentenEliminationParams},
    grid::SimGrid,
    solver::SolverSpec
)
    validate(spec.params)
    validate(grid)
    validate(solver)

    p = pk_param_tuple(spec)
    modifiers = spec.absorption_modifiers

    a0_add, dose_times, dose_amounts = _normalize_doses_with_absorption_modifiers(
        spec.doses, modifiers, grid.t0, grid.t1
    )

    tspan = (grid.t0, grid.t1)
    u0 = [a0_add]

    # Use Rosenbrock23 for stiff nonlinear elimination by default
    alg = solver.alg == :Tsit5 ? Rosenbrock23() : _solver_alg(solver.alg)

    function ode!(du, u, params, t)
        pk_ode!(du, u, params, t, MichaelisMentenElimination())
    end

    prob = ODEProblem(ode!, u0, tspan, p)

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx !== nothing
                integrator.u[1] += dose_amounts[idx]
            end
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    tstops = filter(t -> t > grid.t0 && t < grid.t1, dose_times)

    sol = solve(
        prob,
        alg;
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
        tstops=tstops,
    )

    A = [u[1] for u in sol.u]
    C = [a / p.V for a in A]

    states = Dict(:A_central => A)
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "MichaelisMentenElimination",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount, d.duration) for d in spec.doses],
        "absorption_modifiers" => Dict(
            "ALAG" => modifiers.alag,
            "F" => modifiers.bioavailability
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

function validate(params::MichaelisMentenEliminationParams)
    _require_positive("Vmax", params.Vmax)
    _require_positive("Km", params.Km)
    _require_positive("V", params.V)
    return nothing
end

function pk_param_tuple(spec::ModelSpecWithModifiers{MichaelisMentenElimination,MichaelisMentenEliminationParams})
    return (Vmax=spec.params.Vmax, Km=spec.params.Km, V=spec.params.V)
end

# ============================================================================
# CustomODE Simulation
# ============================================================================

"""
    simulate(spec::ModelSpec{<:CustomODE, CustomODEParams}, grid::SimGrid, solver::SolverSpec)

Simulate a custom ODE model.

This is a generic simulation method that uses the CustomODE interface functions:
- pk_param_tuple: Get parameters as NamedTuple
- pk_u0: Get initial conditions
- pk_ode!: ODE function
- pk_conc: Compute concentration from state
- pk_dose_target_index: Which state receives doses

Supports both bolus and infusion dosing.
"""
function simulate(
    spec::ModelSpec{<:CustomODE, CustomODEParams},
    grid::SimGrid,
    solver::SolverSpec
)
    validate(spec)
    validate(grid)
    validate(solver)

    kind = spec.kind
    p = pk_param_tuple(spec)
    tspan = (grid.t0, grid.t1)
    dose_target_idx = pk_dose_target_index(kind)

    # Check for infusion doses
    if _has_any_infusion(spec.doses)
        # Use infusion-aware simulation
        schedule = build_infusion_schedule(spec.doses, grid.t0, grid.t1)
        u0 = pk_u0(spec, grid)
        u0[dose_target_idx] += schedule.initial_amount

        rate_fn = make_infusion_rate_function(schedule)
        function ode_with_infusion!(du, u, params, t)
            pk_ode_with_infusion!(du, u, params, t, kind, rate_fn(t))
        end

        prob = ODEProblem(ode_with_infusion!, u0, tspan, p)
        tstops = schedule.rate_change_times

        sol = solve(
            prob,
            Rosenbrock23();
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            tstops=tstops,
        )
    else
        # Standard bolus simulation
        _, dose_times, dose_amounts = normalize_doses_for_sim(spec.doses, grid.t0, grid.t1)
        u0 = pk_u0(spec, grid)

        function ode!(du, u, params, t)
            pk_ode!(du, u, params, t, kind)
        end

        prob = ODEProblem(ode!, u0, tspan, p)

        cb = nothing
        if !isempty(dose_times)
            function affect!(integrator)
                idx = findfirst(==(integrator.t), dose_times)
                if idx !== nothing
                    integrator.u[dose_target_idx] += dose_amounts[idx]
                end
            end
            cb = PresetTimeCallback(dose_times, affect!)
        end

        tstops = filter(t -> t > grid.t0 && t < grid.t1, dose_times)

        sol = solve(
            prob,
            Rosenbrock23();
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=grid.saveat,
            callback=cb,
            tstops=tstops,
        )
    end

    # Build state dictionary using state names from CustomODE
    state_names = pk_state_symbols(kind)
    states = Dict{Symbol, Vector{Float64}}()
    for (i, name) in enumerate(state_names)
        states[name] = [u[i] for u in sol.u]
    end

    # Compute concentration using output function
    C = [pk_conc(u, p, kind) for u in sol.u]
    observations = Dict(:conc => C)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "CustomODE",
        "model_description" => kind.description,
        "n_states" => kind.n_states,
        "param_names" => [String(s) for s in kind.param_names],
        "state_names" => [String(s) for s in state_names],
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount, d.duration) for d in spec.doses],
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end
