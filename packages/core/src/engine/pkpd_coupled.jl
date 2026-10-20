export simulate_pkpd_coupled

# Union type for all IRM models that can use coupled simulation
const IRMModelType = Union{
    IndirectResponseTurnover,  # IRM-III (alias)
    IndirectResponseIRM1,
    IndirectResponseIRM2,
    IndirectResponseIRM4,
}

function _preset_dose_callback(
    doses::Vector{DoseEvent}, t0::Float64, t1::Float64, target_index::Int
)
    _, dose_times, dose_amounts = normalize_doses_for_sim(doses, t0, t1)

    if isempty(dose_times)
        return nothing
    end

    function affect!(integrator)
        idx = findfirst(==(integrator.t), dose_times)
        if idx === nothing
            error("Internal error: dose time not found for t=$(integrator.t)")
        end
        integrator.u[target_index] += dose_amounts[idx]
    end

    return PresetTimeCallback(dose_times, affect!)
end

"""
General coupled PKPD simulation for:
- PK: any supported ModelSpec{K,P}
- PD: IndirectResponseTurnover

Coupled states:
- first: PK states
- last:  R response

Observations:
- :conc computed by pk_conc
- pd_spec.output_observation maps to response state R
"""
function simulate_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{IndirectResponseTurnover,IndirectResponseTurnoverParams},
    grid::SimGrid,
    solver::SolverSpec,
) where {K<:ModelKind,P}
    pk_validate(pk_spec)
    validate(pd_spec)
    validate(grid)
    validate(solver)

    kind = pk_spec.kind

    # PK parameter tuple is model-specific
    pkp = pk_param_tuple(pk_spec)

    # PD parameters
    Kin = pd_spec.params.Kin
    Kout = pd_spec.params.Kout
    R0 = pd_spec.params.R0
    Imax = pd_spec.params.Imax
    IC50 = pd_spec.params.IC50

    # Build initial state vector: [PK states..., R]
    a0_add, _, _ = normalize_doses_for_sim(pk_spec.doses, grid.t0, grid.t1)

    u0_pk = pk_u0(pk_spec, grid)

    # Apply the t0 dose add to the PK dose target index
    target_index = pk_dose_target_index(pk_spec.kind)
    u0_pk[target_index] += a0_add

    u0 = vcat(u0_pk, [R0])

    # Combined parameter tuple (NamedTuple), still pure data
    p = merge(pkp, (Kin=Kin, Kout=Kout, Imax=Imax, IC50=IC50))

    n_pk = length(u0_pk)
    r_index = n_pk + 1

    function ode!(du, u, p, t)
        # Split views
        u_pk = @view u[1:n_pk]
        du_pk = @view du[1:n_pk]

        # PK dynamics
        pk_ode!(du_pk, u_pk, p, t, kind)

        # PD dynamics
        R = u[r_index]
        C = pk_conc(u_pk, p, kind)
        I = inhibition(C, p.Imax, p.IC50)
        du[r_index] = p.Kin - p.Kout * (1.0 - I) * R

        return nothing
    end

    prob = ODEProblem(ode!, u0, (grid.t0, grid.t1), p)

    target_index = pk_dose_target_index(kind)
    cb = _preset_dose_callback(pk_spec.doses, grid.t0, grid.t1, target_index)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    # Build state outputs
    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}()

    for (j, sym) in enumerate(state_syms)
        states[sym] = [u[j] for u in sol.u]
    end

    R = [u[r_index] for u in sol.u]
    states[:R] = R

    # Observations
    C = [pk_conc(@view(u[1:n_pk]), p, kind) for u in sol.u]

    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => C, outkey => R)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "IndirectResponseTurnover",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pk_dose_schedule" => [(d.time, d.amount) for d in pk_spec.doses],
        "pd_output_observation" => String(outkey),
        "pd_params" =>
            Dict("Kin" => Kin, "Kout" => Kout, "R0" => R0, "Imax" => Imax, "IC50" => IC50),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

"""
IRM-I: Inhibition of Kin (production).
dR/dt = Kin × (1 - I(C)) - Kout × R

Drug inhibits production → Response decreases below baseline.
"""
function simulate_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{IndirectResponseIRM1,IndirectResponseIRM1Params},
    grid::SimGrid,
    solver::SolverSpec,
) where {K<:ModelKind,P}
    pk_validate(pk_spec)
    validate(pd_spec)
    validate(grid)
    validate(solver)

    kind = pk_spec.kind
    pkp = pk_param_tuple(pk_spec)

    # PD parameters
    Kin = pd_spec.params.Kin
    Kout = pd_spec.params.Kout
    R0 = pd_spec.params.R0
    Imax = pd_spec.params.Imax
    IC50 = pd_spec.params.IC50

    # Build initial state vector: [PK states..., R]
    a0_add, _, _ = normalize_doses_for_sim(pk_spec.doses, grid.t0, grid.t1)
    u0_pk = pk_u0(pk_spec, grid)

    target_index = pk_dose_target_index(pk_spec.kind)
    u0_pk[target_index] += a0_add

    u0 = vcat(u0_pk, [R0])

    # Combined parameter tuple
    p = merge(pkp, (Kin=Kin, Kout=Kout, Imax=Imax, IC50=IC50))

    n_pk = length(u0_pk)
    r_index = n_pk + 1

    function ode!(du, u, p, t)
        u_pk = @view u[1:n_pk]
        du_pk = @view du[1:n_pk]

        pk_ode!(du_pk, u_pk, p, t, kind)

        R = u[r_index]
        C = pk_conc(u_pk, p, kind)
        I = inhibition(C, p.Imax, p.IC50)

        # IRM-I: Inhibition of Kin
        du[r_index] = p.Kin * (1.0 - I) - p.Kout * R

        return nothing
    end

    prob = ODEProblem(ode!, u0, (grid.t0, grid.t1), p)

    target_index = pk_dose_target_index(kind)
    cb = _preset_dose_callback(pk_spec.doses, grid.t0, grid.t1, target_index)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    # Build outputs
    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}()

    for (j, sym) in enumerate(state_syms)
        states[sym] = [u[j] for u in sol.u]
    end

    R = [u[r_index] for u in sol.u]
    states[:R] = R

    C = [pk_conc(@view(u[1:n_pk]), p, kind) for u in sol.u]
    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => C, outkey => R)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "IndirectResponseIRM1",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pk_dose_schedule" => [(d.time, d.amount) for d in pk_spec.doses],
        "pd_output_observation" => String(outkey),
        "pd_params" =>
            Dict("Kin" => Kin, "Kout" => Kout, "R0" => R0, "Imax" => Imax, "IC50" => IC50),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

"""
IRM-II: Stimulation of Kin (production).
dR/dt = Kin × (1 + S(C)) - Kout × R

Drug stimulates production → Response increases above baseline.
"""
function simulate_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{IndirectResponseIRM2,IndirectResponseIRM2Params},
    grid::SimGrid,
    solver::SolverSpec,
) where {K<:ModelKind,P}
    pk_validate(pk_spec)
    validate(pd_spec)
    validate(grid)
    validate(solver)

    kind = pk_spec.kind
    pkp = pk_param_tuple(pk_spec)

    # PD parameters
    Kin = pd_spec.params.Kin
    Kout = pd_spec.params.Kout
    R0 = pd_spec.params.R0
    Smax = pd_spec.params.Smax
    SC50 = pd_spec.params.SC50

    # Build initial state vector
    a0_add, _, _ = normalize_doses_for_sim(pk_spec.doses, grid.t0, grid.t1)
    u0_pk = pk_u0(pk_spec, grid)

    target_index = pk_dose_target_index(pk_spec.kind)
    u0_pk[target_index] += a0_add

    u0 = vcat(u0_pk, [R0])

    p = merge(pkp, (Kin=Kin, Kout=Kout, Smax=Smax, SC50=SC50))

    n_pk = length(u0_pk)
    r_index = n_pk + 1

    function ode!(du, u, p, t)
        u_pk = @view u[1:n_pk]
        du_pk = @view du[1:n_pk]

        pk_ode!(du_pk, u_pk, p, t, kind)

        R = u[r_index]
        C = pk_conc(u_pk, p, kind)
        S = stimulation(C, p.Smax, p.SC50)

        # IRM-II: Stimulation of Kin
        du[r_index] = p.Kin * (1.0 + S) - p.Kout * R

        return nothing
    end

    prob = ODEProblem(ode!, u0, (grid.t0, grid.t1), p)

    target_index = pk_dose_target_index(kind)
    cb = _preset_dose_callback(pk_spec.doses, grid.t0, grid.t1, target_index)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}()

    for (j, sym) in enumerate(state_syms)
        states[sym] = [u[j] for u in sol.u]
    end

    R = [u[r_index] for u in sol.u]
    states[:R] = R

    C = [pk_conc(@view(u[1:n_pk]), p, kind) for u in sol.u]
    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => C, outkey => R)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "IndirectResponseIRM2",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pk_dose_schedule" => [(d.time, d.amount) for d in pk_spec.doses],
        "pd_output_observation" => String(outkey),
        "pd_params" =>
            Dict("Kin" => Kin, "Kout" => Kout, "R0" => R0, "Smax" => Smax, "SC50" => SC50),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

"""
IRM-IV: Stimulation of Kout (elimination).
dR/dt = Kin - Kout × (1 + S(C)) × R

Drug stimulates elimination → Response decreases below baseline.
"""
function simulate_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{IndirectResponseIRM4,IndirectResponseIRM4Params},
    grid::SimGrid,
    solver::SolverSpec,
) where {K<:ModelKind,P}
    pk_validate(pk_spec)
    validate(pd_spec)
    validate(grid)
    validate(solver)

    kind = pk_spec.kind
    pkp = pk_param_tuple(pk_spec)

    # PD parameters
    Kin = pd_spec.params.Kin
    Kout = pd_spec.params.Kout
    R0 = pd_spec.params.R0
    Smax = pd_spec.params.Smax
    SC50 = pd_spec.params.SC50

    # Build initial state vector
    a0_add, _, _ = normalize_doses_for_sim(pk_spec.doses, grid.t0, grid.t1)
    u0_pk = pk_u0(pk_spec, grid)

    target_index = pk_dose_target_index(pk_spec.kind)
    u0_pk[target_index] += a0_add

    u0 = vcat(u0_pk, [R0])

    p = merge(pkp, (Kin=Kin, Kout=Kout, Smax=Smax, SC50=SC50))

    n_pk = length(u0_pk)
    r_index = n_pk + 1

    function ode!(du, u, p, t)
        u_pk = @view u[1:n_pk]
        du_pk = @view du[1:n_pk]

        pk_ode!(du_pk, u_pk, p, t, kind)

        R = u[r_index]
        C = pk_conc(u_pk, p, kind)
        S = stimulation(C, p.Smax, p.SC50)

        # IRM-IV: Stimulation of Kout
        du[r_index] = p.Kin - p.Kout * (1.0 + S) * R

        return nothing
    end

    prob = ODEProblem(ode!, u0, (grid.t0, grid.t1), p)

    target_index = pk_dose_target_index(kind)
    cb = _preset_dose_callback(pk_spec.doses, grid.t0, grid.t1, target_index)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}()

    for (j, sym) in enumerate(state_syms)
        states[sym] = [u[j] for u in sol.u]
    end

    R = [u[r_index] for u in sol.u]
    states[:R] = R

    C = [pk_conc(@view(u[1:n_pk]), p, kind) for u in sol.u]
    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => C, outkey => R)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "IndirectResponseIRM4",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pk_dose_schedule" => [(d.time, d.amount) for d in pk_spec.doses],
        "pd_output_observation" => String(outkey),
        "pd_params" =>
            Dict("Kin" => Kin, "Kout" => Kout, "R0" => R0, "Smax" => Smax, "SC50" => SC50),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

"""
Transit Compartment PD: Signal transduction delay model.
dA1/dt = ktr × (Signal(C) - A1)
dAi/dt = ktr × (A(i-1) - Ai)  for i = 2..N
Effect = AN

Signal(C) = E0 + Emax × C^γ / (EC50^γ + C^γ)
MTT = (N + 1) / ktr
"""
function simulate_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{TransitCompartmentPD,TransitCompartmentPDParams},
    grid::SimGrid,
    solver::SolverSpec,
) where {K<:ModelKind,P}
    pk_validate(pk_spec)
    validate(pd_spec)
    validate(grid)
    validate(solver)

    kind = pk_spec.kind
    pkp = pk_param_tuple(pk_spec)

    # PD parameters
    N = pd_spec.params.N
    ktr = pd_spec.params.ktr
    E0 = pd_spec.params.E0
    Emax = pd_spec.params.Emax
    EC50 = pd_spec.params.EC50
    gamma = pd_spec.params.gamma

    # Build initial state vector: [PK states..., A1, A2, ..., AN]
    a0_add, _, _ = normalize_doses_for_sim(pk_spec.doses, grid.t0, grid.t1)
    u0_pk = pk_u0(pk_spec, grid)

    target_index = pk_dose_target_index(pk_spec.kind)
    u0_pk[target_index] += a0_add

    # Initialize transit compartments at baseline (steady state when C=0)
    u0_transit = fill(E0, N)
    u0 = vcat(u0_pk, u0_transit)

    # Combined parameter tuple
    p = merge(pkp, (N=N, ktr=ktr, E0=E0, Emax=Emax, EC50=EC50, gamma=gamma))

    n_pk = length(u0_pk)

    function ode!(du, u, p, t)
        u_pk = @view u[1:n_pk]
        du_pk = @view du[1:n_pk]

        # PK dynamics
        pk_ode!(du_pk, u_pk, p, t, kind)

        # Calculate drug concentration
        C = pk_conc(u_pk, p, kind)

        # Calculate input signal using sigmoid Emax model
        signal = sigmoid_emax_effect(C, p.E0, p.Emax, p.EC50, p.gamma)

        # Transit compartment dynamics
        # First compartment receives signal
        A1 = u[n_pk + 1]
        du[n_pk + 1] = p.ktr * (signal - A1)

        # Subsequent compartments
        for i in 2:p.N
            A_prev = u[n_pk + i - 1]
            A_curr = u[n_pk + i]
            du[n_pk + i] = p.ktr * (A_prev - A_curr)
        end

        return nothing
    end

    prob = ODEProblem(ode!, u0, (grid.t0, grid.t1), p)

    target_index = pk_dose_target_index(kind)
    cb = _preset_dose_callback(pk_spec.doses, grid.t0, grid.t1, target_index)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    # Build state outputs
    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}()

    for (j, sym) in enumerate(state_syms)
        states[sym] = [u[j] for u in sol.u]
    end

    # Add transit compartment states
    for i in 1:N
        states[Symbol("A$i")] = [u[n_pk + i] for u in sol.u]
    end

    # Effect is the last transit compartment
    effect_idx = n_pk + N
    effect = [u[effect_idx] for u in sol.u]

    C = [pk_conc(@view(u[1:n_pk]), p, kind) for u in sol.u]
    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => C, outkey => effect)

    # Calculate MTT for metadata
    mtt = (N + 1) / ktr

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "TransitCompartmentPD",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pk_dose_schedule" => [(d.time, d.amount) for d in pk_spec.doses],
        "pd_output_observation" => String(outkey),
        "pd_params" => Dict(
            "N" => N, "ktr" => ktr, "E0" => E0,
            "Emax" => Emax, "EC50" => EC50, "gamma" => gamma,
            "MTT" => mtt
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

"""
Disease Progression PD: Tumor growth with drug effect.
Growth models: Linear, Asymptotic, Gompertz, Logistic, Exponential
Drug effect: Cell-kill model (kdrug × C × S)
"""
function simulate_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{DiseaseProgressionPD,DiseaseProgressionPDParams},
    grid::SimGrid,
    solver::SolverSpec,
) where {K<:ModelKind,P}
    pk_validate(pk_spec)
    validate(pd_spec)
    validate(grid)
    validate(solver)

    kind = pk_spec.kind
    pkp = pk_param_tuple(pk_spec)

    # PD parameters
    growth_model = pd_spec.kind.growth_model
    S0 = pd_spec.params.S0
    kgrow = pd_spec.params.kgrow
    Smax = pd_spec.params.Smax
    alpha = pd_spec.params.alpha
    kdrug = pd_spec.params.kdrug

    # Build initial state vector: [PK states..., S]
    a0_add, _, _ = normalize_doses_for_sim(pk_spec.doses, grid.t0, grid.t1)
    u0_pk = pk_u0(pk_spec, grid)

    target_index = pk_dose_target_index(pk_spec.kind)
    u0_pk[target_index] += a0_add

    u0 = vcat(u0_pk, [S0])

    # Combined parameter tuple
    p = merge(pkp, (growth_model=growth_model, kgrow=kgrow, Smax=Smax, alpha=alpha, kdrug=kdrug))

    n_pk = length(u0_pk)
    s_index = n_pk + 1

    function ode!(du, u, p, t)
        u_pk = @view u[1:n_pk]
        du_pk = @view du[1:n_pk]

        # PK dynamics
        pk_ode!(du_pk, u_pk, p, t, kind)

        # Tumor dynamics
        S = u[s_index]
        C = pk_conc(u_pk, p, kind)

        # Natural growth rate
        gr = growth_rate(S, p.growth_model, p.kgrow, p.Smax, p.alpha)

        # Drug-induced cell kill
        de = drug_effect(C, S, p.kdrug)

        # Tumor ODE: dS/dt = growth - drug_effect
        # Ensure tumor size doesn't go negative
        du[s_index] = max(gr - de, -S / 0.1)  # Soft constraint to prevent negative

        return nothing
    end

    prob = ODEProblem(ode!, u0, (grid.t0, grid.t1), p)

    target_index = pk_dose_target_index(kind)
    cb = _preset_dose_callback(pk_spec.doses, grid.t0, grid.t1, target_index)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    # Build state outputs
    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}()

    for (j, sym) in enumerate(state_syms)
        states[sym] = [u[j] for u in sol.u]
    end

    # Tumor size state
    tumor_size = [max(u[s_index], 0.0) for u in sol.u]  # Ensure non-negative
    states[:S] = tumor_size

    C = [pk_conc(@view(u[1:n_pk]), p, kind) for u in sol.u]
    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => C, outkey => tumor_size)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "DiseaseProgressionPD",
        "growth_model" => string(growth_model),
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pk_dose_schedule" => [(d.time, d.amount) for d in pk_spec.doses],
        "pd_output_observation" => String(outkey),
        "pd_params" => Dict(
            "S0" => S0, "kgrow" => kgrow, "Smax" => Smax,
            "alpha" => alpha, "kdrug" => kdrug
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

"""
Tolerance Counter-Regulation PD: Feedback-based tolerance model.

States:
- M: Moderator/feedback variable (starts at 0)

Dynamics:
- dM/dt = kin_mod × E_drug - kout_mod × M
- E_drug = Emax × C^γ / (EC50^γ + C^γ)
- E_net = E0 + E_drug - alpha × M

The moderator M accumulates with drug exposure and reduces net effect.
"""
function simulate_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{ToleranceCounterRegulation,ToleranceCounterRegulationParams},
    grid::SimGrid,
    solver::SolverSpec,
) where {K<:ModelKind,P}
    pk_validate(pk_spec)
    validate(pd_spec)
    validate(grid)
    validate(solver)

    kind = pk_spec.kind
    pkp = pk_param_tuple(pk_spec)

    # PD parameters
    E0 = pd_spec.params.E0
    Emax = pd_spec.params.Emax
    EC50 = pd_spec.params.EC50
    gamma = pd_spec.params.gamma
    kin_mod = pd_spec.params.kin_mod
    kout_mod = pd_spec.params.kout_mod
    alpha = pd_spec.params.alpha

    # Build initial state vector: [PK states..., M]
    a0_add, _, _ = normalize_doses_for_sim(pk_spec.doses, grid.t0, grid.t1)
    u0_pk = pk_u0(pk_spec, grid)

    target_index = pk_dose_target_index(pk_spec.kind)
    u0_pk[target_index] += a0_add

    # Moderator starts at 0 (no tolerance at baseline)
    M0 = 0.0
    u0 = vcat(u0_pk, [M0])

    # Combined parameter tuple
    p = merge(pkp, (E0=E0, Emax=Emax, EC50=EC50, gamma=gamma, kin_mod=kin_mod, kout_mod=kout_mod, alpha=alpha))

    n_pk = length(u0_pk)
    m_index = n_pk + 1

    function ode!(du, u, p, t)
        u_pk = @view u[1:n_pk]
        du_pk = @view du[1:n_pk]

        # PK dynamics
        pk_ode!(du_pk, u_pk, p, t, kind)

        # Drug effect calculation
        C = pk_conc(u_pk, p, kind)
        E_drug = sigmoid_emax_effect(C, 0.0, p.Emax, p.EC50, p.gamma)  # E0=0 for drug component

        # Moderator dynamics
        M = u[m_index]
        du[m_index] = p.kin_mod * E_drug - p.kout_mod * M

        return nothing
    end

    prob = ODEProblem(ode!, u0, (grid.t0, grid.t1), p)

    target_index = pk_dose_target_index(kind)
    cb = _preset_dose_callback(pk_spec.doses, grid.t0, grid.t1, target_index)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    # Build state outputs
    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}()

    for (j, sym) in enumerate(state_syms)
        states[sym] = [u[j] for u in sol.u]
    end

    # Moderator state
    moderator = [u[m_index] for u in sol.u]
    states[:M] = moderator

    # Calculate net effect: E_net = E0 + E_drug - alpha × M
    C = [pk_conc(@view(u[1:n_pk]), p, kind) for u in sol.u]
    E_drug = [sigmoid_emax_effect(c, 0.0, p.Emax, p.EC50, p.gamma) for c in C]
    effect = [p.E0 + ed - p.alpha * m for (ed, m) in zip(E_drug, moderator)]

    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => C, outkey => effect, :E_drug => E_drug)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "ToleranceCounterRegulation",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pk_dose_schedule" => [(d.time, d.amount) for d in pk_spec.doses],
        "pd_output_observation" => String(outkey),
        "pd_params" => Dict(
            "E0" => E0, "Emax" => Emax, "EC50" => EC50, "gamma" => gamma,
            "kin_mod" => kin_mod, "kout_mod" => kout_mod, "alpha" => alpha
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end

"""
Receptor Regulation PD: Tolerance through receptor up/down-regulation.

States:
- R: Receptor density (normalized, baseline = R_baseline, typically 1.0)

Dynamics (down-regulation):
- dR/dt = kreg × (R_baseline - R) - kchange × f(C) × R

Dynamics (up-regulation):
- dR/dt = kreg × (R_baseline - R) + kchange × f(C) × (Rmax - R)

Where f(C) = C^gamma / (EC50^gamma + C^gamma) (fractional drug effect)

Net effect:
- E_net = E0 + R × Emax × f(C)
"""
function simulate_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{ReceptorRegulation,ReceptorRegulationParams},
    grid::SimGrid,
    solver::SolverSpec,
) where {K<:ModelKind,P}
    pk_validate(pk_spec)
    validate(pd_spec)
    validate(grid)
    validate(solver)

    kind = pk_spec.kind
    pkp = pk_param_tuple(pk_spec)

    # PD parameters
    E0 = pd_spec.params.E0
    Emax = pd_spec.params.Emax
    EC50 = pd_spec.params.EC50
    gamma = pd_spec.params.gamma
    R_baseline = pd_spec.params.R_baseline
    kreg = pd_spec.params.kreg
    Rmax = pd_spec.params.Rmax
    kchange = pd_spec.params.kchange
    direction = pd_spec.params.direction

    # Build initial state vector: [PK states..., R]
    a0_add, _, _ = normalize_doses_for_sim(pk_spec.doses, grid.t0, grid.t1)
    u0_pk = pk_u0(pk_spec, grid)

    target_index = pk_dose_target_index(pk_spec.kind)
    u0_pk[target_index] += a0_add

    # Receptor starts at baseline
    u0 = vcat(u0_pk, [R_baseline])

    # Combined parameter tuple
    p = merge(pkp, (E0=E0, Emax=Emax, EC50=EC50, gamma=gamma, R_baseline=R_baseline, kreg=kreg, Rmax=Rmax, kchange=kchange, direction=direction))

    n_pk = length(u0_pk)
    r_index = n_pk + 1

    function ode!(du, u, p, t)
        u_pk = @view u[1:n_pk]
        du_pk = @view du[1:n_pk]

        # PK dynamics
        pk_ode!(du_pk, u_pk, p, t, kind)

        # Fractional drug effect
        C = pk_conc(u_pk, p, kind)
        if C <= 0.0
            f_drug = 0.0
        else
            ratio = (p.EC50 / C)^p.gamma
            f_drug = 1.0 / (1.0 + ratio)
        end

        # Receptor dynamics
        R = u[r_index]
        baseline_return = p.kreg * (p.R_baseline - R)

        if p.direction == :down
            # Down-regulation: drug exposure decreases receptor density
            regulation_effect = -p.kchange * f_drug * R
        else  # :up
            # Up-regulation: drug exposure increases receptor density (toward Rmax)
            regulation_effect = p.kchange * f_drug * (p.Rmax - R)
        end

        du[r_index] = baseline_return + regulation_effect

        return nothing
    end

    prob = ODEProblem(ode!, u0, (grid.t0, grid.t1), p)

    target_index = pk_dose_target_index(kind)
    cb = _preset_dose_callback(pk_spec.doses, grid.t0, grid.t1, target_index)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    # Build state outputs
    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}()

    for (j, sym) in enumerate(state_syms)
        states[sym] = [u[j] for u in sol.u]
    end

    # Receptor state
    receptor = [max(u[r_index], 0.0) for u in sol.u]  # Ensure non-negative
    states[:R] = receptor

    # Calculate net effect: E_net = E0 + R × Emax × f(C)
    C = [pk_conc(@view(u[1:n_pk]), p, kind) for u in sol.u]

    effect = Vector{Float64}(undef, length(C))
    for i in eachindex(C)
        c = C[i]
        if c <= 0.0
            f_drug = 0.0
        else
            ratio = (p.EC50 / c)^p.gamma
            f_drug = 1.0 / (1.0 + ratio)
        end
        effect[i] = p.E0 + receptor[i] * p.Emax * f_drug
    end

    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => C, outkey => effect, :receptor_density => receptor)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "ReceptorRegulation",
        "direction" => String(direction),
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pk_dose_schedule" => [(d.time, d.amount) for d in pk_spec.doses],
        "pd_output_observation" => String(outkey),
        "pd_params" => Dict(
            "E0" => E0, "Emax" => Emax, "EC50" => EC50, "gamma" => gamma,
            "R_baseline" => R_baseline, "kreg" => kreg, "Rmax" => Rmax, "kchange" => kchange
        ),
        "deterministic_output_grid" => true,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return SimResult(Vector{Float64}(sol.t), states, observations, metadata)
end
