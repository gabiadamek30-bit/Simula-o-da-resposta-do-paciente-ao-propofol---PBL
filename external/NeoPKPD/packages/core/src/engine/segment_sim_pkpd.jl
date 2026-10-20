export simulate_segmented_pkpd_coupled

"""
simulate_segmented_pkpd_coupled runs a coupled PKPD model across multiple time segments,
allowing PK parameters to change at segment boundaries while preserving
the global saveat grid and continuous PD response dynamics.

Inputs:
- pk_spec: ModelSpec with base PK params and dose schedule
- pd_spec: PDSpec for IndirectResponseTurnover
- grid: SimGrid with full [t0, t1] and saveat points
- solver: SolverSpec
- segment_starts: Vector of times including t0 and any segment boundaries
- pk_params_per_segment: Vector of typed PK params aligned with segments

Semantics:
- Doses are handled by event semantics per segment
- State continuity: final state (both PK and PD) of segment becomes initial state of next
- PD parameters remain constant across segments (only PK params vary with IOV)
- Returned times exactly equal grid.saveat
"""
function simulate_segmented_pkpd_coupled(
    pk_spec::ModelSpec{K,P},
    pd_spec::PDSpec{IndirectResponseTurnover,IndirectResponseTurnoverParams},
    grid::SimGrid,
    solver::SolverSpec,
    segment_starts::Vector{Float64},
    pk_params_per_segment::Vector{P},
) where {K<:ModelKind,P}
    if isempty(segment_starts)
        error("segment_starts must be non-empty")
    end
    if length(pk_params_per_segment) != length(segment_starts)
        error("pk_params_per_segment must align with segment_starts")
    end
    if !issorted(segment_starts)
        error("segment_starts must be sorted")
    end
    if segment_starts[1] != grid.t0
        error("First segment must start at grid.t0")
    end

    saveat = grid.saveat
    kind = pk_spec.kind

    # PD parameters (constant across segments)
    Kin = pd_spec.params.Kin
    Kout = pd_spec.params.Kout
    R0 = pd_spec.params.R0
    Imax = pd_spec.params.Imax
    IC50 = pd_spec.params.IC50

    # Initial PK state
    seg0 = SimGrid(segment_starts[1], grid.t1, saveat)
    u0_pk = pk_u0(pk_spec, seg0)

    # Apply t0 dosing only for the first segment
    a0_add, _, _ = normalize_doses_for_sim(pk_spec.doses, grid.t0, grid.t1)
    target_index = pk_dose_target_index(kind)
    u0_pk[target_index] += a0_add

    # Combined initial state: [PK states..., R]
    u0 = vcat(u0_pk, [R0])

    n_pk = length(u0_pk)
    r_index = n_pk + 1

    # Prepare outputs aligned with global saveat
    n_t = length(saveat)
    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}(
        sym => Vector{Float64}(undef, n_t) for sym in state_syms
    )
    states[:R] = Vector{Float64}(undef, n_t)
    conc = Vector{Float64}(undef, n_t)
    response = Vector{Float64}(undef, n_t)

    n_seg = length(segment_starts)

    for si in 1:n_seg
        t0 = segment_starts[si]
        t1 = (si < n_seg) ? segment_starts[si + 1] : grid.t1

        if t1 < t0
            error("Invalid segment boundary: t1 < t0 at segment $(si)")
        end

        # Indices of global saveat that fall in this segment
        idx = Int[]
        for (j, t) in enumerate(saveat)
            if t >= t0 && (t < t1 || (si == n_seg && t <= t1))
                push!(idx, j)
            end
        end

        if isempty(idx)
            continue
        end

        seg_saveat = saveat[idx]

        # Build typed params for this segment
        spec_seg = ModelSpec(kind, pk_spec.name, pk_params_per_segment[si], pk_spec.doses)
        pkp = pk_param_tuple(spec_seg)

        # Combined parameter tuple
        p = merge(pkp, (Kin=Kin, Kout=Kout, Imax=Imax, IC50=IC50))

        function ode!(du, u, p, t)
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

        prob = ODEProblem(ode!, u0, (t0, t1), p)

        # Doses within (t0, t1]
        cb = preset_dose_callback(pk_spec.doses, t0, t1, target_index)

        sol = solve(
            prob,
            _solver_alg(solver.alg);
            reltol=solver.reltol,
            abstol=solver.abstol,
            maxiters=solver.maxiters,
            saveat=seg_saveat,
            callback=cb,
        )

        # Fill outputs into global arrays
        for (local_i, global_j) in enumerate(idx)
            u = sol.u[local_i]
            for (k, sym) in enumerate(state_syms)
                states[sym][global_j] = u[k]
            end
            states[:R][global_j] = u[r_index]
            conc[global_j] = pk_conc(@view(u[1:n_pk]), p, kind)
            response[global_j] = u[r_index]
        end

        # Carry final state into next segment
        u0 = copy(sol.u[end])
    end

    outkey = pd_spec.output_observation
    observations = Dict{Symbol,Vector{Float64}}(:conc => conc, outkey => response)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Segmented Coupled PKPD",
        "pk_kind" => string(typeof(kind)),
        "pd_kind" => "IndirectResponseTurnover",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "pd_output_observation" => String(outkey),
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
        "deterministic_output_grid" => true,
    )

    return SimResult(saveat, states, observations, metadata)
end
