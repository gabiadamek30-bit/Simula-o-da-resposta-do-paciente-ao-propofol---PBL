export simulate_segmented_pk

"""
simulate_segmented_pk runs a PK model across multiple time segments,
allowing parameters to change at segment boundaries while preserving
the global saveat grid.

Inputs:
- spec_base: ModelSpec with base params and dose schedule
- grid: SimGrid with full [t0, t1] and saveat points
- solver: SolverSpec
- segment_starts: Vector of times including t0 and any segment boundaries
- params_per_segment: Vector of typed params aligned with segments

Semantics:
- Doses are handled by event semantics per segment using the same normalize_doses_for_sim
- State continuity: final state of segment becomes initial state of next segment
- Returned times exactly equal grid.saveat
"""
function simulate_segmented_pk(
    spec_base::ModelSpec{K,P},
    grid::SimGrid,
    solver::SolverSpec,
    segment_starts::Vector{Float64},
    params_per_segment::Vector{P},
) where {K<:ModelKind,P}
    if isempty(segment_starts)
        error("segment_starts must be non-empty")
    end
    if length(params_per_segment) != length(segment_starts)
        error("params_per_segment must align with segment_starts")
    end
    if !issorted(segment_starts)
        error("segment_starts must be sorted")
    end
    if segment_starts[1] != grid.t0
        error("First segment must start at grid.t0")
    end

    # Precompute global saveat (assumed validated + sorted by SimGrid validate)
    saveat = grid.saveat

    kind = spec_base.kind

    # Initial state for first segment
    seg0 = SimGrid(segment_starts[1], grid.t1, saveat)
    u0 = pk_u0(spec_base, seg0)

    # Apply t0 dosing only for the first segment
    a0_add, _, _ = normalize_doses_for_sim(spec_base.doses, grid.t0, grid.t1)
    u0[pk_dose_target_index(kind)] += a0_add

    # Prepare outputs aligned with global saveat
    n_t = length(saveat)
    state_syms = pk_state_symbols(kind)
    states = Dict{Symbol,Vector{Float64}}(
        sym => Vector{Float64}(undef, n_t) for sym in state_syms
    )
    conc = Vector{Float64}(undef, n_t)

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
        seg_grid = SimGrid(t0, t1, seg_saveat)

        # Params for this segment (typed)
        spec_seg = ModelSpec(kind, spec_base.name, params_per_segment[si], spec_base.doses)
        p_tuple = pk_param_tuple(spec_seg)

        function ode!(du, u, p, t)
            pk_ode!(du, u, p, t, kind)
            return nothing
        end

        prob = ODEProblem(ode!, u0, (t0, t1), p_tuple)

        # Doses within (t0, t1]
        target_index = pk_dose_target_index(kind)
        cb = preset_dose_callback(spec_base.doses, t0, t1, target_index)

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
            conc[global_j] = pk_conc(@view(u[1:length(state_syms)]), p_tuple, kind)
        end

        # Carry final state into next segment
        u0 = copy(sol.u[end])
    end

    observations = Dict{Symbol,Vector{Float64}}(:conc => conc)

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "model" => "Segmented PK",
        "pk_kind" => string(typeof(kind)),
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
        "deterministic_output_grid" => true,
    )

    return SimResult(saveat, states, observations, metadata)
end
