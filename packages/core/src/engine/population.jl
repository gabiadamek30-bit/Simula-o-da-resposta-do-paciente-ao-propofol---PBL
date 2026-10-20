using StableRNGs
using Distributions

export PopulationResult, PopulationSummary, simulate_population

"""
Population summary statistics for a single observation series.

Each vector is aligned with time grid t.
"""
struct PopulationSummary
    observation::Symbol
    probs::Vector{Float64}
    mean::Vector{Float64}
    median::Vector{Float64}
    quantiles::Dict{Float64,Vector{Float64}}
end

"""
Population simulation output.

individuals:
- Vector of SimResult, one per individual

params:
- Vector of Dicts describing realized individual parameters

metadata:
- includes seed, n, iiv kind, and omega map
"""
struct PopulationResult
    individuals::Vector{SimResult}
    params::Vector{Dict{Symbol,Float64}}
    summaries::Dict{Symbol,PopulationSummary}
    metadata::Dict{String,Any}
end

function validate(iiv::IIVSpec{LogNormalIIV})
    if iiv.n < 1
        error("IIVSpec n must be >= 1")
    end
    for (k, ω) in iiv.omegas
        _require_positive("omega for $(k)", ω)
    end
    return nothing
end

"""
Sample eta values for IIV using diagonal (independent) omegas.
Backward compatible - samples from independent Normal distributions.
"""
function _sample_etas_log_normal(base_params, omegas::Dict{Symbol,Float64}, rng)
    T = typeof(base_params)
    fn = fieldnames(T)

    etas = Dict{Symbol,Float64}()
    for f in fn
        if haskey(omegas, f)
            ω = omegas[f]
            etas[f] = rand(rng, Normal(0.0, ω))
        end
    end
    return etas
end

"""
Sample correlated eta values using full covariance matrix.

Uses Cholesky decomposition: η = L * z where z ~ N(0, I)
This properly samples from the multivariate normal distribution.
"""
function _sample_etas_correlated(omega_matrix::OmegaMatrix, rng)
    n = length(omega_matrix.param_names)

    # Sample from standard normal
    z = randn(rng, n)

    # Transform to correlated etas using Cholesky factor
    eta_vec = omega_matrix.cholesky_L * z

    # Convert to Dict
    etas = Dict{Symbol,Float64}()
    for (i, p) in enumerate(omega_matrix.param_names)
        etas[p] = eta_vec[i]
    end

    return etas
end

"""
Sample etas from IIVSpec - uses correlated sampling if available.
"""
function _sample_etas_from_iiv(iiv::IIVSpec{LogNormalIIV}, base_params, rng)
    if iiv.omega_matrix !== nothing
        # Use correlated sampling with full covariance matrix
        return _sample_etas_correlated(iiv.omega_matrix, rng)
    else
        # Fall back to independent sampling
        return _sample_etas_log_normal(base_params, iiv.omegas, rng)
    end
end

"""
Sample etas as a vector (ordered by param_names) for estimation.
Returns (eta_vector, param_names).
"""
function sample_eta_vector(omega_matrix::OmegaMatrix, rng)::Tuple{Vector{Float64}, Vector{Symbol}}
    n = length(omega_matrix.param_names)
    z = randn(rng, n)
    eta_vec = omega_matrix.cholesky_L * z
    return (eta_vec, omega_matrix.param_names)
end

"""
Apply pre-sampled etas to parameters.
theta_i = theta * exp(eta)
"""
function _apply_etas_log_normal(params, etas::Dict{Symbol,Float64})
    T = typeof(params)
    fn = fieldnames(T)

    vals = Dict{Symbol,Float64}()
    for f in fn
        θ = Float64(getfield(params, f))
        if haskey(etas, f)
            θ = θ * exp(etas[f])
        end
        vals[f] = θ
    end

    return T((vals[f] for f in fn)...), vals
end

# Keep the old function for backward compatibility with existing code paths
function _realize_params_log_normal(base_params, omegas::Dict{Symbol,Float64}, rng)
    etas = _sample_etas_log_normal(base_params, omegas, rng)
    return _apply_etas_log_normal(base_params, etas)
end

"""
Realize parameters with full IIV support (including correlations).
"""
function _realize_params_with_iiv(base_params, iiv::IIVSpec{LogNormalIIV}, rng)
    etas = _sample_etas_from_iiv(iiv, base_params, rng)
    return _apply_etas_log_normal(base_params, etas)
end

function _mean(xs::Vector{Float64})
    s = 0.0
    for x in xs
        s += x
    end
    return s / length(xs)
end

function _median_sorted(xs_sorted::Vector{Float64})
    n = length(xs_sorted)
    if isodd(n)
        return xs_sorted[(n + 1) ÷ 2]
    end
    return 0.5 * (xs_sorted[n ÷ 2] + xs_sorted[n ÷ 2 + 1])
end

function _quantile_sorted(xs_sorted::Vector{Float64}, p::Float64)
    n = length(xs_sorted)
    if n == 1
        return xs_sorted[1]
    end
    if p <= 0.0
        return xs_sorted[1]
    end
    if p >= 1.0
        return xs_sorted[end]
    end

    h = 1.0 + (n - 1) * p
    lo = floor(Int, h)
    hi = ceil(Int, h)

    if lo == hi
        return xs_sorted[lo]
    end

    w = h - lo
    return (1.0 - w) * xs_sorted[lo] + w * xs_sorted[hi]
end

function compute_population_summary(
    individuals::Vector{SimResult}, observation::Symbol; probs::Vector{Float64}=[0.05, 0.95]
)
    if isempty(individuals)
        error("Cannot compute summary for empty individuals")
    end

    t = individuals[1].t
    n_t = length(t)

    for (i, ind) in enumerate(individuals)
        if ind.t != t
            error(
                "All individuals must share identical time grid for summaries. Mismatch at index $(i)",
            )
        end
        if !haskey(ind.observations, observation)
            error("Observation $(observation) missing for individual $(i)")
        end
    end

    mean_v = Vector{Float64}(undef, n_t)
    median_v = Vector{Float64}(undef, n_t)
    qmap = Dict{Float64,Vector{Float64}}(p => Vector{Float64}(undef, n_t) for p in probs)

    scratch = Vector{Float64}(undef, length(individuals))

    for ti in 1:n_t
        for (j, ind) in enumerate(individuals)
            scratch[j] = ind.observations[observation][ti]
        end

        sort!(scratch)

        mean_v[ti] = _mean(scratch)
        median_v[ti] = _median_sorted(scratch)

        for p in probs
            qmap[p][ti] = _quantile_sorted(scratch, p)
        end
    end

    return PopulationSummary(observation, probs, mean_v, median_v, qmap)
end

"""
Build unified segment boundaries from IOV and time-varying covariates.
Returns sorted unique times within [t0, t1].
"""
function _build_segment_starts(
    grid::SimGrid,
    doses::Vector{DoseEvent},
    iov::Union{Nothing,IOVSpec},
    covs_i::Union{Nothing,IndividualCovariates},
    has_covariate_model::Bool,
)
    segment_starts = Float64[grid.t0]

    # Add IOV occasion boundaries
    if iov !== nothing
        occ = derive_occasions(doses, grid, iov.occasion_def)
        append!(segment_starts, occ)
    end

    # Add time-varying covariate boundaries
    if has_covariate_model && covs_i !== nothing && covs_i.time_varying !== nothing
        tv_times = covariate_boundary_times(covs_i.time_varying)
        append!(segment_starts, tv_times)
    end

    sort!(segment_starts)

    # unique and within bounds
    uniq = Float64[]
    for t in segment_starts
        if t >= grid.t0 && t <= grid.t1
            if isempty(uniq) || t != uniq[end]
                push!(uniq, t)
            end
        end
    end

    return uniq
end

"""
Check if segmentation is needed for an individual.
"""
function _needs_segmentation(
    iov::Union{Nothing,IOVSpec},
    covs_i::Union{Nothing,IndividualCovariates},
    has_covariate_model::Bool,
)
    if iov !== nothing
        return true
    end
    if has_covariate_model && covs_i !== nothing && covs_i.time_varying !== nothing
        return true
    end
    return false
end

function simulate_population(
    pop::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    pd_spec::Union{Nothing,PDSpec}=nothing,
)
    base = pop.base_model_spec

    # Determine population size and RNG setup for IIV
    if pop.iiv === nothing
        if pop.covariate_model !== nothing
            n = length(pop.covariates)
            if n < 1
                error("covariates must be non-empty when covariate_model is provided")
            end
        else
            n = 1
        end

        seed = UInt64(0)
        iiv_kind = "none"
        omegas = Dict{Symbol,Float64}()
        rng = StableRNG(0)
    else
        validate(pop.iiv)
        n = pop.iiv.n
        seed = pop.iiv.seed
        iiv_kind = string(typeof(pop.iiv.kind))
        omegas = pop.iiv.omegas
        rng = StableRNG(seed)
    end

    if pop.iov !== nothing
        validate(pop.iov)
    end

    # Covariates length checks
    if pop.covariate_model !== nothing
        if length(pop.covariates) != n
            error("covariates must be length n when covariate_model is provided")
        end
    else
        if !isempty(pop.covariates) && length(pop.covariates) != n
            error("covariates must be empty or length n")
        end
    end

    individuals = Vector{SimResult}(undef, n)
    realized_params = Vector{Dict{Symbol,Float64}}(undef, n)

    has_covariate_model = pop.covariate_model !== nothing

    for i in 1:n
        # Get individual covariates if available
        covs_i = if !isempty(pop.covariates)
            pop.covariates[i]
        else
            nothing
        end

        # Sample etas once per individual (constant across segments)
        # Use _sample_etas_from_iiv to support both diagonal and correlated sampling
        etas = if pop.iiv === nothing
            Dict{Symbol,Float64}()
        else
            _sample_etas_from_iiv(pop.iiv, base.params, rng)
        end

        # Precompute IOV occasions and kappas for this individual
        occ_starts = if pop.iov === nothing
            Float64[]
        else
            derive_occasions(base.doses, grid, pop.iov.occasion_def)
        end
        kappas = if pop.iov === nothing
            Vector{Dict{Symbol,Float64}}()
        else
            sample_iov_kappas(pop.iov.pis, length(occ_starts), pop.iov.seed + UInt64(i))
        end

        # Check if we need segmentation
        needs_seg = _needs_segmentation(pop.iov, covs_i, has_covariate_model)

        if !needs_seg
            # Simple path: no segmentation needed
            # Apply static covariates if any
            params_i = base.params
            snap = Dict{Symbol,Float64}()

            if has_covariate_model && covs_i !== nothing
                params_i, snap = apply_covariates(base.params, pop.covariate_model, covs_i)
            else
                T = typeof(base.params)
                for f in fieldnames(T)
                    snap[f] = Float64(getfield(base.params, f))
                end
            end

            # Apply IIV
            if pop.iiv !== nothing
                params_i, vals_iiv = _apply_etas_log_normal(params_i, etas)
                realized = copy(snap)
                for (k, v) in vals_iiv
                    realized[k] = v
                end
                realized_params[i] = realized
            else
                realized_params[i] = snap
            end

            spec_i = ModelSpec(base.kind, base.name * "_i$(i)", params_i, base.doses)

            if pd_spec === nothing
                individuals[i] = simulate(spec_i, grid, solver)
            else
                individuals[i] = simulate_pkpd_coupled(spec_i, pd_spec, grid, solver)
            end
        else
            # Segmentation path: IOV and/or time-varying covariates
            segment_starts = _build_segment_starts(
                grid, base.doses, pop.iov, covs_i, has_covariate_model
            )

            params_per_segment = Vector{typeof(base.params)}(undef, length(segment_starts))
            final_realized = Dict{Symbol,Float64}()

            for (si, tseg) in enumerate(segment_starts)
                # 1) Apply covariates at this time point
                params_cov = base.params
                snap = Dict{Symbol,Float64}()

                if has_covariate_model && covs_i !== nothing
                    params_cov, snap = apply_covariates_at_time(
                        base.params, pop.covariate_model, covs_i, tseg
                    )
                else
                    T = typeof(base.params)
                    for f in fieldnames(T)
                        snap[f] = Float64(getfield(base.params, f))
                    end
                end

                # 2) Apply IIV (etas are constant across segments)
                params_iiv, vals_iiv = _apply_etas_log_normal(params_cov, etas)

                # 3) Apply IOV if applicable
                params_final = params_iiv
                if pop.iov !== nothing && !isempty(occ_starts)
                    oi = occasion_index_at_time(occ_starts, tseg)
                    if oi >= 1 && oi <= length(kappas)
                        params_final, _ = apply_iov(params_iiv, kappas[oi])
                    end
                end

                params_per_segment[si] = params_final

                # Store final realized params from first segment for metadata
                if si == 1
                    final_realized = vals_iiv
                end
            end

            realized_params[i] = final_realized

            spec_i = ModelSpec(base.kind, base.name * "_i$(i)", base.params, base.doses)

            if pd_spec === nothing
                individuals[i] = simulate_segmented_pk(
                    spec_i, grid, solver, segment_starts, params_per_segment
                )
            else
                individuals[i] = simulate_segmented_pkpd_coupled(
                    spec_i, pd_spec, grid, solver, segment_starts, params_per_segment
                )
            end
        end
    end

    summaries = Dict{Symbol,PopulationSummary}()
    if haskey(individuals[1].observations, :conc)
        summaries[:conc] = compute_population_summary(
            individuals, :conc; probs=[0.05, 0.95]
        )
    end

    # Add PD observation summary if pd_spec is provided
    if pd_spec !== nothing
        pd_obs = pd_spec.output_observation
        if haskey(individuals[1].observations, pd_obs)
            summaries[pd_obs] = compute_population_summary(
                individuals, pd_obs; probs=[0.05, 0.95]
            )
        end
    end

    metadata = Dict{String,Any}(
        "n" => n,
        "seed" => seed,
        "iiv_kind" => iiv_kind,
        "omegas" => Dict(String(k) => v for (k, v) in omegas),
        "engine_version" => "0.1.1",
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
        "covariate_model" =>
            pop.covariate_model === nothing ? "none" : pop.covariate_model.name,
        "iov_kind" => pop.iov === nothing ? "none" : string(typeof(pop.iov.kind)),
        "iov_seed" => pop.iov === nothing ? 0 : Int(pop.iov.seed),
        "iov_pis" => if pop.iov === nothing
            Dict{String,Any}()
        else
            Dict(String(k) => v for (k, v) in pop.iov.pis)
        end,
        "occasion_def_mode" =>
            pop.iov === nothing ? "none" : String(pop.iov.occasion_def.mode),
    )

    return PopulationResult(individuals, realized_params, summaries, metadata)
end
