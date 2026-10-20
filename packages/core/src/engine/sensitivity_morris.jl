# Morris Elementary Effects Sensitivity Analysis Implementation
# A computationally efficient screening method for identifying important parameters

using Statistics: mean, std
using StableRNGs: StableRNG

export MorrisIndex, MorrisResult, PopulationMorrisResult
export run_morris_sensitivity, run_population_morris_sensitivity
export compute_morris_indices
export rank_parameters, identify_important_parameters

# ============================================================================
# Result Types
# ============================================================================

"""
    MorrisIndex

Morris Elementary Effects indices for a single parameter.

# Fields
- `mu::Float64`: Mean elementary effect (signed)
- `mu_star::Float64`: Mean of absolute elementary effects (importance measure)
- `sigma::Float64`: Standard deviation of elementary effects

# Interpretation
- μ* large: Parameter is important (affects output significantly)
- μ ≈ 0 but μ* large: Parameter has non-monotonic effect or interactions
- σ large: Parameter has nonlinear effects or interactions with others
- σ/μ* ratio: Indicates degree of nonlinearity/interactions

# Classification (Campolongo et al. 2007)
- μ* < threshold: Negligible effect
- σ/μ* ≈ 0: Mostly linear, additive effect
- σ/μ* ≈ 1: Nonlinear or interacting
"""
struct MorrisIndex
    mu::Float64
    mu_star::Float64
    sigma::Float64
end

"""
    MorrisResult

Complete result from Morris Elementary Effects analysis.

# Fields
- `spec::GlobalSensitivitySpec{MorrisMethod}`: Analysis specification
- `params::Vector{Symbol}`: Parameter names in order
- `indices::Dict{Symbol,MorrisIndex}`: Morris indices per parameter
- `elementary_effects::Dict{Symbol,Vector{Float64}}`: Raw EEs for each parameter
- `n_evaluations::Int`: Total number of model evaluations
- `computation_time::Float64`: Wall-clock time in seconds
- `metadata::Dict{String,Any}`: Additional metadata

# Usage
Morris screening is typically used before expensive variance-based methods (Sobol').
Parameters with low μ* can be fixed, reducing dimensionality for subsequent analysis.
"""
struct MorrisResult
    spec::GlobalSensitivitySpec{MorrisMethod}
    params::Vector{Symbol}
    indices::Dict{Symbol,MorrisIndex}
    elementary_effects::Dict{Symbol,Vector{Float64}}
    n_evaluations::Int
    computation_time::Float64
    metadata::Dict{String,Any}
end

"""
    PopulationMorrisResult

Morris analysis result for population models.
"""
struct PopulationMorrisResult
    spec::GlobalSensitivitySpec{MorrisMethod}
    params::Vector{Symbol}
    indices::Dict{Symbol,MorrisIndex}
    elementary_effects::Dict{Symbol,Vector{Float64}}
    n_evaluations::Int
    computation_time::Float64
    metadata::Dict{String,Any}
end

# ============================================================================
# Main Analysis Functions
# ============================================================================

"""
    run_morris_sensitivity(spec, grid, solver, gsa_spec; parallel_config, aggregator)

Run Morris Elementary Effects screening analysis on a single model.

# Arguments
- `spec::ModelSpec`: Model specification with nominal parameters
- `grid::SimGrid`: Simulation time grid
- `solver::SolverSpec`: ODE solver specification
- `gsa_spec::GlobalSensitivitySpec{MorrisMethod}`: GSA specification with bounds

# Keyword Arguments
- `parallel_config::ParallelConfig`: Parallel execution configuration
- `aggregator::Function`: Function to aggregate time series to scalar (default: mean)

# Returns
- `MorrisResult`: μ, μ*, and σ for each parameter

# Example
```julia
bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
gsa_spec = GlobalSensitivitySpec(MorrisMethod(n_trajectories=20), bounds)
result = run_morris_sensitivity(model_spec, grid, solver, gsa_spec)
println("CL importance: μ*=\$(result.indices[:CL].mu_star)")
println("V importance:  μ*=\$(result.indices[:V].mu_star)")
```

# Computational Cost
Total evaluations = r × (d + 1) where:
- r = n_trajectories
- d = number of parameters
Much cheaper than Sobol' (which requires N × (d + 2))
"""
function run_morris_sensitivity(
    spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    gsa_spec::GlobalSensitivitySpec{MorrisMethod};
    parallel_config::ParallelConfig=ParallelConfig(SerialBackend()),
    aggregator::Function=mean
)::MorrisResult
    start_time = time()
    method = gsa_spec.method
    bounds = gsa_spec.bounds
    d = length(bounds)
    r = method.n_trajectories
    delta = method.delta

    # Initialize RNG
    rng = StableRNG(gsa_spec.seed)

    # Generate Morris trajectories
    trajectories = generate_morris_trajectories(bounds, r, method.n_levels, delta, rng)

    # Flatten all trajectory points for batch evaluation
    all_points = Vector{Vector{Float64}}()
    point_to_traj_step = Vector{Tuple{Int,Int}}()  # (trajectory_idx, step_idx)

    for (traj_idx, trajectory) in enumerate(trajectories)
        n_points = size(trajectory, 1)  # d+1 points per trajectory
        for step_idx in 1:n_points
            push!(all_points, trajectory[step_idx, :])
            push!(point_to_traj_step, (traj_idx, step_idx))
        end
    end

    n_evaluations = length(all_points)

    # Evaluate all points
    outputs = evaluate_morris_parameter_sets(
        spec, grid, solver, all_points, bounds.params,
        gsa_spec.observation, aggregator, parallel_config
    )

    # Reorganize outputs by trajectory
    traj_outputs = [Vector{Float64}() for _ in 1:r]
    for (i, (traj_idx, _)) in enumerate(point_to_traj_step)
        push!(traj_outputs[traj_idx], outputs[i])
    end

    # Compute elementary effects
    elementary_effects = compute_elementary_effects(
        trajectories, traj_outputs, bounds, delta
    )

    # Compute indices
    indices = Dict{Symbol,MorrisIndex}()
    for param in bounds.params
        ees = elementary_effects[param]
        if length(ees) > 0
            mu = mean(ees)
            mu_star = mean(abs.(ees))
            sigma = length(ees) > 1 ? std(ees) : 0.0
            indices[param] = MorrisIndex(mu, mu_star, sigma)
        else
            indices[param] = MorrisIndex(0.0, 0.0, 0.0)
        end
    end

    elapsed = time() - start_time

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "method" => "Morris",
        "n_trajectories" => r,
        "n_levels" => method.n_levels,
        "delta" => delta,
        "n_parameters" => d,
        "seed" => gsa_spec.seed,
        "observation" => String(gsa_spec.observation),
        "aggregator" => string(aggregator),
    )

    return MorrisResult(
        gsa_spec,
        bounds.params,
        indices,
        elementary_effects,
        n_evaluations,
        elapsed,
        metadata
    )
end

"""
    run_population_morris_sensitivity(pop, grid, solver, gsa_spec; parallel_config, aggregator)

Run Morris screening analysis on a population model.

Computes elementary effects on population mean output.
"""
function run_population_morris_sensitivity(
    pop::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    gsa_spec::GlobalSensitivitySpec{MorrisMethod};
    parallel_config::ParallelConfig=ParallelConfig(SerialBackend()),
    aggregator::Function=mean
)::PopulationMorrisResult
    start_time = time()
    method = gsa_spec.method
    bounds = gsa_spec.bounds
    d = length(bounds)
    r = method.n_trajectories
    delta = method.delta

    rng = StableRNG(gsa_spec.seed)

    # Generate trajectories
    trajectories = generate_morris_trajectories(bounds, r, method.n_levels, delta, rng)

    # Flatten all trajectory points
    all_points = Vector{Vector{Float64}}()
    point_to_traj_step = Vector{Tuple{Int,Int}}()

    for (traj_idx, trajectory) in enumerate(trajectories)
        n_points = size(trajectory, 1)
        for step_idx in 1:n_points
            push!(all_points, trajectory[step_idx, :])
            push!(point_to_traj_step, (traj_idx, step_idx))
        end
    end

    n_evaluations = length(all_points)

    # Evaluate all points on population
    outputs = evaluate_morris_population_parameter_sets(
        pop, grid, solver, all_points, bounds.params,
        gsa_spec.observation, aggregator, parallel_config
    )

    # Reorganize outputs by trajectory
    traj_outputs = [Vector{Float64}() for _ in 1:r]
    for (i, (traj_idx, _)) in enumerate(point_to_traj_step)
        push!(traj_outputs[traj_idx], outputs[i])
    end

    # Compute elementary effects
    elementary_effects = compute_elementary_effects(
        trajectories, traj_outputs, bounds, delta
    )

    # Compute indices
    indices = Dict{Symbol,MorrisIndex}()
    for param in bounds.params
        ees = elementary_effects[param]
        if length(ees) > 0
            mu = mean(ees)
            mu_star = mean(abs.(ees))
            sigma = length(ees) > 1 ? std(ees) : 0.0
            indices[param] = MorrisIndex(mu, mu_star, sigma)
        else
            indices[param] = MorrisIndex(0.0, 0.0, 0.0)
        end
    end

    elapsed = time() - start_time

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "method" => "Morris_Population",
        "n_trajectories" => r,
        "n_levels" => method.n_levels,
        "delta" => delta,
        "n_parameters" => d,
        "population_n" => pop.iiv.n,
        "seed" => gsa_spec.seed,
    )

    return PopulationMorrisResult(
        gsa_spec,
        bounds.params,
        indices,
        elementary_effects,
        n_evaluations,
        elapsed,
        metadata
    )
end

# ============================================================================
# Elementary Effects Computation
# ============================================================================

"""
    compute_elementary_effects(trajectories, traj_outputs, bounds, delta)

Compute elementary effects from trajectory evaluations.

# Arguments
- `trajectories::Vector{Matrix{Float64}}`: Morris trajectories
- `traj_outputs::Vector{Vector{Float64}}`: Model outputs for each trajectory point
- `bounds::ParameterBounds`: Parameter bounds
- `delta::Float64`: Step size in normalized space

# Returns
- `Dict{Symbol,Vector{Float64}}`: Elementary effects for each parameter

# Algorithm
For each consecutive pair of points in a trajectory, identify which parameter
changed and compute:
  EE_i = (Y(x + Δe_i) - Y(x)) / Δ
where Δ is the normalized step size.
"""
function compute_elementary_effects(
    trajectories::Vector{Matrix{Float64}},
    traj_outputs::Vector{Vector{Float64}},
    bounds::ParameterBounds,
    delta::Float64
)::Dict{Symbol,Vector{Float64}}
    d = length(bounds)

    # Initialize storage for elementary effects
    effects = Dict{Symbol,Vector{Float64}}()
    for param in bounds.params
        effects[param] = Float64[]
    end

    for (traj_idx, trajectory) in enumerate(trajectories)
        outputs = traj_outputs[traj_idx]
        n_steps = size(trajectory, 1) - 1  # d steps

        for step in 1:n_steps
            # Find which parameter changed
            prev_point = trajectory[step, :]
            curr_point = trajectory[step + 1, :]

            param_idx = find_changed_parameter(prev_point, curr_point, bounds)

            if param_idx > 0 && param_idx <= d
                param = bounds.params[param_idx]

                # Compute elementary effect
                # The step size in parameter space
                range_j = bounds.upper[param_idx] - bounds.lower[param_idx]
                delta_actual = (curr_point[param_idx] - prev_point[param_idx]) / range_j

                if abs(delta_actual) > 1e-10
                    # EE = ΔY / Δx (normalized)
                    ee = (outputs[step + 1] - outputs[step]) / delta_actual
                    push!(effects[param], ee)
                end
            end
        end
    end

    return effects
end

"""
    find_changed_parameter(prev_point, curr_point, bounds)

Find which parameter changed between two trajectory points.

Returns the 1-based index of the changed parameter, or 0 if none found.
"""
function find_changed_parameter(
    prev_point::Vector{Float64},
    curr_point::Vector{Float64},
    bounds::ParameterBounds
)::Int
    d = length(bounds)

    max_change = 0.0
    changed_idx = 0

    for j in 1:d
        range_j = bounds.upper[j] - bounds.lower[j]
        normalized_diff = abs(curr_point[j] - prev_point[j]) / range_j

        if normalized_diff > max_change
            max_change = normalized_diff
            changed_idx = j
        end
    end

    # Only return if the change is significant
    if max_change > 0.01
        return changed_idx
    else
        return 0
    end
end

# ============================================================================
# Helper Functions for Model Evaluation
# ============================================================================

"""
    evaluate_morris_parameter_sets(spec, grid, solver, param_sets, param_names, observation, aggregator, parallel_config)

Evaluate model for Morris trajectory points.
"""
function evaluate_morris_parameter_sets(
    spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    param_sets::Vector{Vector{Float64}},
    param_names::Vector{Symbol},
    observation::Symbol,
    aggregator::Function,
    parallel_config::ParallelConfig
)::Vector{Float64}
    function eval_single(params::Vector{Float64})
        new_params = build_morris_params_from_vector(spec.params, param_names, params)
        new_spec = ModelSpec(spec.kind, spec.name, new_params, spec.doses)

        result = simulate(new_spec, grid, solver)

        if haskey(result.observations, observation)
            series = result.observations[observation]
            return aggregator(series)
        else
            return NaN
        end
    end

    if is_parallel(parallel_config)
        return parallel_map(eval_single, param_sets, parallel_config)
    else
        return [eval_single(p) for p in param_sets]
    end
end

"""
    evaluate_morris_population_parameter_sets(pop, grid, solver, param_sets, param_names, observation, aggregator, parallel_config)

Evaluate population model for Morris trajectory points.
"""
function evaluate_morris_population_parameter_sets(
    pop::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    param_sets::Vector{Vector{Float64}},
    param_names::Vector{Symbol},
    observation::Symbol,
    aggregator::Function,
    parallel_config::ParallelConfig
)::Vector{Float64}
    function eval_single(params::Vector{Float64})
        new_params = build_morris_params_from_vector(pop.base_model_spec.params, param_names, params)
        new_base = ModelSpec(
            pop.base_model_spec.kind,
            pop.base_model_spec.name,
            new_params,
            pop.base_model_spec.doses
        )

        new_pop = PopulationSpec(
            new_base,
            pop.iiv,
            pop.iov,
            pop.covariate_model,
            pop.covariates
        )

        result = simulate_population(new_pop, grid, solver)

        if haskey(result.summaries, observation)
            series_mean = result.summaries[observation].mean
            return aggregator(series_mean)
        else
            return NaN
        end
    end

    if is_parallel(parallel_config)
        return parallel_map(eval_single, param_sets, parallel_config)
    else
        return [eval_single(p) for p in param_sets]
    end
end

"""
    build_morris_params_from_vector(base_params, param_names, values)

Build a parameters struct by updating specified parameters with new values.
"""
function build_morris_params_from_vector(base_params, param_names::Vector{Symbol}, values::Vector{Float64})
    T = typeof(base_params)
    field_names = fieldnames(T)

    new_values = Float64[]
    for fn in field_names
        if fn in param_names
            idx = findfirst(==(fn), param_names)
            push!(new_values, values[idx])
        else
            push!(new_values, Float64(getfield(base_params, fn)))
        end
    end

    return T(new_values...)
end

# ============================================================================
# Public API for Index Computation
# ============================================================================

"""
    compute_morris_indices(elementary_effects)

Compute Morris indices from pre-computed elementary effects.

# Arguments
- `elementary_effects::Dict{Symbol,Vector{Float64}}`: Elementary effects for each parameter

# Returns
- `Dict{Symbol,MorrisIndex}`: Indices for each parameter
"""
function compute_morris_indices(
    elementary_effects::Dict{Symbol,Vector{Float64}}
)::Dict{Symbol,MorrisIndex}
    indices = Dict{Symbol,MorrisIndex}()

    for (param, ees) in elementary_effects
        if length(ees) > 0
            mu = mean(ees)
            mu_star = mean(abs.(ees))
            sigma = length(ees) > 1 ? std(ees) : 0.0
            indices[param] = MorrisIndex(mu, mu_star, sigma)
        else
            indices[param] = MorrisIndex(0.0, 0.0, 0.0)
        end
    end

    return indices
end

# ============================================================================
# Screening Utilities
# ============================================================================

"""
    rank_parameters(result::MorrisResult; by=:mu_star)

Rank parameters by importance metric.

# Arguments
- `result::MorrisResult`: Morris analysis result
- `by::Symbol`: Metric to rank by (:mu_star, :mu, :sigma)

# Returns
- `Vector{Tuple{Symbol,Float64}}`: Parameters sorted by metric (descending)
"""
function rank_parameters(result::MorrisResult; by::Symbol=:mu_star)::Vector{Tuple{Symbol,Float64}}
    rankings = Tuple{Symbol,Float64}[]

    for (param, idx) in result.indices
        val = if by == :mu_star
            idx.mu_star
        elseif by == :mu
            abs(idx.mu)
        elseif by == :sigma
            idx.sigma
        else
            error("Unknown ranking metric: $by")
        end
        push!(rankings, (param, val))
    end

    sort!(rankings, by=x -> x[2], rev=true)
    return rankings
end

"""
    identify_important_parameters(result::MorrisResult; threshold=0.1)

Identify parameters with μ* above a threshold (relative to max μ*).

# Arguments
- `result::MorrisResult`: Morris analysis result
- `threshold::Float64`: Relative threshold (0-1)

# Returns
- `Vector{Symbol}`: Important parameter names
"""
function identify_important_parameters(
    result::MorrisResult;
    threshold::Float64=0.1
)::Vector{Symbol}
    rankings = rank_parameters(result; by=:mu_star)

    if isempty(rankings)
        return Symbol[]
    end

    max_mu_star = rankings[1][2]
    if max_mu_star < 1e-10
        return Symbol[]
    end

    cutoff = threshold * max_mu_star
    return [param for (param, val) in rankings if val >= cutoff]
end
