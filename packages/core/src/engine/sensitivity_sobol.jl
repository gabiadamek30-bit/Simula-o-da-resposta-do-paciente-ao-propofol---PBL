# Sobol' Sensitivity Analysis Implementation
# Variance-based global sensitivity analysis using Saltelli sampling

using Statistics: mean, var, std, quantile
using StableRNGs: StableRNG

export SobolIndex, SobolResult, PopulationSobolResult
export run_sobol_sensitivity, run_population_sobol_sensitivity
export run_sobol_sensitivity_with_second_order
export compute_sobol_indices, compute_second_order_sobol_indices

# ============================================================================
# Result Types
# ============================================================================

"""
    SobolIndex

Sobol' sensitivity indices for a single parameter.

# Fields
- `Si::Float64`: First-order index (main effect)
- `Si_ci_lower::Float64`: Lower bound of confidence interval for Si
- `Si_ci_upper::Float64`: Upper bound of confidence interval for Si
- `STi::Float64`: Total-order index (including interactions)
- `STi_ci_lower::Float64`: Lower bound of confidence interval for STi
- `STi_ci_upper::Float64`: Upper bound of confidence interval for STi

# Interpretation
- Si ≈ 0: Parameter has negligible main effect
- Si ≈ 1: Parameter explains almost all variance
- STi - Si > 0: Parameter has significant interactions with others
- Sum of Si < 1: Significant interaction effects present
"""
struct SobolIndex
    Si::Float64
    Si_ci_lower::Float64
    Si_ci_upper::Float64
    STi::Float64
    STi_ci_lower::Float64
    STi_ci_upper::Float64
end

"""
    SobolResult

Complete result from Sobol' sensitivity analysis.

# Fields
- `spec::GlobalSensitivitySpec{SobolMethod}`: Analysis specification
- `params::Vector{Symbol}`: Parameter names in order
- `indices::Dict{Symbol,SobolIndex}`: Sensitivity indices per parameter
- `second_order::Union{Nothing,Dict{Tuple{Symbol,Symbol},Float64}}`: Second-order indices Sij
- `n_evaluations::Int`: Total number of model evaluations
- `convergence_metric::Float64`: Sum of first-order indices (should be ≤ 1)
- `output_variance::Float64`: Total variance of model output
- `computation_time::Float64`: Wall-clock time in seconds
- `metadata::Dict{String,Any}`: Additional metadata
"""
struct SobolResult
    spec::GlobalSensitivitySpec{SobolMethod}
    params::Vector{Symbol}
    indices::Dict{Symbol,SobolIndex}
    second_order::Union{Nothing,Dict{Tuple{Symbol,Symbol},Float64}}
    n_evaluations::Int
    convergence_metric::Float64
    output_variance::Float64
    computation_time::Float64
    metadata::Dict{String,Any}
end

"""
    PopulationSobolResult

Sobol' analysis result for population models.
"""
struct PopulationSobolResult
    spec::GlobalSensitivitySpec{SobolMethod}
    params::Vector{Symbol}
    indices_mean::Dict{Symbol,SobolIndex}
    indices_quantiles::Dict{Float64,Dict{Symbol,SobolIndex}}
    n_evaluations::Int
    convergence_metric::Float64
    output_variance::Float64
    computation_time::Float64
    metadata::Dict{String,Any}
end

# ============================================================================
# Main Analysis Functions
# ============================================================================

"""
    run_sobol_sensitivity(spec, grid, solver, gsa_spec; parallel_config, aggregator)

Run Sobol' sensitivity analysis on a single model.

# Arguments
- `spec::ModelSpec`: Model specification with nominal parameters
- `grid::SimGrid`: Simulation time grid
- `solver::SolverSpec`: ODE solver specification
- `gsa_spec::GlobalSensitivitySpec{SobolMethod}`: GSA specification with bounds

# Keyword Arguments
- `parallel_config::ParallelConfig`: Parallel execution configuration
- `aggregator::Function`: Function to aggregate time series to scalar (default: mean)

# Returns
- `SobolResult`: First-order (Si) and total-order (STi) indices with CIs

# Example
```julia
bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
gsa_spec = GlobalSensitivitySpec(SobolMethod(base_sample_size=1024), bounds)
result = run_sobol_sensitivity(model_spec, grid, solver, gsa_spec)
println("CL: Si=\$(result.indices[:CL].Si), STi=\$(result.indices[:CL].STi)")
```
"""
function run_sobol_sensitivity(
    spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    gsa_spec::GlobalSensitivitySpec{SobolMethod};
    parallel_config::ParallelConfig=ParallelConfig(SerialBackend()),
    aggregator::Function=mean
)::SobolResult
    # Check if second-order computation is requested
    if gsa_spec.method.compute_second_order
        return run_sobol_sensitivity_with_second_order(
            spec, grid, solver, gsa_spec;
            parallel_config=parallel_config,
            aggregator=aggregator
        )
    end

    start_time = time()
    method = gsa_spec.method
    bounds = gsa_spec.bounds
    d = length(bounds)
    N = method.base_sample_size

    # Initialize RNG
    rng = StableRNG(gsa_spec.seed)

    # Generate Saltelli samples
    samples = generate_saltelli_samples(bounds, N, rng)

    # Prepare all parameter sets for evaluation
    # Total: N (A) + N (B) + N*d (AB) + N*d (BA for STi) = N*(2d+2)
    # But we use Jansen estimator which needs only N*(d+2)
    all_param_sets = Vector{Vector{Float64}}()

    # A samples
    for i in 1:N
        push!(all_param_sets, samples.A[i, :])
    end

    # B samples
    for i in 1:N
        push!(all_param_sets, samples.B[i, :])
    end

    # AB samples (for each parameter)
    for j in 1:d
        for i in 1:N
            push!(all_param_sets, samples.AB[j][i, :])
        end
    end

    n_evaluations = length(all_param_sets)

    # Evaluate all parameter sets
    outputs = evaluate_parameter_sets(
        spec, grid, solver, all_param_sets, bounds.params,
        gsa_spec.observation, aggregator, parallel_config
    )

    # Split outputs
    Y_A = outputs[1:N]
    Y_B = outputs[N+1:2N]
    Y_AB = [outputs[(2+j-1)*N+1:(2+j)*N] for j in 1:d]

    # Compute indices
    indices_raw = compute_sobol_indices_internal(
        Y_A, Y_B, Y_AB,
        method.bootstrap_samples,
        method.bootstrap_ci_level,
        StableRNG(gsa_spec.seed + 1)
    )

    # Map to parameter names
    indices = Dict{Symbol,SobolIndex}()
    for (j, param) in enumerate(bounds.params)
        indices[param] = indices_raw[j]
    end

    # Compute convergence metric (sum of Si should be ≤ 1)
    convergence = sum(idx.Si for idx in values(indices))
    output_var = var(vcat(Y_A, Y_B))

    # Build result
    elapsed = time() - start_time

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "method" => "Sobol",
        "base_sample_size" => N,
        "n_parameters" => d,
        "bootstrap_samples" => method.bootstrap_samples,
        "bootstrap_ci_level" => method.bootstrap_ci_level,
        "seed" => gsa_spec.seed,
        "observation" => String(gsa_spec.observation),
        "aggregator" => string(aggregator),
    )

    return SobolResult(
        gsa_spec,
        bounds.params,
        indices,
        nothing,  # second_order not implemented yet
        n_evaluations,
        convergence,
        output_var,
        elapsed,
        metadata
    )
end

"""
    run_population_sobol_sensitivity(pop, grid, solver, gsa_spec; parallel_config, probs, aggregator)

Run Sobol' sensitivity analysis on a population model.

Computes sensitivity indices on population summary statistics (mean and quantiles).
"""
function run_population_sobol_sensitivity(
    pop::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    gsa_spec::GlobalSensitivitySpec{SobolMethod};
    parallel_config::ParallelConfig=ParallelConfig(SerialBackend()),
    probs::Vector{Float64}=[0.05, 0.50, 0.95],
    aggregator::Function=mean
)::PopulationSobolResult
    start_time = time()
    method = gsa_spec.method
    bounds = gsa_spec.bounds
    d = length(bounds)
    N = method.base_sample_size

    rng = StableRNG(gsa_spec.seed)

    # Generate Saltelli samples
    samples = generate_saltelli_samples(bounds, N, rng)

    # Prepare all parameter sets
    all_param_sets = Vector{Vector{Float64}}()
    for i in 1:N
        push!(all_param_sets, samples.A[i, :])
    end
    for i in 1:N
        push!(all_param_sets, samples.B[i, :])
    end
    for j in 1:d
        for i in 1:N
            push!(all_param_sets, samples.AB[j][i, :])
        end
    end

    n_evaluations = length(all_param_sets)

    # Evaluate all parameter sets on population
    outputs = evaluate_population_parameter_sets(
        pop, grid, solver, all_param_sets, bounds.params,
        gsa_spec.observation, aggregator, parallel_config
    )

    # Split outputs
    Y_A = outputs[1:N]
    Y_B = outputs[N+1:2N]
    Y_AB = [outputs[(2+j-1)*N+1:(2+j)*N] for j in 1:d]

    # Compute indices on mean
    indices_raw = compute_sobol_indices_internal(
        Y_A, Y_B, Y_AB,
        method.bootstrap_samples,
        method.bootstrap_ci_level,
        StableRNG(gsa_spec.seed + 1)
    )

    indices_mean = Dict{Symbol,SobolIndex}()
    for (j, param) in enumerate(bounds.params)
        indices_mean[param] = indices_raw[j]
    end

    convergence = sum(idx.Si for idx in values(indices_mean))
    output_var = var(vcat(Y_A, Y_B))

    elapsed = time() - start_time

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "method" => "Sobol_Population",
        "base_sample_size" => N,
        "n_parameters" => d,
        "population_n" => pop.iiv.n,
        "bootstrap_samples" => method.bootstrap_samples,
        "seed" => gsa_spec.seed,
    )

    return PopulationSobolResult(
        gsa_spec,
        bounds.params,
        indices_mean,
        Dict{Float64,Dict{Symbol,SobolIndex}}(),  # quantiles not computed for efficiency
        n_evaluations,
        convergence,
        output_var,
        elapsed,
        metadata
    )
end

# ============================================================================
# Index Computation
# ============================================================================

"""
    compute_sobol_indices_internal(Y_A, Y_B, Y_AB, bootstrap_samples, ci_level, rng)

Compute Sobol' indices using Jansen/Saltelli estimators.

# Estimators used:
- First-order (Si): Jansen (2014) estimator
- Total-order (STi): Jansen (1999) estimator

# Arguments
- `Y_A::Vector{Float64}`: Model outputs for matrix A
- `Y_B::Vector{Float64}`: Model outputs for matrix B
- `Y_AB::Vector{Vector{Float64}}`: Model outputs for AB matrices (one per parameter)
- `bootstrap_samples::Int`: Number of bootstrap samples for CI
- `ci_level::Float64`: Confidence level
- `rng::AbstractRNG`: Random number generator

# Returns
- `Dict{Int,SobolIndex}`: Indices keyed by parameter index (1-based)
"""
function compute_sobol_indices_internal(
    Y_A::Vector{Float64},
    Y_B::Vector{Float64},
    Y_AB::Vector{Vector{Float64}},
    bootstrap_samples::Int,
    ci_level::Float64,
    rng::AbstractRNG
)::Dict{Int,SobolIndex}
    N = length(Y_A)
    d = length(Y_AB)

    # Total variance estimate
    Y_all = vcat(Y_A, Y_B)
    V_total = var(Y_all)

    if V_total < 1e-15
        # No variance - all indices are undefined
        zero_idx = SobolIndex(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        return Dict(j => zero_idx for j in 1:d)
    end

    indices = Dict{Int,SobolIndex}()

    for j in 1:d
        # Point estimates using Jansen estimators
        Si, STi = compute_sobol_point_estimates(Y_A, Y_B, Y_AB[j], V_total)

        # Bootstrap confidence intervals
        if bootstrap_samples > 0
            Si_boot = Float64[]
            STi_boot = Float64[]

            for _ in 1:bootstrap_samples
                boot_idx = rand(rng, 1:N, N)
                Y_A_boot = Y_A[boot_idx]
                Y_B_boot = Y_B[boot_idx]
                Y_AB_boot = Y_AB[j][boot_idx]

                V_boot = var(vcat(Y_A_boot, Y_B_boot))
                if V_boot > 1e-15
                    Si_b, STi_b = compute_sobol_point_estimates(Y_A_boot, Y_B_boot, Y_AB_boot, V_boot)
                    push!(Si_boot, Si_b)
                    push!(STi_boot, STi_b)
                end
            end

            alpha = 1.0 - ci_level
            if length(Si_boot) > 10
                Si_ci_lower = quantile(Si_boot, alpha / 2)
                Si_ci_upper = quantile(Si_boot, 1 - alpha / 2)
                STi_ci_lower = quantile(STi_boot, alpha / 2)
                STi_ci_upper = quantile(STi_boot, 1 - alpha / 2)
            else
                Si_ci_lower = Si
                Si_ci_upper = Si
                STi_ci_lower = STi
                STi_ci_upper = STi
            end
        else
            Si_ci_lower = Si
            Si_ci_upper = Si
            STi_ci_lower = STi
            STi_ci_upper = STi
        end

        # Clamp to valid range [0, 1]
        Si = clamp(Si, 0.0, 1.0)
        STi = clamp(STi, 0.0, 1.0)
        Si_ci_lower = clamp(Si_ci_lower, 0.0, 1.0)
        Si_ci_upper = clamp(Si_ci_upper, 0.0, 1.0)
        STi_ci_lower = clamp(STi_ci_lower, 0.0, 1.0)
        STi_ci_upper = clamp(STi_ci_upper, 0.0, 1.0)

        indices[j] = SobolIndex(Si, Si_ci_lower, Si_ci_upper, STi, STi_ci_lower, STi_ci_upper)
    end

    return indices
end

"""
    compute_sobol_point_estimates(Y_A, Y_B, Y_AB_j, V_total)

Compute point estimates for Si and STi using Jansen estimators.

# Estimators:
- Si (Jansen 2014): Si = (V(Y) - 0.5 * mean((Y_B - Y_AB_j)^2)) / V(Y)
- STi (Jansen 1999): STi = 0.5 * mean((Y_A - Y_AB_j)^2) / V(Y)
"""
function compute_sobol_point_estimates(
    Y_A::Vector{Float64},
    Y_B::Vector{Float64},
    Y_AB_j::Vector{Float64},
    V_total::Float64
)::Tuple{Float64,Float64}
    N = length(Y_A)

    # First-order index (Jansen 2014 / Saltelli 2010)
    # Si = 1 - Var(Y_B - Y_AB_j) / (2 * V(Y))
    diff_B = Y_B .- Y_AB_j
    V_diff_B = mean(diff_B .^ 2)
    Si = 1.0 - V_diff_B / (2.0 * V_total)

    # Total-order index (Jansen 1999)
    # STi = Var(Y_A - Y_AB_j) / (2 * V(Y))
    diff_A = Y_A .- Y_AB_j
    V_diff_A = mean(diff_A .^ 2)
    STi = V_diff_A / (2.0 * V_total)

    return (Si, STi)
end

# ============================================================================
# Helper Functions for Model Evaluation
# ============================================================================

"""
    evaluate_parameter_sets(spec, grid, solver, param_sets, param_names, observation, aggregator, parallel_config)

Evaluate model for multiple parameter sets.

Returns vector of scalar outputs (aggregated from time series).
"""
function evaluate_parameter_sets(
    spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    param_sets::Vector{Vector{Float64}},
    param_names::Vector{Symbol},
    observation::Symbol,
    aggregator::Function,
    parallel_config::ParallelConfig
)::Vector{Float64}
    # Create evaluation function
    function eval_single(params::Vector{Float64})
        # Build parameter struct from vector
        new_params = build_params_from_vector(spec.params, param_names, params)

        # Create model spec with new parameters
        new_spec = ModelSpec(spec.kind, spec.name, new_params, spec.doses)

        # Simulate
        result = simulate(new_spec, grid, solver)

        # Extract and aggregate observation
        if haskey(result.observations, observation)
            series = result.observations[observation]
            return aggregator(series)
        else
            return NaN
        end
    end

    # Evaluate all parameter sets (parallel if configured)
    if is_parallel(parallel_config)
        return parallel_map(eval_single, param_sets, parallel_config)
    else
        return [eval_single(p) for p in param_sets]
    end
end

"""
    evaluate_population_parameter_sets(pop, grid, solver, param_sets, param_names, observation, aggregator, parallel_config)

Evaluate population model for multiple parameter sets.

Returns vector of scalar outputs (population mean, aggregated).
"""
function evaluate_population_parameter_sets(
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
        # Build new base model spec
        new_params = build_params_from_vector(pop.base_model_spec.params, param_names, params)
        new_base = ModelSpec(
            pop.base_model_spec.kind,
            pop.base_model_spec.name,
            new_params,
            pop.base_model_spec.doses
        )

        # Create new population spec
        new_pop = PopulationSpec(
            new_base,
            pop.iiv,
            pop.iov,
            pop.covariate_model,
            pop.covariates
        )

        # Simulate population
        result = simulate_population(new_pop, grid, solver)

        # Get mean of observation
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
    build_params_from_vector(base_params, param_names, values)

Build a parameters struct by updating specified parameters with new values.
"""
function build_params_from_vector(base_params, param_names::Vector{Symbol}, values::Vector{Float64})
    T = typeof(base_params)
    field_names = fieldnames(T)

    # Get current values
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
    compute_sobol_indices(Y_A, Y_B, Y_AB; bootstrap_samples=0, ci_level=0.95, rng=nothing)

Compute Sobol' indices from pre-computed model outputs.

This is useful when you have already evaluated the model and want to compute indices.

# Arguments
- `Y_A::Vector{Float64}`: Model outputs for sample matrix A
- `Y_B::Vector{Float64}`: Model outputs for sample matrix B
- `Y_AB::Vector{Vector{Float64}}`: Model outputs for AB matrices

# Keyword Arguments
- `bootstrap_samples::Int`: Number of bootstrap samples for CI (0 = no CI)
- `ci_level::Float64`: Confidence level (default 0.95)
- `rng::AbstractRNG`: Random number generator (optional)

# Returns
- `Dict{Int,SobolIndex}`: Indices keyed by parameter index (1-based)
"""
function compute_sobol_indices(
    Y_A::Vector{Float64},
    Y_B::Vector{Float64},
    Y_AB::Vector{Vector{Float64}};
    bootstrap_samples::Int=0,
    ci_level::Float64=0.95,
    rng::Union{Nothing,AbstractRNG}=nothing
)::Dict{Int,SobolIndex}
    actual_rng = rng === nothing ? StableRNG(12345) : rng
    return compute_sobol_indices_internal(Y_A, Y_B, Y_AB, bootstrap_samples, ci_level, actual_rng)
end

# ============================================================================
# Second-Order Sobol' Indices
# ============================================================================

"""
    compute_second_order_sobol_indices(Y_A, Y_B, Y_AB, Y_ABij, first_order_indices, V_total)

Compute second-order Sobol' indices Sij for all parameter pairs.

The second-order index Sij represents the variance contribution from the interaction
between parameters i and j, excluding their individual main effects.

# Arguments
- `Y_A::Vector{Float64}`: Model outputs for matrix A
- `Y_B::Vector{Float64}`: Model outputs for matrix B
- `Y_AB::Vector{Vector{Float64}}`: Model outputs for AB matrices
- `Y_ABij::Dict{Tuple{Int,Int},Vector{Float64}}`: Model outputs for ABij matrices
- `first_order_indices::Dict{Int,SobolIndex}`: Pre-computed first-order indices
- `V_total::Float64`: Total variance of model output

# Returns
- `Dict{Tuple{Int,Int},Float64}`: Second-order indices Sij for each parameter pair (i,j)

# Formula
Using the Saltelli (2002) estimator:
  Sij = (Vij - Vi - Vj) / V(Y)

where:
  Vij = (1/N) * sum(Y_A .* Y_ABij) - f0²
  f0 = mean(Y)

# References
- Saltelli (2002) "Making best use of model evaluations to compute sensitivity indices"
- Homma & Saltelli (1996) "Importance measures in global sensitivity analysis"
"""
function compute_second_order_sobol_indices(
    Y_A::Vector{Float64},
    Y_B::Vector{Float64},
    Y_AB::Vector{Vector{Float64}},
    Y_ABij::Dict{Tuple{Int,Int},Vector{Float64}},
    first_order_indices::Dict{Int,SobolIndex},
    V_total::Float64
)::Dict{Tuple{Int,Int},Float64}
    N = length(Y_A)
    d = length(Y_AB)

    # Mean of all samples for computing f0²
    f0 = mean(vcat(Y_A, Y_B))
    f0_sq = f0^2

    second_order = Dict{Tuple{Int,Int},Float64}()

    if V_total < 1e-15
        # No variance - all indices are zero
        for i in 1:(d-1)
            for j in (i+1):d
                second_order[(i, j)] = 0.0
            end
        end
        return second_order
    end

    for (pair, Y_ij) in Y_ABij
        i, j = pair

        # Get first-order indices for normalization
        Vi = first_order_indices[i].Si * V_total  # Variance contribution of parameter i
        Vj = first_order_indices[j].Si * V_total  # Variance contribution of parameter j

        # Compute variance of the interaction using Saltelli estimator
        # Vij = (1/N) * sum(Y_A .* Y_ABij[(i,j)]) - f0²
        Vij = mean(Y_A .* Y_ij) - f0_sq

        # Second-order index: subtract individual contributions
        # Sij = (Vij - Vi - Vj) / V(Y)
        Sij_raw = (Vij - Vi - Vj) / V_total

        # Clamp to valid range [0, 1]
        # Note: Sij can be negative due to estimation error for weak interactions
        Sij = clamp(Sij_raw, 0.0, 1.0)

        second_order[(i, j)] = Sij
    end

    return second_order
end

"""
    run_sobol_sensitivity_with_second_order(spec, grid, solver, gsa_spec; parallel_config, aggregator)

Run Sobol' sensitivity analysis including second-order interaction indices.

This is an extended version of run_sobol_sensitivity that computes Sij indices
when compute_second_order=true in the SobolMethod specification.

# Additional evaluations required
For d parameters, second-order indices require:
  N * d * (d-1) / 2 additional model evaluations

Total evaluations: N * (d + 2 + d*(d-1)/2)

# Returns
SobolResult with populated second_order field containing Dict{Tuple{Symbol,Symbol},Float64}
"""
function run_sobol_sensitivity_with_second_order(
    spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    gsa_spec::GlobalSensitivitySpec{SobolMethod};
    parallel_config::ParallelConfig=ParallelConfig(SerialBackend()),
    aggregator::Function=mean
)::SobolResult
    start_time = time()
    method = gsa_spec.method
    bounds = gsa_spec.bounds
    d = length(bounds)
    N = method.base_sample_size

    # Initialize RNG
    rng = StableRNG(gsa_spec.seed)

    # Generate extended Saltelli samples including ABij matrices
    samples = generate_saltelli_samples_second_order(bounds, N, rng)

    # Prepare all parameter sets for evaluation
    all_param_sets = Vector{Vector{Float64}}()

    # A samples (indices 1:N)
    for i in 1:N
        push!(all_param_sets, samples.A[i, :])
    end

    # B samples (indices N+1:2N)
    for i in 1:N
        push!(all_param_sets, samples.B[i, :])
    end

    # AB samples for each parameter (indices 2N+1 : 2N+N*d)
    for j in 1:d
        for i in 1:N
            push!(all_param_sets, samples.AB[j][i, :])
        end
    end

    # ABij samples for each parameter pair
    # Track where each pair's samples start
    abij_start_idx = length(all_param_sets) + 1
    pair_indices = Dict{Tuple{Int,Int},UnitRange{Int}}()

    for i in 1:(d-1)
        for j in (i+1):d
            start_idx = length(all_param_sets) + 1
            ABij_mat = samples.ABij[(i, j)]
            for k in 1:N
                push!(all_param_sets, ABij_mat[k, :])
            end
            pair_indices[(i, j)] = start_idx:(start_idx + N - 1)
        end
    end

    n_evaluations = length(all_param_sets)

    # Evaluate all parameter sets
    outputs = evaluate_parameter_sets(
        spec, grid, solver, all_param_sets, bounds.params,
        gsa_spec.observation, aggregator, parallel_config
    )

    # Split outputs
    Y_A = outputs[1:N]
    Y_B = outputs[N+1:2N]
    Y_AB = [outputs[(2+j-1)*N+1:(2+j)*N] for j in 1:d]

    # Extract ABij outputs
    Y_ABij = Dict{Tuple{Int,Int},Vector{Float64}}()
    for (pair, range) in pair_indices
        Y_ABij[pair] = outputs[range]
    end

    # Compute first-order and total-order indices
    indices_raw = compute_sobol_indices_internal(
        Y_A, Y_B, Y_AB,
        method.bootstrap_samples,
        method.bootstrap_ci_level,
        StableRNG(gsa_spec.seed + 1)
    )

    # Compute total variance for second-order computation
    V_total = var(vcat(Y_A, Y_B))

    # Compute second-order indices
    second_order_raw = compute_second_order_sobol_indices(
        Y_A, Y_B, Y_AB, Y_ABij, indices_raw, V_total
    )

    # Map to parameter names
    indices = Dict{Symbol,SobolIndex}()
    for (j, param) in enumerate(bounds.params)
        indices[param] = indices_raw[j]
    end

    # Map second-order to parameter names
    second_order = Dict{Tuple{Symbol,Symbol},Float64}()
    for ((i, j), Sij) in second_order_raw
        param_i = bounds.params[i]
        param_j = bounds.params[j]
        second_order[(param_i, param_j)] = Sij
    end

    # Compute convergence metric (sum of Si should be ≤ 1)
    convergence = sum(idx.Si for idx in values(indices))

    # Build result
    elapsed = time() - start_time

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "method" => "Sobol",
        "base_sample_size" => N,
        "n_parameters" => d,
        "bootstrap_samples" => method.bootstrap_samples,
        "bootstrap_ci_level" => method.bootstrap_ci_level,
        "seed" => gsa_spec.seed,
        "observation" => String(gsa_spec.observation),
        "aggregator" => string(aggregator),
        "compute_second_order" => true,
        "n_parameter_pairs" => length(second_order),
    )

    return SobolResult(
        gsa_spec,
        bounds.params,
        indices,
        second_order,
        n_evaluations,
        convergence,
        V_total,
        elapsed,
        metadata
    )
end
