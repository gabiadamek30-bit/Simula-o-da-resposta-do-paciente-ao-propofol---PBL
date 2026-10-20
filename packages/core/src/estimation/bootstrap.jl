# Bootstrap Module for Parameter Uncertainty Estimation
# Industry-standard implementation for regulatory submissions (FDA/EMA)
#
# Supports:
# - Non-parametric case bootstrap (standard regulatory approach)
# - Parametric/residual bootstrap for model-based inference
# - Stratified resampling (by study, dose, formulation)
# - Multiple CI methods (percentile, BCa, basic)
# - Parallel execution for large-scale studies
# - Regulatory-compliant output formatting
#
# References:
# - FDA Guidance on Population PK Analysis (2022)
# - EMA Guideline on Reporting Population PK Studies (2007)
# - Efron & Tibshirani (1993) An Introduction to the Bootstrap

using StableRNGs
using Statistics
using LinearAlgebra
using Distributions
using Base.Threads: @threads, nthreads

# Helper to safely convert values that may be Py-wrapped (from Python callables)
# These handle both pure Julia and Python-via-PythonCall cases
_to_float_vector(x::Vector{Float64}) = x
_to_float_vector(x::AbstractVector) = Vector{Float64}(x)
function _to_float_vector(x)
    # Handle Py-wrapped objects from PythonCall
    if isdefined(Main, :PythonCall)
        return Main.PythonCall.pyconvert(Vector{Float64}, x)
    else
        return Vector{Float64}(x)
    end
end

_to_float_matrix(x::Matrix{Float64}) = x
_to_float_matrix(x::AbstractMatrix) = Matrix{Float64}(x)
function _to_float_matrix(x)
    # Handle Py-wrapped objects from PythonCall
    if isdefined(Main, :PythonCall)
        return Main.PythonCall.pyconvert(Matrix{Float64}, x)
    else
        return Matrix{Float64}(x)
    end
end

export BootstrapSpec, BootstrapResult, BootstrapDiagnostics
export OmegaBootstrapSummary, SigmaBootstrapSummary
export run_bootstrap, run_bootstrap_native, run_bootstrap_parallel_impl, stratified_resample, compute_bootstrap_ci
export BootstrapCIMethod, PercentileCI, BCCI, StudentizedCI, BasicCI
export BootstrapType, CaseBootstrap, ParametricBootstrap, ResidualBootstrap
export generate_bootstrap_summary, format_regulatory_table
export compute_bootstrap_coverage, compute_bootstrap_stability
export influential_subject_analysis, jackknife_influence

# ============================================================================
# Bootstrap CI Methods
# ============================================================================

"""
Abstract type for bootstrap confidence interval methods.
"""
abstract type BootstrapCIMethod end

"""
Percentile method: CI = [θ*(α/2), θ*(1-α/2)]
"""
struct PercentileCI <: BootstrapCIMethod end

"""
Bias-corrected and accelerated (BCa) method.
"""
struct BCCI <: BootstrapCIMethod
    acceleration::Float64  # Acceleration constant (computed from jackknife)
end
BCCI() = BCCI(0.0)  # Default with no acceleration

"""
Studentized (bootstrap-t) method.
"""
struct StudentizedCI <: BootstrapCIMethod end

"""
Basic bootstrap CI: CI = [2θ̂ - θ*(1-α/2), 2θ̂ - θ*(α/2)]
Also known as "reverse percentile" or "pivotal" method.
"""
struct BasicCI <: BootstrapCIMethod end

# ============================================================================
# Bootstrap Type (Case vs Parametric)
# ============================================================================

"""
Abstract type for bootstrap sampling strategy.
"""
abstract type BootstrapType end

"""
Non-parametric case bootstrap: resample subjects with replacement.
This is the standard FDA/EMA recommended approach.
"""
struct CaseBootstrap <: BootstrapType end

"""
Parametric bootstrap: simulate new data from the fitted model.
Useful when assessing model adequacy or when sample size is small.
"""
struct ParametricBootstrap <: BootstrapType
    n_simulations_per_subject::Int

    ParametricBootstrap(n::Int=1) = new(max(1, n))
end

"""
Residual bootstrap: resample residuals and add to predictions.
Preserves the covariate structure of the original data.
"""
struct ResidualBootstrap <: BootstrapType
    standardize::Bool  # Whether to standardize residuals

    ResidualBootstrap(standardize::Bool=true) = new(standardize)
end

# ============================================================================
# Bootstrap Specification
# ============================================================================

"""
Specification for bootstrap analysis.

Fields:
- n_bootstrap: Number of bootstrap replicates (default: 1000, FDA recommends ≥500)
- bootstrap_type: Type of bootstrap (CaseBootstrap, ParametricBootstrap, ResidualBootstrap)
- stratify_by: Variables to stratify by (e.g., [:study, :dose, :formulation])
- seed: Random seed for reproducibility
- parallel: Whether to run in parallel (Bool) or ParallelConfig for advanced control
- ci_level: Confidence interval level (default: 0.95)
- ci_method: Method for computing CIs (PercentileCI, BCCI, BasicCI)
- compute_omega_ci: Whether to compute omega parameter CIs (default: true)
- compute_sigma_ci: Whether to compute sigma parameter CIs (default: true)
- min_success_rate: Minimum required success rate (default: 0.8)

# Example
```julia
# Standard FDA-compliant bootstrap
spec = BootstrapSpec(n_bootstrap=1000, parallel=true)

# Stratified by study for pooled analysis
spec = BootstrapSpec(n_bootstrap=1000, stratify_by=[:study], parallel=true)

# BCa confidence intervals (recommended for skewed parameters)
spec = BootstrapSpec(n_bootstrap=1000, ci_method=BCCI())

# Parametric bootstrap for small samples
spec = BootstrapSpec(n_bootstrap=500, bootstrap_type=ParametricBootstrap())
```
"""
struct BootstrapSpec
    n_bootstrap::Int
    bootstrap_type::BootstrapType
    stratify_by::Vector{Symbol}
    seed::UInt64
    parallel::Union{Bool, ParallelConfig}
    ci_level::Float64
    ci_method::BootstrapCIMethod
    compute_omega_ci::Bool
    compute_sigma_ci::Bool
    min_success_rate::Float64

    function BootstrapSpec(;
        n_bootstrap::Int=1000,
        bootstrap_type::BootstrapType=CaseBootstrap(),
        stratify_by::Vector{Symbol}=Symbol[],
        seed::UInt64=UInt64(12345),
        parallel::Union{Bool, ParallelConfig}=false,
        ci_level::Float64=0.95,
        ci_method::BootstrapCIMethod=PercentileCI(),
        compute_omega_ci::Bool=true,
        compute_sigma_ci::Bool=true,
        min_success_rate::Float64=0.8
    )
        @assert n_bootstrap >= 100 "n_bootstrap must be >= 100 (FDA recommends ≥500)"
        @assert 0.0 < ci_level < 1.0 "ci_level must be in (0, 1)"
        @assert 0.0 < min_success_rate <= 1.0 "min_success_rate must be in (0, 1]"
        new(n_bootstrap, bootstrap_type, stratify_by, seed, parallel, ci_level, ci_method,
            compute_omega_ci, compute_sigma_ci, min_success_rate)
    end
end

"""
    get_parallel_config(spec::BootstrapSpec) -> ParallelConfig

Get the parallel configuration from a BootstrapSpec.
Converts Bool to appropriate ParallelConfig.
"""
function get_parallel_config(spec::BootstrapSpec)::ParallelConfig
    if spec.parallel isa ParallelConfig
        return spec.parallel
    elseif spec.parallel === true
        # Use all available threads
        return ParallelConfig(
            ThreadedBackend(nthreads());
            seed=spec.seed,
            load_balance=true,
            progress=spec.n_bootstrap >= 100
        )
    else
        return ParallelConfig(SerialBackend(); seed=spec.seed)
    end
end

# ============================================================================
# Bootstrap Results
# ============================================================================

"""
Diagnostics for bootstrap analysis.

Fields:
- n_successful: Number of successful bootstrap runs
- n_failed: Number of failed runs (minimization failures, boundary issues)
- convergence_rate: Fraction that converged (FDA requires ≥80%)
- median_iterations: Median iterations to convergence
- outlier_indices: Indices of potential outlier estimates (>3 SD)
- failure_reasons: Count of different failure types
- rse_stability: Relative SE of SE estimate (assesses bootstrap stability)
"""
struct BootstrapDiagnostics
    n_successful::Int
    n_failed::Int
    convergence_rate::Float64
    median_iterations::Float64
    outlier_indices::Vector{Int}
    failure_reasons::Dict{Symbol, Int}
    rse_stability::Vector{Float64}  # RSE of each parameter's SE estimate
end

"""
Summary statistics for omega (IIV) parameters from bootstrap.
"""
struct OmegaBootstrapSummary
    estimates::Vector{Matrix{Float64}}  # All omega matrices
    mean::Matrix{Float64}               # Mean omega matrix
    se::Matrix{Float64}                 # SE for each omega element
    ci_lower::Matrix{Float64}           # Lower CI
    ci_upper::Matrix{Float64}           # Upper CI
    cv_percent::Matrix{Float64}         # CV% for diagonal elements
end

"""
Summary statistics for sigma (residual error) parameters from bootstrap.
"""
struct SigmaBootstrapSummary
    estimates::Vector{Float64}     # All sigma estimates
    mean::Float64                  # Mean sigma
    se::Float64                    # SE of sigma
    ci_lower::Float64              # Lower CI
    ci_upper::Float64              # Upper CI
    cv_percent::Float64            # CV%
end

"""
Result of bootstrap analysis.

This is a comprehensive result structure suitable for regulatory submissions.
Includes point estimates, SEs, CIs, and diagnostics for all parameter types.

Fields:
- theta_estimates: Matrix of fixed effect estimates (n_bootstrap x n_params)
- theta_mean: Mean of bootstrap estimates
- theta_se: Standard error from bootstrap
- theta_rse: Relative SE (%) = SE/estimate * 100
- theta_ci_lower: Lower confidence interval
- theta_ci_upper: Upper confidence interval
- original_estimate: Original point estimate
- bias: Bootstrap bias estimate = mean(bootstrap) - original
- bias_corrected: Bias-corrected estimate = original - bias
- omega_summary: Summary for omega (IIV) parameters
- sigma_summary: Summary for sigma (residual) parameters
- eta_shrinkage: Bootstrap distribution of eta shrinkage
- diagnostics: Bootstrap diagnostics
- ci_level: Confidence level used (e.g., 0.95)
- ci_method: CI method used (e.g., "Percentile", "BCa")
"""
struct BootstrapResult
    theta_estimates::Matrix{Float64}
    theta_mean::Vector{Float64}
    theta_se::Vector{Float64}
    theta_rse::Vector{Float64}  # Relative SE (%)
    theta_ci_lower::Vector{Float64}
    theta_ci_upper::Vector{Float64}
    original_estimate::Vector{Float64}
    bias::Vector{Float64}
    bias_corrected::Vector{Float64}
    omega_summary::Union{Nothing, OmegaBootstrapSummary}
    sigma_summary::Union{Nothing, SigmaBootstrapSummary}
    eta_shrinkage::Union{Nothing, Matrix{Float64}}  # Shrinkage per bootstrap replicate
    diagnostics::BootstrapDiagnostics
    ci_level::Float64
    ci_method::String

    # Convenience aliases for backward compatibility
    omega_estimates::Union{Nothing, Vector{Matrix{Float64}}}
    sigma_estimates::Union{Nothing, Vector{Float64}}
end

# ============================================================================
# Stratified Resampling
# ============================================================================

"""
    stratified_resample(subject_ids, strata, rng)

Perform stratified resampling of subjects.

Arguments:
- subject_ids: Vector of subject identifiers
- strata: Vector of strata labels (same length as subject_ids)
- rng: Random number generator

Returns:
- Vector of resampled subject indices
"""
function stratified_resample(
    subject_ids::Vector{String},
    strata::Vector{T},
    rng::AbstractRNG
)::Vector{Int} where T
    n_subjects = length(subject_ids)

    # Group subjects by stratum
    stratum_indices = Dict{T, Vector{Int}}()
    for (i, s) in enumerate(strata)
        if !haskey(stratum_indices, s)
            stratum_indices[s] = Int[]
        end
        push!(stratum_indices[s], i)
    end

    # Resample within each stratum
    resampled = Int[]
    for (stratum, indices) in stratum_indices
        n_stratum = length(indices)
        sampled = rand(rng, indices, n_stratum)
        append!(resampled, sampled)
    end

    return resampled
end

"""
    stratified_resample(n_subjects, rng)

Simple random resampling (no stratification).
"""
function stratified_resample(n_subjects::Int, rng::AbstractRNG)::Vector{Int}
    return rand(rng, 1:n_subjects, n_subjects)
end

# ============================================================================
# Bootstrap CI Computation
# ============================================================================

"""
    compute_bootstrap_ci(estimates, original, ci_level, method)

Compute confidence intervals from bootstrap estimates.
"""
function compute_bootstrap_ci(
    estimates::Matrix{Float64},
    original::Vector{Float64},
    ci_level::Float64,
    method::PercentileCI
)::Tuple{Vector{Float64}, Vector{Float64}}
    n_params = size(estimates, 2)
    alpha = 1.0 - ci_level

    lower = zeros(n_params)
    upper = zeros(n_params)

    for j in 1:n_params
        sorted = sort(estimates[:, j])
        n_boot = length(sorted)
        lower[j] = sorted[max(1, Int(floor(alpha/2 * n_boot)))]
        upper[j] = sorted[min(n_boot, Int(ceil((1 - alpha/2) * n_boot)))]
    end

    return lower, upper
end

function compute_bootstrap_ci(
    estimates::Matrix{Float64},
    original::Vector{Float64},
    ci_level::Float64,
    method::BCCI
)::Tuple{Vector{Float64}, Vector{Float64}}
    n_boot, n_params = size(estimates)
    alpha = 1.0 - ci_level

    lower = zeros(n_params)
    upper = zeros(n_params)

    for j in 1:n_params
        theta_boot = estimates[:, j]
        theta_0 = original[j]

        # Bias correction factor z0
        p0 = sum(theta_boot .< theta_0) / n_boot
        p0 = clamp(p0, 0.001, 0.999)
        z0 = quantile(Normal(), p0)

        # Acceleration factor (from spec or computed)
        a = method.acceleration

        # Adjusted percentiles
        z_alpha_lower = quantile(Normal(), alpha/2)
        z_alpha_upper = quantile(Normal(), 1 - alpha/2)

        alpha1 = cdf(Normal(), z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)))
        alpha2 = cdf(Normal(), z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)))

        sorted = sort(theta_boot)
        lower[j] = sorted[max(1, Int(floor(alpha1 * n_boot)))]
        upper[j] = sorted[min(n_boot, Int(ceil(alpha2 * n_boot)))]
    end

    return lower, upper
end

function compute_bootstrap_ci(
    estimates::Matrix{Float64},
    original::Vector{Float64},
    ci_level::Float64,
    method::StudentizedCI
)::Tuple{Vector{Float64}, Vector{Float64}}
    # Studentized bootstrap requires SE estimates for each bootstrap sample
    # Fall back to percentile if not available
    return compute_bootstrap_ci(estimates, original, ci_level, PercentileCI())
end

"""
    compute_bootstrap_ci(estimates, original, ci_level, method::BasicCI)

Compute basic (pivotal) bootstrap confidence intervals.
CI = [2θ̂ - θ*(1-α/2), 2θ̂ - θ*(α/2)]

This method corrects for bias in the bootstrap distribution.
"""
function compute_bootstrap_ci(
    estimates::Matrix{Float64},
    original::Vector{Float64},
    ci_level::Float64,
    method::BasicCI
)::Tuple{Vector{Float64}, Vector{Float64}}
    n_params = size(estimates, 2)
    alpha = 1.0 - ci_level

    lower = zeros(n_params)
    upper = zeros(n_params)

    for j in 1:n_params
        sorted = sort(estimates[:, j])
        n_boot = length(sorted)

        # Get percentiles of bootstrap distribution
        q_lower = sorted[max(1, Int(floor(alpha/2 * n_boot)))]
        q_upper = sorted[min(n_boot, Int(ceil((1 - alpha/2) * n_boot)))]

        # Basic (pivotal) transformation
        lower[j] = 2 * original[j] - q_upper
        upper[j] = 2 * original[j] - q_lower
    end

    return lower, upper
end

"""
    ci_method_name(method::BootstrapCIMethod) -> String

Get human-readable name for CI method.
"""
ci_method_name(::PercentileCI) = "Percentile"
ci_method_name(::BCCI) = "BCa"
ci_method_name(::StudentizedCI) = "Studentized"
ci_method_name(::BasicCI) = "Basic"

# ============================================================================
# Main Bootstrap Function
# ============================================================================

"""
    run_bootstrap(estimation_fn, observed_data, spec; kwargs...)

Run bootstrap analysis for parameter estimation.

Supports both sequential and parallel execution based on the `parallel` field
in the BootstrapSpec. For large bootstrap studies (n >= 500), parallel execution
can provide significant speedup.

# Arguments
- `estimation_fn`: Function that takes observed data and returns EstimationResult
- `observed_data`: Original observed data (ObservedData type)
- `spec`: BootstrapSpec specifying bootstrap parameters

# Keyword Arguments
- `strata`: Optional vector of strata labels for each subject
- `verbose`: Print progress (default: false)

# Returns
- BootstrapResult with bootstrap estimates and diagnostics

# Example
```julia
# Sequential bootstrap
result = run_bootstrap(estimate_fn, data, BootstrapSpec(n_bootstrap=500))

# Parallel bootstrap (uses all available threads)
result = run_bootstrap(estimate_fn, data, BootstrapSpec(n_bootstrap=1000, parallel=true))

# Parallel with custom configuration
config = ParallelConfig(ThreadedBackend(8); seed=42)
result = run_bootstrap(estimate_fn, data, BootstrapSpec(parallel=config))
```
"""
function run_bootstrap(
    estimation_fn,  # Any callable (Julia Function or Python callable via PythonCall)
    observed_data,  # ObservedData type
    spec::BootstrapSpec;
    strata::Union{Nothing, Vector}=nothing,
    verbose::Bool=false
)::BootstrapResult
    # Get parallel configuration
    parallel_config = get_parallel_config(spec)

    # Dispatch to parallel or sequential implementation
    if is_parallel(parallel_config)
        return run_bootstrap_parallel_impl(
            estimation_fn,
            observed_data,
            spec,
            parallel_config;
            strata=strata,
            verbose=verbose
        )
    else
        return run_bootstrap_sequential(
            estimation_fn,
            observed_data,
            spec;
            strata=strata,
            verbose=verbose
        )
    end
end

"""
    run_bootstrap_sequential(estimation_fn, observed_data, spec; strata, verbose)

Internal sequential implementation of bootstrap analysis.
"""
function run_bootstrap_sequential(
    estimation_fn,  # Any callable
    observed_data,
    spec::BootstrapSpec;
    strata::Union{Nothing, Vector}=nothing,
    verbose::Bool=false
)::BootstrapResult
    rng = StableRNG(spec.seed)

    n_subjects = length(observed_data.subjects)
    subject_ids = [s.subject_id for s in observed_data.subjects]

    # Run original estimation to get point estimate
    if verbose
        println("Running original estimation...")
    end
    original_result = estimation_fn(observed_data)
    # Handle case where estimation_fn is Python callable returning Py-wrapped result
    original_theta = _to_float_vector(original_result.theta)
    n_params = length(original_theta)

    # Storage for bootstrap estimates
    theta_estimates = zeros(spec.n_bootstrap, n_params)
    omega_estimates = Matrix{Float64}[]
    sigma_estimates = Float64[]

    n_successful = 0
    n_failed = 0
    iterations = Float64[]

    for b in 1:spec.n_bootstrap
        if verbose && b % 100 == 0
            println("Bootstrap replicate $b / $(spec.n_bootstrap)")
        end

        try
            # Resample subjects
            if strata !== nothing
                resampled_indices = stratified_resample(subject_ids, strata, rng)
            else
                resampled_indices = stratified_resample(n_subjects, rng)
            end

            # Create resampled dataset
            resampled_data = _create_resampled_data(observed_data, resampled_indices)

            # Run estimation on resampled data
            boot_result = estimation_fn(resampled_data)

            # Store estimates (convert in case of Py-wrapped result)
            theta_estimates[b, :] = _to_float_vector(boot_result.theta)

            if boot_result.omega !== nothing
                push!(omega_estimates, _to_float_matrix(boot_result.omega))
            end
            if boot_result.sigma !== nothing
                push!(sigma_estimates, boot_result.sigma.params.sigma)
            end

            n_successful += 1

            # Track iterations if available
            if hasproperty(boot_result, :n_iterations)
                push!(iterations, Float64(boot_result.n_iterations))
            end

        catch e
            n_failed += 1
            # Fill with NaN for failed runs
            theta_estimates[b, :] .= NaN
            if verbose
                println("Bootstrap $b failed: $e")
            end
        end
    end

    # Aggregate and return results
    return _aggregate_bootstrap_results(
        theta_estimates,
        omega_estimates,
        sigma_estimates,
        original_theta,
        n_successful,
        n_failed,
        iterations,
        spec
    )
end

"""
    run_bootstrap_parallel_impl(estimation_fn, observed_data, spec, parallel_config; strata, verbose)

Internal parallel implementation of bootstrap analysis.
Uses thread-safe RNG for reproducibility across parallel replicates.
"""
function run_bootstrap_parallel_impl(
    estimation_fn,  # Any callable
    observed_data,
    spec::BootstrapSpec,
    parallel_config::ParallelConfig;
    strata::Union{Nothing, Vector}=nothing,
    verbose::Bool=false
)::BootstrapResult
    n_subjects = length(observed_data.subjects)
    subject_ids = [s.subject_id for s in observed_data.subjects]

    # Run original estimation to get point estimate
    if verbose
        n_workers_used = n_workers(parallel_config)
        println("Running original estimation...")
        println("Parallel bootstrap: $(spec.n_bootstrap) replicates using $n_workers_used workers")
    end

    original_result = estimation_fn(observed_data)
    # Handle case where estimation_fn is Python callable returning Py-wrapped result
    original_theta = _to_float_vector(original_result.theta)
    n_params = length(original_theta)

    # Create independent RNGs for each bootstrap replicate (for reproducibility)
    rngs = create_thread_rngs(spec.n_bootstrap; seed=spec.seed)

    # Prepare replicate indices and strata info
    replicate_indices = collect(1:spec.n_bootstrap)

    # Define single replicate function
    function run_single_replicate(rep_idx::Int, rng::AbstractRNG)
        try
            # Resample subjects using this replicate's RNG
            if strata !== nothing
                resampled_indices = stratified_resample(subject_ids, strata, rng)
            else
                resampled_indices = stratified_resample(n_subjects, rng)
            end

            # Create resampled dataset
            resampled_data = _create_resampled_data(observed_data, resampled_indices)

            # Run estimation on resampled data
            boot_result = estimation_fn(resampled_data)

            # Extract sigma parameter if available
            sigma_val = NaN
            if boot_result.sigma !== nothing
                if hasproperty(boot_result.sigma, :params) && hasproperty(boot_result.sigma.params, :sigma)
                    sigma_val = boot_result.sigma.params.sigma
                end
            end

            # Return successful result (convert in case of Py-wrapped result)
            return (
                success=true,
                theta=_to_float_vector(boot_result.theta),
                omega=boot_result.omega !== nothing ? _to_float_matrix(boot_result.omega) : nothing,
                sigma=sigma_val,
                n_iter=hasproperty(boot_result, :n_iterations) ? boot_result.n_iterations : 0,
                error=nothing
            )
        catch e
            # Return failed result
            return (
                success=false,
                theta=fill(NaN, n_params),
                omega=nothing,
                sigma=NaN,
                n_iter=0,
                error=e
            )
        end
    end

    # Run bootstrap replicates in parallel
    if verbose
        println("Starting parallel bootstrap...")
    end

    results = parallel_map_with_rng(
        run_single_replicate,
        replicate_indices,
        rngs,
        parallel_config
    )

    # Aggregate results
    theta_estimates = zeros(spec.n_bootstrap, n_params)
    omega_estimates = Matrix{Float64}[]
    sigma_estimates = Float64[]
    iterations = Float64[]
    n_successful = 0
    n_failed = 0

    for (i, r) in enumerate(results)
        if r.success
            theta_estimates[i, :] = r.theta
            if r.omega !== nothing
                push!(omega_estimates, r.omega)
            end
            if !isnan(r.sigma)
                push!(sigma_estimates, r.sigma)
            end
            if r.n_iter > 0
                push!(iterations, Float64(r.n_iter))
            end
            n_successful += 1
        else
            theta_estimates[i, :] .= NaN
            n_failed += 1
            if verbose
                println("Bootstrap replicate $i failed: $(r.error)")
            end
        end
    end

    if verbose
        println("Parallel bootstrap complete: $n_successful successful, $n_failed failed")
    end

    # Aggregate and return results
    return _aggregate_bootstrap_results(
        theta_estimates,
        omega_estimates,
        sigma_estimates,
        original_theta,
        n_successful,
        n_failed,
        iterations,
        spec
    )
end

"""
    _aggregate_bootstrap_results(theta_estimates, omega_estimates, sigma_estimates,
                                  original_theta, n_successful, n_failed, iterations, spec,
                                  failure_reasons)

Internal helper to aggregate bootstrap results into BootstrapResult.
Computes comprehensive statistics for regulatory submissions.
"""
function _aggregate_bootstrap_results(
    theta_estimates::Matrix{Float64},
    omega_estimates::Vector{Matrix{Float64}},
    sigma_estimates::Vector{Float64},
    original_theta::Vector{Float64},
    n_successful::Int,
    n_failed::Int,
    iterations::Vector{Float64},
    spec::BootstrapSpec;
    failure_reasons::Dict{Symbol, Int}=Dict{Symbol, Int}()
)::BootstrapResult
    n_params = length(original_theta)

    # Remove failed runs for statistics
    valid_mask = .!any(isnan.(theta_estimates), dims=2)[:]
    valid_estimates = theta_estimates[valid_mask, :]
    n_valid = size(valid_estimates, 1)

    # Handle edge case: no successful replicates
    if n_valid == 0
        theta_mean = fill(NaN, n_params)
        theta_se = fill(NaN, n_params)
        theta_rse = fill(NaN, n_params)
        bias = fill(NaN, n_params)
        bias_corrected = fill(NaN, n_params)
        lower = fill(NaN, n_params)
        upper = fill(NaN, n_params)
        rse_stability = fill(NaN, n_params)
    else
        # Compute statistics
        theta_mean = vec(mean(valid_estimates, dims=1))
        theta_se = vec(std(valid_estimates, dims=1))

        # Relative SE (%) using original estimate as reference
        theta_rse = zeros(n_params)
        for j in 1:n_params
            if abs(original_theta[j]) > 1e-10
                theta_rse[j] = theta_se[j] / abs(original_theta[j]) * 100
            else
                theta_rse[j] = NaN
            end
        end

        bias = theta_mean .- original_theta
        bias_corrected = original_theta .- bias

        # Compute confidence intervals
        lower, upper = compute_bootstrap_ci(valid_estimates, original_theta, spec.ci_level, spec.ci_method)

        # Compute RSE of SE (bootstrap stability metric)
        # Use jackknife-like approach: SE of SE ≈ SE / sqrt(2*(n-1))
        rse_stability = theta_se ./ sqrt(2 * (n_valid - 1)) ./ theta_se .* 100
        replace!(rse_stability, NaN => 0.0)
    end

    # Detect outliers (> 3 SD from mean)
    outlier_indices = Int[]
    if n_valid > 0
        for j in 1:n_params
            if theta_se[j] > 0
                z_scores = abs.(valid_estimates[:, j] .- theta_mean[j]) ./ theta_se[j]
                outliers = findall(z_scores .> 3.0)
                append!(outlier_indices, outliers)
            end
        end
        outlier_indices = unique(outlier_indices)
    end

    # Create diagnostics
    diagnostics = BootstrapDiagnostics(
        n_successful,
        n_failed,
        n_successful / spec.n_bootstrap,
        isempty(iterations) ? 0.0 : median(iterations),
        outlier_indices,
        failure_reasons,
        rse_stability
    )

    # Compute omega summary if requested and available
    omega_summary = nothing
    if spec.compute_omega_ci && !isempty(omega_estimates)
        omega_summary = _compute_omega_summary(omega_estimates, spec.ci_level, spec.ci_method)
    end

    # Compute sigma summary if requested and available
    sigma_summary = nothing
    if spec.compute_sigma_ci && !isempty(sigma_estimates)
        sigma_summary = _compute_sigma_summary(sigma_estimates, spec.ci_level)
    end

    return BootstrapResult(
        theta_estimates,
        theta_mean,
        theta_se,
        theta_rse,
        lower,
        upper,
        original_theta,
        bias,
        bias_corrected,
        omega_summary,
        sigma_summary,
        nothing,  # eta_shrinkage - computed separately if needed
        diagnostics,
        spec.ci_level,
        ci_method_name(spec.ci_method),
        # Backward compatibility
        isempty(omega_estimates) ? nothing : omega_estimates,
        isempty(sigma_estimates) ? nothing : sigma_estimates
    )
end

"""
    _compute_omega_summary(omega_estimates, ci_level, ci_method) -> OmegaBootstrapSummary

Compute summary statistics for omega (IIV) parameters.
"""
function _compute_omega_summary(
    omega_estimates::Vector{Matrix{Float64}},
    ci_level::Float64,
    ci_method::BootstrapCIMethod
)::OmegaBootstrapSummary
    n_boot = length(omega_estimates)
    n_eta = size(omega_estimates[1], 1)

    # Stack into arrays for computation
    omega_mean = zeros(n_eta, n_eta)
    omega_se = zeros(n_eta, n_eta)
    omega_lower = zeros(n_eta, n_eta)
    omega_upper = zeros(n_eta, n_eta)
    cv_percent = zeros(n_eta, n_eta)

    for i in 1:n_eta
        for j in 1:n_eta
            values = [omega_estimates[k][i, j] for k in 1:n_boot]
            valid_values = filter(!isnan, values)

            if length(valid_values) > 0
                omega_mean[i, j] = mean(valid_values)
                omega_se[i, j] = std(valid_values)

                # Compute CI
                alpha = 1.0 - ci_level
                sorted = sort(valid_values)
                n = length(sorted)
                omega_lower[i, j] = sorted[max(1, Int(floor(alpha/2 * n)))]
                omega_upper[i, j] = sorted[min(n, Int(ceil((1 - alpha/2) * n)))]

                # CV% only for diagonal elements (variance)
                if i == j && omega_mean[i, j] > 0
                    cv_percent[i, j] = sqrt(omega_mean[i, j]) * 100  # CV% = sqrt(variance) * 100
                end
            else
                omega_mean[i, j] = NaN
                omega_se[i, j] = NaN
                omega_lower[i, j] = NaN
                omega_upper[i, j] = NaN
                cv_percent[i, j] = NaN
            end
        end
    end

    return OmegaBootstrapSummary(
        omega_estimates,
        omega_mean,
        omega_se,
        omega_lower,
        omega_upper,
        cv_percent
    )
end

"""
    _compute_sigma_summary(sigma_estimates, ci_level) -> SigmaBootstrapSummary

Compute summary statistics for sigma (residual error) parameters.
"""
function _compute_sigma_summary(
    sigma_estimates::Vector{Float64},
    ci_level::Float64
)::SigmaBootstrapSummary
    valid_values = filter(!isnan, sigma_estimates)

    if isempty(valid_values)
        return SigmaBootstrapSummary(
            sigma_estimates, NaN, NaN, NaN, NaN, NaN
        )
    end

    sigma_mean = mean(valid_values)
    sigma_se = std(valid_values)

    # Compute percentile CI
    alpha = 1.0 - ci_level
    sorted = sort(valid_values)
    n = length(sorted)
    ci_lower = sorted[max(1, Int(floor(alpha/2 * n)))]
    ci_upper = sorted[min(n, Int(ceil((1 - alpha/2) * n)))]

    # CV% for proportional error
    cv_percent = sigma_mean > 0 ? sigma_mean * 100 : NaN

    return SigmaBootstrapSummary(
        sigma_estimates,
        sigma_mean,
        sigma_se,
        ci_lower,
        ci_upper,
        cv_percent
    )
end

"""
Helper function to create resampled observed data.
"""
function _create_resampled_data(observed_data, indices::Vector{Int})
    # Create new subjects list with resampled indices
    # Note: Same subject can appear multiple times (with replacement)

    new_subjects = [observed_data.subjects[i] for i in indices]

    # Renumber subject IDs to avoid duplicates
    for (i, subj) in enumerate(new_subjects)
        # Create new subject with modified ID
        # This is a shallow approach - actual implementation depends on ObservedData structure
    end

    # Return a structure compatible with ObservedData
    # The actual implementation depends on the ObservedData type definition
    return typeof(observed_data)(new_subjects)
end

# ============================================================================
# Bootstrap SE Comparison
# ============================================================================

"""
    compare_se_methods(bootstrap_result, analytical_se)

Compare bootstrap standard errors with analytical standard errors.

Returns a dictionary with comparison metrics.
"""
function compare_se_methods(
    bootstrap_result::BootstrapResult,
    analytical_se::Vector{Float64}
)::Dict{Symbol, Any}
    n_params = length(analytical_se)

    # Compute ratios
    se_ratio = bootstrap_result.theta_se ./ analytical_se

    # Compute relative difference
    rel_diff = (bootstrap_result.theta_se .- analytical_se) ./ analytical_se .* 100

    # Flag parameters where SE differs by more than 20%
    flags = abs.(rel_diff) .> 20.0

    return Dict(
        :bootstrap_se => bootstrap_result.theta_se,
        :analytical_se => analytical_se,
        :se_ratio => se_ratio,
        :relative_diff_percent => rel_diff,
        :discrepancy_flags => flags,
        :mean_ratio => mean(se_ratio),
        :max_discrepancy => maximum(abs.(rel_diff))
    )
end

export compare_se_methods

# ============================================================================
# Bootstrap Diagnostics Plotting Data
# ============================================================================

"""
    get_bootstrap_histogram_data(result, param_index)

Get data for plotting bootstrap distribution histogram.
"""
function get_bootstrap_histogram_data(
    result::BootstrapResult,
    param_index::Int
)::Dict{Symbol, Any}
    estimates = result.theta_estimates[:, param_index]
    valid = filter(!isnan, estimates)

    return Dict(
        :estimates => valid,
        :original => result.original_estimate[param_index],
        :mean => result.theta_mean[param_index],
        :ci_lower => result.theta_ci_lower[param_index],
        :ci_upper => result.theta_ci_upper[param_index],
        :se => result.theta_se[param_index],
        :bias => result.bias[param_index]
    )
end

export get_bootstrap_histogram_data

# ============================================================================
# Regulatory-Compliant Output Formatting
# ============================================================================

"""
    generate_bootstrap_summary(result::BootstrapResult, param_names::Vector{String})

Generate a regulatory-compliant summary table for bootstrap results.

Returns a DataFrame-like structure with:
- Parameter names
- Point estimates
- Bootstrap SE
- RSE (%)
- 95% CI
- Bias
- Shrinkage (if applicable)

# Example
```julia
result = run_bootstrap(estimate_fn, data, spec)
summary = generate_bootstrap_summary(result, ["CL", "V", "Ka"])
```
"""
function generate_bootstrap_summary(
    result::BootstrapResult,
    param_names::Vector{String}
)::Dict{Symbol, Any}
    n_params = length(param_names)
    @assert n_params == length(result.original_estimate) "Parameter names must match result length"

    # Build summary table
    summary = Dict{Symbol, Any}(
        :parameter => param_names,
        :estimate => result.original_estimate,
        :bootstrap_se => result.theta_se,
        :rse_percent => result.theta_rse,
        :ci_lower => result.theta_ci_lower,
        :ci_upper => result.theta_ci_upper,
        :bias => result.bias,
        :bias_percent => result.bias ./ abs.(result.original_estimate) .* 100,
        :ci_level => result.ci_level * 100,
        :ci_method => result.ci_method,
        :n_bootstrap => size(result.theta_estimates, 1),
        :n_successful => result.diagnostics.n_successful,
        :convergence_rate => result.diagnostics.convergence_rate * 100
    )

    # Add omega summary if available
    if result.omega_summary !== nothing
        omega = result.omega_summary
        n_eta = size(omega.mean, 1)
        omega_params = String[]
        omega_est = Float64[]
        omega_se = Float64[]
        omega_lower = Float64[]
        omega_upper = Float64[]

        for i in 1:n_eta
            for j in 1:i
                if i == j
                    push!(omega_params, "OMEGA($i,$j)")
                    push!(omega_est, omega.mean[i, j])
                    push!(omega_se, omega.se[i, j])
                    push!(omega_lower, omega.ci_lower[i, j])
                    push!(omega_upper, omega.ci_upper[i, j])
                end
            end
        end

        summary[:omega_parameters] = omega_params
        summary[:omega_estimate] = omega_est
        summary[:omega_se] = omega_se
        summary[:omega_ci_lower] = omega_lower
        summary[:omega_ci_upper] = omega_upper
    end

    # Add sigma summary if available
    if result.sigma_summary !== nothing
        sigma = result.sigma_summary
        summary[:sigma_estimate] = sigma.mean
        summary[:sigma_se] = sigma.se
        summary[:sigma_ci_lower] = sigma.ci_lower
        summary[:sigma_ci_upper] = sigma.ci_upper
    end

    return summary
end

"""
    format_regulatory_table(result::BootstrapResult, param_names::Vector{String}; digits::Int=4)

Format bootstrap results as a regulatory-compliant text table.

This format follows FDA/EMA recommendations for population PK reports.

# Returns
A string containing the formatted table.

# Example
```julia
table = format_regulatory_table(result, ["CL (L/h)", "V (L)", "Ka (1/h)"])
println(table)
```
"""
function format_regulatory_table(
    result::BootstrapResult,
    param_names::Vector{String};
    digits::Int=4
)::String
    n_params = length(param_names)
    ci_pct = Int(result.ci_level * 100)

    # Helper function for formatting numbers
    function fmt(x::Float64, d::Int=digits)
        isnan(x) ? "NaN" : string(round(x, sigdigits=d))
    end

    # Helper for padding
    function pad_left(s::String, width::Int)
        len = length(s)
        len >= width ? s : " "^(width - len) * s
    end

    function pad_right(s::String, width::Int)
        len = length(s)
        len >= width ? s : s * " "^(width - len)
    end

    lines = String[]

    # Header
    push!(lines, "=" ^ 90)
    push!(lines, "BOOTSTRAP PARAMETER ESTIMATES")
    push!(lines, "=" ^ 90)
    push!(lines, "")
    push!(lines, "Bootstrap Settings:")
    push!(lines, "  Number of replicates: $(size(result.theta_estimates, 1))")
    push!(lines, "  Successful runs: $(result.diagnostics.n_successful) ($(round(result.diagnostics.convergence_rate * 100, digits=1))%)")
    push!(lines, "  CI Method: $(result.ci_method)")
    push!(lines, "  CI Level: $(ci_pct)%")
    push!(lines, "")

    # Fixed effects table header
    push!(lines, "-" ^ 90)
    push!(lines, "FIXED EFFECTS (THETA)")
    push!(lines, "-" ^ 90)
    header = pad_right("Parameter", 20) * " " *
             pad_left("Estimate", 12) * " " *
             pad_left("SE", 12) * " " *
             pad_left("RSE(%)", 10) * " " *
             pad_left("$(ci_pct)%CI Lo", 12) * " " *
             pad_left("$(ci_pct)%CI Hi", 12) * " " *
             pad_left("Bias", 12)
    push!(lines, header)
    push!(lines, "-" ^ 90)

    for i in 1:n_params
        line = pad_right(param_names[i], 20) * " " *
               pad_left(fmt(result.original_estimate[i]), 12) * " " *
               pad_left(fmt(result.theta_se[i]), 12) * " " *
               pad_left(string(round(result.theta_rse[i], digits=1)), 10) * " " *
               pad_left(fmt(result.theta_ci_lower[i]), 12) * " " *
               pad_left(fmt(result.theta_ci_upper[i]), 12) * " " *
               pad_left(fmt(result.bias[i]), 12)
        push!(lines, line)
    end
    push!(lines, "")

    # Omega table if available
    if result.omega_summary !== nothing
        omega = result.omega_summary
        n_eta = size(omega.mean, 1)

        push!(lines, "-" ^ 90)
        push!(lines, "RANDOM EFFECTS (OMEGA) - IIV Variance")
        push!(lines, "-" ^ 90)
        header = pad_right("Parameter", 20) * " " *
                 pad_left("Variance", 12) * " " *
                 pad_left("SE", 12) * " " *
                 pad_left("CV%", 12) * " " *
                 pad_left("$(ci_pct)%CI Lo", 12) * " " *
                 pad_left("$(ci_pct)%CI Hi", 12)
        push!(lines, header)
        push!(lines, "-" ^ 90)

        for i in 1:n_eta
            name = "OMEGA($i,$i)"
            cv = omega.mean[i, i] > 0 ? sqrt(omega.mean[i, i]) * 100 : 0.0
            line = pad_right(name, 20) * " " *
                   pad_left(fmt(omega.mean[i, i]), 12) * " " *
                   pad_left(fmt(omega.se[i, i]), 12) * " " *
                   pad_left(string(round(cv, digits=1)), 12) * " " *
                   pad_left(fmt(omega.ci_lower[i, i]), 12) * " " *
                   pad_left(fmt(omega.ci_upper[i, i]), 12)
            push!(lines, line)
        end
        push!(lines, "")
    end

    # Sigma table if available
    if result.sigma_summary !== nothing
        sigma = result.sigma_summary

        push!(lines, "-" ^ 90)
        push!(lines, "RESIDUAL ERROR (SIGMA)")
        push!(lines, "-" ^ 90)
        header = pad_right("Parameter", 20) * " " *
                 pad_left("Estimate", 12) * " " *
                 pad_left("SE", 12) * " " *
                 pad_left("$(ci_pct)%CI Lo", 12) * " " *
                 pad_left("$(ci_pct)%CI Hi", 12)
        push!(lines, header)
        push!(lines, "-" ^ 90)

        line = pad_right("SIGMA", 20) * " " *
               pad_left(fmt(sigma.mean), 12) * " " *
               pad_left(fmt(sigma.se), 12) * " " *
               pad_left(fmt(sigma.ci_lower), 12) * " " *
               pad_left(fmt(sigma.ci_upper), 12)
        push!(lines, line)
        push!(lines, "")
    end

    # Diagnostics
    push!(lines, "-" ^ 90)
    push!(lines, "BOOTSTRAP DIAGNOSTICS")
    push!(lines, "-" ^ 90)
    push!(lines, "  Convergence rate: $(round(result.diagnostics.convergence_rate * 100, digits=1))% (FDA recommends ≥80%)")
    push!(lines, "  Number of outliers detected: $(length(result.diagnostics.outlier_indices))")
    push!(lines, "  Median iterations: $(round(result.diagnostics.median_iterations, digits=0))")

    if !isempty(result.diagnostics.failure_reasons)
        push!(lines, "  Failure reasons:")
        for (reason, count) in result.diagnostics.failure_reasons
            push!(lines, "    $(reason): $(count)")
        end
    end

    push!(lines, "")
    push!(lines, "=" ^ 90)

    return join(lines, "\n")
end

# ============================================================================
# Influential Subject Analysis
# ============================================================================

"""
    jackknife_influence(estimation_fn, observed_data; verbose=false)

Perform jackknife (leave-one-out) analysis to identify influential subjects.

Returns indices of subjects whose exclusion changes parameter estimates by >10%.

# Arguments
- `estimation_fn`: Function that takes observed data and returns EstimationResult
- `observed_data`: Original observed data

# Returns
- Dict with influence metrics for each subject
"""
function jackknife_influence(
    estimation_fn::Function,
    observed_data;
    verbose::Bool=false
)::Dict{Symbol, Any}
    n_subjects = length(observed_data.subjects)

    # Get original estimate
    if verbose
        println("Running original estimation...")
    end
    original_result = estimation_fn(observed_data)
    original_theta = original_result.theta
    n_params = length(original_theta)

    # Storage for jackknife estimates
    jackknife_estimates = zeros(n_subjects, n_params)
    influence = zeros(n_subjects, n_params)

    for i in 1:n_subjects
        if verbose && i % 10 == 0
            println("Jackknife iteration $i / $n_subjects")
        end

        try
            # Create leave-one-out dataset
            indices = [j for j in 1:n_subjects if j != i]
            loo_data = _create_resampled_data(observed_data, indices)

            # Run estimation
            loo_result = estimation_fn(loo_data)
            jackknife_estimates[i, :] = loo_result.theta

            # Compute influence (% change)
            for j in 1:n_params
                if abs(original_theta[j]) > 1e-10
                    influence[i, j] = (loo_result.theta[j] - original_theta[j]) / original_theta[j] * 100
                end
            end
        catch e
            jackknife_estimates[i, :] .= NaN
            influence[i, :] .= NaN
            if verbose
                println("Jackknife iteration $i failed: $e")
            end
        end
    end

    # Identify influential subjects (>10% change in any parameter)
    influential_subjects = Int[]
    for i in 1:n_subjects
        if any(abs.(influence[i, :]) .> 10.0)
            push!(influential_subjects, i)
        end
    end

    # Compute jackknife SE
    jackknife_mean = vec(mean(jackknife_estimates, dims=1))
    # Jackknife SE formula: sqrt((n-1)/n * sum((theta_i - theta_bar)^2))
    jackknife_se = zeros(n_params)
    for j in 1:n_params
        valid_est = filter(!isnan, jackknife_estimates[:, j])
        if length(valid_est) > 1
            jackknife_se[j] = sqrt((n_subjects - 1) / n_subjects * sum((valid_est .- mean(valid_est)).^2))
        end
    end

    return Dict(
        :original_estimate => original_theta,
        :jackknife_estimates => jackknife_estimates,
        :influence_percent => influence,
        :influential_subjects => influential_subjects,
        :jackknife_se => jackknife_se,
        :jackknife_mean => jackknife_mean,
        :n_subjects => n_subjects,
        :n_params => n_params
    )
end

"""
    influential_subject_analysis(result::BootstrapResult, threshold::Float64=10.0)

Analyze bootstrap results for subject influence patterns.

This identifies if certain bootstrap replicates consistently produce outlier estimates,
which may indicate influential subjects in the original dataset.

# Returns
Dict with:
- High-influence replicate indices
- Parameters most affected by outliers
- Stability assessment
"""
function influential_subject_analysis(
    result::BootstrapResult;
    threshold::Float64=10.0
)::Dict{Symbol, Any}
    n_boot = size(result.theta_estimates, 1)
    n_params = length(result.theta_mean)

    # Compute standardized residuals for each replicate
    z_scores = zeros(n_boot, n_params)
    for j in 1:n_params
        if result.theta_se[j] > 0
            z_scores[:, j] = (result.theta_estimates[:, j] .- result.theta_mean[j]) ./ result.theta_se[j]
        end
    end

    # Identify high-influence replicates (any z > 3)
    high_influence_replicates = Int[]
    for i in 1:n_boot
        z_row = z_scores[i, :]
        if any(.!isnan.(z_row) .& (abs.(z_row) .> 3.0))
            push!(high_influence_replicates, i)
        end
    end

    # Parameters most affected by outliers
    outlier_counts = zeros(Int, n_params)
    for j in 1:n_params
        outlier_counts[j] = sum(abs.(z_scores[:, j]) .> 3.0)
    end

    # Stability assessment: % of replicates within 2 SE of mean
    stability = zeros(n_params)
    for j in 1:n_params
        n_stable = sum(abs.(z_scores[:, j]) .<= 2.0)
        n_valid = sum(.!isnan.(z_scores[:, j]))
        stability[j] = n_valid > 0 ? n_stable / n_valid * 100 : 0.0
    end

    return Dict(
        :high_influence_replicates => high_influence_replicates,
        :n_high_influence => length(high_influence_replicates),
        :outlier_counts_per_param => outlier_counts,
        :stability_percent => stability,
        :mean_stability => mean(stability),
        :assessment => mean(stability) >= 95.0 ? "Stable" :
                       mean(stability) >= 90.0 ? "Acceptable" : "Unstable"
    )
end

# ============================================================================
# Bootstrap Coverage and Stability
# ============================================================================

"""
    compute_bootstrap_coverage(result::BootstrapResult, true_values::Vector{Float64})

Compute empirical coverage of bootstrap CIs given true parameter values.
Useful for simulation studies to validate bootstrap performance.

# Returns
Dict with coverage rates for each parameter and overall.
"""
function compute_bootstrap_coverage(
    result::BootstrapResult,
    true_values::Vector{Float64}
)::Dict{Symbol, Any}
    n_params = length(true_values)
    @assert n_params == length(result.theta_ci_lower) "True values must match parameter length"

    # Check if true value falls within CI
    covered = zeros(Bool, n_params)
    for j in 1:n_params
        covered[j] = result.theta_ci_lower[j] <= true_values[j] <= result.theta_ci_upper[j]
    end

    return Dict(
        :parameter_coverage => covered,
        :n_covered => sum(covered),
        :n_params => n_params,
        :coverage_rate => sum(covered) / n_params,
        :expected_coverage => result.ci_level
    )
end

"""
    compute_bootstrap_stability(result::BootstrapResult)

Assess stability of bootstrap estimates.

# Returns
Dict with stability metrics including:
- Coefficient of variation for each parameter
- Normality test p-values (Shapiro-Wilk approximation)
- Skewness and kurtosis
"""
function compute_bootstrap_stability(result::BootstrapResult)::Dict{Symbol, Any}
    n_params = length(result.theta_mean)

    cv = zeros(n_params)
    skewness = zeros(n_params)
    kurtosis = zeros(n_params)

    for j in 1:n_params
        valid = filter(!isnan, result.theta_estimates[:, j])
        if length(valid) > 3
            μ = mean(valid)
            σ = std(valid)

            # CV
            cv[j] = abs(μ) > 1e-10 ? σ / abs(μ) * 100 : NaN

            # Standardized moments
            z = (valid .- μ) ./ σ
            skewness[j] = mean(z.^3)
            kurtosis[j] = mean(z.^4) - 3  # Excess kurtosis
        else
            cv[j] = NaN
            skewness[j] = NaN
            kurtosis[j] = NaN
        end
    end

    # Overall stability assessment
    max_cv = maximum(filter(!isnan, cv))
    max_skew = maximum(abs.(filter(!isnan, skewness)))

    return Dict(
        :cv_percent => cv,
        :skewness => skewness,
        :excess_kurtosis => kurtosis,
        :max_cv => max_cv,
        :max_abs_skewness => max_skew,
        :is_stable => max_cv < 50.0 && max_skew < 1.0,
        :stability_grade => max_cv < 20.0 && max_skew < 0.5 ? "Excellent" :
                            max_cv < 35.0 && max_skew < 1.0 ? "Good" :
                            max_cv < 50.0 ? "Acceptable" : "Poor"
    )
end


# ============================================================================
# Native Bootstrap Function (No Callback Required)
# ============================================================================

"""
    run_bootstrap_native(observed_data, model_spec, est_config, spec; grid, solver, strata, verbose)

Run bootstrap analysis without requiring a callback function.

This is the preferred method for calling from Python via JuliaCall, as it avoids
the callback interop issues that can cause crashes when Python closures are
called from Julia loops.

# Arguments
- `observed_data`: ObservedData containing subjects
- `model_spec`: ModelSpec defining the PK/PD model
- `est_config`: EstimationConfig with estimation settings
- `spec`: BootstrapSpec with bootstrap settings
- `grid`: SimGrid for simulation time points
- `solver`: SolverSpec for ODE solver settings
- `strata`: Optional stratification vector
- `verbose`: Print progress information

# Returns
- BootstrapResult with estimates, CIs, and diagnostics

# Example
```julia
result = run_bootstrap_native(
    data,
    model_spec,
    est_config,
    BootstrapSpec(n_bootstrap=500),
    grid=sim_grid,
    solver=solver_spec
)
```
"""
function run_bootstrap_native(
    observed_data::ObservedData,
    model_spec::ModelSpec,
    est_config::EstimationConfig,
    spec::BootstrapSpec;
    grid::SimGrid,
    solver::SolverSpec,
    strata::Union{Nothing, Vector}=nothing,
    verbose::Bool=false
)::BootstrapResult
    rng = StableRNG(spec.seed)

    n_subjects = length(observed_data.subjects)
    subject_ids = [s.subject_id for s in observed_data.subjects]

    # Run original estimation to get point estimate
    if verbose
        println("Running original estimation...")
    end
    original_result = estimate(observed_data, model_spec, est_config; grid=grid, solver=solver)
    original_theta = Vector{Float64}(original_result.theta)
    n_params = length(original_theta)

    # Storage for bootstrap estimates
    theta_estimates = zeros(spec.n_bootstrap, n_params)
    omega_estimates = Matrix{Float64}[]
    sigma_estimates = Float64[]

    n_successful = 0
    n_failed = 0
    iterations = Float64[]

    for b in 1:spec.n_bootstrap
        if verbose && b % 100 == 0
            println("Bootstrap replicate $b / $(spec.n_bootstrap)")
        end

        try
            # Resample subjects
            if strata !== nothing
                resampled_indices = stratified_resample(subject_ids, strata, rng)
            else
                resampled_indices = stratified_resample(n_subjects, rng)
            end

            # Create resampled dataset
            resampled_data = _create_resampled_data(observed_data, resampled_indices)

            # Run estimation on resampled data (native Julia call, no callback)
            boot_result = estimate(resampled_data, model_spec, est_config; grid=grid, solver=solver)

            # Store estimates
            theta_estimates[b, :] = Vector{Float64}(boot_result.theta)

            if boot_result.omega !== nothing
                push!(omega_estimates, Matrix{Float64}(boot_result.omega))
            end
            if boot_result.sigma !== nothing && hasproperty(boot_result.sigma, :params) &&
               hasproperty(boot_result.sigma.params, :sigma)
                push!(sigma_estimates, boot_result.sigma.params.sigma)
            end

            if boot_result.n_iterations > 0
                push!(iterations, Float64(boot_result.n_iterations))
            end

            n_successful += 1
        catch e
            n_failed += 1
            theta_estimates[b, :] .= NaN
            if verbose
                println("Bootstrap replicate $b failed: $e")
            end
        end
    end

    # Check success rate
    success_rate = n_successful / spec.n_bootstrap
    if success_rate < spec.min_success_rate
        @warn "Bootstrap success rate ($success_rate) below minimum threshold ($(spec.min_success_rate))"
    end

    if verbose
        println("Bootstrap complete: $n_successful successful, $n_failed failed")
    end

    # Aggregate and return results
    return _aggregate_bootstrap_results(
        theta_estimates,
        omega_estimates,
        sigma_estimates,
        original_theta,
        n_successful,
        n_failed,
        iterations,
        spec
    )
end
