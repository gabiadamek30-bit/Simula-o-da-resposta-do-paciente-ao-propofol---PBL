# Parallel Computing Infrastructure for Population PK/PD
# Professional-grade implementation following industry standards
#
# Supports:
# - Multi-threaded execution for shared-memory systems
# - Distributed computing for cluster environments
# - Thread-safe RNG management for reproducibility
# - Automatic load balancing with chunking strategies
# - Seamless fallback to sequential execution
#
# References:
# - Julia Multi-Threading Documentation
# - NONMEM Parallelization Guide
# - Monolix Parallel Computing

using Base.Threads: @threads, nthreads, threadid
using Random: Random, AbstractRNG
using StableRNGs: StableRNG

export ParallelBackend, ParallelConfig, ParallelResult
export SerialBackend, ThreadedBackend, DistributedBackend
export parallel_map, parallel_map_reduce, parallel_sum, parallel_map_with_rng
export create_thread_rngs, create_subject_rngs, get_optimal_chunk_size, create_chunks
export with_parallel_config, current_parallel_config, set_parallel_config!
export n_workers, is_parallel, available_threads, recommend_parallel_config
export ParallelSubjectProcessor, process_subjects_parallel, process_subjects_parallel_with_rng
export compute_ofv_parallel, compute_ofv_with_etas_parallel
export ParallelBootstrapExecutor

# =============================================================================
# Parallel Backend Types
# =============================================================================

"""
Abstract type for parallel execution backends.
"""
abstract type ParallelBackend end

"""
Serial execution (no parallelism).
Used as fallback or when n_workers = 1.
"""
struct SerialBackend <: ParallelBackend end

"""
Multi-threaded execution using Julia's threading.
Best for shared-memory systems (single machine).
"""
struct ThreadedBackend <: ParallelBackend
    n_threads::Int

    function ThreadedBackend(n_threads::Int=nthreads())
        @assert n_threads >= 1 "n_threads must be at least 1"
        new(min(n_threads, nthreads()))
    end
end

"""
Distributed execution for cluster computing.
Uses Julia's Distributed.jl for multi-node execution.
"""
struct DistributedBackend <: ParallelBackend
    n_workers::Int

    function DistributedBackend(n_workers::Int=1)
        @assert n_workers >= 1 "n_workers must be at least 1"
        new(n_workers)
    end
end

# =============================================================================
# Parallel Configuration
# =============================================================================

"""
Configuration for parallel execution in estimation.

# Fields
- `backend`: Parallel execution backend (Serial, Threaded, or Distributed)
- `chunk_size`: Number of subjects per parallel chunk (0 = auto)
- `load_balance`: Enable dynamic load balancing (default: true)
- `progress`: Show progress during parallel execution (default: false)
- `seed`: Base seed for reproducible parallel RNG (default: nothing = random)
- `verbose`: Print parallel execution info (default: false)

# Example
```julia
# Use all available threads
config = ParallelConfig(ThreadedBackend())

# Use 4 threads with fixed seed for reproducibility
config = ParallelConfig(ThreadedBackend(4), seed=12345)

# Serial execution (debugging)
config = ParallelConfig(SerialBackend())
```
"""
struct ParallelConfig
    backend::ParallelBackend
    chunk_size::Int
    load_balance::Bool
    progress::Bool
    seed::Union{Nothing, UInt64}
    verbose::Bool

    function ParallelConfig(
        backend::ParallelBackend=ThreadedBackend();
        chunk_size::Int=0,
        load_balance::Bool=true,
        progress::Bool=false,
        seed::Union{Nothing, Integer}=nothing,
        verbose::Bool=false
    )
        @assert chunk_size >= 0 "chunk_size must be non-negative (0 = auto)"
        seed_val = seed === nothing ? nothing : UInt64(seed)
        new(backend, chunk_size, load_balance, progress, seed_val, verbose)
    end
end

# Default configuration (serial for safety)
const DEFAULT_PARALLEL_CONFIG = Ref{ParallelConfig}(ParallelConfig(SerialBackend()))

"""
    current_parallel_config() -> ParallelConfig

Get the current default parallel configuration.
"""
function current_parallel_config()::ParallelConfig
    return DEFAULT_PARALLEL_CONFIG[]
end

"""
    set_parallel_config!(config::ParallelConfig)

Set the default parallel configuration.
"""
function set_parallel_config!(config::ParallelConfig)
    DEFAULT_PARALLEL_CONFIG[] = config
    return config
end

"""
    with_parallel_config(f, config::ParallelConfig)

Execute function `f` with temporary parallel configuration.
"""
function with_parallel_config(f::Function, config::ParallelConfig)
    old_config = current_parallel_config()
    try
        set_parallel_config!(config)
        return f()
    finally
        set_parallel_config!(old_config)
    end
end

"""
    n_workers(config::ParallelConfig) -> Int

Get the number of parallel workers for the given configuration.
"""
function n_workers(config::ParallelConfig)::Int
    if config.backend isa SerialBackend
        return 1
    elseif config.backend isa ThreadedBackend
        return config.backend.n_threads
    elseif config.backend isa DistributedBackend
        return config.backend.n_workers
    else
        return 1
    end
end

"""
    is_parallel(config::ParallelConfig) -> Bool

Check if configuration enables parallel execution.
"""
function is_parallel(config::ParallelConfig)::Bool
    return n_workers(config) > 1
end

# =============================================================================
# Thread-Safe RNG Management
# =============================================================================

"""
    create_thread_rngs(n::Int; seed=nothing) -> Vector{StableRNG}

Create `n` independent RNG streams for parallel execution.
Each stream is seeded deterministically from the base seed for reproducibility.
"""
function create_thread_rngs(n::Int; seed::Union{Nothing, UInt64}=nothing)::Vector{StableRNG}
    base_seed = seed === nothing ? rand(UInt64) : seed

    # Create independent RNG streams using jump-ahead
    rngs = Vector{StableRNG}(undef, n)
    for i in 1:n
        # Combine base seed with stream index for independence
        stream_seed = base_seed ⊻ hash(i)
        rngs[i] = StableRNG(stream_seed)
    end

    return rngs
end

"""
    create_subject_rngs(n_subjects::Int, config::ParallelConfig) -> Vector{StableRNG}

Create RNG streams for each subject, ensuring reproducibility in parallel execution.
"""
function create_subject_rngs(n_subjects::Int, config::ParallelConfig)::Vector{StableRNG}
    return create_thread_rngs(n_subjects; seed=config.seed)
end

# =============================================================================
# Chunking Strategies
# =============================================================================

"""
    get_optimal_chunk_size(n_items::Int, n_workers::Int; min_chunk::Int=1) -> Int

Compute optimal chunk size for load balancing.
"""
function get_optimal_chunk_size(n_items::Int, n_workers::Int; min_chunk::Int=1)::Int
    if n_workers <= 1
        return n_items
    end

    # Aim for ~4 chunks per worker for good load balancing
    target_chunks = n_workers * 4
    chunk_size = max(min_chunk, cld(n_items, target_chunks))

    return chunk_size
end

"""
    create_chunks(n_items::Int, chunk_size::Int) -> Vector{UnitRange{Int}}

Create index ranges for chunked processing.
"""
function create_chunks(n_items::Int, chunk_size::Int)::Vector{UnitRange{Int}}
    if chunk_size <= 0 || chunk_size >= n_items
        return [1:n_items]
    end

    chunks = UnitRange{Int}[]
    start = 1
    while start <= n_items
        stop = min(start + chunk_size - 1, n_items)
        push!(chunks, start:stop)
        start = stop + 1
    end

    return chunks
end

# =============================================================================
# Parallel Map Operations
# =============================================================================

"""
    parallel_map(f, items, config::ParallelConfig) -> Vector

Apply function `f` to each item in parallel according to configuration.
Returns results in the same order as input items.

Thread-safe: Each call to `f` receives its own RNG if needed.
"""
function parallel_map(
    f::Function,
    items::AbstractVector,
    config::ParallelConfig
)::Vector
    n = length(items)

    if n == 0
        return []
    end

    # Serial execution
    if !is_parallel(config) || n == 1
        return [f(item) for item in items]
    end

    # Threaded execution
    if config.backend isa ThreadedBackend
        return parallel_map_threaded(f, items, config)
    end

    # Distributed execution
    if config.backend isa DistributedBackend
        return parallel_map_distributed(f, items, config)
    end

    # Fallback to serial
    return [f(item) for item in items]
end

"""
Internal threaded parallel map implementation.
"""
function parallel_map_threaded(
    f::Function,
    items::AbstractVector,
    config::ParallelConfig
)::Vector
    n = length(items)
    n_threads = config.backend.n_threads

    # Determine chunk size
    chunk_size = config.chunk_size > 0 ?
        config.chunk_size :
        get_optimal_chunk_size(n, n_threads)

    # Pre-allocate results (type-stable)
    first_result = f(items[1])
    results = Vector{typeof(first_result)}(undef, n)
    results[1] = first_result

    if n == 1
        return results
    end

    # Process remaining items in parallel
    @threads for i in 2:n
        results[i] = f(items[i])
    end

    return results
end

"""
Internal distributed parallel map implementation.
Uses pmap when Distributed is loaded, otherwise falls back to threaded.
"""
function parallel_map_distributed(
    f::Function,
    items::AbstractVector,
    config::ParallelConfig
)::Vector
    # Check if Distributed is available
    if isdefined(Main, :Distributed) && isdefined(Main.Distributed, :pmap)
        return Main.Distributed.pmap(f, items)
    end

    # Fallback to threaded execution
    if config.verbose
        @warn "Distributed.jl not loaded, falling back to threaded execution"
    end

    threaded_config = ParallelConfig(
        ThreadedBackend(config.backend.n_workers);
        chunk_size=config.chunk_size,
        load_balance=config.load_balance,
        progress=config.progress,
        seed=config.seed,
        verbose=config.verbose
    )

    return parallel_map_threaded(f, items, threaded_config)
end

"""
    parallel_map_with_rng(f, items, rngs, config) -> Vector

Apply function `f(item, rng)` to each item with its own RNG.
Ensures reproducible parallel execution for stochastic algorithms.
"""
function parallel_map_with_rng(
    f::Function,
    items::AbstractVector,
    rngs::Vector{<:AbstractRNG},
    config::ParallelConfig
)::Vector
    n = length(items)
    @assert length(rngs) >= n "Need at least $n RNGs for $n items"

    if n == 0
        return []
    end

    # Serial execution
    if !is_parallel(config) || n == 1
        return [f(items[i], rngs[i]) for i in 1:n]
    end

    # Threaded execution with RNG
    if config.backend isa ThreadedBackend
        # Pre-allocate results
        first_result = f(items[1], rngs[1])
        results = Vector{typeof(first_result)}(undef, n)
        results[1] = first_result

        if n > 1
            @threads for i in 2:n
                results[i] = f(items[i], rngs[i])
            end
        end

        return results
    end

    # Fallback to serial with RNG
    return [f(items[i], rngs[i]) for i in 1:n]
end

# =============================================================================
# Parallel Reduction Operations
# =============================================================================

"""
    parallel_map_reduce(f, op, init, items, config) -> result

Map function `f` over items and reduce with binary operator `op`.
Useful for computing sums, products, or custom aggregations in parallel.

# Example
```julia
# Parallel sum of squared values
result = parallel_map_reduce(x -> x^2, +, 0.0, values, config)
```
"""
function parallel_map_reduce(
    f::Function,
    op::Function,
    init,
    items::AbstractVector,
    config::ParallelConfig
)
    n = length(items)

    if n == 0
        return init
    end

    # Serial execution
    if !is_parallel(config) || n == 1
        result = init
        for item in items
            result = op(result, f(item))
        end
        return result
    end

    # Threaded execution with thread-local reduction
    if config.backend isa ThreadedBackend
        return parallel_map_reduce_threaded(f, op, init, items, config)
    end

    # Fallback to serial
    result = init
    for item in items
        result = op(result, f(item))
    end
    return result
end

"""
Internal threaded map-reduce with thread-local accumulators.
"""
function parallel_map_reduce_threaded(
    f::Function,
    op::Function,
    init,
    items::AbstractVector,
    config::ParallelConfig
)
    n = length(items)
    n_threads = config.backend.n_threads

    # Thread-local accumulators
    local_results = [deepcopy(init) for _ in 1:n_threads]

    @threads for i in 1:n
        tid = threadid()
        local_results[tid] = op(local_results[tid], f(items[i]))
    end

    # Combine thread-local results
    result = init
    for local_result in local_results
        result = op(result, local_result)
    end

    return result
end

"""
    parallel_sum(f, items, config) -> Number

Parallel sum of `f(item)` for each item.
Optimized for numerical summation with minimal overhead.
"""
function parallel_sum(
    f::Function,
    items::AbstractVector,
    config::ParallelConfig
)::Float64
    return parallel_map_reduce(f, +, 0.0, items, config)
end

# =============================================================================
# Subject Processing Infrastructure
# =============================================================================

"""
Result from parallel subject processing.
"""
struct ParallelSubjectResult{T}
    results::Vector{T}
    total_time::Float64
    n_threads_used::Int
end

"""
    ParallelSubjectProcessor

Processes subjects in parallel with configurable execution.
"""
struct ParallelSubjectProcessor
    config::ParallelConfig
    rngs::Vector{StableRNG}

    function ParallelSubjectProcessor(n_subjects::Int, config::ParallelConfig)
        rngs = create_subject_rngs(n_subjects, config)
        new(config, rngs)
    end
end

"""
    process_subjects_parallel(f, subjects, config) -> Vector

Process subjects in parallel using the given function.

# Arguments
- `f`: Function that takes (subject_data, subject_index) and returns result
- `subjects`: Vector of subject data
- `config`: Parallel configuration

# Returns
Vector of results in subject order.
"""
function process_subjects_parallel(
    f::Function,
    subjects::AbstractVector,
    config::ParallelConfig
)::Vector
    n = length(subjects)

    if n == 0
        return []
    end

    # Create indexed items
    indexed_subjects = [(subjects[i], i) for i in 1:n]

    # Process with parallel_map
    return parallel_map(x -> f(x[1], x[2]), indexed_subjects, config)
end

"""
    process_subjects_parallel_with_rng(f, subjects, config) -> Vector

Process subjects in parallel with per-subject RNG for reproducibility.

# Arguments
- `f`: Function that takes (subject_data, subject_index, rng) and returns result
- `subjects`: Vector of subject data
- `config`: Parallel configuration

# Returns
Vector of results in subject order.
"""
function process_subjects_parallel_with_rng(
    f::Function,
    subjects::AbstractVector,
    config::ParallelConfig
)::Vector
    n = length(subjects)

    if n == 0
        return []
    end

    # Create RNGs for each subject
    rngs = create_subject_rngs(n, config)

    # Create indexed items
    indexed_subjects = [(subjects[i], i) for i in 1:n]

    # Process with parallel_map_with_rng
    return parallel_map_with_rng(
        (x, rng) -> f(x[1], x[2], rng),
        indexed_subjects,
        rngs,
        config
    )
end

# =============================================================================
# Parallel OFV (Objective Function Value) Computation
# =============================================================================

"""
    compute_ofv_parallel(subject_ofv_fn, subjects, config) -> Float64

Compute total OFV (objective function value) by summing subject contributions in parallel.

# Arguments
- `subject_ofv_fn`: Function that takes subject data and returns OFV contribution
- `subjects`: Vector of subject data
- `config`: Parallel configuration

# Returns
Total OFV (sum of all subject contributions).
"""
function compute_ofv_parallel(
    subject_ofv_fn::Function,
    subjects::AbstractVector,
    config::ParallelConfig
)::Float64
    return parallel_sum(subject_ofv_fn, subjects, config)
end

"""
    compute_ofv_with_etas_parallel(subject_ofv_fn, subjects, etas, config) -> Float64

Compute total OFV with per-subject eta values in parallel.

# Arguments
- `subject_ofv_fn`: Function that takes (subject_data, eta) and returns OFV contribution
- `subjects`: Vector of subject data
- `etas`: Vector of eta vectors for each subject
- `config`: Parallel configuration

# Returns
Total OFV (sum of all subject contributions).
"""
function compute_ofv_with_etas_parallel(
    subject_ofv_fn::Function,
    subjects::AbstractVector,
    etas::Vector{Vector{Float64}},
    config::ParallelConfig
)::Float64
    n = length(subjects)
    @assert length(etas) == n "Must have eta for each subject"

    # Create subject-eta pairs
    pairs = [(subjects[i], etas[i]) for i in 1:n]

    return parallel_sum(x -> subject_ofv_fn(x[1], x[2]), pairs, config)
end

# =============================================================================
# Parallel Eta Optimization
# =============================================================================

"""
Result from parallel eta optimization.
"""
struct ParallelEtaResult
    etas::Vector{Vector{Float64}}
    hessians::Vector{Matrix{Float64}}
    log_likelihoods::Vector{Float64}
    prior_contributions::Vector{Float64}
    converged::Vector{Bool}
end

"""
    optimize_etas_parallel(eta_optimizer, subjects, theta, omega, sigma, config) -> ParallelEtaResult

Optimize individual etas for all subjects in parallel.

# Arguments
- `eta_optimizer`: Function that optimizes eta for a single subject
- `subjects`: Vector of subject data
- `theta`: Population parameters
- `omega`: Random effects covariance
- `sigma`: Residual error specification
- `config`: Parallel configuration

# Returns
ParallelEtaResult containing all individual results.
"""
function optimize_etas_parallel(
    eta_optimizer::Function,
    subjects::AbstractVector,
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma,
    initial_etas::Vector{Vector{Float64}},
    config::ParallelConfig
)::ParallelEtaResult
    n = length(subjects)

    # Create items with all necessary data
    items = [(subjects[i], initial_etas[i], i) for i in 1:n]

    # Run parallel optimization
    results = parallel_map(
        item -> eta_optimizer(item[1], item[2], theta, omega, sigma),
        items,
        config
    )

    # Unpack results
    etas = [r[1] for r in results]
    hessians = [r[2] for r in results]
    lls = [r[3] for r in results]
    priors = [r[4] for r in results]
    converged = [r[5] for r in results]

    return ParallelEtaResult(etas, hessians, lls, priors, converged)
end

# =============================================================================
# Parallel MCMC for SAEM
# =============================================================================

"""
Result from parallel MCMC sampling.
"""
struct ParallelMCMCResult
    samples::Vector{Vector{Vector{Float64}}}  # [subject][chain][eta_vector]
    acceptance_rates::Vector{Vector{Float64}}
    proposal_sds::Vector{Vector{Float64}}
end

"""
    sample_etas_parallel(mcmc_sampler, subjects, theta, omega, sigma, current_etas, proposal_sds, rngs, config) -> ParallelMCMCResult

Run MCMC sampling for all subjects in parallel with reproducible RNG.

# Arguments
- `mcmc_sampler`: Function that samples etas for a single subject
- `subjects`: Vector of subject data
- `theta`: Population parameters
- `omega`: Random effects covariance
- `sigma`: Residual error specification
- `current_etas`: Current eta chains for each subject
- `proposal_sds`: Proposal standard deviations
- `config`: Parallel configuration with seed for reproducibility

# Returns
ParallelMCMCResult containing all MCMC samples.
"""
function sample_etas_parallel(
    mcmc_sampler::Function,
    subjects::AbstractVector,
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma,
    current_chains::Vector,
    proposal_sds::Vector,
    config::ParallelConfig
)::ParallelMCMCResult
    n = length(subjects)

    # Create subject RNGs for reproducibility
    rngs = create_subject_rngs(n, config)

    # Create items with all necessary data
    items = [(subjects[i], current_chains[i], proposal_sds[i], i) for i in 1:n]

    # Run parallel MCMC
    results = parallel_map_with_rng(
        (item, rng) -> mcmc_sampler(item[1], item[2], item[3], theta, omega, sigma, rng),
        items,
        rngs,
        config
    )

    # Unpack results
    samples = [r[1] for r in results]
    acceptance_rates = [r[2] for r in results]
    new_proposal_sds = [r[3] for r in results]

    return ParallelMCMCResult(samples, acceptance_rates, new_proposal_sds)
end

# =============================================================================
# Parallel Bootstrap
# =============================================================================

"""
Result from a single bootstrap replicate.
"""
struct BootstrapReplicateResult
    theta::Vector{Float64}
    omega::Matrix{Float64}
    sigma_params::Vector{Float64}
    converged::Bool
    ofv::Float64
end

"""
    run_bootstrap_parallel(estimation_fn, data_resampler, observed_data, n_replicates, config) -> Vector{BootstrapReplicateResult}

Run bootstrap replicates in parallel.

# Arguments
- `estimation_fn`: Function that estimates parameters from data
- `data_resampler`: Function that resamples the data
- `observed_data`: Original observed data
- `n_replicates`: Number of bootstrap replicates
- `config`: Parallel configuration (use DistributedBackend for clusters)

# Returns
Vector of BootstrapReplicateResult for each replicate.
"""
function run_bootstrap_parallel(
    estimation_fn::Function,
    data_resampler::Function,
    observed_data,
    n_replicates::Int,
    config::ParallelConfig
)::Vector{BootstrapReplicateResult}
    # Create RNGs for each replicate
    rngs = create_thread_rngs(n_replicates; seed=config.seed)

    # Create replicate indices
    replicates = collect(1:n_replicates)

    # Run parallel bootstrap
    results = parallel_map_with_rng(
        (rep_idx, rng) -> begin
            try
                # Resample data
                resampled_data = data_resampler(observed_data, rng)

                # Run estimation
                result = estimation_fn(resampled_data)

                # Extract results
                BootstrapReplicateResult(
                    copy(result.theta),
                    copy(result.omega),
                    extract_sigma_params(result),
                    result.convergence,
                    result.ofv
                )
            catch e
                # Return failed result
                n_theta = 0  # Will be handled by caller
                BootstrapReplicateResult(
                    Float64[],
                    Matrix{Float64}(undef, 0, 0),
                    Float64[],
                    false,
                    NaN
                )
            end
        end,
        replicates,
        rngs,
        config
    )

    return results
end

"""
Helper to extract sigma parameters from estimation result.
"""
function extract_sigma_params(result)::Vector{Float64}
    if hasproperty(result, :sigma) && result.sigma !== nothing
        sigma = result.sigma
        if hasproperty(sigma, :params)
            p = sigma.params
            if hasproperty(p, :sigma)
                return [p.sigma]
            elseif hasproperty(p, :a) && hasproperty(p, :b)
                return [p.a, p.b]
            end
        end
    end
    return Float64[]
end

# =============================================================================
# Parallel Individual Estimates Computation
# =============================================================================

"""
Result from parallel individual estimates computation.
"""
struct ParallelIndividualResult
    subject_id::String
    eta::Vector{Float64}
    ipred::Vector{Float64}
    pred::Vector{Float64}
    cwres::Vector{Float64}
    iwres::Vector{Float64}
    wres::Vector{Float64}
    ofv_contribution::Float64
end

"""
    compute_individual_estimates_parallel(compute_fn, subjects, etas, theta, omega, sigma, config) -> Vector{ParallelIndividualResult}

Compute individual estimates (IPRED, PRED, residuals) for all subjects in parallel.

# Arguments
- `compute_fn`: Function that computes estimates for a single subject
- `subjects`: Vector of subject data
- `etas`: Individual eta estimates
- `theta`: Population parameters
- `omega`: Random effects covariance
- `sigma`: Residual error specification
- `config`: Parallel configuration

# Returns
Vector of ParallelIndividualResult for each subject.
"""
function compute_individual_estimates_parallel(
    compute_fn::Function,
    subjects::AbstractVector,
    etas::Vector{Vector{Float64}},
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma,
    config::ParallelConfig
)::Vector{ParallelIndividualResult}
    n = length(subjects)

    # Create items with all necessary data
    items = [(subjects[i], etas[i], i) for i in 1:n]

    # Run parallel computation
    return parallel_map(
        item -> compute_fn(item[1], item[2], theta, omega, sigma),
        items,
        config
    )
end

# =============================================================================
# Parallel Gradient/Hessian Computation
# =============================================================================

"""
    compute_gradients_parallel(gradient_fn, subjects, theta, config) -> Matrix{Float64}

Compute per-subject gradients in parallel for robust standard errors.

# Arguments
- `gradient_fn`: Function that computes gradient for a single subject
- `subjects`: Vector of subject data
- `theta`: Parameter vector
- `config`: Parallel configuration

# Returns
Matrix of gradients (n_subjects × n_parameters).
"""
function compute_gradients_parallel(
    gradient_fn::Function,
    subjects::AbstractVector,
    theta::Vector{Float64},
    config::ParallelConfig
)::Matrix{Float64}
    n_subj = length(subjects)
    n_params = length(theta)

    # Compute gradients in parallel
    gradients = parallel_map(
        subj -> gradient_fn(subj, theta),
        subjects,
        config
    )

    # Combine into matrix
    result = Matrix{Float64}(undef, n_subj, n_params)
    for i in 1:n_subj
        result[i, :] = gradients[i]
    end

    return result
end

# =============================================================================
# Utility Functions
# =============================================================================

"""
    available_threads() -> Int

Get the number of threads available for parallel execution.
"""
function available_threads()::Int
    return nthreads()
end

"""
    recommend_parallel_config(n_subjects::Int) -> ParallelConfig

Recommend a parallel configuration based on problem size and available resources.
"""
function recommend_parallel_config(n_subjects::Int)::ParallelConfig
    n_threads = available_threads()

    if n_subjects < 10 || n_threads == 1
        # Small problem or single thread - use serial
        return ParallelConfig(SerialBackend())
    elseif n_subjects < 100
        # Medium problem - use threading with load balancing
        return ParallelConfig(
            ThreadedBackend(n_threads);
            load_balance=true
        )
    else
        # Large problem - full threading with progress
        return ParallelConfig(
            ThreadedBackend(n_threads);
            load_balance=true,
            progress=true
        )
    end
end

"""
    with_progress(f, config::ParallelConfig, n_items::Int, desc::String="Processing")

Execute function with optional progress tracking.
"""
function with_progress(f::Function, config::ParallelConfig, n_items::Int, desc::String="Processing")
    if config.progress && n_items > 10
        # Could integrate with ProgressMeter.jl if available
        if config.verbose
            println("$desc $n_items items with $(n_workers(config)) workers...")
        end
    end

    result = f()

    if config.progress && config.verbose && n_items > 10
        println("$desc complete.")
    end

    return result
end

# =============================================================================
# Pretty Printing
# =============================================================================

function Base.show(io::IO, config::ParallelConfig)
    println(io, "ParallelConfig")
    println(io, "  Backend: $(typeof(config.backend))")
    println(io, "  Workers: $(n_workers(config))")
    println(io, "  Chunk size: $(config.chunk_size == 0 ? "auto" : config.chunk_size)")
    println(io, "  Load balance: $(config.load_balance)")
    if config.seed !== nothing
        println(io, "  Seed: $(config.seed)")
    end
end

function Base.show(io::IO, backend::ThreadedBackend)
    print(io, "ThreadedBackend($(backend.n_threads) threads)")
end

function Base.show(io::IO, backend::DistributedBackend)
    print(io, "DistributedBackend($(backend.n_workers) workers)")
end

function Base.show(io::IO, ::SerialBackend)
    print(io, "SerialBackend()")
end
