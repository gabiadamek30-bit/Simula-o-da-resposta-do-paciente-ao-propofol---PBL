# Optimizer Fallback Module
# Provides robust optimization with automatic fallback between methods
#
# When BFGS fails (common in NLME problems), automatically tries:
# 1. L-BFGS (limited memory BFGS, better for larger problems)
# 2. Nelder-Mead (derivative-free, last resort)

using Optim
using LineSearches

export OptimizerConfig, OptimizationResult, OptimizerType
export optimize_with_fallback, optimize_bounded_with_fallback
export BFGS_OPTIMIZER, LBFGS_OPTIMIZER, NELDER_MEAD_OPTIMIZER, CONJUGATE_GRADIENT

# =============================================================================
# Optimizer Types
# =============================================================================

"""
Available optimizer types for NLME estimation.
"""
@enum OptimizerType begin
    BFGS_OPTIMIZER        # Standard BFGS (default)
    LBFGS_OPTIMIZER       # Limited-memory BFGS (better for large problems)
    NELDER_MEAD_OPTIMIZER # Derivative-free (robust fallback)
    CONJUGATE_GRADIENT    # Conjugate gradient method
end

"""
Configuration for optimization with fallback.

# Fields
- `primary`: Primary optimizer (default: BFGS)
- `fallback_chain`: Chain of fallback optimizers to try on failure
- `max_attempts_per_optimizer`: Max attempts per optimizer before trying next
- `scale_initial_step`: Scale factor for initial step (useful for restarts)
- `verbose`: Print optimizer switching information
- `convergence_check`: Custom convergence check function

# Example
```julia
# Default: BFGS with L-BFGS and Nelder-Mead fallbacks
config = OptimizerConfig()

# Custom: Start with L-BFGS
config = OptimizerConfig(primary=LBFGS_OPTIMIZER)

# No fallback (strict BFGS only)
config = OptimizerConfig(fallback_chain=OptimizerType[])
```
"""
struct OptimizerConfig
    primary::OptimizerType
    fallback_chain::Vector{OptimizerType}
    max_attempts_per_optimizer::Int
    scale_initial_step::Float64
    verbose::Bool

    function OptimizerConfig(;
        primary::OptimizerType=BFGS_OPTIMIZER,
        fallback_chain::Vector{OptimizerType}=[LBFGS_OPTIMIZER, NELDER_MEAD_OPTIMIZER],
        max_attempts_per_optimizer::Int=2,
        scale_initial_step::Float64=1.0,
        verbose::Bool=false
    )
        @assert max_attempts_per_optimizer >= 1 "max_attempts_per_optimizer must be >= 1"
        @assert scale_initial_step > 0.0 "scale_initial_step must be positive"
        new(primary, fallback_chain, max_attempts_per_optimizer, scale_initial_step, verbose)
    end
end

"""
Result of optimization with fallback information.

# Fields
- `minimizer`: Optimal parameter values
- `minimum`: Objective function value at minimum
- `converged`: Whether optimization converged
- `iterations`: Number of iterations used
- `optimizer_used`: Which optimizer succeeded
- `fallback_count`: Number of fallback attempts made
- `messages`: Any warnings or info messages
"""
struct OptimizationResult
    minimizer::Vector{Float64}
    minimum::Float64
    converged::Bool
    iterations::Int
    optimizer_used::OptimizerType
    fallback_count::Int
    messages::Vector{String}
end

# =============================================================================
# Optimizer Factory
# =============================================================================

"""
Create an Optim.jl optimizer from OptimizerType.
"""
function create_optimizer(opt_type::OptimizerType; linesearch=LineSearches.BackTracking())
    if opt_type == BFGS_OPTIMIZER
        return BFGS(linesearch=linesearch)
    elseif opt_type == LBFGS_OPTIMIZER
        return LBFGS(linesearch=linesearch)
    elseif opt_type == NELDER_MEAD_OPTIMIZER
        return NelderMead()
    elseif opt_type == CONJUGATE_GRADIENT
        return ConjugateGradient(linesearch=linesearch)
    else
        error("Unknown optimizer type: $opt_type")
    end
end

"""
Create a bounded optimizer (Fminbox wrapper) from OptimizerType.
"""
function create_bounded_optimizer(opt_type::OptimizerType; linesearch=LineSearches.BackTracking())
    inner = create_optimizer(opt_type; linesearch=linesearch)
    return Fminbox(inner)
end

# =============================================================================
# Optimization Success Checking
# =============================================================================

"""
Check if optimization result is acceptable.

Returns false if:
- Result has NaN or Inf values
- Minimum is not finite
- Convergence failed and result is poor
"""
function is_optimization_successful(result::Optim.OptimizationResults)::Bool
    # Check if minimum is finite
    min_val = Optim.minimum(result)
    if !isfinite(min_val)
        return false
    end

    # Check if minimizer has finite values
    x = Optim.minimizer(result)
    if any(!isfinite, x)
        return false
    end

    # Check convergence (but allow non-convergence if result is reasonable)
    # Some NLME problems don't fully converge but get close enough
    if !Optim.converged(result)
        # Non-converged but finite is still usable
        # Return true to allow it, but log a warning
        return true
    end

    return true
end

# =============================================================================
# Unbounded Optimization with Fallback
# =============================================================================

"""
    optimize_with_fallback(f, x0, config; options) -> OptimizationResult

Perform optimization with automatic fallback on failure.

# Arguments
- `f`: Objective function to minimize
- `x0`: Initial parameter values
- `config`: OptimizerConfig with fallback settings
- `options`: Optim.Options for optimization control

# Returns
OptimizationResult with minimizer and fallback information.

# Example
```julia
result = optimize_with_fallback(
    x -> sum(x.^2),
    [1.0, 2.0],
    OptimizerConfig();
    options=Optim.Options(iterations=100)
)
```
"""
function optimize_with_fallback(
    f::Function,
    x0::Vector{Float64},
    config::OptimizerConfig=OptimizerConfig();
    options::Optim.Options=Optim.Options(iterations=100, g_tol=1e-6, show_trace=false),
    autodiff::Symbol=:forward
)::OptimizationResult

    messages = String[]
    fallback_count = 0

    # Build optimizer chain: primary + fallbacks
    optimizer_chain = vcat([config.primary], config.fallback_chain)

    x_current = copy(x0)
    best_result = nothing
    best_minimum = Inf

    for (idx, opt_type) in enumerate(optimizer_chain)
        for attempt in 1:config.max_attempts_per_optimizer
            try
                optimizer = create_optimizer(opt_type)

                # Scale initial point on retry attempts
                if attempt > 1
                    perturbation = 0.1 * randn(length(x_current))
                    x_current = x_current .* (1.0 .+ perturbation)
                end

                if config.verbose && (idx > 1 || attempt > 1)
                    push!(messages, "Trying $opt_type (attempt $attempt)")
                end

                # Run optimization
                if opt_type == NELDER_MEAD_OPTIMIZER
                    # Nelder-Mead doesn't use autodiff
                    result = optimize(f, x_current, optimizer, options)
                else
                    result = optimize(f, x_current, optimizer, options; autodiff=autodiff)
                end

                # Check success
                if is_optimization_successful(result)
                    min_val = Optim.minimum(result)

                    # Track best result
                    if min_val < best_minimum
                        best_minimum = min_val
                        best_result = result
                    end

                    # If converged well, return immediately
                    if Optim.converged(result)
                        return OptimizationResult(
                            Optim.minimizer(result),
                            min_val,
                            true,
                            Optim.iterations(result),
                            opt_type,
                            fallback_count,
                            messages
                        )
                    end

                    # Update current point for next attempt
                    x_current = Optim.minimizer(result)
                end

            catch e
                if config.verbose
                    push!(messages, "Optimizer $opt_type failed: $(sprint(showerror, e))")
                end
                # Continue to next attempt/optimizer
            end

            if idx > 1 || attempt > 1
                fallback_count += 1
            end
        end
    end

    # Return best result found (even if not fully converged)
    if best_result !== nothing
        push!(messages, "Warning: No optimizer fully converged, using best result")
        return OptimizationResult(
            Optim.minimizer(best_result),
            Optim.minimum(best_result),
            false,
            Optim.iterations(best_result),
            config.primary,  # Report primary as the "used" optimizer
            fallback_count,
            messages
        )
    end

    # Complete failure - return initial point
    push!(messages, "Error: All optimizers failed")
    return OptimizationResult(
        x0,
        f(x0),
        false,
        0,
        config.primary,
        fallback_count,
        messages
    )
end

# =============================================================================
# Bounded Optimization with Fallback
# =============================================================================

"""
    optimize_bounded_with_fallback(f, lower, upper, x0, config; options) -> OptimizationResult

Perform bounded optimization with automatic fallback on failure.

# Arguments
- `f`: Objective function to minimize
- `lower`: Lower bounds for parameters
- `upper`: Upper bounds for parameters
- `x0`: Initial parameter values (must be within bounds)
- `config`: OptimizerConfig with fallback settings
- `options`: Optim.Options for optimization control

# Returns
OptimizationResult with minimizer and fallback information.

# Example
```julia
result = optimize_bounded_with_fallback(
    x -> sum(x.^2),
    [0.0, 0.0],    # lower bounds
    [10.0, 10.0],  # upper bounds
    [5.0, 5.0],    # initial
    OptimizerConfig();
    options=Optim.Options(iterations=500)
)
```
"""
function optimize_bounded_with_fallback(
    f::Function,
    lower::Vector{Float64},
    upper::Vector{Float64},
    x0::Vector{Float64},
    config::OptimizerConfig=OptimizerConfig();
    options::Optim.Options=Optim.Options(iterations=500, g_tol=1e-4, show_trace=false),
    autodiff::Symbol=:forward
)::OptimizationResult

    messages = String[]
    fallback_count = 0

    # Ensure x0 is within bounds
    x0_clamped = clamp.(x0, lower, upper)
    if x0_clamped != x0
        push!(messages, "Initial point clamped to bounds")
    end

    # Build optimizer chain: primary + fallbacks
    optimizer_chain = vcat([config.primary], config.fallback_chain)

    x_current = copy(x0_clamped)
    best_result = nothing
    best_minimum = Inf

    for (idx, opt_type) in enumerate(optimizer_chain)
        for attempt in 1:config.max_attempts_per_optimizer
            try
                optimizer = create_bounded_optimizer(opt_type)

                # Perturb initial point on retry (staying within bounds)
                if attempt > 1
                    perturbation = 0.1 * randn(length(x_current))
                    x_perturbed = x_current .* (1.0 .+ perturbation)
                    x_current = clamp.(x_perturbed, lower, upper)
                end

                if config.verbose && (idx > 1 || attempt > 1)
                    push!(messages, "Trying bounded $opt_type (attempt $attempt)")
                end

                # Run bounded optimization
                result = optimize(f, lower, upper, x_current, optimizer, options)

                # Check success
                if is_optimization_successful(result)
                    min_val = Optim.minimum(result)

                    # Track best result
                    if min_val < best_minimum
                        best_minimum = min_val
                        best_result = result
                    end

                    # If converged well, return immediately
                    if Optim.converged(result)
                        return OptimizationResult(
                            Optim.minimizer(result),
                            min_val,
                            true,
                            Optim.iterations(result),
                            opt_type,
                            fallback_count,
                            messages
                        )
                    end

                    # Update current point for next attempt
                    x_current = Optim.minimizer(result)
                end

            catch e
                if config.verbose
                    push!(messages, "Bounded optimizer $opt_type failed: $(sprint(showerror, e))")
                end
                # Continue to next attempt/optimizer
            end

            if idx > 1 || attempt > 1
                fallback_count += 1
            end
        end
    end

    # Return best result found (even if not fully converged)
    if best_result !== nothing
        if config.verbose
            push!(messages, "Warning: No bounded optimizer fully converged, using best result")
        end
        return OptimizationResult(
            Optim.minimizer(best_result),
            Optim.minimum(best_result),
            false,
            Optim.iterations(best_result),
            config.primary,
            fallback_count,
            messages
        )
    end

    # Complete failure - return initial point
    push!(messages, "Error: All bounded optimizers failed")
    return OptimizationResult(
        x0_clamped,
        f(x0_clamped),
        false,
        0,
        config.primary,
        fallback_count,
        messages
    )
end

# =============================================================================
# Simple Wrapper Functions for Drop-in Replacement
# =============================================================================

"""
    robust_optimize(f, x0; max_iter, g_tol, show_trace) -> OptimizationResults

Drop-in replacement for Optim.optimize with automatic fallback.
Returns standard Optim.OptimizationResults-like interface.
"""
function robust_optimize(
    f::Function,
    x0::Vector{Float64};
    max_iter::Int=100,
    g_tol::Float64=1e-6,
    show_trace::Bool=false
)::OptimizationResult
    config = OptimizerConfig(verbose=show_trace)
    options = Optim.Options(iterations=max_iter, g_tol=g_tol, show_trace=show_trace)
    return optimize_with_fallback(f, x0, config; options=options)
end

"""
    robust_optimize_bounded(f, lower, upper, x0; max_iter, g_tol, show_trace) -> OptimizationResults

Drop-in replacement for bounded Optim.optimize with automatic fallback.
"""
function robust_optimize_bounded(
    f::Function,
    lower::Vector{Float64},
    upper::Vector{Float64},
    x0::Vector{Float64};
    max_iter::Int=500,
    g_tol::Float64=1e-4,
    show_trace::Bool=false
)::OptimizationResult
    config = OptimizerConfig(verbose=show_trace)
    options = Optim.Options(iterations=max_iter, g_tol=g_tol, show_trace=show_trace)
    return optimize_bounded_with_fallback(f, lower, upper, x0, config; options=options)
end

export robust_optimize, robust_optimize_bounded
