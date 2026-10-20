# Gradient Computation Utilities for Parameter Estimation
# Provides finite difference, ForwardDiff, and ReverseDiff support
# with automatic mode selection based on problem dimensions

using ForwardDiff
using ReverseDiff
using LinearAlgebra

export compute_gradient, compute_hessian, gradient_fd, hessian_fd
export ADMode, ForwardMode, ReverseMode, AutoMode
export gradient_with_mode, hessian_with_mode

# ============================================================================
# AD Mode Selection
# ============================================================================

"""
Abstract type for automatic differentiation mode selection.
"""
abstract type ADMode end

"""
Use forward-mode AD (ForwardDiff).
Best for: small number of parameters (n < ~100), computing Hessians.
Complexity: O(n) for gradient, O(n²) for Hessian.
"""
struct ForwardMode <: ADMode end

"""
Use reverse-mode AD (ReverseDiff).
Best for: large number of parameters, when only gradient is needed.
Complexity: O(1) for gradient (in terms of n), O(n) for Hessian.
"""
struct ReverseMode <: ADMode end

"""
Automatically select AD mode based on problem dimensions.
- Forward mode for n < threshold (default 50)
- Reverse mode for n >= threshold
"""
struct AutoMode <: ADMode
    threshold::Int
    AutoMode(threshold::Int=50) = new(threshold)
end

"""
Select AD mode based on problem dimensions.
"""
function select_mode(mode::AutoMode, n_params::Int)::ADMode
    if n_params < mode.threshold
        return ForwardMode()
    else
        return ReverseMode()
    end
end

select_mode(mode::ForwardMode, ::Int) = mode
select_mode(mode::ReverseMode, ::Int) = mode

# ============================================================================
# Finite Difference Methods (fallback)
# ============================================================================

"""
Compute gradient of a function using finite differences.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute gradient
- h: Step size for finite differences (default: 1e-6)

Returns:
- Vector of partial derivatives
"""
function gradient_fd(
    f::Function,
    x::Vector{Float64};
    h::Float64=1e-6
)::Vector{Float64}
    n = length(x)
    grad = zeros(n)

    for i in 1:n
        x_plus = copy(x)
        x_minus = copy(x)
        x_plus[i] += h
        x_minus[i] -= h

        # Central difference
        grad[i] = (f(x_plus) - f(x_minus)) / (2 * h)
    end

    return grad
end

"""
Compute Hessian of a function using finite differences.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute Hessian
- h: Step size for finite differences (default: 1e-5)

Returns:
- Symmetric Hessian matrix
"""
function hessian_fd(
    f::Function,
    x::Vector{Float64};
    h::Float64=1e-5
)::Matrix{Float64}
    n = length(x)
    hess = zeros(n, n)
    f_center = f(x)

    for i in 1:n
        for j in i:n
            if i == j
                # Diagonal: second derivative
                x_plus = copy(x)
                x_minus = copy(x)
                x_plus[i] += h
                x_minus[i] -= h

                hess[i, i] = (f(x_plus) - 2 * f_center + f(x_minus)) / h^2
            else
                # Off-diagonal: mixed partial
                x_pp = copy(x)
                x_pm = copy(x)
                x_mp = copy(x)
                x_mm = copy(x)

                x_pp[i] += h; x_pp[j] += h
                x_pm[i] += h; x_pm[j] -= h
                x_mp[i] -= h; x_mp[j] += h
                x_mm[i] -= h; x_mm[j] -= h

                hess[i, j] = (f(x_pp) - f(x_pm) - f(x_mp) + f(x_mm)) / (4 * h^2)
                hess[j, i] = hess[i, j]
            end
        end
    end

    return hess
end

# ============================================================================
# ForwardDiff Methods
# ============================================================================

"""
Compute gradient using ForwardDiff automatic differentiation.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute gradient

Returns:
- Vector of partial derivatives
"""
function gradient_forward(
    f::Function,
    x::Vector{Float64}
)::Vector{Float64}
    return ForwardDiff.gradient(f, x)
end

"""
Compute Hessian using ForwardDiff automatic differentiation.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute Hessian

Returns:
- Hessian matrix
"""
function hessian_forward(
    f::Function,
    x::Vector{Float64}
)::Matrix{Float64}
    return ForwardDiff.hessian(f, x)
end

# ============================================================================
# ReverseDiff Methods
# ============================================================================

"""
Compute gradient using ReverseDiff automatic differentiation.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute gradient

Returns:
- Vector of partial derivatives

Note: ReverseDiff is more efficient for high-dimensional problems
where n_params >> n_outputs.
"""
function gradient_reverse(
    f::Function,
    x::Vector{Float64}
)::Vector{Float64}
    return ReverseDiff.gradient(f, x)
end

"""
Compute Hessian using ReverseDiff automatic differentiation.

For Hessian computation, we use forward-over-reverse mode,
which is often more efficient than pure reverse mode.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute Hessian

Returns:
- Hessian matrix
"""
function hessian_reverse(
    f::Function,
    x::Vector{Float64}
)::Matrix{Float64}
    # For Hessian, use forward-over-reverse (most efficient)
    # grad_f returns the gradient, then we differentiate that with ForwardDiff
    grad_f = x -> ReverseDiff.gradient(f, x)
    return ForwardDiff.jacobian(grad_f, x)
end

# ============================================================================
# Mode-Aware Interface
# ============================================================================

"""
Compute gradient with specified AD mode.

Arguments:
- f: Objective function
- x: Point at which to compute gradient
- mode: AD mode (ForwardMode, ReverseMode, or AutoMode)

Returns:
- Vector of partial derivatives
"""
function gradient_with_mode(
    f::Function,
    x::Vector{Float64},
    mode::ADMode=AutoMode()
)::Vector{Float64}
    effective_mode = select_mode(mode, length(x))

    if effective_mode isa ForwardMode
        return gradient_forward(f, x)
    else
        return gradient_reverse(f, x)
    end
end

"""
Compute Hessian with specified AD mode.

Arguments:
- f: Objective function
- x: Point at which to compute Hessian
- mode: AD mode (ForwardMode, ReverseMode, or AutoMode)

Returns:
- Hessian matrix

Note: For Hessian computation, forward mode is typically used
regardless of the mode setting, as it's more efficient for
second derivatives.
"""
function hessian_with_mode(
    f::Function,
    x::Vector{Float64},
    mode::ADMode=AutoMode()
)::Matrix{Float64}
    effective_mode = select_mode(mode, length(x))

    if effective_mode isa ForwardMode
        return hessian_forward(f, x)
    else
        # Use forward-over-reverse for Hessian with large n
        return hessian_reverse(f, x)
    end
end

# ============================================================================
# Default Interface (backwards compatible)
# ============================================================================

"""
Compute gradient using ForwardDiff automatic differentiation.
Default method - uses ForwardDiff for compatibility.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute gradient

Returns:
- Vector of partial derivatives
"""
function compute_gradient(
    f::Function,
    x::Vector{Float64}
)::Vector{Float64}
    return ForwardDiff.gradient(f, x)
end

"""
Compute Hessian using ForwardDiff automatic differentiation.
Default method - uses ForwardDiff for compatibility.

Arguments:
- f: Objective function f(x) -> Float64
- x: Point at which to compute Hessian

Returns:
- Hessian matrix
"""
function compute_hessian(
    f::Function,
    x::Vector{Float64}
)::Matrix{Float64}
    return ForwardDiff.hessian(f, x)
end

# ============================================================================
# Gradient Checking and Utilities
# ============================================================================

"""
Check gradient accuracy by comparing AD and FD.
Useful for debugging.

Returns maximum relative difference between AD and FD gradients.
"""
function check_gradient(
    f::Function,
    x::Vector{Float64};
    h::Float64=1e-6,
    verbose::Bool=false,
    mode::ADMode=ForwardMode()
)::Float64
    grad_ad = gradient_with_mode(f, x, mode)
    grad_fd_val = gradient_fd(f, x; h=h)

    if verbose
        println("AD gradient:  ", grad_ad)
        println("FD gradient:  ", grad_fd_val)
        println("Difference:   ", grad_ad .- grad_fd_val)
    end

    # Relative difference
    rel_diff = abs.(grad_ad .- grad_fd_val) ./ max.(abs.(grad_ad), 1.0)
    return maximum(rel_diff)
end

export check_gradient

"""
Compare forward and reverse mode gradients.
Useful for verifying correctness.
"""
function compare_ad_modes(
    f::Function,
    x::Vector{Float64};
    verbose::Bool=false
)::Float64
    grad_forward = gradient_forward(f, x)
    grad_reverse = gradient_reverse(f, x)

    if verbose
        println("Forward gradient: ", grad_forward)
        println("Reverse gradient: ", grad_reverse)
        println("Difference:       ", grad_forward .- grad_reverse)
    end

    rel_diff = abs.(grad_forward .- grad_reverse) ./ max.(abs.(grad_forward), 1.0)
    return maximum(rel_diff)
end

export compare_ad_modes

"""
Compute Jacobian of a vector-valued function.

Arguments:
- f: Function f(x) -> Vector{Float64}
- x: Point at which to compute Jacobian

Returns:
- Jacobian matrix (m x n where m = length(f(x)), n = length(x))
"""
function compute_jacobian(
    f::Function,
    x::Vector{Float64}
)::Matrix{Float64}
    return ForwardDiff.jacobian(f, x)
end

export compute_jacobian

"""
Compute numerical condition number of a matrix.
"""
function condition_number(A::Matrix{Float64})::Float64
    eigenvalues = eigvals(A)
    real_eigenvalues = real.(eigenvalues)

    max_eig = maximum(abs.(real_eigenvalues))
    min_eig = minimum(abs.(real_eigenvalues))

    if min_eig < 1e-12
        return Inf
    end

    return max_eig / min_eig
end

export condition_number

"""
Check if matrix is positive definite.
"""
function is_positive_definite(A::Matrix{Float64})::Bool
    try
        cholesky(Symmetric(A))
        return true
    catch
        return false
    end
end

export is_positive_definite

# ============================================================================
# Precompiled Tape Support for ReverseDiff
# ============================================================================

"""
Create a compiled gradient tape for ReverseDiff.
Provides significant speedup for repeated gradient evaluations.

Arguments:
- f: Objective function
- x_template: Template input for compilation

Returns:
- Compiled gradient tape
"""
function compile_gradient_tape(f::Function, x_template::Vector{Float64})
    tape = ReverseDiff.GradientTape(f, x_template)
    return ReverseDiff.compile(tape)
end

export compile_gradient_tape

"""
Compute gradient using precompiled tape.

Arguments:
- compiled_tape: Precompiled ReverseDiff tape
- x: Point at which to compute gradient
- result: Pre-allocated result vector (optional)

Returns:
- Vector of partial derivatives
"""
function gradient_compiled!(
    result::Vector{Float64},
    compiled_tape,
    x::Vector{Float64}
)
    ReverseDiff.gradient!(result, compiled_tape, x)
    return result
end

function gradient_compiled(compiled_tape, x::Vector{Float64})::Vector{Float64}
    result = similar(x)
    return gradient_compiled!(result, compiled_tape, x)
end

export gradient_compiled!, gradient_compiled

# ============================================================================
# Benchmarking Utilities
# ============================================================================

"""
Benchmark different gradient computation methods.

Arguments:
- f: Objective function
- x: Test point
- n_runs: Number of benchmark runs

Returns:
- Dict with timing results for each method
"""
function benchmark_gradients(
    f::Function,
    x::Vector{Float64};
    n_runs::Int=100
)::Dict{String, Float64}
    results = Dict{String, Float64}()

    # Warmup
    _ = gradient_fd(f, x)
    _ = gradient_forward(f, x)
    _ = gradient_reverse(f, x)

    # Finite differences
    t_fd = @elapsed for _ in 1:n_runs
        gradient_fd(f, x)
    end
    results["finite_diff"] = t_fd / n_runs

    # ForwardDiff
    t_forward = @elapsed for _ in 1:n_runs
        gradient_forward(f, x)
    end
    results["forward_diff"] = t_forward / n_runs

    # ReverseDiff
    t_reverse = @elapsed for _ in 1:n_runs
        gradient_reverse(f, x)
    end
    results["reverse_diff"] = t_reverse / n_runs

    # ReverseDiff compiled
    compiled = compile_gradient_tape(f, x)
    t_compiled = @elapsed for _ in 1:n_runs
        gradient_compiled(compiled, x)
    end
    results["reverse_compiled"] = t_compiled / n_runs

    return results
end

export benchmark_gradients
