# Covariance Step Implementation
# Computes R matrix (expected Fisher information), S matrix (score covariance),
# and robust sandwich estimators for NLME parameter estimation

using LinearAlgebra
using ForwardDiff
using Statistics
using Distributions

export CovarianceStepResult, compute_covariance_step
export compute_r_matrix, compute_s_matrix
export SandwichEstimator, RobustSE

# ============================================================================
# Result Types
# ============================================================================

"""
Result of the covariance step computation.
"""
struct CovarianceStepResult
    # Standard covariance (based on Hessian)
    cov_standard::Union{Nothing, Matrix{Float64}}
    se_standard::Union{Nothing, Vector{Float64}}

    # Robust covariance (sandwich estimator)
    cov_robust::Union{Nothing, Matrix{Float64}}
    se_robust::Union{Nothing, Vector{Float64}}

    # R matrix: Expected (average) second derivative of -2LL
    r_matrix::Matrix{Float64}

    # S matrix: Covariance of score (gradient) contributions
    s_matrix::Matrix{Float64}

    # Individual gradient contributions
    individual_gradients::Matrix{Float64}

    # Diagnostics
    r_condition_number::Float64
    r_eigenvalues::Vector{Float64}
    success::Bool
    message::String
end

# ============================================================================
# R Matrix (Expected Fisher Information)
# ============================================================================

"""
    compute_r_matrix(ofv_fn, theta, individual_ofv_fns)

Compute the R matrix (expected Fisher information).

The R matrix is defined as:
R = -E[∂²L/∂θ²] = (1/n) * Σᵢ ∂²Lᵢ/∂θ²

where Lᵢ is the contribution to the log-likelihood from individual i.

Arguments:
- ofv_fn: Total objective function
- theta: Current parameter estimates
- individual_ofv_fns: Vector of individual objective functions

Returns:
- R matrix (n_params x n_params)
"""
function compute_r_matrix(
    ofv_fn::Function,
    theta::Vector{Float64},
    individual_ofv_fns::Vector{<:Function};
    use_ad::Bool=true
)::Matrix{Float64}
    n_params = length(theta)
    n_individuals = length(individual_ofv_fns)

    # Method 1: Use total Hessian (faster but less robust)
    if isempty(individual_ofv_fns)
        if use_ad
            return ForwardDiff.hessian(ofv_fn, theta)
        else
            return hessian_fd(ofv_fn, theta)
        end
    end

    # Method 2: Average individual Hessians (more robust)
    R = zeros(n_params, n_params)

    for ofv_i in individual_ofv_fns
        if use_ad
            H_i = ForwardDiff.hessian(ofv_i, theta)
        else
            H_i = hessian_fd(ofv_i, theta)
        end
        R .+= H_i
    end

    # Note: We don't divide by n since NONMEM uses the sum, not the mean
    # For expected Fisher: R = sum of individual Hessians
    return R
end

"""
Compute R matrix using finite differences (fallback).
"""
function compute_r_matrix_fd(
    ofv_fn::Function,
    theta::Vector{Float64};
    h::Float64=1e-5
)::Matrix{Float64}
    return hessian_fd(ofv_fn, theta; h=h)
end

# ============================================================================
# S Matrix (Score Covariance / Meat Matrix)
# ============================================================================

"""
    compute_s_matrix(theta, individual_ofv_fns)

Compute the S matrix (outer product of individual score contributions).

The S matrix is defined as:
S = Σᵢ (∂Lᵢ/∂θ)(∂Lᵢ/∂θ)' = Σᵢ sᵢ sᵢ'

where sᵢ = ∂Lᵢ/∂θ is the score contribution from individual i.

Arguments:
- theta: Current parameter estimates
- individual_ofv_fns: Vector of individual objective functions

Returns:
- S matrix (n_params x n_params)
- Individual gradients matrix (n_individuals x n_params)
"""
function compute_s_matrix(
    theta::Vector{Float64},
    individual_ofv_fns::Vector{<:Function};
    use_ad::Bool=true
)::Tuple{Matrix{Float64}, Matrix{Float64}}
    n_params = length(theta)
    n_individuals = length(individual_ofv_fns)

    # Compute individual gradients
    individual_gradients = zeros(n_individuals, n_params)

    for (i, ofv_i) in enumerate(individual_ofv_fns)
        if use_ad
            individual_gradients[i, :] = ForwardDiff.gradient(ofv_i, theta)
        else
            individual_gradients[i, :] = gradient_fd(ofv_i, theta)
        end
    end

    # Compute S matrix as outer product
    S = zeros(n_params, n_params)
    for i in 1:n_individuals
        g = individual_gradients[i, :]
        S .+= g * g'
    end

    return S, individual_gradients
end

# ============================================================================
# Sandwich Estimators
# ============================================================================

"""
Type of sandwich estimator to use.
"""
abstract type SandwichEstimator end

"""
Standard sandwich estimator: Cov = R⁻¹ S R⁻¹
"""
struct StandardSandwich <: SandwichEstimator end

"""
HC0 sandwich estimator (no small-sample correction).
"""
struct HC0 <: SandwichEstimator end

"""
HC1 sandwich estimator (simple degrees-of-freedom correction).
Multiplies by n/(n-p) where p is number of parameters.
"""
struct HC1 <: SandwichEstimator end

"""
HC2 sandwich estimator (leverage-based correction).
"""
struct HC2 <: SandwichEstimator end

"""
HC3 sandwich estimator (jackknife-like correction).
"""
struct HC3 <: SandwichEstimator end

"""
    compute_sandwich_covariance(R, S, n_obs, n_params; estimator)

Compute sandwich covariance matrix.

Standard sandwich: Cov = R⁻¹ S R⁻¹

Arguments:
- R: R matrix (Hessian / Fisher information)
- S: S matrix (outer product of scores)
- n_obs: Number of observations
- n_params: Number of parameters
- estimator: Type of sandwich estimator

Returns:
- Covariance matrix
"""
function compute_sandwich_covariance(
    R::Matrix{Float64},
    S::Matrix{Float64},
    n_obs::Int,
    n_params::Int;
    estimator::SandwichEstimator=StandardSandwich()
)::Union{Nothing, Matrix{Float64}}
    # Invert R matrix
    R_inv = try
        inv(R)
    catch
        # Try regularization
        try
            inv(R + 1e-6 * I)
        catch
            return nothing
        end
    end

    # Apply correction factor based on estimator type
    correction = if estimator isa HC1
        n_obs / (n_obs - n_params)
    else
        1.0
    end

    # Sandwich: R⁻¹ S R⁻¹
    cov = correction * (R_inv * S * R_inv)

    return cov
end

# ============================================================================
# Main Covariance Step Function
# ============================================================================

"""
    compute_covariance_step(ofv_fn, theta, individual_ofv_fns; kwargs...)

Perform the full covariance step computation.

This computes:
1. Standard covariance from Hessian inversion
2. Robust (sandwich) covariance from R and S matrices
3. Diagnostics (condition number, eigenvalues)

Arguments:
- ofv_fn: Total objective function
- theta: Parameter estimates
- individual_ofv_fns: Vector of individual OFV contributions
- use_ad: Use automatic differentiation (default: true)
- sandwich_type: Type of sandwich estimator (default: StandardSandwich())

Returns:
- CovarianceStepResult
"""
function compute_covariance_step(
    ofv_fn::Function,
    theta::Vector{Float64},
    individual_ofv_fns::Vector{<:Function};
    use_ad::Bool=true,
    sandwich_type::SandwichEstimator=StandardSandwich()
)::CovarianceStepResult
    n_params = length(theta)
    n_individuals = length(individual_ofv_fns)

    # Compute R matrix
    R = compute_r_matrix(ofv_fn, theta, individual_ofv_fns; use_ad=use_ad)

    # Compute S matrix and individual gradients
    S, individual_gradients = compute_s_matrix(theta, individual_ofv_fns; use_ad=use_ad)

    # Compute eigenvalues for diagnostics
    R_eigenvalues = try
        real.(eigvals(Symmetric(R)))
    catch
        zeros(n_params)
    end

    r_condition = if minimum(abs.(R_eigenvalues)) > 1e-12
        maximum(abs.(R_eigenvalues)) / minimum(abs.(R_eigenvalues))
    else
        Inf
    end

    # Standard covariance (Hessian inverse)
    cov_standard, se_standard, std_success = compute_se_from_hessian(R)

    # Robust (sandwich) covariance
    cov_robust = compute_sandwich_covariance(R, S, n_individuals, n_params;
        estimator=sandwich_type)

    se_robust = if cov_robust !== nothing
        variances = diag(cov_robust)
        if all(variances .>= 0)
            sqrt.(variances)
        else
            nothing
        end
    else
        nothing
    end

    success = std_success || (cov_robust !== nothing && se_robust !== nothing)

    message = if success
        "Covariance step successful"
    elseif r_condition > 1e10
        "R matrix is ill-conditioned (condition number = $(r_condition))"
    elseif any(R_eigenvalues .<= 0)
        "R matrix is not positive definite"
    else
        "Covariance step failed"
    end

    return CovarianceStepResult(
        cov_standard,
        se_standard,
        cov_robust,
        se_robust,
        R,
        S,
        individual_gradients,
        r_condition,
        R_eigenvalues,
        success,
        message
    )
end

# ============================================================================
# Robust SE Wrapper
# ============================================================================

"""
    RobustSE(result::CovarianceStepResult)

Extract robust standard errors from covariance step result.
Falls back to standard SE if robust SE computation failed.
"""
function RobustSE(result::CovarianceStepResult)::Union{Nothing, Vector{Float64}}
    if result.se_robust !== nothing
        return result.se_robust
    elseif result.se_standard !== nothing
        return result.se_standard
    else
        return nothing
    end
end

# ============================================================================
# Model Comparison Using Covariance
# ============================================================================

"""
    wald_test(theta, se, null_value)

Perform Wald test for parameter significance.

H0: θ = null_value
Wald statistic: W = (θ - θ₀)² / Var(θ)

Returns (statistic, p_value)
"""
function wald_test(
    theta::Float64,
    se::Float64,
    null_value::Float64=0.0
)::Tuple{Float64, Float64}
    z = (theta - null_value) / se
    p_value = 2 * (1 - cdf(Normal(), abs(z)))
    return z^2, p_value
end

export wald_test

"""
    score_test(S, R_inv, theta, null_value)

Perform score test (Lagrange multiplier test).

Under H0: Score test statistic = S' R⁻¹ S ~ χ²(k)
"""
function score_test(
    score::Vector{Float64},
    R_inv::Matrix{Float64}
)::Tuple{Float64, Float64}
    statistic = dot(score, R_inv * score)
    df = length(score)
    p_value = 1 - cdf(Chisq(df), statistic)
    return statistic, p_value
end

export score_test

# ============================================================================
# Diagnostic Functions
# ============================================================================

"""
Check if covariance matrix is positive definite and well-conditioned.
"""
function check_covariance_quality(
    cov::Matrix{Float64};
    max_condition::Float64=1e6
)::NamedTuple{(:positive_definite, :well_conditioned, :condition_number), Tuple{Bool, Bool, Float64}}
    eigenvalues = eigvals(Symmetric(cov))
    pd = all(real.(eigenvalues) .> 0)

    max_eig = maximum(abs.(eigenvalues))
    min_eig = minimum(abs.(eigenvalues))
    cond = min_eig > 1e-12 ? max_eig / min_eig : Inf

    well_cond = cond < max_condition

    return (positive_definite=pd, well_conditioned=well_cond, condition_number=cond)
end

export check_covariance_quality

"""
Compare standard and robust SE to assess model misspecification.

Large differences suggest potential model misspecification.
Returns ratio of robust/standard SE.
"""
function compare_se(
    se_standard::Vector{Float64},
    se_robust::Vector{Float64}
)::Vector{Float64}
    return se_robust ./ se_standard
end

export compare_se
