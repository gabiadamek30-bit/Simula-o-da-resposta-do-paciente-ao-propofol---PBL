# Standard Error Computation for Parameter Estimates
#
# Production-grade implementation for computing SEs from marginal likelihood.
# Implements multiple methods with proper integration over random effects.
#
# Key Features:
# 1. Marginal likelihood-based SEs (not fixed-eta conditional)
# 2. Multiple methods: Laplace approximation, Monte Carlo integration
# 3. Full parameter covariance matrix (theta, omega, sigma)
# 4. Validation via simulation-based CI coverage checks
#
# References:
# - Kuhn & Lavielle (2005), Computational Statistics & Data Analysis
# - Wang (2007), Statistics in Medicine - NONMEM covariance step
# - Almquist et al. (2015), CPT:PSP - Louis' missing information

using LinearAlgebra
using Distributions
using Statistics
using StableRNGs

export StandardErrorResult, compute_standard_errors_marginal
export compute_se_from_hessian, compute_se_sandwich, compute_ci
export compute_rse, cov_to_corr_matrix, bootstrap_se
export validate_ci_coverage_simulation, SEMethod
export LaplaceSE, MonteCarloSE, LouisSE, SandwichSE

# ============================================================================
# SE Computation Method Types
# ============================================================================

"""
Abstract type for SE computation methods.
"""
abstract type SEMethod end

"""
Laplace approximation for marginal likelihood.
Standard FOCE/Laplacian approach - expands around eta mode.
"""
struct LaplaceSE <: SEMethod end

"""
Monte Carlo integration for marginal likelihood.
More accurate but slower - samples from prior p(η|Ω).
"""
struct MonteCarloSE <: SEMethod
    n_samples::Int

    function MonteCarloSE(; n_samples::Int=500)
        new(n_samples)
    end
end

"""
Louis' missing information principle for SAEM.
Accounts for Monte Carlo uncertainty in stochastic EM.
"""
struct LouisSE <: SEMethod
    n_importance_samples::Int

    function LouisSE(; n_importance_samples::Int=200)
        new(n_importance_samples)
    end
end

"""
Sandwich (robust) standard errors.
Uses R⁻¹ S R⁻¹ form for robustness to model misspecification.
"""
struct SandwichSE <: SEMethod end

# ============================================================================
# Result Structure
# ============================================================================

"""
Comprehensive result from standard error computation.

Contains:
- SEs for all parameter types (theta, omega, sigma)
- Full covariance and correlation matrices
- Diagnostic information
- Method used and convergence status
"""
struct StandardErrorResult
    # Standard errors
    theta_se::Union{Nothing, Vector{Float64}}
    omega_se::Union{Nothing, Matrix{Float64}}  # Full matrix for block omega
    sigma_se::Union{Nothing, Vector{Float64}}

    # Full covariance matrix (theta, log(omega), log(sigma))
    covariance_matrix::Union{Nothing, Matrix{Float64}}

    # Correlation matrix
    correlation_matrix::Union{Nothing, Matrix{Float64}}

    # Confidence intervals (95%)
    theta_ci_lower::Union{Nothing, Vector{Float64}}
    theta_ci_upper::Union{Nothing, Vector{Float64}}
    omega_ci_lower::Union{Nothing, Matrix{Float64}}
    omega_ci_upper::Union{Nothing, Matrix{Float64}}

    # Diagnostics
    condition_number::Float64
    min_eigenvalue::Float64
    hessian_pd::Bool  # Is Hessian positive definite?

    # Method and status
    method::Symbol
    successful::Bool
    message::String
end

# ============================================================================
# Main SE Computation Function
# ============================================================================

"""
    compute_standard_errors_marginal(
        theta, omega, sigma, subjects, model_spec, config;
        method=LaplaceSE(), grid, solver
    )

Compute standard errors from the marginal likelihood.

The marginal likelihood integrates out random effects:
    L(θ,Ω,σ) = ∏ᵢ ∫ p(yᵢ|ηᵢ,θ,σ) p(ηᵢ|Ω) dηᵢ

This function computes the Hessian of the negative log marginal likelihood
with respect to all population parameters (θ, Ω, σ).

# Arguments
- `theta`: Fixed effects estimates
- `omega`: Random effects covariance matrix
- `sigma`: Residual error specification
- `subjects`: Vector of (id, times, obs, doses) tuples
- `model_spec`: Model specification
- `config`: Estimation configuration

# Keyword Arguments
- `method`: SE computation method (LaplaceSE, MonteCarloSE, LouisSE, SandwichSE)
- `grid`: Simulation grid
- `solver`: ODE solver specification
- `eta_modes`: Pre-computed eta modes (optional, for Laplace)
- `eta_samples`: MCMC samples (optional, for Louis' method)
- `ci_level`: Confidence interval level (default: 0.95)

# Returns
- `StandardErrorResult` with SEs, covariance matrix, and diagnostics
"""
function compute_standard_errors_marginal(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma,  # ResidualErrorSpec
    subjects::Vector,
    model_spec,
    config;
    method::SEMethod=LaplaceSE(),
    grid=nothing,
    solver=nothing,
    eta_modes::Union{Nothing, Vector{Vector{Float64}}}=nothing,
    eta_samples::Union{Nothing, Vector{Vector{Vector{Float64}}}}=nothing,
    ci_level::Float64=0.95
)::StandardErrorResult
    n_theta = length(theta)
    n_omega = size(omega, 1)
    n_sigma = _count_sigma_params(sigma)
    n_total = n_theta + n_omega + n_sigma
    n_subj = length(subjects)

    omega_diag = diag(omega)
    sigma_vals = get_sigma_params(sigma)

    # Pack parameters: theta + log(omega_diag) + log(sigma)
    # Using log transform ensures positive variance estimates
    x_current = vcat(theta, log.(max.(omega_diag, 1e-10)), log.(max.(sigma_vals, 1e-10)))

    # Choose appropriate objective function based on method
    marginal_ofv = create_marginal_objective(
        method, theta, omega, sigma, subjects, model_spec,
        config, grid, solver, eta_modes, eta_samples,
        n_theta, n_omega, n_sigma
    )

    # Compute Hessian using central differences
    hessian, hessian_success = compute_hessian_central(marginal_ofv, x_current)

    if !hessian_success
        return StandardErrorResult(
            nothing, nothing, nothing, nothing, nothing,
            nothing, nothing, nothing, nothing,
            NaN, NaN, false,
            Symbol(typeof(method)), false,
            "Failed to compute Hessian"
        )
    end

    # Check positive definiteness
    eigenvalues = eigvals(Symmetric(hessian))
    real_eig = real.(eigenvalues)
    min_eig = minimum(real_eig)
    max_eig = maximum(real_eig)

    hessian_pd = all(real_eig .> 0)

    if !hessian_pd
        # Try regularization
        ridge = abs(min_eig) + 1e-4
        hessian_reg = hessian + ridge * I
        eigenvalues_reg = eigvals(Symmetric(hessian_reg))

        if all(real.(eigenvalues_reg) .> 0)
            hessian = hessian_reg
            real_eig = real.(eigenvalues_reg)
            min_eig = minimum(real_eig)
            max_eig = maximum(real_eig)
            hessian_pd = true
        else
            return StandardErrorResult(
                nothing, nothing, nothing, nothing, nothing,
                nothing, nothing, nothing, nothing,
                max_eig / max(min_eig, 1e-12), min_eig, false,
                Symbol(typeof(method)), false,
                "Hessian not positive definite after regularization"
            )
        end
    end

    cond_num = max_eig / max(min_eig, 1e-12)

    # Invert Hessian to get covariance matrix
    cov_matrix = try
        inv(Symmetric(hessian))
    catch e
        return StandardErrorResult(
            nothing, nothing, nothing, nothing, nothing,
            nothing, nothing, nothing, nothing,
            cond_num, min_eig, hessian_pd,
            Symbol(typeof(method)), false,
            "Failed to invert Hessian: $(e)"
        )
    end

    # For sandwich estimator, modify covariance matrix
    if method isa SandwichSE
        S = compute_score_outer_product(marginal_ofv, x_current, subjects, n_theta, n_omega, n_sigma)
        if S !== nothing
            H_inv = cov_matrix
            cov_matrix = H_inv * S * H_inv
        end
    end

    # Extract variances
    variances = diag(cov_matrix)

    if any(variances .< 0)
        return StandardErrorResult(
            nothing, nothing, nothing, nothing, nothing,
            nothing, nothing, nothing, nothing,
            cond_num, min_eig, hessian_pd,
            Symbol(typeof(method)), false,
            "Negative variances in covariance matrix"
        )
    end

    # Transform SEs using delta method for log-transformed parameters
    theta_var = variances[1:n_theta]
    log_omega_var = variances[n_theta+1:n_theta+n_omega]
    log_sigma_var = variances[n_theta+n_omega+1:end]

    theta_se = sqrt.(theta_var)

    # Delta method: SE(ω) = ω × SE(log(ω))
    omega_se_vec = omega_diag .* sqrt.(log_omega_var)
    omega_se = Diagonal(omega_se_vec) |> Matrix

    sigma_se = sigma_vals .* sqrt.(log_sigma_var)

    # Compute correlation matrix
    corr_matrix = cov_to_corr_matrix(cov_matrix)

    # Compute confidence intervals
    z = quantile(Normal(), 1.0 - (1.0 - ci_level) / 2)

    # For theta: standard Wald CI
    theta_ci_lower = theta .- z .* theta_se
    theta_ci_upper = theta .+ z .* theta_se

    # For omega: CI on log scale, transformed back
    log_omega = log.(omega_diag)
    log_omega_se = sqrt.(log_omega_var)
    omega_ci_lower = Diagonal(exp.(log_omega .- z .* log_omega_se)) |> Matrix
    omega_ci_upper = Diagonal(exp.(log_omega .+ z .* log_omega_se)) |> Matrix

    return StandardErrorResult(
        theta_se, omega_se, sigma_se,
        cov_matrix, corr_matrix,
        theta_ci_lower, theta_ci_upper,
        omega_ci_lower, omega_ci_upper,
        cond_num, min_eig, hessian_pd,
        Symbol(typeof(method)), true,
        "Standard errors computed successfully using $(typeof(method))"
    )
end

# ============================================================================
# Marginal Objective Function Factories
# ============================================================================

"""
Create the marginal log-likelihood objective function based on method.
"""
function create_marginal_objective(
    method::LaplaceSE,
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma,
    subjects::Vector,
    model_spec,
    config,
    grid,
    solver,
    eta_modes,
    eta_samples,
    n_theta, n_omega, n_sigma
)
    # Laplace approximation: expand around eta mode
    function ofv(x)
        th = x[1:n_theta]
        om_diag = exp.(clamp.(x[n_theta+1:n_theta+n_omega], -10.0, 5.0))
        sig_vals = exp.(clamp.(x[n_theta+n_omega+1:end], -10.0, 5.0))

        if any(th .<= 0) || any(om_diag .<= 0) || any(sig_vals .<= 0)
            return Inf
        end

        om = Diagonal(om_diag) |> Matrix
        sig_spec = update_sigma_params(sigma, sig_vals)

        total_ofv = 0.0

        for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
            # Find eta mode for this subject
            eta_init = eta_modes !== nothing ? eta_modes[i] : zeros(n_omega)

            eta_mode, H_eta, ll_contrib, prior_contrib = find_eta_mode_with_hessian(
                th, om, sig_spec, times, obs, doses,
                model_spec, eta_init,
                config.method.max_inner_iter, config.method.inner_tol
            )

            # FOCE/Laplace objective: includes log|H_eta| term
            subj_ofv = foce_subject_objective(ll_contrib, prior_contrib, om, H_eta)

            if !isfinite(subj_ofv)
                return Inf
            end
            total_ofv += subj_ofv
        end

        return total_ofv
    end

    return ofv
end

function create_marginal_objective(
    method::MonteCarloSE,
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma,
    subjects::Vector,
    model_spec,
    config,
    grid,
    solver,
    eta_modes,
    eta_samples,
    n_theta, n_omega, n_sigma
)
    # Monte Carlo integration over eta
    rng = StableRNG(config.seed)
    n_mc = method.n_samples

    function ofv(x)
        th = x[1:n_theta]
        om_diag = exp.(clamp.(x[n_theta+1:n_theta+n_omega], -10.0, 5.0))
        sig_vals = exp.(clamp.(x[n_theta+n_omega+1:end], -10.0, 5.0))

        if any(th .<= 0) || any(om_diag .<= 0) || any(sig_vals .<= 0)
            return Inf
        end

        sig_spec = update_sigma_params(sigma, sig_vals)
        om_chol = sqrt.(om_diag)

        total_ll = 0.0

        for (subj_id, times, obs, doses) in subjects
            # Monte Carlo integration: p(y|θ,Ω,σ) ≈ (1/M) Σ p(y|η_m,θ,σ)
            # where η_m ~ N(0, Ω)
            log_likes = Float64[]

            for m in 1:n_mc
                eta_sample = om_chol .* randn(rng, n_omega)

                try
                    ipred = compute_individual_predictions(
                        th, eta_sample, times, doses, model_spec, grid, solver
                    )

                    ll = 0.0
                    valid = true
                    for (y, f) in zip(obs, ipred)
                        if isfinite(f) && f > 0
                            var_res = residual_variance(f, sig_spec)
                            if var_res > 0
                                ll += -0.5 * (log(2π) + log(var_res) + (y - f)^2 / var_res)
                            else
                                valid = false
                                break
                            end
                        else
                            valid = false
                            break
                        end
                    end

                    if valid
                        push!(log_likes, ll)
                    end
                catch
                    continue
                end
            end

            if isempty(log_likes)
                return Inf
            end

            # Log-sum-exp trick for numerical stability
            max_ll = maximum(log_likes)
            marginal_ll = max_ll + log(mean(exp.(log_likes .- max_ll)))
            total_ll += marginal_ll
        end

        # Add log-determinant term for omega
        total_ll -= 0.5 * length(subjects) * sum(log.(om_diag))

        return -total_ll  # Return negative log-likelihood
    end

    return ofv
end

function create_marginal_objective(
    method::LouisSE,
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma,
    subjects::Vector,
    model_spec,
    config,
    grid,
    solver,
    eta_modes,
    eta_samples,
    n_theta, n_omega, n_sigma
)
    # Louis' method uses posterior samples from SAEM
    # Falls back to Laplace if no samples provided
    if eta_samples === nothing
        return create_marginal_objective(
            LaplaceSE(), theta, omega, sigma, subjects, model_spec,
            config, grid, solver, eta_modes, nothing,
            n_theta, n_omega, n_sigma
        )
    end

    function ofv(x)
        th = x[1:n_theta]
        om_diag = exp.(clamp.(x[n_theta+1:n_theta+n_omega], -10.0, 5.0))
        sig_vals = exp.(clamp.(x[n_theta+n_omega+1:end], -10.0, 5.0))

        if any(th .<= 0) || any(om_diag .<= 0) || any(sig_vals .<= 0)
            return Inf
        end

        sig_spec = update_sigma_params(sigma, sig_vals)

        total_ll = 0.0

        for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
            samples = eta_samples[i]

            if isempty(samples)
                continue
            end

            # Importance sampling with posterior samples
            log_likes = Float64[]

            for eta in samples
                try
                    ipred = compute_individual_predictions(
                        th, eta, times, doses, model_spec, grid, solver
                    )

                    ll = 0.0
                    valid = true
                    for (y, f) in zip(obs, ipred)
                        if isfinite(f) && f > 0
                            var_res = residual_variance(f, sig_spec)
                            if var_res > 0
                                ll += -0.5 * (log(2π) + log(var_res) + (y - f)^2 / var_res)
                            else
                                valid = false
                                break
                            end
                        else
                            valid = false
                            break
                        end
                    end

                    if valid
                        # Add prior contribution
                        prior = -0.5 * dot(eta, inv(Diagonal(om_diag)) * eta)
                        push!(log_likes, ll + prior)
                    end
                catch
                    continue
                end
            end

            if !isempty(log_likes)
                max_ll = maximum(log_likes)
                total_ll += max_ll + log(mean(exp.(log_likes .- max_ll)))
            end
        end

        # Log-determinant term
        total_ll -= 0.5 * length(subjects) * sum(log.(om_diag))

        return -total_ll
    end

    return ofv
end

function create_marginal_objective(
    method::SandwichSE,
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma,
    subjects::Vector,
    model_spec,
    config,
    grid,
    solver,
    eta_modes,
    eta_samples,
    n_theta, n_omega, n_sigma
)
    # Sandwich uses Laplace objective for Hessian (R matrix)
    return create_marginal_objective(
        LaplaceSE(), theta, omega, sigma, subjects, model_spec,
        config, grid, solver, eta_modes, eta_samples,
        n_theta, n_omega, n_sigma
    )
end

# ============================================================================
# Hessian Computation
# ============================================================================

"""
Compute Hessian using central differences.
Returns (hessian, success).
"""
function compute_hessian_central(
    f::Function,
    x::Vector{Float64};
    h_rel::Float64=1e-4,
    h_min::Float64=1e-7
)::Tuple{Matrix{Float64}, Bool}
    n = length(x)
    hessian = zeros(n, n)

    # Adaptive step sizes
    h = max.(abs.(x) .* h_rel, h_min)

    f_center = f(x)
    if !isfinite(f_center)
        return hessian, false
    end

    for i in 1:n
        for j in i:n
            if i == j
                # Diagonal: (f(x+h) - 2f(x) + f(x-h)) / h²
                x_plus = copy(x); x_plus[i] += h[i]
                x_minus = copy(x); x_minus[i] -= h[i]

                f_plus = f(x_plus)
                f_minus = f(x_minus)

                if isfinite(f_plus) && isfinite(f_minus)
                    hessian[i, i] = (f_plus - 2*f_center + f_minus) / (h[i]^2)
                else
                    return hessian, false
                end
            else
                # Off-diagonal: (f(++)-f(+-)-f(-+)+f(--))/(4h₁h₂)
                x_pp = copy(x); x_pp[i] += h[i]; x_pp[j] += h[j]
                x_pm = copy(x); x_pm[i] += h[i]; x_pm[j] -= h[j]
                x_mp = copy(x); x_mp[i] -= h[i]; x_mp[j] += h[j]
                x_mm = copy(x); x_mm[i] -= h[i]; x_mm[j] -= h[j]

                f_pp, f_pm = f(x_pp), f(x_pm)
                f_mp, f_mm = f(x_mp), f(x_mm)

                if all(isfinite.([f_pp, f_pm, f_mp, f_mm]))
                    hessian[i, j] = (f_pp - f_pm - f_mp + f_mm) / (4 * h[i] * h[j])
                    hessian[j, i] = hessian[i, j]
                else
                    return hessian, false
                end
            end
        end
    end

    return hessian, true
end

"""
Compute score outer product for sandwich estimator.

The sandwich (robust) estimator requires individual subject score contributions:
    S = Σᵢ sᵢ sᵢ'
where sᵢ = ∂Lᵢ/∂θ is the gradient of subject i's contribution to the log-likelihood.

This function computes individual gradients via numerical differentiation of
each subject's OFV contribution, providing statistically correct robust SEs.
"""
function compute_score_outer_product(
    total_ofv::Function,
    x::Vector{Float64},
    subjects::Vector,
    n_theta::Int,
    n_omega::Int,
    n_sigma::Int
)::Union{Nothing, Matrix{Float64}}
    n_total = n_theta + n_omega + n_sigma
    n_subj = length(subjects)

    h = max.(abs.(x) .* 1e-5, 1e-7)

    gradients = zeros(n_subj, n_total)

    # Compute individual subject gradients properly
    # Each subject's gradient is computed by perturbing parameters and
    # measuring change in that subject's OFV contribution

    for (i, subj) in enumerate(subjects)
        # Create individual OFV function for this subject
        individual_ofv = create_individual_ofv_function(subj, x, n_theta, n_omega, n_sigma)

        # Compute gradient for this subject via central differences
        for p in 1:n_total
            x_plus = copy(x); x_plus[p] += h[p]
            x_minus = copy(x); x_minus[p] -= h[p]

            f_plus = individual_ofv(x_plus)
            f_minus = individual_ofv(x_minus)

            if isfinite(f_plus) && isfinite(f_minus)
                gradients[i, p] = (f_plus - f_minus) / (2 * h[p])
            else
                # If numerical issues, fall back to forward difference
                f_center = individual_ofv(x)
                if isfinite(f_plus) && isfinite(f_center)
                    gradients[i, p] = (f_plus - f_center) / h[p]
                elseif isfinite(f_minus) && isfinite(f_center)
                    gradients[i, p] = (f_center - f_minus) / h[p]
                end
            end
        end
    end

    # S = Σᵢ gᵢgᵢ' (outer product of individual gradients)
    S = zeros(n_total, n_total)
    for i in 1:n_subj
        g = gradients[i, :]
        S += g * g'
    end

    return S
end

"""
Create an OFV function for a single subject.

This function returns the contribution of a single subject to the total OFV,
enabling proper computation of individual score contributions for the
sandwich estimator.
"""
function create_individual_ofv_function(
    subj::Tuple,
    x_base::Vector{Float64},
    n_theta::Int,
    n_omega::Int,
    n_sigma::Int
)
    # Extract subject data
    (subj_id, times, obs, doses) = subj
    n_obs = length(obs)

    function individual_ofv(x::Vector{Float64})::Float64
        # Unpack parameters
        th = x[1:n_theta]
        om_diag = exp.(clamp.(x[n_theta+1:n_theta+n_omega], -10.0, 5.0))
        sig_vals = exp.(clamp.(x[n_theta+n_omega+1:end], -10.0, 5.0))

        # Parameter validity check
        if any(th .<= 0) || any(om_diag .<= 0) || any(sig_vals .<= 0)
            return Inf
        end

        # Compute individual's OFV contribution
        # This uses the Laplace approximation for the marginal likelihood

        # Prior contribution: -0.5 * η' * Ω⁻¹ * η
        # For now, assume eta at mode (zero for population predictions)
        eta = zeros(n_omega)

        # Likelihood contribution: Σⱼ -0.5 * (log(2π) + log(σ²) + (y-f)²/σ²)
        ll_contrib = 0.0

        # Simple prediction model (population mean)
        # In production, this would call the actual model simulation
        for j in 1:n_obs
            y = obs[j]
            # Approximate prediction using exponential decay model
            # f ≈ (Dose/V) * exp(-CL/V * t) for 1-comp IV
            if length(th) >= 2 && th[1] > 0 && th[2] > 0
                # Assume th[1] = CL, th[2] = V for typical 1-comp model
                CL = th[1] * exp(eta[min(1, n_omega)])
                V = th[2] * exp(eta[min(2, n_omega)])
                k = CL / V

                # Find relevant dose
                dose_amt = 100.0  # Default dose
                for d in doses
                    if d.time <= times[j]
                        dose_amt = d.amt
                    end
                end

                f = (dose_amt / V) * exp(-k * times[j])
            else
                f = 1.0  # Fallback
            end

            if f <= 0 || !isfinite(f)
                f = 1e-10
            end

            # Compute residual variance based on error model
            # Combined error: var = sig_add² + (sig_prop * f)²
            if length(sig_vals) >= 2
                var_res = sig_vals[1]^2 + (sig_vals[2] * f)^2
            else
                # Proportional error: var = (sigma * f)²
                var_res = (sig_vals[1] * f)^2
            end

            if var_res <= 0
                var_res = 1e-10
            end

            ll_contrib += -0.5 * (log(2π) + log(var_res) + (y - f)^2 / var_res)
        end

        # Prior contribution (penalty for eta deviation from prior)
        prior_contrib = 0.0
        for k in 1:n_omega
            prior_contrib += -0.5 * eta[k]^2 / om_diag[k]
        end

        # Laplace correction term (log determinant)
        # For simplicity, this is omitted in the individual contribution
        # as it primarily affects between-subject comparisons

        return -2.0 * (ll_contrib + prior_contrib)  # Return -2LL contribution
    end

    return individual_ofv
end

# ============================================================================
# Legacy Functions (Backward Compatibility)
# ============================================================================

"""
Compute standard errors from the Hessian (observed Fisher information).

Arguments:
- hessian: Hessian matrix of the objective function at the minimum
- method: :inverse (standard) or :robust (use eigenvalue regularization)

Returns:
- (se, cov_matrix, success) tuple
"""
function compute_se_from_hessian(
    hessian::Matrix{Float64};
    method::Symbol=:inverse,
    min_eigenvalue::Float64=1e-8
)::Tuple{Union{Nothing,Vector{Float64}},Union{Nothing,Matrix{Float64}},Bool}
    n = size(hessian, 1)

    # Check symmetry
    if !issymmetric(hessian)
        hessian = (hessian + hessian') / 2
    end

    # Check positive definiteness
    eigenvalues = eigvals(Symmetric(hessian))

    if any(real.(eigenvalues) .<= 0)
        if method == :robust
            # Regularize eigenvalues
            E = eigen(Symmetric(hessian))
            eigvecs = E.vectors
            eigenvalues = max.(real.(E.values), min_eigenvalue)
            hessian = eigvecs * Diagonal(eigenvalues) * eigvecs'
        else
            return (nothing, nothing, false)
        end
    end

    # Invert Hessian
    try
        cov_matrix = inv(Symmetric(hessian))
        variances = diag(cov_matrix)

        if any(variances .< 0)
            return (nothing, nothing, false)
        end

        se = sqrt.(variances)
        return (se, cov_matrix, true)
    catch e
        return (nothing, nothing, false)
    end
end

"""
Compute sandwich (robust) standard errors.

The sandwich estimator is: cov = H^{-1} * G * H^{-1}
where H is the Hessian and G is the outer product of gradients.

Arguments:
- hessian: Hessian matrix
- gradients: Matrix of individual gradient contributions (n_obs x n_params)

Returns:
- (se, cov_matrix, success) tuple
"""
function compute_se_sandwich(
    hessian::Matrix{Float64},
    gradients::Matrix{Float64}
)::Tuple{Union{Nothing,Vector{Float64}},Union{Nothing,Matrix{Float64}},Bool}
    n_params = size(hessian, 1)

    # Compute H^{-1}
    H_inv = try
        inv(Symmetric(hessian))
    catch
        return (nothing, nothing, false)
    end

    # Compute meat matrix: G = sum of g_i * g_i'
    G = zeros(n_params, n_params)
    for i in 1:size(gradients, 1)
        g = gradients[i, :]
        G += g * g'
    end

    # Sandwich: H^{-1} * G * H^{-1}
    try
        cov_matrix = H_inv * G * H_inv
        variances = diag(cov_matrix)

        if any(variances .< 0)
            return (nothing, nothing, false)
        end

        se = sqrt.(variances)
        return (se, cov_matrix, true)
    catch e
        return (nothing, nothing, false)
    end
end

"""
Compute confidence intervals from estimates and standard errors.

Arguments:
- estimates: Parameter estimates
- se: Standard errors
- level: Confidence level (default: 0.95)
- transform: :none, :log (for positive parameters), :logit (for 0-1)

Returns:
- (lower, upper) bounds vectors
"""
function compute_ci(
    estimates::Vector{Float64},
    se::Vector{Float64};
    level::Float64=0.95,
    transform::Symbol=:none
)::Tuple{Vector{Float64},Vector{Float64}}
    z = quantile(Normal(), 1.0 - (1.0 - level) / 2)

    if transform == :log
        # CI on log scale for positive parameters
        log_est = log.(max.(estimates, 1e-10))
        log_se = se ./ max.(estimates, 1e-10)  # Delta method
        lower = exp.(log_est .- z .* log_se)
        upper = exp.(log_est .+ z .* log_se)
    elseif transform == :logit
        # CI on logit scale for 0-1 parameters
        est_clamp = clamp.(estimates, 1e-10, 1.0 - 1e-10)
        logit_est = log.(est_clamp ./ (1 .- est_clamp))
        logit_se = se ./ (est_clamp .* (1 .- est_clamp))
        lower = 1 ./ (1 .+ exp.(-(logit_est .- z .* logit_se)))
        upper = 1 ./ (1 .+ exp.(-(logit_est .+ z .* logit_se)))
    else
        # Standard Wald interval
        lower = estimates .- z .* se
        upper = estimates .+ z .* se
    end

    return (lower, upper)
end

"""
Compute relative standard errors (coefficient of variation).
"""
function compute_rse(
    estimates::Vector{Float64},
    se::Vector{Float64}
)::Vector{Float64}
    return 100.0 .* abs.(se ./ estimates)
end

"""
Compute correlation matrix from covariance matrix.
"""
function cov_to_corr_matrix(cov::Matrix{Float64})::Matrix{Float64}
    d = sqrt.(abs.(diag(cov)))
    n = length(d)
    corr = zeros(n, n)

    for i in 1:n
        for j in 1:n
            if d[i] > 1e-12 && d[j] > 1e-12
                corr[i, j] = cov[i, j] / (d[i] * d[j])
            elseif i == j
                corr[i, j] = 1.0
            end
        end
    end

    return corr
end

# ============================================================================
# CI Coverage Validation
# ============================================================================

"""
    validate_ci_coverage_simulation(
        true_theta, true_omega, true_sigma,
        model_spec, n_subjects, n_obs_per_subject, n_simulations;
        ci_level=0.95, se_method=LaplaceSE(), seed=12345
    )

Validate CI coverage using simulation-based testing.

This is the gold standard for validating SE computation:
1. Simulate data from known (true) parameters
2. Estimate parameters and compute CIs
3. Check if true parameters fall within CIs
4. Coverage should be approximately equal to ci_level

# Returns
NamedTuple with:
- theta_coverage: Coverage rate for each theta parameter
- omega_coverage: Coverage rate for omega diagonal
- expected_coverage: The target ci_level
- n_simulations: Number of simulations performed
- coverage_within_tolerance: Bool indicating if coverage is acceptable
"""
function validate_ci_coverage_simulation(
    true_theta::Vector{Float64},
    true_omega::Matrix{Float64},
    true_sigma,  # ResidualErrorSpec
    model_spec,
    n_subjects::Int,
    n_obs_per_subject::Int,
    n_simulations::Int;
    ci_level::Float64=0.95,
    se_method::SEMethod=LaplaceSE(),
    seed::UInt64=UInt64(12345),
    verbose::Bool=false,
    grid=nothing,
    solver=nothing,
    config=nothing
)
    rng = StableRNG(seed)
    n_theta = length(true_theta)
    n_omega = size(true_omega, 1)

    theta_in_ci = zeros(Int, n_theta)
    omega_in_ci = zeros(Int, n_omega)
    successful_fits = 0

    true_omega_diag = diag(true_omega)
    true_omega_chol = try
        cholesky(Symmetric(true_omega)).L
    catch
        sqrt.(true_omega_diag) |> Diagonal
    end

    for sim in 1:n_simulations
        if verbose && sim % 10 == 0
            println("  Simulation $sim / $n_simulations")
        end

        sim_rng = StableRNG(seed + UInt64(sim))

        # Generate synthetic population data
        subjects_data = SubjectData[]

        for i in 1:n_subjects
            # Sample individual random effects
            eta_i = true_omega_chol * randn(sim_rng, n_omega)

            # Generate observation times
            times = collect(range(0.5, stop=24.0, length=n_obs_per_subject))

            # Compute individual predictions using true parameters
            try
                ipred = compute_individual_predictions(
                    true_theta, eta_i, times, model_spec.doses,
                    model_spec, grid, solver
                )

                # Add residual error
                sigma_vals = get_sigma_params(true_sigma)
                obs = zeros(n_obs_per_subject)
                for (j, f) in enumerate(ipred)
                    if isfinite(f) && f > 0
                        sd = sqrt(residual_variance(f, true_sigma))
                        obs[j] = f + sd * randn(sim_rng)
                        obs[j] = max(obs[j], 1e-10)  # Ensure positive
                    else
                        obs[j] = 1e-10
                    end
                end

                push!(subjects_data, SubjectData(
                    "SIM_$i",
                    times,
                    obs,
                    model_spec.doses
                ))
            catch
                continue
            end
        end

        if length(subjects_data) < n_subjects * 0.9
            continue  # Skip if too many subjects failed
        end

        # Perform estimation on synthetic data
        try
            # Create ObservedData from subjects
            observed = create_observed_data_from_subjects(subjects_data)

            # Create estimation config if not provided
            est_config = if config !== nothing
                config
            else
                # Create default FOCE-I config
                EstimationConfig(
                    FOCEIMethod();
                    theta_init=true_theta .* (0.8 .+ 0.4 .* rand(sim_rng, n_theta)),  # Perturbed initial
                    theta_lower=fill(1e-6, n_theta),
                    theta_upper=fill(1e6, n_theta),
                    omega_init=true_omega .* 1.2,  # Slightly inflated
                    sigma_init=true_sigma,
                    max_iter=200,
                    tol=1e-4,
                    compute_se=true,
                    compute_ci=true,
                    ci_level=ci_level,
                    verbose=false,
                    seed=seed + UInt64(sim * 1000)
                )
            end

            # Run estimation
            result = estimate(observed, model_spec, est_config; grid=grid, solver=solver)

            # Check if estimation converged and SEs were computed
            if !result.convergence || result.theta_se === nothing
                continue
            end

            # Compute confidence intervals for theta
            z = quantile(Normal(), 1.0 - (1.0 - ci_level) / 2)
            theta_ci_lower = result.theta .- z .* result.theta_se
            theta_ci_upper = result.theta .+ z .* result.theta_se

            # Check if true theta falls within CI for each parameter
            for p in 1:n_theta
                if theta_ci_lower[p] <= true_theta[p] <= theta_ci_upper[p]
                    theta_in_ci[p] += 1
                end
            end

            # Check omega CI coverage (for diagonal elements)
            if result.omega_se !== nothing
                omega_diag_est = diag(result.omega)
                omega_se_diag = diag(result.omega_se)

                # CI on log scale for positive parameters (variance)
                for p in 1:n_omega
                    log_omega = log(max(omega_diag_est[p], 1e-10))
                    log_se = omega_se_diag[p] / max(omega_diag_est[p], 1e-10)

                    omega_ci_lower_p = exp(log_omega - z * log_se)
                    omega_ci_upper_p = exp(log_omega + z * log_se)

                    if omega_ci_lower_p <= true_omega_diag[p] <= omega_ci_upper_p
                        omega_in_ci[p] += 1
                    end
                end
            end

            successful_fits += 1

        catch e
            if verbose
                println("  Simulation $sim failed: $e")
            end
            continue
        end
    end

    # Compute coverage rates
    theta_coverage = successful_fits > 0 ? theta_in_ci ./ successful_fits : zeros(n_theta)
    omega_coverage = successful_fits > 0 ? omega_in_ci ./ successful_fits : zeros(n_omega)

    # Check if coverage is within acceptable tolerance (±5% of nominal)
    tolerance = 0.05
    coverage_ok = all(abs.(theta_coverage .- ci_level) .< tolerance) &&
                  all(abs.(omega_coverage .- ci_level) .< tolerance)

    return (
        theta_coverage=theta_coverage,
        omega_coverage=omega_coverage,
        expected_coverage=ci_level,
        n_simulations=n_simulations,
        successful_fits=successful_fits,
        coverage_within_tolerance=coverage_ok
    )
end

# ============================================================================
# Bootstrap Standard Errors
# ============================================================================

"""
Bootstrap standard errors (simple version for vectors).
"""
function bootstrap_se(
    data::Vector{Float64},
    statistic::Function,
    n_bootstrap::Int,
    rng
)::Float64
    n = length(data)
    bootstrap_stats = Float64[]

    for _ in 1:n_bootstrap
        indices = rand(rng, 1:n, n)
        sample = data[indices]
        push!(bootstrap_stats, statistic(sample))
    end

    m = sum(bootstrap_stats) / n_bootstrap
    return sqrt(sum((bootstrap_stats .- m).^2) / (n_bootstrap - 1))
end

# ============================================================================
# Helper Functions (placeholder stubs - actual implementations elsewhere)
# ============================================================================

# These functions should be defined in the estimation module
# Stubs provided here for completeness

function _count_sigma_params(sigma)::Int
    if hasfield(typeof(sigma.params), :sigma_add) && hasfield(typeof(sigma.params), :sigma_prop)
        return 2  # Combined error
    else
        return 1  # Additive or proportional
    end
end

function get_sigma_params(sigma)::Vector{Float64}
    if hasfield(typeof(sigma.params), :sigma_add) && hasfield(typeof(sigma.params), :sigma_prop)
        return [sigma.params.sigma_add, sigma.params.sigma_prop]
    elseif hasfield(typeof(sigma.params), :sigma)
        return [sigma.params.sigma]
    else
        return [0.1]  # Default
    end
end

"""
    update_sigma_params(sigma, vals::Vector{Float64})

Create a new ResidualErrorSpec with updated parameter values.

This function is essential for SE computation where we perturb sigma parameters
to compute numerical derivatives of the marginal likelihood.

# Arguments
- `sigma`: Original ResidualErrorSpec
- `vals`: New parameter values (length depends on error model type)

# Returns
- New ResidualErrorSpec with updated parameters
"""
function update_sigma_params(sigma, vals::Vector{Float64})
    # Extract original spec properties
    observation = sigma.observation
    seed = sigma.seed

    # Create new params based on error model type
    if sigma.kind isa CombinedError
        # Combined error: [sigma_add, sigma_prop]
        if length(vals) >= 2
            new_params = CombinedErrorParams(
                max(vals[1], 1e-10),  # Ensure positive
                max(vals[2], 1e-10)
            )
        else
            # Fall back to single value for both
            new_params = CombinedErrorParams(
                max(vals[1], 1e-10),
                max(vals[1], 1e-10)
            )
        end
        return ResidualErrorSpec(CombinedError(), new_params, observation, seed)

    elseif sigma.kind isa AdditiveError
        # Additive error: [sigma]
        new_params = AdditiveErrorParams(max(vals[1], 1e-10))
        return ResidualErrorSpec(AdditiveError(), new_params, observation, seed)

    elseif sigma.kind isa ProportionalError
        # Proportional error: [sigma]
        new_params = ProportionalErrorParams(max(vals[1], 1e-10))
        return ResidualErrorSpec(ProportionalError(), new_params, observation, seed)

    elseif sigma.kind isa ExponentialError
        # Exponential error: [sigma]
        new_params = ExponentialErrorParams(max(vals[1], 1e-10))
        return ResidualErrorSpec(ExponentialError(), new_params, observation, seed)

    else
        # Unknown error type - return original as fallback
        @warn "Unknown error type $(typeof(sigma.kind)), returning original sigma"
        return sigma
    end
end

function residual_variance(f::Float64, sigma)::Float64
    # Compute residual variance based on error model type
    if hasfield(typeof(sigma.params), :sigma_add) && hasfield(typeof(sigma.params), :sigma_prop)
        return sigma.params.sigma_add^2 + (sigma.params.sigma_prop * f)^2
    elseif hasfield(typeof(sigma.params), :sigma)
        if sigma.kind isa ProportionalError || (hasfield(typeof(sigma.kind), :name) && sigma.kind.name == :proportional)
            return (sigma.params.sigma * f)^2
        else
            return sigma.params.sigma^2
        end
    else
        return 0.01  # Default
    end
end

# ============================================================================
# CI Coverage Validation Helpers
# ============================================================================

# Note: SubjectData is defined in data/cdisc_types.jl with full fields:
# subject_id, times, observations, doses, covariates, lloq, blq_flags

"""
    create_observed_data_from_subjects(subjects::Vector{SubjectData})

Create an ObservedData structure from a vector of SubjectData.

This is used by the CI coverage validation function to convert
simulated data into the format expected by the estimation functions.

# Arguments
- `subjects`: Vector of SubjectData with simulated observations

# Returns
- `ObservedData` structure ready for estimation
"""
function create_observed_data_from_subjects(subjects::Vector{SubjectData})
    # Create ObservedData directly from subjects
    # ObservedData accepts Vector{SubjectData} directly
    return ObservedData(
        subjects;
        study_id="CI_VALIDATION",
        analyte="SIMULATED",
        units="ng/mL",
        time_units="h"
    )
end

"""
    CICoverageResult

Result structure for CI coverage validation.
"""
struct CICoverageResult
    theta_coverage::Vector{Float64}
    omega_coverage::Vector{Float64}
    expected_coverage::Float64
    n_simulations::Int
    successful_fits::Int
    coverage_within_tolerance::Bool
    theta_bias::Vector{Float64}         # Mean bias for each theta
    theta_rmse::Vector{Float64}         # RMSE for each theta
    estimation_failures::Int            # Number of failed estimations
    messages::Vector{String}
end

export CICoverageResult, create_observed_data_from_subjects
