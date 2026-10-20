# Laplacian Estimation Method
# Simple approximation using mode of posterior for random effects

using Optim
using LinearAlgebra

export laplacian_estimate

"""
    laplacian_estimate(observed, model_spec, config, grid, solver, rng) -> EstimationResult

Estimate parameters using the Laplacian approximation method.

The Laplacian method approximates the marginal likelihood by:
1. Finding the mode of the posterior for eta (random effects)
2. Using the Laplace approximation for the integral over eta

This is computationally efficient but less accurate than FOCE-I for
complex models with high IIV.
"""
function laplacian_estimate(
    observed::ObservedData,
    model_spec::ModelSpec,
    config::EstimationConfig{LaplacianMethod},
    grid::SimGrid,
    solver::SolverSpec,
    rng;
    parallel_config::ParallelConfig=ParallelConfig(SerialBackend())
)::EstimationResult
    start_time = time()
    n_theta = length(config.theta_init)
    n_eta = size(config.omega_init, 1)
    n_subj = n_subjects(observed)
    n_obs = n_observations(observed)

    subjects = extract_subject_data(observed)

    # Extract BLQ information for each subject
    blq_config = config.blq_config
    subject_blq_flags = Vector{Vector{Bool}}(undef, n_subj)
    subject_lloq = Vector{Float64}(undef, n_subj)
    for (i, s) in enumerate(observed.subjects)
        subject_blq_flags[i] = get_blq_flags_for_subject(s, blq_config)
        subject_lloq[i] = get_lloq_for_subject(s, blq_config)
    end

    # Pack parameters for optimization: [theta; omega_chol_diag; sigma_params]
    # For simplicity, we'll optimize theta and omega diagonal, fix sigma initially

    # Initialize
    theta_current = copy(config.theta_init)
    omega_current = copy(config.omega_init)
    sigma_current = config.sigma_init

    # Storage for individual etas
    eta_per_subject = [zeros(n_eta) for _ in 1:n_subj]

    # Optimization history
    ofv_history = Float64[]
    converged = false
    n_iter = 0
    final_gradient_norm = Inf

    # Main optimization loop
    for iter in 1:config.max_iter
        n_iter = iter

        # E-step: Find mode of posterior for each subject's eta (with BLQ support)
        # Uses parallel processing if enabled
        if is_parallel(parallel_config)
            subject_indices = collect(1:n_subj)
            initial_etas = [copy(eta_per_subject[i]) for i in 1:n_subj]

            eta_results = parallel_map(
                idx -> begin
                    subj_id, times, obs, doses = subjects[idx]
                    find_eta_mode_laplacian(
                        theta_current, omega_current, sigma_current,
                        times, obs, doses, model_spec, grid, solver,
                        initial_etas[idx], config.method.max_inner_iter, config.method.inner_tol;
                        blq_config=blq_config, blq_flags=subject_blq_flags[idx], lloq=subject_lloq[idx]
                    )
                end,
                subject_indices,
                parallel_config
            )

            for (i, eta) in enumerate(eta_results)
                eta_per_subject[i] = eta
            end
        else
            for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
                eta_per_subject[i] = find_eta_mode_laplacian(
                    theta_current, omega_current, sigma_current,
                    times, obs, doses, model_spec, grid, solver,
                    eta_per_subject[i], config.method.max_inner_iter, config.method.inner_tol;
                    blq_config=blq_config, blq_flags=subject_blq_flags[i], lloq=subject_lloq[i]
                )
            end
        end

        # M-step: Update theta and omega given eta modes (with BLQ support)
        theta_new, omega_new, ofv = m_step_laplacian(
            theta_current, omega_current, sigma_current,
            eta_per_subject, subjects, model_spec, grid, solver, config;
            blq_config=blq_config, subject_blq_flags=subject_blq_flags, subject_lloq=subject_lloq
        )

        push!(ofv_history, ofv)

        # Check convergence
        theta_change = norm(theta_new - theta_current) / max(1.0, norm(theta_current))
        omega_change = norm(omega_new - omega_current) / max(1.0, norm(omega_current))

        if config.verbose && iter % 10 == 0
            println("  Iter $iter: OFV=$(round(ofv, digits=3)), theta_change=$(round(theta_change, sigdigits=3))")
        end

        if theta_change < config.tol && omega_change < config.tol
            converged = true
            theta_current = theta_new
            omega_current = omega_new
            break
        end

        theta_current = theta_new
        omega_current = omega_new
    end

    # Final OFV computation with BLQ support
    final_ofv = compute_laplacian_ofv(
        theta_current, omega_current, sigma_current,
        eta_per_subject, subjects, model_spec, grid, solver;
        blq_config=blq_config, subject_blq_flags=subject_blq_flags, subject_lloq=subject_lloq
    )

    # Compute individual estimates with BLQ support
    individual_estimates = compute_individual_estimates_laplacian(
        theta_current, omega_current, sigma_current,
        eta_per_subject, subjects, model_spec, grid, solver;
        blq_config=blq_config, subject_blq_flags=subject_blq_flags, subject_lloq=subject_lloq
    )

    # Compute standard errors if requested
    theta_se = nothing
    theta_rse = nothing
    omega_se = nothing
    sigma_se = nothing
    covariance_successful = false
    condition_num = NaN
    eigenvalue_ratio = NaN

    if config.compute_se
        theta_se, omega_se, sigma_se, covariance_successful, condition_num, eigenvalue_ratio =
            compute_standard_errors_laplacian(
                theta_current, omega_current, sigma_current,
                eta_per_subject, subjects, model_spec, grid, solver, config;
                blq_config=blq_config, subject_blq_flags=subject_blq_flags, subject_lloq=subject_lloq
            )
        if theta_se !== nothing
            theta_rse = 100.0 .* abs.(theta_se ./ theta_current)
        end
    end

    # Confidence intervals
    theta_ci_lower = nothing
    theta_ci_upper = nothing
    if config.compute_ci && theta_se !== nothing
        z = quantile(Normal(), 1.0 - (1.0 - config.ci_level) / 2)
        theta_ci_lower = theta_current .- z .* theta_se
        theta_ci_upper = theta_current .+ z .* theta_se
    end

    # Compute correlation matrix from omega
    omega_corr = cov_to_corr(omega_current)

    # Model fit statistics
    n_params = n_theta + _count_omega_params(config.omega_structure, n_eta) +
               _count_sigma_params(sigma_current)
    aic = compute_aic(final_ofv, n_params)
    bic = compute_bic(final_ofv, n_params, n_obs)

    elapsed = time() - start_time

    # Build messages
    messages = String["Laplacian estimation"]

    # Compute BLQ summary if BLQ handling is enabled
    blq_summary = nothing
    if blq_config !== nothing && blq_config.report_blq_summary
        blq_summary = compute_blq_summary(observed, blq_config)
        push!(messages, "BLQ handling: $(blq_config.method) with $(blq_summary.blq_observations)/$(blq_summary.total_observations) BLQ observations ($(round(blq_summary.blq_percentage, digits=1))%)")
        for warning in blq_summary.warnings
            push!(messages, "BLQ Warning: $warning")
        end
    end

    return EstimationResult(
        config,
        theta_current, theta_se, nothing,  # theta_se_robust = nothing for Laplacian
        theta_rse, theta_ci_lower, theta_ci_upper,
        omega_current, omega_se, omega_corr,
        sigma_current, sigma_se,
        final_ofv, aic, bic,
        individual_estimates,
        converged, n_iter, final_gradient_norm, condition_num, eigenvalue_ratio,
        covariance_successful,
        messages,
        elapsed,
        blq_summary
    )
end

"""
Find the mode of the posterior for eta given current parameters.
Supports BLQ/censoring handling via optional parameters.
"""
function find_eta_mode_laplacian(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    times::Vector{Float64},
    obs::Vector{Float64},
    doses::Vector{DoseEvent},
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    eta_init::Vector{Float64},
    max_iter::Int,
    tol::Float64;
    blq_config::Union{Nothing,BLQConfig}=nothing,
    blq_flags::Vector{Bool}=Bool[],
    lloq::Float64=0.0
)::Vector{Float64}
    n_eta = length(eta_init)
    omega_inv = inv(omega)

    # Check if BLQ handling is enabled
    has_blq = blq_config !== nothing && !isempty(blq_flags)

    # Objective function: -2 * log(p(eta|y,theta))
    # = -2LL(y|theta,eta) + eta' * omega_inv * eta + log(det(omega)) + const
    function eta_objective(eta)
        # Prior contribution
        prior_contrib = eta' * omega_inv * eta

        # Likelihood contribution - wrap in try-catch for invalid parameters
        try
            ipred = compute_individual_predictions(
                theta, eta, times, doses, model_spec, grid, solver
            )

            ll_contrib = 0.0
            for (i, (y, f)) in enumerate(zip(obs, ipred))
                # Check if this observation is BLQ
                is_blq = has_blq && i <= length(blq_flags) && blq_flags[i]

                if is_blq && blq_config !== nothing
                    ll_contrib += observation_log_likelihood_blq(y, f, sigma, blq_config, true, lloq)
                else
                    ll_contrib += observation_log_likelihood(y, f, sigma)
                end
            end

            return ll_contrib + prior_contrib
        catch e
            # Invalid parameters (e.g., negative CL or V) - return large value
            return 1e10 + prior_contrib
        end
    end

    # Optimize using fallback optimizer chain (BFGS -> L-BFGS -> Nelder-Mead)
    opt_config = OptimizerConfig(
        max_attempts_per_optimizer=1,  # Keep inner optimization fast
        verbose=false
    )
    opt_options = Optim.Options(iterations=max_iter, g_tol=tol, show_trace=false)

    opt_result = optimize_with_fallback(
        eta_objective,
        eta_init,
        opt_config;
        options=opt_options
    )

    return opt_result.minimizer
end

"""
M-step: Update theta and omega given eta modes.
Supports BLQ/censoring handling via optional parameters.
"""
function m_step_laplacian(
    theta_current::Vector{Float64},
    omega_current::Matrix{Float64},
    sigma::ResidualErrorSpec,
    etas::Vector{Vector{Float64}},
    subjects::Vector,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    config::EstimationConfig;
    blq_config::Union{Nothing,BLQConfig}=nothing,
    subject_blq_flags::Vector{Vector{Bool}}=Vector{Vector{Bool}}(),
    subject_lloq::Vector{Float64}=Float64[]
)::Tuple{Vector{Float64},Matrix{Float64},Float64}
    n_subj = length(subjects)
    n_eta = size(omega_current, 1)

    # Update omega: empirical covariance of eta modes
    eta_matrix = hcat(etas...)  # n_eta x n_subj
    omega_new = (eta_matrix * eta_matrix') / n_subj

    # Ensure positive definiteness
    omega_new = ensure_positive_definite(omega_new)

    # If diagonal structure, zero out off-diagonals
    if config.omega_structure isa DiagonalOmega
        omega_new = Diagonal(diag(omega_new)) |> Matrix
    end

    # Check if BLQ handling is enabled
    has_blq = blq_config !== nothing && !isempty(subject_blq_flags)

    # Update theta: minimize objective given fixed etas
    function theta_objective(theta)
        # Check bounds
        if any(theta .< config.theta_lower) || any(theta .> config.theta_upper)
            return Inf
        end

        ofv = 0.0
        omega_inv = inv(omega_new)

        for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
            eta = etas[i]

            # Prior contribution
            ofv += eta' * omega_inv * eta

            # Get BLQ info for this subject
            blq_flags = has_blq && i <= length(subject_blq_flags) ? subject_blq_flags[i] : Bool[]
            lloq = has_blq && i <= length(subject_lloq) ? subject_lloq[i] : 0.0
            subj_has_blq = has_blq && !isempty(blq_flags)

            # Likelihood contribution
            ipred = compute_individual_predictions(
                theta, eta, times, doses, model_spec, grid, solver
            )

            for (j, (y, f)) in enumerate(zip(obs, ipred))
                is_blq = subj_has_blq && j <= length(blq_flags) && blq_flags[j]

                if is_blq && blq_config !== nothing
                    ofv += observation_log_likelihood_blq(y, f, sigma, blq_config, true, lloq)
                else
                    ofv += observation_log_likelihood(y, f, sigma)
                end
            end
        end

        # Add log determinant term
        ofv += n_subj * logdet(omega_new)

        return ofv
    end

    # Optimize theta with fallback optimizer
    opt_config = OptimizerConfig(verbose=false)
    opt_options = Optim.Options(iterations=100, g_tol=1e-5, show_trace=false)

    opt_result = optimize_bounded_with_fallback(
        theta_objective,
        config.theta_lower,
        config.theta_upper,
        theta_current,
        opt_config;
        options=opt_options
    )

    theta_new = opt_result.minimizer
    final_ofv = opt_result.minimum

    return theta_new, omega_new, final_ofv
end

"""
Compute Laplacian objective function value.
Supports BLQ/censoring handling via optional parameters.
"""
function compute_laplacian_ofv(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    etas::Vector{Vector{Float64}},
    subjects::Vector,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec;
    blq_config::Union{Nothing,BLQConfig}=nothing,
    subject_blq_flags::Vector{Vector{Bool}}=Vector{Vector{Bool}}(),
    subject_lloq::Vector{Float64}=Float64[]
)::Float64
    n_subj = length(subjects)
    omega_inv = inv(omega)
    ofv = 0.0

    # Check if BLQ handling is enabled
    has_blq = blq_config !== nothing && !isempty(subject_blq_flags)

    for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
        eta = etas[i]

        # Prior contribution
        ofv += eta' * omega_inv * eta

        # Get BLQ info for this subject
        blq_flags = has_blq && i <= length(subject_blq_flags) ? subject_blq_flags[i] : Bool[]
        lloq = has_blq && i <= length(subject_lloq) ? subject_lloq[i] : 0.0
        subj_has_blq = has_blq && !isempty(blq_flags)

        # Likelihood contribution
        ipred = compute_individual_predictions(
            theta, eta, times, doses, model_spec, grid, solver
        )

        for (j, (y, f)) in enumerate(zip(obs, ipred))
            is_blq = subj_has_blq && j <= length(blq_flags) && blq_flags[j]

            if is_blq && blq_config !== nothing
                ofv += observation_log_likelihood_blq(y, f, sigma, blq_config, true, lloq)
            else
                ofv += observation_log_likelihood(y, f, sigma)
            end
        end
    end

    # Add log determinant term
    ofv += n_subj * logdet(omega)

    return ofv
end

"""
Compute individual estimates for all subjects.
Supports BLQ/censoring handling via optional parameters.
"""
function compute_individual_estimates_laplacian(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    etas::Vector{Vector{Float64}},
    subjects::Vector,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec;
    blq_config::Union{Nothing,BLQConfig}=nothing,
    subject_blq_flags::Vector{Vector{Bool}}=Vector{Vector{Bool}}(),
    subject_lloq::Vector{Float64}=Float64[]
)::Vector{IndividualEstimate}
    estimates = IndividualEstimate[]
    omega_inv = inv(omega)

    # Check if BLQ handling is enabled
    has_blq = blq_config !== nothing && !isempty(subject_blq_flags)

    for (i, (subj_id, times, obs, doses)) in enumerate(subjects)
        eta = etas[i]

        # Get BLQ info for this subject
        blq_flags = has_blq && i <= length(subject_blq_flags) ? subject_blq_flags[i] : Bool[]
        lloq = has_blq && i <= length(subject_lloq) ? subject_lloq[i] : 0.0
        subj_has_blq = has_blq && !isempty(blq_flags)

        # Compute IPRED
        ipred = compute_individual_predictions(
            theta, eta, times, doses, model_spec, grid, solver
        )

        # Compute PRED (population prediction with eta=0)
        pred = compute_individual_predictions(
            theta, zeros(length(eta)), times, doses, model_spec, grid, solver
        )

        # Compute residuals (CWRES = 0 for BLQ with M1/M3 following NONMEM convention)
        cwres = compute_cwres(obs, ipred, sigma)
        iwres = compute_iwres(obs, ipred, sigma)
        wres = compute_wres(obs, pred, sigma)

        # Set CWRES/IWRES to 0 for BLQ observations with M1 or M3 method
        if subj_has_blq && blq_config !== nothing
            if blq_config.method == BLQ_M1_DISCARD || blq_config.method == BLQ_M3_LIKELIHOOD
                for j in 1:length(blq_flags)
                    if j <= length(cwres) && blq_flags[j]
                        cwres[j] = 0.0
                        iwres[j] = 0.0
                    end
                end
            end
        end

        # OFV contribution with BLQ handling
        ofv_contrib = eta' * omega_inv * eta
        for (j, (y, f)) in enumerate(zip(obs, ipred))
            is_blq = subj_has_blq && j <= length(blq_flags) && blq_flags[j]

            if is_blq && blq_config !== nothing
                ofv_contrib += observation_log_likelihood_blq(y, f, sigma, blq_config, true, lloq)
            else
                ofv_contrib += observation_log_likelihood(y, f, sigma)
            end
        end

        push!(estimates, IndividualEstimate(
            subj_id,
            eta,
            nothing,  # eta_se computed separately if needed
            ipred,
            pred,
            cwres,
            iwres,
            wres,
            ofv_contrib
        ))
    end

    return estimates
end

"""
Compute standard errors using the observed Fisher information.
Supports BLQ/censoring handling via optional parameters.
"""
function compute_standard_errors_laplacian(
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    etas::Vector{Vector{Float64}},
    subjects::Vector,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    config::EstimationConfig;
    blq_config::Union{Nothing,BLQConfig}=nothing,
    subject_blq_flags::Vector{Vector{Bool}}=Vector{Vector{Bool}}(),
    subject_lloq::Vector{Float64}=Float64[]
)
    # Compute numerical Hessian of OFV w.r.t. theta
    n_theta = length(theta)
    h = 1e-5  # Step size for finite differences

    function ofv_wrapper(th)
        return compute_laplacian_ofv(
            th, omega, sigma, etas, subjects, model_spec, grid, solver;
            blq_config=blq_config, subject_blq_flags=subject_blq_flags, subject_lloq=subject_lloq
        )
    end

    # Compute Hessian using central differences
    hessian = zeros(n_theta, n_theta)

    for i in 1:n_theta
        for j in i:n_theta
            if i == j
                # Diagonal element
                theta_plus = copy(theta)
                theta_minus = copy(theta)
                theta_plus[i] += h
                theta_minus[i] -= h

                f_plus = ofv_wrapper(theta_plus)
                f_center = ofv_wrapper(theta)
                f_minus = ofv_wrapper(theta_minus)

                hessian[i, i] = (f_plus - 2*f_center + f_minus) / h^2
            else
                # Off-diagonal element
                theta_pp = copy(theta)
                theta_pm = copy(theta)
                theta_mp = copy(theta)
                theta_mm = copy(theta)

                theta_pp[i] += h; theta_pp[j] += h
                theta_pm[i] += h; theta_pm[j] -= h
                theta_mp[i] -= h; theta_mp[j] += h
                theta_mm[i] -= h; theta_mm[j] -= h

                f_pp = ofv_wrapper(theta_pp)
                f_pm = ofv_wrapper(theta_pm)
                f_mp = ofv_wrapper(theta_mp)
                f_mm = ofv_wrapper(theta_mm)

                hessian[i, j] = (f_pp - f_pm - f_mp + f_mm) / (4 * h^2)
                hessian[j, i] = hessian[i, j]
            end
        end
    end

    # Compute condition number and eigenvalue ratio
    eigenvalues = eigvals(hessian)
    real_eigenvalues = real.(eigenvalues)

    if any(real_eigenvalues .<= 0)
        # Hessian not positive definite
        return nothing, nothing, nothing, false, NaN, NaN
    end

    condition_num = maximum(real_eigenvalues) / minimum(real_eigenvalues)
    eigenvalue_ratio = minimum(real_eigenvalues) / maximum(real_eigenvalues)

    # Invert Hessian to get covariance matrix
    try
        cov_matrix = inv(hessian)
        variances = diag(cov_matrix)

        if any(variances .< 0)
            return nothing, nothing, nothing, false, condition_num, eigenvalue_ratio
        end

        theta_se = sqrt.(variances)

        # For omega and sigma SEs, we would need to expand the parameter vector
        # For now, return nothing for these
        return theta_se, nothing, nothing, true, condition_num, eigenvalue_ratio
    catch e
        return nothing, nothing, nothing, false, condition_num, eigenvalue_ratio
    end
end

"""
Ensure matrix is positive definite by adjusting small eigenvalues.
"""
function ensure_positive_definite(A::Matrix{Float64}; min_eigenvalue::Float64=1e-6)::Matrix{Float64}
    eigenvals, eigenvecs = eigen(Symmetric(A))
    eigenvals = max.(eigenvals, min_eigenvalue)
    return eigenvecs * Diagonal(eigenvals) * eigenvecs'
end

"""
Convert covariance matrix to correlation matrix.
"""
function cov_to_corr(cov::Matrix{Float64})::Matrix{Float64}
    d = sqrt.(diag(cov))
    corr = cov ./ (d * d')
    # Fix diagonal to exactly 1
    for i in 1:size(corr, 1)
        corr[i, i] = 1.0
    end
    return corr
end

export ensure_positive_definite, cov_to_corr
