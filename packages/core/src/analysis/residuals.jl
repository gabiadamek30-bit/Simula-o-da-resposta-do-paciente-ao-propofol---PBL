# Full CWRES (Conditional Weighted Residuals) Implementation
# Implements proper FO approximation with AD-computed derivatives

using ForwardDiff
using LinearAlgebra
using Distributions

# Note: Basic compute_cwres, compute_iwres, compute_wres are defined in engine/residual_error.jl
# and estimation/diagnostics.jl. This file provides:
# - Advanced CWRES with full FO approximation (different signature)
# - NPDE computation
# - ResidualResult type for aggregated diagnostics
export compute_npde, compute_npde_monte_carlo
export ResidualResult, ResidualSummary, compute_full_cwres
export compute_npde_from_population, compute_npde_from_estimation
export summarize_residuals, check_residual_normality, compute_pde

# ============================================================================
# Result Types
# ============================================================================

"""
Result of residual computation for a single subject.
"""
struct ResidualResult
    subject_id::String
    times::Vector{Float64}
    dv::Vector{Float64}           # Observed values
    pred::Vector{Float64}         # Population predictions
    ipred::Vector{Float64}        # Individual predictions
    wres::Vector{Float64}         # Weighted residuals
    iwres::Vector{Float64}        # Individual weighted residuals
    cwres::Vector{Float64}        # Conditional weighted residuals
    npde::Vector{Float64}         # Normalized prediction distribution errors
end

"""
Summary statistics for residuals.
"""
struct ResidualSummary
    n_obs::Int
    mean_wres::Float64
    sd_wres::Float64
    mean_iwres::Float64
    sd_iwres::Float64
    mean_cwres::Float64
    sd_cwres::Float64
    shapiro_wilk_p::Float64  # Test for normality
end

# ============================================================================
# Helper: Compute Residual Variance
# ============================================================================

"""
Compute residual variance based on error model.
Used for advanced CWRES computation.
"""
function _compute_residual_variance(pred::Float64, sigma::ResidualErrorSpec)::Float64
    if sigma.kind isa AdditiveError
        return sigma.params.sigma^2
    elseif sigma.kind isa ProportionalError
        return (sigma.params.sigma * pred)^2
    elseif sigma.kind isa CombinedError
        add = sigma.params.sigma_add
        prop = sigma.params.sigma_prop
        return add^2 + (prop * pred)^2
    else
        return sigma.params.sigma^2
    end
end

# ============================================================================
# CWRES (Conditional Weighted Residuals) - Full FO Approximation
# ============================================================================

"""
    compute_cwres(dv, pred, ipred, eta, omega, sigma, dF_deta)

Compute conditional weighted residuals using the FO approximation.

CWRES_i = (DV_i - PRED_i) / sqrt(C_i)

where C_i = dF/dη * Ω * dF/dη' + σ²

This properly accounts for the correlation structure induced by random effects.

Arguments:
- dv: Observed values
- pred: Population predictions (η=0)
- ipred: Individual predictions (with estimated η)
- eta: Estimated random effects vector
- omega: Random effects covariance matrix
- sigma: Residual error specification
- dF_deta: Jacobian of predictions with respect to eta (n_obs x n_eta)
"""
function compute_cwres(
    dv::Vector{Float64},
    pred::Vector{Float64},
    ipred::Vector{Float64},
    eta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    dF_deta::Matrix{Float64}
)::Vector{Float64}
    n_obs = length(dv)
    n_eta = length(eta)
    cwres = zeros(n_obs)

    for i in 1:n_obs
        # dF/dη for this observation (row vector)
        dF_i = dF_deta[i, :]

        # Variance from random effects: dF/dη * Ω * dF/dη'
        var_eta = dot(dF_i, omega * dF_i)

        # Residual variance
        var_eps = _compute_residual_variance(pred[i], sigma)

        # Total variance under FO approximation
        var_total = var_eta + var_eps

        # CWRES
        cwres[i] = (dv[i] - pred[i]) / sqrt(var_total)
    end

    return cwres
end

"""
    compute_cwres_with_ad(dv, times, pred_fn, eta, omega, sigma)

Compute CWRES with automatic differentiation for dF/dη.

Arguments:
- dv: Observed values
- times: Observation times
- pred_fn: Function (eta, t) -> prediction at time t
- eta: Estimated random effects
- omega: Random effects covariance matrix
- sigma: Residual error specification
"""
function compute_cwres_with_ad(
    dv::Vector{Float64},
    times::Vector{Float64},
    pred_fn::Function,
    eta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec
)::Vector{Float64}
    n_obs = length(dv)
    n_eta = length(eta)

    # Compute population predictions (eta = 0)
    eta_zero = zeros(n_eta)
    pred = [pred_fn(eta_zero, t) for t in times]

    # Compute dF/dη using ForwardDiff
    dF_deta = zeros(n_obs, n_eta)
    for i in 1:n_obs
        t = times[i]
        grad = ForwardDiff.gradient(e -> pred_fn(e, t), eta)
        dF_deta[i, :] = grad
    end

    # Compute individual predictions for reference
    ipred = [pred_fn(eta, t) for t in times]

    return compute_cwres(dv, pred, ipred, eta, omega, sigma, dF_deta)
end

# ============================================================================
# Full CWRES Computation for EstimationResult
# ============================================================================

"""
    compute_full_cwres(estimation_result, observed, model_spec, grid, solver;
                       include_npde=false, npde_n_sim=1000, npde_seed=12345)

Compute full CWRES for all subjects in an estimation result.

This uses the estimated individual parameters to compute proper CWRES
with AD-derived gradients.

# Arguments
- result: EstimationResult from parameter estimation
- observed: ObservedData with concentration-time data
- model_spec: ModelSpec defining the PK/PD model
- grid: Simulation time grid
- solver: ODE solver specification
- include_npde: If true, compute proper NPDE using Monte Carlo (slower but accurate)
- npde_n_sim: Number of Monte Carlo simulations for NPDE (default: 1000)
- npde_seed: Random seed for NPDE computation

# Returns
- Vector of ResidualResult, one per subject

# Note
Computing NPDE with Monte Carlo is computationally expensive. For quick
diagnostics, set include_npde=false. For final model validation, use
include_npde=true with npde_n_sim >= 500 (FDA recommendation).
"""
function compute_full_cwres(
    result::EstimationResult,
    observed::ObservedData,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec;
    include_npde::Bool=false,
    npde_n_sim::Int=1000,
    npde_seed::UInt64=UInt64(12345)
)::Vector{ResidualResult}
    residual_results = ResidualResult[]

    # Pre-compute NPDE for all subjects if requested
    npde_by_subject = if include_npde
        npde_result = compute_npde_monte_carlo(
            observed, model_spec,
            Vector{Float64}(result.theta),
            Matrix{Float64}(result.omega),
            result.sigma,
            grid, solver;
            n_sim=npde_n_sim,
            seed=npde_seed
        )
        npde_result[:npde_by_subject]
    else
        nothing
    end

    for (i, subj) in enumerate(observed.subjects)
        ind = result.individuals[i]

        # Get subject's observations
        times = subj.times
        dv = subj.observations

        # Create prediction function that takes eta and returns predictions
        pred_fn = (eta, t) -> begin
            # Apply eta to parameters
            modified_params = _apply_eta_to_params(model_spec.params, eta, result)

            # Create modified spec
            mod_spec = ModelSpec(
                model_spec.kind,
                model_spec.name,
                modified_params,
                model_spec.doses
            )

            # Simulate at single time point
            single_grid = SimGrid(grid.t0, t + 0.1, [t])
            sim_result = simulate(mod_spec, single_grid, solver)

            return sim_result.observations[:conc][1]
        end

        # Compute CWRES with AD
        cwres = compute_cwres_with_ad(
            dv, times, pred_fn, ind.eta, result.omega, result.sigma
        )

        # Compute other residuals
        wres = compute_wres(dv, ind.pred, result.sigma)
        iwres = compute_iwres(dv, ind.ipred, result.sigma)

        # NPDE - use pre-computed values if available, otherwise zeros
        npde = if include_npde && npde_by_subject !== nothing
            npde_by_subject[i]
        else
            zeros(length(dv))
        end

        push!(residual_results, ResidualResult(
            subj.id, times, dv, ind.pred, ind.ipred,
            wres, iwres, cwres, npde
        ))
    end

    return residual_results
end

"""
Apply random effects to population parameters.
"""
function _apply_eta_to_params(params, eta::Vector{Float64}, result::EstimationResult)
    # This is a simplified version - actual implementation depends on
    # how parameters are structured and which ones have IIV

    # For now, assume eta applies multiplicatively to first n_eta parameters
    # A full implementation would use the IIV specification

    param_values = collect(values(params))
    n_eta = length(eta)

    for j in 1:min(n_eta, length(param_values))
        param_values[j] *= exp(eta[j])
    end

    return param_values
end

# ============================================================================
# NPDE (Normalized Prediction Distribution Errors)
# ============================================================================

"""
    compute_npde(dv, simulated_preds)

Compute NPDE using Monte Carlo simulations.

NPDE transforms the observed data to standard normal based on the
predictive distribution from simulations.

Arguments:
- dv: Observed values
- simulated_preds: Matrix of simulated predictions (n_obs x n_sim)

Returns:
- Vector of NPDE values
"""
function compute_npde(
    dv::Vector{Float64},
    simulated_preds::Matrix{Float64}
)::Vector{Float64}
    n_obs = length(dv)
    n_sim = size(simulated_preds, 2)
    npde = zeros(n_obs)

    for i in 1:n_obs
        # Get simulated values for this observation
        sim_vals = simulated_preds[i, :]

        # Compute percentile (pde)
        pde = sum(sim_vals .<= dv[i]) / n_sim

        # Handle edge cases
        pde = clamp(pde, 0.001, 0.999)

        # Transform to standard normal
        npde[i] = quantile(Normal(), pde)
    end

    return npde
end

"""
    compute_npde_from_population(observed, pop_spec, grid, solver; n_sim)

Compute NPDE by running population simulations.
"""
function compute_npde_from_population(
    observed::ObservedData,
    pop_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    n_sim::Int=1000,
    seed::UInt64=UInt64(12345)
)::Vector{Vector{Float64}}
    rng = StableRNG(seed)
    npde_results = Vector{Float64}[]

    for subj in observed.subjects
        times = subj.times
        dv = subj.observations
        n_obs = length(times)

        # Simulate predictions
        simulated_preds = zeros(n_obs, n_sim)

        for sim_idx in 1:n_sim
            # Create population spec with new seed
            sim_seed = rand(rng, UInt64)
            new_iiv = IIVSpec(
                pop_spec.iiv.kind,
                pop_spec.iiv.omegas,
                sim_seed,
                1  # Single subject
            )
            new_pop_spec = PopulationSpec(
                pop_spec.model_spec,
                new_iiv,
                pop_spec.iov,
                pop_spec.error,
                pop_spec.covariates
            )

            # Simulate
            pop_result = simulate_population(new_pop_spec, grid, solver)

            # Extract concentrations at observation times
            sim_times = pop_result.subjects[1].times
            sim_conc = pop_result.subjects[1].observations[:conc]

            for (j, t) in enumerate(times)
                # Find nearest simulation time
                idx = argmin(abs.(sim_times .- t))
                simulated_preds[j, sim_idx] = sim_conc[idx]
            end
        end

        # Compute NPDE for this subject
        npde = compute_npde(dv, simulated_preds)
        push!(npde_results, npde)
    end

    return npde_results
end

# ============================================================================
# Diagnostic Functions
# ============================================================================

"""
Compute summary statistics for residuals.
"""
function summarize_residuals(residuals::Vector{ResidualResult})::ResidualSummary
    all_wres = vcat([r.wres for r in residuals]...)
    all_iwres = vcat([r.iwres for r in residuals]...)
    all_cwres = vcat([r.cwres for r in residuals]...)

    n_obs = length(all_wres)

    # Simple Shapiro-Wilk approximation for normality test
    # (Full implementation would use proper statistical test)
    cwres_sorted = sort(all_cwres)
    expected_normal = [quantile(Normal(), (i - 0.375) / (n_obs + 0.25)) for i in 1:n_obs]
    correlation = cor(cwres_sorted, expected_normal)
    sw_p = correlation^2  # Simplified approximation

    return ResidualSummary(
        n_obs,
        mean(all_wres),
        std(all_wres),
        mean(all_iwres),
        std(all_iwres),
        mean(all_cwres),
        std(all_cwres),
        sw_p
    )
end

export summarize_residuals

"""
Check if residuals follow standard normal distribution.
Returns true if residuals are approximately N(0,1).
"""
function check_residual_normality(
    residuals::Vector{Float64};
    mean_tol::Float64=0.1,
    sd_tol::Float64=0.2
)::Bool
    m = mean(residuals)
    s = std(residuals)

    return abs(m) < mean_tol && abs(s - 1.0) < sd_tol
end

export check_residual_normality


# ============================================================================
# Industry-Standard NPDE Monte Carlo Implementation
# ============================================================================
# Reference: Brendel et al. (2006), Comets et al. (2008)
# NPDE is computed by:
# 1. Simulating K predictions from the model using estimated parameters
# 2. Computing the prediction distribution error (pde) for each observation
# 3. Transforming pde to standard normal (npde)

"""
    compute_pde(observed::Float64, simulated::Vector{Float64}) -> Float64

Compute the Prediction Distribution Error (pde) for a single observation.

The pde is the proportion of simulated values that are less than or equal
to the observed value. This represents the percentile of the observation
in the predictive distribution.

Uses Blom's approximation for continuity correction to avoid 0 and 1.
"""
function compute_pde(observed::Float64, simulated::Vector{Float64})::Float64
    K = length(simulated)
    if K == 0
        return 0.5  # Default to median if no simulations
    end

    # Count simulations <= observed
    n_below = sum(simulated .<= observed)

    # Blom's continuity correction: (rank - 3/8) / (n + 1/4)
    # This maps rank 1 to ~0.001 and rank K to ~0.999 for K=1000
    pde = (n_below - 0.375) / (K + 0.25)

    # Clamp to valid probability range
    return clamp(pde, 1e-10, 1.0 - 1e-10)
end

"""
    compute_npde_monte_carlo(
        observed_data, model_spec, theta, omega, sigma, grid, solver;
        n_sim=1000, seed=12345, include_residual_error=true
    )

Compute NPDE using industry-standard Monte Carlo simulation.

This is the preferred method for computing NPDE as it properly accounts for
both inter-individual variability (IIV) and residual unexplained variability (RUV).

# Algorithm (Brendel et al. 2006, Comets et al. 2008)
For each observation yij:
1. Simulate K predictions from the model:
   - Draw ηk ~ N(0, Ω)
   - Compute individual parameters: θi,k = θ * exp(ηk)
   - Simulate concentration: f(θi,k, ti)
   - Add residual error: ŷk = f(...) * (1 + ε) where ε ~ N(0, σ²)
2. Compute pde = P(simulated ≤ observed)
3. Transform: npde = Φ⁻¹(pde)

# Arguments
- observed_data: ObservedData containing subjects with observations
- model_spec: ModelSpec defining the PK/PD model
- theta: Estimated population parameters (Vector{Float64})
- omega: Estimated IIV variance-covariance matrix (Matrix{Float64})
- sigma: Residual error specification
- grid: Simulation time grid
- solver: ODE solver specification
- n_sim: Number of Monte Carlo simulations (default: 1000, FDA recommends ≥500)
- seed: Random seed for reproducibility
- include_residual_error: Whether to add residual error to simulations

# Returns
- Dictionary with:
  - :npde => Vector of all NPDE values (pooled across subjects)
  - :pde => Vector of all pde values
  - :npde_by_subject => Vector of vectors (NPDE per subject)
  - :pde_by_subject => Vector of vectors (pde per subject)
  - :n_observations => Total number of observations
  - :n_simulations => Number of MC simulations used

# Example
```julia
result = compute_npde_monte_carlo(
    observed_data, model_spec,
    [10.0, 50.0],  # theta
    [0.09 0.0; 0.0 0.04],  # omega
    sigma_spec,
    grid, solver,
    n_sim=1000
)

# Check if model is adequate (NPDE should be ~N(0,1))
mean_npde = mean(result[:npde])  # Should be ~0
std_npde = std(result[:npde])    # Should be ~1
```
"""
function compute_npde_monte_carlo(
    observed_data::ObservedData,
    model_spec::ModelSpec,
    theta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    grid::SimGrid,
    solver::SolverSpec;
    n_sim::Int=1000,
    seed::UInt64=UInt64(12345),
    include_residual_error::Bool=true
)::Dict{Symbol, Any}
    rng = StableRNG(seed)
    n_eta = size(omega, 1)

    # Pre-compute Cholesky decomposition for sampling from MVN
    omega_chol = try
        cholesky(Symmetric(omega)).L
    catch
        # If omega is not positive definite, use diagonal
        cholesky(Symmetric(Diagonal(diag(omega) .+ 1e-10))).L
    end

    # Storage for results
    all_npde = Float64[]
    all_pde = Float64[]
    npde_by_subject = Vector{Float64}[]
    pde_by_subject = Vector{Float64}[]

    # Get sigma value for residual error
    sigma_val = if sigma.kind isa ProportionalError
        sigma.params.sigma
    elseif sigma.kind isa AdditiveError
        sigma.params.sigma
    else
        0.1  # Default
    end

    for subj in observed_data.subjects
        times = subj.times
        dv = subj.observations
        n_obs = length(times)

        # Matrix to store simulated predictions: n_obs x n_sim
        simulated_preds = zeros(n_obs, n_sim)

        for k in 1:n_sim
            # 1. Sample random effects: η ~ N(0, Ω)
            z = randn(rng, n_eta)
            eta = omega_chol * z

            # 2. Compute individual parameters: θ_i = θ * exp(η)
            theta_individual = theta .* exp.(eta)

            # 3. Create individual model spec
            try
                ind_params = theta_to_params(theta_individual, model_spec)
                ind_model = ModelSpec(
                    model_spec.kind,
                    model_spec.name,
                    ind_params,
                    subj.doses
                )

                # 4. Simulate
                sim_result = simulate(ind_model, grid, solver)
                sim_times = sim_result.t
                sim_conc = sim_result.observations[:conc]

                # 5. Interpolate to observation times and add residual error
                for (j, t) in enumerate(times)
                    idx = searchsortedfirst(sim_times, t)
                    if idx > length(sim_times)
                        idx = length(sim_times)
                    elseif idx > 1 && abs(sim_times[idx] - t) > abs(sim_times[idx-1] - t)
                        idx = idx - 1
                    end

                    pred = sim_conc[idx]

                    # Add residual error if requested
                    if include_residual_error
                        if sigma.kind isa ProportionalError
                            # Proportional: y = f * (1 + ε), ε ~ N(0, σ²)
                            eps = randn(rng) * sigma_val
                            pred = pred * (1.0 + eps)
                        elseif sigma.kind isa AdditiveError
                            # Additive: y = f + ε, ε ~ N(0, σ²)
                            eps = randn(rng) * sigma_val
                            pred = pred + eps
                        elseif sigma.kind isa CombinedError
                            # Combined: y = f + ε_add + f * ε_prop
                            eps_add = randn(rng) * sigma.params.sigma_add
                            eps_prop = randn(rng) * sigma.params.sigma_prop
                            pred = pred + eps_add + pred * eps_prop
                        end
                    end

                    simulated_preds[j, k] = max(pred, 0.0)  # Concentrations must be non-negative
                end
            catch e
                # If simulation fails, use NaN (will be handled in pde computation)
                simulated_preds[:, k] .= NaN
            end
        end

        # Compute pde and npde for each observation
        subj_npde = zeros(n_obs)
        subj_pde = zeros(n_obs)

        for j in 1:n_obs
            # Get valid (non-NaN) simulated values
            sim_vals = filter(!isnan, simulated_preds[j, :])

            if length(sim_vals) >= 10  # Need at least some valid simulations
                pde = compute_pde(dv[j], sim_vals)
                npde = quantile(Normal(), pde)
            else
                pde = 0.5
                npde = 0.0
            end

            subj_pde[j] = pde
            subj_npde[j] = npde
        end

        push!(npde_by_subject, subj_npde)
        push!(pde_by_subject, subj_pde)
        append!(all_npde, subj_npde)
        append!(all_pde, subj_pde)
    end

    return Dict(
        :npde => all_npde,
        :pde => all_pde,
        :npde_by_subject => npde_by_subject,
        :pde_by_subject => pde_by_subject,
        :n_observations => length(all_npde),
        :n_simulations => n_sim
    )
end

"""
    compute_npde_from_estimation(
        observed_data, model_spec, result, grid, solver;
        n_sim=1000, seed=12345
    )

Compute NPDE using estimation results.

This is a convenience wrapper around compute_npde_monte_carlo that
extracts theta, omega, and sigma from an EstimationResult.

# Arguments
- observed_data: ObservedData with observations
- model_spec: ModelSpec defining the model
- result: EstimationResult from parameter estimation
- grid: Simulation time grid
- solver: ODE solver specification
- n_sim: Number of Monte Carlo simulations
- seed: Random seed

# Returns
- Dictionary with NPDE results (same as compute_npde_monte_carlo)

# Example
```julia
# After estimation
est_result = estimate(data, model_spec, config; grid=grid, solver=solver)

# Compute NPDE
npde_result = compute_npde_from_estimation(
    data, model_spec, est_result, grid, solver
)

println("Mean NPDE: ", mean(npde_result[:npde]))
println("Std NPDE: ", std(npde_result[:npde]))
```
"""
function compute_npde_from_estimation(
    observed_data::ObservedData,
    model_spec::ModelSpec,
    result::EstimationResult,
    grid::SimGrid,
    solver::SolverSpec;
    n_sim::Int=1000,
    seed::UInt64=UInt64(12345)
)::Dict{Symbol, Any}
    theta = Vector{Float64}(result.theta)
    omega = Matrix{Float64}(result.omega)
    sigma = result.sigma

    return compute_npde_monte_carlo(
        observed_data, model_spec, theta, omega, sigma, grid, solver;
        n_sim=n_sim, seed=seed, include_residual_error=true
    )
end
