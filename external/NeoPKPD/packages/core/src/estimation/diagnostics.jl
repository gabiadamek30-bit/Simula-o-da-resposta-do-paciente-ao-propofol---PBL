# Diagnostics for Parameter Estimation
# Implements additional residual computations and model evaluation metrics

using LinearAlgebra

# Note: compute_cwres and compute_iwres are already defined in residual_error.jl
# We just export additional diagnostic functions here
export compute_wres, compute_npde
export shrinkage_eta, shrinkage_epsilon, vpc_check

"""
Compute Weighted Residuals (WRES).

Uses population predictions (eta = 0).

Arguments:
- obs: Observed values
- pred: Population predictions (eta = 0)
- sigma: Residual error specification

Returns:
- Vector of WRES
"""
function compute_wres(
    obs::Vector{Float64},
    pred::Vector{Float64},
    sigma::ResidualErrorSpec
)::Vector{Float64}
    # WRES uses the same computation as IWRES/CWRES
    return compute_iwres(obs, pred, sigma)
end

"""
Compute Normalized Prediction Distribution Errors (NPDE).

NPDE are simulation-based diagnostics that should follow N(0,1).

Arguments:
- obs: Observed value
- simulated: Matrix of simulated values (n_sim x 1)
- rng: Random number generator

Returns:
- NPDE value
"""
function compute_npde_single(
    obs::Float64,
    simulated::Vector{Float64},
    rng
)::Float64
    n_sim = length(simulated)

    # Compute cumulative distribution
    pde = sum(simulated .<= obs) / n_sim

    # Add small random jitter to avoid exact 0 or 1
    pde = pde + rand(rng) / n_sim

    # Clamp to valid range
    pde = clamp(pde, 1e-10, 1.0 - 1e-10)

    # Transform to normal
    return quantile(Normal(), pde)
end

function compute_npde(
    obs::Vector{Float64},
    simulated::Matrix{Float64},  # n_obs x n_sim
    rng
)::Vector{Float64}
    n_obs = length(obs)
    npde = zeros(n_obs)

    for i in 1:n_obs
        npde[i] = compute_npde_single(obs[i], simulated[i, :], rng)
    end

    return npde
end

"""
Compute eta shrinkage.

Shrinkage indicates how much information the data provides about
individual random effects. High shrinkage (>30%) suggests the
individual estimates are regressing toward the population mean.

shrinkage = 1 - SD(eta_hat) / omega

Arguments:
- etas: Matrix of estimated etas (n_subj x n_eta)
- omega: Omega matrix (variance of random effects)

Returns:
- Vector of shrinkage values (one per eta)
"""
function shrinkage_eta(
    etas::Vector{Vector{Float64}},
    omega::Matrix{Float64}
)::Vector{Float64}
    n_subj = length(etas)
    n_eta = length(etas[1])

    # Compute empirical SD of each eta
    eta_matrix = hcat(etas...)'  # n_subj x n_eta

    shrinkage = zeros(n_eta)
    for j in 1:n_eta
        eta_j = eta_matrix[:, j]
        sd_empirical = sqrt(var_sample(eta_j))
        sd_theoretical = sqrt(omega[j, j])

        if sd_theoretical > 0
            shrinkage[j] = 1.0 - sd_empirical / sd_theoretical
        else
            shrinkage[j] = NaN
        end
    end

    return shrinkage
end

"""
Compute epsilon shrinkage.

Epsilon shrinkage measures the information content for residual error.

shrinkage = 1 - SD(IWRES)

Arguments:
- iwres_all: All IWRES values from all subjects

Returns:
- Epsilon shrinkage value
"""
function shrinkage_epsilon(
    iwres_all::Vector{Float64}
)::Float64
    # Filter out NaN values
    valid_iwres = filter(!isnan, iwres_all)

    if length(valid_iwres) < 2
        return NaN
    end

    sd_iwres = sqrt(var_sample(valid_iwres))
    return 1.0 - sd_iwres
end

"""
Sample variance computation.
"""
function var_sample(x::Vector{Float64})::Float64
    n = length(x)
    if n < 2
        return 0.0
    end
    m = sum(x) / n
    return sum((x .- m).^2) / (n - 1)
end

export shrinkage_eta, shrinkage_epsilon

"""
VPC-based model check statistics.

Computes the percentage of observed data falling within
the prediction intervals.

Arguments:
- obs: Observed values
- pi_lower: Lower prediction interval
- pi_upper: Upper prediction interval

Returns:
- Percentage within PI
"""
function vpc_check(
    obs::Vector{Float64},
    pi_lower::Vector{Float64},
    pi_upper::Vector{Float64}
)::Float64
    n = length(obs)
    within = sum((obs .>= pi_lower) .& (obs .<= pi_upper))
    return 100.0 * within / n
end

"""
Compute objective function value for model comparison.
"""
function compute_ofv(
    result::EstimationResult
)::Float64
    return result.ofv
end

export compute_ofv

"""
Likelihood ratio test between two nested models.

Arguments:
- ofv_full: OFV of full model
- ofv_reduced: OFV of reduced model
- df: Degrees of freedom difference

Returns:
- (chi_sq, p_value) tuple
"""
function likelihood_ratio_test(
    ofv_full::Float64,
    ofv_reduced::Float64,
    df::Int
)::Tuple{Float64,Float64}
    chi_sq = ofv_reduced - ofv_full  # Should be positive if full is better

    if chi_sq < 0
        chi_sq = 0.0
    end

    p_value = 1.0 - cdf(Chisq(df), chi_sq)

    return (chi_sq, p_value)
end

export likelihood_ratio_test

"""
Compute log-likelihood from OFV.
"""
function loglikelihood_from_ofv(ofv::Float64)::Float64
    return -ofv / 2
end

export loglikelihood_from_ofv

"""
Summary statistics for residuals.
"""
function residual_summary(residuals::Vector{Float64})::Dict{Symbol,Float64}
    valid = filter(!isnan, residuals)
    n = length(valid)

    if n == 0
        return Dict{Symbol,Float64}()
    end

    m = sum(valid) / n
    v = n > 1 ? sum((valid .- m).^2) / (n - 1) : 0.0
    sd = sqrt(v)

    sorted = sort(valid)
    median_val = n % 2 == 1 ? sorted[(n+1)÷2] : (sorted[n÷2] + sorted[n÷2+1]) / 2

    return Dict(
        :mean => m,
        :sd => sd,
        :median => median_val,
        :min => minimum(valid),
        :max => maximum(valid),
        :n => Float64(n)
    )
end

export residual_summary
