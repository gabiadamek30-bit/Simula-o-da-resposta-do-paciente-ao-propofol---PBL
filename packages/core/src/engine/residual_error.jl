# Residual Error Implementation
# Functions for applying error, computing likelihoods, and diagnostics

using Random
using Distributions
using StableRNGs

export apply_residual_error, residual_log_likelihood, compute_iwres, compute_cwres
export residual_variance, residual_sd

# -------------------------
# Variance/SD Computation
# -------------------------

"""
Compute the residual variance for a given prediction and error model.
"""
function residual_variance(F::Float64, spec::ResidualErrorSpec{AdditiveError})::Float64
    return spec.params.sigma^2
end

function residual_variance(F::Float64, spec::ResidualErrorSpec{ProportionalError})::Float64
    return (spec.params.sigma * F)^2
end

function residual_variance(F::Float64, spec::ResidualErrorSpec{CombinedError})::Float64
    return spec.params.sigma_add^2 + (spec.params.sigma_prop * F)^2
end

function residual_variance(F::Float64, spec::ResidualErrorSpec{ExponentialError})::Float64
    # For exponential model, variance is on log scale
    return spec.params.sigma^2
end

"""
Compute the residual standard deviation for a given prediction.
"""
function residual_sd(F::Float64, spec::ResidualErrorSpec)::Float64
    return sqrt(residual_variance(F, spec))
end

# -------------------------
# Error Application
# -------------------------

"""
Apply residual error to a vector of predictions.

Returns a new vector with error applied. The original predictions are not modified.

Arguments:
- F: Vector of predictions (IPRED)
- spec: Residual error specification
- rng: Random number generator (optional, uses seeded RNG from spec if not provided)

Returns:
- Vector of observations with error applied (simulated DV)
"""
function apply_residual_error(
    F::Vector{Float64},
    spec::ResidualErrorSpec{AdditiveError};
    rng::Union{Nothing,AbstractRNG}=nothing
)::Vector{Float64}
    rng = isnothing(rng) ? StableRNG(spec.seed) : rng
    σ = spec.params.sigma
    return F .+ randn(rng, length(F)) .* σ
end

function apply_residual_error(
    F::Vector{Float64},
    spec::ResidualErrorSpec{ProportionalError};
    rng::Union{Nothing,AbstractRNG}=nothing
)::Vector{Float64}
    rng = isnothing(rng) ? StableRNG(spec.seed) : rng
    σ = spec.params.sigma
    return F .* (1.0 .+ randn(rng, length(F)) .* σ)
end

function apply_residual_error(
    F::Vector{Float64},
    spec::ResidualErrorSpec{CombinedError};
    rng::Union{Nothing,AbstractRNG}=nothing
)::Vector{Float64}
    rng = isnothing(rng) ? StableRNG(spec.seed) : rng
    σ_add = spec.params.sigma_add
    σ_prop = spec.params.sigma_prop

    n = length(F)
    Y = similar(F)
    eps = randn(rng, n)

    for i in 1:n
        sd = sqrt(σ_add^2 + (σ_prop * F[i])^2)
        Y[i] = F[i] + sd * eps[i]
    end

    return Y
end

function apply_residual_error(
    F::Vector{Float64},
    spec::ResidualErrorSpec{ExponentialError};
    rng::Union{Nothing,AbstractRNG}=nothing
)::Vector{Float64}
    rng = isnothing(rng) ? StableRNG(spec.seed) : rng
    σ = spec.params.sigma

    # Ensure predictions are positive for log transform
    F_pos = max.(F, eps(Float64))

    return F_pos .* exp.(randn(rng, length(F)) .* σ)
end

# -------------------------
# Log-Likelihood Computation
# -------------------------

"""
Compute the log-likelihood contribution for observed vs predicted.

This is used in parameter estimation (e.g., FOCE, Laplacian).

Arguments:
- Y: Observed values
- F: Predicted values (IPRED)
- spec: Residual error specification

Returns:
- Total log-likelihood (sum over all observations)
"""
function residual_log_likelihood(
    Y::Vector{Float64},
    F::Vector{Float64},
    spec::ResidualErrorSpec{AdditiveError}
)::Float64
    length(Y) == length(F) || error("Y and F must have same length")

    σ = spec.params.sigma
    n = length(Y)

    ll = 0.0
    for i in 1:n
        ll += logpdf(Normal(F[i], σ), Y[i])
    end

    return ll
end

function residual_log_likelihood(
    Y::Vector{Float64},
    F::Vector{Float64},
    spec::ResidualErrorSpec{ProportionalError}
)::Float64
    length(Y) == length(F) || error("Y and F must have same length")

    σ = spec.params.sigma
    n = length(Y)

    ll = 0.0
    for i in 1:n
        sd = abs(σ * F[i])
        if sd > 0
            ll += logpdf(Normal(F[i], sd), Y[i])
        else
            # Handle F[i] = 0 case
            ll += Y[i] == 0 ? 0.0 : -Inf
        end
    end

    return ll
end

function residual_log_likelihood(
    Y::Vector{Float64},
    F::Vector{Float64},
    spec::ResidualErrorSpec{CombinedError}
)::Float64
    length(Y) == length(F) || error("Y and F must have same length")

    σ_add = spec.params.sigma_add
    σ_prop = spec.params.sigma_prop
    n = length(Y)

    ll = 0.0
    for i in 1:n
        sd = sqrt(σ_add^2 + (σ_prop * F[i])^2)
        ll += logpdf(Normal(F[i], sd), Y[i])
    end

    return ll
end

function residual_log_likelihood(
    Y::Vector{Float64},
    F::Vector{Float64},
    spec::ResidualErrorSpec{ExponentialError}
)::Float64
    length(Y) == length(F) || error("Y and F must have same length")

    σ = spec.params.sigma
    n = length(Y)

    ll = 0.0
    for i in 1:n
        if Y[i] > 0 && F[i] > 0
            # Log-normal: log(Y) ~ Normal(log(F), σ)
            # Include Jacobian: -log(Y)
            ll += logpdf(Normal(log(F[i]), σ), log(Y[i])) - log(Y[i])
        else
            ll += -Inf
        end
    end

    return ll
end

# -------------------------
# Residual Diagnostics
# -------------------------

"""
Compute Individual Weighted Residuals (IWRES).

IWRES = (Y - F) / SD(Y|F)

These are standardized residuals that should be approximately N(0,1)
if the model is correct.
"""
function compute_iwres(
    Y::Vector{Float64},
    F::Vector{Float64},
    spec::ResidualErrorSpec
)::Vector{Float64}
    length(Y) == length(F) || error("Y and F must have same length")

    n = length(Y)
    iwres = similar(Y)

    for i in 1:n
        sd = residual_sd(F[i], spec)
        iwres[i] = sd > 0 ? (Y[i] - F[i]) / sd : 0.0
    end

    return iwres
end

"""
Compute Individual Weighted Residuals for exponential error model.

For exponential error: IWRES = (log(Y) - log(F)) / σ
"""
function compute_iwres(
    Y::Vector{Float64},
    F::Vector{Float64},
    spec::ResidualErrorSpec{ExponentialError}
)::Vector{Float64}
    length(Y) == length(F) || error("Y and F must have same length")

    σ = spec.params.sigma
    n = length(Y)
    iwres = similar(Y)

    for i in 1:n
        if Y[i] > 0 && F[i] > 0
            iwres[i] = (log(Y[i]) - log(F[i])) / σ
        else
            iwres[i] = 0.0
        end
    end

    return iwres
end

"""
Compute Conditional Weighted Residuals (CWRES).

CWRES are a more sophisticated residual metric that accounts for
the correlation structure in mixed-effects models. For population
models, CWRES ≈ IWRES when η is at its mode.

For simple models without random effects, CWRES = IWRES.

Note: Full CWRES computation requires the gradient of predictions
with respect to η, which is implemented in the estimation module.
This function provides a simplified version for diagnostics.
"""
function compute_cwres(
    Y::Vector{Float64},
    F::Vector{Float64},
    spec::ResidualErrorSpec
)::Vector{Float64}
    # Simplified CWRES (equals IWRES for models without random effects)
    return compute_iwres(Y, F, spec)
end
