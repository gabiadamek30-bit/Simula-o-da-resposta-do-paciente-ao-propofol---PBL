# MCP-Mod Implementation
# Multiple Comparison Procedure - Modelling for Phase II Dose-Response Studies
# FDA-qualified method for establishing proof-of-concept and dose selection

export MCPModConfig, MCPModResult, MCPStepResult, ModStepResult
export DoseResponseModel, LinearModel, QuadraticModel, EmaxModel, SigmoidEmaxModel
export ExponentialModel, BetaModel, LogLinearModel
export run_mcpmod, mcp_step, mod_step
export compute_optimal_contrasts, compute_critical_value
export estimate_target_dose, estimate_model_parameters
export model_mean, model_mean_vector, create_mcpmod_models, default_mcpmod_config

using LinearAlgebra
using Statistics
using Distributions
using Random

# ============================================================================
# Dose-Response Model Types
# ============================================================================

"""
Abstract type for dose-response models used in MCP-Mod.
"""
abstract type DoseResponseModel end

"""
    LinearModel

Linear dose-response: f(d) = e0 + δ * d

# Fields
- `name::Symbol`: Model identifier
- `bounds::NamedTuple`: Parameter bounds for δ
"""
struct LinearModel <: DoseResponseModel
    name::Symbol
    bounds::NamedTuple{(:delta_lower, :delta_upper), Tuple{Float64, Float64}}

    function LinearModel(;
            bounds::NamedTuple{(:delta_lower, :delta_upper), Tuple{Float64, Float64}} =
                (delta_lower = -Inf, delta_upper = Inf))
        new(:linear, bounds)
    end
end

"""
    QuadraticModel

Quadratic dose-response: f(d) = e0 + β₁ * d + β₂ * d²

# Fields
- `name::Symbol`: Model identifier
- `delta::Float64`: Location of maximum/minimum (β₁ / (-2β₂))
"""
struct QuadraticModel <: DoseResponseModel
    name::Symbol
    delta::Float64

    function QuadraticModel(; delta::Float64 = 1.0)
        @assert delta != 0 "Delta parameter must be non-zero"
        new(:quadratic, delta)
    end
end

"""
    EmaxModel

Emax dose-response: f(d) = e0 + emax * d / (ed50 + d)

# Fields
- `name::Symbol`: Model identifier
- `ed50::Float64`: Dose producing 50% of maximum effect
"""
struct EmaxModel <: DoseResponseModel
    name::Symbol
    ed50::Float64

    function EmaxModel(; ed50::Float64)
        @assert ed50 > 0 "ED50 must be positive"
        new(:emax, ed50)
    end
end

"""
    SigmoidEmaxModel

Sigmoid Emax dose-response: f(d) = e0 + emax * d^h / (ed50^h + d^h)

# Fields
- `name::Symbol`: Model identifier
- `ed50::Float64`: Dose producing 50% of maximum effect
- `h::Float64`: Hill coefficient (steepness parameter)
"""
struct SigmoidEmaxModel <: DoseResponseModel
    name::Symbol
    ed50::Float64
    h::Float64

    function SigmoidEmaxModel(; ed50::Float64, h::Float64 = 1.0)
        @assert ed50 > 0 "ED50 must be positive"
        @assert h > 0 "Hill coefficient must be positive"
        new(:sigmoid_emax, ed50, h)
    end
end

"""
    ExponentialModel

Exponential dose-response: f(d) = e0 + e1 * (exp(d/δ) - 1)

# Fields
- `name::Symbol`: Model identifier
- `delta::Float64`: Rate parameter
"""
struct ExponentialModel <: DoseResponseModel
    name::Symbol
    delta::Float64

    function ExponentialModel(; delta::Float64)
        @assert delta > 0 "Delta must be positive"
        new(:exponential, delta)
    end
end

"""
    BetaModel

Beta dose-response: f(d) = e0 + emax * B(δ₁, δ₂) * (d/scal)^δ₁ * (1 - d/scal)^δ₂

where B(δ₁, δ₂) is a scaling factor ensuring maximum effect equals emax.

# Fields
- `name::Symbol`: Model identifier
- `delta1::Float64`: First shape parameter
- `delta2::Float64`: Second shape parameter
- `scal::Float64`: Scaling parameter (max dose for normalization)
"""
struct BetaModel <: DoseResponseModel
    name::Symbol
    delta1::Float64
    delta2::Float64
    scal::Float64

    function BetaModel(; delta1::Float64, delta2::Float64, scal::Float64)
        @assert delta1 > 0 "Delta1 must be positive"
        @assert delta2 > 0 "Delta2 must be positive"
        @assert scal > 0 "Scale must be positive"
        new(:beta, delta1, delta2, scal)
    end
end

"""
    LogLinearModel

Log-linear dose-response: f(d) = e0 + δ * log(d + off)

# Fields
- `name::Symbol`: Model identifier
- `off::Float64`: Offset parameter to ensure log is defined at d=0
"""
struct LogLinearModel <: DoseResponseModel
    name::Symbol
    off::Float64

    function LogLinearModel(; off::Float64 = 1.0)
        @assert off > 0 "Offset must be positive"
        new(:loglinear, off)
    end
end

# ============================================================================
# Model Mean Functions
# ============================================================================

"""
    model_mean(model::DoseResponseModel, dose::Float64, params::Dict{Symbol, Float64})

Compute the mean response at a given dose for the specified model.

# Arguments
- `model`: Dose-response model
- `dose`: Dose level
- `params`: Model parameters (e0, emax, etc.)

# Returns
- Mean response value
"""
function model_mean(model::LinearModel, dose::Float64, params::Dict{Symbol, Float64})
    e0 = get(params, :e0, 0.0)
    delta = get(params, :delta, 1.0)
    return e0 + delta * dose
end

function model_mean(model::QuadraticModel, dose::Float64, params::Dict{Symbol, Float64})
    e0 = get(params, :e0, 0.0)
    b1 = get(params, :b1, 1.0)
    b2 = get(params, :b2, -1.0 / (2 * model.delta))
    return e0 + b1 * dose + b2 * dose^2
end

function model_mean(model::EmaxModel, dose::Float64, params::Dict{Symbol, Float64})
    e0 = get(params, :e0, 0.0)
    emax = get(params, :emax, 1.0)
    ed50 = get(params, :ed50, model.ed50)
    return e0 + emax * dose / (ed50 + dose)
end

function model_mean(model::SigmoidEmaxModel, dose::Float64, params::Dict{Symbol, Float64})
    e0 = get(params, :e0, 0.0)
    emax = get(params, :emax, 1.0)
    ed50 = get(params, :ed50, model.ed50)
    h = get(params, :h, model.h)
    return e0 + emax * dose^h / (ed50^h + dose^h)
end

function model_mean(model::ExponentialModel, dose::Float64, params::Dict{Symbol, Float64})
    e0 = get(params, :e0, 0.0)
    e1 = get(params, :e1, 1.0)
    delta = get(params, :delta, model.delta)
    return e0 + e1 * (exp(dose / delta) - 1.0)
end

function model_mean(model::BetaModel, dose::Float64, params::Dict{Symbol, Float64})
    e0 = get(params, :e0, 0.0)
    emax = get(params, :emax, 1.0)
    d1, d2 = model.delta1, model.delta2
    scal = model.scal

    # Normalized dose
    x = dose / scal

    # Scaling factor for beta model
    mode = d1 / (d1 + d2)
    B = mode^d1 * (1 - mode)^d2

    if dose <= 0 || dose >= scal
        return e0
    end

    return e0 + emax * (x^d1 * (1 - x)^d2) / B
end

function model_mean(model::LogLinearModel, dose::Float64, params::Dict{Symbol, Float64})
    e0 = get(params, :e0, 0.0)
    delta = get(params, :delta, 1.0)
    return e0 + delta * log(dose + model.off)
end

"""
    model_mean_vector(model::DoseResponseModel, doses::Vector{Float64}, params::Dict{Symbol, Float64})

Compute mean responses for a vector of doses.
"""
function model_mean_vector(model::DoseResponseModel, doses::Vector{Float64}, params::Dict{Symbol, Float64})
    return [model_mean(model, d, params) for d in doses]
end

# ============================================================================
# Standardized Model Functions (for contrast computation)
# ============================================================================

"""
    standardized_model(model::DoseResponseModel, dose::Float64)

Compute the standardized model function value (e0=0, emax=1 equivalent).
Used for optimal contrast computation.
"""
function standardized_model(model::LinearModel, dose::Float64, dose_max::Float64)
    return dose / dose_max
end

function standardized_model(model::QuadraticModel, dose::Float64, dose_max::Float64)
    delta = model.delta
    # Standardize so that f(0)=0 and range is [0,1]
    f = dose - dose^2 / (2 * delta)
    f_max = dose_max - dose_max^2 / (2 * delta)
    return f / f_max
end

function standardized_model(model::EmaxModel, dose::Float64, dose_max::Float64)
    ed50 = model.ed50
    f = dose / (ed50 + dose)
    f_max = dose_max / (ed50 + dose_max)
    return f / f_max
end

function standardized_model(model::SigmoidEmaxModel, dose::Float64, dose_max::Float64)
    ed50, h = model.ed50, model.h
    f = dose^h / (ed50^h + dose^h)
    f_max = dose_max^h / (ed50^h + dose_max^h)
    return f / f_max
end

function standardized_model(model::ExponentialModel, dose::Float64, dose_max::Float64)
    delta = model.delta
    f = exp(dose / delta) - 1.0
    f_max = exp(dose_max / delta) - 1.0
    return f / f_max
end

function standardized_model(model::BetaModel, dose::Float64, dose_max::Float64)
    d1, d2 = model.delta1, model.delta2
    scal = model.scal

    x = dose / scal
    mode = d1 / (d1 + d2)
    B = mode^d1 * (1 - mode)^d2

    if dose <= 0 || dose >= scal
        return 0.0
    end

    f = (x^d1 * (1 - x)^d2) / B

    # Compute max value
    x_max = dose_max / scal
    if dose_max <= 0 || dose_max >= scal
        f_max = 1.0
    else
        f_max = (x_max^d1 * (1 - x_max)^d2) / B
    end

    return f / max(f_max, 1e-10)
end

function standardized_model(model::LogLinearModel, dose::Float64, dose_max::Float64)
    off = model.off
    f = log(dose + off) - log(off)
    f_max = log(dose_max + off) - log(off)
    return f / f_max
end

# ============================================================================
# Configuration and Result Types
# ============================================================================

"""
    MCPModConfig

Configuration for MCP-Mod analysis.

# Fields
- `doses::Vector{Float64}`: Dose levels including placebo (0)
- `models::Vector{DoseResponseModel}`: Candidate dose-response models
- `alpha::Float64`: Overall significance level (default 0.025 one-sided)
- `delta::Float64`: Clinically relevant effect size for target dose estimation
- `direction::Symbol`: Direction of beneficial effect (:increasing or :decreasing)
- `adjust_method::Symbol`: Method for p-value adjustment (:dunnett, :bonferroni)
- `n_per_dose::Union{Vector{Int}, Nothing}`: Sample sizes per dose (for unbalanced designs)
"""
struct MCPModConfig
    doses::Vector{Float64}
    models::Vector{DoseResponseModel}
    alpha::Float64
    delta::Float64
    direction::Symbol
    adjust_method::Symbol
    n_per_dose::Union{Vector{Int}, Nothing}

    function MCPModConfig(;
            doses::Vector{Float64},
            models::Vector{DoseResponseModel},
            alpha::Float64 = 0.025,
            delta::Float64 = 0.0,
            direction::Symbol = :increasing,
            adjust_method::Symbol = :dunnett,
            n_per_dose::Union{Vector{Int}, Nothing} = nothing)
        @assert 0.0 in doses "Doses must include placebo (0)"
        @assert !isempty(models) "At least one candidate model required"
        @assert 0 < alpha < 1 "Alpha must be between 0 and 1"
        @assert direction in [:increasing, :decreasing] "Direction must be :increasing or :decreasing"
        @assert adjust_method in [:dunnett, :bonferroni] "Adjustment method must be :dunnett or :bonferroni"

        if n_per_dose !== nothing
            @assert length(n_per_dose) == length(doses) "n_per_dose must match doses length"
        end

        new(sort(doses), models, alpha, delta, direction, adjust_method, n_per_dose)
    end
end

"""
    MCPStepResult

Results from the MCP (multiple comparison) step.

# Fields
- `test_statistics::Dict{Symbol, Float64}`: Test statistic for each model
- `p_values::Dict{Symbol, Float64}`: Adjusted p-values for each model
- `critical_value::Float64`: Critical value for the test
- `contrasts::Dict{Symbol, Vector{Float64}}`: Optimal contrasts for each model
- `significant_models::Vector{Symbol}`: Models with significant dose-response
- `max_test_statistic::Float64`: Maximum test statistic across models
- `proof_of_concept::Bool`: Whether PoC was established
"""
struct MCPStepResult
    test_statistics::Dict{Symbol, Float64}
    p_values::Dict{Symbol, Float64}
    critical_value::Float64
    contrasts::Dict{Symbol, Vector{Float64}}
    significant_models::Vector{Symbol}
    max_test_statistic::Float64
    proof_of_concept::Bool
end

"""
    ModStepResult

Results from the Mod (modeling) step.

# Fields
- `fitted_models::Dict{Symbol, Dict{Symbol, Float64}}`: Fitted parameters for each model
- `model_selection::Symbol`: Selected best model
- `aic_values::Dict{Symbol, Float64}`: AIC for each model
- `bic_values::Dict{Symbol, Float64}`: BIC for each model
- `gaic_values::Dict{Symbol, Float64}`: Generalized AIC for each model
- `target_dose::Float64`: Estimated target dose for clinically relevant effect
- `dose_response_curve::Dict{Symbol, Vector{Float64}}`: Predicted responses at each dose
- `residuals::Dict{Symbol, Vector{Float64}}`: Model residuals
"""
struct ModStepResult
    fitted_models::Dict{Symbol, Dict{Symbol, Float64}}
    model_selection::Symbol
    aic_values::Dict{Symbol, Float64}
    bic_values::Dict{Symbol, Float64}
    gaic_values::Dict{Symbol, Float64}
    target_dose::Float64
    dose_response_curve::Dict{Symbol, Vector{Float64}}
    residuals::Dict{Symbol, Vector{Float64}}
end

"""
    MCPModResult

Complete results from MCP-Mod analysis.

# Fields
- `config::MCPModConfig`: Configuration used
- `mcp_result::MCPStepResult`: Results from MCP step
- `mod_result::Union{Nothing, ModStepResult}`: Results from Mod step (if PoC established)
- `observed_means::Vector{Float64}`: Observed means at each dose
- `observed_se::Vector{Float64}`: Standard errors at each dose
- `sample_sizes::Vector{Int}`: Sample sizes at each dose
"""
struct MCPModResult
    config::MCPModConfig
    mcp_result::MCPStepResult
    mod_result::Union{Nothing, ModStepResult}
    observed_means::Vector{Float64}
    observed_se::Vector{Float64}
    sample_sizes::Vector{Int}
end

# ============================================================================
# Optimal Contrast Computation
# ============================================================================

"""
    compute_optimal_contrasts(config::MCPModConfig;
                               n_per_dose::Union{Vector{Int}, Nothing} = nothing)

Compute optimal contrasts for each candidate model.

The optimal contrast c* for testing H0: θ = 0 vs H1: θ > 0 is:
c* = S⁻¹(μ - μ̄1) / ||S⁻¹(μ - μ̄1)||

where μ is the model mean vector, S is the covariance matrix, and 1 is a vector of ones.

# Arguments
- `config`: MCP-Mod configuration
- `n_per_dose`: Sample sizes per dose level (default: balanced)

# Returns
- Dictionary mapping model names to contrast vectors
"""
function compute_optimal_contrasts(config::MCPModConfig;
                                    n_per_dose::Union{Vector{Int}, Nothing} = nothing)
    doses = config.doses
    K = length(doses)
    dose_max = maximum(doses)

    # Use provided n_per_dose or assume balanced
    if n_per_dose === nothing
        n_per_dose = config.n_per_dose
    end
    if n_per_dose === nothing
        n_per_dose = fill(10, K)  # Assume balanced
    end

    # Weight matrix (diagonal with 1/n)
    W = diagm(1.0 ./ n_per_dose)

    contrasts = Dict{Symbol, Vector{Float64}}()

    for model in config.models
        # Compute standardized mean vector
        μ = [standardized_model(model, d, dose_max) for d in doses]

        # Center: subtract weighted mean
        μ_bar = sum(μ .* n_per_dose) / sum(n_per_dose)
        μ_centered = μ .- μ_bar

        # Optimal contrast: c* = S⁻¹(μ - μ̄1) normalized
        # For independent observations: c* ∝ (μ - μ̄1) / n
        c = μ_centered .* n_per_dose

        # Normalize to unit length under S
        norm_val = sqrt(sum(c.^2 ./ n_per_dose))
        if norm_val > 1e-10
            c = c ./ norm_val
        end

        # Ensure contrasts sum to zero
        c = c .- mean(c)

        # Adjust sign based on direction
        if config.direction == :decreasing
            c = -c
        end

        contrasts[model.name] = c
    end

    return contrasts
end

# ============================================================================
# Critical Value Computation
# ============================================================================

"""
    compute_critical_value(config::MCPModConfig, df::Int;
                           n_sim::Int = 100000, seed::Int = 12345)

Compute the critical value for the MCP test using simulation.

For Dunnett-type adjustment, simulates the distribution of the maximum
of correlated t-statistics under the null hypothesis.

# Arguments
- `config`: MCP-Mod configuration
- `df`: Degrees of freedom
- `n_sim`: Number of simulations
- `seed`: Random seed

# Returns
- Critical value for the test
"""
function compute_critical_value(config::MCPModConfig, df::Int;
                                n_sim::Int = 100000, seed::Int = 12345)
    if config.adjust_method == :bonferroni
        # Bonferroni correction
        return quantile(TDist(df), 1 - config.alpha / length(config.models))
    end

    # Dunnett-type adjustment via simulation
    rng = MersenneTwister(seed)

    # Get contrasts
    contrasts = compute_optimal_contrasts(config)
    C = hcat([contrasts[m.name] for m in config.models]...)'  # M x K matrix

    # Compute correlation matrix of test statistics
    # Cor(T_i, T_j) = C_i' S^{-1} C_j / sqrt(C_i' S^{-1} C_i * C_j' S^{-1} C_j)
    K = length(config.doses)
    n_models = length(config.models)

    # Assuming balanced design
    n = config.n_per_dose !== nothing ? config.n_per_dose : fill(10, K)
    S_inv = diagm(n .* 1.0)

    # Correlation matrix
    R = zeros(n_models, n_models)
    for i in 1:n_models
        for j in 1:n_models
            ci = C[i, :]
            cj = C[j, :]
            num = ci' * S_inv * cj
            denom = sqrt((ci' * S_inv * ci) * (cj' * S_inv * cj))
            R[i, j] = num / denom
        end
    end

    # Simulate maximum of correlated t-statistics
    # T = Z / sqrt(χ²/df) where Z ~ MVN(0, R)
    try
        mvn = MvNormal(zeros(n_models), R)
        chi2_dist = Chisq(df)

        max_t = zeros(n_sim)
        for i in 1:n_sim
            z = rand(rng, mvn)
            s = sqrt(rand(rng, chi2_dist) / df)
            t = z ./ s
            max_t[i] = maximum(t)
        end

        return quantile(max_t, 1 - config.alpha)
    catch
        # Fallback to Bonferroni if correlation matrix is singular
        return quantile(TDist(df), 1 - config.alpha / n_models)
    end
end

# ============================================================================
# MCP Step
# ============================================================================

"""
    mcp_step(means::Vector{Float64}, se::Vector{Float64}, n::Vector{Int},
             config::MCPModConfig)

Perform the MCP (multiple comparison procedure) step.

Tests whether there is a significant dose-response relationship for any
of the candidate models.

# Arguments
- `means`: Observed mean responses at each dose
- `se`: Standard errors at each dose
- `n`: Sample sizes at each dose
- `config`: MCP-Mod configuration

# Returns
- `MCPStepResult` with test statistics, p-values, and PoC decision
"""
function mcp_step(means::Vector{Float64}, se::Vector{Float64}, n::Vector{Int},
                  config::MCPModConfig)
    K = length(config.doses)
    @assert length(means) == K "means length must match doses"
    @assert length(se) == K "se length must match doses"
    @assert length(n) == K "n length must match doses"

    # Degrees of freedom
    df = sum(n) - K

    # Pooled variance estimate (weighted by df)
    pooled_var = sum((n .- 1) .* se.^2) / df

    # Compute optimal contrasts
    contrasts = compute_optimal_contrasts(config; n_per_dose = n)

    # Compute critical value
    crit_val = compute_critical_value(config, df; n_sim = 50000)

    # Compute test statistics for each model
    test_stats = Dict{Symbol, Float64}()
    p_values = Dict{Symbol, Float64}()

    for model in config.models
        c = contrasts[model.name]

        # Test statistic: T = c'ȳ / sqrt(σ² * c'Sc)
        # Where S = diag(1/n)
        numerator = sum(c .* means)
        denominator = sqrt(pooled_var * sum(c.^2 ./ n))

        t_stat = numerator / denominator
        test_stats[model.name] = t_stat

        # P-value (one-sided)
        p_val = 1 - cdf(TDist(df), t_stat)
        p_values[model.name] = p_val
    end

    # Apply multiplicity adjustment
    if config.adjust_method == :bonferroni
        for model in config.models
            p_values[model.name] = min(1.0, p_values[model.name] * length(config.models))
        end
    end
    # For Dunnett, p-values would need simulation-based adjustment
    # For simplicity, using critical value comparison instead

    # Identify significant models
    max_t = maximum(values(test_stats))
    significant = Symbol[]
    for model in config.models
        if test_stats[model.name] >= crit_val
            push!(significant, model.name)
        end
    end

    poc = !isempty(significant)

    return MCPStepResult(
        test_stats,
        p_values,
        crit_val,
        contrasts,
        significant,
        max_t,
        poc
    )
end

# ============================================================================
# Model Parameter Estimation
# ============================================================================

"""
    estimate_model_parameters(model::DoseResponseModel, doses::Vector{Float64},
                               means::Vector{Float64}, weights::Vector{Float64};
                               max_iter::Int = 1000, tol::Float64 = 1e-8)

Estimate model parameters using weighted nonlinear least squares.

Uses Gauss-Newton optimization with Levenberg-Marquardt damping.

# Arguments
- `model`: Dose-response model
- `doses`: Dose levels
- `means`: Observed means
- `weights`: Weights (typically 1/variance)
- `max_iter`: Maximum iterations
- `tol`: Convergence tolerance

# Returns
- Dictionary of estimated parameters
"""
function estimate_model_parameters(model::LinearModel, doses::Vector{Float64},
                                   means::Vector{Float64}, weights::Vector{Float64};
                                   max_iter::Int = 1000, tol::Float64 = 1e-8)
    # Weighted linear regression: y = e0 + delta * d
    n = length(doses)
    W = diagm(weights)
    X = hcat(ones(n), doses)

    # WLS: β = (X'WX)^{-1} X'Wy
    beta = (X' * W * X) \ (X' * W * means)

    return Dict{Symbol, Float64}(:e0 => beta[1], :delta => beta[2])
end

function estimate_model_parameters(model::QuadraticModel, doses::Vector{Float64},
                                   means::Vector{Float64}, weights::Vector{Float64};
                                   max_iter::Int = 1000, tol::Float64 = 1e-8)
    # Weighted polynomial regression: y = e0 + b1*d + b2*d²
    n = length(doses)
    W = diagm(weights)
    X = hcat(ones(n), doses, doses.^2)

    beta = (X' * W * X) \ (X' * W * means)

    return Dict{Symbol, Float64}(:e0 => beta[1], :b1 => beta[2], :b2 => beta[3])
end

function estimate_model_parameters(model::EmaxModel, doses::Vector{Float64},
                                   means::Vector{Float64}, weights::Vector{Float64};
                                   max_iter::Int = 1000, tol::Float64 = 1e-8)
    # Emax model: y = e0 + emax * d / (ed50 + d)
    # Use grid search + refinement for ed50, linear for e0, emax

    dose_max = maximum(doses)
    best_sse = Inf
    best_params = Dict{Symbol, Float64}()

    # Grid search for ed50
    ed50_grid = range(dose_max * 0.01, dose_max * 2.0, length = 50)

    for ed50 in ed50_grid
        # For fixed ed50, this is a linear problem in (e0, emax)
        f = doses ./ (ed50 .+ doses)
        X = hcat(ones(length(doses)), f)
        W = diagm(weights)

        try
            beta = (X' * W * X) \ (X' * W * means)
            e0, emax = beta[1], beta[2]

            pred = e0 .+ emax .* f
            sse = sum(weights .* (means .- pred).^2)

            if sse < best_sse
                best_sse = sse
                best_params = Dict(:e0 => e0, :emax => emax, :ed50 => ed50)
            end
        catch
            continue
        end
    end

    # Refinement with Gauss-Newton
    if !isempty(best_params)
        best_params = refine_nonlinear(model, doses, means, weights, best_params;
                                        max_iter = max_iter, tol = tol)
    else
        best_params = Dict(:e0 => mean(means), :emax => 0.0, :ed50 => model.ed50)
    end

    return best_params
end

function estimate_model_parameters(model::SigmoidEmaxModel, doses::Vector{Float64},
                                   means::Vector{Float64}, weights::Vector{Float64};
                                   max_iter::Int = 1000, tol::Float64 = 1e-8)
    # Sigmoid Emax: y = e0 + emax * d^h / (ed50^h + d^h)
    dose_max = maximum(doses)
    best_sse = Inf
    best_params = Dict{Symbol, Float64}()

    # Grid search for ed50 and h
    ed50_grid = range(dose_max * 0.01, dose_max * 2.0, length = 20)
    h_grid = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]

    for ed50 in ed50_grid, h in h_grid
        f = (doses.^h) ./ (ed50^h .+ doses.^h)
        X = hcat(ones(length(doses)), f)
        W = diagm(weights)

        try
            beta = (X' * W * X) \ (X' * W * means)
            e0, emax = beta[1], beta[2]

            pred = e0 .+ emax .* f
            sse = sum(weights .* (means .- pred).^2)

            if sse < best_sse
                best_sse = sse
                best_params = Dict(:e0 => e0, :emax => emax, :ed50 => ed50, :h => h)
            end
        catch
            continue
        end
    end

    if !isempty(best_params)
        best_params = refine_nonlinear(model, doses, means, weights, best_params;
                                        max_iter = max_iter, tol = tol)
    else
        best_params = Dict(:e0 => mean(means), :emax => 0.0, :ed50 => model.ed50, :h => model.h)
    end

    return best_params
end

function estimate_model_parameters(model::ExponentialModel, doses::Vector{Float64},
                                   means::Vector{Float64}, weights::Vector{Float64};
                                   max_iter::Int = 1000, tol::Float64 = 1e-8)
    # Exponential: y = e0 + e1 * (exp(d/delta) - 1)
    dose_max = maximum(doses)
    best_sse = Inf
    best_params = Dict{Symbol, Float64}()

    # Grid search for delta
    delta_grid = range(dose_max * 0.1, dose_max * 5.0, length = 50)

    for delta in delta_grid
        f = exp.(doses ./ delta) .- 1.0
        X = hcat(ones(length(doses)), f)
        W = diagm(weights)

        try
            beta = (X' * W * X) \ (X' * W * means)
            e0, e1 = beta[1], beta[2]

            pred = e0 .+ e1 .* f
            sse = sum(weights .* (means .- pred).^2)

            if sse < best_sse
                best_sse = sse
                best_params = Dict(:e0 => e0, :e1 => e1, :delta => delta)
            end
        catch
            continue
        end
    end

    if isempty(best_params)
        best_params = Dict(:e0 => mean(means), :e1 => 0.0, :delta => model.delta)
    end

    return best_params
end

function estimate_model_parameters(model::BetaModel, doses::Vector{Float64},
                                   means::Vector{Float64}, weights::Vector{Float64};
                                   max_iter::Int = 1000, tol::Float64 = 1e-8)
    # Beta model with fixed shape parameters - estimate e0 and emax
    d1, d2 = model.delta1, model.delta2
    scal = model.scal

    mode = d1 / (d1 + d2)
    B = mode^d1 * (1 - mode)^d2

    f = zeros(length(doses))
    for (i, d) in enumerate(doses)
        x = d / scal
        if 0 < x < 1
            f[i] = (x^d1 * (1 - x)^d2) / B
        end
    end

    X = hcat(ones(length(doses)), f)
    W = diagm(weights)

    try
        beta = (X' * W * X) \ (X' * W * means)
        return Dict(:e0 => beta[1], :emax => beta[2])
    catch
        return Dict(:e0 => mean(means), :emax => 0.0)
    end
end

function estimate_model_parameters(model::LogLinearModel, doses::Vector{Float64},
                                   means::Vector{Float64}, weights::Vector{Float64};
                                   max_iter::Int = 1000, tol::Float64 = 1e-8)
    # Log-linear: y = e0 + delta * log(d + off)
    f = log.(doses .+ model.off)
    X = hcat(ones(length(doses)), f)
    W = diagm(weights)

    try
        beta = (X' * W * X) \ (X' * W * means)
        return Dict(:e0 => beta[1], :delta => beta[2])
    catch
        return Dict(:e0 => mean(means), :delta => 0.0)
    end
end

"""
    refine_nonlinear(model, doses, means, weights, init_params; max_iter, tol)

Refine parameter estimates using Gauss-Newton iteration.
"""
function refine_nonlinear(model::DoseResponseModel, doses::Vector{Float64},
                          means::Vector{Float64}, weights::Vector{Float64},
                          init_params::Dict{Symbol, Float64};
                          max_iter::Int = 100, tol::Float64 = 1e-8)
    params = copy(init_params)
    prev_sse = Inf

    for iter in 1:max_iter
        # Compute predictions
        pred = model_mean_vector(model, doses, params)
        residuals = means .- pred

        # Check convergence
        sse = sum(weights .* residuals.^2)
        if iter > 1 && abs(sse - prev_sse) / max(1e-10, prev_sse) < tol
            break
        end
        prev_sse = sse

        # Numerical Jacobian
        J = compute_jacobian(model, doses, params)

        # Gauss-Newton update with damping
        W = diagm(weights)
        H = J' * W * J
        g = J' * W * residuals

        # Levenberg-Marquardt damping
        lambda = 1e-4 * maximum(diag(H))
        H_damped = H + lambda * I

        try
            delta = H_damped \ g
            update_params!(params, delta, model)
        catch
            break
        end
    end

    return params
end

"""
    compute_jacobian(model, doses, params)

Compute numerical Jacobian of the model.
"""
function compute_jacobian(model::DoseResponseModel, doses::Vector{Float64},
                          params::Dict{Symbol, Float64})
    eps = 1e-7
    n = length(doses)
    param_names = collect(keys(params))
    p = length(param_names)

    J = zeros(n, p)
    f0 = model_mean_vector(model, doses, params)

    for (j, name) in enumerate(param_names)
        params_plus = copy(params)
        params_plus[name] += eps
        f_plus = model_mean_vector(model, doses, params_plus)
        J[:, j] = (f_plus .- f0) ./ eps
    end

    return J
end

"""
    update_params!(params, delta, model)

Update parameter dictionary with Gauss-Newton step.
"""
function update_params!(params::Dict{Symbol, Float64}, delta::Vector{Float64},
                        model::DoseResponseModel)
    param_names = collect(keys(params))
    for (i, name) in enumerate(param_names)
        params[name] += delta[i]
    end

    # Ensure positivity constraints
    for name in [:ed50, :delta, :h]
        if haskey(params, name) && params[name] <= 0
            params[name] = 1e-6
        end
    end
end

# ============================================================================
# Mod Step
# ============================================================================

"""
    mod_step(means::Vector{Float64}, se::Vector{Float64}, n::Vector{Int},
             config::MCPModConfig, significant_models::Vector{Symbol})

Perform the Mod (modeling) step.

Fits significant models to the data, selects the best model, and estimates
the target dose for achieving a clinically relevant effect.

# Arguments
- `means`: Observed mean responses at each dose
- `se`: Standard errors at each dose
- `n`: Sample sizes at each dose
- `config`: MCP-Mod configuration
- `significant_models`: Models that passed the MCP step

# Returns
- `ModStepResult` with fitted models and target dose estimation
"""
function mod_step(means::Vector{Float64}, se::Vector{Float64}, n::Vector{Int},
                  config::MCPModConfig, significant_models::Vector{Symbol})
    doses = config.doses
    K = length(doses)

    # Weights for fitting (inverse variance)
    weights = 1.0 ./ (se.^2 .+ 1e-10)

    # Fit each significant model
    fitted_models = Dict{Symbol, Dict{Symbol, Float64}}()
    aic_values = Dict{Symbol, Float64}()
    bic_values = Dict{Symbol, Float64}()
    gaic_values = Dict{Symbol, Float64}()
    residuals = Dict{Symbol, Vector{Float64}}()
    dose_response = Dict{Symbol, Vector{Float64}}()

    for model in config.models
        if !(model.name in significant_models)
            continue
        end

        # Estimate parameters
        params = estimate_model_parameters(model, doses, means, weights)
        fitted_models[model.name] = params

        # Compute predictions and residuals
        pred = model_mean_vector(model, doses, params)
        res = means .- pred
        residuals[model.name] = res
        dose_response[model.name] = pred

        # Compute model selection criteria
        sse = sum(weights .* res.^2)
        n_total = sum(n)
        p = length(params)

        # Log-likelihood (weighted)
        ll = -0.5 * n_total * log(2π * sse / n_total) - 0.5 * n_total

        # AIC, BIC, gAIC
        aic_values[model.name] = -2 * ll + 2 * p
        bic_values[model.name] = -2 * ll + log(n_total) * p
        gaic_values[model.name] = -2 * ll + 3 * p  # Generalized AIC
    end

    # Select best model (minimum AIC)
    if isempty(aic_values)
        # No significant models - use first model
        model = config.models[1]
        params = estimate_model_parameters(model, doses, means, weights)
        fitted_models[model.name] = params
        pred = model_mean_vector(model, doses, params)
        dose_response[model.name] = pred
        residuals[model.name] = means .- pred
        best_model = model.name
    else
        best_model = argmin(aic_values)
    end

    # Estimate target dose
    target_dose = estimate_target_dose(
        config.models[findfirst(m -> m.name == best_model, config.models)],
        fitted_models[best_model],
        config.delta,
        config.doses,
        config.direction
    )

    return ModStepResult(
        fitted_models,
        best_model,
        aic_values,
        bic_values,
        gaic_values,
        target_dose,
        dose_response,
        residuals
    )
end

"""
    estimate_target_dose(model::DoseResponseModel, params::Dict{Symbol, Float64},
                         delta::Float64, doses::Vector{Float64},
                         direction::Symbol)

Estimate the target dose achieving a clinically relevant effect delta.

# Arguments
- `model`: Fitted dose-response model
- `params`: Fitted parameters
- `delta`: Clinically relevant effect size (relative to placebo)
- `doses`: Available dose levels
- `direction`: Direction of beneficial effect

# Returns
- Estimated target dose
"""
function estimate_target_dose(model::DoseResponseModel, params::Dict{Symbol, Float64},
                              delta::Float64, doses::Vector{Float64},
                              direction::Symbol)
    if delta <= 0
        # If no target delta specified, return dose at 80% of Emax
        e0 = get(params, :e0, 0.0)
        dose_max = maximum(doses)

        # Find max effect
        max_effect = model_mean(model, dose_max, params) - e0
        target_effect = 0.8 * max_effect
        delta = target_effect
    end

    e0 = get(params, :e0, 0.0)
    dose_max = maximum(doses)

    # Binary search for target dose
    target = e0 + (direction == :increasing ? delta : -delta)

    lo, hi = 0.0, dose_max * 2.0
    for _ in 1:50
        mid = (lo + hi) / 2
        effect = model_mean(model, mid, params)

        if direction == :increasing
            if effect < target
                lo = mid
            else
                hi = mid
            end
        else
            if effect > target
                lo = mid
            else
                hi = mid
            end
        end

        if abs(hi - lo) < dose_max * 1e-6
            break
        end
    end

    target_dose = (lo + hi) / 2

    # Clamp to dose range
    target_dose = clamp(target_dose, minimum(doses), maximum(doses))

    return target_dose
end

# ============================================================================
# Main MCP-Mod Function
# ============================================================================

"""
    run_mcpmod(data::Dict{Float64, Vector{Float64}}, config::MCPModConfig)

Run complete MCP-Mod analysis on dose-response data.

# Arguments
- `data`: Dictionary mapping doses to vectors of observations
- `config`: MCP-Mod configuration

# Returns
- `MCPModResult` with complete analysis results

# Example
```julia
# Define candidate models
models = [
    LinearModel(),
    EmaxModel(ed50 = 25.0),
    SigmoidEmaxModel(ed50 = 25.0, h = 2.0),
    ExponentialModel(delta = 50.0)
]

# Configure analysis
config = MCPModConfig(
    doses = [0.0, 10.0, 25.0, 50.0, 100.0],
    models = models,
    alpha = 0.025,
    delta = 0.4,  # Clinically relevant effect
    direction = :increasing
)

# Data: dose => observations
data = Dict(
    0.0 => [0.1, -0.2, 0.3, ...],
    10.0 => [0.5, 0.4, 0.6, ...],
    25.0 => [1.0, 0.8, 1.2, ...],
    50.0 => [1.5, 1.3, 1.6, ...],
    100.0 => [1.8, 1.7, 1.9, ...]
)

result = run_mcpmod(data, config)
```
"""
function run_mcpmod(data::Dict{Float64, Vector{Float64}}, config::MCPModConfig)
    # Extract summary statistics
    doses = config.doses
    K = length(doses)

    means = zeros(K)
    se = zeros(K)
    n = zeros(Int, K)

    for (i, dose) in enumerate(doses)
        obs = get(data, dose, Float64[])
        @assert !isempty(obs) "No observations for dose $dose"

        n[i] = length(obs)
        means[i] = mean(obs)
        se[i] = std(obs) / sqrt(n[i])
    end

    return run_mcpmod(means, se, n, config)
end

"""
    run_mcpmod(means::Vector{Float64}, se::Vector{Float64}, n::Vector{Int},
               config::MCPModConfig)

Run MCP-Mod analysis from summary statistics.

# Arguments
- `means`: Mean responses at each dose
- `se`: Standard errors at each dose
- `n`: Sample sizes at each dose
- `config`: MCP-Mod configuration

# Returns
- `MCPModResult`
"""
function run_mcpmod(means::Vector{Float64}, se::Vector{Float64}, n::Vector{Int},
                    config::MCPModConfig)
    # Step 1: MCP (Multiple Comparison Procedure)
    mcp_result = mcp_step(means, se, n, config)

    # Step 2: Mod (Modeling) - only if PoC established
    mod_result = nothing
    if mcp_result.proof_of_concept
        mod_result = mod_step(means, se, n, config, mcp_result.significant_models)
    end

    return MCPModResult(config, mcp_result, mod_result, means, se, n)
end

# ============================================================================
# Convenience Functions
# ============================================================================

"""
    create_mcpmod_models(doses::Vector{Float64};
                          include_linear::Bool = true,
                          include_quadratic::Bool = true,
                          include_emax::Bool = true,
                          include_sigemax::Bool = true,
                          include_exponential::Bool = true,
                          include_beta::Bool = false,
                          include_loglinear::Bool = false)

Create a set of candidate dose-response models.

Automatically sets reasonable parameter values based on the dose range.

# Arguments
- `doses`: Dose levels
- `include_*`: Whether to include each model type

# Returns
- Vector of `DoseResponseModel`
"""
function create_mcpmod_models(doses::Vector{Float64};
                               include_linear::Bool = true,
                               include_quadratic::Bool = true,
                               include_emax::Bool = true,
                               include_sigemax::Bool = true,
                               include_exponential::Bool = true,
                               include_beta::Bool = false,
                               include_loglinear::Bool = false)
    dose_max = maximum(doses)
    dose_mid = dose_max / 2

    models = DoseResponseModel[]

    if include_linear
        push!(models, LinearModel())
    end

    if include_quadratic
        push!(models, QuadraticModel(delta = dose_max))
    end

    if include_emax
        # Multiple ED50 guesses
        push!(models, EmaxModel(ed50 = dose_mid * 0.5))
        push!(models, EmaxModel(ed50 = dose_mid))
    end

    if include_sigemax
        push!(models, SigmoidEmaxModel(ed50 = dose_mid, h = 2.0))
        push!(models, SigmoidEmaxModel(ed50 = dose_mid, h = 4.0))
    end

    if include_exponential
        push!(models, ExponentialModel(delta = dose_max))
    end

    if include_beta
        push!(models, BetaModel(delta1 = 1.0, delta2 = 1.0, scal = dose_max * 1.2))
    end

    if include_loglinear
        push!(models, LogLinearModel(off = dose_max * 0.1))
    end

    return models
end

"""
    default_mcpmod_config(doses::Vector{Float64}; kwargs...)

Create a default MCP-Mod configuration with standard candidate models.

# Arguments
- `doses`: Dose levels (must include 0 for placebo)
- `kwargs`: Additional configuration options

# Returns
- `MCPModConfig`
"""
function default_mcpmod_config(doses::Vector{Float64};
                                alpha::Float64 = 0.025,
                                delta::Float64 = 0.0,
                                direction::Symbol = :increasing)
    models = create_mcpmod_models(doses)

    return MCPModConfig(
        doses = doses,
        models = models,
        alpha = alpha,
        delta = delta,
        direction = direction
    )
end
