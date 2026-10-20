# Exposure-Response Analysis Module
# Links PK exposure metrics to clinical outcomes for dose selection
# Integrates with NCA module for exposure calculations

export ExposureResponseConfig, ExposureResponseResult
export ExposureMetric, CmaxMetric, AUCMetric, CtroughMetric, CavgMetric, CustomExposureMetric
export ResponseType, BinaryResponse, ContinuousResponse, TimeToEventResponse, CountResponse
export ERModel, LogisticERModel, LinearERModel, EmaxERModel, LogLinearERModel
export run_exposure_response, fit_er_model, predict_response
export compute_exposure_quantiles, bootstrap_er_model, plot_er_data
export predict_response_at_dose, find_target_exposure, extract_exposure

using Statistics
using LinearAlgebra
using Random
using Distributions

# ============================================================================
# Exposure Metric Types
# ============================================================================

"""
Abstract type for exposure metrics.
"""
abstract type ExposureMetric end

"""
    CmaxMetric

Use Cmax (maximum concentration) as exposure metric.
"""
struct CmaxMetric <: ExposureMetric
    log_transform::Bool
    CmaxMetric(; log_transform::Bool = false) = new(log_transform)
end

"""
    AUCMetric

Use AUC (area under the curve) as exposure metric.

# Fields
- `type::Symbol`: `:auc_0_inf`, `:auc_0_t`, `:auc_0_tau`
- `log_transform::Bool`: Whether to log-transform the metric
"""
struct AUCMetric <: ExposureMetric
    type::Symbol
    log_transform::Bool

    function AUCMetric(; type::Symbol = :auc_0_inf, log_transform::Bool = false)
        @assert type in [:auc_0_inf, :auc_0_t, :auc_0_tau] "Invalid AUC type"
        new(type, log_transform)
    end
end

"""
    CtroughMetric

Use Ctrough (trough concentration) as exposure metric.
"""
struct CtroughMetric <: ExposureMetric
    log_transform::Bool
    CtroughMetric(; log_transform::Bool = false) = new(log_transform)
end

"""
    CavgMetric

Use Cavg (average concentration) as exposure metric.
"""
struct CavgMetric <: ExposureMetric
    log_transform::Bool
    CavgMetric(; log_transform::Bool = false) = new(log_transform)
end

"""
    CustomExposureMetric

Custom exposure metric with user-defined extraction function.

# Fields
- `name::Symbol`: Metric name
- `extractor::Function`: Function to extract metric from NCAResult
- `log_transform::Bool`: Whether to log-transform
"""
struct CustomExposureMetric <: ExposureMetric
    name::Symbol
    extractor::Function
    log_transform::Bool

    function CustomExposureMetric(name::Symbol, extractor::Function;
                                   log_transform::Bool = false)
        new(name, extractor, log_transform)
    end
end

# ============================================================================
# Response Type Definitions
# ============================================================================

"""
Abstract type for response outcomes.
"""
abstract type ResponseType end

"""
    BinaryResponse

Binary (0/1) response outcome (e.g., responder/non-responder).
"""
struct BinaryResponse <: ResponseType
    success_label::Symbol
    BinaryResponse(; success_label::Symbol = :responder) = new(success_label)
end

"""
    ContinuousResponse

Continuous response outcome (e.g., change in biomarker).
"""
struct ContinuousResponse <: ResponseType
    baseline_corrected::Bool
    ContinuousResponse(; baseline_corrected::Bool = true) = new(baseline_corrected)
end

"""
    TimeToEventResponse

Time-to-event response outcome (e.g., survival, progression).
"""
struct TimeToEventResponse <: ResponseType
    censoring_symbol::Symbol
    TimeToEventResponse(; censoring_symbol::Symbol = :censored) = new(censoring_symbol)
end

"""
    CountResponse

Count response outcome (e.g., number of events).
"""
struct CountResponse <: ResponseType
    offset_variable::Union{Symbol, Nothing}
    CountResponse(; offset_variable::Union{Symbol, Nothing} = nothing) = new(offset_variable)
end

# ============================================================================
# Exposure-Response Models
# ============================================================================

"""
Abstract type for exposure-response models.
"""
abstract type ERModel end

"""
    LogisticERModel

Logistic regression model for binary outcomes.
P(response) = 1 / (1 + exp(-(β0 + β1*exposure)))

# Fields
- `intercept::Union{Float64, Nothing}`: Fitted intercept
- `slope::Union{Float64, Nothing}`: Fitted slope
- `covariates::Vector{Symbol}`: Additional covariates
"""
mutable struct LogisticERModel <: ERModel
    intercept::Union{Float64, Nothing}
    slope::Union{Float64, Nothing}
    covariates::Vector{Symbol}
    covariate_coefficients::Dict{Symbol, Float64}

    function LogisticERModel(; covariates::Vector{Symbol} = Symbol[])
        new(nothing, nothing, covariates, Dict{Symbol, Float64}())
    end
end

"""
    LinearERModel

Linear regression model for continuous outcomes.
E[Y] = β0 + β1*exposure

# Fields
- `intercept::Union{Float64, Nothing}`: Fitted intercept
- `slope::Union{Float64, Nothing}`: Fitted slope
- `residual_variance::Union{Float64, Nothing}`: Residual variance
- `covariates::Vector{Symbol}`: Additional covariates
"""
mutable struct LinearERModel <: ERModel
    intercept::Union{Float64, Nothing}
    slope::Union{Float64, Nothing}
    residual_variance::Union{Float64, Nothing}
    covariates::Vector{Symbol}
    covariate_coefficients::Dict{Symbol, Float64}

    function LinearERModel(; covariates::Vector{Symbol} = Symbol[])
        new(nothing, nothing, nothing, covariates, Dict{Symbol, Float64}())
    end
end

"""
    EmaxERModel

Emax model for exposure-response relationships.
E[Y] = E0 + Emax * exposure / (EC50 + exposure)

# Fields
- `e0::Union{Float64, Nothing}`: Baseline effect
- `emax::Union{Float64, Nothing}`: Maximum effect
- `ec50::Union{Float64, Nothing}`: Exposure at 50% max effect
"""
mutable struct EmaxERModel <: ERModel
    e0::Union{Float64, Nothing}
    emax::Union{Float64, Nothing}
    ec50::Union{Float64, Nothing}

    function EmaxERModel()
        new(nothing, nothing, nothing)
    end
end

"""
    LogLinearERModel

Log-linear model for exposure-response.
E[Y] = β0 + β1*log(exposure)
"""
mutable struct LogLinearERModel <: ERModel
    intercept::Union{Float64, Nothing}
    slope::Union{Float64, Nothing}

    function LogLinearERModel()
        new(nothing, nothing)
    end
end

# ============================================================================
# Configuration and Result Types
# ============================================================================

"""
    ExposureResponseConfig

Configuration for exposure-response analysis.

# Fields
- `exposure_metric::ExposureMetric`: Metric for exposure
- `response_type::ResponseType`: Type of outcome
- `model::ERModel`: E-R model type
- `confidence_level::Float64`: CI level (default 0.95)
- `n_bootstrap::Int`: Number of bootstrap samples for CI
- `seed::Int`: Random seed
"""
struct ExposureResponseConfig
    exposure_metric::ExposureMetric
    response_type::ResponseType
    model::ERModel
    confidence_level::Float64
    n_bootstrap::Int
    seed::Int

    function ExposureResponseConfig(;
            exposure_metric::ExposureMetric = AUCMetric(),
            response_type::ResponseType = BinaryResponse(),
            model::ERModel = LogisticERModel(),
            confidence_level::Float64 = 0.95,
            n_bootstrap::Int = 1000,
            seed::Int = 12345)
        @assert 0 < confidence_level < 1 "confidence_level must be between 0 and 1"
        new(exposure_metric, response_type, model, confidence_level, n_bootstrap, seed)
    end
end

"""
    ExposureResponseResult

Results from exposure-response analysis.

# Fields
- `config::ExposureResponseConfig`: Configuration used
- `fitted_model::ERModel`: Fitted model with estimated parameters
- `n_subjects::Int`: Number of subjects
- `exposure_values::Vector{Float64}`: Exposure values used
- `response_values::Vector{Float64}`: Response values (0/1 for binary)
- `predicted_responses::Vector{Float64}`: Model predictions
- `parameter_estimates::Dict{Symbol, Float64}`: Parameter estimates
- `parameter_se::Dict{Symbol, Float64}`: Standard errors
- `parameter_ci::Dict{Symbol, Tuple{Float64, Float64}}`: Confidence intervals
- `odds_ratios::Dict{Symbol, Float64}`: Odds ratios (logistic only)
- `model_fit::Dict{Symbol, Float64}`: Goodness of fit metrics
- `exposure_quantiles::Dict{Float64, Float64}`: Exposure quantiles
- `response_at_quantiles::Dict{Float64, Float64}`: Predicted response at quantiles
"""
struct ExposureResponseResult
    config::ExposureResponseConfig
    fitted_model::ERModel
    n_subjects::Int
    exposure_values::Vector{Float64}
    response_values::Vector{Float64}
    predicted_responses::Vector{Float64}
    parameter_estimates::Dict{Symbol, Float64}
    parameter_se::Dict{Symbol, Float64}
    parameter_ci::Dict{Symbol, Tuple{Float64, Float64}}
    odds_ratios::Dict{Symbol, Float64}
    model_fit::Dict{Symbol, Float64}
    exposure_quantiles::Dict{Float64, Float64}
    response_at_quantiles::Dict{Float64, Float64}
end

# ============================================================================
# Exposure Extraction Functions
# ============================================================================

"""
    extract_exposure(nca_result::NCAResult, metric::ExposureMetric)

Extract exposure value from NCA result.
"""
function extract_exposure(nca_result, metric::CmaxMetric)
    val = nca_result.cmax
    return metric.log_transform ? log(val) : val
end

function extract_exposure(nca_result, metric::AUCMetric)
    val = if metric.type == :auc_0_inf
        nca_result.auc_0_inf
    elseif metric.type == :auc_0_t
        nca_result.auc_0_t
    elseif metric.type == :auc_0_tau
        nca_result.auc_0_tau
    end

    if val === nothing
        return NaN
    end

    return metric.log_transform ? log(val) : val
end

function extract_exposure(nca_result, metric::CtroughMetric)
    val = nca_result.cmin
    if val === nothing
        return NaN
    end
    return metric.log_transform ? log(val) : val
end

function extract_exposure(nca_result, metric::CavgMetric)
    val = nca_result.cavg
    if val === nothing
        return NaN
    end
    return metric.log_transform ? log(val) : val
end

function extract_exposure(nca_result, metric::CustomExposureMetric)
    val = metric.extractor(nca_result)
    return metric.log_transform ? log(val) : val
end

# ============================================================================
# Model Fitting Functions
# ============================================================================

"""
    fit_logistic(x::Vector{Float64}, y::Vector{Float64};
                 max_iter::Int = 100, tol::Float64 = 1e-8)

Fit logistic regression using iteratively reweighted least squares (IRLS).

Returns (intercept, slope, converged)
"""
function fit_logistic(x::Vector{Float64}, y::Vector{Float64};
                      max_iter::Int = 100, tol::Float64 = 1e-8)
    n = length(x)
    @assert length(y) == n "x and y must have same length"
    @assert all(yi -> yi == 0 || yi == 1, y) "y must be binary (0 or 1)"

    # Initialize with linear approximation
    p_init = clamp.(y .+ 0.5, 0.25, 0.75)
    logit_p = log.(p_init ./ (1 .- p_init))

    X = hcat(ones(n), x)
    beta = X \ logit_p

    # IRLS iteration
    for iter in 1:max_iter
        # Predicted probabilities
        eta = X * beta
        mu = 1.0 ./ (1.0 .+ exp.(-eta))
        mu = clamp.(mu, 1e-10, 1 - 1e-10)

        # Weights
        W = diagm(mu .* (1 .- mu))

        # Working response
        z = eta .+ (y .- mu) ./ (mu .* (1 .- mu))

        # Update
        try
            beta_new = (X' * W * X) \ (X' * W * z)
            delta = maximum(abs.(beta_new .- beta))
            beta = beta_new

            if delta < tol
                return beta[1], beta[2], true
            end
        catch
            return beta[1], beta[2], false
        end
    end

    return beta[1], beta[2], false
end

"""
    fit_linear(x::Vector{Float64}, y::Vector{Float64})

Fit linear regression using ordinary least squares.

Returns (intercept, slope, residual_variance)
"""
function fit_linear(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    X = hcat(ones(n), x)

    beta = (X' * X) \ (X' * y)
    residuals = y .- X * beta
    residual_var = sum(residuals.^2) / (n - 2)

    return beta[1], beta[2], residual_var
end

"""
    fit_emax(x::Vector{Float64}, y::Vector{Float64};
             max_iter::Int = 100, tol::Float64 = 1e-8)

Fit Emax model using Gauss-Newton optimization.
E[Y] = E0 + Emax * x / (EC50 + x)

Returns (e0, emax, ec50)
"""
function fit_emax(x::Vector{Float64}, y::Vector{Float64};
                  max_iter::Int = 100, tol::Float64 = 1e-8)
    n = length(x)

    # Grid search for initial EC50
    x_max = maximum(x)
    best_sse = Inf
    best_params = (mean(y), 0.0, x_max / 2)

    ec50_grid = range(0.1 * x_max, 2.0 * x_max, length = 20)

    for ec50 in ec50_grid
        f = x ./ (ec50 .+ x)
        X = hcat(ones(n), f)

        try
            beta = (X' * X) \ (X' * y)
            pred = X * beta
            sse = sum((y .- pred).^2)

            if sse < best_sse
                best_sse = sse
                best_params = (beta[1], beta[2], ec50)
            end
        catch
            continue
        end
    end

    e0, emax, ec50 = best_params

    # Gauss-Newton refinement
    params = [e0, emax, ec50]
    prev_sse = Inf

    for iter in 1:max_iter
        pred = params[1] .+ params[2] .* x ./ (params[3] .+ x)
        residuals = y .- pred

        sse = sum(residuals.^2)
        if iter > 1 && abs(sse - prev_sse) / max(1e-10, prev_sse) < tol
            break
        end
        prev_sse = sse

        # Jacobian
        J = zeros(n, 3)
        J[:, 1] .= 1.0  # dE/dE0
        J[:, 2] = x ./ (params[3] .+ x)  # dE/dEmax
        J[:, 3] = -params[2] .* x ./ (params[3] .+ x).^2  # dE/dEC50

        # Update
        try
            delta = (J' * J + 1e-6 * I) \ (J' * residuals)
            params .+= delta
            params[3] = max(1e-6, params[3])  # Ensure EC50 > 0
        catch
            break
        end
    end

    return params[1], params[2], params[3]
end

# ============================================================================
# Model Fitting Wrapper
# ============================================================================

"""
    fit_er_model(exposure::Vector{Float64}, response::Vector{Float64},
                 model::ERModel; covariates::Matrix{Float64} = zeros(0, 0))

Fit exposure-response model to data.
"""
function fit_er_model(exposure::Vector{Float64}, response::Vector{Float64},
                      model::LogisticERModel; covariates::Matrix{Float64} = zeros(0, 0))
    intercept, slope, converged = fit_logistic(exposure, response)

    model.intercept = intercept
    model.slope = slope

    return model
end

function fit_er_model(exposure::Vector{Float64}, response::Vector{Float64},
                      model::LinearERModel; covariates::Matrix{Float64} = zeros(0, 0))
    intercept, slope, residual_var = fit_linear(exposure, response)

    model.intercept = intercept
    model.slope = slope
    model.residual_variance = residual_var

    return model
end

function fit_er_model(exposure::Vector{Float64}, response::Vector{Float64},
                      model::EmaxERModel; covariates::Matrix{Float64} = zeros(0, 0))
    e0, emax, ec50 = fit_emax(exposure, response)

    model.e0 = e0
    model.emax = emax
    model.ec50 = ec50

    return model
end

function fit_er_model(exposure::Vector{Float64}, response::Vector{Float64},
                      model::LogLinearERModel; covariates::Matrix{Float64} = zeros(0, 0))
    # Log-transform exposure
    log_exp = log.(max.(exposure, 1e-10))
    intercept, slope, _ = fit_linear(log_exp, response)

    model.intercept = intercept
    model.slope = slope

    return model
end

# ============================================================================
# Prediction Functions
# ============================================================================

"""
    predict_response(model::ERModel, exposure::Vector{Float64})

Predict response from fitted model.
"""
function predict_response(model::LogisticERModel, exposure::Vector{Float64})
    @assert model.intercept !== nothing "Model not fitted"
    eta = model.intercept .+ model.slope .* exposure
    return 1.0 ./ (1.0 .+ exp.(-eta))
end

function predict_response(model::LinearERModel, exposure::Vector{Float64})
    @assert model.intercept !== nothing "Model not fitted"
    return model.intercept .+ model.slope .* exposure
end

function predict_response(model::EmaxERModel, exposure::Vector{Float64})
    @assert model.e0 !== nothing "Model not fitted"
    return model.e0 .+ model.emax .* exposure ./ (model.ec50 .+ exposure)
end

function predict_response(model::LogLinearERModel, exposure::Vector{Float64})
    @assert model.intercept !== nothing "Model not fitted"
    log_exp = log.(max.(exposure, 1e-10))
    return model.intercept .+ model.slope .* log_exp
end

# ============================================================================
# Bootstrap for Confidence Intervals
# ============================================================================

"""
    bootstrap_er_model(exposure::Vector{Float64}, response::Vector{Float64},
                       model_type::Type{<:ERModel}, n_bootstrap::Int;
                       seed::Int = 12345, alpha::Float64 = 0.05)

Bootstrap exposure-response model for parameter uncertainty.

Returns dictionary of parameter CIs.
"""
function bootstrap_er_model(exposure::Vector{Float64}, response::Vector{Float64},
                            model_type::Type{<:ERModel}, n_bootstrap::Int;
                            seed::Int = 12345, alpha::Float64 = 0.05)
    rng = MersenneTwister(seed)
    n = length(exposure)

    # Storage for bootstrap estimates
    params_list = []

    for b in 1:n_bootstrap
        # Bootstrap sample
        idx = rand(rng, 1:n, n)
        exp_b = exposure[idx]
        resp_b = response[idx]

        # Fit model
        try
            model = model_type()
            fitted = fit_er_model(exp_b, resp_b, model)

            if model_type == LogisticERModel
                push!(params_list, (intercept = fitted.intercept, slope = fitted.slope))
            elseif model_type == LinearERModel
                push!(params_list, (intercept = fitted.intercept, slope = fitted.slope))
            elseif model_type == EmaxERModel
                push!(params_list, (e0 = fitted.e0, emax = fitted.emax, ec50 = fitted.ec50))
            elseif model_type == LogLinearERModel
                push!(params_list, (intercept = fitted.intercept, slope = fitted.slope))
            end
        catch
            continue
        end
    end

    if isempty(params_list)
        return Dict{Symbol, Tuple{Float64, Float64}}()
    end

    # Compute CIs
    ci_dict = Dict{Symbol, Tuple{Float64, Float64}}()
    param_names = keys(params_list[1])

    for param in param_names
        values = [p[param] for p in params_list if !isnothing(p[param]) && !isnan(p[param])]
        if length(values) >= 10
            lower = quantile(values, alpha / 2)
            upper = quantile(values, 1 - alpha / 2)
            ci_dict[param] = (lower, upper)
        end
    end

    return ci_dict
end

# ============================================================================
# Main Analysis Function
# ============================================================================

"""
    run_exposure_response(exposure::Vector{Float64}, response::Vector{Float64},
                          config::ExposureResponseConfig)

Run complete exposure-response analysis.

# Arguments
- `exposure`: Exposure values
- `response`: Response values (0/1 for binary, continuous for others)
- `config`: E-R configuration

# Returns
- `ExposureResponseResult`: Complete analysis results

# Example
```julia
# Binary response example
exposure = [10.0, 20.0, 30.0, 40.0, 50.0, ...]
response = [0, 0, 1, 0, 1, ...]  # 0=non-responder, 1=responder

config = ExposureResponseConfig(
    exposure_metric = AUCMetric(type = :auc_0_inf),
    response_type = BinaryResponse(),
    model = LogisticERModel()
)

result = run_exposure_response(exposure, response, config)
```
"""
function run_exposure_response(exposure::Vector{Float64}, response::Vector{Float64},
                               config::ExposureResponseConfig)
    n = length(exposure)
    @assert length(response) == n "exposure and response must have same length"

    # Remove missing values
    valid_idx = findall(i -> !isnan(exposure[i]) && !isnan(response[i]), 1:n)
    exp_clean = exposure[valid_idx]
    resp_clean = response[valid_idx]
    n_valid = length(exp_clean)

    # Fit model
    fitted_model = fit_er_model(exp_clean, resp_clean, config.model)

    # Predictions
    predicted = predict_response(fitted_model, exp_clean)

    # Parameter estimates
    param_estimates = Dict{Symbol, Float64}()
    param_se = Dict{Symbol, Float64}()

    if config.model isa LogisticERModel
        param_estimates[:intercept] = fitted_model.intercept
        param_estimates[:slope] = fitted_model.slope

        # Compute standard errors using Fisher information
        eta = fitted_model.intercept .+ fitted_model.slope .* exp_clean
        mu = 1.0 ./ (1.0 .+ exp.(-eta))
        mu = clamp.(mu, 1e-10, 1 - 1e-10)
        W = diagm(mu .* (1 .- mu))
        X = hcat(ones(n_valid), exp_clean)
        info = X' * W * X
        try
            cov = inv(info)
            param_se[:intercept] = sqrt(cov[1, 1])
            param_se[:slope] = sqrt(cov[2, 2])
        catch
            param_se[:intercept] = NaN
            param_se[:slope] = NaN
        end

    elseif config.model isa LinearERModel
        param_estimates[:intercept] = fitted_model.intercept
        param_estimates[:slope] = fitted_model.slope
        param_estimates[:residual_variance] = fitted_model.residual_variance

        # Standard errors
        X = hcat(ones(n_valid), exp_clean)
        residual_var = fitted_model.residual_variance
        try
            cov = residual_var * inv(X' * X)
            param_se[:intercept] = sqrt(cov[1, 1])
            param_se[:slope] = sqrt(cov[2, 2])
        catch
            param_se[:intercept] = NaN
            param_se[:slope] = NaN
        end

    elseif config.model isa EmaxERModel
        param_estimates[:e0] = fitted_model.e0
        param_estimates[:emax] = fitted_model.emax
        param_estimates[:ec50] = fitted_model.ec50

    elseif config.model isa LogLinearERModel
        param_estimates[:intercept] = fitted_model.intercept
        param_estimates[:slope] = fitted_model.slope
    end

    # Bootstrap CIs
    model_type = typeof(config.model)
    param_ci = bootstrap_er_model(exp_clean, resp_clean, model_type,
                                   config.n_bootstrap; seed = config.seed,
                                   alpha = 1 - config.confidence_level)

    # Odds ratios (logistic only)
    odds_ratios = Dict{Symbol, Float64}()
    if config.model isa LogisticERModel
        # OR for 1 unit increase in exposure
        odds_ratios[:per_unit] = exp(fitted_model.slope)

        # OR for interquartile range
        iqr = quantile(exp_clean, 0.75) - quantile(exp_clean, 0.25)
        odds_ratios[:per_iqr] = exp(fitted_model.slope * iqr)
    end

    # Model fit statistics
    model_fit = Dict{Symbol, Float64}()

    if config.model isa LogisticERModel || config.response_type isa BinaryResponse
        # Log-likelihood
        mu = clamp.(predicted, 1e-10, 1 - 1e-10)
        ll = sum(resp_clean .* log.(mu) .+ (1 .- resp_clean) .* log.(1 .- mu))
        model_fit[:log_likelihood] = ll

        # Null log-likelihood
        p0 = mean(resp_clean)
        ll_null = n_valid * (p0 * log(p0) + (1 - p0) * log(1 - p0))
        model_fit[:null_log_likelihood] = ll_null

        # McFadden's R²
        model_fit[:mcfadden_r2] = 1 - ll / ll_null

        # AIC, BIC
        k = 2  # Number of parameters
        model_fit[:aic] = -2 * ll + 2 * k
        model_fit[:bic] = -2 * ll + log(n_valid) * k

        # Hosmer-Lemeshow (simplified)
        # Brier score
        model_fit[:brier_score] = mean((predicted .- resp_clean).^2)

    elseif config.model isa LinearERModel
        # R²
        ss_res = sum((resp_clean .- predicted).^2)
        ss_tot = sum((resp_clean .- mean(resp_clean)).^2)
        model_fit[:r_squared] = 1 - ss_res / ss_tot

        # Adjusted R²
        model_fit[:adjusted_r_squared] = 1 - (1 - model_fit[:r_squared]) * (n_valid - 1) / (n_valid - 2)

        # RMSE
        model_fit[:rmse] = sqrt(mean((predicted .- resp_clean).^2))

        # AIC, BIC
        ll = -n_valid / 2 * log(2π * fitted_model.residual_variance) -
             ss_res / (2 * fitted_model.residual_variance)
        k = 3  # intercept, slope, variance
        model_fit[:aic] = -2 * ll + 2 * k
        model_fit[:bic] = -2 * ll + log(n_valid) * k
    end

    # Exposure quantiles
    quantile_levels = [0.1, 0.25, 0.5, 0.75, 0.9]
    exposure_quantiles = Dict{Float64, Float64}()
    response_at_quantiles = Dict{Float64, Float64}()

    for q in quantile_levels
        exp_q = quantile(exp_clean, q)
        exposure_quantiles[q] = exp_q
        response_at_quantiles[q] = predict_response(fitted_model, [exp_q])[1]
    end

    return ExposureResponseResult(
        config,
        fitted_model,
        n_valid,
        exp_clean,
        resp_clean,
        predicted,
        param_estimates,
        param_se,
        param_ci,
        odds_ratios,
        model_fit,
        exposure_quantiles,
        response_at_quantiles
    )
end

"""
    run_exposure_response_from_nca(nca_results::Vector, responses::Vector{Float64},
                                    config::ExposureResponseConfig)

Run exposure-response analysis using NCA results.

# Arguments
- `nca_results`: Vector of NCAResult objects
- `responses`: Response values for each subject
- `config`: E-R configuration

# Returns
- `ExposureResponseResult`
"""
function run_exposure_response_from_nca(nca_results::Vector, responses::Vector{Float64},
                                         config::ExposureResponseConfig)
    n = length(nca_results)
    @assert length(responses) == n "nca_results and responses must have same length"

    # Extract exposure from NCA results
    exposure = Float64[]
    for nca in nca_results
        try
            exp_val = extract_exposure(nca, config.exposure_metric)
            push!(exposure, exp_val)
        catch
            push!(exposure, NaN)
        end
    end

    return run_exposure_response(exposure, responses, config)
end

# ============================================================================
# Exposure Quantile Analysis
# ============================================================================

"""
    compute_exposure_quantiles(exposure::Vector{Float64}, response::Vector{Float64};
                                n_groups::Int = 4)

Compute response rates by exposure quartile (or other grouping).

# Arguments
- `exposure`: Exposure values
- `response`: Response values (binary)
- `n_groups`: Number of exposure groups (default: 4 for quartiles)

# Returns
- Dictionary with group bounds and response rates
"""
function compute_exposure_quantiles(exposure::Vector{Float64}, response::Vector{Float64};
                                     n_groups::Int = 4)
    n = length(exposure)

    # Remove missing
    valid_idx = findall(i -> !isnan(exposure[i]) && !isnan(response[i]), 1:n)
    exp_clean = exposure[valid_idx]
    resp_clean = response[valid_idx]

    # Compute quantile boundaries
    quantile_probs = range(0, 1, length = n_groups + 1)
    boundaries = quantile(exp_clean, quantile_probs)

    results = Dict{Symbol, Any}()
    results[:n_groups] = n_groups
    results[:boundaries] = boundaries
    results[:group_stats] = []

    for g in 1:n_groups
        lower = boundaries[g]
        upper = boundaries[g + 1]

        # Include lower boundary, exclude upper (except for last group)
        if g == n_groups
            in_group = findall(x -> lower <= x <= upper, exp_clean)
        else
            in_group = findall(x -> lower <= x < upper, exp_clean)
        end

        n_in_group = length(in_group)
        mean_exposure = n_in_group > 0 ? mean(exp_clean[in_group]) : NaN
        response_rate = n_in_group > 0 ? mean(resp_clean[in_group]) : NaN
        n_responders = n_in_group > 0 ? sum(resp_clean[in_group]) : 0

        push!(results[:group_stats], (
            group = g,
            lower = lower,
            upper = upper,
            n = n_in_group,
            mean_exposure = mean_exposure,
            response_rate = response_rate,
            n_responders = n_responders
        ))
    end

    return results
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
    predict_response_at_dose(er_result::ExposureResponseResult,
                             exposure_value::Float64)

Predict response probability/value at a specific exposure level.
"""
function predict_response_at_dose(er_result::ExposureResponseResult,
                                  exposure_value::Float64)
    return predict_response(er_result.fitted_model, [exposure_value])[1]
end

"""
    find_target_exposure(er_result::ExposureResponseResult,
                         target_response::Float64)

Find exposure level needed to achieve target response (inverse prediction).

For logistic models, finds exposure where P(response) = target_response.
"""
function find_target_exposure(er_result::ExposureResponseResult,
                              target_response::Float64)
    model = er_result.fitted_model

    if model isa LogisticERModel
        @assert 0 < target_response < 1 "Target response must be between 0 and 1"

        # Inverse logistic: exposure = (logit(p) - intercept) / slope
        logit_p = log(target_response / (1 - target_response))
        return (logit_p - model.intercept) / model.slope

    elseif model isa LinearERModel
        # Inverse linear: exposure = (target - intercept) / slope
        return (target_response - model.intercept) / model.slope

    elseif model isa EmaxERModel
        # Inverse Emax: exposure = EC50 * (target - E0) / (Emax - target + E0)
        @assert target_response > model.e0 "Target must be above E0"
        @assert target_response < model.e0 + model.emax "Target must be below E0 + Emax"

        return model.ec50 * (target_response - model.e0) / (model.emax - (target_response - model.e0))

    else
        error("Inverse prediction not implemented for $(typeof(model))")
    end
end
