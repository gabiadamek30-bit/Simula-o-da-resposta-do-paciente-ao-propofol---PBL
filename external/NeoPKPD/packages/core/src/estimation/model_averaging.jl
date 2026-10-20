# Model Averaging Module
# AIC/BIC weighted model averaging for predictions and parameter estimates

using Statistics
using LinearAlgebra

export ModelAveragingSpec, ModelAveragingResult, CandidateModel
export compute_model_weights, average_predictions, average_parameters
export AICWeighting, BICWeighting, StackedWeighting
export model_selection_table, best_model

# ============================================================================
# Weighting Methods
# ============================================================================

"""
Abstract type for model weighting methods.
"""
abstract type WeightingMethod end

"""
Akaike Information Criterion (AIC) weighting.
w_i = exp(-0.5 * Δ_i) / Σ exp(-0.5 * Δ_j)
where Δ_i = AIC_i - AIC_min
"""
struct AICWeighting <: WeightingMethod end

"""
Bayesian Information Criterion (BIC) weighting.
w_i = exp(-0.5 * Δ_i) / Σ exp(-0.5 * Δ_j)
where Δ_i = BIC_i - BIC_min
"""
struct BICWeighting <: WeightingMethod end

"""
Stacked regression weighting (Breiman 1996).
Weights determined by minimizing cross-validated prediction error.
"""
struct StackedWeighting <: WeightingMethod
    n_folds::Int
    StackedWeighting(n_folds::Int=5) = new(n_folds)
end

# ============================================================================
# Model Specification
# ============================================================================

"""
A candidate model for model averaging.

Fields:
- name: Model identifier
- n_params: Number of parameters (k)
- log_likelihood: Log-likelihood at MLE
- aic: Akaike Information Criterion
- bic: Bayesian Information Criterion
- aicc: Corrected AIC (for small samples)
- theta: Parameter estimates
- predictions: Model predictions (optional)
- estimation_result: Full estimation result (optional)
"""
struct CandidateModel
    name::String
    n_params::Int
    n_obs::Int
    log_likelihood::Float64
    aic::Float64
    bic::Float64
    aicc::Float64
    theta::Vector{Float64}
    theta_se::Union{Nothing, Vector{Float64}}
    predictions::Union{Nothing, Vector{Float64}}
    estimation_result::Any  # EstimationResult type

    function CandidateModel(
        name::String,
        n_params::Int,
        n_obs::Int,
        log_likelihood::Float64;
        theta::Vector{Float64}=Float64[],
        theta_se::Union{Nothing, Vector{Float64}}=nothing,
        predictions::Union{Nothing, Vector{Float64}}=nothing,
        estimation_result=nothing
    )
        k = n_params
        n = n_obs

        # AIC = -2 * LL + 2k
        aic = -2.0 * log_likelihood + 2.0 * k

        # BIC = -2 * LL + k * ln(n)
        bic = -2.0 * log_likelihood + k * log(n)

        # AICc = AIC + 2k(k+1)/(n-k-1)  (corrected for small samples)
        if n > k + 1
            aicc = aic + 2.0 * k * (k + 1) / (n - k - 1)
        else
            aicc = Inf  # Not defined for small samples
        end

        new(name, n_params, n_obs, log_likelihood, aic, bic, aicc,
            theta, theta_se, predictions, estimation_result)
    end
end

# ============================================================================
# Model Averaging Specification
# ============================================================================

"""
Specification for model averaging.
"""
struct ModelAveragingSpec
    method::WeightingMethod
    use_aicc::Bool  # Use AICc instead of AIC for AICWeighting
    min_weight::Float64  # Minimum weight threshold (models below this are excluded)

    function ModelAveragingSpec(;
        method::WeightingMethod=AICWeighting(),
        use_aicc::Bool=true,
        min_weight::Float64=0.0
    )
        @assert 0.0 <= min_weight < 1.0 "min_weight must be in [0, 1)"
        new(method, use_aicc, min_weight)
    end
end

# ============================================================================
# Model Averaging Results
# ============================================================================

"""
Result of model averaging.
"""
struct ModelAveragingResult
    models::Vector{CandidateModel}
    weights::Vector{Float64}
    averaged_predictions::Union{Nothing, Vector{Float64}}
    prediction_uncertainty::Union{Nothing, Vector{Float64}}
    best_model_index::Int
    evidence_ratios::Vector{Float64}  # w_best / w_i
    cumulative_weight::Vector{Float64}  # Cumulative weights for confidence set
end

# ============================================================================
# Weight Computation
# ============================================================================

"""
    compute_model_weights(models, method)

Compute model weights using the specified weighting method.
"""
function compute_model_weights(
    models::Vector{CandidateModel},
    method::AICWeighting;
    use_aicc::Bool=true
)::Vector{Float64}
    n_models = length(models)

    # Get information criterion values
    if use_aicc
        ic_values = [m.aicc for m in models]
    else
        ic_values = [m.aic for m in models]
    end

    # Handle Inf values
    valid_mask = isfinite.(ic_values)
    if !any(valid_mask)
        return ones(n_models) / n_models
    end

    # Compute deltas relative to minimum
    ic_min = minimum(ic_values[valid_mask])
    deltas = ic_values .- ic_min

    # Compute Akaike weights
    # w_i = exp(-0.5 * Δ_i) / Σ exp(-0.5 * Δ_j)
    raw_weights = exp.(-0.5 .* deltas)
    raw_weights[.!valid_mask] .= 0.0

    # Normalize
    total = sum(raw_weights)
    if total > 0
        return raw_weights ./ total
    else
        return ones(n_models) / n_models
    end
end

function compute_model_weights(
    models::Vector{CandidateModel},
    method::BICWeighting;
    use_aicc::Bool=false  # Ignored for BIC
)::Vector{Float64}
    n_models = length(models)

    bic_values = [m.bic for m in models]

    # Handle Inf values
    valid_mask = isfinite.(bic_values)
    if !any(valid_mask)
        return ones(n_models) / n_models
    end

    # Compute deltas relative to minimum
    bic_min = minimum(bic_values[valid_mask])
    deltas = bic_values .- bic_min

    # Compute BIC weights (same formula as AIC weights)
    raw_weights = exp.(-0.5 .* deltas)
    raw_weights[.!valid_mask] .= 0.0

    # Normalize
    total = sum(raw_weights)
    if total > 0
        return raw_weights ./ total
    else
        return ones(n_models) / n_models
    end
end

function compute_model_weights(
    models::Vector{CandidateModel},
    method::StackedWeighting;
    use_aicc::Bool=false  # Ignored for stacking
)::Vector{Float64}
    # Stacked regression (Breiman 1996): find weights that minimize
    # cross-validated prediction error using NNLS optimization

    n_models = length(models)

    if n_models == 1
        return [1.0]
    end

    # Check if predictions are available
    has_preds = all(m -> m.predictions !== nothing, models)

    if !has_preds
        # Fall back to AIC weighting if no predictions
        return compute_model_weights(models, AICWeighting(); use_aicc=true)
    end

    # Try to extract observed values from estimation results
    y = _extract_observations(models)

    if y === nothing
        # Fall back to AIC weighting if no observations available
        return compute_model_weights(models, AICWeighting(); use_aicc=true)
    end

    n_obs = length(y)

    # Build prediction matrix X (n_obs × n_models)
    X = zeros(n_obs, n_models)
    for (i, model) in enumerate(models)
        if length(model.predictions) != n_obs
            # Inconsistent prediction lengths - fall back to AIC
            return compute_model_weights(models, AICWeighting(); use_aicc=true)
        end
        X[:, i] = model.predictions
    end

    # Solve: min ||Xw - y||² subject to w ≥ 0, sum(w) = 1
    weights = _nnls_simplex(X, y)

    return weights
end

"""
    _extract_observations(models)

Extract observed values from estimation results.
Returns nothing if observations cannot be extracted.
"""
function _extract_observations(models::Vector{CandidateModel})::Union{Nothing, Vector{Float64}}
    # Try to extract observations from the first model's estimation result
    if isempty(models)
        return nothing
    end

    est = models[1].estimation_result
    if est === nothing
        return nothing
    end

    # Try different ways to access observations
    # EstimationResult typically stores population data
    if hasproperty(est, :observations) && est.observations !== nothing
        return Vector{Float64}(est.observations)
    elseif hasproperty(est, :population) && est.population !== nothing
        pop = est.population
        if hasproperty(pop, :observations)
            return Vector{Float64}(pop.observations)
        elseif hasproperty(pop, :dv) && pop.dv !== nothing
            return Vector{Float64}(pop.dv)
        end
    elseif hasproperty(est, :data) && est.data !== nothing
        data = est.data
        if hasproperty(data, :dv)
            return Vector{Float64}(data.dv)
        end
    end

    return nothing
end

"""
    _nnls_simplex(X, y; max_iter=1000, tol=1e-8)

Solve non-negative least squares with simplex constraint:
    min ||Xw - y||² subject to w ≥ 0, sum(w) = 1

Uses projected gradient descent with adaptive step size.
"""
function _nnls_simplex(
    X::Matrix{Float64},
    y::Vector{Float64};
    max_iter::Int=1000,
    tol::Float64=1e-8
)::Vector{Float64}
    n_obs, n_models = size(X)

    # Initialize with equal weights
    w = ones(n_models) / n_models

    # Precompute X'X and X'y for efficiency
    XtX = X' * X
    Xty = X' * y

    # Lipschitz constant for step size
    L = opnorm(XtX)
    if L < 1e-10
        return w  # Degenerate case
    end
    step_size = 1.0 / L

    prev_obj = Inf

    for iter in 1:max_iter
        # Gradient: ∇f = 2(X'Xw - X'y) = 2*X'*(Xw - y)
        grad = 2.0 * (XtX * w - Xty)

        # Gradient step
        w_new = w - step_size * grad

        # Project onto probability simplex
        w_new = _project_to_simplex(w_new)

        # Compute objective for convergence check
        residual = X * w_new - y
        obj = dot(residual, residual)

        # Check convergence
        if abs(prev_obj - obj) < tol * max(1.0, prev_obj)
            w = w_new
            break
        end

        prev_obj = obj
        w = w_new
    end

    # Ensure strict non-negativity and normalization
    w = max.(w, 0.0)
    w_sum = sum(w)
    if w_sum > 0
        w ./= w_sum
    else
        # Fallback to equal weights if something went wrong
        w = ones(n_models) / n_models
    end

    return w
end

"""
    _project_to_simplex(v)

Project vector v onto the probability simplex {w : w ≥ 0, sum(w) = 1}.
Uses the algorithm from Duchi et al. (2008).
"""
function _project_to_simplex(v::Vector{Float64})::Vector{Float64}
    n = length(v)

    # Sort in descending order
    u = sort(v, rev=true)

    # Find the threshold
    cssv = cumsum(u)
    rho = 0
    for j in 1:n
        if u[j] + (1.0 - cssv[j]) / j > 0
            rho = j
        end
    end

    # Compute threshold
    theta = (cssv[rho] - 1.0) / rho

    # Apply thresholding
    return max.(v .- theta, 0.0)
end

# ============================================================================
# Prediction Averaging
# ============================================================================

"""
    average_predictions(models, weights)

Compute weighted average of model predictions.

Returns (averaged_predictions, prediction_uncertainty)
"""
function average_predictions(
    models::Vector{CandidateModel},
    weights::Vector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}}
    # Check if all models have predictions
    has_preds = all(m -> m.predictions !== nothing, models)
    if !has_preds
        return (Float64[], Float64[])
    end

    n_obs = length(models[1].predictions)
    n_models = length(models)

    # Weighted average: E[Y] = Σ w_i * pred_i
    avg_pred = zeros(n_obs)
    for (i, model) in enumerate(models)
        avg_pred .+= weights[i] .* model.predictions
    end

    # Prediction uncertainty includes:
    # 1. Within-model variance (if available)
    # 2. Between-model variance

    # Between-model variance: Var = Σ w_i * (pred_i - avg)^2
    between_var = zeros(n_obs)
    for (i, model) in enumerate(models)
        between_var .+= weights[i] .* (model.predictions .- avg_pred).^2
    end

    pred_uncertainty = sqrt.(between_var)

    return (avg_pred, pred_uncertainty)
end

# ============================================================================
# Parameter Averaging
# ============================================================================

"""
    average_parameters(models, weights, param_names)

Compute weighted average of parameters across models.

Only averages parameters that appear in all models with the same interpretation.
Returns Dict(param_name => (avg_value, avg_se))
"""
function average_parameters(
    models::Vector{CandidateModel},
    weights::Vector{Float64},
    param_names::Vector{Symbol}
)::Dict{Symbol, Tuple{Float64, Float64}}
    result = Dict{Symbol, Tuple{Float64, Float64}}()

    n_models = length(models)

    for (j, pname) in enumerate(param_names)
        # Check if all models have this parameter
        values = Float64[]
        se_values = Float64[]

        for model in models
            if j <= length(model.theta)
                push!(values, model.theta[j])
                if model.theta_se !== nothing && j <= length(model.theta_se)
                    push!(se_values, model.theta_se[j])
                end
            end
        end

        if length(values) == n_models
            # Weighted average of parameter
            avg_value = sum(weights .* values)

            # Averaged SE includes:
            # 1. Average of within-model SEs (conditional SE)
            # 2. Between-model variance (unconditional component)

            if length(se_values) == n_models
                # Within-model average SE
                avg_within_var = sum(weights .* (se_values.^2))

                # Between-model variance
                between_var = sum(weights .* ((values .- avg_value).^2))

                # Total variance (unconditional)
                total_var = avg_within_var + between_var
                avg_se = sqrt(total_var)
            else
                avg_se = NaN
            end

            result[pname] = (avg_value, avg_se)
        end
    end

    return result
end

# ============================================================================
# Model Selection Table
# ============================================================================

"""
    model_selection_table(models, spec)

Generate a model selection table with ranks, weights, and evidence ratios.
"""
function model_selection_table(
    models::Vector{CandidateModel},
    spec::ModelAveragingSpec
)::Dict{Symbol, Any}
    weights = compute_model_weights(models, spec.method; use_aicc=spec.use_aicc)

    # Sort by weight (descending)
    sorted_idx = sortperm(weights, rev=true)

    # Evidence ratios (relative to best model)
    best_weight = maximum(weights)
    evidence_ratios = best_weight ./ max.(weights, 1e-10)

    # Cumulative weights (for confidence set)
    sorted_weights = weights[sorted_idx]
    cumulative = cumsum(sorted_weights)

    return Dict(
        :model_names => [models[i].name for i in sorted_idx],
        :n_params => [models[i].n_params for i in sorted_idx],
        :log_likelihood => [models[i].log_likelihood for i in sorted_idx],
        :aic => [models[i].aic for i in sorted_idx],
        :aicc => [models[i].aicc for i in sorted_idx],
        :bic => [models[i].bic for i in sorted_idx],
        :delta_aic => [models[i].aic - minimum(m.aic for m in models) for i in sorted_idx],
        :delta_bic => [models[i].bic - minimum(m.bic for m in models) for i in sorted_idx],
        :weights => sorted_weights,
        :evidence_ratios => evidence_ratios[sorted_idx],
        :cumulative_weight => cumulative,
        :ranks => 1:length(models)
    )
end

"""
    best_model(models, spec)

Return the best model according to the weighting method.
"""
function best_model(
    models::Vector{CandidateModel},
    spec::ModelAveragingSpec
)::CandidateModel
    weights = compute_model_weights(models, spec.method; use_aicc=spec.use_aicc)
    best_idx = argmax(weights)
    return models[best_idx]
end

# ============================================================================
# Full Model Averaging
# ============================================================================

"""
    model_averaging(models, spec)

Perform full model averaging analysis.
"""
function model_averaging(
    models::Vector{CandidateModel},
    spec::ModelAveragingSpec
)::ModelAveragingResult
    weights = compute_model_weights(models, spec.method; use_aicc=spec.use_aicc)

    # Apply minimum weight threshold
    weights[weights .< spec.min_weight] .= 0.0
    if sum(weights) > 0
        weights ./= sum(weights)
    end

    # Average predictions
    avg_pred, pred_unc = average_predictions(models, weights)

    # Best model
    best_idx = argmax(weights)
    best_weight = weights[best_idx]

    # Evidence ratios
    evidence_ratios = best_weight ./ max.(weights, 1e-10)

    # Cumulative weights (sorted)
    sorted_idx = sortperm(weights, rev=true)
    cumulative = cumsum(weights[sorted_idx])

    return ModelAveragingResult(
        models,
        weights,
        isempty(avg_pred) ? nothing : avg_pred,
        isempty(pred_unc) ? nothing : pred_unc,
        best_idx,
        evidence_ratios,
        cumulative
    )
end

export model_averaging

# ============================================================================
# Confidence Set
# ============================================================================

"""
    confidence_set(result; level=0.95)

Get the set of models that contain the best model with specified confidence.

Returns indices of models in the confidence set.
"""
function confidence_set(
    result::ModelAveragingResult;
    level::Float64=0.95
)::Vector{Int}
    sorted_idx = sortperm(result.weights, rev=true)
    cumulative = cumsum(result.weights[sorted_idx])

    # Find models until cumulative weight >= level
    n_in_set = findfirst(c -> c >= level, cumulative)
    if n_in_set === nothing
        n_in_set = length(sorted_idx)
    end

    return sorted_idx[1:n_in_set]
end

export confidence_set

# ============================================================================
# Relative Importance
# ============================================================================

"""
    parameter_importance(models, weights, param_name)

Compute the relative importance of a parameter.

Importance = sum of weights for models containing the parameter.
"""
function parameter_importance(
    models::Vector{CandidateModel},
    weights::Vector{Float64},
    param_index::Int
)::Float64
    importance = 0.0
    for (i, model) in enumerate(models)
        if param_index <= length(model.theta)
            importance += weights[i]
        end
    end
    return importance
end

export parameter_importance
