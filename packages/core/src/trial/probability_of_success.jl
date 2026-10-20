# Probability of Success (PoS) Module
# Bayesian framework for Go/No-Go decision support in clinical development
# Uses Bayesian assurance and predictive probability methods

export PoSConfig, PoSResult, AssuranceResult
export PoSPrior, NormalPoSPrior, BetaPoSPrior, HistoricalDataPrior
export SuccessCriterion, SuperioritySuccess, NonInferioritySuccess, EquivalenceSuccess
export compute_pos, compute_assurance, compute_predictive_probability
export sensitivity_analysis_pos, pos_decision_boundary
export create_pos_prior_from_trials, summarize_pos_result
export to_normal_prior, prior_to_distribution, compute_power

using Statistics
using LinearAlgebra
using Random
using Distributions

# ============================================================================
# Prior Specification Types
# ============================================================================

"""
Abstract type for Probability of Success priors.
"""
abstract type PoSPrior end

"""
    NormalPoSPrior

Normal prior for treatment effect.

# Fields
- `mean::Float64`: Prior mean for treatment effect
- `variance::Float64`: Prior variance
- `source::Symbol`: Source of prior (:meta_analysis, :previous_trial, :expert, :historical)
"""
struct NormalPoSPrior <: PoSPrior
    mean::Float64
    variance::Float64
    source::Symbol

    function NormalPoSPrior(mean::Float64, variance::Float64;
                            source::Symbol = :expert)
        @assert variance > 0 "Variance must be positive"
        new(mean, variance, source)
    end
end

"""
    BetaPoSPrior

Beta prior for response rate (binary outcomes).

# Fields
- `alpha::Float64`: Alpha parameter (pseudo successes)
- `beta::Float64`: Beta parameter (pseudo failures)
- `source::Symbol`: Source of prior
"""
struct BetaPoSPrior <: PoSPrior
    alpha::Float64
    beta::Float64
    source::Symbol

    function BetaPoSPrior(alpha::Float64, beta::Float64;
                          source::Symbol = :expert)
        @assert alpha > 0 && beta > 0 "Alpha and beta must be positive"
        new(alpha, beta, source)
    end
end

"""
    HistoricalDataPrior

Prior derived from historical trial data.

# Fields
- `effect_estimates::Vector{Float64}`: Observed effects from historical trials
- `standard_errors::Vector{Float64}`: SEs of effect estimates
- `discount_factor::Float64`: Discount for prior variance (0-1, 1=full weight)
- `source::Symbol`: :historical
"""
struct HistoricalDataPrior <: PoSPrior
    effect_estimates::Vector{Float64}
    standard_errors::Vector{Float64}
    discount_factor::Float64
    source::Symbol

    function HistoricalDataPrior(effect_estimates::Vector{Float64},
                                  standard_errors::Vector{Float64};
                                  discount_factor::Float64 = 0.5)
        @assert length(effect_estimates) == length(standard_errors)
        @assert all(se -> se > 0, standard_errors)
        @assert 0 < discount_factor <= 1
        new(effect_estimates, standard_errors, discount_factor, :historical)
    end
end

"""
Convert HistoricalDataPrior to Normal prior using meta-analysis.
"""
function to_normal_prior(prior::HistoricalDataPrior)::NormalPoSPrior
    n = length(prior.effect_estimates)

    # Fixed-effects meta-analysis
    weights = 1.0 ./ (prior.standard_errors.^2)
    pooled_mean = sum(weights .* prior.effect_estimates) / sum(weights)
    pooled_var = 1.0 / sum(weights)

    # Apply discount factor (inflate variance)
    discounted_var = pooled_var / prior.discount_factor

    return NormalPoSPrior(pooled_mean, discounted_var; source = :historical)
end

# ============================================================================
# Success Criterion Types
# ============================================================================

"""
Abstract type for success criteria.
"""
abstract type SuccessCriterion end

"""
    SuperioritySuccess

Superiority success criterion.
Success if: effect > margin with probability > threshold

# Fields
- `margin::Float64`: Superiority margin (typically 0 for strict superiority)
- `one_sided::Bool`: One-sided test
- `alpha::Float64`: Significance level
"""
struct SuperioritySuccess <: SuccessCriterion
    margin::Float64
    one_sided::Bool
    alpha::Float64

    function SuperioritySuccess(; margin::Float64 = 0.0,
                                  one_sided::Bool = true,
                                  alpha::Float64 = 0.025)
        new(margin, one_sided, alpha)
    end
end

"""
    NonInferioritySuccess

Non-inferiority success criterion.
Success if: effect > -margin

# Fields
- `margin::Float64`: Non-inferiority margin (positive value, test vs -margin)
- `alpha::Float64`: Significance level
"""
struct NonInferioritySuccess <: SuccessCriterion
    margin::Float64
    alpha::Float64

    function NonInferioritySuccess(; margin::Float64 = 0.1,
                                     alpha::Float64 = 0.025)
        @assert margin > 0 "Margin must be positive"
        new(margin, alpha)
    end
end

"""
    EquivalenceSuccess

Equivalence (bioequivalence) success criterion.
Success if: -margin < effect < margin

# Fields
- `lower_margin::Float64`: Lower equivalence margin
- `upper_margin::Float64`: Upper equivalence margin
- `alpha::Float64`: Significance level
"""
struct EquivalenceSuccess <: SuccessCriterion
    lower_margin::Float64
    upper_margin::Float64
    alpha::Float64

    function EquivalenceSuccess(; lower_margin::Float64 = -0.2,
                                  upper_margin::Float64 = 0.2,
                                  alpha::Float64 = 0.05)
        @assert lower_margin < upper_margin
        new(lower_margin, upper_margin, alpha)
    end
end

# ============================================================================
# Configuration and Result Types
# ============================================================================

"""
    PoSConfig

Configuration for Probability of Success calculation.

# Fields
- `prior::PoSPrior`: Prior distribution for treatment effect
- `criterion::SuccessCriterion`: Success criterion
- `n_samples::Int`: Number of Monte Carlo samples
- `n_subjects_planned::Int`: Planned sample size for future trial
- `include_uncertainty::Bool`: Include parameter uncertainty
- `seed::Int`: Random seed
"""
struct PoSConfig
    prior::PoSPrior
    criterion::SuccessCriterion
    n_samples::Int
    n_subjects_planned::Int
    include_uncertainty::Bool
    seed::Int

    function PoSConfig(;
            prior::PoSPrior,
            criterion::SuccessCriterion = SuperioritySuccess(),
            n_samples::Int = 10000,
            n_subjects_planned::Int = 100,
            include_uncertainty::Bool = true,
            seed::Int = 12345)
        new(prior, criterion, n_samples, n_subjects_planned, include_uncertainty, seed)
    end
end

"""
    AssuranceResult

Result from assurance calculation (Bayesian predictive probability of success).

# Fields
- `assurance::Float64`: Probability of trial success (0-1)
- `assurance_ci::Tuple{Float64, Float64}`: 95% CI for assurance
- `predictive_mean::Float64`: Predictive mean for treatment effect
- `predictive_sd::Float64`: Predictive SD
- `predictive_quantiles::Dict{Float64, Float64}`: Quantiles of predictive distribution
- `probability_above_threshold::Float64`: P(effect > threshold)
- `expected_power::Float64`: Expected power across prior
- `n_samples::Int`: Number of MC samples used
"""
struct AssuranceResult
    assurance::Float64
    assurance_ci::Tuple{Float64, Float64}
    predictive_mean::Float64
    predictive_sd::Float64
    predictive_quantiles::Dict{Float64, Float64}
    probability_above_threshold::Float64
    expected_power::Float64
    n_samples::Int
end

"""
    PoSResult

Complete Probability of Success analysis result.

# Fields
- `config::PoSConfig`: Configuration used
- `assurance::AssuranceResult`: Assurance calculation result
- `posterior_samples::Vector{Float64}`: Posterior samples of treatment effect
- `power_curve::Dict{Float64, Float64}`: Power as function of true effect
- `go_nogo_recommendation::Symbol`: :go, :nogo, or :uncertain
- `decision_rationale::String`: Explanation of recommendation
- `sensitivity_results::Union{Nothing, Dict{String, AssuranceResult}}`: Sensitivity analysis
"""
struct PoSResult
    config::PoSConfig
    assurance::AssuranceResult
    posterior_samples::Vector{Float64}
    power_curve::Dict{Float64, Float64}
    go_nogo_recommendation::Symbol
    decision_rationale::String
    sensitivity_results::Union{Nothing, Dict{String, AssuranceResult}}
end

# ============================================================================
# Power Calculation Functions
# ============================================================================

"""
    compute_power(effect::Float64, se::Float64, criterion::SuccessCriterion)

Compute power for a given true effect size.
"""
function compute_power(effect::Float64, se::Float64, criterion::SuperioritySuccess)
    # Z = (effect - margin) / se
    z_crit = quantile(Normal(), 1 - criterion.alpha)
    z = (effect - criterion.margin) / se

    return cdf(Normal(), z - z_crit)
end

function compute_power(effect::Float64, se::Float64, criterion::NonInferioritySuccess)
    # Test: H0: effect <= -margin vs H1: effect > -margin
    z_crit = quantile(Normal(), 1 - criterion.alpha)
    z = (effect + criterion.margin) / se

    return cdf(Normal(), z - z_crit)
end

function compute_power(effect::Float64, se::Float64, criterion::EquivalenceSuccess)
    # Two one-sided tests (TOST)
    z_crit = quantile(Normal(), 1 - criterion.alpha)

    # Lower bound test
    z_lower = (effect - criterion.lower_margin) / se
    power_lower = cdf(Normal(), z_lower - z_crit)

    # Upper bound test
    z_upper = (criterion.upper_margin - effect) / se
    power_upper = cdf(Normal(), z_upper - z_crit)

    # Both must pass
    return power_lower * power_upper
end

# ============================================================================
# Assurance Calculation
# ============================================================================

"""
    compute_assurance(config::PoSConfig)

Compute Bayesian assurance (average power weighted by prior).

Assurance = E_prior[Power(θ)]
         = ∫ Power(θ) × π(θ) dθ

# Arguments
- `config`: PoS configuration with prior and criterion

# Returns
- `AssuranceResult`
"""
function compute_assurance(config::PoSConfig)
    rng = MersenneTwister(config.seed)

    # Get prior distribution
    prior_dist = prior_to_distribution(config.prior)

    # Sample from prior
    n_samples = config.n_samples
    prior_samples = rand(rng, prior_dist, n_samples)

    # Expected SE for planned trial
    n_arm = config.n_subjects_planned ÷ 2
    expected_se = estimate_expected_se(config.prior, n_arm)

    # Compute power for each prior sample
    powers = zeros(n_samples)
    successes = zeros(n_samples)

    for i in 1:n_samples
        effect = prior_samples[i]
        powers[i] = compute_power(effect, expected_se, config.criterion)

        if config.include_uncertainty
            # Also sample from sampling distribution
            observed_effect = effect + expected_se * randn(rng)
            successes[i] = would_succeed(observed_effect, expected_se, config.criterion)
        else
            successes[i] = powers[i]
        end
    end

    # Assurance = mean power across prior
    assurance = mean(successes)
    expected_power = mean(powers)

    # Bootstrap CI for assurance
    n_bootstrap = 1000
    assurance_bootstrap = zeros(n_bootstrap)
    for b in 1:n_bootstrap
        idx = rand(rng, 1:n_samples, n_samples)
        assurance_bootstrap[b] = mean(successes[idx])
    end
    assurance_ci = (quantile(assurance_bootstrap, 0.025),
                    quantile(assurance_bootstrap, 0.975))

    # Predictive distribution statistics
    predictive_mean = mean(prior_samples)
    predictive_sd = std(prior_samples)

    quantile_levels = [0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975]
    predictive_quantiles = Dict{Float64, Float64}()
    for q in quantile_levels
        predictive_quantiles[q] = quantile(prior_samples, q)
    end

    # P(effect > threshold)
    threshold = get_success_threshold(config.criterion)
    prob_above = mean(prior_samples .> threshold)

    return AssuranceResult(
        assurance,
        assurance_ci,
        predictive_mean,
        predictive_sd,
        predictive_quantiles,
        prob_above,
        expected_power,
        n_samples
    )
end

"""
Convert prior to Distribution object.
"""
function prior_to_distribution(prior::NormalPoSPrior)
    return Normal(prior.mean, sqrt(prior.variance))
end

function prior_to_distribution(prior::BetaPoSPrior)
    return Beta(prior.alpha, prior.beta)
end

function prior_to_distribution(prior::HistoricalDataPrior)
    normal_prior = to_normal_prior(prior)
    return Normal(normal_prior.mean, sqrt(normal_prior.variance))
end

"""
Estimate expected standard error for planned trial.
"""
function estimate_expected_se(prior::PoSPrior, n_per_arm::Int)
    # For two-arm trial, assuming equal allocation
    # SE = sqrt(2 * sigma^2 / n) where sigma is residual SD

    # Use prior to estimate plausible variability
    if prior isa NormalPoSPrior
        # Assume effect variance reflects prior uncertainty
        # Use pooled SD estimate
        sigma = sqrt(prior.variance) * 2  # Rough scaling
    elseif prior isa BetaPoSPrior
        # For binary, SE = sqrt(p*(1-p)/n)
        p = prior.alpha / (prior.alpha + prior.beta)
        return sqrt(2 * p * (1 - p) / n_per_arm)
    else
        sigma = 1.0  # Default
    end

    return sqrt(2 * sigma^2 / n_per_arm)
end

"""
Check if observed effect would lead to success.
"""
function would_succeed(observed_effect::Float64, se::Float64,
                       criterion::SuperioritySuccess)
    z_crit = quantile(Normal(), 1 - criterion.alpha)
    z = (observed_effect - criterion.margin) / se
    return z > z_crit ? 1.0 : 0.0
end

function would_succeed(observed_effect::Float64, se::Float64,
                       criterion::NonInferioritySuccess)
    z_crit = quantile(Normal(), 1 - criterion.alpha)
    z = (observed_effect + criterion.margin) / se
    return z > z_crit ? 1.0 : 0.0
end

function would_succeed(observed_effect::Float64, se::Float64,
                       criterion::EquivalenceSuccess)
    z_crit = quantile(Normal(), 1 - criterion.alpha)

    z_lower = (observed_effect - criterion.lower_margin) / se
    z_upper = (criterion.upper_margin - observed_effect) / se

    return (z_lower > z_crit && z_upper > z_crit) ? 1.0 : 0.0
end

"""
Get success threshold from criterion.
"""
function get_success_threshold(criterion::SuperioritySuccess)
    return criterion.margin
end

function get_success_threshold(criterion::NonInferioritySuccess)
    return -criterion.margin
end

function get_success_threshold(criterion::EquivalenceSuccess)
    return 0.0  # Center of equivalence region
end

# ============================================================================
# Main PoS Calculation
# ============================================================================

"""
    compute_pos(config::PoSConfig; run_sensitivity::Bool = false)

Compute complete Probability of Success analysis.

# Arguments
- `config`: PoS configuration
- `run_sensitivity`: Whether to run sensitivity analysis

# Returns
- `PoSResult`

# Example
```julia
# Define prior from Phase 2 data
prior = NormalPoSPrior(0.3, 0.04; source = :previous_trial)

# Define success criterion
criterion = SuperioritySuccess(margin = 0.0, alpha = 0.025)

config = PoSConfig(
    prior = prior,
    criterion = criterion,
    n_subjects_planned = 400
)

result = compute_pos(config)
println("Assurance: \$(result.assurance.assurance * 100)%")
println("Recommendation: \$(result.go_nogo_recommendation)")
```
"""
function compute_pos(config::PoSConfig; run_sensitivity::Bool = false)
    rng = MersenneTwister(config.seed)

    # Compute assurance
    assurance = compute_assurance(config)

    # Generate posterior samples
    prior_dist = prior_to_distribution(config.prior)
    posterior_samples = rand(rng, prior_dist, config.n_samples)

    # Compute power curve
    power_curve = Dict{Float64, Float64}()
    n_arm = config.n_subjects_planned ÷ 2
    expected_se = estimate_expected_se(config.prior, n_arm)

    effect_range = range(-0.5, 1.0, length = 30)
    for effect in effect_range
        power_curve[effect] = compute_power(effect, expected_se, config.criterion)
    end

    # Generate recommendation
    go_nogo, rationale = generate_recommendation(assurance, config)

    # Sensitivity analysis
    sensitivity_results = nothing
    if run_sensitivity
        sensitivity_results = sensitivity_analysis_pos(config)
    end

    return PoSResult(
        config,
        assurance,
        posterior_samples,
        power_curve,
        go_nogo,
        rationale,
        sensitivity_results
    )
end

"""
Generate go/no-go recommendation based on assurance.
"""
function generate_recommendation(assurance::AssuranceResult, config::PoSConfig)
    pos = assurance.assurance

    if pos >= 0.80
        return :go, "Strong assurance ($(round(pos*100, digits=1))%) supports proceeding. " *
                    "Expected power: $(round(assurance.expected_power*100, digits=1))%."

    elseif pos >= 0.60
        return :go, "Moderate assurance ($(round(pos*100, digits=1))%) suggests proceeding with caution. " *
                    "Consider risk mitigation strategies."

    elseif pos >= 0.40
        return :uncertain, "Uncertain outlook ($(round(pos*100, digits=1))% assurance). " *
                          "Additional data or design modifications recommended before decision."

    elseif pos >= 0.25
        return :nogo, "Low assurance ($(round(pos*100, digits=1))%). " *
                      "Significant changes to program needed before proceeding."

    else
        return :nogo, "Very low assurance ($(round(pos*100, digits=1))%). " *
                      "Program unlikely to succeed. Consider discontinuation."
    end
end

# ============================================================================
# Sensitivity Analysis
# ============================================================================

"""
    sensitivity_analysis_pos(config::PoSConfig)

Perform sensitivity analysis on PoS by varying key assumptions.

Varies:
- Prior mean (+/- 20%)
- Prior variance (2x, 0.5x)
- Sample size (+/- 25%)
"""
function sensitivity_analysis_pos(config::PoSConfig)
    results = Dict{String, AssuranceResult}()

    # Get baseline prior parameters
    if config.prior isa NormalPoSPrior
        base_mean = config.prior.mean
        base_var = config.prior.variance

        # Prior mean sensitivity
        for factor in [0.8, 1.2]
            modified_prior = NormalPoSPrior(base_mean * factor, base_var;
                                            source = config.prior.source)
            modified_config = PoSConfig(
                prior = modified_prior,
                criterion = config.criterion,
                n_samples = config.n_samples,
                n_subjects_planned = config.n_subjects_planned,
                seed = config.seed
            )
            label = "Prior mean $(Int(factor*100))%"
            results[label] = compute_assurance(modified_config)
        end

        # Prior variance sensitivity
        for factor in [0.5, 2.0]
            modified_prior = NormalPoSPrior(base_mean, base_var * factor;
                                            source = config.prior.source)
            modified_config = PoSConfig(
                prior = modified_prior,
                criterion = config.criterion,
                n_samples = config.n_samples,
                n_subjects_planned = config.n_subjects_planned,
                seed = config.seed
            )
            label = "Prior variance $(Int(factor*100))%"
            results[label] = compute_assurance(modified_config)
        end
    end

    # Sample size sensitivity
    base_n = config.n_subjects_planned
    for factor in [0.75, 1.25]
        modified_config = PoSConfig(
            prior = config.prior,
            criterion = config.criterion,
            n_samples = config.n_samples,
            n_subjects_planned = Int(round(base_n * factor)),
            seed = config.seed
        )
        label = "Sample size $(Int(factor*100))%"
        results[label] = compute_assurance(modified_config)
    end

    return results
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
    create_pos_prior_from_trials(effects::Vector{Float64},
                                  standard_errors::Vector{Float64};
                                  discount::Float64 = 0.5)

Create PoS prior from historical trial results using meta-analysis.

# Arguments
- `effects`: Estimated treatment effects from historical trials
- `standard_errors`: Standard errors of the estimates
- `discount`: Discount factor for prior (0-1, lower = more skeptical)

# Returns
- `NormalPoSPrior`
"""
function create_pos_prior_from_trials(effects::Vector{Float64},
                                       standard_errors::Vector{Float64};
                                       discount::Float64 = 0.5)
    hist_prior = HistoricalDataPrior(effects, standard_errors;
                                      discount_factor = discount)
    return to_normal_prior(hist_prior)
end

"""
    pos_decision_boundary(config::PoSConfig; target_assurance::Float64 = 0.80)

Find the minimum sample size needed to achieve target assurance.

Returns estimated required sample size.
"""
function pos_decision_boundary(config::PoSConfig;
                               target_assurance::Float64 = 0.80)
    # Binary search for required sample size
    min_n = 20
    max_n = 2000

    for _ in 1:20
        mid_n = (min_n + max_n) ÷ 2

        test_config = PoSConfig(
            prior = config.prior,
            criterion = config.criterion,
            n_samples = config.n_samples ÷ 10,  # Faster for search
            n_subjects_planned = mid_n,
            seed = config.seed
        )

        assurance = compute_assurance(test_config)

        if assurance.assurance >= target_assurance
            max_n = mid_n
        else
            min_n = mid_n
        end

        if max_n - min_n <= 10
            break
        end
    end

    return max_n
end

"""
    summarize_pos_result(result::PoSResult)

Generate summary report of PoS analysis.
"""
function summarize_pos_result(result::PoSResult)
    a = result.assurance

    summary = """
    Probability of Success Analysis Summary
    ========================================

    Prior Information:
    - Type: $(typeof(result.config.prior))
    - Source: $(result.config.prior.source)

    Planned Trial:
    - Sample size: $(result.config.n_subjects_planned)
    - Success criterion: $(typeof(result.config.criterion))

    Results:
    - Assurance: $(round(a.assurance * 100, digits=1))%
    - 95% CI: ($(round(a.assurance_ci[1] * 100, digits=1))%, $(round(a.assurance_ci[2] * 100, digits=1))%)
    - Expected Power: $(round(a.expected_power * 100, digits=1))%

    Predictive Distribution:
    - Mean: $(round(a.predictive_mean, digits=3))
    - SD: $(round(a.predictive_sd, digits=3))
    - 95% PI: ($(round(a.predictive_quantiles[0.025], digits=3)), $(round(a.predictive_quantiles[0.975], digits=3)))

    Recommendation: $(result.go_nogo_recommendation)
    $(result.decision_rationale)
    """

    return summary
end

"""
    compute_predictive_probability(current_data::Dict{Symbol, Float64},
                                   final_n::Int,
                                   config::PoSConfig)

Compute predictive probability of success given interim data.

# Arguments
- `current_data`: Dict with :n_current, :effect_current, :se_current
- `final_n`: Planned final sample size
- `config`: PoS configuration

# Returns
- Float64: Predictive probability of success at final analysis
"""
function compute_predictive_probability(current_data::Dict{Symbol, Float64},
                                        final_n::Int,
                                        config::PoSConfig)
    n_current = Int(current_data[:n_current])
    effect_current = current_data[:effect_current]
    se_current = current_data[:se_current]

    rng = MersenneTwister(config.seed)
    n_samples = config.n_samples

    # Posterior given current data
    prior_dist = prior_to_distribution(config.prior)
    prior_precision = 1.0 / var(prior_dist)
    data_precision = 1.0 / se_current^2

    posterior_precision = prior_precision + data_precision
    posterior_var = 1.0 / posterior_precision
    posterior_mean = posterior_var * (prior_precision * mean(prior_dist) +
                                      data_precision * effect_current)

    posterior = Normal(posterior_mean, sqrt(posterior_var))

    # Simulate to final analysis
    n_remaining = final_n - n_current
    n_arm_remaining = n_remaining ÷ 2

    successes = 0.0
    for _ in 1:n_samples
        # Sample true effect from posterior
        true_effect = rand(rng, posterior)

        # Simulate remaining data
        remaining_se = estimate_expected_se(config.prior, n_arm_remaining)
        remaining_effect = true_effect + remaining_se * randn(rng)

        # Combine with current data
        final_effect = (n_current * effect_current + n_remaining * remaining_effect) / final_n
        final_se = sqrt(se_current^2 * n_current / final_n +
                       remaining_se^2 * n_remaining / final_n)

        # Check success
        if would_succeed(final_effect, final_se, config.criterion) > 0
            successes += 1
        end
    end

    return successes / n_samples
end
