# Adaptive Trial Simulation
# Industry-standard implementation of adaptive clinical trial designs
#
# Features:
# - Interim analyses with alpha spending (O'Brien-Fleming, Pocock, etc.)
# - Response-adaptive randomization (RAR)
# - Information-adaptive designs (sample size re-estimation)
# - Treatment selection (adaptive arm dropping)
# - Biomarker-driven enrichment designs
#
# References:
# - FDA Guidance: Adaptive Designs for Clinical Trials of Drugs and Biologics (2019)
# - Thall PF, Wathen JK (2007). Practical Bayesian adaptive randomisation
# - Mehta CR, Pocock SJ (2011). Adaptive increase in sample size

using StableRNGs
using SpecialFunctions: erf, erfinv
using Distributions: Beta, Normal, pdf, cdf, quantile

export simulate_adaptive_trial, AdaptiveTrialResult, InterimResult
export check_efficacy_stop, check_futility_stop, compute_conditional_power
export get_efficacy_boundary, get_alpha_spending
# Response-adaptive randomization exports
export RARSpec, ThallWathenRAR, DoublyAdaptiveBiasedCoin, compute_rar_probabilities
# Sample size re-estimation exports
export SSRSpec, ConditionalPowerSSR, VarianceBasedSSR, compute_sample_size_reestimation
# Treatment selection exports
export TreatmentSelectionSpec, DropLoserSelection, select_treatments
# Enrichment exports
export EnrichmentSpec, BiomarkerEnrichment, apply_enrichment, EnrichmentResult
# Full adaptive simulation
export simulate_adaptive_trial_full, TreatmentSelectionResult

# ============================================================================
# Result Types
# ============================================================================

"""
    InterimResult

Results from a single interim analysis.

# Fields
- `analysis_number`: Which interim analysis (1, 2, ...)
- `information_fraction`: Fraction of total planned information
- `n_enrolled`: Number of subjects enrolled at this interim
- `primary_effect`: Estimated treatment effect
- `primary_se`: Standard error of effect estimate
- `primary_z_stat`: Z-statistic for primary endpoint
- `p_value`: P-value for primary endpoint
- `efficacy_boundary`: Efficacy stopping boundary at this interim
- `futility_boundary`: Futility stopping boundary
- `conditional_power`: Conditional power estimate
- `stop_for_efficacy`: Whether stopped for efficacy
- `stop_for_futility`: Whether stopped for futility
"""
struct InterimResult
    analysis_number::Int
    information_fraction::Float64
    n_enrolled::Int
    primary_effect::Float64
    primary_se::Float64
    primary_z_stat::Float64
    p_value::Float64
    efficacy_boundary::Float64
    futility_boundary::Union{Nothing, Float64}
    conditional_power::Float64
    stop_for_efficacy::Bool
    stop_for_futility::Bool
end

"""
    AdaptiveTrialResult

Complete results from adaptive trial simulation.

# Fields
- `final_result`: Final endpoint analysis
- `interim_results`: Vector of InterimResult from each interim
- `stopped_early`: Whether trial stopped before planned completion
- `stop_reason`: Reason for early stopping (:efficacy, :futility, :none)
- `final_n`: Final number of subjects enrolled
- `initial_n`: Planned sample size
- `seed`: Random seed used
"""
struct AdaptiveTrialResult
    final_result::Dict{Symbol, Any}
    interim_results::Vector{InterimResult}
    stopped_early::Bool
    stop_reason::Union{Nothing, Symbol}
    final_n::Int
    initial_n::Int
    seed::UInt64
end

# ============================================================================
# Alpha Spending Functions
# ============================================================================

"""
    get_alpha_spending(spending_function::Symbol, info_fraction::Float64, alpha::Float64)

Compute cumulative alpha spent at a given information fraction.

# Arguments
- `spending_function`: Type of spending function (:obrien_fleming, :pocock, :haybittle_peto)
- `info_fraction`: Information fraction (0 to 1)
- `alpha`: Overall alpha level (typically 0.025 one-sided)

# Returns
- Cumulative alpha spent at this information fraction
"""
function get_alpha_spending(spending_function::Symbol, info_fraction::Float64, alpha::Float64)
    t = info_fraction

    if spending_function == :obrien_fleming
        # O'Brien-Fleming: α(t) = 1 - Φ(z_α/√t) for one-sided alpha
        # Very conservative early, spends most alpha at end
        z_alpha = sqrt(2) * erfinv(1 - 2*alpha)
        return 1 - 0.5 * (1 + erf(z_alpha / sqrt(t) / sqrt(2)))

    elseif spending_function == :pocock
        # Pocock: α(t) = α * log(1 + (e-1)*t)
        # Constant boundaries, more aggressive early
        return alpha * log(1 + (exp(1) - 1) * t)

    elseif spending_function == :haybittle_peto
        # Haybittle-Peto: Use z=3 for interim, full alpha at final
        # Very conservative for early analyses
        if t < 1.0
            return 0.001  # Approximately 2*(1-Φ(3))
        else
            return alpha
        end

    elseif spending_function == :linear
        # Linear spending (Kim-DeMets family with rho=1)
        return alpha * t

    elseif spending_function == :none
        # No alpha spending - all at final
        if t >= 1.0
            return alpha
        else
            return 0.0
        end

    else
        # Default to O'Brien-Fleming
        z_alpha = sqrt(2) * erfinv(1 - 2*alpha)
        return 1 - 0.5 * (1 + erf(z_alpha / sqrt(t) / sqrt(2)))
    end
end

"""
    get_efficacy_boundary(spending_function::Symbol, info_fraction::Float64, alpha::Float64)

Get the efficacy stopping boundary (critical z-value) at a given information fraction.

# Returns
- Z-value threshold for efficacy stopping
"""
function get_efficacy_boundary(spending_function::Symbol, info_fraction::Float64, alpha::Float64)
    # Compute incremental alpha for this look
    alpha_spent = get_alpha_spending(spending_function, info_fraction, alpha)

    # Convert to z-value
    if alpha_spent <= 0
        return Inf  # No stopping at this interim
    elseif alpha_spent >= 1
        return 0.0
    else
        # Z-value corresponding to this alpha (one-sided)
        return sqrt(2) * erfinv(1 - 2*alpha_spent)
    end
end

# ============================================================================
# Stopping Rules
# ============================================================================

"""
    check_efficacy_stop(z_stat::Float64, boundary::Float64)

Check if trial should stop for efficacy.

# Arguments
- `z_stat`: Observed z-statistic
- `boundary`: Efficacy boundary (z-value)

# Returns
- Boolean: true if should stop for efficacy
"""
function check_efficacy_stop(z_stat::Float64, boundary::Float64)::Bool
    return z_stat > boundary
end

"""
    check_futility_stop(conditional_power::Float64, threshold::Float64)

Check if trial should stop for futility based on conditional power.

# Arguments
- `conditional_power`: Estimated conditional power
- `threshold`: Futility threshold (typically 0.10 to 0.20)

# Returns
- Boolean: true if should stop for futility
"""
function check_futility_stop(conditional_power::Float64, threshold::Float64)::Bool
    return conditional_power < threshold
end

"""
    compute_conditional_power(effect::Float64, se::Float64, info_frac::Float64,
                               target_effect::Float64; alpha=0.025)

Compute conditional power given interim results.

Conditional power is the probability of reaching significance at the end
of the trial, given the current data and assuming a specific true effect.

# Arguments
- `effect`: Observed treatment effect at interim
- `se`: Standard error of effect at interim
- `info_frac`: Information fraction at interim
- `target_effect`: Assumed true effect for power calculation
- `alpha`: Significance level (one-sided)

# Returns
- Conditional power (0 to 1)
"""
function compute_conditional_power(
    effect::Float64,
    se::Float64,
    info_frac::Float64,
    target_effect::Float64;
    alpha::Float64 = 0.025
)::Float64

    if info_frac >= 1.0 || info_frac <= 0.0 || se <= 0
        return NaN
    end

    # Information at interim
    info_interim = 1 / se^2

    # Projected information at final
    info_final = info_interim / info_frac

    # Information remaining
    info_remaining = info_final - info_interim

    # Z-statistic at interim
    z_interim = effect / se

    # Z-value required for final significance
    z_crit = sqrt(2) * erfinv(1 - 2*alpha)

    # Under the alternative (true effect = target_effect):
    # Final Z ~ N(target_effect * sqrt(info_final), 1)
    # Conditional on interim data, remaining contribution is:
    # target_effect * sqrt(info_remaining) + noise

    # Predicted final z-statistic (weighted combination)
    # z_final = (z_interim * sqrt(info_interim) + z_remaining * sqrt(info_remaining)) / sqrt(info_final)

    # Expected z_remaining under alternative
    expected_z_remaining = target_effect * sqrt(info_remaining)

    # Predicted final z
    predicted_z_final = (z_interim * sqrt(info_interim) + expected_z_remaining) / sqrt(info_final)

    # Variance of remaining contribution (under alternative)
    var_remaining = 1.0

    # Conditional power
    # P(Z_final > z_crit | interim data, theta = target_effect)
    # = P(Z_remaining > required_remaining)
    # where required_remaining = (z_crit * sqrt(info_final) - z_interim * sqrt(info_interim)) / sqrt(info_remaining)

    required_remaining = (z_crit * sqrt(info_final) - z_interim * sqrt(info_interim)) / sqrt(info_remaining)

    # Under H1: Z_remaining ~ N(target_effect * sqrt(info_remaining), 1)
    mean_remaining = target_effect * sqrt(info_remaining)

    # P(Z_remaining > required_remaining) = 1 - Φ((required_remaining - mean_remaining) / 1)
    cp = 1 - 0.5 * (1 + erf((required_remaining - mean_remaining) / sqrt(2)))

    return max(0.0, min(1.0, cp))
end

# ============================================================================
# Adaptive Trial Simulation
# ============================================================================

"""
    simulate_adaptive_trial(trial_spec::TrialSpec, design::AdaptiveDesign;
                            grid=nothing, solver=nothing)

Simulate an adaptive clinical trial with interim analyses and stopping rules.

# Arguments
- `trial_spec`: Trial specification
- `design`: AdaptiveDesign with interim analysis schedule and rules

# Returns
- `AdaptiveTrialResult` with complete results including interim analyses

# Example
```julia
base_design = ParallelDesign(2)
adaptive = AdaptiveDesign(base_design;
    interim_analyses = [0.5],  # One interim at 50%
    alpha_spending = :obrien_fleming
)
result = simulate_adaptive_trial(trial_spec, adaptive)
```
"""
function simulate_adaptive_trial(
    trial_spec::TrialSpec,
    design::AdaptiveDesign;
    grid = nothing,
    solver = nothing,
    target_effect::Float64 = 0.0,  # For conditional power (0 = use observed)
    futility_threshold::Float64 = 0.10
)::AdaptiveTrialResult

    rng = StableRNG(trial_spec.seed)

    # Validate inputs
    if isempty(trial_spec.arms)
        error("Trial must have at least one arm")
    end

    # Total planned sample size
    planned_n = sum(arm.n_subjects for arm in trial_spec.arms)

    # Interim analysis schedule
    interim_fractions = sort(design.interim_analyses)
    push!(interim_fractions, 1.0)  # Add final analysis

    # Get alpha level (default 0.025 one-sided)
    alpha = get(design.adaptation_rules, :alpha, 0.025)

    # Futility settings
    futility_enabled = get(design.adaptation_rules, :futility_enabled, true)

    # Generate full virtual population upfront
    population = generate_virtual_population(trial_spec.virtual_population, planned_n)

    # Cumulative data storage
    cumulative_results = Dict{String, Vector{Any}}()
    cumulative_values = Dict{String, Dict{Symbol, Vector{Float64}}}()

    for arm in trial_spec.arms
        cumulative_results[arm.name] = []
        cumulative_values[arm.name] = Dict{Symbol, Vector{Float64}}()
        for ep in trial_spec.endpoints
            cumulative_values[arm.name][ep.name] = Float64[]
        end
    end

    # Track interim results
    interim_results = InterimResult[]
    stopped_early = false
    stop_reason = nothing
    current_n = 0

    # Process each interim analysis point
    for (analysis_num, target_frac) in enumerate(interim_fractions)
        target_n = round(Int, target_frac * planned_n)
        n_to_enroll = target_n - current_n

        if n_to_enroll <= 0
            continue
        end

        # Enroll subjects up to this target
        for arm in trial_spec.arms
            arm_target = round(Int, arm.n_subjects * target_frac)
            arm_current = length(cumulative_results[arm.name])
            arm_n_to_enroll = arm_target - arm_current

            for i in 1:arm_n_to_enroll
                # Get subject from population
                subj_idx = arm_current + i
                if subj_idx > length(population)
                    break
                end
                subject = population[subj_idx]

                # Simulate individual
                result = simulate_individual(
                    arm, subject, trial_spec.duration_days,
                    trial_spec.pk_sampling_times;
                    rng = rng, grid = grid, solver = solver
                )

                push!(cumulative_results[arm.name], result)

                # Calculate endpoints
                for ep in trial_spec.endpoints
                    value = calculate_individual_endpoint(result, ep, subject)
                    push!(cumulative_values[arm.name][ep.name], value)
                end
            end
        end

        current_n = sum(length(cumulative_results[arm.name]) for arm in trial_spec.arms)

        # Skip interim analysis at final (just complete enrollment)
        if target_frac >= 1.0
            continue
        end

        # Perform interim analysis on primary endpoint
        if !isempty(trial_spec.endpoints)
            primary_ep = trial_spec.endpoints[1]

            # Get arm values
            arm_names = [arm.name for arm in trial_spec.arms]
            if length(arm_names) >= 2
                v1 = cumulative_values[arm_names[1]][primary_ep.name]
                v2 = cumulative_values[arm_names[2]][primary_ep.name]

                # Filter NaN
                v1_valid = filter(!isnan, v1)
                v2_valid = filter(!isnan, v2)

                if length(v1_valid) >= 2 && length(v2_valid) >= 2
                    # Compute effect and SE
                    n1, n2 = length(v1_valid), length(v2_valid)
                    mean1, mean2 = sum(v1_valid)/n1, sum(v2_valid)/n2
                    var1 = sum((v1_valid .- mean1).^2) / (n1 - 1)
                    var2 = sum((v2_valid .- mean2).^2) / (n2 - 1)

                    effect = mean2 - mean1
                    se = sqrt(var1/n1 + var2/n2)
                    z_stat = effect / se

                    # Two-tailed p-value
                    p_value = 2 * (1 - 0.5 * (1 + erf(abs(z_stat) / sqrt(2))))

                    # Get boundaries
                    info_frac = target_frac
                    efficacy_boundary = get_efficacy_boundary(design.alpha_spending, info_frac, alpha)

                    # Conditional power
                    cp_target = target_effect > 0 ? target_effect : effect
                    cp = compute_conditional_power(effect, se, info_frac, cp_target; alpha = alpha)

                    # Check stopping rules
                    stop_eff = check_efficacy_stop(abs(z_stat), efficacy_boundary)
                    stop_fut = futility_enabled && check_futility_stop(cp, futility_threshold)

                    # Record interim result
                    push!(interim_results, InterimResult(
                        analysis_num,
                        info_frac,
                        current_n,
                        effect,
                        se,
                        z_stat,
                        p_value,
                        efficacy_boundary,
                        futility_enabled ? futility_threshold : nothing,
                        cp,
                        stop_eff,
                        stop_fut
                    ))

                    # Execute stopping
                    if stop_eff
                        stopped_early = true
                        stop_reason = :efficacy
                        break
                    elseif stop_fut
                        stopped_early = true
                        stop_reason = :futility
                        break
                    end
                end
            end
        end
    end

    # Build final results
    arm_results = Dict{String, ArmResult}()
    for arm in trial_spec.arms
        n_enrolled = length(cumulative_results[arm.name])

        # Build summary stats
        summary_stats = Dict{Symbol, Dict{Symbol, Float64}}()
        for ep in trial_spec.endpoints
            values = cumulative_values[arm.name][ep.name]
            valid = filter(!isnan, values)
            if !isempty(valid)
                mean_val = sum(valid) / length(valid)
                var_val = length(valid) > 1 ? sum((valid .- mean_val).^2) / (length(valid) - 1) : 0.0
                summary_stats[ep.name] = Dict{Symbol, Float64}(
                    :mean => mean_val,
                    :sd => sqrt(var_val),
                    :n => Float64(length(valid))
                )
            end
        end

        arm_results[arm.name] = ArmResult(
            arm.name,
            n_enrolled,
            n_enrolled,  # Assume no dropout for now
            0,
            cumulative_results[arm.name],
            cumulative_values[arm.name],
            summary_stats
        )
    end

    # Final endpoint analysis
    final_analyses = Dict{Symbol, Dict{Symbol, Any}}()
    for ep in trial_spec.endpoints
        final_analyses[ep.name] = analyze_trial_endpoint(arm_results, ep, design.base_design)
    end

    return AdaptiveTrialResult(
        final_analyses,
        interim_results,
        stopped_early,
        stop_reason,
        current_n,
        planned_n,
        trial_spec.seed
    )
end

# ============================================================================
# Response-Adaptive Randomization (RAR)
# ============================================================================

"""
Abstract type for response-adaptive randomization specifications.
"""
abstract type RARSpec end

"""
    ThallWathenRAR

Thall-Wathen Bayesian response-adaptive randomization.

Uses posterior probabilities of treatment superiority to update allocation ratios.
Reference: Thall PF, Wathen JK (2007). Practical Bayesian adaptive randomisation.

# Fields
- `prior_alpha`: Beta prior alpha parameter (default 1.0)
- `prior_beta`: Beta prior beta parameter (default 1.0)
- `tuning_parameter`: Controls adaptation aggressiveness (default 0.5)
- `min_allocation`: Minimum allocation probability per arm (default 0.1)
- `burn_in_fraction`: Fraction of trial with fixed allocation (default 0.2)
"""
struct ThallWathenRAR <: RARSpec
    prior_alpha::Float64
    prior_beta::Float64
    tuning_parameter::Float64
    min_allocation::Float64
    burn_in_fraction::Float64

    function ThallWathenRAR(;
        prior_alpha::Float64 = 1.0,
        prior_beta::Float64 = 1.0,
        tuning_parameter::Float64 = 0.5,
        min_allocation::Float64 = 0.1,
        burn_in_fraction::Float64 = 0.2
    )
        @assert prior_alpha > 0 "prior_alpha must be positive"
        @assert prior_beta > 0 "prior_beta must be positive"
        @assert 0.0 <= tuning_parameter <= 1.0 "tuning_parameter must be in [0, 1]"
        @assert 0.0 < min_allocation < 0.5 "min_allocation must be in (0, 0.5)"
        @assert 0.0 <= burn_in_fraction < 1.0 "burn_in_fraction must be in [0, 1)"
        new(prior_alpha, prior_beta, tuning_parameter, min_allocation, burn_in_fraction)
    end
end

"""
    DoublyAdaptiveBiasedCoin

Doubly Adaptive Biased Coin Design (DBCD) for response-adaptive randomization.

Balances treatment allocation based on both response rates and allocation imbalance.
Reference: Hu F, Zhang LX (2004). Asymptotic properties of DBCD.

# Fields
- `gamma`: Allocation function parameter (default 2.0, higher = more aggressive)
- `target_ratio`: Target allocation ratio (default equal allocation)
- `min_allocation`: Minimum allocation probability per arm (default 0.1)
"""
struct DoublyAdaptiveBiasedCoin <: RARSpec
    gamma::Float64
    target_ratio::Vector{Float64}
    min_allocation::Float64

    function DoublyAdaptiveBiasedCoin(;
        gamma::Float64 = 2.0,
        target_ratio::Vector{Float64} = Float64[],
        min_allocation::Float64 = 0.1
    )
        @assert gamma >= 0.0 "gamma must be non-negative"
        @assert 0.0 < min_allocation < 0.5 "min_allocation must be in (0, 0.5)"
        new(gamma, target_ratio, min_allocation)
    end
end

"""
    compute_rar_probabilities(rar_spec::ThallWathenRAR, arm_data, n_total, n_planned)

Compute response-adaptive randomization probabilities using Thall-Wathen method.

# Arguments
- `rar_spec`: ThallWathenRAR specification
- `arm_data`: Vector of (n_subjects, n_responders) tuples for each arm
- `n_total`: Total subjects enrolled so far
- `n_planned`: Total planned sample size

# Returns
- Vector of allocation probabilities for each arm
"""
function compute_rar_probabilities(
    rar_spec::ThallWathenRAR,
    arm_data::Vector{Tuple{Int, Int}},
    n_total::Int,
    n_planned::Int
)::Vector{Float64}
    n_arms = length(arm_data)

    # During burn-in, use equal allocation
    if n_total < rar_spec.burn_in_fraction * n_planned
        return fill(1.0 / n_arms, n_arms)
    end

    # Compute posterior parameters for each arm
    posterior_alphas = Float64[]
    posterior_betas = Float64[]

    for (n_subj, n_resp) in arm_data
        push!(posterior_alphas, rar_spec.prior_alpha + n_resp)
        push!(posterior_betas, rar_spec.prior_beta + (n_subj - n_resp))
    end

    # Compute posterior means (expected response rates)
    posterior_means = [a / (a + b) for (a, b) in zip(posterior_alphas, posterior_betas)]

    # Apply tuning transformation: π_k ∝ (posterior_mean_k)^c
    c = rar_spec.tuning_parameter
    raw_probs = [p^c for p in posterior_means]

    # Normalize
    total_prob = sum(raw_probs)
    probs = raw_probs ./ total_prob

    # Apply minimum allocation constraint
    probs = clamp.(probs, rar_spec.min_allocation, 1.0 - rar_spec.min_allocation * (n_arms - 1))

    # Renormalize
    probs ./= sum(probs)

    return probs
end

"""
    compute_rar_probabilities(rar_spec::DoublyAdaptiveBiasedCoin, arm_data, n_total, n_planned)

Compute response-adaptive randomization probabilities using DBCD method.
"""
function compute_rar_probabilities(
    rar_spec::DoublyAdaptiveBiasedCoin,
    arm_data::Vector{Tuple{Int, Int}},
    n_total::Int,
    n_planned::Int
)::Vector{Float64}
    n_arms = length(arm_data)

    # Default to equal allocation if no target specified
    target = isempty(rar_spec.target_ratio) ? fill(1.0 / n_arms, n_arms) : rar_spec.target_ratio

    # Current allocation proportions
    n_subjects = [d[1] for d in arm_data]
    current_props = n_total > 0 ? n_subjects ./ n_total : fill(1.0 / n_arms, n_arms)

    # Response rates
    response_rates = [n > 0 ? r / n : 0.5 for (n, r) in arm_data]

    # DBCD allocation function
    # g(ρ, N) = ρ^γ / Σ_k ρ_k^γ where ρ is Neyman allocation based on response
    # Neyman allocation: ρ_k ∝ sqrt(p_k * (1 - p_k))
    neyman_allocation = [sqrt(p * (1 - p) + 1e-6) for p in response_rates]
    neyman_allocation ./= sum(neyman_allocation)

    # Apply gamma transformation
    gamma = rar_spec.gamma
    if gamma > 0
        transformed = neyman_allocation .^ gamma
        probs = transformed ./ sum(transformed)
    else
        probs = fill(1.0 / n_arms, n_arms)
    end

    # Apply minimum allocation constraint
    probs = clamp.(probs, rar_spec.min_allocation, 1.0 - rar_spec.min_allocation * (n_arms - 1))
    probs ./= sum(probs)

    return probs
end

# ============================================================================
# Sample Size Re-estimation (SSR)
# ============================================================================

"""
Abstract type for sample size re-estimation specifications.
"""
abstract type SSRSpec end

"""
    ConditionalPowerSSR

Sample size re-estimation based on conditional power.

Increases sample size to maintain target conditional power at interim analysis.
Reference: Mehta CR, Pocock SJ (2011). Adaptive increase in sample size.

# Fields
- `target_power`: Target conditional power (default 0.80)
- `max_increase_factor`: Maximum allowed sample size increase (default 2.0)
- `min_conditional_power`: Below this CP, allow re-estimation (default 0.50)
- `promising_zone`: CP range for re-estimation [lower, upper] (default [0.36, 0.80])
"""
struct ConditionalPowerSSR <: SSRSpec
    target_power::Float64
    max_increase_factor::Float64
    min_conditional_power::Float64
    promising_zone::Tuple{Float64, Float64}

    function ConditionalPowerSSR(;
        target_power::Float64 = 0.80,
        max_increase_factor::Float64 = 2.0,
        min_conditional_power::Float64 = 0.50,
        promising_zone::Tuple{Float64, Float64} = (0.36, 0.80)
    )
        @assert 0.5 <= target_power < 1.0 "target_power must be in [0.5, 1.0)"
        @assert max_increase_factor >= 1.0 "max_increase_factor must be >= 1.0"
        @assert 0.0 <= min_conditional_power < target_power "min_conditional_power must be < target_power"
        @assert promising_zone[1] < promising_zone[2] "promising_zone must be ordered"
        new(target_power, max_increase_factor, min_conditional_power, promising_zone)
    end
end

"""
    VarianceBasedSSR

Sample size re-estimation based on observed variance (blinded or unblinded).

Adjusts sample size when interim variance estimate differs from planned.
Reference: Gould AL (2001). Sample size re-estimation.

# Fields
- `variance_ratio_trigger`: Re-estimate if observed/planned variance ratio exceeds this (default 1.2)
- `max_increase_factor`: Maximum allowed sample size increase (default 1.5)
- `blinded`: Use blinded variance estimate (default true for regulatory acceptance)
"""
struct VarianceBasedSSR <: SSRSpec
    variance_ratio_trigger::Float64
    max_increase_factor::Float64
    blinded::Bool

    function VarianceBasedSSR(;
        variance_ratio_trigger::Float64 = 1.2,
        max_increase_factor::Float64 = 1.5,
        blinded::Bool = true
    )
        @assert variance_ratio_trigger > 1.0 "variance_ratio_trigger must be > 1.0"
        @assert max_increase_factor >= 1.0 "max_increase_factor must be >= 1.0"
        new(variance_ratio_trigger, max_increase_factor, blinded)
    end
end

"""
    compute_sample_size_reestimation(ssr_spec::ConditionalPowerSSR, interim_result, planned_n, effect_size)

Compute new sample size based on conditional power at interim.

# Arguments
- `ssr_spec`: ConditionalPowerSSR specification
- `interim_result`: InterimResult from interim analysis
- `planned_n`: Originally planned sample size
- `effect_size`: Target effect size for power calculations

# Returns
- (new_n::Int, reason::Symbol) - New sample size and reason for change
"""
function compute_sample_size_reestimation(
    ssr_spec::ConditionalPowerSSR,
    interim_result::InterimResult,
    planned_n::Int,
    effect_size::Float64;
    alpha::Float64 = 0.025
)::Tuple{Int, Symbol}
    cp = interim_result.conditional_power
    info_frac = interim_result.information_fraction

    # Check if in promising zone
    if cp < ssr_spec.promising_zone[1]
        # Below promising zone - don't increase (likely futile)
        return (planned_n, :futile)
    elseif cp > ssr_spec.promising_zone[2]
        # Above promising zone - on track, no increase needed
        return (planned_n, :on_track)
    end

    # In promising zone - compute sample size increase
    # Target: CP = target_power after increase

    # Current information
    se_current = interim_result.primary_se
    info_current = 1.0 / se_current^2

    # Information needed for target power
    z_alpha = sqrt(2) * erfinv(1 - 2*alpha)
    z_beta = sqrt(2) * erfinv(2*ssr_spec.target_power - 1)

    # Required total information for target power
    info_required = ((z_alpha + z_beta) / effect_size)^2

    # Additional information needed
    info_additional = max(0.0, info_required - info_current)

    # Convert to sample size (assuming equal allocation, 2-arm trial)
    # n_per_arm = 2 * sigma^2 / effect_size^2 for power
    # Approximate: n_new = n_planned * (info_required / info_total_planned)
    info_total_planned = info_current / info_frac
    increase_factor = info_required / info_total_planned

    # Apply constraints
    increase_factor = clamp(increase_factor, 1.0, ssr_spec.max_increase_factor)

    new_n = ceil(Int, planned_n * increase_factor)

    return (new_n, :promising)
end

"""
    compute_sample_size_reestimation(ssr_spec::VarianceBasedSSR, observed_variance, planned_variance, planned_n)

Compute new sample size based on observed variance.
"""
function compute_sample_size_reestimation(
    ssr_spec::VarianceBasedSSR,
    observed_variance::Float64,
    planned_variance::Float64,
    planned_n::Int
)::Tuple{Int, Symbol}
    variance_ratio = observed_variance / planned_variance

    if variance_ratio < ssr_spec.variance_ratio_trigger
        # Variance not inflated enough to trigger re-estimation
        return (planned_n, :no_change)
    end

    # Increase sample size proportionally to variance inflation
    increase_factor = min(variance_ratio, ssr_spec.max_increase_factor)
    new_n = ceil(Int, planned_n * increase_factor)

    return (new_n, :variance_inflated)
end

# ============================================================================
# Treatment Selection (Adaptive Arm Dropping)
# ============================================================================

"""
Abstract type for treatment selection specifications.
"""
abstract type TreatmentSelectionSpec end

"""
    DropLoserSelection

Drop inferior treatment arms based on interim efficacy assessment.

Reference: Stallard N, Todd S (2003). Sequential designs for Phase II/III trials.

# Fields
- `selection_criterion`: Criterion for dropping (:posterior_probability, :frequentist_pvalue, :effect_size)
- `threshold`: Threshold for dropping (interpretation depends on criterion)
- `min_arms`: Minimum number of arms to keep (default 2, control + best treatment)
- `selection_fraction`: Information fraction at which to select (default 0.5)
- `control_arm`: Index of control arm (never dropped, default 1)
"""
struct DropLoserSelection <: TreatmentSelectionSpec
    selection_criterion::Symbol
    threshold::Float64
    min_arms::Int
    selection_fraction::Float64
    control_arm::Int

    function DropLoserSelection(;
        selection_criterion::Symbol = :posterior_probability,
        threshold::Float64 = 0.10,
        min_arms::Int = 2,
        selection_fraction::Float64 = 0.5,
        control_arm::Int = 1
    )
        @assert selection_criterion in [:posterior_probability, :frequentist_pvalue, :effect_size]
        @assert threshold > 0 "threshold must be positive"
        @assert min_arms >= 2 "min_arms must be at least 2"
        @assert 0.0 < selection_fraction < 1.0 "selection_fraction must be in (0, 1)"
        new(selection_criterion, threshold, min_arms, selection_fraction, control_arm)
    end
end

"""
    TreatmentSelectionResult

Result of treatment selection analysis.
"""
struct TreatmentSelectionResult
    arms_to_continue::Vector{Int}
    arms_to_drop::Vector{Int}
    selection_scores::Dict{Int, Float64}
    reason::String
end

"""
    select_treatments(selection_spec::DropLoserSelection, arm_results, n_total, n_planned)

Select which treatment arms to continue based on interim results.

# Arguments
- `selection_spec`: DropLoserSelection specification
- `arm_results`: Dict mapping arm index to (mean_effect, se, n_subjects)
- `n_total`: Total subjects enrolled so far
- `n_planned`: Total planned sample size

# Returns
- TreatmentSelectionResult with arms to continue and drop
"""
function select_treatments(
    selection_spec::DropLoserSelection,
    arm_results::Dict{Int, Tuple{Float64, Float64, Int}},
    n_total::Int,
    n_planned::Int
)::TreatmentSelectionResult
    info_frac = n_total / n_planned

    # Check if at selection point
    if info_frac < selection_spec.selection_fraction - 0.05
        # Not yet at selection point - keep all arms
        all_arms = collect(keys(arm_results))
        return TreatmentSelectionResult(all_arms, Int[], Dict{Int, Float64}(), "Not at selection point")
    end

    n_arms = length(arm_results)
    control_idx = selection_spec.control_arm

    # Get control arm data
    if !haskey(arm_results, control_idx)
        control_idx = 1  # Default to first arm
    end
    control_mean, control_se, control_n = arm_results[control_idx]

    # Compute selection scores for each treatment arm
    selection_scores = Dict{Int, Float64}()

    for (arm_idx, (mean_eff, se, n)) in arm_results
        if arm_idx == control_idx
            selection_scores[arm_idx] = Inf  # Control never dropped
            continue
        end

        if selection_spec.selection_criterion == :posterior_probability
            # Bayesian posterior probability of being better than control
            # P(effect_arm > effect_control | data)
            effect_diff = mean_eff - control_mean
            se_diff = sqrt(se^2 + control_se^2)
            z = effect_diff / se_diff
            # P(diff > 0) using normal approximation
            prob_better = 0.5 * (1 + erf(z / sqrt(2)))
            selection_scores[arm_idx] = prob_better

        elseif selection_spec.selection_criterion == :frequentist_pvalue
            # P-value for superiority vs control
            effect_diff = mean_eff - control_mean
            se_diff = sqrt(se^2 + control_se^2)
            z = effect_diff / se_diff
            # Two-sided p-value
            p_value = 2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2))))
            selection_scores[arm_idx] = 1.0 - p_value  # Higher is better

        else  # :effect_size
            effect_diff = mean_eff - control_mean
            selection_scores[arm_idx] = effect_diff
        end
    end

    # Sort arms by score (descending, best first)
    sorted_arms = sort(collect(keys(selection_scores)), by=k -> -selection_scores[k])

    # Determine which arms to drop
    arms_to_continue = Int[]
    arms_to_drop = Int[]

    for (rank, arm_idx) in enumerate(sorted_arms)
        score = selection_scores[arm_idx]

        # Always keep control
        if arm_idx == control_idx
            push!(arms_to_continue, arm_idx)
            continue
        end

        # Keep if meets threshold and haven't reached min_arms
        meets_threshold = if selection_spec.selection_criterion == :posterior_probability
            score >= selection_spec.threshold
        elseif selection_spec.selection_criterion == :frequentist_pvalue
            score >= (1.0 - selection_spec.threshold)
        else
            score > 0  # Effect size > 0
        end

        if meets_threshold || length(arms_to_continue) < selection_spec.min_arms
            push!(arms_to_continue, arm_idx)
        else
            push!(arms_to_drop, arm_idx)
        end
    end

    reason = "Selected $(length(arms_to_continue)) arms at $(round(info_frac*100, digits=1))% information"

    return TreatmentSelectionResult(arms_to_continue, arms_to_drop, selection_scores, reason)
end

# ============================================================================
# Biomarker Enrichment Designs
# ============================================================================

"""
Abstract type for enrichment design specifications.
"""
abstract type EnrichmentSpec end

"""
    BiomarkerEnrichment

Biomarker-driven enrichment design that can adapt enrollment criteria based on interim results.

Reference: Simon R, Wang SJ (2006). Use of genomic signatures in therapeutics development.

# Fields
- `biomarker_name`: Name of the biomarker covariate
- `initial_threshold`: Initial threshold for biomarker-positive classification
- `adapt_threshold`: Whether to adapt threshold based on interim results
- `enrichment_fraction`: Information fraction at which to apply enrichment (default 0.5)
- `min_biomarker_positive_rate`: Minimum rate of biomarker-positive to trigger enrichment
- `treatment_effect_threshold`: Min treatment effect in biomarker+ to restrict enrollment
"""
struct BiomarkerEnrichment <: EnrichmentSpec
    biomarker_name::Symbol
    initial_threshold::Float64
    adapt_threshold::Bool
    enrichment_fraction::Float64
    min_biomarker_positive_rate::Float64
    treatment_effect_threshold::Float64

    function BiomarkerEnrichment(;
        biomarker_name::Symbol = :biomarker,
        initial_threshold::Float64 = 0.0,
        adapt_threshold::Bool = true,
        enrichment_fraction::Float64 = 0.5,
        min_biomarker_positive_rate::Float64 = 0.3,
        treatment_effect_threshold::Float64 = 0.0
    )
        @assert 0.0 < enrichment_fraction < 1.0 "enrichment_fraction must be in (0, 1)"
        @assert 0.0 < min_biomarker_positive_rate < 1.0 "min_biomarker_positive_rate must be in (0, 1)"
        new(biomarker_name, initial_threshold, adapt_threshold, enrichment_fraction,
            min_biomarker_positive_rate, treatment_effect_threshold)
    end
end

"""
    EnrichmentResult

Result of enrichment analysis.
"""
struct EnrichmentResult
    should_enrich::Bool
    new_threshold::Float64
    biomarker_positive_rate::Float64
    effect_overall::Float64
    effect_biomarker_positive::Float64
    effect_biomarker_negative::Float64
    reason::String
end

"""
    apply_enrichment(enrichment_spec::BiomarkerEnrichment, subject_data, arm_effects, n_total, n_planned)

Determine whether to apply biomarker enrichment based on interim results.

# Arguments
- `enrichment_spec`: BiomarkerEnrichment specification
- `subject_data`: Vector of subject biomarker values and outcomes
- `arm_effects`: Dict of treatment effects by arm for biomarker+/- subgroups
- `n_total`: Total subjects enrolled so far
- `n_planned`: Total planned sample size

# Returns
- EnrichmentResult with enrichment decision and analysis results
"""
function apply_enrichment(
    enrichment_spec::BiomarkerEnrichment,
    subject_data::Vector{Tuple{Float64, Float64, Int}},  # (biomarker, outcome, arm)
    n_total::Int,
    n_planned::Int
)::EnrichmentResult
    info_frac = n_total / n_planned

    # Default to no enrichment
    default_result = EnrichmentResult(
        false, enrichment_spec.initial_threshold,
        0.0, 0.0, 0.0, 0.0, "Not at enrichment point"
    )

    if info_frac < enrichment_spec.enrichment_fraction - 0.05
        return default_result
    end

    if isempty(subject_data)
        return default_result
    end

    # Classify subjects by biomarker status
    threshold = enrichment_spec.initial_threshold
    biomarker_positive = [(bm, outcome, arm) for (bm, outcome, arm) in subject_data if bm >= threshold]
    biomarker_negative = [(bm, outcome, arm) for (bm, outcome, arm) in subject_data if bm < threshold]

    n_positive = length(biomarker_positive)
    n_negative = length(biomarker_negative)
    n_subjects = n_positive + n_negative

    biomarker_positive_rate = n_positive / n_subjects

    # Compute effects in subgroups (assuming arm 1 = control, arm 2 = treatment)
    function compute_effect(data::Vector{Tuple{Float64, Float64, Int}})
        control = [outcome for (_, outcome, arm) in data if arm == 1]
        treatment = [outcome for (_, outcome, arm) in data if arm == 2]

        if isempty(control) || isempty(treatment)
            return 0.0, Inf
        end

        mean_control = sum(control) / length(control)
        mean_treatment = sum(treatment) / length(treatment)
        effect = mean_treatment - mean_control

        var_control = length(control) > 1 ? sum((control .- mean_control).^2) / (length(control) - 1) : 1.0
        var_treatment = length(treatment) > 1 ? sum((treatment .- mean_treatment).^2) / (length(treatment) - 1) : 1.0
        se = sqrt(var_control/length(control) + var_treatment/length(treatment))

        return effect, se
    end

    effect_overall, _ = compute_effect(subject_data)
    effect_positive, se_positive = n_positive >= 4 ? compute_effect(biomarker_positive) : (0.0, Inf)
    effect_negative, se_negative = n_negative >= 4 ? compute_effect(biomarker_negative) : (0.0, Inf)

    # Decision logic for enrichment
    should_enrich = false
    reason = ""
    new_threshold = threshold

    # Check if treatment effect is differential
    if effect_positive > enrichment_spec.treatment_effect_threshold &&
       effect_positive > 2.0 * effect_overall &&
       biomarker_positive_rate >= enrichment_spec.min_biomarker_positive_rate

        should_enrich = true
        reason = "Treatment effect concentrated in biomarker-positive subgroup"

        # Optionally adapt threshold to maximize effect
        if enrichment_spec.adapt_threshold
            # Simple approach: use current threshold
            # Advanced: could search for optimal threshold
            new_threshold = threshold
        end
    elseif effect_overall <= 0 && effect_positive > 0
        should_enrich = true
        reason = "Overall effect null but positive in biomarker-positive subgroup"
    else
        reason = "No evidence for differential treatment effect by biomarker status"
    end

    return EnrichmentResult(
        should_enrich,
        new_threshold,
        biomarker_positive_rate,
        effect_overall,
        effect_positive,
        effect_negative,
        reason
    )
end

# ============================================================================
# Enhanced Adaptive Trial Simulation with All Features
# ============================================================================

"""
    simulate_adaptive_trial_full(trial_spec::TrialSpec, design::AdaptiveDesign;
                                  rar_spec=nothing, ssr_spec=nothing,
                                  selection_spec=nothing, enrichment_spec=nothing, ...)

Full adaptive trial simulation supporting all adaptive design features.

# Arguments
- `trial_spec`: Trial specification
- `design`: AdaptiveDesign with interim analysis schedule
- `rar_spec`: Optional RARSpec for response-adaptive randomization
- `ssr_spec`: Optional SSRSpec for sample size re-estimation
- `selection_spec`: Optional TreatmentSelectionSpec for arm dropping
- `enrichment_spec`: Optional EnrichmentSpec for biomarker enrichment

# Returns
- Enhanced AdaptiveTrialResult with full adaptation history
"""
function simulate_adaptive_trial_full(
    trial_spec::TrialSpec,
    design::AdaptiveDesign;
    grid = nothing,
    solver = nothing,
    target_effect::Float64 = 0.0,
    futility_threshold::Float64 = 0.10,
    rar_spec::Union{Nothing, RARSpec} = nothing,
    ssr_spec::Union{Nothing, SSRSpec} = nothing,
    selection_spec::Union{Nothing, TreatmentSelectionSpec} = nothing,
    enrichment_spec::Union{Nothing, EnrichmentSpec} = nothing
)::AdaptiveTrialResult

    rng = StableRNG(trial_spec.seed)

    # Validate inputs
    if isempty(trial_spec.arms)
        error("Trial must have at least one arm")
    end

    # Total planned sample size (may be updated via SSR)
    planned_n = sum(arm.n_subjects for arm in trial_spec.arms)
    current_planned_n = planned_n

    # Track active arms (may be modified via treatment selection)
    active_arms = collect(1:length(trial_spec.arms))

    # Interim analysis schedule
    interim_fractions = sort(design.interim_analyses)
    push!(interim_fractions, 1.0)

    # Get alpha level
    alpha = get(design.adaptation_rules, :alpha, 0.025)
    futility_enabled = get(design.adaptation_rules, :futility_enabled, true)

    # Generate full virtual population upfront
    max_n = ceil(Int, planned_n * (ssr_spec !== nothing ? 2.0 : 1.0))  # Allow for SSR increase
    population = generate_virtual_population(trial_spec.virtual_population, max_n)

    # Allocation probabilities (may be updated via RAR)
    n_arms = length(trial_spec.arms)
    allocation_probs = fill(1.0 / n_arms, n_arms)

    # Cumulative data storage
    cumulative_results = Dict{String, Vector{Any}}()
    cumulative_values = Dict{String, Dict{Symbol, Vector{Float64}}}()
    arm_responders = Dict{Int, Int}()  # Track responders per arm for RAR

    for (i, arm) in enumerate(trial_spec.arms)
        cumulative_results[arm.name] = []
        cumulative_values[arm.name] = Dict{Symbol, Vector{Float64}}()
        for ep in trial_spec.endpoints
            cumulative_values[arm.name][ep.name] = Float64[]
        end
        arm_responders[i] = 0
    end

    # Track interim results
    interim_results = InterimResult[]
    stopped_early = false
    stop_reason = nothing
    current_n = 0
    enrichment_applied = false

    # Process each interim analysis point
    for (analysis_num, target_frac) in enumerate(interim_fractions)
        target_n = round(Int, target_frac * current_planned_n)
        n_to_enroll = target_n - current_n

        if n_to_enroll <= 0
            continue
        end

        # Enroll subjects up to this target with current allocation
        for subj_idx in 1:n_to_enroll
            pop_idx = current_n + subj_idx
            if pop_idx > length(population)
                break
            end
            subject = population[pop_idx]

            # Apply enrichment filter if active
            if enrichment_applied && enrichment_spec !== nothing
                biomarker_val = get(subject.covariates, enrichment_spec.biomarker_name, 0.0)
                if biomarker_val < enrichment_spec.initial_threshold
                    continue  # Skip biomarker-negative subjects
                end
            end

            # Determine arm assignment using allocation probabilities
            # Only assign to active arms
            active_probs = [i in active_arms ? allocation_probs[i] : 0.0 for i in 1:n_arms]
            active_probs ./= sum(active_probs)

            r = rand(rng)
            cumsum_prob = 0.0
            assigned_arm_idx = 1
            for (i, p) in enumerate(active_probs)
                cumsum_prob += p
                if r <= cumsum_prob
                    assigned_arm_idx = i
                    break
                end
            end

            arm = trial_spec.arms[assigned_arm_idx]

            # Simulate individual
            result = simulate_individual(
                arm, subject, trial_spec.duration_days,
                trial_spec.pk_sampling_times;
                rng = rng, grid = grid, solver = solver
            )

            push!(cumulative_results[arm.name], result)

            # Calculate endpoints
            for ep in trial_spec.endpoints
                value = calculate_individual_endpoint(result, ep, subject)
                push!(cumulative_values[arm.name][ep.name], value)

                # Track responders (for RAR)
                if ep.name == trial_spec.endpoints[1].name && value > 0
                    arm_responders[assigned_arm_idx] += 1
                end
            end
        end

        current_n = sum(length(cumulative_results[arm.name]) for arm in trial_spec.arms)

        # Skip detailed analysis at final
        if target_frac >= 1.0
            continue
        end

        # Update RAR probabilities if enabled
        if rar_spec !== nothing
            arm_data = [(length(cumulative_results[trial_spec.arms[i].name]),
                        arm_responders[i]) for i in 1:n_arms]
            allocation_probs = compute_rar_probabilities(rar_spec, arm_data, current_n, current_planned_n)
        end

        # Perform interim analysis on primary endpoint
        if !isempty(trial_spec.endpoints)
            primary_ep = trial_spec.endpoints[1]
            arm_names = [arm.name for arm in trial_spec.arms]

            if length(arm_names) >= 2
                v1 = cumulative_values[arm_names[1]][primary_ep.name]
                v2 = cumulative_values[arm_names[2]][primary_ep.name]

                v1_valid = filter(!isnan, v1)
                v2_valid = filter(!isnan, v2)

                if length(v1_valid) >= 2 && length(v2_valid) >= 2
                    n1, n2 = length(v1_valid), length(v2_valid)
                    mean1, mean2 = sum(v1_valid)/n1, sum(v2_valid)/n2
                    var1 = sum((v1_valid .- mean1).^2) / (n1 - 1)
                    var2 = sum((v2_valid .- mean2).^2) / (n2 - 1)

                    effect = mean2 - mean1
                    se = sqrt(var1/n1 + var2/n2)
                    z_stat = effect / se
                    p_value = 2 * (1 - 0.5 * (1 + erf(abs(z_stat) / sqrt(2))))

                    info_frac = target_frac
                    efficacy_boundary = get_efficacy_boundary(design.alpha_spending, info_frac, alpha)

                    cp_target = target_effect > 0 ? target_effect : effect
                    cp = compute_conditional_power(effect, se, info_frac, cp_target; alpha = alpha)

                    stop_eff = check_efficacy_stop(abs(z_stat), efficacy_boundary)
                    stop_fut = futility_enabled && check_futility_stop(cp, futility_threshold)

                    # Record interim result
                    interim_res = InterimResult(
                        analysis_num, info_frac, current_n, effect, se, z_stat, p_value,
                        efficacy_boundary, futility_enabled ? futility_threshold : nothing,
                        cp, stop_eff, stop_fut
                    )
                    push!(interim_results, interim_res)

                    # Sample size re-estimation
                    if ssr_spec !== nothing && !stop_eff && !stop_fut
                        new_n, ssr_reason = compute_sample_size_reestimation(
                            ssr_spec, interim_res, current_planned_n, cp_target; alpha=alpha
                        )
                        if new_n > current_planned_n
                            current_planned_n = new_n
                        end
                    end

                    # Treatment selection (arm dropping)
                    if selection_spec !== nothing && length(active_arms) > 2
                        arm_results_dict = Dict{Int, Tuple{Float64, Float64, Int}}()
                        for (i, arm) in enumerate(trial_spec.arms)
                            vals = filter(!isnan, cumulative_values[arm.name][primary_ep.name])
                            if !isempty(vals)
                                m = sum(vals) / length(vals)
                                v = length(vals) > 1 ? sum((vals .- m).^2) / (length(vals) - 1) : 1.0
                                arm_results_dict[i] = (m, sqrt(v / length(vals)), length(vals))
                            end
                        end

                        selection_result = select_treatments(selection_spec, arm_results_dict,
                                                            current_n, current_planned_n)
                        active_arms = selection_result.arms_to_continue
                    end

                    # Check stopping rules
                    if stop_eff
                        stopped_early = true
                        stop_reason = :efficacy
                        break
                    elseif stop_fut
                        stopped_early = true
                        stop_reason = :futility
                        break
                    end
                end
            end
        end

        # Enrichment check
        if enrichment_spec !== nothing && !enrichment_applied
            # Gather biomarker data
            subject_biomarker_data = Tuple{Float64, Float64, Int}[]
            for (arm_idx, arm) in enumerate(trial_spec.arms)
                for (j, result) in enumerate(cumulative_results[arm.name])
                    bm_val = 0.0  # Would need to track biomarker from subject
                    outcome_vals = cumulative_values[arm.name][trial_spec.endpoints[1].name]
                    if j <= length(outcome_vals)
                        push!(subject_biomarker_data, (bm_val, outcome_vals[j], arm_idx))
                    end
                end
            end

            enrichment_result = apply_enrichment(enrichment_spec, subject_biomarker_data,
                                                 current_n, current_planned_n)
            if enrichment_result.should_enrich
                enrichment_applied = true
            end
        end
    end

    # Build final results (same as original)
    arm_results = Dict{String, ArmResult}()
    for arm in trial_spec.arms
        n_enrolled = length(cumulative_results[arm.name])
        summary_stats = Dict{Symbol, Dict{Symbol, Float64}}()
        for ep in trial_spec.endpoints
            values = cumulative_values[arm.name][ep.name]
            valid = filter(!isnan, values)
            if !isempty(valid)
                mean_val = sum(valid) / length(valid)
                var_val = length(valid) > 1 ? sum((valid .- mean_val).^2) / (length(valid) - 1) : 0.0
                summary_stats[ep.name] = Dict{Symbol, Float64}(
                    :mean => mean_val, :sd => sqrt(var_val), :n => Float64(length(valid))
                )
            end
        end
        arm_results[arm.name] = ArmResult(arm.name, n_enrolled, n_enrolled, 0,
                                          cumulative_results[arm.name],
                                          cumulative_values[arm.name], summary_stats)
    end

    final_analyses = Dict{Symbol, Dict{Symbol, Any}}()
    for ep in trial_spec.endpoints
        final_analyses[ep.name] = analyze_trial_endpoint(arm_results, ep, design.base_design)
    end

    return AdaptiveTrialResult(final_analyses, interim_results, stopped_early, stop_reason,
                               current_n, planned_n, trial_spec.seed)
end
