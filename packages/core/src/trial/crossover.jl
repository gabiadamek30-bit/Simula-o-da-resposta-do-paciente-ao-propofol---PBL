# Crossover Trial Analysis
# Specialized analysis for crossover designs including bioequivalence

using SpecialFunctions: erf

export analyze_crossover, CrossoverAnalysis, CrossoverBEResult
export test_period_effect, test_sequence_effect, assess_bioequivalence
export compute_within_subject_cv

"""
    CrossoverAnalysis

Results from crossover trial analysis.

# Fields
- `treatment_effect`: Paired comparison result between treatments
- `period_effect`: Test for period effect
- `sequence_effect`: Test for carryover/sequence effect
- `within_subject_cv`: Within-subject coefficient of variation
- `be_assessment`: Bioequivalence assessment (if applicable)
"""
struct CrossoverAnalysis
    treatment_effect::Dict{Symbol, Any}
    period_effect::Dict{Symbol, Any}
    sequence_effect::Dict{Symbol, Any}
    within_subject_cv::Float64
    be_assessment::Union{Nothing, Dict{Symbol, Any}}
end

"""
    CrossoverBEResult

Results from bioequivalence assessment in crossover trials using TOST.

# Fields
- `geometric_mean_ratio`: Point estimate of GMR
- `ci_90`: 90% confidence interval for GMR
- `is_bioequivalent`: Whether BE criteria are met
- `theta_limits`: BE acceptance limits used

Note: This is the trial module's BE result type. For NCA BE analysis,
see `BioequivalenceResult` in the NCA module.
"""
struct CrossoverBEResult
    geometric_mean_ratio::Float64
    ci_90::Tuple{Float64, Float64}
    is_bioequivalent::Bool
    theta_limits::Tuple{Float64, Float64}
end

"""
    analyze_crossover(arm_results, endpoint, design::CrossoverDesign)

Perform complete crossover analysis for an endpoint.

# Arguments
- `arm_results`: Dict of ArmResult by arm name
- `endpoint`: EndpointSpec for the parameter to analyze
- `design`: CrossoverDesign specification

# Returns
- `CrossoverAnalysis` with treatment effect, period effect, sequence effect, and BE assessment

# Example
```julia
result = analyze_crossover(arm_results, PKEndpoint(:auc_0_inf), design)
```
"""
function analyze_crossover(
    arm_results::Dict{String, ArmResult},
    endpoint::EndpointSpec,
    design::CrossoverDesign
)::CrossoverAnalysis

    # Extract data by arm
    arm_names = collect(keys(arm_results))
    if length(arm_names) != 2
        error("Crossover analysis requires exactly 2 arms, got $(length(arm_names))")
    end

    ref_name, test_name = arm_names[1], arm_names[2]

    ref_values = if haskey(arm_results[ref_name].endpoint_values, endpoint.name)
        arm_results[ref_name].endpoint_values[endpoint.name]
    else
        Float64[]
    end

    test_values = if haskey(arm_results[test_name].endpoint_values, endpoint.name)
        arm_results[test_name].endpoint_values[endpoint.name]
    else
        Float64[]
    end

    # Filter NaN values (must maintain pairing)
    valid_pairs = [(r, t) for (r, t) in zip(ref_values, test_values)
                   if !isnan(r) && !isnan(t)]

    if length(valid_pairs) < 3
        return CrossoverAnalysis(
            Dict{Symbol, Any}(:error => "Insufficient paired data"),
            Dict{Symbol, Any}(),
            Dict{Symbol, Any}(),
            NaN,
            nothing
        )
    end

    ref_valid = [p[1] for p in valid_pairs]
    test_valid = [p[2] for p in valid_pairs]

    # Primary analysis: paired comparison
    treatment_effect = if endpoint isa PKEndpoint && endpoint.log_transform
        # Use log-transformed analysis for PK
        compare_arms(ref_valid, test_valid; test = :ratio, paired = true)
    else
        # Use difference-based analysis
        compare_arms(ref_valid, test_valid; test = :ttest, paired = true)
    end

    # Period effect analysis
    period_effect = test_period_effect(ref_valid, test_valid, design)

    # Sequence effect analysis (potential carryover)
    sequence_effect = test_sequence_effect(ref_valid, test_valid, design)

    # Within-subject CV
    within_cv = compute_within_subject_cv(ref_valid, test_valid)

    # Bioequivalence assessment (for PK endpoints with log_transform)
    be_assessment = if endpoint isa PKEndpoint && endpoint.log_transform
        assess_bioequivalence(ref_valid, test_valid)
    else
        nothing
    end

    return CrossoverAnalysis(
        treatment_effect,
        period_effect,
        sequence_effect,
        within_cv,
        be_assessment
    )
end

"""
    test_period_effect(values_period1, values_period2, design)

Test for period effect in crossover design.

The period effect tests whether there's a systematic difference between periods
independent of treatment. A significant period effect suggests environmental
or learning effects that could bias the treatment comparison.

# Returns
- Dict with :effect_size, :p_value, :significant
"""
function test_period_effect(
    values_period1::Vector{Float64},
    values_period2::Vector{Float64},
    design::CrossoverDesign
)::Dict{Symbol, Any}

    n = length(values_period1)
    if n != length(values_period2) || n < 3
        return Dict{Symbol, Any}(:error => "Insufficient data for period effect test")
    end

    # Period effect: average of sums within subject
    # For 2x2 crossover: period effect = mean((Y_ij1 + Y_ij2)/2) across subjects
    # where j indicates treatment
    sums = values_period1 .+ values_period2
    mean_sum = sum(sums) / n

    # Half the sum represents the period-average response
    period_avg = mean_sum / 2

    # Calculate difference (simplified period effect test)
    # In a 2x2 crossover, period effect estimated from sequence-by-period interaction
    # Simplified: test if mean(period1) ≠ mean(period2) using paired test
    differences = values_period2 .- values_period1
    mean_diff = sum(differences) / n
    var_diff = sum((differences .- mean_diff).^2) / (n - 1)
    se_diff = sqrt(var_diff / n)

    # Period effect is half the mean difference (for equal sequences)
    period_effect_est = mean_diff / 2

    result = Dict{Symbol, Any}(
        :period_effect_estimate => period_effect_est,
        :mean_period1 => sum(values_period1) / n,
        :mean_period2 => sum(values_period2) / n,
        :n_subjects => n
    )

    # Significance test
    if se_diff > 0
        t_stat = mean_diff / se_diff / 2  # Period effect is half the sequence difference
        df = n - 1
        p_value = _approximate_t_pvalue(abs(t_stat), Float64(df))
        result[:t_statistic] = t_stat
        result[:p_value] = p_value
        result[:significant] = p_value < 0.05
    else
        result[:significant] = false
    end

    return result
end

"""
    test_sequence_effect(values_period1, values_period2, design)

Test for sequence (carryover) effect in crossover design.

A significant sequence effect suggests carryover from treatment in period 1
to response in period 2. If detected, the crossover analysis may be biased.

# Returns
- Dict with :effect_size, :p_value, :significant
"""
function test_sequence_effect(
    values_period1::Vector{Float64},
    values_period2::Vector{Float64},
    design::CrossoverDesign
)::Dict{Symbol, Any}

    n = length(values_period1)
    if n != length(values_period2) || n < 4
        return Dict{Symbol, Any}(:error => "Insufficient data for sequence effect test")
    end

    # For a 2x2 crossover with sequences RT and TR:
    # Sequence effect tested by comparing subject sums between sequences
    # Without actual sequence assignment, we use a simplified test

    # Sum for each subject (should be similar across sequences if no carryover)
    subject_sums = values_period1 .+ values_period2

    mean_sum = sum(subject_sums) / n
    var_sum = sum((subject_sums .- mean_sum).^2) / (n - 1)
    sd_sum = sqrt(var_sum)

    # CV of sums (high CV suggests sequence effect or high variability)
    cv_sum = mean_sum > 0 ? sd_sum / mean_sum * 100 : 0.0

    result = Dict{Symbol, Any}(
        :mean_subject_sum => mean_sum,
        :sd_subject_sum => sd_sum,
        :cv_sum_percent => cv_sum,
        :n_subjects => n
    )

    # For proper sequence effect, we'd need sequence assignments
    # Without them, flag potential issue if CV is very high
    if design.n_sequences == 2 && n >= 6
        # Use first half vs second half as proxy for sequences
        # (This is a placeholder - real implementation needs sequence data)
        n_half = n ÷ 2
        seq1_sums = subject_sums[1:n_half]
        seq2_sums = subject_sums[(n_half+1):end]

        if length(seq1_sums) >= 2 && length(seq2_sums) >= 2
            comparison = compare_arms(seq1_sums, seq2_sums; test = :ttest, paired = false)
            if haskey(comparison, :p_value)
                result[:sequence_comparison] = comparison
                result[:p_value] = comparison[:p_value]
                result[:significant] = comparison[:p_value] < 0.10  # Conservative threshold
            end
        end
    end

    if !haskey(result, :significant)
        result[:significant] = false
        result[:note] = "Full sequence effect test requires sequence assignment data"
    end

    return result
end

"""
    compute_within_subject_cv(values1, values2)

Compute within-subject coefficient of variation from paired data.

This is the intra-subject CV, which represents the variability of
repeated measurements within the same subject.

# Returns
- Within-subject CV as a percentage
"""
function compute_within_subject_cv(
    values1::Vector{Float64},
    values2::Vector{Float64}
)::Float64

    n = length(values1)
    if n != length(values2) || n < 2
        return NaN
    end

    # For log-normal data (typical for PK), use log-transformed approach
    # Filter positive values
    valid_pairs = [(v1, v2) for (v1, v2) in zip(values1, values2)
                   if v1 > 0 && v2 > 0]

    if length(valid_pairs) < 2
        return NaN
    end

    log_v1 = [log(p[1]) for p in valid_pairs]
    log_v2 = [log(p[2]) for p in valid_pairs]

    # Within-subject variance from log differences
    log_diff = log_v2 .- log_v1
    n_valid = length(log_diff)
    mean_log_diff = sum(log_diff) / n_valid

    # Residual variance (within-subject)
    ss_residual = sum((log_diff .- mean_log_diff).^2)
    mse = ss_residual / (n_valid - 1)

    # Within-subject CV = sqrt(exp(MSE) - 1) * 100
    cv_within = sqrt(exp(mse) - 1) * 100

    return cv_within
end

"""
    assess_bioequivalence(reference_values, test_values;
                           theta1=0.80, theta2=1.25, alpha=0.05)

Perform bioequivalence assessment using Two One-Sided Tests (TOST).

The TOST procedure tests whether the 90% CI for the geometric mean ratio
(Test/Reference) lies entirely within the acceptance limits [theta1, theta2].

# Arguments
- `reference_values`: PK parameter values for reference formulation
- `test_values`: PK parameter values for test formulation
- `theta1`: Lower BE limit (default 0.80)
- `theta2`: Upper BE limit (default 1.25)
- `alpha`: Significance level for each one-sided test (default 0.05)

# Returns
- Dict with :geometric_mean_ratio, :ci_90, :is_bioequivalent, :conclusion
"""
function assess_bioequivalence(
    reference_values::Vector{Float64},
    test_values::Vector{Float64};
    theta1::Float64 = 0.80,
    theta2::Float64 = 1.25,
    alpha::Float64 = 0.05
)::Dict{Symbol, Any}

    # Filter positive values (required for log transformation)
    valid_pairs = [(r, t) for (r, t) in zip(reference_values, test_values)
                   if r > 0 && t > 0 && !isnan(r) && !isnan(t)]

    n = length(valid_pairs)

    if n < 3
        return Dict{Symbol, Any}(:error => "Insufficient data for BE assessment")
    end

    log_ref = [log(p[1]) for p in valid_pairs]
    log_test = [log(p[2]) for p in valid_pairs]

    # Calculate log ratio (within-subject)
    log_ratios = log_test .- log_ref

    mean_log_ratio = sum(log_ratios) / n
    var_log_ratio = sum((log_ratios .- mean_log_ratio).^2) / (n - 1)
    se_log = sqrt(var_log_ratio / n)

    # Geometric mean ratio
    gmr = exp(mean_log_ratio)

    # 90% CI (using t-distribution with n-1 df)
    t_crit = _t_critical(n - 1, alpha)  # One-sided alpha for 90% CI
    ci_lower = exp(mean_log_ratio - t_crit * se_log)
    ci_upper = exp(mean_log_ratio + t_crit * se_log)

    # BE conclusion: 90% CI must be entirely within [theta1, theta2]
    is_be = (ci_lower >= theta1) && (ci_upper <= theta2)

    # TOST p-values (for completeness)
    # H01: mu_T - mu_R < log(theta1)  -> reject if t1 > t_crit
    # H02: mu_T - mu_R > log(theta2)  -> reject if t2 < -t_crit
    t1 = (mean_log_ratio - log(theta1)) / se_log
    t2 = (mean_log_ratio - log(theta2)) / se_log

    # P-values for one-sided tests
    p_value_lower = 1 - 0.5 * (1 + erf(t1 * sqrt((n-2)/n) / sqrt(2)))
    p_value_upper = 0.5 * (1 + erf(t2 * sqrt((n-2)/n) / sqrt(2)))

    # Overall TOST p-value is the maximum of the two
    p_value_tost = max(p_value_lower, p_value_upper)

    return Dict{Symbol, Any}(
        :geometric_mean_ratio => gmr,
        :ci_90_lower => ci_lower,
        :ci_90_upper => ci_upper,
        :is_bioequivalent => is_be,
        :theta_limits => (theta1, theta2),
        :n_subjects => n,
        :intra_subject_cv => sqrt(exp(var_log_ratio) - 1) * 100,
        :point_estimate_percent => gmr * 100,
        :p_value_tost => p_value_tost,
        :conclusion => is_be ?
            "Bioequivalence demonstrated (90% CI within $(theta1*100)%-$(theta2*100)%)" :
            "Bioequivalence NOT demonstrated (90% CI: $(round(ci_lower*100, digits=1))%-$(round(ci_upper*100, digits=1))%)"
    )
end

"""
    analyze_crossover_endpoint(arm_results, endpoint, design)

Wrapper to analyze a single endpoint in crossover design context.
Returns analysis dict compatible with standard endpoint analysis format.
"""
function analyze_crossover_endpoint(
    arm_results::Dict{String, ArmResult},
    endpoint::EndpointSpec,
    design::CrossoverDesign
)::Dict{Symbol, Any}

    crossover_result = analyze_crossover(arm_results, endpoint, design)

    result = Dict{Symbol, Any}(
        :crossover_analysis => crossover_result,
        :treatment_comparison => crossover_result.treatment_effect,
        :is_paired => true,
        :test_type => get(crossover_result.treatment_effect, :test_type, :paired_ttest),
        :within_subject_cv => crossover_result.within_subject_cv
    )

    # Add BE results if available
    if crossover_result.be_assessment !== nothing
        result[:be_assessment] = crossover_result.be_assessment
    end

    # Add warnings for significant period/sequence effects
    warnings = String[]
    if get(crossover_result.period_effect, :significant, false)
        push!(warnings, "Significant period effect detected (p < 0.05)")
    end
    if get(crossover_result.sequence_effect, :significant, false)
        push!(warnings, "Potential sequence/carryover effect detected")
    end
    if !isempty(warnings)
        result[:warnings] = warnings
    end

    return result
end
