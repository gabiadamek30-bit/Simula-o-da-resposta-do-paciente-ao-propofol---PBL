# Trial Endpoints
# Endpoint calculation and analysis

using SpecialFunctions: erf, erfinv

export calculate_endpoint, analyze_endpoint, compare_arms
export extract_pk_endpoint, extract_pd_endpoint, evaluate_safety_endpoint
export responder_analysis

"""
    extract_pk_endpoint(simulation_result, endpoint::PKEndpoint)

Extract PK endpoint value from a simulation result.

# Arguments
- `simulation_result`: Individual simulation result
- `endpoint::PKEndpoint`: PK endpoint specification

# Returns
- `Float64`: Endpoint value
"""
function extract_pk_endpoint(simulation_result, endpoint::PKEndpoint)
    t = simulation_result["t"]
    obs = simulation_result["observations"]

    # Find concentration data
    conc_key = "conc"
    for key in keys(obs)
        if contains(string(key), "conc") || contains(string(key), "Concentration")
            conc_key = key
            break
        end
    end

    if !haskey(obs, conc_key)
        return NaN
    end

    c = obs[conc_key]

    if endpoint.metric == :cmax
        return maximum(c)
    elseif endpoint.metric == :tmax
        idx = argmax(c)
        return t[idx]
    elseif endpoint.metric == :auc_0_t
        # Linear trapezoidal AUC
        auc = 0.0
        for i in 1:(length(t)-1)
            auc += (t[i+1] - t[i]) * (c[i] + c[i+1]) / 2
        end
        return auc
    elseif endpoint.metric == :auc_0_inf
        # AUC with extrapolation (simplified)
        auc_0_t = 0.0
        for i in 1:(length(t)-1)
            auc_0_t += (t[i+1] - t[i]) * (c[i] + c[i+1]) / 2
        end
        # Estimate lambda_z from terminal phase
        if length(t) >= 3 && c[end] > 0
            lambda_z = log(c[end-1] / c[end]) / (t[end] - t[end-1])
            if lambda_z > 0
                auc_extra = c[end] / lambda_z
                return auc_0_t + auc_extra
            end
        end
        return auc_0_t
    elseif endpoint.metric == :t_half
        # Estimate from terminal phase
        if length(t) >= 3 && c[end] > 0 && c[end-1] > c[end]
            lambda_z = log(c[end-1] / c[end]) / (t[end] - t[end-1])
            if lambda_z > 0
                return log(2) / lambda_z
            end
        end
        return NaN
    else
        return NaN
    end
end


"""
    extract_pd_endpoint(simulation_result, endpoint::PDEndpoint, baseline_value=nothing)

Extract PD endpoint value from a simulation result.

# Arguments
- `simulation_result`: Individual simulation result
- `endpoint::PDEndpoint`: PD endpoint specification
- `baseline_value`: Optional baseline value

# Returns
- `Float64`: Endpoint value
"""
function extract_pd_endpoint(simulation_result, endpoint::PDEndpoint,
                              baseline_value::Union{Nothing, Float64} = nothing)
    t = simulation_result["t"]
    obs = simulation_result["observations"]

    # Find effect data
    effect_key = "effect"
    for key in keys(obs)
        if contains(string(key), "effect") || contains(string(key), "response")
            effect_key = key
            break
        end
    end

    if !haskey(obs, effect_key)
        return NaN
    end

    e = obs[effect_key]

    if endpoint.metric == :emax
        return maximum(e)
    elseif endpoint.metric == :change_from_baseline
        baseline = baseline_value !== nothing ? baseline_value : e[1]
        # Find value at last assessment time
        last_time = maximum(endpoint.assessment_times)
        idx = argmin(abs.(t .- last_time))
        return e[idx] - baseline
    elseif endpoint.metric == :percent_change
        baseline = baseline_value !== nothing ? baseline_value : e[1]
        if baseline == 0
            return NaN
        end
        last_time = maximum(endpoint.assessment_times)
        idx = argmin(abs.(t .- last_time))
        return 100 * (e[idx] - baseline) / baseline
    elseif endpoint.metric == :time_to_response
        if endpoint.responder_threshold !== nothing
            for i in eachindex(t)
                if e[i] >= endpoint.responder_threshold
                    return t[i]
                end
            end
        end
        return NaN  # No response
    else
        return NaN
    end
end


"""
    evaluate_safety_endpoint(adverse_events, endpoint::SafetyEndpoint, n_subjects::Int)

Evaluate a safety endpoint.

# Arguments
- `adverse_events::Vector{Tuple{Int, Float64, Float64}}`: (subject_id, time, severity)
- `endpoint::SafetyEndpoint`: Safety endpoint specification
- `n_subjects::Int`: Total number of subjects

# Returns
- `Dict{Symbol, Any}`: Safety endpoint analysis
"""
function evaluate_safety_endpoint(adverse_events::Vector{Tuple{Int, Float64, Float64}},
                                   endpoint::SafetyEndpoint,
                                   n_subjects::Int)
    # Count events by severity
    severity_map = Dict(1.0 => :mild, 2.0 => :moderate, 3.0 => :severe)

    event_counts = Dict{Symbol, Int}()
    subjects_with_event = Set{Int}()

    for (subj_id, _, severity) in adverse_events
        sev_symbol = get(severity_map, severity, :unknown)
        if sev_symbol in endpoint.severity_levels
            event_counts[sev_symbol] = get(event_counts, sev_symbol, 0) + 1
            push!(subjects_with_event, subj_id)
        end
    end

    total_events = sum(values(event_counts); init=0)
    subjects_with_any = length(subjects_with_event)

    result = Dict{Symbol, Any}(
        :total_events => total_events,
        :subjects_with_event => subjects_with_any,
        :event_rate => subjects_with_any / n_subjects,
        :event_counts => event_counts
    )

    # Check threshold
    if endpoint.threshold !== nothing
        result[:exceeds_threshold] = subjects_with_any / n_subjects > endpoint.threshold
    end

    return result
end


"""
    calculate_endpoint(results::Vector, endpoint::EndpointSpec; kwargs...)

Calculate endpoint values for multiple subjects.

# Arguments
- `results::Vector`: Individual simulation results
- `endpoint::EndpointSpec`: Endpoint specification

# Keyword Arguments
- `baseline_values::Vector{Float64}`: Baseline values for PD endpoints

# Returns
- `Vector{Float64}`: Endpoint values
"""
function calculate_endpoint(results::Vector, endpoint::PKEndpoint;
                            baseline_values::Union{Nothing, Vector{Float64}} = nothing)
    return [extract_pk_endpoint(r, endpoint) for r in results]
end

function calculate_endpoint(results::Vector, endpoint::PDEndpoint;
                            baseline_values::Union{Nothing, Vector{Float64}} = nothing)
    if baseline_values !== nothing
        return [extract_pd_endpoint(r, endpoint, baseline_values[i])
                for (i, r) in enumerate(results)]
    else
        return [extract_pd_endpoint(r, endpoint) for r in results]
    end
end


"""
    analyze_endpoint(values::Vector{Float64}, endpoint::EndpointSpec)

Perform statistical analysis on endpoint values.

# Arguments
- `values::Vector{Float64}`: Endpoint values
- `endpoint::EndpointSpec`: Endpoint specification

# Returns
- `Dict{Symbol, Any}`: Analysis results
"""
function analyze_endpoint(values::Vector{Float64}, endpoint::EndpointSpec)
    # Filter NaN values
    valid_values = filter(!isnan, values)
    n = length(valid_values)

    if n == 0
        return Dict{Symbol, Any}(:n => 0, :mean => NaN, :sd => NaN)
    end

    mean_val = sum(valid_values) / n
    variance = n > 1 ? sum((valid_values .- mean_val).^2) / (n - 1) : 0.0
    sd_val = sqrt(variance)

    sorted = sort(valid_values)
    median_val = n % 2 == 1 ? sorted[(n+1)÷2] : (sorted[n÷2] + sorted[n÷2+1]) / 2

    result = Dict{Symbol, Any}(
        :n => n,
        :mean => mean_val,
        :sd => sd_val,
        :median => median_val,
        :min => minimum(valid_values),
        :max => maximum(valid_values),
        :q25 => sorted[max(1, round(Int, 0.25 * n))],
        :q75 => sorted[min(n, round(Int, 0.75 * n))]
    )

    # Log-transformed statistics for PK endpoints
    if endpoint isa PKEndpoint && endpoint.log_transform
        pos_values = filter(v -> v > 0, valid_values)
        if length(pos_values) > 0
            log_values = log.(pos_values)
            log_mean = sum(log_values) / length(pos_values)
            log_var = length(pos_values) > 1 ?
                sum((log_values .- log_mean).^2) / (length(pos_values) - 1) : 0.0

            result[:geometric_mean] = exp(log_mean)
            result[:geometric_cv] = sqrt(exp(log_var) - 1) * 100  # CV%
        end
    end

    return result
end


"""
    compare_arms(values1::Vector{Float64}, values2::Vector{Float64};
                 test::Symbol=:ttest, paired::Bool=false)

Compare endpoint values between two arms.

# Arguments
- `values1::Vector{Float64}`: Values from arm 1 (or period 1 for crossover)
- `values2::Vector{Float64}`: Values from arm 2 (or period 2 for crossover)

# Keyword Arguments
- `test::Symbol`: Statistical test (:ttest, :wilcoxon, :ratio)
- `paired::Bool`: Whether comparison is paired (crossover). When true, uses paired t-test.

# Returns
- `Dict{Symbol, Any}`: Comparison results including:
  - :test_type - :paired_ttest, :welch_ttest, or :ratio
  - :t_statistic, :df, :p_value - test statistics
  - :ci_lower, :ci_upper - 95% confidence interval
  - :difference or :geometric_mean_ratio - effect estimate

# Notes
For paired analysis, values1 and values2 must have equal length and be matched by subject.
"""
function compare_arms(values1::Vector{Float64}, values2::Vector{Float64};
                      test::Symbol = :ttest, paired::Bool = false)
    v1 = filter(!isnan, values1)
    v2 = filter(!isnan, values2)

    n1, n2 = length(v1), length(v2)

    if n1 < 2 || n2 < 2
        return Dict{Symbol, Any}(:error => "Insufficient data")
    end

    # For paired analysis, we need matched data
    if paired && n1 != n2
        return Dict{Symbol, Any}(:error => "Paired analysis requires equal sample sizes. Got $(n1) vs $(n2)")
    end

    mean1, mean2 = sum(v1)/n1, sum(v2)/n2
    var1 = sum((v1 .- mean1).^2) / (n1 - 1)
    var2 = sum((v2 .- mean2).^2) / (n2 - 1)

    result = Dict{Symbol, Any}(
        :n1 => n1, :n2 => n2,
        :mean1 => mean1, :mean2 => mean2,
        :sd1 => sqrt(var1), :sd2 => sqrt(var2),
        :difference => mean2 - mean1
    )

    if test == :ttest
        if paired
            # PAIRED T-TEST
            # Calculate paired differences
            differences = v2 .- v1
            n = length(differences)
            mean_diff = sum(differences) / n
            var_diff = sum((differences .- mean_diff).^2) / (n - 1)
            sd_diff = sqrt(var_diff)
            se_diff = sd_diff / sqrt(n)

            # T-statistic for paired test
            t_stat = mean_diff / se_diff
            df = n - 1

            # Two-tailed p-value using t-distribution approximation
            # Using the approximation: p ≈ 2 * (1 - Φ(|t| * √(df/(df-2)))) for large df
            # For small df, use more accurate approximation
            p_value = _approximate_t_pvalue(abs(t_stat), df)

            # 95% CI for mean difference
            t_crit = _trial_t_critical(df, 0.025)  # Two-tailed 95% CI
            ci_lower = mean_diff - t_crit * se_diff
            ci_upper = mean_diff + t_crit * se_diff

            result[:test_type] = :paired_ttest
            result[:t_statistic] = t_stat
            result[:df] = df
            result[:se] = se_diff
            result[:p_value] = p_value
            result[:ci_lower] = ci_lower
            result[:ci_upper] = ci_upper
            result[:mean_diff] = mean_diff
            result[:sd_diff] = sd_diff

            # Correlation between paired samples
            if var1 > 0 && var2 > 0
                covariance = sum((v1 .- mean1) .* (v2 .- mean2)) / (n - 1)
                correlation = covariance / (sqrt(var1) * sqrt(var2))
                result[:correlation] = correlation
            end
        else
            # UNPAIRED (WELCH'S) T-TEST
            se = sqrt(var1/n1 + var2/n2)
            t_stat = (mean2 - mean1) / se

            # Degrees of freedom (Welch-Satterthwaite)
            df = (var1/n1 + var2/n2)^2 /
                 ((var1/n1)^2/(n1-1) + (var2/n2)^2/(n2-1))

            # Two-tailed p-value
            p_value = _approximate_t_pvalue(abs(t_stat), df)

            # 95% CI for difference
            t_crit = _trial_t_critical(df, 0.025)
            ci_lower = (mean2 - mean1) - t_crit * se
            ci_upper = (mean2 - mean1) + t_crit * se

            result[:test_type] = :welch_ttest
            result[:t_statistic] = t_stat
            result[:df] = df
            result[:se] = se
            result[:p_value] = p_value
            result[:ci_lower] = ci_lower
            result[:ci_upper] = ci_upper
        end

    elseif test == :ratio
        # Geometric mean ratio (for log-normal data)
        pos_v1 = filter(v -> v > 0, v1)
        pos_v2 = filter(v -> v > 0, v2)

        if length(pos_v1) >= 2 && length(pos_v2) >= 2
            if paired && length(pos_v1) == length(pos_v2)
                # Paired ratio analysis (log difference)
                log_ratios = log.(pos_v2) .- log.(pos_v1)
                n = length(log_ratios)
                mean_log_ratio = sum(log_ratios) / n
                var_log_ratio = sum((log_ratios .- mean_log_ratio).^2) / (n - 1)
                se_log = sqrt(var_log_ratio / n)

                gmr = exp(mean_log_ratio)
                t_crit = _trial_t_critical(n - 1, 0.05)  # 90% CI for BE

                result[:test_type] = :paired_ratio
                result[:geometric_mean_ratio] = gmr
                result[:ci_lower_ratio] = exp(mean_log_ratio - t_crit * se_log)
                result[:ci_upper_ratio] = exp(mean_log_ratio + t_crit * se_log)
            else
                # Unpaired ratio analysis
                log_mean1 = sum(log.(pos_v1)) / length(pos_v1)
                log_mean2 = sum(log.(pos_v2)) / length(pos_v2)

                log_var1 = sum((log.(pos_v1) .- log_mean1).^2) / (length(pos_v1) - 1)
                log_var2 = sum((log.(pos_v2) .- log_mean2).^2) / (length(pos_v2) - 1)

                gmr = exp(log_mean2 - log_mean1)
                se_log = sqrt(log_var1/length(pos_v1) + log_var2/length(pos_v2))

                result[:test_type] = :unpaired_ratio
                result[:geometric_mean_ratio] = gmr
                result[:ci_lower_ratio] = exp(log(gmr) - 1.96 * se_log)
                result[:ci_upper_ratio] = exp(log(gmr) + 1.96 * se_log)
            end
        end
    end

    return result
end

"""
Approximate t-distribution critical value using Wilson-Hilferty transformation.
For a given df and alpha (one-tailed), returns the critical t-value.
Used internally by trial endpoint analysis (separate from NCA bioequivalence).
"""
function _trial_t_critical(df::Float64, alpha::Float64)
    # For large df, use normal approximation
    if df > 100
        # Standard normal critical values
        if alpha ≈ 0.025
            return 1.96
        elseif alpha ≈ 0.05
            return 1.645
        else
            # Use inverse normal approximation
            return sqrt(2) * erfinv(1 - 2*alpha)
        end
    end

    # For smaller df, use approximation
    # Based on Gleason (1988) approximation
    z = sqrt(2) * erfinv(1 - 2*alpha)  # Normal quantile
    g1 = (z^3 + z) / 4
    g2 = (5*z^5 + 16*z^3 + 3*z) / 96
    g3 = (3*z^7 + 19*z^5 + 17*z^3 - 15*z) / 384

    return z + g1/df + g2/df^2 + g3/df^3
end

function _trial_t_critical(df::Int, alpha::Float64)
    return _trial_t_critical(Float64(df), alpha)
end

"""
Approximate two-tailed p-value for t-distribution.
"""
function _approximate_t_pvalue(t_abs::Float64, df::Int)
    return _approximate_t_pvalue(t_abs, Float64(df))
end

function _approximate_t_pvalue(t_abs::Float64, df::Float64)
    # For large df, use normal approximation
    if df > 100
        # Using standard normal CDF approximation
        z = t_abs
        # Approximation of 2*(1 - Φ(z))
        return 2 * (1 - 0.5 * (1 + erf(z / sqrt(2))))
    end

    # For smaller df, use approximation based on incomplete beta function
    # p = 2 * (1 - F_t(|t|, df)) where F_t is t-distribution CDF
    # Using the relationship with incomplete beta: F_t = 1 - 0.5*I_x(df/2, 1/2) where x = df/(df + t^2)
    x = df / (df + t_abs^2)

    # Approximation of incomplete beta function for this case
    # For t >> 1, p ≈ 2 * df^(df/2) * |t|^(-df) / B(df/2, 1/2)
    if t_abs > 3 && df > 2
        # Use asymptotic approximation
        log_p = (df/2) * log(df / (df + t_abs^2)) + 0.5 * log(df / (df + t_abs^2))
        return 2 * exp(log_p) / sqrt(π * df)
    else
        # Use normal approximation with finite df correction
        z = t_abs * sqrt((df - 2) / df)
        return 2 * (1 - 0.5 * (1 + erf(z / sqrt(2))))
    end
end


"""
    responder_analysis(values::Vector{Float64}, threshold::Float64;
                       direction::Symbol=:greater)

Perform responder analysis.

# Arguments
- `values::Vector{Float64}`: Endpoint values
- `threshold::Float64`: Responder threshold

# Keyword Arguments
- `direction::Symbol`: :greater or :less

# Returns
- `Dict{Symbol, Any}`: Responder analysis results
"""
function responder_analysis(values::Vector{Float64}, threshold::Float64;
                            direction::Symbol = :greater)
    valid = filter(!isnan, values)
    n = length(valid)

    if n == 0
        return Dict{Symbol, Any}(:n => 0, :n_responders => 0, :response_rate => NaN)
    end

    n_responders = if direction == :greater
        count(v -> v >= threshold, valid)
    else
        count(v -> v <= threshold, valid)
    end

    rate = n_responders / n

    # 95% CI for proportion (Wilson score interval)
    z = 1.96
    denominator = 1 + z^2 / n
    center = (rate + z^2 / (2*n)) / denominator
    margin = z * sqrt(rate * (1 - rate) / n + z^2 / (4*n^2)) / denominator

    return Dict{Symbol, Any}(
        :n => n,
        :n_responders => n_responders,
        :response_rate => rate,
        :ci_lower => max(0, center - margin),
        :ci_upper => min(1, center + margin)
    )
end
