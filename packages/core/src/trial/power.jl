# Power Analysis
# Power calculation, sample size estimation, and alpha spending

using SpecialFunctions: erf, erfinv

export estimate_power, estimate_sample_size, alpha_spending_function
export estimate_power_analytical, incremental_alpha
export PowerResult, SampleSizeResult

using StableRNGs

"""
    PowerResult

Result of power analysis.

# Fields
- `power::Float64`: Estimated power
- `n_replicates::Int`: Number of simulation replicates
- `n_significant::Int`: Number of significant results
- `alpha::Float64`: Significance level used
- `effect_size::Float64`: Effect size used
- `ci_lower::Float64`: Lower 95% CI for power
- `ci_upper::Float64`: Upper 95% CI for power
"""
struct PowerResult
    power::Float64
    n_replicates::Int
    n_significant::Int
    alpha::Float64
    effect_size::Float64
    ci_lower::Float64
    ci_upper::Float64
end


"""
    SampleSizeResult

Result of sample size estimation.

# Fields
- `n_per_arm::Int`: Sample size per arm
- `total_n::Int`: Total sample size
- `target_power::Float64`: Target power
- `achieved_power::Float64`: Achieved power at recommended N
- `effect_size::Float64`: Effect size used
- `alpha::Float64`: Significance level
"""
struct SampleSizeResult
    n_per_arm::Int
    total_n::Int
    target_power::Float64
    achieved_power::Float64
    effect_size::Float64
    alpha::Float64
end


"""
    estimate_power(simulate_func, n_per_arm::Int, effect_size::Float64;
                   n_replicates::Int=1000, alpha::Float64=0.05, kwargs...)

Estimate power through simulation.

# Arguments
- `simulate_func`: Function that runs one trial and returns p-value
- `n_per_arm::Int`: Sample size per arm
- `effect_size::Float64`: True effect size

# Keyword Arguments
- `n_replicates::Int`: Number of simulation replicates
- `alpha::Float64`: Significance level
- Additional kwargs passed to simulate_func

# Returns
- `PowerResult`: Power analysis result

# Example
```julia
# Define a simulation function that returns a p-value
function sim_two_sample_ttest(n, effect; rng=nothing, sd=1.0)
    rng = isnothing(rng) ? Random.default_rng() : rng
    # Generate two groups
    group1 = randn(rng, n)
    group2 = randn(rng, n) .+ effect  # Shift by effect size
    # Calculate t-test p-value
    m1, m2 = mean(group1), mean(group2)
    s1, s2 = std(group1), std(group2)
    se = sqrt(s1^2/n + s2^2/n)
    t_stat = (m2 - m1) / se
    df = 2*n - 2
    # Two-sided p-value (approximation)
    p_value = 2 * (1 - cdf(TDist(df), abs(t_stat)))
    return p_value
end
result = estimate_power(sim_two_sample_ttest, 50, 0.5; n_replicates=1000)
```
"""
function estimate_power(simulate_func::Function, n_per_arm::Int, effect_size::Float64;
                        n_replicates::Int = 1000, alpha::Float64 = 0.05,
                        seed::UInt64 = UInt64(12345), kwargs...)

    rng = StableRNG(seed)
    n_significant = 0

    for _ in 1:n_replicates
        p_value = simulate_func(n_per_arm, effect_size; rng = rng, kwargs...)
        if p_value < alpha
            n_significant += 1
        end
    end

    power = n_significant / n_replicates

    # 95% CI for proportion
    z = 1.96
    se = sqrt(power * (1 - power) / n_replicates)
    ci_lower = max(0, power - z * se)
    ci_upper = min(1, power + z * se)

    return PowerResult(power, n_replicates, n_significant, alpha, effect_size,
                       ci_lower, ci_upper)
end


"""
    estimate_power_analytical(n_per_arm::Int, effect_size::Float64,
                              sd::Float64; alpha::Float64=0.05, n_arms::Int=2)

Estimate power analytically for a two-sample t-test.

# Arguments
- `n_per_arm::Int`: Sample size per arm
- `effect_size::Float64`: Expected difference in means
- `sd::Float64`: Standard deviation (assumed equal)

# Keyword Arguments
- `alpha::Float64`: Significance level
- `n_arms::Int`: Number of arms (for adjustment)

# Returns
- `Float64`: Power estimate
"""
function estimate_power_analytical(n_per_arm::Int, effect_size::Float64,
                                    sd::Float64; alpha::Float64 = 0.05, n_arms::Int = 2)
    # Two-sample t-test power calculation
    se = sd * sqrt(2 / n_per_arm)
    ncp = abs(effect_size) / se  # Non-centrality parameter

    # Critical value (approximate)
    z_alpha = if alpha == 0.05
        1.96
    elseif alpha == 0.025
        2.24
    elseif alpha == 0.01
        2.58
    else
        1.96
    end

    # Power = P(Z > z_alpha - ncp)
    # Approximate using normal distribution
    power = 0.5 * (1 + erf((ncp - z_alpha) / sqrt(2)))

    return power
end


"""
    estimate_sample_size(target_power::Float64, effect_size::Float64,
                         sd::Float64; alpha::Float64=0.05, n_arms::Int=2)

Estimate sample size needed to achieve target power.

# Arguments
- `target_power::Float64`: Target power (e.g., 0.80)
- `effect_size::Float64`: Expected difference in means
- `sd::Float64`: Standard deviation

# Keyword Arguments
- `alpha::Float64`: Significance level
- `n_arms::Int`: Number of arms

# Returns
- `SampleSizeResult`: Sample size estimation result

# Example
```julia
result = estimate_sample_size(0.80, 0.5, 1.0)
println("Need \$(result.n_per_arm) subjects per arm")
```
"""
function estimate_sample_size(target_power::Float64, effect_size::Float64,
                               sd::Float64; alpha::Float64 = 0.05, n_arms::Int = 2)

    # Binary search for sample size
    n_low = 4
    n_high = 10000

    achieved_power = 0.0
    best_n = n_high

    while n_high - n_low > 1
        n_mid = (n_low + n_high) ÷ 2
        power = estimate_power_analytical(n_mid, effect_size, sd; alpha = alpha)

        if power >= target_power
            best_n = n_mid
            achieved_power = power
            n_high = n_mid
        else
            n_low = n_mid
        end
    end

    # Final check at n_high
    power_high = estimate_power_analytical(n_high, effect_size, sd; alpha = alpha)
    if power_high >= target_power
        best_n = n_high
        achieved_power = power_high
    end

    return SampleSizeResult(best_n, best_n * n_arms, target_power, achieved_power,
                            effect_size, alpha)
end


"""
    estimate_sample_size_simulation(simulate_func, target_power::Float64,
                                    effect_size::Float64; kwargs...)

Estimate sample size through simulation.

# Arguments
- `simulate_func`: Function that simulates trial and returns p-value
- `target_power::Float64`: Target power
- `effect_size::Float64`: Effect size

# Keyword Arguments
- `n_range::Tuple{Int, Int}`: Range of sample sizes to search
- `n_replicates::Int`: Replicates per sample size evaluation
- `alpha::Float64`: Significance level

# Returns
- `SampleSizeResult`: Sample size estimation result
"""
function estimate_sample_size_simulation(simulate_func::Function,
                                          target_power::Float64,
                                          effect_size::Float64;
                                          n_range::Tuple{Int, Int} = (10, 500),
                                          n_replicates::Int = 500,
                                          alpha::Float64 = 0.05,
                                          seed::UInt64 = UInt64(12345))

    n_low, n_high = n_range
    best_n = n_high
    achieved_power = 0.0

    # Binary search
    while n_high - n_low > 5
        n_mid = (n_low + n_high) ÷ 2

        result = estimate_power(simulate_func, n_mid, effect_size;
                                n_replicates = n_replicates, alpha = alpha,
                                seed = seed)

        if result.power >= target_power
            best_n = n_mid
            achieved_power = result.power
            n_high = n_mid
        else
            n_low = n_mid
        end
    end

    # Fine-tune search
    for n in n_low:n_high
        result = estimate_power(simulate_func, n, effect_size;
                                n_replicates = n_replicates, alpha = alpha,
                                seed = seed)
        if result.power >= target_power
            best_n = n
            achieved_power = result.power
            break
        end
    end

    return SampleSizeResult(best_n, best_n * 2, target_power, achieved_power,
                            effect_size, alpha)
end


"""
    alpha_spending_function(information_fraction::Float64,
                            total_alpha::Float64,
                            spending_type::Symbol)

Calculate cumulative alpha spent at a given information fraction.

# Arguments
- `information_fraction::Float64`: Fraction of total information (0 to 1)
- `total_alpha::Float64`: Total type I error
- `spending_type::Symbol`: Spending function type

# Spending Types
- `:obrien_fleming`: O'Brien-Fleming spending (conservative early)
- `:pocock`: Pocock spending (linear)
- `:haybittle_peto`: Haybittle-Peto (very conservative early)
- `:linear`: Linear spending

# Returns
- `Float64`: Cumulative alpha spent

# Example
```julia
# Alpha spent at 50% interim analysis
alpha_spent = alpha_spending_function(0.5, 0.05, :obrien_fleming)
```
"""
function alpha_spending_function(information_fraction::Float64,
                                  total_alpha::Float64,
                                  spending_type::Symbol)

    t = clamp(information_fraction, 0.0, 1.0)

    if spending_type == :obrien_fleming
        # O'Brien-Fleming: α(t) = 2 - 2Φ(z_α/2 / √t)
        if t <= 0
            return 0.0
        end
        z_alpha2 = 1.96  # For α = 0.05 two-sided
        return 2 * (1 - 0.5 * (1 + erf(z_alpha2 / sqrt(t) / sqrt(2))))

    elseif spending_type == :pocock
        # Pocock: α(t) = α * log(1 + (e-1)*t)
        return total_alpha * log(1 + (exp(1) - 1) * t)

    elseif spending_type == :haybittle_peto
        # Haybittle-Peto: Very small alpha until final
        if t < 1.0
            return 0.001  # Nominal 0.001 at interim
        else
            return total_alpha
        end

    elseif spending_type == :linear
        # Linear spending
        return total_alpha * t

    else
        # Default to O'Brien-Fleming
        if t <= 0
            return 0.0
        end
        z_alpha2 = 1.96
        return 2 * (1 - 0.5 * (1 + erf(z_alpha2 / sqrt(t) / sqrt(2))))
    end
end


"""
    incremental_alpha(information_fractions::Vector{Float64},
                      total_alpha::Float64, spending_type::Symbol)

Calculate incremental alpha for each analysis.

# Arguments
- `information_fractions::Vector{Float64}`: Information fractions for each analysis
- `total_alpha::Float64`: Total type I error
- `spending_type::Symbol`: Spending function type

# Returns
- `Vector{Float64}`: Incremental alpha for each analysis
"""
function incremental_alpha(information_fractions::Vector{Float64},
                            total_alpha::Float64, spending_type::Symbol)

    # Add final analysis if not included
    fractions = copy(information_fractions)
    if fractions[end] < 1.0
        push!(fractions, 1.0)
    end

    cumulative = [alpha_spending_function(t, total_alpha, spending_type) for t in fractions]

    # Incremental = cumulative[i] - cumulative[i-1]
    incremental = zeros(length(cumulative))
    incremental[1] = cumulative[1]
    for i in 2:length(cumulative)
        incremental[i] = cumulative[i] - cumulative[i-1]
    end

    return incremental
end


"""
    futility_boundary(information_fraction::Float64, target_power::Float64;
                      boundary_type::Symbol=:conditional_power)

Calculate futility boundary value.

# Arguments
- `information_fraction::Float64`: Current information fraction
- `target_power::Float64`: Target power for the trial

# Keyword Arguments
- `boundary_type::Symbol`: Type of futility boundary
  - `:conditional_power`: Based on conditional power
  - `:predictive_probability`: Based on predictive probability

# Returns
- `Float64`: Futility boundary (stop if below this value)
"""
function futility_boundary(information_fraction::Float64, target_power::Float64;
                            boundary_type::Symbol = :conditional_power)

    if boundary_type == :conditional_power
        # Typically stop if conditional power < 0.10 or 0.20
        # Return z-score boundary
        remaining_info = 1 - information_fraction
        if remaining_info <= 0
            return -Inf
        end

        # Approximate: need z-score that gives conditional power = futility threshold
        futility_threshold = 0.10
        # This is a simplified approximation
        return -0.5 * sqrt(information_fraction / remaining_info)

    elseif boundary_type == :predictive_probability
        # Stop if predictive probability of success is very low
        return -1.0 * sqrt(information_fraction)

    else
        return -1.0
    end
end


"""
    conditional_power(z_current::Float64, information_fraction::Float64,
                      z_final::Float64; effect_size::Float64=0.0)

Calculate conditional power given current results.

# Arguments
- `z_current::Float64`: Current z-statistic
- `information_fraction::Float64`: Current information fraction
- `z_final::Float64`: Critical value for final analysis

# Keyword Arguments
- `effect_size::Float64`: Assumed effect size (0 = current estimate)

# Returns
- `Float64`: Conditional power
"""
function conditional_power(z_current::Float64, information_fraction::Float64,
                            z_final::Float64; effect_size::Float64 = 0.0)

    remaining = 1 - information_fraction

    if remaining <= 0
        return z_current >= z_final ? 1.0 : 0.0
    end

    # Under current trend (effect_size = 0), project z-statistic
    projected_z = z_current + effect_size * sqrt(remaining)

    # Conditional power = P(Z_final > z_final | Z_current)
    mean_z_final = z_current * sqrt(information_fraction) +
                   effect_size * sqrt(remaining)
    sd_z_final = sqrt(remaining)

    # P(Z > z_final)
    cp = 0.5 * (1 + erf((mean_z_final - z_final) / (sd_z_final * sqrt(2))))

    return clamp(cp, 0.0, 1.0)
end
