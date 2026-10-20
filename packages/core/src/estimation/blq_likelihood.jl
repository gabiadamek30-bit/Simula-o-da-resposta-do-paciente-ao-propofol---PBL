# BLQ/Censoring Likelihood Functions
# Implements M1, M2, M3 methods following FDA/NONMEM conventions
#
# References:
# - FDA Guidance: Population PK Analysis (2022)
# - NONMEM User Guide: LAPLACIAN with M3 method
# - Beal SL (2001). Ways to fit a PK model with some data below the quantification limit
# - Bergstrand M, Karlsson MO (2009). Handling data below the limit of quantification

using SpecialFunctions: erf

export observation_log_likelihood_blq, censored_log_likelihood, log_phi_stable
export compute_blq_summary, generate_blq_warnings, maximum_consecutive_true

"""
    observation_log_likelihood_blq(y, f, spec, blq_config, is_censored, lloq)

Compute -2LL contribution for a single observation, handling BLQ/censoring.

# Arguments
- `y::Float64`: Observed value
- `f::Float64`: Predicted value (IPRED)
- `spec::ResidualErrorSpec`: Residual error specification
- `blq_config::BLQConfig`: BLQ handling configuration
- `is_censored::Bool`: Whether this observation is BLQ
- `lloq::Float64`: Lower limit of quantification

# Returns
- `Float64`: -2LL contribution

# Methods
- M1 (Discard): Returns 0.0 (no contribution to likelihood)
- M2 (Impute): Uses standard Gaussian likelihood with imputed value
- M3 (Censored): Uses censored likelihood P(Y < LLOQ)
"""
function observation_log_likelihood_blq(
    y::Float64,
    f::Float64,
    spec::ResidualErrorSpec,
    blq_config::BLQConfig,
    is_censored::Bool,
    lloq::Float64
)::Float64

    if !is_censored
        # Normal observation - standard Gaussian likelihood
        return observation_log_likelihood(y, f, spec)
    end

    # Handle BLQ observation based on method
    if blq_config.method == BLQ_M1_DISCARD
        # M1: Return 0 contribution (effectively discards)
        return 0.0

    elseif blq_config.method == BLQ_M2_IMPUTE_HALF
        # M2a: Impute with LLOQ/2
        y_imputed = lloq / 2.0
        return observation_log_likelihood(y_imputed, f, spec)

    elseif blq_config.method == BLQ_M2_IMPUTE_ZERO
        # M2b: Impute with 0
        return observation_log_likelihood(0.0, f, spec)

    elseif blq_config.method == BLQ_M3_LIKELIHOOD
        # M3: Censored likelihood - P(Y < LLOQ)
        return censored_log_likelihood(f, lloq, spec)
    end

    # Fallback to standard likelihood (should not reach here)
    return observation_log_likelihood(y, f, spec)
end

"""
    censored_log_likelihood(ipred, lloq, spec) -> Float64

Compute -2LL for a left-censored observation using M3 method.

The M3 method treats BLQ observations as censored data:
    L = P(Y < LLOQ | IPRED, sigma)
      = Phi((LLOQ - IPRED) / sigma)

    -2LL = -2 * log(Phi((LLOQ - IPRED) / sigma))

# Arguments
- `ipred::Float64`: Individual prediction
- `lloq::Float64`: Lower limit of quantification
- `spec::ResidualErrorSpec`: Residual error specification

# Returns
- `Float64`: -2LL contribution for censored observation

Uses numerically stable log_phi for small probabilities.
"""
function censored_log_likelihood(
    ipred::Float64,
    lloq::Float64,
    spec::ResidualErrorSpec
)::Float64
    # Get residual SD from error model
    var_res = residual_variance(ipred, spec)

    if var_res <= 0.0 || !isfinite(var_res)
        return 1e10  # Large penalty for invalid variance
    end

    sigma = sqrt(var_res)

    # Standardized distance: z = (LLOQ - IPRED) / sigma
    z = (lloq - ipred) / sigma

    # log(Phi(z)) with numerical stability
    log_phi = log_phi_stable(z)

    # Return -2LL contribution
    return -2.0 * log_phi
end

"""
    log_phi_stable(z) -> Float64

Compute log(Phi(z)) with numerical stability for extreme z values.

Phi(z) is the CDF of the standard normal distribution.
This function handles extreme values of z to avoid underflow/overflow.

# Regions:
- z > 6: Phi(z) ≈ 1, so log(Phi(z)) ≈ 0
- -38 < z ≤ 6: Standard computation is numerically stable
- z ≤ -38: Use asymptotic expansion to avoid underflow

The asymptotic expansion for very negative z:
    log(Phi(z)) ≈ -z²/2 - log(-z) - log(sqrt(2π)) + O(1/z²)
"""
function log_phi_stable(z::Float64)::Float64
    if z > 6.0
        # Phi(z) ≈ 1, log(Phi(z)) ≈ 0
        # More precise: -exp(-z²/2)/(z*sqrt(2π)) but effectively 0
        return 0.0
    elseif z > -38.0
        # Standard computation is stable in this range
        # Phi(z) = 0.5 * (1 + erf(z / sqrt(2)))
        phi = 0.5 * (1.0 + erf(z / sqrt(2.0)))
        return log(max(phi, 1e-300))
    else
        # Very negative z: use asymptotic expansion to avoid underflow
        # For z → -∞: log(Phi(z)) ≈ -z²/2 - log(-z) - 0.5*log(2π)
        # This comes from: Phi(z) ≈ phi(z)/(-z) for large negative z
        # where phi(z) = exp(-z²/2)/sqrt(2π) is the PDF
        return -0.5 * z^2 - log(-z) - 0.5 * log(2π)
    end
end

"""
    compute_blq_summary(observed, blq_config) -> BLQSummary

Compute summary statistics for BLQ observations in the dataset.

# Arguments
- `observed::ObservedData`: Observed data with subjects
- `blq_config::BLQConfig`: BLQ configuration

# Returns
- `BLQSummary`: Summary of BLQ observations including counts, percentages, and warnings
"""
function compute_blq_summary(observed::ObservedData, blq_config::BLQConfig)::BLQSummary
    total_obs = 0
    total_blq = 0
    blq_by_subject = Dict{String,Int}()

    for subject in observed.subjects
        n_obs = length(subject.observations)
        total_obs += n_obs

        # Count BLQ observations for this subject
        n_blq = if !isempty(subject.blq_flags)
            sum(subject.blq_flags)
        elseif blq_config.lloq > 0.0
            # Use global LLOQ if no per-subject flags
            sum(subject.observations .< blq_config.lloq)
        else
            0
        end

        total_blq += n_blq
        blq_by_subject[subject.subject_id] = n_blq
    end

    blq_pct = total_obs > 0 ? 100.0 * total_blq / total_obs : 0.0

    # Generate warnings
    warnings = generate_blq_warnings(observed, blq_config, total_obs, total_blq, blq_pct)

    return BLQSummary(
        total_obs,
        total_blq,
        blq_pct,
        blq_by_subject,
        blq_config.method,
        warnings
    )
end

"""
    generate_blq_warnings(observed, blq_config, total_obs, total_blq, blq_pct) -> Vector{String}

Generate warnings about BLQ handling based on the data characteristics.

# Warnings generated:
- High BLQ percentage (>30%)
- Consecutive BLQ observations exceeding threshold
- Subjects with >50% BLQ
- Using M1/M2 with high BLQ percentage
"""
function generate_blq_warnings(
    observed::ObservedData,
    blq_config::BLQConfig,
    total_obs::Int,
    total_blq::Int,
    blq_pct::Float64
)::Vector{String}
    warnings = String[]

    # Warning for high overall BLQ percentage
    if blq_pct > 30.0
        if blq_config.method == BLQ_M1_DISCARD
            push!(warnings, "High BLQ percentage ($(round(blq_pct, digits=1))%) with M1 method. " *
                           "Consider using M3 (censored likelihood) to reduce bias.")
        elseif blq_config.method == BLQ_M2_IMPUTE_HALF || blq_config.method == BLQ_M2_IMPUTE_ZERO
            push!(warnings, "High BLQ percentage ($(round(blq_pct, digits=1))%) with M2 method. " *
                           "Consider using M3 (censored likelihood) for more accurate estimates.")
        else
            push!(warnings, "High BLQ percentage ($(round(blq_pct, digits=1))%). " *
                           "Estimates may be sensitive to BLQ handling assumptions.")
        end
    end

    # Check for consecutive BLQ observations
    for subject in observed.subjects
        blq_flags = if !isempty(subject.blq_flags)
            subject.blq_flags
        elseif blq_config.lloq > 0.0
            subject.observations .< blq_config.lloq
        else
            Bool[]
        end

        if !isempty(blq_flags)
            max_consecutive = maximum_consecutive_true(blq_flags)
            if max_consecutive > blq_config.max_consecutive_blq
                push!(warnings, "Subject $(subject.subject_id) has $max_consecutive consecutive BLQ " *
                               "observations (threshold: $(blq_config.max_consecutive_blq)).")
            end

            # Warning for subjects with >50% BLQ
            subject_blq_pct = 100.0 * sum(blq_flags) / length(blq_flags)
            if subject_blq_pct > 50.0
                push!(warnings, "Subject $(subject.subject_id) has $(round(subject_blq_pct, digits=1))% " *
                               "BLQ observations.")
            end
        end
    end

    return warnings
end

"""
    maximum_consecutive_true(flags::Vector{Bool}) -> Int

Find the maximum number of consecutive true values in a boolean vector.

# Example
```julia
flags = [false, true, true, true, false, true]
maximum_consecutive_true(flags)  # Returns 3
```
"""
function maximum_consecutive_true(flags::Vector{Bool})::Int
    if isempty(flags)
        return 0
    end

    max_count = 0
    current_count = 0

    for flag in flags
        if flag
            current_count += 1
            max_count = max(max_count, current_count)
        else
            current_count = 0
        end
    end

    return max_count
end

"""
    get_blq_flags_for_subject(subject, blq_config) -> Vector{Bool}

Get BLQ flags for a subject, using either stored flags or computing from LLOQ.

# Arguments
- `subject::SubjectData`: Subject data
- `blq_config::BLQConfig`: BLQ configuration

# Returns
- `Vector{Bool}`: BLQ flags for each observation
"""
function get_blq_flags_for_subject(subject::SubjectData, blq_config::Union{Nothing,BLQConfig})::Vector{Bool}
    if blq_config === nothing
        return Bool[]
    end

    # Use stored BLQ flags if available
    if !isempty(subject.blq_flags)
        return subject.blq_flags
    end

    # Otherwise compute from LLOQ
    lloq = subject.lloq > 0.0 ? subject.lloq : blq_config.lloq
    if lloq > 0.0
        return subject.observations .< lloq
    end

    return falses(length(subject.observations))
end

"""
    get_lloq_for_subject(subject, blq_config) -> Float64

Get the LLOQ for a subject, using subject-specific or global value.

# Arguments
- `subject::SubjectData`: Subject data
- `blq_config::BLQConfig`: BLQ configuration

# Returns
- `Float64`: LLOQ value to use
"""
function get_lloq_for_subject(subject::SubjectData, blq_config::Union{Nothing,BLQConfig})::Float64
    if blq_config === nothing
        return 0.0
    end

    # Use subject-specific LLOQ if available, otherwise use global
    return subject.lloq > 0.0 ? subject.lloq : blq_config.lloq
end

export get_blq_flags_for_subject, get_lloq_for_subject
