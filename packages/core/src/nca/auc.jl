# AUC Calculations
# FDA/EMA compliant Area Under the Curve calculations

export auc_0_t, auc_0_inf, auc_0_tau, aumc_0_t, aumc_0_inf
export auc_interval, auc_linear, auc_log_linear, auc_partial

# =============================================================================
# Core AUC Functions
# =============================================================================

"""
    auc_0_t(t, c, config)

Calculate AUC from time 0 to the last measurable concentration (AUC0-t).

Uses the method specified in config (linear, log-linear, or lin-log mixed).

# Arguments
- `t::Vector{Float64}`: Time points (sorted, ascending)
- `c::Vector{Float64}`: Concentration values
- `config::NCAConfig`: NCA configuration

# Returns
- `Float64`: AUC from 0 to tlast
"""
function auc_0_t(t::Vector{Float64}, c::Vector{Float64}, config::NCAConfig)
    n = length(t)
    @assert n == length(c) "Time and concentration vectors must have same length"
    @assert n >= 2 "Need at least 2 points for AUC calculation"

    auc = 0.0

    for i in 2:n
        # Skip if both concentrations are zero or negative
        if c[i-1] <= 0.0 && c[i] <= 0.0
            continue
        end

        dt = t[i] - t[i-1]

        if dt <= 0.0
            continue
        end

        auc += _auc_interval(t[i-1], t[i], c[i-1], c[i], config.method)
    end

    return auc
end

"""
    auc_0_inf(t, c, lambda_z, clast, config)

Calculate AUC from time 0 to infinity (AUC0-∞).

AUC0-∞ = AUC0-t + Clast/λz

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `lambda_z::Float64`: Terminal elimination rate constant
- `clast::Float64`: Last measurable concentration
- `config::NCAConfig`: NCA configuration

# Returns
- `Tuple{Float64, Float64}`: (AUC0-inf, percent_extrapolated)
"""
function auc_0_inf(
    t::Vector{Float64},
    c::Vector{Float64},
    lambda_z::Float64,
    clast::Float64,
    config::NCAConfig
)
    @assert lambda_z > 0.0 "lambda_z must be positive"
    @assert clast >= 0.0 "clast must be non-negative"

    auc_0t = auc_0_t(t, c, config)
    auc_extra = clast / lambda_z
    auc_inf = auc_0t + auc_extra

    pct_extra = (auc_extra / auc_inf) * 100.0

    return (auc_inf, pct_extra)
end

"""
    auc_0_tau(t, c, tau, config)

Calculate AUC over a dosing interval (AUC0-τ) for steady-state analysis.

# Arguments
- `t::Vector{Float64}`: Time points within dosing interval
- `c::Vector{Float64}`: Concentration values
- `tau::Float64`: Dosing interval duration
- `config::NCAConfig`: NCA configuration

# Returns
- `Float64`: AUC over the dosing interval
"""
function auc_0_tau(
    t::Vector{Float64},
    c::Vector{Float64},
    tau::Float64,
    config::NCAConfig
)
    @assert tau > 0.0 "Dosing interval tau must be positive"

    # Filter points within the dosing interval
    valid_idx = findall(ti -> 0.0 <= ti <= tau, t)

    if length(valid_idx) < 2
        error("Insufficient points within dosing interval [0, $tau]")
    end

    t_tau = t[valid_idx]
    c_tau = c[valid_idx]

    return auc_0_t(t_tau, c_tau, config)
end

"""
    auc_interval(t1, t2, c1, c2, method)

Calculate AUC for a single interval using specified method.

# Arguments
- `t1, t2::Float64`: Start and end times
- `c1, c2::Float64`: Concentrations at t1 and t2
- `method::NCAMethod`: Calculation method

# Returns
- `Float64`: AUC for the interval
"""
function auc_interval(t1::Float64, t2::Float64, c1::Float64, c2::Float64, method::NCAMethod)
    return _auc_interval(t1, t2, c1, c2, method)
end

# =============================================================================
# Internal AUC Methods
# =============================================================================

function _auc_interval(t1::Float64, t2::Float64, c1::Float64, c2::Float64, ::LinearMethod)
    return _auc_linear(t1, t2, c1, c2)
end

function _auc_interval(t1::Float64, t2::Float64, c1::Float64, c2::Float64, ::LogLinearMethod)
    return _auc_log_linear(t1, t2, c1, c2)
end

function _auc_interval(t1::Float64, t2::Float64, c1::Float64, c2::Float64, ::LinLogMixedMethod)
    # Lin-up/Log-down: Use linear if ascending, log-linear if descending
    if c2 > c1
        # Ascending (absorption phase) - use linear
        return _auc_linear(t1, t2, c1, c2)
    elseif c2 < c1 && c1 > 0.0 && c2 > 0.0
        # Descending (elimination phase) with positive concentrations - use log-linear
        return _auc_log_linear(t1, t2, c1, c2)
    else
        # Equal or one is zero - use linear
        return _auc_linear(t1, t2, c1, c2)
    end
end

"""
    _auc_linear(t1, t2, c1, c2)

Linear trapezoidal rule: AUC = (t2 - t1) * (c1 + c2) / 2
"""
function _auc_linear(t1::Float64, t2::Float64, c1::Float64, c2::Float64)
    dt = t2 - t1
    return 0.5 * dt * (c1 + c2)
end

"""
    _auc_log_linear(t1, t2, c1, c2)

Log-linear trapezoidal rule: AUC = (t2 - t1) * (c1 - c2) / ln(c1/c2)

For exponential decay: C(t) = C1 * exp(-k*(t-t1))
AUC = ∫C dt = (C1 - C2) / k = (C1 - C2) * (t2 - t1) / ln(C1/C2)
"""
function _auc_log_linear(t1::Float64, t2::Float64, c1::Float64, c2::Float64)
    # Handle edge cases
    if c1 <= 0.0 || c2 <= 0.0
        return _auc_linear(t1, t2, max(c1, 0.0), max(c2, 0.0))
    end

    if abs(c1 - c2) < 1e-12 * max(c1, c2)
        # Concentrations nearly equal - use linear to avoid division by ~0
        return _auc_linear(t1, t2, c1, c2)
    end

    dt = t2 - t1
    log_ratio = log(c1 / c2)

    return dt * (c1 - c2) / log_ratio
end

# =============================================================================
# AUMC (Area Under First Moment Curve)
# =============================================================================

"""
    aumc_0_t(t, c, config)

Calculate AUMC from time 0 to the last measurable concentration.

AUMC = ∫ t * C(t) dt

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `config::NCAConfig`: NCA configuration

# Returns
- `Float64`: AUMC from 0 to tlast
"""
function aumc_0_t(t::Vector{Float64}, c::Vector{Float64}, config::NCAConfig)
    n = length(t)
    @assert n == length(c) "Time and concentration vectors must have same length"
    @assert n >= 2 "Need at least 2 points for AUMC calculation"

    aumc = 0.0

    for i in 2:n
        if c[i-1] <= 0.0 && c[i] <= 0.0
            continue
        end

        dt = t[i] - t[i-1]

        if dt <= 0.0
            continue
        end

        aumc += _aumc_interval(t[i-1], t[i], c[i-1], c[i], config.method)
    end

    return aumc
end

"""
    aumc_0_inf(t, c, lambda_z, clast, tlast, config)

Calculate AUMC from time 0 to infinity.

AUMC0-∞ = AUMC0-t + (tlast * Clast / λz) + (Clast / λz²)

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `lambda_z::Float64`: Terminal elimination rate constant
- `clast::Float64`: Last measurable concentration
- `tlast::Float64`: Time of last measurable concentration
- `config::NCAConfig`: NCA configuration

# Returns
- `Float64`: AUMC from 0 to infinity
"""
function aumc_0_inf(
    t::Vector{Float64},
    c::Vector{Float64},
    lambda_z::Float64,
    clast::Float64,
    tlast::Float64,
    config::NCAConfig
)
    @assert lambda_z > 0.0 "lambda_z must be positive"

    aumc_0t = aumc_0_t(t, c, config)

    # Extrapolation: AUMC_extra = (tlast * Clast / λz) + (Clast / λz²)
    aumc_extra = (tlast * clast / lambda_z) + (clast / (lambda_z^2))

    return aumc_0t + aumc_extra
end

# =============================================================================
# Internal AUMC Methods
# =============================================================================

function _aumc_interval(t1::Float64, t2::Float64, c1::Float64, c2::Float64, ::LinearMethod)
    return _aumc_linear(t1, t2, c1, c2)
end

function _aumc_interval(t1::Float64, t2::Float64, c1::Float64, c2::Float64, ::LogLinearMethod)
    return _aumc_log_linear(t1, t2, c1, c2)
end

function _aumc_interval(t1::Float64, t2::Float64, c1::Float64, c2::Float64, ::LinLogMixedMethod)
    if c2 > c1
        return _aumc_linear(t1, t2, c1, c2)
    elseif c2 < c1 && c1 > 0.0 && c2 > 0.0
        return _aumc_log_linear(t1, t2, c1, c2)
    else
        return _aumc_linear(t1, t2, c1, c2)
    end
end

"""
    _aumc_linear(t1, t2, c1, c2)

Linear trapezoidal AUMC:
AUMC = (t2 - t1) * (t1*c1 + t2*c2) / 2 + (t2 - t1)² * (c2 - c1) / 6
"""
function _aumc_linear(t1::Float64, t2::Float64, c1::Float64, c2::Float64)
    dt = t2 - t1
    return dt * (t1 * c1 + t2 * c2) / 2.0 + (dt^2) * (c2 - c1) / 6.0
end

"""
    _aumc_log_linear(t1, t2, c1, c2)

Log-linear trapezoidal AUMC for exponential decay.
"""
function _aumc_log_linear(t1::Float64, t2::Float64, c1::Float64, c2::Float64)
    if c1 <= 0.0 || c2 <= 0.0
        return _aumc_linear(t1, t2, max(c1, 0.0), max(c2, 0.0))
    end

    if abs(c1 - c2) < 1e-12 * max(c1, c2)
        return _aumc_linear(t1, t2, c1, c2)
    end

    dt = t2 - t1
    log_ratio = log(c1 / c2)

    # AUMC = dt * (t2*c2 - t1*c1) / ln(c1/c2) + dt² * (c1 - c2) / ln(c1/c2)²
    term1 = dt * (t2 * c2 - t1 * c1) / log_ratio
    term2 = (dt^2) * (c1 - c2) / (log_ratio^2)

    return term1 + term2
end

# =============================================================================
# Partial AUC
# =============================================================================

"""
    auc_partial(t, c, t_start, t_end, config)

Calculate partial AUC between specified time points.

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `t_start::Float64`: Start time for partial AUC
- `t_end::Float64`: End time for partial AUC
- `config::NCAConfig`: NCA configuration

# Returns
- `Float64`: Partial AUC between t_start and t_end
"""
function auc_partial(
    t::Vector{Float64},
    c::Vector{Float64},
    t_start::Float64,
    t_end::Float64,
    config::NCAConfig
)
    @assert t_start < t_end "t_start must be less than t_end"

    # Find indices within range
    idx_start = findfirst(ti -> ti >= t_start, t)
    idx_end = findlast(ti -> ti <= t_end, t)

    if idx_start === nothing || idx_end === nothing || idx_start > idx_end
        return 0.0
    end

    # Interpolate at boundaries if needed
    t_partial = Float64[]
    c_partial = Float64[]

    # Interpolate start point if necessary - respect config.method
    if t[idx_start] > t_start && idx_start > 1
        c_interp = _interpolate_concentration(
            t[idx_start-1], t[idx_start],
            c[idx_start-1], c[idx_start],
            t_start, config.method
        )
        push!(t_partial, t_start)
        push!(c_partial, c_interp)
    end

    # Add points within range
    for i in idx_start:idx_end
        push!(t_partial, t[i])
        push!(c_partial, c[i])
    end

    # Interpolate end point if necessary - respect config.method
    if t[idx_end] < t_end && idx_end < length(t)
        c_interp = _interpolate_concentration(
            t[idx_end], t[idx_end+1],
            c[idx_end], c[idx_end+1],
            t_end, config.method
        )
        push!(t_partial, t_end)
        push!(c_partial, c_interp)
    end

    return auc_0_t(t_partial, c_partial, config)
end

"""
    _interpolate_concentration(t1, t2, c1, c2, t_target, method)

Interpolate concentration at target time using specified method.

# Arguments
- `t1, t2::Float64`: Bounding time points
- `c1, c2::Float64`: Concentrations at t1 and t2
- `t_target::Float64`: Target time for interpolation
- `method::NCAMethod`: Interpolation method to use

# Returns
- `Float64`: Interpolated concentration at t_target
"""
function _interpolate_concentration(
    t1::Float64, t2::Float64,
    c1::Float64, c2::Float64,
    t_target::Float64,
    method::NCAMethod = LinLogMixedMethod()
)
    return _interpolate_concentration_dispatch(t1, t2, c1, c2, t_target, method)
end

# Linear method: always use linear interpolation
function _interpolate_concentration_dispatch(
    t1::Float64, t2::Float64,
    c1::Float64, c2::Float64,
    t_target::Float64,
    ::LinearMethod
)
    slope = (c2 - c1) / (t2 - t1)
    return c1 + slope * (t_target - t1)
end

# Log-linear method: use log-linear if both positive, else linear
function _interpolate_concentration_dispatch(
    t1::Float64, t2::Float64,
    c1::Float64, c2::Float64,
    t_target::Float64,
    ::LogLinearMethod
)
    if c1 > 0.0 && c2 > 0.0
        # Log-linear interpolation
        k = log(c1 / c2) / (t2 - t1)
        return c1 * exp(-k * (t_target - t1))
    else
        # Fall back to linear if not both positive
        slope = (c2 - c1) / (t2 - t1)
        return c1 + slope * (t_target - t1)
    end
end

# Lin-Log Mixed method: use log-linear for descending, linear for ascending
function _interpolate_concentration_dispatch(
    t1::Float64, t2::Float64,
    c1::Float64, c2::Float64,
    t_target::Float64,
    ::LinLogMixedMethod
)
    if c1 > 0.0 && c2 > 0.0 && c2 < c1
        # Descending and both positive: use log-linear interpolation
        k = log(c1 / c2) / (t2 - t1)
        return c1 * exp(-k * (t_target - t1))
    else
        # Ascending or zero/negative: use linear interpolation
        slope = (c2 - c1) / (t2 - t1)
        return c1 + slope * (t_target - t1)
    end
end
