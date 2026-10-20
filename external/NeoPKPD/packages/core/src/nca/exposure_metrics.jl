# Exposure Metrics
# Primary PK exposure parameters for NCA

export nca_cmax, nca_tmax, nca_cmin, nca_clast, nca_tlast, nca_cavg
export find_cmax, find_tmax, find_clast, find_tlast
export nca_ctrough, nca_c_at_time, time_above_concentration, time_to_cmax

# =============================================================================
# Primary Exposure Metrics
# =============================================================================

"""
    nca_cmax(c)

Find maximum observed concentration (Cmax).

# Arguments
- `c::Vector{Float64}`: Concentration values

# Returns
- `Float64`: Maximum concentration
"""
function nca_cmax(c::Vector{Float64})
    @assert !isempty(c) "Concentration vector cannot be empty"
    return maximum(c)
end

"""
    nca_tmax(t, c)

Find time of maximum observed concentration (Tmax).

If multiple timepoints have the same Cmax, returns the first occurrence.

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values

# Returns
- `Float64`: Time of maximum concentration
"""
function nca_tmax(t::Vector{Float64}, c::Vector{Float64})
    @assert length(t) == length(c) "Time and concentration vectors must have same length"
    @assert !isempty(c) "Vectors cannot be empty"

    idx = argmax(c)
    return t[idx]
end

"""
    find_cmax(t, c)

Find Cmax and its index.

# Returns
- `Tuple{Float64, Int}`: (cmax, index)
"""
function find_cmax(t::Vector{Float64}, c::Vector{Float64})
    idx = argmax(c)
    return (c[idx], idx)
end

"""
    find_tmax(t, c)

Find Tmax and its index.

# Returns
- `Tuple{Float64, Int}`: (tmax, index)
"""
function find_tmax(t::Vector{Float64}, c::Vector{Float64})
    idx = argmax(c)
    return (t[idx], idx)
end

"""
    nca_cmin(c)

Find minimum observed concentration (Cmin).

For single dose studies, this is typically at the end.
For multiple dose studies, this is the trough concentration.

# Arguments
- `c::Vector{Float64}`: Concentration values

# Returns
- `Float64`: Minimum concentration
"""
function nca_cmin(c::Vector{Float64})
    @assert !isempty(c) "Concentration vector cannot be empty"
    return minimum(c)
end

"""
    nca_clast(t, c; lloq=0.0)

Find last measurable concentration (Clast).

Clast is the last concentration above the LLOQ (or > 0 if LLOQ not specified).

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `lloq::Float64`: Lower limit of quantification (default: 0.0)

# Returns
- `Float64`: Last measurable concentration
"""
function nca_clast(t::Vector{Float64}, c::Vector{Float64}; lloq::Float64 = 0.0)
    @assert length(t) == length(c) "Time and concentration vectors must have same length"

    # Find last index where concentration > lloq
    for i in length(c):-1:1
        if c[i] > lloq
            return c[i]
        end
    end

    # If no concentration above LLOQ, return 0
    return 0.0
end

"""
    nca_tlast(t, c; lloq=0.0)

Find time of last measurable concentration (Tlast).

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `lloq::Float64`: Lower limit of quantification (default: 0.0)

# Returns
- `Float64`: Time of last measurable concentration
"""
function nca_tlast(t::Vector{Float64}, c::Vector{Float64}; lloq::Float64 = 0.0)
    @assert length(t) == length(c) "Time and concentration vectors must have same length"

    # Find last index where concentration > lloq
    for i in length(c):-1:1
        if c[i] > lloq
            return t[i]
        end
    end

    # If no concentration above LLOQ, return last time
    return t[end]
end

"""
    find_clast(t, c; lloq=0.0)

Find Clast and its index.

# Returns
- `Tuple{Float64, Float64, Int}`: (clast, tlast, index)
"""
function find_clast(t::Vector{Float64}, c::Vector{Float64}; lloq::Float64 = 0.0)
    for i in length(c):-1:1
        if c[i] > lloq
            return (c[i], t[i], i)
        end
    end
    return (0.0, t[end], length(t))
end

"""
    find_tlast(t, c; lloq=0.0)

Find Tlast and its index.

# Returns
- `Tuple{Float64, Int}`: (tlast, index)
"""
function find_tlast(t::Vector{Float64}, c::Vector{Float64}; lloq::Float64 = 0.0)
    for i in length(c):-1:1
        if c[i] > lloq
            return (t[i], i)
        end
    end
    return (t[end], length(t))
end

"""
    nca_cavg(t, c, tau, config)

Calculate average concentration over a dosing interval (Cavg).

Cavg = AUC0-τ / τ

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `tau::Float64`: Dosing interval
- `config::NCAConfig`: NCA configuration

# Returns
- `Float64`: Average concentration
"""
function nca_cavg(t::Vector{Float64}, c::Vector{Float64}, tau::Float64, config::NCAConfig)
    @assert tau > 0.0 "Dosing interval must be positive"

    auc_tau = auc_0_tau(t, c, tau, config)
    return auc_tau / tau
end

# =============================================================================
# Additional Exposure Metrics
# =============================================================================

"""
    nca_ctrough(t, c, tau)

Find trough concentration at end of dosing interval.

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `tau::Float64`: Dosing interval

# Returns
- `Float64`: Trough concentration
"""
function nca_ctrough(t::Vector{Float64}, c::Vector{Float64}, tau::Float64)
    # Find concentration at time closest to tau
    idx = argmin(abs.(t .- tau))
    return c[idx]
end

"""
    nca_c_at_time(t, c, target_time)

Get concentration at a specific time point (with interpolation if needed).

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `target_time::Float64`: Target time

# Returns
- `Float64`: Concentration at target time
"""
function nca_c_at_time(t::Vector{Float64}, c::Vector{Float64}, target_time::Float64)
    n = length(t)

    # Check bounds
    if target_time <= t[1]
        return c[1]
    elseif target_time >= t[end]
        return c[end]
    end

    # Find bracketing indices
    idx_upper = findfirst(ti -> ti >= target_time, t)

    if idx_upper === nothing
        return c[end]
    end

    if t[idx_upper] == target_time
        return c[idx_upper]
    end

    idx_lower = idx_upper - 1

    # Interpolate
    t1, t2 = t[idx_lower], t[idx_upper]
    c1, c2 = c[idx_lower], c[idx_upper]

    # Use log-linear interpolation if descending with positive concentrations
    if c1 > 0.0 && c2 > 0.0 && c2 < c1
        k = log(c1 / c2) / (t2 - t1)
        return c1 * exp(-k * (target_time - t1))
    else
        # Linear interpolation
        slope = (c2 - c1) / (t2 - t1)
        return c1 + slope * (target_time - t1)
    end
end

"""
    time_above_concentration(t, c, threshold)

Calculate time above a concentration threshold.

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `threshold::Float64`: Concentration threshold

# Returns
- `Float64`: Total time where concentration > threshold
"""
function time_above_concentration(t::Vector{Float64}, c::Vector{Float64}, threshold::Float64)
    n = length(t)
    time_above = 0.0

    for i in 2:n
        dt = t[i] - t[i-1]

        # Both above threshold
        if c[i-1] > threshold && c[i] > threshold
            time_above += dt
        # Crossing from above to below
        elseif c[i-1] > threshold && c[i] <= threshold
            # Interpolate crossing time
            t_cross = t[i-1] + (threshold - c[i-1]) * dt / (c[i] - c[i-1])
            time_above += t_cross - t[i-1]
        # Crossing from below to above
        elseif c[i-1] <= threshold && c[i] > threshold
            t_cross = t[i-1] + (threshold - c[i-1]) * dt / (c[i] - c[i-1])
            time_above += t[i] - t_cross
        end
    end

    return time_above
end

"""
    time_to_cmax(t, c)

Find time to reach maximum concentration (equivalent to Tmax).

# Returns
- `Float64`: Time to Cmax
"""
function time_to_cmax(t::Vector{Float64}, c::Vector{Float64})
    return nca_tmax(t, c)
end
