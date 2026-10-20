export emin, time_below, auc_above_baseline

"""
Minimum value of a response curve.
"""
function emin(t::Vector{Float64}, y::Vector{Float64})
    if length(t) != length(y)
        error("Length mismatch in emin")
    end
    m = Inf
    for v in y
        if v < m
            m = v
        end
    end
    return m
end

"""
Total time where response y is below a threshold.
Assumes t is sorted, uses left-constant rule per interval [t[i-1], t[i]] with y[i-1].
"""
function time_below(t::Vector{Float64}, y::Vector{Float64}, thr::Float64)
    if length(t) != length(y)
        error("Length mismatch in time_below")
    end
    if length(t) < 2
        return 0.0
    end
    s = 0.0
    for i in 2:length(t)
        dt = t[i] - t[i - 1]
        if dt < 0
            error("t must be sorted for time_below")
        end
        if y[i - 1] < thr
            s += dt
        end
    end
    return s
end

"""
AUC of (baseline - y) where y is below baseline. This measures suppression burden.
Uses trapezoid on the transformed curve max(0, baseline - y).
"""
function auc_above_baseline(t::Vector{Float64}, y::Vector{Float64}, baseline::Float64)
    if length(t) != length(y)
        error("Length mismatch in auc_above_baseline")
    end
    if length(t) < 2
        return 0.0
    end
    s = 0.0
    for i in 2:length(t)
        dt = t[i] - t[i - 1]
        if dt < 0
            error("t must be sorted for auc_above_baseline")
        end
        a0 = max(0.0, baseline - y[i - 1])
        a1 = max(0.0, baseline - y[i])
        s += 0.5 * dt * (a0 + a1)
    end
    return s
end
