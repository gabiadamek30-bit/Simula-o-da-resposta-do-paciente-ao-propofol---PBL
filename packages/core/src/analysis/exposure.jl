export cmax, auc_trapezoid

"""
Maximum concentration.
"""
function cmax(t::Vector{Float64}, c::Vector{Float64})
    if length(t) != length(c)
        error("Length mismatch in cmax")
    end
    m = -Inf
    for x in c
        if x > m
            m = x
        end
    end
    return m
end

"""
Trapezoidal AUC on the provided grid.
Assumes t is sorted and matches c length.
"""
function auc_trapezoid(t::Vector{Float64}, c::Vector{Float64})
    if length(t) != length(c)
        error("Length mismatch in auc_trapezoid")
    end
    if length(t) < 2
        return 0.0
    end
    s = 0.0
    for i in 2:length(t)
        dt = t[i] - t[i - 1]
        if dt < 0
            error("t must be sorted for auc_trapezoid")
        end
        s += 0.5 * dt * (c[i] + c[i - 1])
    end
    return s
end
