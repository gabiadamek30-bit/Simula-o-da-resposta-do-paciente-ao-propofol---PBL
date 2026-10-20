export covariate_value_at, covariate_boundary_times

function _require_sorted_unique(xs::Vector{Float64}, name::String)
    if !issorted(xs)
        error("$(name) must be sorted")
    end
    for i in 2:length(xs)
        if xs[i] == xs[i - 1]
            error("$(name) must have unique times")
        end
    end
    return nothing
end

function covariate_value_at(s::TimeCovariateSeries{StepTimeCovariate}, t::Float64)
    ts = s.times
    vs = s.values
    _require_sorted_unique(ts, "TimeCovariateSeries.times")
    if length(ts) != length(vs)
        error("TimeCovariateSeries values length mismatch")
    end

    if t <= ts[1]
        return vs[1]
    end
    if t >= ts[end]
        return vs[end]
    end

    # rightmost ts[k] <= t
    k = searchsortedlast(ts, t)
    return vs[k]
end

function covariate_value_at(s::TimeCovariateSeries{LinearTimeCovariate}, t::Float64)
    ts = s.times
    vs = s.values
    _require_sorted_unique(ts, "TimeCovariateSeries.times")
    if length(ts) != length(vs)
        error("TimeCovariateSeries values length mismatch")
    end

    if t <= ts[1]
        return vs[1]
    end
    if t >= ts[end]
        return vs[end]
    end

    hi = searchsortedfirst(ts, t)
    lo = hi - 1

    t0 = ts[lo]
    t1 = ts[hi]
    v0 = vs[lo]
    v1 = vs[hi]

    w = (t - t0) / (t1 - t0)
    return (1.0 - w) * v0 + w * v1
end

"""
Return boundary times where covariate parameters may change.
For Step: every knot time is a boundary.
For Linear: every knot time is a boundary (piecewise linear).
"""
function covariate_boundary_times(tvc::TimeVaryingCovariates)
    out = Float64[]
    for (_, s_any) in tvc.series
        s = s_any
        append!(out, s.times)
    end
    sort!(out)

    # unique
    uniq = Float64[]
    for t in out
        if isempty(uniq) || t != uniq[end]
            push!(uniq, t)
        end
    end
    return uniq
end
