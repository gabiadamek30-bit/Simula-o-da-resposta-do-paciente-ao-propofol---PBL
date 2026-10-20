using StableRNGs
using Distributions

export derive_occasions, sample_iov_kappas, apply_iov, occasion_index_at_time

function validate(iov::IOVSpec{LogNormalIIV})
    for (k, π) in iov.pis
        _require_positive("pi for $(k)", π)
    end
    if iov.occasion_def.mode != :dose_times
        error(
            "Unsupported occasion_def.mode: $(iov.occasion_def.mode). Supported: :dose_times",
        )
    end
    return nothing
end

"""
derive_occasions returns sorted unique occasion boundary times.

Semantics v1:
- t0 is the start of occasion 1
- each unique dose time strictly greater than t0 starts a new occasion
- boundaries returned include t0 and subsequent occasion starts within (t0, t1]
"""
function derive_occasions(doses::Vector{DoseEvent}, grid::SimGrid, def::OccasionDefinition)
    if def.mode != :dose_times
        error("Unsupported occasion definition mode: $(def.mode)")
    end

    t0 = grid.t0
    t1 = grid.t1

    ts = Float64[t0]
    for d in doses
        if d.time > t0 && d.time <= t1
            push!(ts, d.time)
        end
    end

    # unique and sorted
    sort!(ts)
    unique_ts = Float64[]
    for t in ts
        if isempty(unique_ts) || t != unique_ts[end]
            push!(unique_ts, t)
        end
    end

    return unique_ts
end

"""
sample_iov_kappas returns a vector of Dicts, one per occasion, with kappa values.

kappa ~ Normal(0, pi^2) per parameter.
"""
function sample_iov_kappas(pis::Dict{Symbol,Float64}, n_occ::Int, seed::UInt64)
    rng = StableRNG(seed)
    kappas = Vector{Dict{Symbol,Float64}}(undef, n_occ)

    for occ in 1:n_occ
        d = Dict{Symbol,Float64}()
        for (param, π) in pis
            d[param] = rand(rng, Normal(0.0, π))
        end
        kappas[occ] = d
    end

    return kappas
end

"""
Apply IOV kappas to typed params:
theta_occ = theta_base * exp(kappa) for specified parameters.
"""
function apply_iov(params, kappa::Dict{Symbol,Float64})
    T = typeof(params)
    fn = fieldnames(T)

    vals = Dict{Symbol,Float64}()
    for f in fn
        θ = Float64(getfield(params, f))
        if haskey(kappa, f)
            θ = θ * exp(kappa[f])
        end
        vals[f] = θ
    end

    return T((vals[f] for f in fn)...), vals
end

"""
Return the occasion index (1-based) for a given time t.
Returns the index of the rightmost occasion start <= t.
"""
function occasion_index_at_time(occ_starts::Vector{Float64}, t::Float64)
    # rightmost occ_start <= t
    return searchsortedlast(occ_starts, t)
end
