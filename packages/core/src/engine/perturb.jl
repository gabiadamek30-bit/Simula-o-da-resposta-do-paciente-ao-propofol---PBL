export apply_perturbation, apply_plan

function _require_finite(name::String, x::Float64)
    if !isfinite(x)
        error("Expected finite value for $(name), got $(x)")
    end
    return nothing
end

function apply_perturbation(params, p::Perturbation{RelativePerturbation})
    T = typeof(params)
    fn = fieldnames(T)
    if !(p.param in fn)
        error("Parameter $(p.param) not found in $(T)")
    end

    vals = Vector{Float64}(undef, length(fn))
    for (i, f) in enumerate(fn)
        v = Float64(getfield(params, f))
        if f == p.param
            v = v * (1.0 + p.delta)
        end
        _require_finite(String(f), v)
        vals[i] = v
    end
    return T(vals...)
end

function apply_perturbation(params, p::Perturbation{AbsolutePerturbation})
    T = typeof(params)
    fn = fieldnames(T)
    if !(p.param in fn)
        error("Parameter $(p.param) not found in $(T)")
    end

    vals = Vector{Float64}(undef, length(fn))
    for (i, f) in enumerate(fn)
        v = Float64(getfield(params, f))
        if f == p.param
            v = v + p.delta
        end
        _require_finite(String(f), v)
        vals[i] = v
    end
    return T(vals...)
end

function apply_plan(params, plan::PerturbationPlan)
    out = params
    for p in plan.perturbations
        out = apply_perturbation(out, p)
    end
    return out
end
