export validate, evaluate

function validate(spec::PDSpec{DirectEmax,DirectEmaxParams})
    E0 = spec.params.E0
    Emax = spec.params.Emax
    EC50 = spec.params.EC50

    # E0 can be any real number depending on baseline definition
    _require_positive("Emax", Emax)
    _require_positive("EC50", EC50)

    return nothing
end

function evaluate(spec::PDSpec{DirectEmax,DirectEmaxParams}, input_series::Vector{Float64})
    validate(spec)

    E0 = spec.params.E0
    Emax = spec.params.Emax
    EC50 = spec.params.EC50

    out = Vector{Float64}(undef, length(input_series))
    for i in eachindex(input_series)
        C = input_series[i]
        out[i] = E0 + (Emax * C) / (EC50 + C)
    end
    return out
end
