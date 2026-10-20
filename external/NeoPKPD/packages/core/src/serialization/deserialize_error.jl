# Deserialization for Residual Error Models

export deserialize_error_spec

const _ERROR_KIND_MAP = Dict{String,Function}(
    "AdditiveError" => () -> AdditiveError(),
    "ProportionalError" => () -> ProportionalError(),
    "CombinedError" => () -> CombinedError(),
    "ExponentialError" => () -> ExponentialError(),
)

"""
Deserialize a Dict to a ResidualErrorSpec.
"""
function deserialize_error_spec(d::Dict)::ResidualErrorSpec
    kind_str = String(d["kind"])

    if !haskey(_ERROR_KIND_MAP, kind_str)
        error("Unsupported error kind: $kind_str")
    end

    kind = _ERROR_KIND_MAP[kind_str]()
    params_d = d["params"]
    observation = Symbol(String(d["observation"]))
    seed = UInt64(d["seed"])

    if kind isa AdditiveError
        sigma = Float64(params_d["sigma"])
        params = AdditiveErrorParams(sigma)
        return ResidualErrorSpec(kind, params, observation, seed)
    elseif kind isa ProportionalError
        sigma = Float64(params_d["sigma"])
        params = ProportionalErrorParams(sigma)
        return ResidualErrorSpec(kind, params, observation, seed)
    elseif kind isa CombinedError
        sigma_add = Float64(params_d["sigma_add"])
        sigma_prop = Float64(params_d["sigma_prop"])
        params = CombinedErrorParams(sigma_add, sigma_prop)
        return ResidualErrorSpec(kind, params, observation, seed)
    elseif kind isa ExponentialError
        sigma = Float64(params_d["sigma"])
        params = ExponentialErrorParams(sigma)
        return ResidualErrorSpec(kind, params, observation, seed)
    end

    error("Internal error: error kind parsed but not handled")
end
