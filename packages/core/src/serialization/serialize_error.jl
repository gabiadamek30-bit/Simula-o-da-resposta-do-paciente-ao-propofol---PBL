# Serialization for Residual Error Models

export serialize_error_spec

"""
Serialize a ResidualErrorSpec to a Dict.
"""
function serialize_error_spec(spec::ResidualErrorSpec{AdditiveError})
    return Dict(
        "kind" => "AdditiveError",
        "params" => Dict("sigma" => spec.params.sigma),
        "observation" => String(spec.observation),
        "seed" => spec.seed
    )
end

function serialize_error_spec(spec::ResidualErrorSpec{ProportionalError})
    return Dict(
        "kind" => "ProportionalError",
        "params" => Dict("sigma" => spec.params.sigma),
        "observation" => String(spec.observation),
        "seed" => spec.seed
    )
end

function serialize_error_spec(spec::ResidualErrorSpec{CombinedError})
    return Dict(
        "kind" => "CombinedError",
        "params" => Dict(
            "sigma_add" => spec.params.sigma_add,
            "sigma_prop" => spec.params.sigma_prop
        ),
        "observation" => String(spec.observation),
        "seed" => spec.seed
    )
end

function serialize_error_spec(spec::ResidualErrorSpec{ExponentialError})
    return Dict(
        "kind" => "ExponentialError",
        "params" => Dict("sigma" => spec.params.sigma),
        "observation" => String(spec.observation),
        "seed" => spec.seed
    )
end
