# GSA Deserialization - Sobol' and Morris results
# This file must be loaded after sensitivity_sobol.jl and sensitivity_morris.jl

export deserialize_sobol_result, deserialize_morris_result

function _parse_parameter_bounds(d::Dict)::ParameterBounds
    params = [Symbol(p) for p in d["params"]]
    lower = [Float64(x) for x in d["lower"]]
    upper = [Float64(x) for x in d["upper"]]
    return ParameterBounds(params, lower, upper)
end

function _parse_sobol_method(d::Dict)::SobolMethod
    return SobolMethod(
        base_sample_size=Int(d["base_sample_size"]),
        compute_second_order=Bool(d["compute_second_order"]),
        bootstrap_samples=Int(d["bootstrap_samples"]),
        bootstrap_ci_level=Float64(d["bootstrap_ci_level"]),
    )
end

function _parse_morris_method(d::Dict)::MorrisMethod
    return MorrisMethod(
        n_trajectories=Int(d["n_trajectories"]),
        n_levels=Int(d["n_levels"]),
        delta=Float64(d["delta"]),
    )
end

function _parse_gsa_spec_sobol(d::Dict)::GlobalSensitivitySpec{SobolMethod}
    method = _parse_sobol_method(_to_dict(d["method"]))
    bounds = _parse_parameter_bounds(_to_dict(d["bounds"]))
    observation = Symbol(String(d["observation"]))
    seed = UInt64(d["seed"])
    return GlobalSensitivitySpec(method, bounds; observation=observation, seed=seed)
end

function _parse_gsa_spec_morris(d::Dict)::GlobalSensitivitySpec{MorrisMethod}
    method = _parse_morris_method(_to_dict(d["method"]))
    bounds = _parse_parameter_bounds(_to_dict(d["bounds"]))
    observation = Symbol(String(d["observation"]))
    seed = UInt64(d["seed"])
    return GlobalSensitivitySpec(method, bounds; observation=observation, seed=seed)
end

function _parse_sobol_index(d::Dict)::SobolIndex
    return SobolIndex(
        Float64(d["Si"]),
        Float64(d["Si_ci_lower"]),
        Float64(d["Si_ci_upper"]),
        Float64(d["STi"]),
        Float64(d["STi_ci_lower"]),
        Float64(d["STi_ci_upper"]),
    )
end

function _parse_morris_index(d::Dict)::MorrisIndex
    return MorrisIndex(
        Float64(d["mu"]),
        Float64(d["mu_star"]),
        Float64(d["sigma"]),
    )
end

"""
    deserialize_sobol_result(artifact::Dict)

Deserialize a Sobol' sensitivity analysis result from a dictionary.
"""
function deserialize_sobol_result(artifact::Dict)::SobolResult
    if haskey(artifact, "artifact_type")
        if String(artifact["artifact_type"]) != "sobol_sensitivity"
            error("Expected artifact_type=sobol_sensitivity")
        end
    end

    schema = String(artifact["artifact_schema_version"])
    if schema != ARTIFACT_SCHEMA_VERSION
        error("Unsupported artifact schema version: $(schema). Expected: $(ARTIFACT_SCHEMA_VERSION)")
    end

    gsa_spec = _parse_gsa_spec_sobol(_to_dict(artifact["gsa_spec"]))
    params = [Symbol(p) for p in artifact["params"]]

    indices = Dict{Symbol,SobolIndex}()
    indices_dict = _to_dict(artifact["indices"])
    for (param_str, idx_dict) in indices_dict
        param = Symbol(param_str)
        indices[param] = _parse_sobol_index(_to_dict(idx_dict))
    end

    second_order = nothing
    if haskey(artifact, "second_order") && artifact["second_order"] !== nothing
        second_order = Dict{Tuple{Symbol,Symbol},Float64}()
        so_dict = _to_dict(artifact["second_order"])
        for (key, val) in so_dict
            parts = split(key, ",")
            p1 = Symbol(parts[1])
            p2 = Symbol(parts[2])
            second_order[(p1, p2)] = Float64(val)
        end
    end

    metadata = Dict{String,Any}()
    if haskey(artifact, "metadata")
        for (k, v) in artifact["metadata"]
            metadata[String(k)] = v
        end
    end

    return SobolResult(
        gsa_spec,
        params,
        indices,
        second_order,
        Int(artifact["n_evaluations"]),
        Float64(artifact["convergence_metric"]),
        Float64(artifact["output_variance"]),
        Float64(artifact["computation_time"]),
        metadata,
    )
end

"""
    deserialize_morris_result(artifact::Dict)

Deserialize a Morris sensitivity analysis result from a dictionary.
"""
function deserialize_morris_result(artifact::Dict)::MorrisResult
    if haskey(artifact, "artifact_type")
        if String(artifact["artifact_type"]) != "morris_sensitivity"
            error("Expected artifact_type=morris_sensitivity")
        end
    end

    schema = String(artifact["artifact_schema_version"])
    if schema != ARTIFACT_SCHEMA_VERSION
        error("Unsupported artifact schema version: $(schema). Expected: $(ARTIFACT_SCHEMA_VERSION)")
    end

    gsa_spec = _parse_gsa_spec_morris(_to_dict(artifact["gsa_spec"]))
    params = [Symbol(p) for p in artifact["params"]]

    indices = Dict{Symbol,MorrisIndex}()
    indices_dict = _to_dict(artifact["indices"])
    for (param_str, idx_dict) in indices_dict
        param = Symbol(param_str)
        indices[param] = _parse_morris_index(_to_dict(idx_dict))
    end

    elementary_effects = Dict{Symbol,Vector{Float64}}()
    ee_dict = _to_dict(artifact["elementary_effects"])
    for (param_str, ees) in ee_dict
        param = Symbol(param_str)
        elementary_effects[param] = [Float64(x) for x in ees]
    end

    metadata = Dict{String,Any}()
    if haskey(artifact, "metadata")
        for (k, v) in artifact["metadata"]
            metadata[String(k)] = v
        end
    end

    return MorrisResult(
        gsa_spec,
        params,
        indices,
        elementary_effects,
        Int(artifact["n_evaluations"]),
        Float64(artifact["computation_time"]),
        metadata,
    )
end
