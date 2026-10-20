# GSA Serialization - Sobol' and Morris results
# This file must be loaded after sensitivity_sobol.jl and sensitivity_morris.jl

export serialize_sobol_result, write_sobol_json
export serialize_morris_result, write_morris_json

function _serialize_parameter_bounds(bounds::ParameterBounds)
    return Dict(
        "params" => [String(p) for p in bounds.params],
        "lower" => bounds.lower,
        "upper" => bounds.upper,
    )
end

function _serialize_sobol_method(method::SobolMethod)
    return Dict(
        "method_type" => "Sobol",
        "base_sample_size" => method.base_sample_size,
        "compute_second_order" => method.compute_second_order,
        "bootstrap_samples" => method.bootstrap_samples,
        "bootstrap_ci_level" => method.bootstrap_ci_level,
    )
end

function _serialize_morris_method(method::MorrisMethod)
    return Dict(
        "method_type" => "Morris",
        "n_trajectories" => method.n_trajectories,
        "n_levels" => method.n_levels,
        "delta" => method.delta,
    )
end

function _serialize_gsa_spec(spec::GlobalSensitivitySpec{SobolMethod})
    return Dict(
        "method" => _serialize_sobol_method(spec.method),
        "bounds" => _serialize_parameter_bounds(spec.bounds),
        "observation" => String(spec.observation),
        "seed" => spec.seed,
    )
end

function _serialize_gsa_spec(spec::GlobalSensitivitySpec{MorrisMethod})
    return Dict(
        "method" => _serialize_morris_method(spec.method),
        "bounds" => _serialize_parameter_bounds(spec.bounds),
        "observation" => String(spec.observation),
        "seed" => spec.seed,
    )
end

function _serialize_sobol_index(idx::SobolIndex)
    return Dict(
        "Si" => idx.Si,
        "Si_ci_lower" => idx.Si_ci_lower,
        "Si_ci_upper" => idx.Si_ci_upper,
        "STi" => idx.STi,
        "STi_ci_lower" => idx.STi_ci_lower,
        "STi_ci_upper" => idx.STi_ci_upper,
    )
end

function _serialize_morris_index(idx::MorrisIndex)
    return Dict(
        "mu" => idx.mu,
        "mu_star" => idx.mu_star,
        "sigma" => idx.sigma,
    )
end

"""
    serialize_sobol_result(result::SobolResult)

Serialize a Sobol' sensitivity analysis result to a dictionary.
"""
function serialize_sobol_result(result::SobolResult)
    indices_dict = Dict{String,Any}()
    for (param, idx) in result.indices
        indices_dict[String(param)] = _serialize_sobol_index(idx)
    end

    second_order_dict = nothing
    if result.second_order !== nothing
        second_order_dict = Dict{String,Float64}()
        for ((p1, p2), val) in result.second_order
            key = "$(String(p1)),$(String(p2))"
            second_order_dict[key] = val
        end
    end

    return Dict(
        "artifact_type" => "sobol_sensitivity",
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "semantics_fingerprint" => semantics_fingerprint(),
        "gsa_spec" => _serialize_gsa_spec(result.spec),
        "params" => [String(p) for p in result.params],
        "indices" => indices_dict,
        "second_order" => second_order_dict,
        "n_evaluations" => result.n_evaluations,
        "convergence_metric" => result.convergence_metric,
        "output_variance" => result.output_variance,
        "computation_time" => result.computation_time,
        "metadata" => result.metadata,
    )
end

"""
    write_sobol_json(path::AbstractString, result::SobolResult)

Write a Sobol' result to a JSON file.
"""
function write_sobol_json(path::AbstractString, result::SobolResult)
    artifact = serialize_sobol_result(result)
    open(path, "w") do io
        JSON.print(io, artifact; indent=2)
    end
    return path
end

"""
    serialize_morris_result(result::MorrisResult)

Serialize a Morris sensitivity analysis result to a dictionary.
"""
function serialize_morris_result(result::MorrisResult)
    indices_dict = Dict{String,Any}()
    for (param, idx) in result.indices
        indices_dict[String(param)] = _serialize_morris_index(idx)
    end

    ee_dict = Dict{String,Vector{Float64}}()
    for (param, ees) in result.elementary_effects
        ee_dict[String(param)] = ees
    end

    return Dict(
        "artifact_type" => "morris_sensitivity",
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "semantics_fingerprint" => semantics_fingerprint(),
        "gsa_spec" => _serialize_gsa_spec(result.spec),
        "params" => [String(p) for p in result.params],
        "indices" => indices_dict,
        "elementary_effects" => ee_dict,
        "n_evaluations" => result.n_evaluations,
        "computation_time" => result.computation_time,
        "metadata" => result.metadata,
    )
end

"""
    write_morris_json(path::AbstractString, result::MorrisResult)

Write a Morris result to a JSON file.
"""
function write_morris_json(path::AbstractString, result::MorrisResult)
    artifact = serialize_morris_result(result)
    open(path, "w") do io
        JSON.print(io, artifact; indent=2)
    end
    return path
end
