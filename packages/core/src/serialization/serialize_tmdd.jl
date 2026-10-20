# =============================================================================
# TMDD Serialization Functions
# =============================================================================
#
# Serialization support for TMDD model artifacts, enabling golden validation.
# =============================================================================

export serialize_tmdd_execution, write_tmdd_execution_json

"""
Serialize TMDD spec to dictionary format.
"""
function _serialize_tmdd_spec(spec::TMDDSpec{K,P}) where {K,P}
    kind_info = Dict{String,Any}(
        "type" => string(K),
    )

    # Add approximation and route if available
    if hasfield(K, :approximation)
        kind_info["approximation"] = string(spec.kind.approximation)
    end
    if hasfield(K, :route)
        kind_info["route"] = string(spec.kind.route)
    end

    return Dict{String,Any}(
        "kind" => kind_info,
        "name" => spec.name,
        "params" => Dict(
            string(k) => getfield(spec.params, k) for k in fieldnames(P)
        ),
        "doses" => _serialize_doses(spec.doses),
        "target_units" => string(spec.target_units),
        "drug_units" => string(spec.drug_units),
    )
end

"""
Serialize TMDD simulation result to dictionary format.
"""
function _serialize_tmdd_result(result::TMDDSimResult)
    return Dict{String,Any}(
        "t" => result.t,
        "states" => Dict(string(k) => v for (k, v) in result.states),
        "observations" => Dict(string(k) => v for (k, v) in result.observations),
        "metadata" => result.metadata,
    )
end

"""
    serialize_tmdd_execution(; tmdd_spec, grid, solver, result)

Create a full TMDD execution artifact as a Dict for JSON serialization.

# Arguments
- `tmdd_spec::TMDDSpec`: TMDD model specification
- `grid::SimGrid`: Simulation time grid
- `solver::SolverSpec`: ODE solver configuration
- `result::TMDDSimResult`: Simulation results

# Returns
Dict suitable for JSON serialization
"""
function serialize_tmdd_execution(;
    tmdd_spec::TMDDSpec,
    grid::SimGrid,
    solver::SolverSpec,
    result::TMDDSimResult,
)
    artifact = Dict{String,Any}(
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "execution_mode" => "tmdd",
        "semantics_fingerprint" => semantics_fingerprint(),
        "tmdd_spec" => _serialize_tmdd_spec(tmdd_spec),
        "grid" => _serialize_grid(grid),
        "solver" => _serialize_solver(solver),
        "result" => _serialize_tmdd_result(result),
    )

    return artifact
end

"""
    write_tmdd_execution_json(path; kwargs...)

Write TMDD execution artifact to a JSON file.
"""
function write_tmdd_execution_json(path::AbstractString; kwargs...)
    artifact = serialize_tmdd_execution(; kwargs...)
    open(path, "w") do io
        write(io, JSON.json(artifact; pretty=true))
    end
    return path
end
