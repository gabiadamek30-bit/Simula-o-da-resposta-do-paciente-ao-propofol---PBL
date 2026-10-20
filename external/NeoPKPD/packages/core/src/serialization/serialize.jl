using JSON

export serialize_execution, write_execution_json

function _serialize_doses(doses::Vector{DoseEvent})
    return [Dict("time" => d.time, "amount" => d.amount, "duration" => d.duration) for d in doses]
end

function _serialize_model_spec(spec::ModelSpec)
    return Dict(
        "kind" => string(typeof(spec.kind)),
        "name" => spec.name,
        "params" => Dict(
            string(k) => getfield(spec.params, k) for k in fieldnames(typeof(spec.params))
        ),
        "doses" => _serialize_doses(spec.doses),
    )
end

function _serialize_pd_spec(spec::PDSpec)
    return Dict(
        "kind" => string(typeof(spec.kind)),
        "name" => spec.name,
        "params" => Dict(
            string(k) => getfield(spec.params, k) for k in fieldnames(typeof(spec.params))
        ),
        "input_observation" => String(spec.input_observation),
        "output_observation" => String(spec.output_observation),
    )
end

function _serialize_solver(solver::SolverSpec)
    return Dict(
        "alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "maxiters" => solver.maxiters,
    )
end

function _serialize_grid(grid::SimGrid)
    return Dict("t0" => grid.t0, "t1" => grid.t1, "saveat" => grid.saveat)
end

function _serialize_results(res::SimResult)
    return Dict(
        "t" => res.t,
        "states" => Dict(String(k) => v for (k, v) in res.states),
        "observations" => Dict(String(k) => v for (k, v) in res.observations),
        "metadata" => res.metadata,
    )
end

"""
Create a full execution artifact as a Dict.
"""
function serialize_execution(;
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec,
    result::SimResult,
    pd_spec::Union{Nothing,PDSpec}=nothing,
)
    artifact = Dict(
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "model_spec" => _serialize_model_spec(model_spec),
        "grid" => _serialize_grid(grid),
        "solver" => _serialize_solver(solver),
        "result" => _serialize_results(result),
    )

    if pd_spec !== nothing
        artifact["pd_spec"] = _serialize_pd_spec(pd_spec)
    end

    mode = "pk"
    if pd_spec !== nothing
        # Direct PD models (pk_then_pd mode)
        if pd_spec.kind isa DirectEmax || pd_spec.kind isa SigmoidEmax || pd_spec.kind isa BiophaseEquilibration
            mode = "pk_then_pd"
        # Coupled ODE PD models (pkpd_coupled mode)
        elseif pd_spec.kind isa IndirectResponseTurnover
            mode = "pkpd_coupled"
        end
    end
    artifact["execution_mode"] = mode
    artifact["semantics_fingerprint"] = semantics_fingerprint()

    return artifact
end

"""
Write execution artifact to a JSON file.
"""
function write_execution_json(path::AbstractString; kwargs...)
    artifact = serialize_execution(; kwargs...)
    open(path, "w") do io
        write(io, JSON.json(artifact; pretty=true))
    end
    return path
end
