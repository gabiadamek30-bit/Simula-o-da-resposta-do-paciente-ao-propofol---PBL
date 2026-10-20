export serialize_sensitivity_execution, write_sensitivity_json
export serialize_population_sensitivity_execution, write_population_sensitivity_json

function _serialize_plan(plan::PerturbationPlan)
    perts = Vector{Any}(undef, length(plan.perturbations))
    for (i, p) in enumerate(plan.perturbations)
        kind = if p.kind isa RelativePerturbation
            "RelativePerturbation"
        else
            "AbsolutePerturbation"
        end
        perts[i] = Dict("kind" => kind, "param" => String(p.param), "delta" => p.delta)
    end
    return Dict("name" => plan.name, "perturbations" => perts)
end

function _serialize_metric(m::SensitivityMetric)
    return Dict(
        "max_abs_delta" => m.max_abs_delta,
        "max_rel_delta" => m.max_rel_delta,
        "l2_norm_delta" => m.l2_norm_delta,
    )
end

function serialize_sensitivity_execution(;
    model_spec::ModelSpec, grid::SimGrid, solver::SolverSpec, result::SensitivityResult
)
    return Dict(
        "artifact_type" => "sensitivity_single",
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "semantics_fingerprint" => semantics_fingerprint(),
        "model_spec" => _serialize_model_spec(model_spec),
        "grid" => _serialize_grid(grid),
        "solver" => _serialize_solver(solver),
        "plan" => _serialize_plan(result.plan),
        "observation" => String(result.observation),
        "base_series" => result.base_metric_series,
        "pert_series" => result.pert_metric_series,
        "metrics" => _serialize_metric(result.metrics),
        "metadata" => result.metadata,
    )
end

function write_sensitivity_json(path::AbstractString; kwargs...)
    artifact = serialize_sensitivity_execution(; kwargs...)
    open(path, "w") do io
        JSON.print(io, artifact; indent=2)
    end
    return path
end

function serialize_population_sensitivity_execution(;
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    result::PopulationSensitivityResult,
)
    return Dict(
        "artifact_type" => "sensitivity_population",
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "semantics_fingerprint" => semantics_fingerprint(),
        "population_spec" => _serialize_population_spec(population_spec),
        "grid" => _serialize_grid(grid),
        "solver" => _serialize_solver(solver),
        "plan" => _serialize_plan(result.plan),
        "observation" => String(result.observation),
        "probs" => result.probs,
        "base_mean" => result.base_summary_mean,
        "pert_mean" => result.pert_summary_mean,
        "base_quantiles" => Dict(string(p) => v for (p, v) in result.base_quantiles),
        "pert_quantiles" => Dict(string(p) => v for (p, v) in result.pert_quantiles),
        "metrics_mean" => _serialize_metric(result.metrics_mean),
        "metadata" => result.metadata,
    )
end

function write_population_sensitivity_json(path::AbstractString; kwargs...)
    artifact = serialize_population_sensitivity_execution(; kwargs...)
    open(path, "w") do io
        JSON.print(io, artifact; indent=2)
    end
    return path
end
