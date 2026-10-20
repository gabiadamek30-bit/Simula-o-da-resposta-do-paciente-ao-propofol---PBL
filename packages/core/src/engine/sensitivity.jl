export SensitivityResult, run_sensitivity

struct SensitivityResult
    plan::PerturbationPlan
    observation::Symbol
    base_metric_series::Vector{Float64}
    pert_metric_series::Vector{Float64}
    metrics::SensitivityMetric
    metadata::Dict{String,Any}
end

function run_sensitivity(
    spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec;
    plan::PerturbationPlan,
    observation::Symbol=:conc,
)
    base_res = simulate(spec, grid, solver)

    if !haskey(base_res.observations, observation)
        error("Base result missing observation $(observation)")
    end

    pert_params = apply_plan(spec.params, plan)
    pert_spec = ModelSpec(
        spec.kind, spec.name * "_sens_" * plan.name, pert_params, spec.doses
    )

    pert_res = simulate(pert_spec, grid, solver)

    base_series = base_res.observations[observation]
    pert_series = pert_res.observations[observation]

    m = compute_metrics(base_series, pert_series)

    md = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "plan_name" => plan.name,
        "observation" => String(observation),
        "model_kind" => string(typeof(spec.kind)),
        "base_name" => spec.name,
        "pert_name" => pert_spec.name,
    )

    return SensitivityResult(plan, observation, base_series, pert_series, m, md)
end
