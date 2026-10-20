export PopulationSensitivityResult, run_population_sensitivity

struct PopulationSensitivityResult
    plan::PerturbationPlan
    observation::Symbol
    probs::Vector{Float64}
    base_summary_mean::Vector{Float64}
    pert_summary_mean::Vector{Float64}
    base_quantiles::Dict{Float64,Vector{Float64}}
    pert_quantiles::Dict{Float64,Vector{Float64}}
    metrics_mean::SensitivityMetric
    metadata::Dict{String,Any}
end

function run_population_sensitivity(
    pop::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    plan::PerturbationPlan,
    observation::Symbol=:conc,
    probs::Vector{Float64}=[0.05, 0.95],
)
    base_res = simulate_population(pop, grid, solver)

    if !haskey(base_res.summaries, observation)
        error("Base population missing summary for $(observation)")
    end

    base_summary = base_res.summaries[observation]

    pert_params = apply_plan(pop.base_model_spec.params, plan)
    pert_base = ModelSpec(
        pop.base_model_spec.kind,
        pop.base_model_spec.name * "_sens_" * plan.name,
        pert_params,
        pop.base_model_spec.doses,
    )
    pert_pop = PopulationSpec(
        pert_base, pop.iiv, pop.iov, pop.covariate_model, pop.covariates
    )

    pert_res = simulate_population(pert_pop, grid, solver)

    if !haskey(pert_res.summaries, observation)
        error("Perturbed population missing summary for $(observation)")
    end

    pert_summary = pert_res.summaries[observation]

    m_mean = compute_metrics(base_summary.mean, pert_summary.mean)

    md = Dict{String,Any}(
        "engine_version" => "0.1.1",
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "plan_name" => plan.name,
        "observation" => String(observation),
        "population_n" => base_res.metadata["n"],
        "seed" => base_res.metadata["seed"],
    )

    return PopulationSensitivityResult(
        plan,
        observation,
        probs,
        base_summary.mean,
        pert_summary.mean,
        base_summary.quantiles,
        pert_summary.quantiles,
        m_mean,
        md,
    )
end
