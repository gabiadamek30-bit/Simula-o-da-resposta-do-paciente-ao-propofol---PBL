export deserialize_sensitivity_execution, replay_sensitivity_execution
export deserialize_population_sensitivity_execution, replay_population_sensitivity_execution

function _parse_plan(d::Dict)::PerturbationPlan
    name = String(d["name"])
    perts_in = d["perturbations"]

    perts = Perturbation[]
    for item in perts_in
        kind_s = String(item["kind"])
        param = Symbol(String(item["param"]))
        delta = Float64(item["delta"])

        if kind_s == "RelativePerturbation"
            push!(perts, Perturbation(RelativePerturbation(), param, delta))
        elseif kind_s == "AbsolutePerturbation"
            push!(perts, Perturbation(AbsolutePerturbation(), param, delta))
        else
            error("Unsupported perturbation kind: $(kind_s)")
        end
    end

    return PerturbationPlan(name, perts)
end

function deserialize_sensitivity_execution(artifact::Dict)
    if haskey(artifact, "artifact_type")
        if String(artifact["artifact_type"]) != "sensitivity_single"
            error("Expected artifact_type=sensitivity_single")
        end
    end

    schema = String(artifact["artifact_schema_version"])
    if schema != ARTIFACT_SCHEMA_VERSION
        error(
            "Unsupported artifact schema version: $(schema). Expected: $(ARTIFACT_SCHEMA_VERSION)",
        )
    end

    spec = _parse_model_spec(_to_dict(artifact["model_spec"]))
    grid = _parse_grid(_to_dict(artifact["grid"]))
    solver = _parse_solver(_to_dict(artifact["solver"]))
    plan = _parse_plan(_to_dict(artifact["plan"]))

    obs = Symbol(String(artifact["observation"]))

    return (model_spec=spec, grid=grid, solver=solver, plan=plan, observation=obs)
end

function replay_sensitivity_execution(artifact::Dict)::SensitivityResult
    parsed = deserialize_sensitivity_execution(artifact)
    return run_sensitivity(
        parsed.model_spec,
        parsed.grid,
        parsed.solver;
        plan=parsed.plan,
        observation=parsed.observation,
    )
end

function deserialize_population_sensitivity_execution(artifact::Dict)
    if haskey(artifact, "artifact_type")
        if String(artifact["artifact_type"]) != "sensitivity_population"
            error("Expected artifact_type=sensitivity_population")
        end
    end

    schema = String(artifact["artifact_schema_version"])
    if schema != ARTIFACT_SCHEMA_VERSION
        error(
            "Unsupported artifact schema version: $(schema). Expected: $(ARTIFACT_SCHEMA_VERSION)",
        )
    end

    pop = _parse_population_spec(_to_dict(artifact["population_spec"]))
    grid = _parse_grid(_to_dict(artifact["grid"]))
    solver = _parse_solver(_to_dict(artifact["solver"]))
    plan = _parse_plan(_to_dict(artifact["plan"]))

    obs = Symbol(String(artifact["observation"]))

    probs = [Float64(x) for x in artifact["probs"]]

    return (
        population_spec=pop,
        grid=grid,
        solver=solver,
        plan=plan,
        observation=obs,
        probs=probs,
    )
end

function replay_population_sensitivity_execution(
    artifact::Dict
)::PopulationSensitivityResult
    parsed = deserialize_population_sensitivity_execution(artifact)
    return run_population_sensitivity(
        parsed.population_spec,
        parsed.grid,
        parsed.solver;
        plan=parsed.plan,
        observation=parsed.observation,
        probs=parsed.probs,
    )
end
