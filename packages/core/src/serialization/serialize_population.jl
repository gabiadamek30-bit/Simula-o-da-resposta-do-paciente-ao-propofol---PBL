export serialize_population_execution, write_population_json

function _serialize_iiv(iiv::Union{Nothing,IIVSpec})
    if iiv === nothing
        return nothing
    end

    result = Dict(
        "kind" => string(typeof(iiv.kind)),
        "omegas" => Dict(String(k) => v for (k, v) in iiv.omegas),
        "seed" => Int(iiv.seed),
        "n" => iiv.n,
    )

    # Serialize full covariance matrix if present
    if iiv.omega_matrix !== nothing
        om = iiv.omega_matrix
        result["omega_matrix"] = Dict(
            "param_names" => [String(p) for p in om.param_names],
            "matrix" => [om.matrix[i, j] for i in 1:size(om.matrix, 1), j in 1:size(om.matrix, 2)] |> vec,
            "n_params" => length(om.param_names),
        )
    end

    return result
end

function _serialize_time_varying(tvc::Union{Nothing,TimeVaryingCovariates})
    if tvc === nothing
        return nothing
    end
    out = Dict{String,Any}()
    for (name, s_any) in tvc.series
        s = s_any
        kind = s.kind isa StepTimeCovariate ? "StepTimeCovariate" : "LinearTimeCovariate"
        out[String(name)] = Dict("kind" => kind, "times" => s.times, "values" => s.values)
    end
    return out
end

function _serialize_individual_covariates(c::IndividualCovariates)
    return Dict(
        "values" => Dict(String(k) => v for (k, v) in c.values),
        "time_varying" => _serialize_time_varying(c.time_varying),
    )
end

function _serialize_population_spec(pop::PopulationSpec)
    return Dict(
        "base_model_spec" => _serialize_model_spec(pop.base_model_spec),
        "iiv" => _serialize_iiv(pop.iiv),
        "covariates" => [_serialize_individual_covariates(c) for c in pop.covariates],
        "covariate_model" => _serialize_covariate_model(pop.covariate_model),
        "iov" => _serialize_iov(pop.iov),
    )
end

function _serialize_population_result(popres::PopulationResult)
    return Dict(
        "metadata" => popres.metadata,
        "params" => [Dict(String(k) => v for (k, v) in d) for d in popres.params],
        "summaries" => Dict(
            String(k) => _serialize_population_summary(v) for (k, v) in popres.summaries
        ),
        "individuals" => [_serialize_results(r) for r in popres.individuals],
    )
end

"""
Create a population execution artifact as a Dict.

Stores:
- population_spec
- grid
- solver
- population_result (including per-individual outputs)
- pd_spec (optional, for coupled PKPD population simulations)
- schema version and semantics fingerprint
"""
function serialize_population_execution(;
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    result::PopulationResult,
    pd_spec::Union{Nothing,PDSpec}=nothing,
)
    artifact = Dict(
        "artifact_type" => "population",
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
        "semantics_fingerprint" => semantics_fingerprint(),
        "population_spec" => _serialize_population_spec(population_spec),
        "grid" => _serialize_grid(grid),
        "solver" => _serialize_solver(solver),
        "population_result" => _serialize_population_result(result),
    )

    if pd_spec !== nothing
        artifact["pd_spec"] = _serialize_pd_spec(pd_spec)
    end

    return artifact
end

function write_population_json(path::AbstractString; kwargs...)
    artifact = serialize_population_execution(; kwargs...)
    open(path, "w") do io
        JSON.print(io, artifact, 2)
    end
    return path
end

function _serialize_population_summary(s::PopulationSummary)
    return Dict(
        "observation" => String(s.observation),
        "probs" => s.probs,
        "mean" => s.mean,
        "median" => s.median,
        "quantiles" => Dict(string(p) => v for (p, v) in s.quantiles),
    )
end

function _serialize_covariate_model(cm::Union{Nothing,CovariateModel})
    if cm === nothing
        return nothing
    end
    effs = Vector{Any}(undef, length(cm.effects))
    for (i, e) in enumerate(cm.effects)
        kind = if e.kind isa LinearCovariate
            "LinearCovariate"
        elseif e.kind isa PowerCovariate
            "PowerCovariate"
        else
            "ExpCovariate"
        end

        effs[i] = Dict(
            "kind" => kind,
            "param" => String(e.param),
            "covariate" => String(e.covariate),
            "beta" => e.beta,
            "ref" => e.ref,
        )
    end
    return Dict("name" => cm.name, "effects" => effs)
end

function _serialize_iov(iov::Union{Nothing,IOVSpec})
    if iov === nothing
        return nothing
    end
    return Dict(
        "kind" => string(typeof(iov.kind)),
        "pis" => Dict(String(k) => v for (k, v) in iov.pis),
        "seed" => Int(iov.seed),
        "occasion_def" => Dict("mode" => String(iov.occasion_def.mode)),
    )
end
