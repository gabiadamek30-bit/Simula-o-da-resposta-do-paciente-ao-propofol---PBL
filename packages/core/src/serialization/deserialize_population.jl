export deserialize_population_execution, replay_population_execution

function _parse_iiv(d)::Union{Nothing,IIVSpec}
    if d === nothing
        return nothing
    end

    kind_str = String(d["kind"])

    # Normalize both qualified and unqualified forms
    if kind_str == "LogNormalIIV" || kind_str == "NeoPKPD.LogNormalIIV"
        kind = LogNormalIIV()
    else
        error("Unsupported IIV kind in artifact: $(kind_str)")
    end

    seed = UInt64(Int(d["seed"]))
    n_ind = Int(d["n"])

    # Check if full covariance matrix is present
    if haskey(d, "omega_matrix") && d["omega_matrix"] !== nothing
        om_data = d["omega_matrix"]
        param_names = [Symbol(String(p)) for p in om_data["param_names"]]
        n_params = Int(om_data["n_params"])

        # Reconstruct matrix from flattened vector
        matrix_vec = [Float64(v) for v in om_data["matrix"]]
        matrix = reshape(matrix_vec, n_params, n_params)

        omega_matrix = OmegaMatrix(param_names, matrix)
        return IIVSpec(kind, omega_matrix, seed, n_ind)
    else
        # Fall back to diagonal-only (backward compatible)
        omegas = Dict{Symbol,Float64}()
        for (k, v) in d["omegas"]
            omegas[Symbol(String(k))] = Float64(v)
        end
        return IIVSpec(kind, omegas, seed, n_ind)
    end
end

function _parse_time_varying(d)::Union{Nothing,TimeVaryingCovariates}
    if d === nothing
        return nothing
    end
    series = Dict{Symbol,Any}()
    for (k, v) in d
        kind_s = String(v["kind"])
        times = [Float64(x) for x in v["times"]]
        values = [Float64(x) for x in v["values"]]

        if kind_s == "StepTimeCovariate"
            series[Symbol(String(k))] = TimeCovariateSeries(
                StepTimeCovariate(), times, values
            )
        elseif kind_s == "LinearTimeCovariate"
            series[Symbol(String(k))] = TimeCovariateSeries(
                LinearTimeCovariate(), times, values
            )
        else
            error("Unsupported time covariate kind: $(kind_s)")
        end
    end
    return TimeVaryingCovariates(series)
end

function _parse_covariates(arr)::Vector{IndividualCovariates}
    covs = IndividualCovariates[]
    for item in arr
        # Handle both old format (flat dict) and new format (with values/time_varying)
        if haskey(item, "values")
            # New format
            values_in = item["values"]
            d = Dict{Symbol,Float64}()
            for (k, v) in values_in
                d[Symbol(String(k))] = Float64(v)
            end
            tv = if haskey(item, "time_varying")
                _parse_time_varying(item["time_varying"])
            else
                nothing
            end
            push!(covs, IndividualCovariates(d, tv))
        else
            # Old format (flat dict, no time_varying)
            d = Dict{Symbol,Float64}()
            for (k, v) in item
                d[Symbol(String(k))] = Float64(v)
            end
            push!(covs, IndividualCovariates(d, nothing))
        end
    end
    return covs
end

function _parse_population_spec(d)::PopulationSpec
    base = _parse_model_spec(_to_dict(d["base_model_spec"]))
    iiv = _parse_iiv(d["iiv"] === nothing ? nothing : _to_dict(d["iiv"]))
    iov = _parse_iov(d["iov"] === nothing ? nothing : _to_dict(d["iov"]))
    cm = _parse_covariate_model(
        d["covariate_model"] === nothing ? nothing : _to_dict(d["covariate_model"])
    )
    covs = _parse_covariates(d["covariates"])
    return PopulationSpec(base, iiv, iov, cm, covs)
end

"""
Deserialize a population execution artifact into objects.
Returns a NamedTuple:
- population_spec
- grid
- solver
- pd_spec (may be nothing)
"""
function deserialize_population_execution(artifact::Dict)
    schema = get(artifact, "artifact_schema_version", "1.0.0")
    schema_str = String(schema)

    # Handle backward compatibility via migration
    if schema_str != ARTIFACT_SCHEMA_VERSION
        if can_migrate(schema_str, ARTIFACT_SCHEMA_VERSION)
            migration_result = migrate_artifact(artifact)
            if migration_result.success
                artifact = migration_result.migrated_artifact
            else
                error(
                    "Failed to migrate artifact from $(schema_str) to $(ARTIFACT_SCHEMA_VERSION): $(join(migration_result.errors, "; "))"
                )
            end
        else
            error(
                "Unsupported artifact schema version: $(schema_str). Expected: $(ARTIFACT_SCHEMA_VERSION). No migration path available.",
            )
        end
    end

    if haskey(artifact, "artifact_type")
        if String(artifact["artifact_type"]) != "population"
            error("Expected artifact_type=population")
        end
    end

    pop_spec = _parse_population_spec(_to_dict(artifact["population_spec"]))
    grid = _parse_grid(_to_dict(artifact["grid"]))
    solver = _parse_solver(_to_dict(artifact["solver"]))

    pd_spec = nothing
    if haskey(artifact, "pd_spec")
        pd_spec = _parse_pd_spec(_to_dict(artifact["pd_spec"]))
    end

    return (population_spec=pop_spec, grid=grid, solver=solver, pd_spec=pd_spec)
end

"""
Replay a population artifact by rerunning simulate_population.
Returns PopulationResult.
"""
function replay_population_execution(artifact::Dict)::PopulationResult
    parsed = deserialize_population_execution(artifact)
    return simulate_population(
        parsed.population_spec, parsed.grid, parsed.solver; pd_spec=parsed.pd_spec
    )
end

function _parse_covariate_model(d)::Union{Nothing,CovariateModel}
    if d === nothing
        return nothing
    end
    name = String(d["name"])
    effs_in = d["effects"]
    effs = CovariateEffect[]
    for item in effs_in
        kind_s = String(item["kind"])
        param = Symbol(String(item["param"]))
        cov = Symbol(String(item["covariate"]))
        beta = Float64(item["beta"])
        ref = Float64(item["ref"])

        if kind_s == "LinearCovariate"
            push!(effs, CovariateEffect(LinearCovariate(), param, cov, beta, ref))
        elseif kind_s == "PowerCovariate"
            push!(effs, CovariateEffect(PowerCovariate(), param, cov, beta, ref))
        elseif kind_s == "ExpCovariate"
            push!(effs, CovariateEffect(ExpCovariate(), param, cov, beta, ref))
        else
            error("Unsupported covariate kind: $(kind_s)")
        end
    end
    return CovariateModel(name, effs)
end

function _parse_iov(d)::Union{Nothing,IOVSpec}
    if d === nothing
        return nothing
    end

    kind_str = String(d["kind"])
    # Normalize both qualified and unqualified forms
    if kind_str != "NeoPKPD.LogNormalIIV" && kind_str != "LogNormalIIV"
        error("Unsupported IOV kind in artifact: $(kind_str)")
    end

    pis = Dict{Symbol,Float64}()
    for (k, v) in d["pis"]
        pis[Symbol(String(k))] = Float64(v)
    end

    seed = UInt64(Int(d["seed"]))
    mode = Symbol(String(d["occasion_def"]["mode"]))

    return IOVSpec(LogNormalIIV(), pis, seed, OccasionDefinition(mode))
end
