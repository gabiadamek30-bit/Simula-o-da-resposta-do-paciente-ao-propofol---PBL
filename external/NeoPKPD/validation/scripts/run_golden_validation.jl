using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD

function _load_json(path::String)
    return read_execution_json(path)
end

function _require_equal(a, b, label::String)
    if a != b
        error("Mismatch for $(label). Expected $(a), got $(b)")
    end
end

# Tolerance for cross-platform validation: 1e-10 allows for typical floating-point
# variations across different platforms, Julia versions, and BLAS implementations
# while still ensuring ~10 significant figures of agreement.
function _compare_vectors(a::Vector{Float64}, b::Vector{Float64}, label::String; rtol=1e-10, atol=1e-10)
    if length(a) != length(b)
        error("Length mismatch for $(label). Expected $(length(a)), got $(length(b))")
    end
    for i in eachindex(a)
        if !isapprox(a[i], b[i]; rtol=rtol, atol=atol)
            diff = abs(a[i] - b[i])
            rel_diff = a[i] != 0.0 ? diff / abs(a[i]) : diff
            error("Value mismatch for $(label) at index $(i). Expected $(a[i]), got $(b[i]) (diff=$(diff), rel=$(rel_diff))")
        end
    end
end

function _compare_dict_of_series(expected::Dict{String, Any}, actual::Dict{Symbol, Vector{Float64}}, label::String)
    exp_keys = sort(collect(keys(expected)))
    act_keys = sort([String(k) for k in keys(actual)])

    _require_equal(exp_keys, act_keys, "$(label) keys")

    for k in exp_keys
        exp_v = [Float64(x) for x in expected[k]]
        act_v = actual[Symbol(k)]
        _compare_vectors(exp_v, act_v, "$(label).$(k)")
    end
end

function validate_tmdd(path::String, artifact::Dict)
    # Check semantics fingerprint if present
    if haskey(artifact, "semantics_fingerprint")
        current = semantics_fingerprint()
        stored = Dict{String, Any}(artifact["semantics_fingerprint"])

        for (k, v) in current
            if !haskey(stored, k)
                error("Stored semantics fingerprint missing key $(k) in $(path)")
            end
            if String(stored[k]) != String(v)
                error(
                    "Semantics version mismatch for $(k) in $(path). " *
                    "Stored=$(stored[k]), Current=$(v). " *
                    "Golden update requires intentional semantics bump."
                )
            end
        end
    end

    # Check schema version
    schema = String(artifact["artifact_schema_version"])
    if schema != ARTIFACT_SCHEMA_VERSION
        error("Artifact schema version mismatch in $(path). Expected $(ARTIFACT_SCHEMA_VERSION), got $(schema)")
    end

    # Replay TMDD simulation
    replay = replay_tmdd_execution(artifact)

    # Get stored result
    stored = Dict{String, Any}(artifact["result"])
    stored_t = [Float64(x) for x in stored["t"]]
    stored_states = Dict{String, Any}(stored["states"])
    stored_obs = Dict{String, Any}(stored["observations"])

    # Compare time grid
    if replay.t != stored_t
        error("Time grid mismatch in $(path). Stored t does not match replay t")
    end

    # Compare states
    _compare_dict_of_series(stored_states, replay.states, "tmdd.states")

    # Compare observations
    _compare_dict_of_series(stored_obs, replay.observations, "tmdd.observations")

    return true
end

function validate_one(path::String)
    artifact = _load_json(path)

    atype = "single"
    if haskey(artifact, "artifact_type")
        atype = String(artifact["artifact_type"])
    end

    # Check execution_mode for TMDD
    exec_mode = get(artifact, "execution_mode", "pk")
    if exec_mode == "tmdd"
        return validate_tmdd(path, artifact)
    end

    if atype == "population"
        # Fingerprint enforcement if present
        if haskey(artifact, "semantics_fingerprint")
            current = semantics_fingerprint()
            storedfp = Dict{String, Any}(artifact["semantics_fingerprint"])
            for (k, v) in current
                if !haskey(storedfp, k)
                    error("Stored semantics fingerprint missing key $(k) in $(path)")
                end
                if String(storedfp[k]) != String(v)
                    error(
                        "Semantics version mismatch for $(k) in $(path). " *
                        "Stored=$(storedfp[k]), Current=$(v). Golden update requires intentional semantics bump."
                    )
                end
            end
        end

        replay = replay_population_execution(artifact)

        stored = artifact["population_result"]
        stored_meta = Dict{String, Any}(stored["metadata"])
        stored_params = stored["params"]
        stored_inds = stored["individuals"]

        if length(stored_inds) != length(replay.individuals)
            error("Population size mismatch in $(path)")
        end

        # Compare each individual result series
        for i in 1:length(replay.individuals)
            s_i = stored_inds[i]
            stored_t = [Float64(x) for x in s_i["t"]]
            if replay.individuals[i].t != stored_t
                error("Time grid mismatch for individual $(i) in $(path)")
            end

            stored_states = Dict{String, Any}(s_i["states"])
            stored_obs = Dict{String, Any}(s_i["observations"])

            _compare_dict_of_series(stored_states, replay.individuals[i].states, "pop.states[$i]")
            _compare_dict_of_series(stored_obs, replay.individuals[i].observations, "pop.observations[$i]")
        end

        if !haskey(stored, "summaries")
            error("Stored population artifact missing summaries in $(path)")
        end

        stored_summaries = Dict{String, Any}(stored["summaries"])

        # Replay summaries exist on PopulationResult
        replay_summaries = replay.summaries

        exp_keys = sort(collect(keys(stored_summaries)))
        act_keys = sort([String(k) for k in keys(replay_summaries)])
        _require_equal(exp_keys, act_keys, "population summaries keys")

        for k in exp_keys
            s_exp = stored_summaries[k]
            s_act = replay_summaries[Symbol(k)]

            # Compare mean and median series
            exp_mean = [Float64(x) for x in s_exp["mean"]]
            exp_median = [Float64(x) for x in s_exp["median"]]

            _compare_vectors(exp_mean, s_act.mean, "summary.$(k).mean")
            _compare_vectors(exp_median, s_act.median, "summary.$(k).median")

            # Compare quantiles
            exp_q = Dict{String, Any}(s_exp["quantiles"])
            act_q = s_act.quantiles

            exp_q_keys = sort(collect(keys(exp_q)))
            act_q_keys = sort([string(p) for p in keys(act_q)])
            _require_equal(exp_q_keys, act_q_keys, "summary.$(k).quantile keys")

            for pk in exp_q_keys
                p = parse(Float64, pk)
                exp_v = [Float64(x) for x in exp_q[pk]]
                act_v = act_q[p]
                _compare_vectors(exp_v, act_v, "summary.$(k).q$(pk)")
            end
        end

        # Enforce semantics version keys in metadata
        if !haskey(stored_meta, "event_semantics_version")
            error("Stored population artifact missing event_semantics_version in $(path)")
        end
        if !haskey(stored_meta, "solver_semantics_version")
            error("Stored population artifact missing solver_semantics_version in $(path)")
        end

        _require_equal(String(stored_meta["event_semantics_version"]), EVENT_SEMANTICS_VERSION, "event_semantics_version")
        _require_equal(String(stored_meta["solver_semantics_version"]), SOLVER_SEMANTICS_VERSION, "solver_semantics_version")

        return true
    end

        if atype == "sensitivity_single"
        if haskey(artifact, "semantics_fingerprint")
            current = semantics_fingerprint()
            storedfp = Dict{String, Any}(artifact["semantics_fingerprint"])
            for (k, v) in current
                if !haskey(storedfp, k)
                    error("Stored semantics fingerprint missing key $(k) in $(path)")
                end
                if String(storedfp[k]) != String(v)
                    error("Semantics version mismatch for $(k) in $(path). Stored=$(storedfp[k]), Current=$(v)")
                end
            end
        end

        replay = replay_sensitivity_execution(artifact)

        base_series = [Float64(x) for x in artifact["base_series"]]
        pert_series = [Float64(x) for x in artifact["pert_series"]]

        _compare_vectors(base_series, replay.base_metric_series, "sensitivity.base_series")
        _compare_vectors(pert_series, replay.pert_metric_series, "sensitivity.pert_series")

        return true
    end

    if atype == "sensitivity_population"
        if haskey(artifact, "semantics_fingerprint")
            current = semantics_fingerprint()
            storedfp = Dict{String, Any}(artifact["semantics_fingerprint"])
            for (k, v) in current
                if !haskey(storedfp, k)
                    error("Stored semantics fingerprint missing key $(k) in $(path)")
                end
                if String(storedfp[k]) != String(v)
                    error("Semantics version mismatch for $(k) in $(path). Stored=$(storedfp[k]), Current=$(v)")
                end
            end
        end

        replay = replay_population_sensitivity_execution(artifact)

        base_mean = [Float64(x) for x in artifact["base_mean"]]
        pert_mean = [Float64(x) for x in artifact["pert_mean"]]

        _compare_vectors(base_mean, replay.base_summary_mean, "pop_sens.base_mean")
        _compare_vectors(pert_mean, replay.pert_summary_mean, "pop_sens.pert_mean")

        return true
    end


    if haskey(artifact, "semantics_fingerprint")
        current = semantics_fingerprint()
        stored = Dict{String, Any}(artifact["semantics_fingerprint"])

        for (k, v) in current
            if !haskey(stored, k)
                error("Stored semantics fingerprint missing key $(k) in $(path)")
            end
            if String(stored[k]) != String(v)
                error(
                    "Semantics version mismatch for $(k) in $(path). " *
                    "Stored=$(stored[k]), Current=$(v). " *
                    "Golden update requires intentional semantics bump."
                )
            end
        end
    end


    schema = String(artifact["artifact_schema_version"])
    if schema != ARTIFACT_SCHEMA_VERSION
        error("Artifact schema version mismatch in $(path). Expected $(ARTIFACT_SCHEMA_VERSION), got $(schema)")
    end

    stored = artifact["result"]
    stored_t = [Float64(x) for x in stored["t"]]
    stored_states = Dict{String, Any}(stored["states"])
    stored_obs = Dict{String, Any}(stored["observations"])
    stored_meta = Dict{String, Any}(stored["metadata"])

    replay = replay_execution(artifact)

    if replay.t != stored_t
        error("Time grid mismatch in $(path). Stored t does not match replay t")
    end

    _compare_dict_of_series(stored_states, replay.states, "states")
    _compare_dict_of_series(stored_obs, replay.observations, "observations")

    if !haskey(stored_meta, "event_semantics_version")
        error("Stored artifact missing event_semantics_version in $(path)")
    end
    if !haskey(stored_meta, "solver_semantics_version")
        error("Stored artifact missing solver_semantics_version in $(path)")
    end

    _require_equal(String(stored_meta["event_semantics_version"]), EVENT_SEMANTICS_VERSION, "event_semantics_version")
    _require_equal(String(stored_meta["solver_semantics_version"]), SOLVER_SEMANTICS_VERSION, "solver_semantics_version")

    return true
end

function main()
    golden_dir = "validation/golden"
    files = filter(f -> endswith(f, ".json"), readdir(golden_dir; join=true))

    if isempty(files)
        error("No golden artifacts found in $(golden_dir)")
    end

    println("Running golden validation on $(length(files)) artifact(s)")

    for f in sort(files)
        println("Validating: $(f)")
        validate_one(f)
    end

    println("Golden validation passed")
end

main()
