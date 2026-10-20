using JSON

export read_execution_json, deserialize_execution, replay_execution

# -------------------------
# JSON normalization
# -------------------------

_to_dict(x::Dict) = x
_to_dict(x::JSON.Object) = Dict{String,Any}(x)

# -------------------------
# Kind resolution
# -------------------------

const _MODEL_KIND_MAP = Dict{String,Function}(
    "NeoPKPD.OneCompIVBolus" => () -> OneCompIVBolus(),
    "NeoPKPD.OneCompOralFirstOrder" => () -> OneCompOralFirstOrder(),
    "NeoPKPD.TwoCompIVBolus" => () -> TwoCompIVBolus(),
    "NeoPKPD.TwoCompOral" => () -> TwoCompOral(),
    "NeoPKPD.ThreeCompIVBolus" => () -> ThreeCompIVBolus(),
    "NeoPKPD.TransitAbsorption" => () -> TransitAbsorption(),
    "NeoPKPD.MichaelisMentenElimination" => () -> MichaelisMentenElimination(),
)

const _PD_KIND_MAP = Dict{String,Function}(
    "NeoPKPD.DirectEmax" => () -> DirectEmax(),
    "NeoPKPD.SigmoidEmax" => () -> SigmoidEmax(),
    "NeoPKPD.IndirectResponseTurnover" => () -> IndirectResponseTurnover(),
    "NeoPKPD.BiophaseEquilibration" => () -> BiophaseEquilibration(),
)

function _require_key(d::Dict, k::String)
    if !haskey(d, k)
        error("Missing required key: $(k)")
    end
    return d[k]
end

# -------------------------
# Low-level parsing
# -------------------------

function _parse_doses(arr)::Vector{DoseEvent}
    doses = DoseEvent[]
    for item in arr
        t = Float64(item["time"])
        a = Float64(item["amount"])
        # Support duration field (backward compatible - defaults to 0.0 for bolus)
        d = haskey(item, "duration") ? Float64(item["duration"]) : 0.0
        push!(doses, DoseEvent(t, a, d))
    end
    return doses
end

function _parse_solver(d::Dict)::SolverSpec
    alg = Symbol(String(d["alg"]))
    reltol = Float64(d["reltol"])
    abstol = Float64(d["abstol"])
    maxiters = Int(d["maxiters"])
    return SolverSpec(alg, reltol, abstol, maxiters)
end

function _parse_grid(d::Dict)::SimGrid
    t0 = Float64(d["t0"])
    t1 = Float64(d["t1"])
    saveat = [Float64(x) for x in d["saveat"]]
    return SimGrid(t0, t1, saveat)
end

function _parse_model_spec(d::Dict)::ModelSpec
    raw_kind = String(d["kind"])
    kind_str = occursin(".", raw_kind) ? raw_kind : "NeoPKPD.$raw_kind"

    if !haskey(_MODEL_KIND_MAP, kind_str)
        error("Unsupported model kind in artifact: $(kind_str)")
    end
    kind = _MODEL_KIND_MAP[kind_str]()

    name = String(d["name"])
    params_d = d["params"]
    doses = _parse_doses(d["doses"])

    if kind isa OneCompIVBolus
        CL = Float64(params_d["CL"])
        V = Float64(params_d["V"])
        params = OneCompIVBolusParams(CL, V)
        return ModelSpec(kind, name, params, doses)
    elseif kind isa OneCompOralFirstOrder
        Ka = Float64(params_d["Ka"])
        CL = Float64(params_d["CL"])
        V = Float64(params_d["V"])
        params = OneCompOralFirstOrderParams(Ka, CL, V)
        return ModelSpec(kind, name, params, doses)
    elseif kind isa TwoCompIVBolus
        CL = Float64(params_d["CL"])
        V1 = Float64(params_d["V1"])
        Q = Float64(params_d["Q"])
        V2 = Float64(params_d["V2"])
        params = TwoCompIVBolusParams(CL, V1, Q, V2)
        return ModelSpec(kind, name, params, doses)
    elseif kind isa TwoCompOral
        Ka = Float64(params_d["Ka"])
        CL = Float64(params_d["CL"])
        V1 = Float64(params_d["V1"])
        Q = Float64(params_d["Q"])
        V2 = Float64(params_d["V2"])
        params = TwoCompOralParams(Ka, CL, V1, Q, V2)
        return ModelSpec(kind, name, params, doses)
    elseif kind isa ThreeCompIVBolus
        CL = Float64(params_d["CL"])
        V1 = Float64(params_d["V1"])
        Q2 = Float64(params_d["Q2"])
        V2 = Float64(params_d["V2"])
        Q3 = Float64(params_d["Q3"])
        V3 = Float64(params_d["V3"])
        params = ThreeCompIVBolusParams(CL, V1, Q2, V2, Q3, V3)
        return ModelSpec(kind, name, params, doses)
    elseif kind isa TransitAbsorption
        N = Int(params_d["N"])
        Ktr = Float64(params_d["Ktr"])
        Ka = Float64(params_d["Ka"])
        CL = Float64(params_d["CL"])
        V = Float64(params_d["V"])
        params = TransitAbsorptionParams(N, Ktr, Ka, CL, V)
        return ModelSpec(kind, name, params, doses)
    elseif kind isa MichaelisMentenElimination
        Vmax = Float64(params_d["Vmax"])
        Km = Float64(params_d["Km"])
        V = Float64(params_d["V"])
        params = MichaelisMentenEliminationParams(Vmax, Km, V)
        return ModelSpec(kind, name, params, doses)
    end

    error("Internal error: model kind parsed but not handled")
end

function _parse_pd_spec(d::Dict)::PDSpec
    raw_kind = String(d["kind"])
    kind_str = occursin(".", raw_kind) ? raw_kind : "NeoPKPD.$raw_kind"

    if !haskey(_PD_KIND_MAP, kind_str)
        error("Unsupported PD kind in artifact: $(kind_str)")
    end
    kind = _PD_KIND_MAP[kind_str]()

    name = String(d["name"])
    params_d = d["params"]

    input_obs = Symbol(String(d["input_observation"]))
    output_obs = Symbol(String(d["output_observation"]))

    if kind isa DirectEmax
        E0 = Float64(params_d["E0"])
        Emax = Float64(params_d["Emax"])
        EC50 = Float64(params_d["EC50"])
        params = DirectEmaxParams(E0, Emax, EC50)
        return PDSpec(kind, name, params, input_obs, output_obs)
    elseif kind isa SigmoidEmax
        E0 = Float64(params_d["E0"])
        Emax = Float64(params_d["Emax"])
        EC50 = Float64(params_d["EC50"])
        gamma = Float64(params_d["gamma"])
        params = SigmoidEmaxParams(E0, Emax, EC50, gamma)
        return PDSpec(kind, name, params, input_obs, output_obs)
    elseif kind isa IndirectResponseTurnover
        Kin = Float64(params_d["Kin"])
        Kout = Float64(params_d["Kout"])
        R0 = Float64(params_d["R0"])
        Imax = Float64(params_d["Imax"])
        IC50 = Float64(params_d["IC50"])
        params = IndirectResponseTurnoverParams(Kin, Kout, R0, Imax, IC50)
        return PDSpec(kind, name, params, input_obs, output_obs)
    elseif kind isa BiophaseEquilibration
        ke0 = Float64(params_d["ke0"])
        E0 = Float64(params_d["E0"])
        Emax = Float64(params_d["Emax"])
        EC50 = Float64(params_d["EC50"])
        params = BiophaseEquilibrationParams(ke0, E0, Emax, EC50)
        return PDSpec(kind, name, params, input_obs, output_obs)
    end

    error("Internal error: PD kind parsed but not handled")
end

# -------------------------
# Public API
# -------------------------

"""
Read an execution artifact JSON file into a Dict.
"""
function read_execution_json(path::AbstractString)::Dict{String,Any}
    obj = JSON.parsefile(path)
    return Dict{String,Any}(obj)
end

"""
Deserialize an execution artifact into core objects.
Returns a NamedTuple containing:
- model_spec
- grid
- solver
- pd_spec (may be nothing)
- execution_mode (String)
"""
function deserialize_execution(artifact::Dict)
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

    model_spec = _parse_model_spec(_to_dict(_require_key(artifact, "model_spec")))
    grid = _parse_grid(_to_dict(_require_key(artifact, "grid")))
    solver = _parse_solver(_to_dict(_require_key(artifact, "solver")))

    pd_spec = nothing
    if haskey(artifact, "pd_spec")
        pd_spec = _parse_pd_spec(_to_dict(artifact["pd_spec"]))
    end

    # Optional field, inferred if absent
    mode = "pk"
    if haskey(artifact, "execution_mode")
        mode = String(artifact["execution_mode"])
    elseif pd_spec !== nothing
        # Direct PD models (pk_then_pd mode)
        if pd_spec.kind isa DirectEmax || pd_spec.kind isa SigmoidEmax || pd_spec.kind isa BiophaseEquilibration
            mode = "pk_then_pd"
        # Coupled ODE PD models (pkpd_coupled mode)
        elseif pd_spec.kind isa IndirectResponseTurnover
            mode = "pkpd_coupled"
        end
    end

    return (
        model_spec=model_spec,
        grid=grid,
        solver=solver,
        pd_spec=pd_spec,
        execution_mode=mode,
    )
end

"""
Replay an artifact by rerunning simulation from deserialized specs.
Supports:
- pk
- pk_then_pd (DirectEmax)
- pkpd_coupled (IndirectResponseTurnover)
"""
function replay_execution(artifact::Dict)::SimResult
    parsed = deserialize_execution(artifact)

    if parsed.execution_mode == "pk"
        return simulate(parsed.model_spec, parsed.grid, parsed.solver)
    end

    if parsed.pd_spec === nothing
        error(
            "execution_mode=$(parsed.execution_mode) requires pd_spec but artifact has none"
        )
    end

    if parsed.execution_mode == "pk_then_pd"
        return simulate_pkpd(parsed.model_spec, parsed.pd_spec, parsed.grid, parsed.solver)
    end

    if parsed.execution_mode == "pkpd_coupled"
        return simulate_pkpd_coupled(
            parsed.model_spec, parsed.pd_spec, parsed.grid, parsed.solver
        )
    end

    error("Unsupported execution_mode: $(parsed.execution_mode)")
end
