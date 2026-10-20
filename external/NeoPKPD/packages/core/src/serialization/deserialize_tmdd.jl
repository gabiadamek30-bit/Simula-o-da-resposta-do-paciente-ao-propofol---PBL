# =============================================================================
# TMDD Deserialization Functions
# =============================================================================
#
# Replay and deserialization support for TMDD model artifacts.
# =============================================================================

export replay_tmdd_execution, read_tmdd_execution_json

"""
    _deserialize_tmdd_kind(kind_dict)

Reconstruct TMDD model kind from serialized dictionary.
"""
function _deserialize_tmdd_kind(kind_dict::Dict)
    type_str = String(kind_dict["type"])

    # Parse approximation
    approx_str = get(kind_dict, "approximation", "QSS")
    approx = if approx_str == "FullTMDD"
        FullTMDD
    elseif approx_str == "QSS"
        QSS
    elseif approx_str == "QE"
        QE
    elseif approx_str == "RapidBinding"
        RapidBinding
    elseif approx_str == "WagnerBinding"
        WagnerBinding
    elseif approx_str == "IrreversibleBinding"
        IrreversibleBinding
    else
        QSS  # Default
    end

    # Parse route
    route_str = get(kind_dict, "route", "IVBolus")
    route = if route_str == "IVBolus"
        IVBolus
    elseif route_str == "IVInfusion"
        IVInfusion
    elseif route_str == "Subcutaneous"
        Subcutaneous
    else
        IVBolus  # Default
    end

    # Construct appropriate type
    if type_str == "OneCptTMDD"
        return OneCptTMDD(approx, route)
    elseif type_str == "TwoCptTMDD"
        return TwoCptTMDD(approx, route)
    elseif type_str == "TwoCptTMDDFcRn"
        return TwoCptTMDDFcRn(approx, route)
    elseif type_str == "TwoCptTMDDADA"
        return TwoCptTMDDADA(approx, route)
    elseif type_str == "TwoCptTMDDSoluble"
        return TwoCptTMDDSoluble(approx, route)
    elseif type_str == "TwoCptTMDDPeripheral"
        return TwoCptTMDDPeripheral(approx, route)
    else
        error("Unknown TMDD model type: $type_str")
    end
end

"""
    _deserialize_tmdd_params(kind, params_dict)

Reconstruct TMDD parameters from serialized dictionary.
"""
function _deserialize_tmdd_params(kind::OneCptTMDD, params_dict::Dict)
    return OneCptTMDDParams(
        Float64(params_dict["CL"]),
        Float64(params_dict["V"]),
        Float64(params_dict["KSS"]),
        Float64(params_dict["kint"]),
        Float64(params_dict["ksyn"]),
        Float64(params_dict["kdeg"]),
        Float64(params_dict["R0"]),
        Float64(get(params_dict, "ka", 0.0)),
        Float64(get(params_dict, "F", 1.0)),
    )
end

function _deserialize_tmdd_params(kind::TwoCptTMDD, params_dict::Dict)
    return TwoCptTMDDParams(
        Float64(params_dict["CL"]),
        Float64(params_dict["V1"]),
        Float64(params_dict["V2"]),
        Float64(params_dict["Q"]),
        Float64(params_dict["KSS"]),
        Float64(params_dict["kint"]),
        Float64(params_dict["ksyn"]),
        Float64(params_dict["kdeg"]),
        Float64(params_dict["R0"]),
        Float64(get(params_dict, "ka", 0.0)),
        Float64(get(params_dict, "F", 1.0)),
        Float64(get(params_dict, "Tlag", 0.0)),
    )
end

function _deserialize_tmdd_params(kind::TwoCptTMDDFcRn, params_dict::Dict)
    return TwoCptTMDDFcRnParams(
        Float64(params_dict["V1"]),
        Float64(params_dict["V2"]),
        Float64(params_dict["Q"]),
        Float64(params_dict["CLup"]),
        Float64(params_dict["FR"]),
        Float64(params_dict["KSS"]),
        Float64(params_dict["kint"]),
        Float64(params_dict["ksyn"]),
        Float64(params_dict["kdeg"]),
        Float64(params_dict["R0"]),
        Float64(get(params_dict, "ka", 0.0)),
        Float64(get(params_dict, "F", 1.0)),
        Float64(get(params_dict, "Tlag", 0.0)),
    )
end

function _deserialize_tmdd_params(kind::TwoCptTMDDADA, params_dict::Dict)
    return TwoCptTMDDADAParams(
        Float64(params_dict["CL"]),
        Float64(params_dict["V1"]),
        Float64(params_dict["V2"]),
        Float64(params_dict["Q"]),
        Float64(params_dict["KSS"]),
        Float64(params_dict["kint"]),
        Float64(params_dict["ksyn"]),
        Float64(params_dict["kdeg"]),
        Float64(params_dict["R0"]),
        Float64(params_dict["kADA_prod"]),
        Float64(params_dict["kADA_deg"]),
        Float64(params_dict["kon_ADA"]),
        Float64(params_dict["koff_ADA"]),
        Float64(params_dict["CL_complex"]),
        Float64(params_dict["T_onset"]),
        Float64(get(params_dict, "ka", 0.0)),
        Float64(get(params_dict, "F", 1.0)),
    )
end

# Generic fallback for other TMDD types
function _deserialize_tmdd_params(kind::TMDDModelKind, params_dict::Dict)
    # Try TwoCptTMDD params as fallback
    return TwoCptTMDDParams(
        Float64(params_dict["CL"]),
        Float64(params_dict["V1"]),
        Float64(params_dict["V2"]),
        Float64(params_dict["Q"]),
        Float64(params_dict["KSS"]),
        Float64(params_dict["kint"]),
        Float64(params_dict["ksyn"]),
        Float64(params_dict["kdeg"]),
        Float64(params_dict["R0"]),
        Float64(get(params_dict, "ka", 0.0)),
        Float64(get(params_dict, "F", 1.0)),
        Float64(get(params_dict, "Tlag", 0.0)),
    )
end

"""
    _deserialize_tmdd_spec(spec_dict)

Reconstruct TMDDSpec from serialized dictionary.
"""
function _deserialize_tmdd_spec(spec_dict::Dict)
    kind = _deserialize_tmdd_kind(Dict{String,Any}(spec_dict["kind"]))
    params = _deserialize_tmdd_params(kind, Dict{String,Any}(spec_dict["params"]))

    doses = DoseEvent[]
    for d in spec_dict["doses"]
        dd = Dict{String,Any}(d)
        push!(doses, DoseEvent(
            Float64(dd["time"]),
            Float64(dd["amount"]),
            Float64(get(dd, "duration", 0.0)),
        ))
    end

    target_units = Symbol(get(spec_dict, "target_units", "nM"))
    drug_units = Symbol(get(spec_dict, "drug_units", "mg_L"))

    return TMDDSpec(kind, String(spec_dict["name"]), params, doses, target_units, drug_units)
end

"""
    replay_tmdd_execution(artifact::Dict) -> TMDDSimResult

Replay a TMDD simulation from a serialized artifact.

# Arguments
- `artifact::Dict`: Serialized TMDD execution artifact

# Returns
- `TMDDSimResult`: Result of replaying the simulation
"""
function replay_tmdd_execution(artifact::Dict)
    # Reconstruct spec
    spec_dict = Dict{String,Any}(artifact["tmdd_spec"])
    spec = _deserialize_tmdd_spec(spec_dict)

    # Reconstruct grid
    grid_dict = Dict{String,Any}(artifact["grid"])
    grid = SimGrid(
        Float64(grid_dict["t0"]),
        Float64(grid_dict["t1"]),
        [Float64(x) for x in grid_dict["saveat"]],
    )

    # Reconstruct solver
    solver_dict = Dict{String,Any}(artifact["solver"])
    solver = SolverSpec(
        Symbol(solver_dict["alg"]),
        Float64(solver_dict["reltol"]),
        Float64(solver_dict["abstol"]),
        Int(solver_dict["maxiters"]),
    )

    # Run simulation
    return solve_tmdd(spec, grid, solver)
end

"""
    read_tmdd_execution_json(path::AbstractString) -> Dict

Read a TMDD execution artifact from a JSON file.
"""
function read_tmdd_execution_json(path::AbstractString)
    return JSON.parsefile(path)
end
