using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD
using JSON

function parse_maybe_float(s::AbstractString)
    ss = strip(replace(String(s), "\"" => ""))
    if isempty(ss) || lowercase(ss) in ["na", "nan", "."]
        return nothing
    end
    return parse(Float64, ss)
end

parse_int(s::AbstractString) = parse(Int, strip(replace(String(s), "\"" => "")))
parse_float(s::AbstractString) = parse(Float64, strip(replace(String(s), "\"" => "")))
parse_str(s::AbstractString) = strip(replace(String(s), "\"" => ""))

function read_warfarin_csv(path::String)
    lines = readlines(path)
    isempty(lines) && error("Empty CSV: " * path)

    header = [replace(strip(h), "\"" => "") for h in split(strip(lines[1]), ",")]
    required = ["id", "time", "amt", "dv", "dvid", "evid", "wt", "age", "sex"]
    header != required && error("CSV header mismatch. Expected $(required) Got: $(header)")

    rows = Any[]
    for (li, ln) in enumerate(lines[2:end])
        s = strip(ln)
        isempty(s) && continue
        f = split(s, ",")
        length(f) != length(required) && error("Bad column count at line $(li+1)")

        push!(rows, (
            id = parse_int(f[1]),
            time = parse_float(f[2]),
            amt = parse_float(f[3]),
            dv = parse_maybe_float(f[4]),
            dvid = parse_str(f[5]),
            evid = parse_int(f[6]),
            wt = parse_float(f[7]),
            age = parse_float(f[8]),
            sex = parse_str(f[9]),
        ))
    end
    return rows
end

function group_by_id(rows)
    d = Dict{Int,Vector{Any}}()
    for r in rows
        if !haskey(d, r.id)
            d[r.id] = Any[]
        end
        push!(d[r.id], r)
    end
    for (_, v) in d
        sort!(v, by = x -> x.time)
    end
    return d
end

# Fixed parameters for deterministic validation, not fitted.
# These are not claims about true warfarin PK/PD.
const PK_KA = 1.2
const PK_CL = 3.0
const PK_V = 35.0

# Indirect response turnover PD params
# Kin, Kout, R0, Emax, EC50
const PD_KIN = 10.0
const PD_KOUT = 0.5
const PD_R0 = PD_KIN / PD_KOUT
const PD_EMAX = 0.8
const PD_EC50 = 2.0

function build_specs(doses::Vector{DoseEvent}, subj_id::Int)
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "warfarin_pk_subj_" * string(subj_id),
        OneCompOralFirstOrderParams(PK_KA, PK_CL, PK_V),
        doses,
    )

    pd = PDSpec(
        IndirectResponseTurnover(),
        "warfarin_pd_turnover",
        IndirectResponseTurnoverParams(PD_KIN, PD_KOUT, PD_R0, PD_EMAX, PD_EC50),
        :conc,
        :response,
    )

    return pk, pd
end

function simulate_subject(rows_for_id)
    subj_id = rows_for_id[1].id

    dose_rows = filter(r -> r.amt > 0.0, rows_for_id)
    length(dose_rows) == 0 && error("No dosing rows for subject $(subj_id)")

    doses = DoseEvent[]
    for r in dose_rows
        push!(doses, DoseEvent(r.time, r.amt))
    end

    obs_rows = filter(r -> r.dv !== nothing, rows_for_id)
    length(obs_rows) == 0 && error("No DV observations for subject $(subj_id)")

    obs_times = unique([r.time for r in obs_rows])
    sort!(obs_times)

    pk, pd = build_specs(doses, subj_id)

    grid = SimGrid(minimum(obs_times), maximum(obs_times), obs_times)
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    # Coupled PKPD simulation
    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    # Map predictions to each observation row by time
    pred_conc = res.observations[:conc]
    pred_resp = res.observations[:response]
    pred_at_conc = Dict{Float64,Float64}()
    pred_at_resp = Dict{Float64,Float64}()

    for (i, tt) in enumerate(res.t)
        pred_at_conc[tt] = pred_conc[i]
        pred_at_resp[tt] = pred_resp[i]
    end

    # Dataset uses dv + dvid to represent endpoints. We compare:
    # dvid "cp" or "PK" -> concentration
    # all others -> response
    sse_pk = 0.0
    n_pk = 0
    sse_pd = 0.0
    n_pd = 0

    for r in obs_rows
        y = Float64(r.dv)
        dvid = lowercase(String(r.dvid))

        if dvid == "cp" || dvid == "pk" || dvid == "dose"
            yhat = pred_at_conc[r.time]
            e = yhat - y
            sse_pk += e * e
            n_pk += 1
        else
            yhat = pred_at_resp[r.time]
            e = yhat - y
            sse_pd += e * e
            n_pd += 1
        end
    end

    rmse_pk = n_pk == 0 ? NaN : sqrt(sse_pk / n_pk)
    rmse_pd = n_pd == 0 ? NaN : sqrt(sse_pd / n_pd)

    res.metadata["dataset"] = "nlmixr2data::warfarin"
    res.metadata["subject_id"] = subj_id
    res.metadata["rmse_pk"] = rmse_pk
    res.metadata["rmse_pd"] = rmse_pd

    metrics = Dict(
        "id" => subj_id,
        "wt" => rows_for_id[1].wt,
        "age" => rows_for_id[1].age,
        "sex" => rows_for_id[1].sex,
        "dose_events" => length(doses),
        "n_pk" => n_pk,
        "n_pd" => n_pd,
        "rmse_pk" => rmse_pk,
        "rmse_pd" => rmse_pd,
    )

    return pk, pd, grid, solver, res, metrics
end

function main()
    base = "docs/examples/real_world_validation"
    out_dir = joinpath(base, "studies/warfarin_pkpd/output")
    mkpath(out_dir)

    data_path = joinpath(base, "datasets/warfarin_nlmixr2data/warfarin.csv")
    rows = read_warfarin_csv(data_path)
    byid = group_by_id(rows)

    metrics = Any[]

    for id in sort(collect(keys(byid)))
        pk, pd, grid, solver, res, m = simulate_subject(byid[id])
        push!(metrics, m)

        out_art = joinpath(out_dir, "subj_$(id).json")
        write_execution_json(out_art; model_spec = pk, grid = grid, solver = solver, result = res, pd_spec = pd)
        println("Wrote artifact: " * out_art)
    end

    metrics_path = joinpath(out_dir, "metrics.json")
    open(metrics_path, "w") do io
        JSON.print(io, metrics, 2)
    end
    println("Wrote metrics: " * metrics_path)
end

main()
