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

function parse_int(s::AbstractString)
    return parse(Int, strip(replace(String(s), "\"" => "")))
end

function parse_float(s::AbstractString)
    return parse(Float64, strip(replace(String(s), "\"" => "")))
end

function read_csv_strict(path::String)
    lines = readlines(path)
    isempty(lines) && error("Empty CSV: " * path)

    header = [replace(strip(h), "\"" => "") for h in split(strip(lines[1]), ",")]
    required = ["ID", "TIME", "DV", "AMT", "EVID", "CMT", "WT"]
    header != required && error("CSV header mismatch. Expected $(required) Got: $(header)")

    rows = Any[]
    for (li, ln) in enumerate(lines[2:end])
        s = strip(ln)
        isempty(s) && continue
        f = split(s, ",")
        length(f) != length(required) && error("Bad column count at line $(li+1)")

        push!(rows, (
            ID = parse_int(f[1]),
            TIME = parse_float(f[2]),
            DV = parse_maybe_float(f[3]),
            AMT = parse_float(f[4]),
            EVID = parse_int(f[5]),
            CMT = parse_int(f[6]),
            WT = parse_float(f[7]),
        ))
    end
    return rows
end

function group_by_id(rows)
    d = Dict{Int,Vector{Any}}()
    for r in rows
        if !haskey(d, r.ID)
            d[r.ID] = Any[]
        end
        push!(d[r.ID], r)
    end
    for (_, v) in d
        sort!(v, by = x -> x.TIME)
    end
    return d
end

# Fixed params for this validation study (no fitting)
const KA = 1.59
const CL_REF = 2.75
const V_REF = 31.8

function dose_mg_from_row(r)
    amt = r.AMT
    wt = r.WT
    if amt < 50.0
        return amt * wt, "AMT_as_mg_per_kg_times_WT"
    end
    return amt, "AMT_as_mg"
end

function simulate_subject(rows_for_id)
    # Dosing events: any record with AMT > 0
    dose_rows = filter(r -> r.AMT > 0.0, rows_for_id)
    length(dose_rows) == 0 && error("No dose rows for subject")

    doses = DoseEvent[]
    rules = String[]
    for r in dose_rows
        dm, rule = dose_mg_from_row(r)
        push!(doses, DoseEvent(r.TIME, dm))
        push!(rules, rule)
    end

    # Observation records: use DV presence as the anchor
    obs_rows = filter(r -> r.DV !== nothing, rows_for_id)
    length(obs_rows) == 0 && error("No DV observations for subject")

    obs_times = unique([r.TIME for r in obs_rows])
    sort!(obs_times)

    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "theo_md_subj_" * string(rows_for_id[1].ID),
        OneCompOralFirstOrderParams(KA, CL_REF, V_REF),
        doses,
    )

    grid = SimGrid(minimum(obs_times), maximum(obs_times), obs_times)
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    # Metadata: dose unit rule summary
    unique_rules = unique(rules)
    dose_rule = length(unique_rules) == 1 ? unique_rules[1] : "MIXED_RULES"

    res.metadata["dose_unit_rule"] = dose_rule
    res.metadata["dose_event_count"] = length(doses)

    # RMSE over DV observations (DV non-missing by construction here)
    pred = res.observations[:conc]
    pred_at = Dict{Float64,Float64}()
    for (i, tt) in enumerate(res.t)
        pred_at[tt] = pred[i]
    end

    sse = 0.0
    n = 0
    for r in obs_rows
        y = Float64(r.DV)
        yhat = pred_at[r.TIME]
        e = yhat - y
        sse += e * e
        n += 1
    end

    rmse = sqrt(sse / n)

    metrics = Dict(
        "id" => rows_for_id[1].ID,
        "wt" => rows_for_id[1].WT,
        "dose_unit_rule" => dose_rule,
        "dose_events" => length(doses),
        "rmse" => rmse,
    )

    return spec, grid, solver, res, metrics
end

function main()
    base = "docs/examples/real_world_validation"
    out_dir = joinpath(base, "studies/theophylline_theo_md/output")
    mkpath(out_dir)

    data_path = joinpath(base, "datasets/theophylline_theo_md/theo_md.csv")
    rows = read_csv_strict(data_path)
    byid = group_by_id(rows)

    metrics = Any[]

    for id in sort(collect(keys(byid)))
        spec, grid, solver, res, m = simulate_subject(byid[id])
        push!(metrics, m)

        out_art = joinpath(out_dir, "subj_$(id).json")
        write_execution_json(out_art; model_spec = spec, grid = grid, solver = solver, result = res)
        println("Wrote artifact: " * out_art)
    end

    metrics_path = joinpath(out_dir, "metrics.json")
    open(metrics_path, "w") do io
        JSON.print(io, metrics, 2)
    end
    println("Wrote metrics: " * metrics_path)
end

main()
