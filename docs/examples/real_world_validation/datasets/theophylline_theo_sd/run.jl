using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD
using JSON

function parse_maybe_float(s::String)
    ss = strip(s)
    if isempty(ss) || lowercase(ss) in ["na", "nan"]
        return nothing
    end
    return parse(Float64, ss)
end

DV = parse_maybe_float(fields[3]),

# Minimal strict CSV reader without extra deps.
function read_theo_csv(path::String)
    lines = readlines(path)
    isempty(lines) && error("Empty CSV: " * path)

    header = [replace(strip(h), "\"" => "") for h in split(strip(lines[1]), ",")]
    required = ["ID","TIME","DV","AMT","EVID","CMT","WT"]
    header != required && error("CSV header mismatch. Expected: $(required) Got: $(header)")

    rows = []
    for (i, ln) in enumerate(lines[2:end])
        s = strip(ln)
        isempty(s) && continue
        fields = split(s, ",")
        length(fields) != length(required) && error("Bad column count at line $(i+1)")

        # remove quotes if present
        fields = map(x -> replace(x, "\"" => ""), fields)

        push!(rows, (
            ID = parse(Int, fields[1]),
            TIME = parse(Float64, fields[2]),
            DV = parse(Float64, fields[3]),
            AMT = parse(Float64, fields[4]),
            EVID = parse(Int, fields[5]),
            CMT = parse(Int, fields[6]),
            WT = parse(Float64, fields[7]),
        ))
    end
    return rows
end

function group_by_id(rows)
    d = Dict{Int, Vector{Any}}()
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

# Fixed parameters from nlmixr example output:
# Ka ~ 1.59 1/hr, CL ~ 2.75 L/hr, V ~ 31.8 L. :contentReference[oaicite:5]{index=5}
const KA = 1.59
const CL = 2.75
const V  = 31.8

function simulate_subject(rows_for_id)
    # one oral dose event at TIME=0 with AMT>0 (EVID != 0)
    dose_rows = filter(r -> r.AMT > 0.0, rows_for_id)
    length(dose_rows) == 0 && error("No dose row found")
    length(dose_rows) > 1 && error("Multiple dose rows found; this study expects single-dose theo_sd day1")

    wt = dose_rows[1].WT
    amt = dose_rows[1].AMT

    # Explicit, deterministic rule:
    # If AMT looks like mg/kg (typically small numbers), convert to mg via WT.
    # Otherwise treat as mg already.
    dose_mg =
        amt < 50.0 ? amt * wt : amt

    dose_rule =
        amt < 50.0 ? "AMT_as_mg_per_kg_times_WT" : "AMT_as_mg"

    dose_time = dose_rows[1].TIME

    obs_rows = filter(r -> r.EVID == 0, rows_for_id)
    obs_times = unique([r.TIME for r in obs_rows])
    sort!(obs_times)

    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "theo_sd_subject_" * string(rows_for_id[1].ID),
        OneCompOralFirstOrderParams(KA, CL, V),
        [DoseEvent(dose_time, dose_amt)],
    )

    grid = SimGrid(minimum(obs_times), maximum(obs_times), obs_times)
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    res.metadata["dose_unit_rule"] = dose_rule
    res.metadata["raw_amt"] = amt
    res.metadata["wt"] = wt


    pred = res.observations[:conc]
    dv = [r.DV for r in obs_rows]
    t  = [r.TIME for r in obs_rows]

    # RMSE only over non-missing DV
    sse = 0.0
    n = 0
    for i in eachindex(dv)
        if dv[i] === nothing
            continue
        end
        e = pred_full[i] - Float64(dv[i])
        sse += e * e
        n += 1
    end

    rmse = n == 0 ? NaN : sqrt(sse / n)

    pred_full = [pred_at[tt] for tt in t]

    # RMSE against observed DV (purely descriptive, no fitting)
    sse = 0.0
    for i in eachindex(dv)
        e = pred_full[i] - dv[i]
        sse += e * e
    end
    rmse = sqrt(sse / length(dv))

    return spec, grid, solver, res, Dict(
        "id" => rows_for_id[1].ID,
        "wt" => rows_for_id[1].WT,
        "dose_mg" => dose_amt,
        "rmse" => rmse,
    )
end

function main()
    base = "docs/examples/real_world_validation"
    out_dir = joinpath(base, "studies/theophylline_theo_sd/output")
    mkpath(out_dir)

    data_path = joinpath(base, "datasets/theophylline_theo_sd/theo_sd.csv")
    rows = read_theo_csv(data_path)
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
