using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD
using JSON

function compare_metrics(expected_path::String, actual_path::String)
    exp = JSON.parsefile(expected_path)
    act = JSON.parsefile(actual_path)

    length(exp) == length(act) || error("Metrics length mismatch")

    sort!(exp, by = x -> Int(x["id"]))
    sort!(act, by = x -> Int(x["id"]))

    for i in eachindex(exp)
        Int(exp[i]["id"]) == Int(act[i]["id"]) || error("ID mismatch at index $(i)")

        for k in ["wt", "age", "rmse_pk", "rmse_pd"]
            e = Float64(exp[i][k])
            a = Float64(act[i][k])
            # Handle NaN comparison
            if isnan(e) && isnan(a)
                continue
            end
            abs(e - a) <= 1e-12 || error("Metric mismatch $(k) for id=$(exp[i]["id"]). Expected=$(e) Actual=$(a)")
        end

        for k in ["sex"]
            String(exp[i][k]) == String(act[i][k]) || error("Metric mismatch $(k) for id=$(exp[i]["id"])")
        end

        for k in ["dose_events", "n_pk", "n_pd"]
            Int(exp[i][k]) == Int(act[i][k]) || error("Metric mismatch $(k) for id=$(exp[i]["id"])")
        end
    end
end

function replay_and_compare(expected_art::String, actual_art::String)
    expA = read_execution_json(expected_art)
    actA = read_execution_json(actual_art)

    expR = replay_execution(expA)
    actR = replay_execution(actA)

    expC = expR.observations[:conc]
    actC = actR.observations[:conc]
    length(expC) == length(actC) || error("Conc length mismatch")
    for i in eachindex(expC)
        abs(expC[i] - actC[i]) <= 1e-12 || error("Conc mismatch at idx $(i)")
    end

    expY = expR.observations[:response]
    actY = actR.observations[:response]
    length(expY) == length(actY) || error("Response length mismatch")
    for i in eachindex(expY)
        abs(expY[i] - actY[i]) <= 1e-12 || error("Response mismatch at idx $(i)")
    end
end

function main()
    base = "docs/examples/real_world_validation/studies/warfarin_pkpd"

    compare_metrics(
        joinpath(base, "expected", "metrics.json"),
        joinpath(base, "output", "metrics.json"),
    )

    exp_metrics = JSON.parsefile(joinpath(base, "expected", "metrics.json"))
    for m in exp_metrics
        id = Int(m["id"])
        replay_and_compare(
            joinpath(base, "expected", "subj_$(id).json"),
            joinpath(base, "output", "subj_$(id).json"),
        )
    end

    println("Warfarin PKPD real-data validation passed")
end

main()
