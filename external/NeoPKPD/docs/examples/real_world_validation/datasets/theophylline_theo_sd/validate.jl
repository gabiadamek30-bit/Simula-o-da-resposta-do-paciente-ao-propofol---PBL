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

        # 208.6: Compare dose_unit_rule
        exp[i]["dose_unit_rule"] == act[i]["dose_unit_rule"] || error(
            "dose_unit_rule mismatch for id=$(exp[i]["id"])"
        )

        # Compare numeric fields with strict tolerance
        for k in ["wt", "dose_mg", "rmse_fixed", "rmse_wt_scaled"]
            e = Float64(exp[i][k])
            a = Float64(act[i][k])
            # Handle NaN comparison
            if isnan(e) && isnan(a)
                continue
            end
            abs(e - a) <= 1e-12 || error(
                "Metric mismatch $(k) for id=$(exp[i]["id"]). Expected=$(e) Actual=$(a)"
            )
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

    length(expC) == length(actC) || error("Length mismatch in conc series")

    for i in eachindex(expC)
        abs(expC[i] - actC[i]) <= 1e-12 || error("Conc mismatch at idx $(i)")
    end
end

function main()
    base = "docs/examples/real_world_validation/studies/theophylline_theo_sd"

    # Compare metrics
    compare_metrics(
        joinpath(base, "expected", "metrics.json"),
        joinpath(base, "output", "metrics.json"),
    )

    # Get subject IDs from expected metrics
    exp_metrics = JSON.parsefile(joinpath(base, "expected", "metrics.json"))

    # 208.6: Compare both scenarios for each subject
    for m in exp_metrics
        id = Int(m["id"])

        # Fixed scenario
        replay_and_compare(
            joinpath(base, "expected", "subj_$(id)_fixed.json"),
            joinpath(base, "output", "subj_$(id)_fixed.json"),
        )

        # WT-scaled scenario
        replay_and_compare(
            joinpath(base, "expected", "subj_$(id)_wt_scaled.json"),
            joinpath(base, "output", "subj_$(id)_wt_scaled.json"),
        )
    end

    println("Theophylline theo_sd real-data validation passed")
end

main()
