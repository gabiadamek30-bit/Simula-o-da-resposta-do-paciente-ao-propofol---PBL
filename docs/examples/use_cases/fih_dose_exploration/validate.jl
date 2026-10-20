using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD
using JSON

function compare_metrics(expected_path::String, actual_path::String)
    exp = JSON.parsefile(expected_path)
    act = JSON.parsefile(actual_path)

    if length(exp) != length(act)
        error("Metrics length mismatch")
    end

    for i in 1:length(exp)
        if Float64(exp[i]["dose_mg"]) != Float64(act[i]["dose_mg"])
            error("Dose mismatch at index $(i)")
        end

        # Very strict because everything is deterministic
        for k in ["cmax_mean", "auc0_24_mean"]
            e = Float64(exp[i][k])
            a = Float64(act[i][k])
            if abs(e - a) > 1e-12
                error("Metric mismatch for $(k) at dose $(exp[i]["dose_mg"]). Expected=$(e) Actual=$(a)")
            end
        end
    end
end

function replay_and_compare_artifact(expected_path::String, actual_path::String)
    exp_art = read_execution_json(expected_path)
    act_art = read_execution_json(actual_path)

    # enforce semantics fingerprint and compare via replay using existing golden runner style
    exp_res = replay_population_execution(exp_art)
    act_res = replay_population_execution(act_art)

    # compare mean concentration summary only, that is what drives decision making
    exp_mean = exp_res.summaries[:conc].mean
    act_mean = act_res.summaries[:conc].mean

    if length(exp_mean) != length(act_mean)
        error("Mean series length mismatch")
    end

    for i in eachindex(exp_mean)
        if abs(exp_mean[i] - act_mean[i]) > 1e-12
            error("Mean conc mismatch at index $(i)")
        end
    end
end

function main()
    base = "docs/examples/use_cases/fih_dose_exploration"

    # Compare metrics
    compare_metrics(
        joinpath(base, "expected", "metrics.json"),
        joinpath(base, "output", "metrics.json"),
    )

    # Compare each artifact by replaying and comparing core summary
    doses = [10, 30, 100, 300]
    for d in doses
        expected_art = joinpath(base, "expected", "pop_iv_$(d)mg.json")
        actual_art = joinpath(base, "output", "pop_iv_$(d)mg.json")
        replay_and_compare_artifact(expected_art, actual_art)
    end

    println("FIH dose exploration validation passed")
end

main()
