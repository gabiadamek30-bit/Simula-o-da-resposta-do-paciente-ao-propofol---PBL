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

    # Sort by scenario and iov flag so ordering is stable
    function key(m)
        return (String(m["scenario"]), Bool(m["with_iov"]))
    end

    sort!(exp, by = key)
    sort!(act, by = key)

    for i in 1:length(exp)
        if key(exp[i]) != key(act[i])
            error("Metric key mismatch at index $(i)")
        end

        for k in ["emin_mean", "time_below_80pct_mean", "suppression_auc_mean"]
            e = Float64(exp[i][k])
            a = Float64(act[i][k])
            if abs(e - a) > 1e-12
                error("Metric mismatch for $(k) in $(exp[i]["scenario"]) with_iov=$(exp[i]["with_iov"]). Expected=$(e) Actual=$(a)")
            end
        end
    end
end

function replay_and_compare_mean_response(expected_path::String, actual_path::String)
    exp_art = read_execution_json(expected_path)
    act_art = read_execution_json(actual_path)

    exp_res = replay_population_execution(exp_art)
    act_res = replay_population_execution(act_art)

    exp_mean = exp_res.summaries[:response].mean
    act_mean = act_res.summaries[:response].mean

    if length(exp_mean) != length(act_mean)
        error("Mean response length mismatch")
    end

    for i in eachindex(exp_mean)
        if abs(exp_mean[i] - act_mean[i]) > 1e-12
            error("Mean response mismatch at index $(i)")
        end
    end
end

function main()
    base = "docs/examples/use_cases/pkpd_biomarker_turnover"

    compare_metrics(
        joinpath(base, "expected", "metrics.json"),
        joinpath(base, "output", "metrics.json"),
    )

    files = [
        "qd_100mg_no_iov.json",
        "qd_100mg_iov.json",
        "bid_50mg_no_iov.json",
        "bid_50mg_iov.json",
    ]

    for f in files
        replay_and_compare_mean_response(
            joinpath(base, "expected", f),
            joinpath(base, "output", f),
        )
    end

    println("PKPD biomarker turnover validation passed")
end

main()
