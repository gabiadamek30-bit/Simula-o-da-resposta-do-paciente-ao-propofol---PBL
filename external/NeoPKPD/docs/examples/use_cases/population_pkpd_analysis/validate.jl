# Validate Population PKPD Analysis Outputs

using JSON

const BASE_DIR = "docs/examples/use_cases/population_pkpd_analysis"

function validate()
    println("Validating population PKPD analysis outputs...")
    all_pass = true

    # Check summary
    summary_path = joinpath(BASE_DIR, "output", "summary.json")
    if isfile(summary_path)
        summary = JSON.parsefile(summary_path)

        # Check subject count
        if summary["n"] == 50
            println("PASS: 50 subjects analyzed")
        else
            println("FAIL: Expected 50 subjects, got $(summary["n"])")
            all_pass = false
        end

        # Check PK Cmax range
        cmax = summary["pk"]["Cmax_mean"]
        if 1.5 <= cmax <= 5.0
            println("PASS: Cmax $(round(cmax, digits=2)) mg/L in expected range")
        else
            println("FAIL: Cmax $cmax outside expected range [1.5-5.0]")
            all_pass = false
        end

        # Check max inhibition
        max_inh = summary["pd"]["max_inhibition_mean"]
        if 30 <= max_inh <= 70
            println("PASS: Max inhibition $(round(max_inh, digits=1))% in expected range")
        else
            println("FAIL: Max inhibition $max_inh outside expected range [30-70%]")
            all_pass = false
        end
    else
        println("FAIL: Summary file not found")
        all_pass = false
    end

    # Check PK results
    pk_path = joinpath(BASE_DIR, "output", "pk_results.json")
    if isfile(pk_path)
        pk = JSON.parsefile(pk_path)
        if length(pk) == 50
            println("PASS: PK results for 50 subjects")
        else
            println("FAIL: Expected 50 PK results, got $(length(pk))")
            all_pass = false
        end
    else
        println("FAIL: PK results file not found")
        all_pass = false
    end

    # Check PD results
    pd_path = joinpath(BASE_DIR, "output", "pd_results.json")
    if isfile(pd_path)
        pd = JSON.parsefile(pd_path)
        if length(pd) == 50
            println("PASS: PD results for 50 subjects")
        else
            println("FAIL: Expected 50 PD results, got $(length(pd))")
            all_pass = false
        end
    else
        println("FAIL: PD results file not found")
        all_pass = false
    end

    println()
    if all_pass
        println("All validations PASSED")
        exit(0)
    else
        println("Some validations FAILED")
        exit(1)
    end
end

validate()
