# Validate Real-World Theophylline Analysis Outputs

using JSON

const BASE_DIR = "docs/examples/use_cases/real_world_theophylline"

function validate()
    println("Validating theophylline analysis outputs...")
    all_pass = true

    # Check NCA results
    nca_path = joinpath(BASE_DIR, "output", "nca_results.json")
    if isfile(nca_path)
        nca = JSON.parsefile(nca_path)

        # Check subject count
        n_subjects = length(nca["individual"])
        if n_subjects != 12
            println("FAIL: Expected 12 subjects, got $n_subjects")
            all_pass = false
        else
            println("PASS: NCA results for 12 subjects")
        end

        # Check Cmax range (expected 6-11 mg/L)
        cmax_mean = nca["summary"]["Cmax_mean"]
        if 6.0 <= cmax_mean <= 11.0
            println("PASS: Cmax mean $(round(cmax_mean, digits=2)) in expected range [6-11] mg/L")
        else
            println("FAIL: Cmax mean $cmax_mean outside expected range")
            all_pass = false
        end
    else
        println("FAIL: NCA results file not found")
        all_pass = false
    end

    # Check population simulation
    pop_path = joinpath(BASE_DIR, "output", "population_simulation.json")
    if isfile(pop_path)
        println("PASS: Population simulation artifact exists")
    else
        println("FAIL: Population simulation file not found")
        all_pass = false
    end

    # Check VPC summary
    vpc_path = joinpath(BASE_DIR, "output", "vpc_summary.json")
    if isfile(vpc_path)
        println("PASS: VPC summary exists")
    else
        println("FAIL: VPC summary file not found")
        all_pass = false
    end

    # Check report
    report_path = joinpath(BASE_DIR, "output", "report.json")
    if isfile(report_path)
        println("PASS: Analysis report exists")
    else
        println("FAIL: Report file not found")
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
