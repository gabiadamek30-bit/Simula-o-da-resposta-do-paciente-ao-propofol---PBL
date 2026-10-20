# Validate Bioequivalence Study Outputs

using JSON

const BASE_DIR = "docs/examples/use_cases/bioequivalence_study"

function validate()
    println("Validating bioequivalence study outputs...")
    all_pass = true

    # Check BE statistics
    stats_path = joinpath(BASE_DIR, "output", "be_statistics.json")
    if isfile(stats_path)
        stats = JSON.parsefile(stats_path)

        # Check subject count
        if stats["n_subjects"] == 24
            println("PASS: 24 subjects evaluated")
        else
            println("FAIL: Expected 24 subjects, got $(stats["n_subjects"])")
            all_pass = false
        end

        # Check GMR Cmax range (expected ~0.93 due to Ka difference)
        gmr_cmax = stats["Cmax"]["GMR"]
        if 0.85 <= gmr_cmax <= 1.00
            println("PASS: Cmax GMR $(round(gmr_cmax, digits=3)) in expected range")
        else
            println("FAIL: Cmax GMR $gmr_cmax outside expected range [0.85-1.00]")
            all_pass = false
        end

        # Check AUC GMR (should be ~1.0 since only Ka differs)
        gmr_auc = stats["AUC"]["GMR"]
        if 0.95 <= gmr_auc <= 1.05
            println("PASS: AUC GMR $(round(gmr_auc, digits=3)) in expected range")
        else
            println("FAIL: AUC GMR $gmr_auc outside expected range [0.95-1.05]")
            all_pass = false
        end

        # Check BE conclusion
        if stats["overall_BE"] == true
            println("PASS: BE conclusion is BIOEQUIVALENT")
        else
            println("WARN: BE conclusion is NOT BIOEQUIVALENT (may vary by seed)")
        end
    else
        println("FAIL: BE statistics file not found")
        all_pass = false
    end

    # Check individual results
    ind_path = joinpath(BASE_DIR, "output", "individual_results.json")
    if isfile(ind_path)
        results = JSON.parsefile(ind_path)
        if length(results) == 48  # 24 subjects * 2 periods
            println("PASS: Individual results for 48 observations (24 subjects x 2 periods)")
        else
            println("FAIL: Expected 48 observations, got $(length(results))")
            all_pass = false
        end
    else
        println("FAIL: Individual results file not found")
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
