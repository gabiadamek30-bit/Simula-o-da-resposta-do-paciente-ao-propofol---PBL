# Validate NONMEM Migration Outputs

using JSON

const BASE_DIR = "docs/examples/use_cases/nonmem_migration"

function validate()
    println("Validating NONMEM migration outputs...")
    all_pass = true

    # Check migration report
    report_path = joinpath(BASE_DIR, "output", "migration_report.json")
    if isfile(report_path)
        report = JSON.parsefile(report_path)

        # Check migration status
        if report["status"] == "MIGRATION_SUCCESSFUL"
            println("PASS: Migration status is SUCCESSFUL")
        else
            println("FAIL: Migration status is $(report["status"])")
            all_pass = false
        end

        # Check model type mapping
        if report["target"]["model_type"] == "OneCompOralFirstOrder"
            println("PASS: Model type correctly mapped (ADVAN2 -> OneCompOralFirstOrder)")
        else
            println("FAIL: Unexpected model type: $(report["target"]["model_type"])")
            all_pass = false
        end

        # Check parameters
        params = report["target"]["parameters"]
        if params["Ka"] == 1.5 && params["CL"] == 2.8 && params["V"] == 35.0
            println("PASS: Parameters correctly extracted (Ka=1.5, CL=2.8, V=35)")
        else
            println("FAIL: Parameter mismatch")
            all_pass = false
        end

        # Check IIV values
        iiv = report["target"]["iiv"]
        if abs(iiv["CL_cv"] - 30.0) < 1.0  # sqrt(0.09) * 100 = 30%
            println("PASS: IIV values correctly converted")
        else
            println("FAIL: IIV conversion error")
            all_pass = false
        end

        # Check validation metrics are reasonable
        val = report["validation"]
        if 5.0 <= val["cmax_mean"] <= 12.0
            println("PASS: Cmax $(val["cmax_mean"]) mg/L in expected range")
        else
            println("FAIL: Cmax outside expected range")
            all_pass = false
        end
    else
        println("FAIL: Migration report not found")
        all_pass = false
    end

    # Check population artifact
    pop_path = joinpath(BASE_DIR, "output", "migrated_population.json")
    if isfile(pop_path)
        println("PASS: Population artifact exists")
    else
        println("FAIL: Population artifact not found")
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
