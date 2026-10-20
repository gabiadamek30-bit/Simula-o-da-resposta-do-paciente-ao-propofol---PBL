# Validation report generation for IQ/OQ/PQ protocols
# Supports FDA 21 CFR Part 11 and GAMP 5 validation requirements

export ValidationReport, ValidationTestResult, ValidationStatus
export VALIDATION_PASSED, VALIDATION_FAILED, VALIDATION_SKIPPED
export generate_validation_report, write_validation_report_markdown

"""
    ValidationStatus

Status of a validation test.
"""
@enum ValidationStatus begin
    VALIDATION_PASSED = 1
    VALIDATION_FAILED = 2
    VALIDATION_SKIPPED = 3
end

"""
    ValidationTestResult

Result of a single validation test.

# Fields
- `test_id::String`: Unique identifier for the test
- `test_name::String`: Human-readable test name
- `test_description::String`: Description of what the test validates
- `status::ValidationStatus`: Pass/Fail/Skip status
- `expected_value::Any`: Expected result
- `actual_value::Any`: Actual result
- `error_message::Union{Nothing,String}`: Error message if failed
- `execution_time_ms::Float64`: Test execution time in milliseconds
"""
struct ValidationTestResult
    test_id::String
    test_name::String
    test_description::String
    status::ValidationStatus
    expected_value::Any
    actual_value::Any
    error_message::Union{Nothing,String}
    execution_time_ms::Float64
end

"""
    ValidationReport

Complete validation report for IQ, OQ, or PQ protocols.

# Fields
- `report_type::String`: "IQ", "OQ", or "PQ"
- `report_id::String`: Unique report identifier (UUID)
- `report_timestamp::String`: ISO 8601 UTC timestamp
- `environment::EnvironmentSnapshot`: Environment at validation time
- `tests::Vector{ValidationTestResult}`: All test results
- `total_tests::Int`: Total number of tests
- `passed_tests::Int`: Number of passed tests
- `failed_tests::Int`: Number of failed tests
- `skipped_tests::Int`: Number of skipped tests
- `overall_status::ValidationStatus`: Overall validation status
- `execution_time_ms::Float64`: Total execution time
- `notes::Vector{String}`: Additional notes or observations
"""
struct ValidationReport
    report_type::String
    report_id::String
    report_timestamp::String
    environment::EnvironmentSnapshot
    tests::Vector{ValidationTestResult}
    total_tests::Int
    passed_tests::Int
    failed_tests::Int
    skipped_tests::Int
    overall_status::ValidationStatus
    execution_time_ms::Float64
    notes::Vector{String}
end

"""
    _status_to_string(status::ValidationStatus) -> String

Convert ValidationStatus to display string.
"""
function _status_to_string(status::ValidationStatus)::String
    if status == VALIDATION_PASSED
        return "PASSED"
    elseif status == VALIDATION_FAILED
        return "FAILED"
    else
        return "SKIPPED"
    end
end

"""
    _status_emoji(status::ValidationStatus) -> String

Get emoji for validation status.
"""
function _status_emoji(status::ValidationStatus)::String
    if status == VALIDATION_PASSED
        return "[PASS]"
    elseif status == VALIDATION_FAILED
        return "[FAIL]"
    else
        return "[SKIP]"
    end
end

"""
    generate_validation_report(;
        report_type::String,
        tests::Vector{ValidationTestResult},
        notes::Vector{String} = String[],
        execution_time_ms::Float64 = 0.0
    ) -> ValidationReport

Generate a validation report from test results.
"""
function generate_validation_report(;
    report_type::String,
    tests::Vector{ValidationTestResult},
    notes::Vector{String}=String[],
    execution_time_ms::Float64=0.0
)::ValidationReport
    passed = count(t -> t.status == VALIDATION_PASSED, tests)
    failed = count(t -> t.status == VALIDATION_FAILED, tests)
    skipped = count(t -> t.status == VALIDATION_SKIPPED, tests)

    overall_status = if failed > 0
        VALIDATION_FAILED
    elseif passed == 0 && skipped == length(tests)
        VALIDATION_SKIPPED
    else
        VALIDATION_PASSED
    end

    return ValidationReport(
        report_type,
        string(uuid4()),
        _utc_timestamp(),
        capture_environment(),
        tests,
        length(tests),
        passed,
        failed,
        skipped,
        overall_status,
        execution_time_ms,
        notes
    )
end

"""
    write_validation_report_markdown(report::ValidationReport, path::String)

Write a validation report to a Markdown file.
"""
function write_validation_report_markdown(report::ValidationReport, path::String)
    open(path, "w") do io
        # Header
        println(io, "# $(report.report_type) Validation Report")
        println(io)
        println(io, "**Report ID:** `$(report.report_id)`")
        println(io)
        println(io, "**Generated:** $(report.report_timestamp)")
        println(io)
        println(io, "**Overall Status:** $(_status_emoji(report.overall_status)) $(_status_to_string(report.overall_status))")
        println(io)

        # Summary
        println(io, "## Summary")
        println(io)
        println(io, "| Metric | Value |")
        println(io, "|--------|-------|")
        println(io, "| Total Tests | $(report.total_tests) |")
        println(io, "| Passed | $(report.passed_tests) |")
        println(io, "| Failed | $(report.failed_tests) |")
        println(io, "| Skipped | $(report.skipped_tests) |")
        println(io, "| Execution Time | $(round(report.execution_time_ms, digits=2)) ms |")
        println(io)

        # Environment
        println(io, "## Environment")
        println(io)
        println(io, "| Component | Version |")
        println(io, "|-----------|---------|")
        println(io, "| Julia | $(report.environment.julia_version) |")
        println(io, "| NeoPKPD | $(report.environment.neopkpd_version) |")
        println(io, "| Schema | $(report.environment.artifact_schema_version) |")
        println(io, "| DifferentialEquations | $(report.environment.differentialequations_version) |")
        println(io, "| OS | $(report.environment.os) $(report.environment.os_version) |")
        println(io, "| Architecture | $(report.environment.architecture) |")
        if report.environment.git_commit !== nothing
            println(io, "| Git Commit | $(report.environment.git_commit[1:min(8, length(report.environment.git_commit))]) |")
        end
        println(io)

        # Test Results
        println(io, "## Test Results")
        println(io)
        for test in report.tests
            println(io, "### $(_status_emoji(test.status)) $(test.test_id): $(test.test_name)")
            println(io)
            println(io, "**Description:** $(test.test_description)")
            println(io)
            println(io, "**Status:** $(_status_to_string(test.status))")
            println(io)
            if test.status == VALIDATION_FAILED && test.error_message !== nothing
                println(io, "**Error:** $(test.error_message)")
                println(io)
            end
            println(io, "**Execution Time:** $(round(test.execution_time_ms, digits=2)) ms")
            println(io)
        end

        # Notes
        if !isempty(report.notes)
            println(io, "## Notes")
            println(io)
            for note in report.notes
                println(io, "- $(note)")
            end
            println(io)
        end

        # Footer
        println(io, "---")
        println(io)
        println(io, "*Generated by NeoPKPD $(report.environment.neopkpd_version) Compliance Module*")
    end

    return path
end

"""
    serialize_validation_report(report::ValidationReport) -> Dict{String, Any}

Serialize a validation report to a dictionary for JSON encoding.
"""
function serialize_validation_report(report::ValidationReport)::Dict{String,Any}
    tests_serialized = [
        Dict{String,Any}(
            "test_id" => t.test_id,
            "test_name" => t.test_name,
            "test_description" => t.test_description,
            "status" => _status_to_string(t.status),
            "expected_value" => t.expected_value,
            "actual_value" => t.actual_value,
            "error_message" => t.error_message,
            "execution_time_ms" => t.execution_time_ms
        )
        for t in report.tests
    ]

    return Dict{String,Any}(
        "report_type" => report.report_type,
        "report_id" => report.report_id,
        "report_timestamp" => report.report_timestamp,
        "environment" => serialize_environment(report.environment),
        "tests" => tests_serialized,
        "total_tests" => report.total_tests,
        "passed_tests" => report.passed_tests,
        "failed_tests" => report.failed_tests,
        "skipped_tests" => report.skipped_tests,
        "overall_status" => _status_to_string(report.overall_status),
        "execution_time_ms" => report.execution_time_ms,
        "notes" => report.notes
    )
end

export serialize_validation_report
