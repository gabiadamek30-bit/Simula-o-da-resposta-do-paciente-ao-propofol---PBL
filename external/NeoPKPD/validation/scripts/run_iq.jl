# Installation Qualification (IQ) Automated Validation Script
# FDA 21 CFR Part 11 Compliant
#
# Usage:
#   julia validation/scripts/run_iq.jl
#   julia validation/scripts/run_iq.jl --output validation/reports/iq_report.md

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "packages", "core"))

using NeoPKPD
using Test

# Parse command line arguments
output_path = "validation/reports/iq_report.md"
for i in 1:length(ARGS)
    if ARGS[i] == "--output" && i < length(ARGS)
        global output_path = ARGS[i+1]
    end
end

println("=" ^ 60)
println("NeoPKPD Installation Qualification (IQ)")
println("=" ^ 60)
println()

# Collect test results
test_results = ValidationTestResult[]
start_time = time()

# IQ-001: Julia Version Verification
println("Running IQ-001: Julia Version Verification")
t1 = time()
try
    julia_version = VERSION
    min_version = v"1.9.0"
    passed = julia_version >= min_version
    push!(test_results, ValidationTestResult(
        "IQ-001",
        "Julia Version Verification",
        "Verify Julia installation meets minimum version requirements (>= 1.9.0)",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        ">= $min_version",
        string(julia_version),
        passed ? nothing : "Julia version $julia_version is below minimum $min_version",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED", " (", julia_version, ")")
catch e
    push!(test_results, ValidationTestResult(
        "IQ-001", "Julia Version Verification",
        "Verify Julia installation meets minimum version requirements",
        VALIDATION_FAILED, ">= 1.9.0", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# IQ-002: NeoPKPD Package Installation
println("Running IQ-002: NeoPKPD Package Installation")
t1 = time()
try
    version = NEOPKPD_VERSION
    passed = !isempty(version)
    push!(test_results, ValidationTestResult(
        "IQ-002",
        "NeoPKPD Package Installation",
        "Verify NeoPKPD Core package is properly installed",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Version string",
        version,
        passed ? nothing : "Empty version string",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED", " (", version, ")")
catch e
    push!(test_results, ValidationTestResult(
        "IQ-002", "NeoPKPD Package Installation",
        "Verify NeoPKPD Core package is properly installed",
        VALIDATION_FAILED, "Version string", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# IQ-003: Dependency Verification
println("Running IQ-003: Dependency Verification")
t1 = time()
try
    deps = Pkg.dependencies()
    required = ["DifferentialEquations", "SciMLBase", "JSON", "SHA", "UUIDs"]
    missing_deps = String[]
    for dep in required
        found = false
        for (uuid, info) in deps
            if info.name == dep
                found = true
                break
            end
        end
        if !found
            push!(missing_deps, dep)
        end
    end
    passed = isempty(missing_deps)
    push!(test_results, ValidationTestResult(
        "IQ-003",
        "Dependency Verification",
        "Verify all required dependencies are installed",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "All dependencies present",
        passed ? "All present" : "Missing: $(join(missing_deps, ", "))",
        passed ? nothing : "Missing dependencies: $(join(missing_deps, ", "))",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "IQ-003", "Dependency Verification",
        "Verify all required dependencies are installed",
        VALIDATION_FAILED, "All dependencies present", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# IQ-004: Schema Version Verification
println("Running IQ-004: Schema Version Verification")
t1 = time()
try
    schema = ARTIFACT_SCHEMA_VERSION
    expected = "1.1.0"
    passed = schema == expected
    push!(test_results, ValidationTestResult(
        "IQ-004",
        "Schema Version Verification",
        "Verify artifact schema version is current (1.1.0)",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        expected,
        schema,
        passed ? nothing : "Schema version mismatch",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED", " (", schema, ")")
catch e
    push!(test_results, ValidationTestResult(
        "IQ-004", "Schema Version Verification",
        "Verify artifact schema version is current",
        VALIDATION_FAILED, "1.1.0", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# IQ-005: Compliance Module Initialization
println("Running IQ-005: Compliance Module Initialization")
t1 = time()
try
    config = get_compliance_config()
    passed = config.level == COMPLIANCE_STANDARD
    push!(test_results, ValidationTestResult(
        "IQ-005",
        "Compliance Module Initialization",
        "Verify compliance module loads with default COMPLIANCE_STANDARD level",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "COMPLIANCE_STANDARD",
        string(config.level),
        passed ? nothing : "Default compliance level incorrect",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED", " (", config.level, ")")
catch e
    push!(test_results, ValidationTestResult(
        "IQ-005", "Compliance Module Initialization",
        "Verify compliance module loads correctly",
        VALIDATION_FAILED, "COMPLIANCE_STANDARD", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# IQ-006: Environment Capture Capability
println("Running IQ-006: Environment Capture Capability")
t1 = time()
try
    env = capture_environment()
    passed = !isempty(env.julia_version) && !isempty(env.neopkpd_version) && !isempty(env.os)
    push!(test_results, ValidationTestResult(
        "IQ-006",
        "Environment Capture Capability",
        "Verify environment snapshot can be captured with valid values",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Valid environment snapshot",
        "Julia=$(env.julia_version), NeoPKPD=$(env.neopkpd_version), OS=$(env.os)",
        passed ? nothing : "Environment snapshot incomplete",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "IQ-006", "Environment Capture Capability",
        "Verify environment snapshot can be captured",
        VALIDATION_FAILED, "Valid environment snapshot", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# IQ-007: Integrity Hash Computation
println("Running IQ-007: Integrity Hash Computation")
t1 = time()
try
    hash = compute_content_hash("test data for hashing")
    passed = length(hash) == 64 && all(c -> c in "0123456789abcdef", hash)
    push!(test_results, ValidationTestResult(
        "IQ-007",
        "Integrity Hash Computation",
        "Verify SHA-256 hashing is operational",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "64-character hexadecimal hash",
        hash,
        passed ? nothing : "Invalid hash format",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "IQ-007", "Integrity Hash Computation",
        "Verify SHA-256 hashing is operational",
        VALIDATION_FAILED, "64-character hex hash", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# IQ-008: Audit Record Creation
println("Running IQ-008: Audit Record Creation")
t1 = time()
try
    record = create_audit_record()
    passed = !isempty(record.execution_id) && !isempty(record.timestamp_utc)
    verified = verify_audit_record(record)
    passed = passed && verified
    push!(test_results, ValidationTestResult(
        "IQ-008",
        "Audit Record Creation",
        "Verify audit records can be created and verified",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Valid audit record with verified hash",
        "ID=$(record.execution_id[1:8])..., verified=$verified",
        passed ? nothing : "Audit record verification failed",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "IQ-008", "Audit Record Creation",
        "Verify audit records can be created and verified",
        VALIDATION_FAILED, "Valid audit record", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

total_time = (time() - start_time) * 1000

# Generate report
println()
println("=" ^ 60)
println("Generating Validation Report...")
println("=" ^ 60)

report = generate_validation_report(
    report_type="IQ",
    tests=test_results,
    notes=[
        "Automated IQ validation executed by run_iq.jl",
        "Part of FDA 21 CFR Part 11 compliance validation suite"
    ],
    execution_time_ms=total_time
)

# Ensure output directory exists
mkpath(dirname(output_path))

# Write report
write_validation_report_markdown(report, output_path)
println("Report written to: ", output_path)

# Print summary
println()
println("=" ^ 60)
println("IQ VALIDATION SUMMARY")
println("=" ^ 60)
println("Total Tests: ", report.total_tests)
println("Passed: ", report.passed_tests)
println("Failed: ", report.failed_tests)
println("Skipped: ", report.skipped_tests)
println("Overall Status: ", report.overall_status == VALIDATION_PASSED ? "PASSED" : "FAILED")
println("Execution Time: ", round(total_time, digits=2), " ms")
println()

# Exit with appropriate code
exit(report.overall_status == VALIDATION_PASSED ? 0 : 1)
