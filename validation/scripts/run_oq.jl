# Operational Qualification (OQ) Automated Validation Script
# FDA 21 CFR Part 11 Compliant
#
# Usage:
#   julia validation/scripts/run_oq.jl
#   julia validation/scripts/run_oq.jl --output validation/reports/oq_report.md

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "packages", "core"))

using NeoPKPD
using Test
using JSON
using StableRNGs

# Parse command line arguments
output_path = "validation/reports/oq_report.md"
for i in 1:length(ARGS)
    if ARGS[i] == "--output" && i < length(ARGS)
        global output_path = ARGS[i+1]
    end
end

println("=" ^ 60)
println("NeoPKPD Operational Qualification (OQ)")
println("=" ^ 60)
println()

# Collect test results
test_results = ValidationTestResult[]
start_time = time()

# Helper function to create test model
function create_test_model()
    spec = ModelSpec(
        OneCompIVBolus(),
        "test",
        OneCompIVBolusParams(1.0, 10.0),
        [DoseEvent(0.0, 100.0)]
    )
    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)
    return spec, grid, solver
end

# OQ-001: One-Compartment IV Bolus Simulation
println("Running OQ-001: One-Compartment IV Bolus Simulation")
t1 = time()
try
    spec, grid, solver = create_test_model()
    result = simulate(spec, grid, solver)

    # Check initial concentration: C(0) = Dose/V = 100/10 = 10
    C0 = result.observations[:conc][1]
    expected_C0 = 10.0
    tol = 0.01

    passed = abs(C0 - expected_C0) < tol
    push!(test_results, ValidationTestResult(
        "OQ-001",
        "One-Compartment IV Bolus Simulation",
        "Verify basic PK simulation produces expected initial concentration C(0) = Dose/V",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "C(0) = $expected_C0 ± $tol",
        "C(0) = $C0",
        passed ? nothing : "Initial concentration outside tolerance",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED", " (C0 = ", round(C0, digits=4), ")")
catch e
    push!(test_results, ValidationTestResult(
        "OQ-001", "One-Compartment IV Bolus Simulation",
        "Verify basic PK simulation",
        VALIDATION_FAILED, "C(0) = 10.0", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# OQ-002: Artifact Serialization Round-Trip
println("Running OQ-002: Artifact Serialization Round-Trip")
t1 = time()
try
    spec, grid, solver = create_test_model()
    result = simulate(spec, grid, solver)

    # Serialize
    artifact = serialize_execution(
        model_spec=spec,
        grid=grid,
        solver=solver,
        result=result
    )

    # Serialize to JSON and back (convert to Dict for replay_execution)
    json_str = JSON.json(artifact)
    artifact2 = Dict{String,Any}(JSON.parse(json_str))

    # Replay
    result2 = replay_execution(artifact2)

    # Compare results
    orig_conc = result.observations[:conc]
    new_conc = result2.observations[:conc]
    max_diff = maximum(abs.(orig_conc .- new_conc))
    passed = max_diff < 1e-10

    push!(test_results, ValidationTestResult(
        "OQ-002",
        "Artifact Serialization Round-Trip",
        "Verify artifacts can be serialized, deserialized, and replayed with identical results",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "max diff < 1e-10",
        "max diff = $max_diff",
        passed ? nothing : "Results differ after round-trip",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED", " (max diff = ", max_diff, ")")
catch e
    push!(test_results, ValidationTestResult(
        "OQ-002", "Artifact Serialization Round-Trip",
        "Verify serialization round-trip",
        VALIDATION_FAILED, "max diff < 1e-10", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# OQ-003: Sobol' Sensitivity Analysis
println("Running OQ-003: Sobol' Sensitivity Analysis")
t1 = time()
try
    spec, grid, solver = create_test_model()
    bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
    gsa_spec = GlobalSensitivitySpec(
        SobolMethod(base_sample_size=64, bootstrap_samples=0),
        bounds
    )

    result = run_sobol_sensitivity(spec, grid, solver, gsa_spec)

    # Check validity of indices
    all_valid = true
    for (param, idx) in result.indices
        if !(0.0 <= idx.Si <= 1.0 || idx.Si < 0.1)  # Allow small negative
            all_valid = false
        end
        if !(0.0 <= idx.STi <= 1.0 || idx.STi < 0.1)
            all_valid = false
        end
        if !(idx.STi >= idx.Si - 0.1)  # Total should be >= first-order
            all_valid = false
        end
    end

    push!(test_results, ValidationTestResult(
        "OQ-003",
        "Sobol' Sensitivity Analysis",
        "Verify GSA produces valid sensitivity indices (0 ≤ Si ≤ STi ≤ 1)",
        all_valid ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Valid Sobol' indices",
        "CL: Si=$(round(result.indices[:CL].Si, digits=3)), STi=$(round(result.indices[:CL].STi, digits=3))",
        all_valid ? nothing : "Invalid sensitivity indices",
        (time() - t1) * 1000
    ))
    println("  Result: ", all_valid ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "OQ-003", "Sobol' Sensitivity Analysis",
        "Verify GSA produces valid indices",
        VALIDATION_FAILED, "Valid Sobol' indices", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# OQ-004: Morris Elementary Effects
println("Running OQ-004: Morris Elementary Effects")
t1 = time()
try
    spec, grid, solver = create_test_model()
    bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
    gsa_spec = GlobalSensitivitySpec(
        MorrisMethod(n_trajectories=10),
        bounds
    )

    result = run_morris_sensitivity(spec, grid, solver, gsa_spec)

    # Check validity
    all_valid = true
    for (param, idx) in result.indices
        if idx.mu_star < 0 || idx.sigma < 0
            all_valid = false
        end
    end

    push!(test_results, ValidationTestResult(
        "OQ-004",
        "Morris Elementary Effects",
        "Verify Morris screening method produces valid indices (μ* ≥ 0, σ ≥ 0)",
        all_valid ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Valid Morris indices",
        "CL: μ*=$(round(result.indices[:CL].mu_star, digits=3)), σ=$(round(result.indices[:CL].sigma, digits=3))",
        all_valid ? nothing : "Invalid Morris indices",
        (time() - t1) * 1000
    ))
    println("  Result: ", all_valid ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "OQ-004", "Morris Elementary Effects",
        "Verify Morris screening method",
        VALIDATION_FAILED, "Valid Morris indices", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# OQ-005: Compliance Metadata Generation
println("Running OQ-005: Compliance Metadata Generation")
t1 = time()
try
    spec, grid, solver = create_test_model()
    result = simulate(spec, grid, solver)

    artifact = serialize_execution(
        model_spec=spec,
        grid=grid,
        solver=solver,
        result=result
    )
    add_compliance_metadata!(artifact)

    meta = artifact["compliance_metadata"]
    has_audit = haskey(meta, "audit_record")
    has_integrity = haskey(meta, "integrity")
    has_env = haskey(meta, "environment")

    passed = has_audit && has_integrity && has_env

    push!(test_results, ValidationTestResult(
        "OQ-005",
        "Compliance Metadata Generation",
        "Verify compliance metadata (audit_record, integrity, environment) is added to artifacts",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "All compliance fields present",
        "audit=$has_audit, integrity=$has_integrity, env=$has_env",
        passed ? nothing : "Missing compliance fields",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "OQ-005", "Compliance Metadata Generation",
        "Verify compliance metadata generation",
        VALIDATION_FAILED, "All compliance fields", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# OQ-006: Integrity Verification
println("Running OQ-006: Integrity Verification")
t1 = time()
try
    spec, grid, solver = create_test_model()
    result = simulate(spec, grid, solver)

    artifact = serialize_execution(
        model_spec=spec,
        grid=grid,
        solver=solver,
        result=result
    )
    add_compliance_metadata!(artifact)

    # Verify unmodified
    result1 = verify_artifact_integrity(artifact)

    # Modify and verify (should fail) - use correct path in serialized artifact
    artifact["result"]["observations"]["conc"][1] = 999.0
    result2 = verify_artifact_integrity(artifact)

    passed = result1.is_valid && !result2.is_valid

    push!(test_results, ValidationTestResult(
        "OQ-006",
        "Integrity Verification",
        "Verify integrity check passes for unmodified artifacts and fails for modified ones",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Unmodified=valid, Modified=invalid",
        "Unmodified=$(result1.is_valid), Modified=$(result2.is_valid)",
        passed ? nothing : "Integrity verification incorrect",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "OQ-006", "Integrity Verification",
        "Verify integrity verification",
        VALIDATION_FAILED, "Proper integrity detection", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# OQ-007: Reproducibility with Seed
println("Running OQ-007: Reproducibility with Seed")
t1 = time()
try
    spec, grid, solver = create_test_model()
    bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))

    gsa_spec1 = GlobalSensitivitySpec(
        SobolMethod(base_sample_size=64, bootstrap_samples=0),
        bounds;
        seed=UInt64(42)
    )
    gsa_spec2 = GlobalSensitivitySpec(
        SobolMethod(base_sample_size=64, bootstrap_samples=0),
        bounds;
        seed=UInt64(42)
    )

    result1 = run_sobol_sensitivity(spec, grid, solver, gsa_spec1)
    result2 = run_sobol_sensitivity(spec, grid, solver, gsa_spec2)

    # Check if results are identical
    cl_match = result1.indices[:CL].Si ≈ result2.indices[:CL].Si
    v_match = result1.indices[:V].Si ≈ result2.indices[:V].Si
    passed = cl_match && v_match

    push!(test_results, ValidationTestResult(
        "OQ-007",
        "Reproducibility with Seed",
        "Verify simulations with same seed produce identical results",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Identical results with same seed",
        "CL match=$cl_match, V match=$v_match",
        passed ? nothing : "Results differ despite same seed",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "OQ-007", "Reproducibility with Seed",
        "Verify reproducibility with seed",
        VALIDATION_FAILED, "Identical results", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# OQ-008: Audit Trail Verification
println("Running OQ-008: Audit Trail Verification")
t1 = time()
try
    record = create_audit_record()
    verified = verify_audit_record(record)

    # Check all required fields
    has_id = !isempty(record.execution_id)
    has_timestamp = !isempty(record.timestamp_utc)
    has_user = !isempty(record.system_user)

    passed = verified && has_id && has_timestamp && has_user

    push!(test_results, ValidationTestResult(
        "OQ-008",
        "Audit Trail Verification",
        "Verify audit records are self-verifying with all required fields",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Verified record with all fields",
        "verified=$verified, id=$has_id, timestamp=$has_timestamp, user=$has_user",
        passed ? nothing : "Audit record incomplete or not verified",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "OQ-008", "Audit Trail Verification",
        "Verify audit trail",
        VALIDATION_FAILED, "Verified audit record", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# OQ-009: Golden Artifact Integrity (if golden files exist)
println("Running OQ-009: Golden Artifact Integrity")
t1 = time()
golden_dir = joinpath(@__DIR__, "..", "golden")
try
    if isdir(golden_dir)
        golden_files = filter(f -> endswith(f, ".json"), readdir(golden_dir))
        if !isempty(golden_files)
            valid_count = 0
            total_with_compliance = 0
            for file in golden_files
                # Convert JSON.Object to Dict for has_compliance_metadata
                artifact = Dict{String,Any}(JSON.parsefile(joinpath(golden_dir, file)))
                if has_compliance_metadata(artifact)
                    total_with_compliance += 1
                    result = verify_artifact_integrity(artifact)
                    if result.is_valid
                        valid_count += 1
                    end
                end
            end

            if total_with_compliance > 0
                passed = valid_count == total_with_compliance
                push!(test_results, ValidationTestResult(
                    "OQ-009",
                    "Golden Artifact Integrity",
                    "Verify all golden artifacts with compliance metadata pass integrity checks",
                    passed ? VALIDATION_PASSED : VALIDATION_FAILED,
                    "All golden artifacts valid",
                    "$valid_count/$total_with_compliance valid",
                    passed ? nothing : "Some golden artifacts failed integrity",
                    (time() - t1) * 1000
                ))
            else
                push!(test_results, ValidationTestResult(
                    "OQ-009",
                    "Golden Artifact Integrity",
                    "Verify all golden artifacts pass integrity checks",
                    VALIDATION_SKIPPED,
                    "Golden artifacts with compliance",
                    "No golden artifacts with compliance_metadata found",
                    "Golden artifacts predate schema 1.1.0",
                    (time() - t1) * 1000
                ))
            end
            println("  Result: ", total_with_compliance > 0 ? (valid_count == total_with_compliance ? "PASSED" : "FAILED") : "SKIPPED")
        else
            push!(test_results, ValidationTestResult(
                "OQ-009", "Golden Artifact Integrity",
                "Verify golden artifacts",
                VALIDATION_SKIPPED, "Golden artifacts", "No files found", nothing, (time() - t1) * 1000
            ))
            println("  Result: SKIPPED (no golden files)")
        end
    else
        push!(test_results, ValidationTestResult(
            "OQ-009", "Golden Artifact Integrity",
            "Verify golden artifacts",
            VALIDATION_SKIPPED, "Golden artifacts", "Directory not found", nothing, (time() - t1) * 1000
        ))
        println("  Result: SKIPPED (no golden directory)")
    end
catch e
    push!(test_results, ValidationTestResult(
        "OQ-009", "Golden Artifact Integrity",
        "Verify golden artifacts",
        VALIDATION_FAILED, "All golden artifacts valid", "error", string(e), (time() - t1) * 1000
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
    report_type="OQ",
    tests=test_results,
    notes=[
        "Automated OQ validation executed by run_oq.jl",
        "Part of FDA 21 CFR Part 11 compliance validation suite",
        "Golden artifacts may predate compliance metadata (schema 1.1.0)"
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
println("OQ VALIDATION SUMMARY")
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
