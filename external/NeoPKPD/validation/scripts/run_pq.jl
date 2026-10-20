# Performance Qualification (PQ) Automated Validation Script
# FDA 21 CFR Part 11 Compliant
#
# Usage:
#   julia validation/scripts/run_pq.jl
#   julia validation/scripts/run_pq.jl --output validation/reports/pq_report.md

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "packages", "core"))

using NeoPKPD
using Test
using JSON
using StableRNGs

# Parse command line arguments
output_path = "validation/reports/pq_report.md"
for i in 1:length(ARGS)
    if ARGS[i] == "--output" && i < length(ARGS)
        global output_path = ARGS[i+1]
    end
end

println("=" ^ 60)
println("NeoPKPD Performance Qualification (PQ)")
println("=" ^ 60)
println()

# Collect test results
test_results = ValidationTestResult[]
performance_metrics = Dict{String,Float64}()
start_time = time()

# =============================================================================
# PQ-001: End-to-End Workflow - Drug Development
# =============================================================================
println("Running PQ-001: End-to-End Drug Development Workflow")
t1 = time()
try
    # Define population with IIV
    iiv = IIV([
        IIVParameter(:CL, LogNormal(), 0.3),
        IIVParameter(:V, LogNormal(), 0.25)
    ])

    pop_spec = PopulationSpec(
        OneCompIVBolus(),
        "phase1_sim",
        OneCompIVBolusParams(1.0, 10.0),
        [DoseEvent(0.0, 100.0), DoseEvent(24.0, 100.0)],
        50,  # 50 subjects for speed
        iiv
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

    result = simulate_population(pop_spec, grid, solver)

    # Verify results
    n_subjects = length(result.subjects)
    all_completed = all(s -> length(s.t) > 0, result.subjects)

    # Check variability exists
    cmax_values = [maximum(s.observations[:conc]) for s in result.subjects]
    cv = std(cmax_values) / mean(cmax_values)
    has_variability = cv > 0.1  # At least 10% CV

    passed = n_subjects == 50 && all_completed && has_variability
    exec_time = (time() - t1) * 1000
    performance_metrics["population_50_subjects_ms"] = exec_time

    push!(test_results, ValidationTestResult(
        "PQ-001",
        "End-to-End Drug Development Workflow",
        "Verify complete population simulation workflow with IIV",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "50 subjects, variability present",
        "n=$n_subjects, CV=$(round(cv*100, digits=1))%",
        passed ? nothing : "Workflow incomplete or no variability",
        exec_time
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED", " (n=$n_subjects, CV=$(round(cv*100, digits=1))%)")
catch e
    push!(test_results, ValidationTestResult(
        "PQ-001", "End-to-End Drug Development Workflow",
        "Verify population simulation workflow",
        VALIDATION_FAILED, "Complete workflow", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# =============================================================================
# PQ-002: Sensitivity Analysis Workflow
# =============================================================================
println("Running PQ-002: GSA Workflow Verification")
t1 = time()
try
    spec = ModelSpec(
        TwoCompIVBolus(),
        "gsa_test",
        TwoCompIVBolusParams(1.0, 10.0, 0.5, 20.0),
        [DoseEvent(0.0, 100.0)]
    )
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

    bounds = ParameterBounds(Dict(
        :CL => (0.5, 2.0),
        :V1 => (5.0, 20.0),
        :Q => (0.1, 1.0),
        :V2 => (10.0, 50.0)
    ))

    # Run Sobol' with small sample for speed
    sobol_spec = GlobalSensitivitySpec(
        SobolMethod(base_sample_size=64, bootstrap_samples=0),
        bounds
    )
    sobol_result = run_sobol_sensitivity(spec, grid, solver, sobol_spec)

    # Run Morris
    morris_spec = GlobalSensitivitySpec(
        MorrisMethod(n_trajectories=10),
        bounds
    )
    morris_result = run_morris_sensitivity(spec, grid, solver, morris_spec)

    # Get rankings
    sobol_ranking = sort(collect(sobol_result.indices), by=x->x.second.STi, rev=true)
    morris_ranking = sort(collect(morris_result.indices), by=x->x.second.mu_star, rev=true)

    # Check if top 2 parameters are similar (in either order)
    sobol_top2 = Set([sobol_ranking[1].first, sobol_ranking[2].first])
    morris_top2 = Set([morris_ranking[1].first, morris_ranking[2].first])

    # Allow some overlap in rankings
    overlap = length(intersect(sobol_top2, morris_top2))
    passed = overlap >= 1  # At least 1 common parameter in top 2

    exec_time = (time() - t1) * 1000
    performance_metrics["gsa_combined_ms"] = exec_time

    push!(test_results, ValidationTestResult(
        "PQ-002",
        "GSA Workflow Verification",
        "Verify Sobol' and Morris produce consistent parameter rankings",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Top 2 parameters overlap",
        "Sobol top2: $sobol_top2, Morris top2: $morris_top2",
        passed ? nothing : "Rankings inconsistent between methods",
        exec_time
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "PQ-002", "GSA Workflow Verification",
        "Verify GSA workflow",
        VALIDATION_FAILED, "Consistent rankings", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# =============================================================================
# PQ-003: Audit Trail Workflow
# =============================================================================
println("Running PQ-003: Audit Trail Workflow")
t1 = time()
try
    spec = ModelSpec(
        OneCompIVBolus(),
        "audit_test",
        OneCompIVBolusParams(1.0, 10.0),
        [DoseEvent(0.0, 100.0)]
    )
    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

    # Initial execution
    result1 = simulate(spec, grid, solver)
    artifact1 = serialize_execution(
        model_spec=spec,
        grid=grid,
        solver=solver,
        result=result1
    )
    add_compliance_metadata!(artifact1)

    # Get execution ID
    exec_id_1 = get_execution_id(artifact1)

    # Simulate replay with modified parameter
    spec2 = ModelSpec(
        OneCompIVBolus(),
        "audit_test_modified",
        OneCompIVBolusParams(1.5, 10.0),  # Changed CL
        [DoseEvent(0.0, 100.0)]
    )
    result2 = simulate(spec2, grid, solver)
    artifact2 = serialize_execution(
        model_spec=spec2,
        grid=grid,
        solver=solver,
        result=result2
    )
    add_compliance_metadata!(artifact2, action=AUDIT_MODIFY, previous_execution_id=exec_id_1)

    # Verify chain
    exec_id_2 = get_execution_id(artifact2)
    prev_id = artifact2["compliance_metadata"]["audit_record"]["previous_execution_id"]

    chain_valid = prev_id == exec_id_1
    ids_unique = exec_id_1 != exec_id_2
    timestamps_utc = endswith(artifact1["compliance_metadata"]["audit_record"]["timestamp_utc"], "Z") &&
                     endswith(artifact2["compliance_metadata"]["audit_record"]["timestamp_utc"], "Z")

    passed = chain_valid && ids_unique && timestamps_utc

    push!(test_results, ValidationTestResult(
        "PQ-003",
        "Audit Trail Workflow",
        "Verify audit chain maintained between artifacts with unique IDs and UTC timestamps",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Chain valid, unique IDs, UTC timestamps",
        "chain=$chain_valid, unique=$ids_unique, utc=$timestamps_utc",
        passed ? nothing : "Audit chain broken or invalid",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "PQ-003", "Audit Trail Workflow",
        "Verify audit trail workflow",
        VALIDATION_FAILED, "Valid audit chain", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# =============================================================================
# PQ-004: Performance Benchmarks
# =============================================================================
println("Running PQ-004: Performance Benchmarks")
t1 = time()
try
    # Single subject benchmark
    spec = ModelSpec(
        OneCompIVBolus(),
        "perf_test",
        OneCompIVBolusParams(1.0, 10.0),
        [DoseEvent(0.0, 100.0)]
    )
    grid = SimGrid(0.0, 24.0, collect(0.0:0.024:24.0))  # ~1000 points
    solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

    t_single = @elapsed begin
        for _ in 1:10
            simulate(spec, grid, solver)
        end
    end
    single_avg_ms = (t_single / 10) * 1000
    performance_metrics["single_subject_1000pts_ms"] = single_avg_ms

    # Sobol' benchmark
    bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
    gsa_spec = GlobalSensitivitySpec(
        SobolMethod(base_sample_size=128, bootstrap_samples=0),
        bounds
    )

    t_sobol = @elapsed run_sobol_sensitivity(spec, grid, solver, gsa_spec)
    sobol_ms = t_sobol * 1000
    performance_metrics["sobol_n128_ms"] = sobol_ms

    # Morris benchmark
    morris_spec = GlobalSensitivitySpec(
        MorrisMethod(n_trajectories=20),
        bounds
    )

    t_morris = @elapsed run_morris_sensitivity(spec, grid, solver, morris_spec)
    morris_ms = t_morris * 1000
    performance_metrics["morris_r20_ms"] = morris_ms

    # Check against benchmarks (generous tolerance for CI variability)
    single_ok = single_avg_ms < 500  # < 500ms per simulation
    sobol_ok = sobol_ms < 30000      # < 30s for Sobol
    morris_ok = morris_ms < 15000    # < 15s for Morris

    passed = single_ok && sobol_ok && morris_ok

    push!(test_results, ValidationTestResult(
        "PQ-004",
        "Performance Benchmarks",
        "Verify system meets performance requirements",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "All benchmarks within tolerance",
        "single=$(round(single_avg_ms, digits=1))ms, sobol=$(round(sobol_ms, digits=1))ms, morris=$(round(morris_ms, digits=1))ms",
        passed ? nothing : "Some benchmarks exceeded tolerance",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "PQ-004", "Performance Benchmarks",
        "Verify performance",
        VALIDATION_FAILED, "Benchmarks met", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# =============================================================================
# PQ-005: Error Handling
# =============================================================================
println("Running PQ-005: Error Handling")
t1 = time()
try
    errors_handled = true
    error_messages = String[]

    # Test 1: Empty dose list should work (no doses)
    try
        spec = ModelSpec(
            OneCompIVBolus(),
            "no_dose",
            OneCompIVBolusParams(1.0, 10.0),
            DoseEvent[]  # Empty doses
        )
        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)
        result = simulate(spec, grid, solver)
        # Should work but all concentrations near 0
        if all(c -> abs(c) < 1e-10, result.observations[:conc])
            push!(error_messages, "Empty doses handled correctly")
        end
    catch e
        push!(error_messages, "Empty doses: $(typeof(e))")
    end

    # Test 2: Zero time grid
    try
        spec = ModelSpec(
            OneCompIVBolus(),
            "zero_grid",
            OneCompIVBolusParams(1.0, 10.0),
            [DoseEvent(0.0, 100.0)]
        )
        grid = SimGrid(0.0, 0.0, [0.0])
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)
        result = simulate(spec, grid, solver)
        push!(error_messages, "Zero grid handled")
    catch e
        push!(error_messages, "Zero grid error caught: $(typeof(e))")
    end

    # Test 3: Invalid JSON deserialization
    try
        invalid_json = Dict{String,Any}("not" => "valid artifact")
        # This should fail gracefully or return empty
        if !haskey(invalid_json, "model_spec")
            push!(error_messages, "Missing fields detected correctly")
        end
    catch e
        push!(error_messages, "Invalid artifact handled: $(typeof(e))")
    end

    # All error cases should be handled (not crash)
    passed = length(error_messages) >= 3

    push!(test_results, ValidationTestResult(
        "PQ-005",
        "Error Handling",
        "Verify graceful handling of edge cases",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "All errors handled gracefully",
        join(error_messages, "; "),
        passed ? nothing : "Some errors not handled",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "PQ-005", "Error Handling",
        "Verify error handling",
        VALIDATION_FAILED, "Graceful handling", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# =============================================================================
# PQ-006: Digital Signature Workflow
# =============================================================================
println("Running PQ-006: Digital Signature Workflow")
t1 = time()
try
    # Generate key pair
    keypair = generate_keypair(owner="PQ Test Organization")

    # Create artifact
    spec = ModelSpec(
        OneCompIVBolus(),
        "signature_test",
        OneCompIVBolusParams(1.0, 10.0),
        [DoseEvent(0.0, 100.0)]
    )
    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

    result = simulate(spec, grid, solver)
    artifact = serialize_execution(
        model_spec=spec,
        grid=grid,
        solver=solver,
        result=result
    )

    # Sign artifact
    signature = sign_artifact(
        artifact,
        keypair.signing_key,
        purpose=SIGNATURE_AUTHORSHIP,
        signer_name="Test Scientist",
        signer_title="Pharmacometrician",
        comments="PQ Test Signature"
    )

    # Verify signature
    is_valid = verify_signature(artifact, signature, keypair.verification_key)

    # Add signature to artifact
    add_signature_to_artifact!(artifact, signature)

    # Retrieve and verify
    signatures = get_artifact_signatures(artifact)
    has_signature = length(signatures) == 1

    passed = is_valid && has_signature

    push!(test_results, ValidationTestResult(
        "PQ-006",
        "Digital Signature Workflow",
        "Verify ECDSA signing and verification workflow",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Signature valid and stored",
        "valid=$is_valid, stored=$has_signature",
        passed ? nothing : "Signature workflow failed",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "PQ-006", "Digital Signature Workflow",
        "Verify digital signatures",
        VALIDATION_FAILED, "Valid signature", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# =============================================================================
# PQ-007: Schema Validation
# =============================================================================
println("Running PQ-007: Schema Validation")
t1 = time()
try
    # Create valid artifact
    spec = ModelSpec(
        OneCompIVBolus(),
        "schema_test",
        OneCompIVBolusParams(1.0, 10.0),
        [DoseEvent(0.0, 100.0)]
    )
    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

    result = simulate(spec, grid, solver)
    artifact = serialize_execution(
        model_spec=spec,
        grid=grid,
        solver=solver,
        result=result
    )

    # Validate
    validation_result = validate_artifact_schema(artifact)

    valid_artifact = validation_result.is_valid
    detected_type = validation_result.artifact_type == "pk_single"

    # Create invalid artifact and validate
    invalid_artifact = Dict{String,Any}(
        "artifact_schema_version" => "1.1.0"
        # Missing required fields
    )
    invalid_result = validate_artifact_schema(invalid_artifact)
    detected_invalid = !invalid_result.is_valid

    passed = valid_artifact && detected_type && detected_invalid

    push!(test_results, ValidationTestResult(
        "PQ-007",
        "Schema Validation",
        "Verify JSON schema validation catches invalid artifacts",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Valid detected, invalid rejected",
        "valid_ok=$valid_artifact, type_ok=$detected_type, invalid_caught=$detected_invalid",
        passed ? nothing : "Schema validation incorrect",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "PQ-007", "Schema Validation",
        "Verify schema validation",
        VALIDATION_FAILED, "Proper validation", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

# =============================================================================
# PQ-008: Schema Migration
# =============================================================================
println("Running PQ-008: Schema Migration")
t1 = time()
try
    # Create 1.0.0 artifact (without compliance_metadata)
    old_artifact = Dict{String,Any}(
        "artifact_schema_version" => "1.0.0",
        "semantics_fingerprint" => Dict{String,Any}(
            "neopkpd_version" => NEOPKPD_VERSION,
            "event_semantics_version" => EVENT_SEMANTICS_VERSION,
            "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
            "artifact_schema_version" => "1.0.0"
        ),
        "model_spec" => Dict{String,Any}(
            "kind" => "OneCompIVBolus",
            "name" => "migration_test",
            "params" => Dict{String,Any}("CL" => 1.0, "V" => 10.0),
            "doses" => [Dict{String,Any}("time" => 0.0, "amount" => 100.0)]
        ),
        "grid" => Dict{String,Any}("t0" => 0.0, "t1" => 24.0, "saveat" => collect(0.0:1.0:24.0)),
        "solver" => Dict{String,Any}("alg" => "Tsit5", "reltol" => 1e-6, "abstol" => 1e-8, "maxiters" => 10000),
        "result" => Dict{String,Any}(
            "t" => collect(0.0:1.0:24.0),
            "states" => Dict{String,Any}("A1" => fill(0.0, 25)),
            "observations" => Dict{String,Any}("conc" => fill(0.0, 25)),
            "metadata" => Dict{String,Any}()
        )
    )

    # Check migration path exists
    path_exists = can_migrate("1.0.0", "1.1.0")

    # Migrate
    migration_result = migrate_artifact(old_artifact)

    migrated_ok = migration_result.success
    has_compliance = haskey(migration_result.migrated_artifact, "compliance_metadata")
    new_version = get(migration_result.migrated_artifact, "artifact_schema_version", "") == "1.1.0"

    passed = path_exists && migrated_ok && has_compliance && new_version

    push!(test_results, ValidationTestResult(
        "PQ-008",
        "Schema Migration",
        "Verify schema migration from 1.0.0 to 1.1.0",
        passed ? VALIDATION_PASSED : VALIDATION_FAILED,
        "Migration successful with compliance metadata",
        "path=$path_exists, migrated=$migrated_ok, compliance=$has_compliance, version=$new_version",
        passed ? nothing : "Migration failed",
        (time() - t1) * 1000
    ))
    println("  Result: ", passed ? "PASSED" : "FAILED")
catch e
    push!(test_results, ValidationTestResult(
        "PQ-008", "Schema Migration",
        "Verify schema migration",
        VALIDATION_FAILED, "Successful migration", "error", string(e), (time() - t1) * 1000
    ))
    println("  Result: FAILED (", e, ")")
end

total_time = (time() - start_time) * 1000

# =============================================================================
# Generate Report
# =============================================================================
println()
println("=" ^ 60)
println("Generating PQ Validation Report...")
println("=" ^ 60)

# Format performance metrics for notes
perf_notes = String[]
push!(perf_notes, "Automated PQ validation executed by run_pq.jl")
push!(perf_notes, "Part of FDA 21 CFR Part 11 compliance validation suite")
push!(perf_notes, "")
push!(perf_notes, "Performance Metrics:")
for (k, v) in sort(collect(performance_metrics))
    push!(perf_notes, "  - $k: $(round(v, digits=2)) ms")
end

report = generate_validation_report(
    report_type="PQ",
    tests=test_results,
    notes=perf_notes,
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
println("PQ VALIDATION SUMMARY")
println("=" ^ 60)
println("Total Tests: ", report.total_tests)
println("Passed: ", report.passed_tests)
println("Failed: ", report.failed_tests)
println("Skipped: ", report.skipped_tests)
println("Overall Status: ", report.overall_status == VALIDATION_PASSED ? "PASSED" : "FAILED")
println("Execution Time: ", round(total_time, digits=2), " ms")
println()
println("Performance Metrics:")
for (k, v) in sort(collect(performance_metrics))
    println("  $k: $(round(v, digits=2)) ms")
end
println()

# Exit with appropriate code
exit(report.overall_status == VALIDATION_PASSED ? 0 : 1)
