using Test
using NeoPKPD
using JSON

@testset "FDA 21 CFR Part 11 Compliance Module" begin

    @testset "ComplianceConfig" begin
        @testset "Default configuration" begin
            config = get_compliance_config()
            @test config.level == COMPLIANCE_STANDARD
            @test config.enable_audit_trail == true
            @test config.enable_hashing == true
            @test config.enable_environment_capture == true
            @test config.enable_validation_reports == false
        end

        @testset "Configuration from level" begin
            disabled = ComplianceConfig(COMPLIANCE_DISABLED)
            @test disabled.enable_audit_trail == false
            @test disabled.enable_hashing == false

            minimal = ComplianceConfig(COMPLIANCE_MINIMAL)
            @test minimal.enable_environment_capture == true
            @test minimal.enable_hashing == false

            standard = ComplianceConfig(COMPLIANCE_STANDARD)
            @test standard.enable_audit_trail == true
            @test standard.enable_hashing == true

            strict = ComplianceConfig(COMPLIANCE_STRICT)
            @test strict.enable_validation_reports == true
        end

        @testset "Set and get configuration" begin
            old_config = get_compliance_config()

            set_compliance_config!(COMPLIANCE_DISABLED)
            @test get_compliance_config().level == COMPLIANCE_DISABLED

            set_compliance_config!(COMPLIANCE_STANDARD)
            @test get_compliance_config().level == COMPLIANCE_STANDARD

            # Restore original
            set_compliance_config!(old_config)
        end

        @testset "with_compliance context" begin
            @test get_compliance_config().level == COMPLIANCE_STANDARD

            with_compliance(COMPLIANCE_DISABLED) do
                @test get_compliance_config().level == COMPLIANCE_DISABLED
            end

            # Should be restored
            @test get_compliance_config().level == COMPLIANCE_STANDARD
        end

        @testset "Helper functions" begin
            set_compliance_config!(COMPLIANCE_STANDARD)
            @test is_audit_trail_enabled() == true
            @test is_hashing_enabled() == true
            @test is_environment_capture_enabled() == true
            @test is_validation_reports_enabled() == false

            with_compliance(COMPLIANCE_DISABLED) do
                @test is_audit_trail_enabled() == false
                @test is_hashing_enabled() == false
            end
        end
    end

    @testset "EnvironmentSnapshot" begin
        @testset "Capture environment" begin
            env = capture_environment()

            @test !isempty(env.julia_version)
            @test !isempty(env.neopkpd_version)
            @test env.neopkpd_version == NEOPKPD_VERSION
            @test env.artifact_schema_version == ARTIFACT_SCHEMA_VERSION
            @test !isempty(env.os)
            @test env.cpu_cores > 0
            @test env.memory_gb > 0
            @test !isempty(env.architecture)
            @test env.word_size in [32, 64]
            @test !isempty(env.capture_timestamp)
        end

        @testset "Serialize environment" begin
            env = capture_environment()
            serialized = serialize_environment(env)

            @test serialized["julia_version"] == env.julia_version
            @test serialized["neopkpd_version"] == env.neopkpd_version
            @test serialized["cpu_cores"] == env.cpu_cores
        end

        @testset "Deserialize environment" begin
            env = capture_environment()
            serialized = serialize_environment(env)
            deserialized = deserialize_environment(serialized)

            @test deserialized.julia_version == env.julia_version
            @test deserialized.neopkpd_version == env.neopkpd_version
            @test deserialized.cpu_cores == env.cpu_cores
        end
    end

    @testset "IntegrityMetadata" begin
        @testset "Compute content hash" begin
            hash1 = compute_content_hash("test data")
            hash2 = compute_content_hash("test data")
            hash3 = compute_content_hash("different data")

            @test length(hash1) == 64
            @test all(c -> c in "0123456789abcdef", hash1)
            @test hash1 == hash2
            @test hash1 != hash3
        end

        @testset "Compute dict hash" begin
            d = Dict("key" => "value", "number" => 42)
            hash = compute_content_hash(d)

            @test length(hash) == 64
        end

        @testset "Canonical JSON determinism" begin
            d1 = Dict("b" => 2, "a" => 1)
            d2 = Dict("a" => 1, "b" => 2)

            hash1 = compute_content_hash(d1)
            hash2 = compute_content_hash(d2)

            @test hash1 == hash2  # Should be same despite insertion order
        end

        @testset "Compute integrity metadata" begin
            artifact = Dict{String,Any}(
                "model_spec" => Dict("CL" => 1.0, "V" => 10.0),
                "result" => Dict("concentration" => [10.0, 5.0, 2.5]),
                "semantics_fingerprint" => Dict("version" => "1.0.0")
            )

            meta = compute_integrity_metadata(artifact)

            @test length(meta.content_hash) == 64
            @test meta.hash_algorithm == "SHA-256"
            @test length(meta.input_hash) == 64
            @test length(meta.output_hash) == 64
            @test length(meta.semantics_hash) == 64
            @test !isempty(meta.hash_computed_at)
        end

        @testset "Serialize integrity" begin
            artifact = Dict{String,Any}("key" => "value")
            meta = compute_integrity_metadata(artifact)
            serialized = serialize_integrity(meta)

            @test serialized["content_hash"] == meta.content_hash
            @test serialized["hash_algorithm"] == "SHA-256"
        end
    end

    @testset "AuditRecord" begin
        @testset "Create audit record" begin
            record = create_audit_record()

            @test !isempty(record.execution_id)
            @test length(split(record.execution_id, "-")) == 5  # UUID format
            @test !isempty(record.timestamp_utc)
            @test endswith(record.timestamp_utc, "Z")  # UTC
            @test !isempty(record.system_user)
            @test !isempty(record.hostname)
            @test record.action == AUDIT_CREATE
            @test record.previous_execution_id === nothing
            @test record.neopkpd_version == NEOPKPD_VERSION
            @test length(record.record_hash) == 64
        end

        @testset "Create replay audit record" begin
            original = create_audit_record()
            replay = create_audit_record(
                action=AUDIT_REPLAY,
                previous_execution_id=original.execution_id
            )

            @test replay.action == AUDIT_REPLAY
            @test replay.previous_execution_id == original.execution_id
            @test replay.execution_id != original.execution_id
        end

        @testset "Verify audit record" begin
            record = create_audit_record()
            @test verify_audit_record(record) == true

            # Tampered record should fail
            tampered = AuditRecord(
                record.execution_id,
                record.timestamp_utc,
                record.timezone_offset,
                "fake_user",  # Tampered
                record.hostname,
                record.execution_context,
                record.action,
                record.previous_execution_id,
                record.neopkpd_version,
                record.record_hash  # Original hash won't match
            )
            @test verify_audit_record(tampered) == false
        end

        @testset "Serialize audit record" begin
            record = create_audit_record()
            serialized = serialize_audit_record(record)

            @test serialized["execution_id"] == record.execution_id
            @test serialized["action"] == "create"
            @test serialized["neopkpd_version"] == NEOPKPD_VERSION
        end

        @testset "Deserialize audit record" begin
            record = create_audit_record()
            serialized = serialize_audit_record(record)
            deserialized = deserialize_audit_record(serialized)

            @test deserialized.execution_id == record.execution_id
            @test deserialized.action == record.action
        end
    end

    @testset "Compliance Metadata Integration" begin
        @testset "Add compliance metadata to artifact" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
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

            add_compliance_metadata!(artifact)

            @test haskey(artifact, "compliance_metadata")
            meta = artifact["compliance_metadata"]
            @test haskey(meta, "audit_record")
            @test haskey(meta, "integrity")
            @test haskey(meta, "environment")
        end

        @testset "Verify artifact integrity" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
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
            add_compliance_metadata!(artifact)

            # Should be valid
            verification = verify_artifact_integrity(artifact)
            @test verification.is_valid == true
            @test verification.content_hash_valid == true

            # Tamper and verify (should fail)
            artifact["result"]["observations"]["conc"][1] = 999.0
            verification2 = verify_artifact_integrity(artifact)
            @test verification2.is_valid == false
        end

        @testset "Has compliance metadata" begin
            artifact_with = Dict{String,Any}("compliance_metadata" => Dict())
            artifact_without = Dict{String,Any}("data" => 123)

            @test has_compliance_metadata(artifact_with) == true
            @test has_compliance_metadata(artifact_without) == false
        end

        @testset "Get execution ID" begin
            artifact = Dict{String,Any}(
                "compliance_metadata" => Dict{String,Any}(
                    "audit_record" => Dict{String,Any}(
                        "execution_id" => "test-id-123"
                    )
                )
            )

            @test get_execution_id(artifact) == "test-id-123"
            @test get_execution_id(Dict{String,Any}()) === nothing
        end

        @testset "Compliance disabled mode" begin
            artifact = Dict{String,Any}("data" => 123)

            with_compliance(COMPLIANCE_DISABLED) do
                add_compliance_metadata!(artifact)
            end

            @test !haskey(artifact, "compliance_metadata")
        end
    end

    @testset "Validation Reports" begin
        @testset "Generate validation report" begin
            tests = [
                ValidationTestResult(
                    "TEST-001", "Test 1", "Description 1",
                    VALIDATION_PASSED, "expected", "actual", nothing, 100.0
                ),
                ValidationTestResult(
                    "TEST-002", "Test 2", "Description 2",
                    VALIDATION_FAILED, "expected", "wrong", "Error message", 50.0
                ),
            ]

            report = generate_validation_report(
                report_type="IQ",
                tests=tests,
                notes=["Test note"],
                execution_time_ms=150.0
            )

            @test report.report_type == "IQ"
            @test report.total_tests == 2
            @test report.passed_tests == 1
            @test report.failed_tests == 1
            @test report.overall_status == VALIDATION_FAILED
        end

        @testset "Write validation report" begin
            tests = [
                ValidationTestResult(
                    "TEST-001", "Test 1", "Description",
                    VALIDATION_PASSED, "expected", "actual", nothing, 100.0
                ),
            ]

            report = generate_validation_report(
                report_type="IQ",
                tests=tests
            )

            temp_path = tempname() * ".md"
            write_validation_report_markdown(report, temp_path)

            @test isfile(temp_path)
            content = read(temp_path, String)
            @test occursin("IQ Validation Report", content)
            @test occursin("TEST-001", content)

            rm(temp_path)
        end
    end

    @testset "Schema Version" begin
        @test ARTIFACT_SCHEMA_VERSION == "1.1.0"
    end

    # =========================================================================
    # Phase 1 Gold-Tier Tests
    # =========================================================================

    @testset "Digital Signatures (Phase 1 Gold-Tier)" begin
        @testset "Key pair generation - ECDSA" begin
            keypair = generate_keypair(algorithm=ECDSA_P256, owner="Test Org")

            @test !isempty(keypair.signing_key.key_id)
            @test keypair.signing_key.key_id == keypair.verification_key.key_id
            @test keypair.signing_key.algorithm == ECDSA_P256
            @test keypair.verification_key.algorithm == ECDSA_P256
            @test keypair.signing_key.owner == "Test Org"
            @test length(keypair.signing_key.key_data) == 32  # ECDSA P-256 private key
            @test length(keypair.verification_key.key_data) == 64  # ECDSA P-256 public key
        end

        @testset "Key pair generation - RSA" begin
            keypair = generate_keypair(algorithm=RSA_2048, owner="RSA Test")

            @test keypair.signing_key.algorithm == RSA_2048
            @test keypair.verification_key.algorithm == RSA_2048
            @test length(keypair.signing_key.key_data) == 256  # RSA 2048-bit key
            @test length(keypair.verification_key.key_data) == 256
        end

        @testset "Sign artifact" begin
            keypair = generate_keypair(owner="Signer Corp")

            artifact = Dict{String,Any}(
                "model_spec" => Dict("CL" => 1.0, "V" => 10.0),
                "result" => Dict("concentration" => [10.0, 5.0])
            )

            signature = sign_artifact(
                artifact,
                keypair.signing_key,
                purpose=SIGNATURE_AUTHORSHIP,
                signer_name="Dr. Test",
                signer_title="Scientist",
                comments="Test signature"
            )

            @test !isempty(signature.signature_id)
            @test signature.key_id == keypair.signing_key.key_id
            @test signature.algorithm == ECDSA_P256
            @test length(signature.signature_bytes) == 64  # ECDSA signature
            @test length(signature.signed_hash) == 64  # SHA-256 hash
            @test signature.purpose == SIGNATURE_AUTHORSHIP
            @test signature.signer_name == "Dr. Test"
            @test signature.signer_title == "Scientist"
            @test !isempty(signature.signed_at)
        end

        @testset "Verify signature - valid" begin
            keypair = generate_keypair(owner="Verifier Corp")

            artifact = Dict{String,Any}(
                "data" => "test content",
                "value" => 42
            )

            signature = sign_artifact(artifact, keypair.signing_key)
            is_valid = verify_signature(artifact, signature, keypair.verification_key)

            @test is_valid == true
        end

        @testset "Verify signature - tampered content" begin
            keypair = generate_keypair(owner="Tamper Test")

            artifact = Dict{String,Any}(
                "data" => "original content"
            )

            signature = sign_artifact(artifact, keypair.signing_key)

            # Tamper with artifact
            artifact["data"] = "modified content"

            is_valid = verify_signature(artifact, signature, keypair.verification_key)
            @test is_valid == false
        end

        @testset "Verify signature - wrong key" begin
            keypair1 = generate_keypair(owner="Owner 1")
            keypair2 = generate_keypair(owner="Owner 2")

            artifact = Dict{String,Any}("data" => "test")

            signature = sign_artifact(artifact, keypair1.signing_key)

            # Try to verify with wrong key
            is_valid = verify_signature(artifact, signature, keypair2.verification_key)
            @test is_valid == false
        end

        @testset "Multiple signature purposes" begin
            keypair = generate_keypair(owner="Multi-Purpose Test")
            artifact = Dict{String,Any}("data" => "test")

            for purpose in [SIGNATURE_AUTHORSHIP, SIGNATURE_APPROVAL, SIGNATURE_REVIEW,
                            SIGNATURE_VERIFICATION, SIGNATURE_RESPONSIBILITY]
                sig = sign_artifact(artifact, keypair.signing_key, purpose=purpose)
                @test sig.purpose == purpose
            end
        end

        @testset "Serialize and deserialize signature" begin
            keypair = generate_keypair(owner="Serialize Test")
            artifact = Dict{String,Any}("data" => "test")

            signature = sign_artifact(
                artifact,
                keypair.signing_key,
                signer_name="Test Signer",
                signer_title="Engineer",
                signer_organization="Tech Corp"
            )

            serialized = serialize_signature(signature)

            @test serialized["signature_id"] == signature.signature_id
            @test serialized["algorithm"] == "ECDSA-P256"
            @test serialized["purpose"] == "authorship"
            @test serialized["signer_name"] == "Test Signer"

            deserialized = deserialize_signature(serialized)

            @test deserialized.signature_id == signature.signature_id
            @test deserialized.algorithm == signature.algorithm
            @test deserialized.purpose == signature.purpose
        end

        @testset "Add signature to artifact" begin
            keypair = generate_keypair(owner="Artifact Test")
            artifact = Dict{String,Any}("data" => "test")

            signature = sign_artifact(artifact, keypair.signing_key)
            add_signature_to_artifact!(artifact, signature)

            @test haskey(artifact, "compliance_metadata")
            @test haskey(artifact["compliance_metadata"], "signatures")
            @test length(artifact["compliance_metadata"]["signatures"]["signatures"]) == 1
        end

        @testset "Get artifact signatures" begin
            keypair = generate_keypair(owner="Get Signatures Test")
            artifact = Dict{String,Any}("data" => "test")

            sig1 = sign_artifact(artifact, keypair.signing_key, purpose=SIGNATURE_AUTHORSHIP)
            sig2 = sign_artifact(artifact, keypair.signing_key, purpose=SIGNATURE_APPROVAL)

            add_signature_to_artifact!(artifact, sig1)
            add_signature_to_artifact!(artifact, sig2)

            signatures = get_artifact_signatures(artifact)
            @test length(signatures) == 2
        end

        @testset "Verify all signatures" begin
            keypair = generate_keypair(owner="Verify All Test")
            artifact = Dict{String,Any}("data" => "test")

            sig1 = sign_artifact(artifact, keypair.signing_key, purpose=SIGNATURE_AUTHORSHIP)
            sig2 = sign_artifact(artifact, keypair.signing_key, purpose=SIGNATURE_REVIEW)

            add_signature_to_artifact!(artifact, sig1)
            add_signature_to_artifact!(artifact, sig2)

            results = verify_all_signatures(artifact, [keypair.verification_key])

            @test length(results) == 2
            @test all(r -> r[2] == true, results)  # All valid
        end

        @testset "Export and import public key" begin
            keypair = generate_keypair(owner="Export Test")

            exported = export_public_key(keypair.verification_key)

            @test exported["key_type"] == "verification"
            @test exported["algorithm"] == "ECDSA-P256"
            @test !isempty(exported["key_data"])

            imported = import_public_key(exported)

            @test imported.key_id == keypair.verification_key.key_id
            @test imported.algorithm == keypair.verification_key.algorithm
            @test imported.owner == keypair.verification_key.owner
        end
    end

    @testset "Schema Validation (Phase 1 Gold-Tier)" begin
        @testset "Validate pk_single artifact" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "validation_test",
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

            validation = validate_artifact_schema(artifact)

            @test validation.is_valid == true
            @test validation.artifact_type == "pk_single"
            @test isempty(validation.errors)
        end

        @testset "Detect missing required fields" begin
            artifact = Dict{String,Any}(
                "artifact_schema_version" => "1.1.0"
                # Missing: semantics_fingerprint, model_spec, grid, solver, result
            )

            validation = validate_artifact_schema(artifact)

            @test validation.is_valid == false
            @test !isempty(validation.errors)

            # Check for specific missing field errors
            missing_fields = [e.path for e in validation.errors if e.error_type == "missing_field"]
            @test "semantics_fingerprint" in missing_fields
        end

        @testset "Detect invalid field types" begin
            artifact = Dict{String,Any}(
                "artifact_schema_version" => 123,  # Should be String
                "semantics_fingerprint" => "not a dict"  # Should be Dict
            )

            validation = validate_artifact_schema(artifact)

            @test validation.is_valid == false

            type_errors = [e for e in validation.errors if e.error_type == "invalid_type"]
            @test !isempty(type_errors)
        end

        @testset "Detect artifact type" begin
            # pk_single
            pk_artifact = Dict{String,Any}(
                "model_spec" => Dict(),
                "result" => Dict()
            )
            @test validate_artifact_schema(pk_artifact).artifact_type == "pk_single"

            # population
            pop_artifact = Dict{String,Any}(
                "population_spec" => Dict()
            )
            @test validate_artifact_schema(pop_artifact).artifact_type == "population"
        end

        @testset "Get schema definition" begin
            pk_schema = get_schema_definition("pk_single")
            @test pk_schema !== nothing
            @test haskey(pk_schema, "required_fields")
            @test "model_spec" in pk_schema["required_fields"]

            unknown_schema = get_schema_definition("unknown_type")
            @test unknown_schema === nothing
        end

        @testset "Schema registry completeness" begin
            @test haskey(SCHEMA_REGISTRY, "common")
            @test haskey(SCHEMA_REGISTRY, "pk_single")
            @test haskey(SCHEMA_REGISTRY, "population")
            @test haskey(SCHEMA_REGISTRY, "sobol_sensitivity")
            @test haskey(SCHEMA_REGISTRY, "morris_sensitivity")
            @test haskey(SCHEMA_REGISTRY, "compliance_metadata")
        end

        @testset "Validate compliance metadata structure" begin
            artifact = Dict{String,Any}(
                "artifact_schema_version" => "1.1.0",
                "semantics_fingerprint" => Dict{String,Any}(),
                "compliance_metadata" => Dict{String,Any}(
                    "audit_record" => Dict{String,Any}(
                        "execution_id" => "test-id",
                        "timestamp_utc" => "2026-01-10T00:00:00Z",
                        "system_user" => "test",
                        "hostname" => "test-host",
                        "action" => "create",
                        "neopkpd_version" => "0.1.0",
                        "record_hash" => "abc123"
                    ),
                    "integrity" => Dict{String,Any}(
                        "content_hash" => "abc",
                        "hash_algorithm" => "SHA-256",
                        "input_hash" => "def",
                        "output_hash" => "ghi",
                        "hash_computed_at" => "2026-01-10T00:00:00Z"
                    ),
                    "environment" => Dict{String,Any}(
                        "julia_version" => "1.10.0",
                        "neopkpd_version" => "0.1.0",
                        "os" => "Darwin",
                        "capture_timestamp" => "2026-01-10T00:00:00Z"
                    )
                )
            )

            validation = validate_artifact_schema(artifact)
            # Should have some warnings about unknown type but no critical errors
            @test validation.artifact_type == "unknown"
        end

        @testset "Serialize validation result" begin
            validation = SchemaValidationResult(
                true,
                "1.1.0",
                "pk_single",
                SchemaValidationError[],
                String[],
                "2026-01-10T00:00:00Z"
            )

            serialized = serialize_validation_result(validation)

            @test serialized["is_valid"] == true
            @test serialized["schema_version"] == "1.1.0"
            @test serialized["artifact_type"] == "pk_single"
        end
    end

    @testset "Schema Migration (Phase 1 Gold-Tier)" begin
        @testset "Migration registry populated" begin
            @test haskey(MIGRATION_REGISTRY, ("1.0.0", "1.1.0"))
            @test haskey(MIGRATION_FUNCTIONS, ("1.0.0", "1.1.0"))
        end

        @testset "Migration step metadata" begin
            step = MIGRATION_REGISTRY[("1.0.0", "1.1.0")]

            @test step.from_version == "1.0.0"
            @test step.to_version == "1.1.0"
            @test !isempty(step.description)
            @test isempty(step.breaking_changes)
            @test "compliance_metadata" in step.added_fields
        end

        @testset "Can migrate check" begin
            @test can_migrate("1.0.0", "1.1.0") == true
            @test can_migrate("1.1.0", "1.1.0") == true  # Same version
            @test can_migrate("0.9.0", "1.1.0") == false  # No path
            @test can_migrate("1.1.0", "1.0.0") == false  # No downgrade path
        end

        @testset "Get migration path" begin
            path = get_migration_path("1.0.0", "1.1.0")

            @test length(path) == 1
            @test path[1].from_version == "1.0.0"
            @test path[1].to_version == "1.1.0"

            # Same version - empty path
            same_path = get_migration_path("1.1.0", "1.1.0")
            @test isempty(same_path)

            # No path exists
            no_path = get_migration_path("0.5.0", "1.1.0")
            @test isempty(no_path)
        end

        @testset "Migrate artifact 1.0.0 to 1.1.0" begin
            old_artifact = Dict{String,Any}(
                "artifact_schema_version" => "1.0.0",
                "semantics_fingerprint" => Dict{String,Any}(
                    "neopkpd_version" => NEOPKPD_VERSION,
                    "artifact_schema_version" => "1.0.0"
                ),
                "model_spec" => Dict{String,Any}(
                    "kind" => "OneCompIVBolus",
                    "name" => "test",
                    "params" => Dict("CL" => 1.0, "V" => 10.0),
                    "doses" => []
                ),
                "grid" => Dict{String,Any}("t0" => 0.0, "t1" => 24.0, "saveat" => [0.0, 1.0]),
                "solver" => Dict{String,Any}("alg" => "Tsit5", "reltol" => 1e-6, "abstol" => 1e-8, "maxiters" => 10000),
                "result" => Dict{String,Any}("t" => [0.0, 1.0], "states" => Dict(), "observations" => Dict(), "metadata" => Dict())
            )

            result = migrate_artifact(old_artifact)

            @test result.success == true
            @test result.source_version == "1.0.0"
            @test result.target_version == "1.1.0"
            @test length(result.steps_applied) == 1
            @test isempty(result.errors)
            @test !isempty(result.migration_id)

            migrated = result.migrated_artifact
            @test migrated["artifact_schema_version"] == "1.1.0"
            @test haskey(migrated, "compliance_metadata")
            @test haskey(migrated["compliance_metadata"], "audit_record")
            @test haskey(migrated["compliance_metadata"], "integrity")
            @test haskey(migrated["compliance_metadata"], "environment")
        end

        @testset "Migration already at target" begin
            artifact = Dict{String,Any}(
                "artifact_schema_version" => "1.1.0"
            )

            result = migrate_artifact(artifact)

            @test result.success == true
            @test isempty(result.steps_applied)
            @test occursin("already at target", lowercase(result.warnings[1]))
        end

        @testset "Migration no path exists" begin
            artifact = Dict{String,Any}(
                "artifact_schema_version" => "0.5.0"  # Unknown version
            )

            result = migrate_artifact(artifact)

            @test result.success == false
            @test !isempty(result.errors)
            @test any(e -> occursin("No migration path", e), result.errors)
        end

        @testset "Serialize migration result" begin
            result = MigrationResult(
                true,
                "1.0.0",
                "1.1.0",
                Dict{String,Any}(),
                [MIGRATION_REGISTRY[("1.0.0", "1.1.0")]],
                String[],
                ["Test warning"],
                "test-id",
                "2026-01-10T00:00:00Z"
            )

            serialized = serialize_migration_result(result)

            @test serialized["success"] == true
            @test serialized["source_version"] == "1.0.0"
            @test serialized["target_version"] == "1.1.0"
            @test length(serialized["steps_applied"]) == 1
        end

        @testset "Get supported versions" begin
            versions = get_supported_versions()

            @test "1.0.0" in versions
            @test "1.1.0" in versions
            @test issorted(versions)
        end

        @testset "Version compatibility matrix" begin
            matrix = get_version_compatibility_matrix()

            @test haskey(matrix, "1.0.0")
            @test "1.1.0" in matrix["1.0.0"]
            @test "1.0.0" in matrix["1.0.0"]  # Can migrate to self
        end
    end

end
