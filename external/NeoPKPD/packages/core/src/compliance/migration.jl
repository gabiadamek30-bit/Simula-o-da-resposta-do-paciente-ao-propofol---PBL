# Schema Migration Framework for FDA 21 CFR Part 11 Compliance
# Handles artifact version upgrades with full audit trail
#
# GAMP 5: "...change control and configuration management..."

export MigrationResult, MigrationStep
export migrate_artifact, can_migrate, get_migration_path
export MIGRATION_REGISTRY, MIGRATION_FUNCTIONS

"""
    MigrationStep

A single step in the migration process.

# Fields
- `from_version::String`: Source schema version
- `to_version::String`: Target schema version
- `description::String`: Human-readable description of changes
- `breaking_changes::Vector{String}`: List of breaking changes
- `added_fields::Vector{String}`: New fields added
- `removed_fields::Vector{String}`: Fields removed
- `modified_fields::Vector{String}`: Fields with changed semantics
"""
struct MigrationStep
    from_version::String
    to_version::String
    description::String
    breaking_changes::Vector{String}
    added_fields::Vector{String}
    removed_fields::Vector{String}
    modified_fields::Vector{String}
end

"""
    MigrationResult

Result of migrating an artifact.

# Fields
- `success::Bool`: Whether migration succeeded
- `source_version::String`: Original schema version
- `target_version::String`: Target schema version
- `migrated_artifact::Union{Nothing, Dict}`: The migrated artifact (if successful)
- `steps_applied::Vector{MigrationStep}`: Migration steps that were applied
- `errors::Vector{String}`: Any errors encountered
- `warnings::Vector{String}`: Non-fatal warnings
- `migration_id::String`: UUID of this migration operation
- `migrated_at::String`: ISO 8601 UTC timestamp
"""
struct MigrationResult
    success::Bool
    source_version::String
    target_version::String
    migrated_artifact::Union{Nothing,Dict}
    steps_applied::Vector{MigrationStep}
    errors::Vector{String}
    warnings::Vector{String}
    migration_id::String
    migrated_at::String
end

# =============================================================================
# Migration Registry
# =============================================================================

"""
Migration functions registry: (from_version, to_version) => migration_function
"""
const MIGRATION_FUNCTIONS = Dict{Tuple{String,String},Function}()

"""
Migration metadata registry
"""
const MIGRATION_REGISTRY = Dict{Tuple{String,String},MigrationStep}()

# -----------------------------------------------------------------------------
# Migration: 1.0.0 -> 1.1.0
# -----------------------------------------------------------------------------

const MIGRATION_1_0_0_TO_1_1_0 = MigrationStep(
    "1.0.0",
    "1.1.0",
    "Added FDA 21 CFR Part 11 compliance metadata block",
    String[],  # No breaking changes
    [
        "compliance_metadata",
        "compliance_metadata.audit_record",
        "compliance_metadata.integrity",
        "compliance_metadata.environment"
    ],
    String[],  # No removed fields
    String[]   # No modified fields
)

"""
    _migrate_1_0_0_to_1_1_0(artifact::Dict) -> Dict

Migrate artifact from schema 1.0.0 to 1.1.0.
Adds compliance_metadata block with audit trail, integrity, and environment.
"""
function _migrate_1_0_0_to_1_1_0(artifact::Dict)::Dict
    # Deep copy to avoid modifying original
    migrated = deepcopy(artifact)

    # Update schema version
    migrated["artifact_schema_version"] = "1.1.0"

    # Add compliance metadata if not present
    if !haskey(migrated, "compliance_metadata")
        # Create audit record for migration
        record = create_audit_record(action=AUDIT_MODIFY)

        # Capture environment
        env = capture_environment()

        # Compute integrity (without compliance_metadata to avoid circularity)
        integrity = compute_integrity_metadata(migrated)

        # Build compliance metadata
        migrated["compliance_metadata"] = Dict{String,Any}(
            "audit_record" => serialize_audit_record(record),
            "integrity" => serialize_integrity(integrity),
            "environment" => serialize_environment(env),
            "migration_note" => "Migrated from schema 1.0.0 to 1.1.0"
        )
    end

    # Update semantics fingerprint if present
    if haskey(migrated, "semantics_fingerprint")
        migrated["semantics_fingerprint"]["artifact_schema_version"] = "1.1.0"
    end

    return migrated
end

# Register the migration
MIGRATION_REGISTRY[("1.0.0", "1.1.0")] = MIGRATION_1_0_0_TO_1_1_0
MIGRATION_FUNCTIONS[("1.0.0", "1.1.0")] = _migrate_1_0_0_to_1_1_0

# -----------------------------------------------------------------------------
# Future migrations can be added here
# -----------------------------------------------------------------------------

# Example: Migration 1.1.0 -> 1.2.0 (placeholder for future)
# const MIGRATION_1_1_0_TO_1_2_0 = MigrationStep(...)
# MIGRATION_REGISTRY[("1.1.0", "1.2.0")] = MIGRATION_1_1_0_TO_1_2_0
# MIGRATION_FUNCTIONS[("1.1.0", "1.2.0")] = _migrate_1_1_0_to_1_2_0

# =============================================================================
# Migration Path Finding
# =============================================================================

"""
    get_migration_path(from_version::String, to_version::String) -> Vector{MigrationStep}

Find the sequence of migration steps to go from source to target version.
Returns empty vector if no path exists.
"""
function get_migration_path(from_version::String, to_version::String)::Vector{MigrationStep}
    if from_version == to_version
        return MigrationStep[]
    end

    # Simple direct lookup first
    if haskey(MIGRATION_REGISTRY, (from_version, to_version))
        return [MIGRATION_REGISTRY[(from_version, to_version)]]
    end

    # BFS for multi-step path
    visited = Set{String}([from_version])
    queue = [(from_version, MigrationStep[])]

    while !isempty(queue)
        current_version, path = popfirst!(queue)

        # Find all migrations from current version
        for ((src, dst), step) in MIGRATION_REGISTRY
            if src == current_version && !(dst in visited)
                new_path = vcat(path, [step])

                if dst == to_version
                    return new_path
                end

                push!(visited, dst)
                push!(queue, (dst, new_path))
            end
        end
    end

    return MigrationStep[]  # No path found
end

"""
    can_migrate(from_version::String, to_version::String) -> Bool

Check if migration is possible between two versions.
"""
function can_migrate(from_version::String, to_version::String)::Bool
    if from_version == to_version
        return true
    end
    return !isempty(get_migration_path(from_version, to_version))
end

# =============================================================================
# Migration Execution
# =============================================================================

"""
    migrate_artifact(
        artifact::Dict;
        target_version::String = ARTIFACT_SCHEMA_VERSION
    ) -> MigrationResult

Migrate an artifact to the target schema version.

# Arguments
- `artifact`: The artifact to migrate
- `target_version`: Target schema version (default: current version)

# Returns
MigrationResult with migrated artifact or errors.

# Example
```julia
# Load old artifact
old_artifact = JSON.parsefile("old_simulation.json")

# Migrate to current version
result = migrate_artifact(old_artifact)

if result.success
    # Use migrated artifact
    new_artifact = result.migrated_artifact
    JSON.print("new_simulation.json", new_artifact)
else
    println("Migration failed: ", result.errors)
end
```
"""
function migrate_artifact(
    artifact::Dict;
    target_version::String=ARTIFACT_SCHEMA_VERSION
)::MigrationResult
    migration_id = string(uuid4())
    migrated_at = _utc_timestamp()
    errors = String[]
    warnings = String[]

    # Get source version
    source_version = get(artifact, "artifact_schema_version", "1.0.0")

    # Already at target?
    if source_version == target_version
        return MigrationResult(
            true,
            source_version,
            target_version,
            artifact,
            MigrationStep[],
            String[],
            ["Artifact already at target version"],
            migration_id,
            migrated_at
        )
    end

    # Find migration path
    path = get_migration_path(source_version, target_version)

    if isempty(path)
        return MigrationResult(
            false,
            source_version,
            target_version,
            nothing,
            MigrationStep[],
            ["No migration path from $source_version to $target_version"],
            String[],
            migration_id,
            migrated_at
        )
    end

    # Execute migrations
    current_artifact = deepcopy(artifact)
    applied_steps = MigrationStep[]

    for step in path
        key = (step.from_version, step.to_version)

        if !haskey(MIGRATION_FUNCTIONS, key)
            push!(errors, "Missing migration function for $(step.from_version) -> $(step.to_version)")
            return MigrationResult(
                false,
                source_version,
                target_version,
                nothing,
                applied_steps,
                errors,
                warnings,
                migration_id,
                migrated_at
            )
        end

        try
            migrate_fn = MIGRATION_FUNCTIONS[key]
            current_artifact = migrate_fn(current_artifact)
            push!(applied_steps, step)

            # Add warning for breaking changes
            for bc in step.breaking_changes
                push!(warnings, "Breaking change in $(step.from_version) -> $(step.to_version): $bc")
            end
        catch e
            push!(errors, "Migration $(step.from_version) -> $(step.to_version) failed: $(string(e))")
            return MigrationResult(
                false,
                source_version,
                target_version,
                nothing,
                applied_steps,
                errors,
                warnings,
                migration_id,
                migrated_at
            )
        end
    end

    # Verify final version
    final_version = get(current_artifact, "artifact_schema_version", "unknown")
    if final_version != target_version
        push!(errors, "Final version $final_version does not match target $target_version")
        return MigrationResult(
            false,
            source_version,
            target_version,
            nothing,
            applied_steps,
            errors,
            warnings,
            migration_id,
            migrated_at
        )
    end

    # Validate migrated artifact
    validation = validate_artifact_schema(current_artifact)
    if !validation.is_valid
        for err in validation.errors
            push!(warnings, "Validation warning: $(err.message) at $(err.path)")
        end
    end

    return MigrationResult(
        true,
        source_version,
        target_version,
        current_artifact,
        applied_steps,
        String[],
        warnings,
        migration_id,
        migrated_at
    )
end

"""
    serialize_migration_result(result::MigrationResult) -> Dict{String, Any}

Serialize a MigrationResult for logging/auditing.
"""
function serialize_migration_result(result::MigrationResult)::Dict{String,Any}
    steps_serialized = [
        Dict{String,Any}(
            "from_version" => s.from_version,
            "to_version" => s.to_version,
            "description" => s.description,
            "breaking_changes" => s.breaking_changes,
            "added_fields" => s.added_fields,
            "removed_fields" => s.removed_fields,
            "modified_fields" => s.modified_fields
        )
        for s in result.steps_applied
    ]

    return Dict{String,Any}(
        "success" => result.success,
        "source_version" => result.source_version,
        "target_version" => result.target_version,
        "steps_applied" => steps_serialized,
        "errors" => result.errors,
        "warnings" => result.warnings,
        "migration_id" => result.migration_id,
        "migrated_at" => result.migrated_at
    )
end

"""
    get_supported_versions() -> Vector{String}

Get list of all schema versions that can be migrated.
"""
function get_supported_versions()::Vector{String}
    versions = Set{String}()
    for (src, dst) in keys(MIGRATION_REGISTRY)
        push!(versions, src)
        push!(versions, dst)
    end
    return sort(collect(versions))
end

"""
    get_version_compatibility_matrix() -> Dict{String, Vector{String}}

Get a compatibility matrix showing which versions can migrate to which.
"""
function get_version_compatibility_matrix()::Dict{String,Vector{String}}
    versions = get_supported_versions()
    matrix = Dict{String,Vector{String}}()

    for src in versions
        targets = String[]
        for dst in versions
            if can_migrate(src, dst)
                push!(targets, dst)
            end
        end
        matrix[src] = targets
    end

    return matrix
end

export serialize_migration_result, get_supported_versions, get_version_compatibility_matrix
