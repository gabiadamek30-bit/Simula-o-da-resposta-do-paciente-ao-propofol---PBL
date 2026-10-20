# Compliance-aware serialization
# Extends base artifact serialization with compliance_metadata block

export add_compliance_metadata, create_compliance_metadata
export ComplianceMetadata, serialize_compliance_metadata

"""
    ComplianceMetadata

Complete compliance metadata block for an artifact.

# Fields
- `audit_record::Union{Nothing, AuditRecord}`: Audit trail record
- `integrity::Union{Nothing, IntegrityMetadata}`: Integrity hashes
- `environment::Union{Nothing, EnvironmentSnapshot}`: Environment snapshot
"""
struct ComplianceMetadata
    audit_record::Union{Nothing,AuditRecord}
    integrity::Union{Nothing,IntegrityMetadata}
    environment::Union{Nothing,EnvironmentSnapshot}
end

"""
    create_compliance_metadata(;
        action::AuditAction = AUDIT_CREATE,
        previous_execution_id::Union{Nothing,String} = nothing
    ) -> ComplianceMetadata

Create compliance metadata based on current configuration.

# Arguments
- `action`: The audit action being performed
- `previous_execution_id`: For replays, the original execution ID
"""
function create_compliance_metadata(;
    action::AuditAction=AUDIT_CREATE,
    previous_execution_id::Union{Nothing,String}=nothing
)::ComplianceMetadata
    config = get_compliance_config()

    audit_record = if config.enable_audit_trail
        create_audit_record(action=action, previous_execution_id=previous_execution_id)
    else
        nothing
    end

    environment = if config.enable_environment_capture
        capture_environment()
    else
        nothing
    end

    # Note: integrity is computed after artifact is created
    # Use compute_and_add_integrity! to add integrity hashes
    return ComplianceMetadata(audit_record, nothing, environment)
end

"""
    serialize_compliance_metadata(meta::ComplianceMetadata) -> Dict{String, Any}

Serialize ComplianceMetadata to a dictionary.
"""
function serialize_compliance_metadata(meta::ComplianceMetadata)::Dict{String,Any}
    result = Dict{String,Any}()

    if meta.audit_record !== nothing
        result["audit_record"] = serialize_audit_record(meta.audit_record)
    end

    if meta.integrity !== nothing
        result["integrity"] = serialize_integrity(meta.integrity)
    end

    if meta.environment !== nothing
        result["environment"] = serialize_environment(meta.environment)
    end

    return result
end

"""
    add_compliance_metadata!(artifact::Dict;
        action::AuditAction = AUDIT_CREATE,
        previous_execution_id::Union{Nothing,String} = nothing
    ) -> Dict

Add compliance_metadata block to an artifact in-place.
Computes integrity hashes based on current artifact content.
Returns the modified artifact.

# Arguments
- `artifact`: The artifact dictionary to modify
- `action`: The audit action being performed
- `previous_execution_id`: For replays, the original execution ID
"""
function add_compliance_metadata!(artifact::Dict;
    action::AuditAction=AUDIT_CREATE,
    previous_execution_id::Union{Nothing,String}=nothing
)::Dict
    config = get_compliance_config()

    if config.level == COMPLIANCE_DISABLED
        return artifact
    end

    meta = create_compliance_metadata(
        action=action,
        previous_execution_id=previous_execution_id
    )

    # Compute integrity if enabled
    integrity = if config.enable_hashing
        compute_integrity_metadata(artifact)
    else
        nothing
    end

    # Create full metadata with integrity
    full_meta = ComplianceMetadata(meta.audit_record, integrity, meta.environment)

    # Add to artifact
    artifact["compliance_metadata"] = serialize_compliance_metadata(full_meta)

    return artifact
end

"""
    add_compliance_metadata(artifact::Dict; kwargs...) -> Dict

Non-mutating version of add_compliance_metadata!.
Returns a new artifact with compliance_metadata added.
"""
function add_compliance_metadata(artifact::Dict; kwargs...)::Dict
    artifact_copy = deepcopy(artifact)
    return add_compliance_metadata!(artifact_copy; kwargs...)
end

export add_compliance_metadata!
