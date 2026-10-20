# Compliance metadata deserialization
# Parses compliance_metadata block from artifacts

export deserialize_compliance_metadata, extract_compliance_metadata
export has_compliance_metadata, get_execution_id

"""
    has_compliance_metadata(artifact::Dict) -> Bool

Check if an artifact has a compliance_metadata block.
"""
function has_compliance_metadata(artifact::Dict)::Bool
    return haskey(artifact, "compliance_metadata")
end

"""
    get_execution_id(artifact::Dict) -> Union{Nothing, String}

Extract the execution ID from an artifact's compliance metadata.
Returns nothing if no compliance metadata or audit record exists.
"""
function get_execution_id(artifact::Dict)::Union{Nothing,String}
    if !has_compliance_metadata(artifact)
        return nothing
    end

    meta = artifact["compliance_metadata"]
    if !haskey(meta, "audit_record")
        return nothing
    end

    return String(meta["audit_record"]["execution_id"])
end

"""
    deserialize_compliance_metadata(d::Dict) -> ComplianceMetadata

Deserialize a compliance_metadata dictionary block.
"""
function deserialize_compliance_metadata(d::Dict)::ComplianceMetadata
    audit_record = if haskey(d, "audit_record")
        deserialize_audit_record(d["audit_record"])
    else
        nothing
    end

    integrity = if haskey(d, "integrity")
        deserialize_integrity(d["integrity"])
    else
        nothing
    end

    environment = if haskey(d, "environment")
        deserialize_environment(d["environment"])
    else
        nothing
    end

    return ComplianceMetadata(audit_record, integrity, environment)
end

"""
    extract_compliance_metadata(artifact::Dict) -> Union{Nothing, ComplianceMetadata}

Extract and deserialize compliance metadata from an artifact.
Returns nothing if no compliance metadata exists.
"""
function extract_compliance_metadata(artifact::Dict)::Union{Nothing,ComplianceMetadata}
    if !has_compliance_metadata(artifact)
        return nothing
    end

    return deserialize_compliance_metadata(artifact["compliance_metadata"])
end

"""
    get_artifact_timestamp(artifact::Dict) -> Union{Nothing, String}

Extract the creation timestamp from an artifact.
Returns the audit record timestamp if available, nothing otherwise.
"""
function get_artifact_timestamp(artifact::Dict)::Union{Nothing,String}
    if !has_compliance_metadata(artifact)
        return nothing
    end

    meta = artifact["compliance_metadata"]
    if !haskey(meta, "audit_record")
        return nothing
    end

    return String(meta["audit_record"]["timestamp_utc"])
end

"""
    get_artifact_user(artifact::Dict) -> Union{Nothing, String}

Extract the system user who created the artifact.
Returns nothing if no compliance metadata exists.
"""
function get_artifact_user(artifact::Dict)::Union{Nothing,String}
    if !has_compliance_metadata(artifact)
        return nothing
    end

    meta = artifact["compliance_metadata"]
    if !haskey(meta, "audit_record")
        return nothing
    end

    return String(meta["audit_record"]["system_user"])
end

"""
    get_artifact_environment(artifact::Dict) -> Union{Nothing, EnvironmentSnapshot}

Extract the environment snapshot from an artifact.
Returns nothing if no compliance metadata or environment exists.
"""
function get_artifact_environment(artifact::Dict)::Union{Nothing,EnvironmentSnapshot}
    if !has_compliance_metadata(artifact)
        return nothing
    end

    meta = artifact["compliance_metadata"]
    if !haskey(meta, "environment")
        return nothing
    end

    return deserialize_environment(meta["environment"])
end

"""
    compare_environments(artifact1::Dict, artifact2::Dict) -> Dict{String, Tuple{Any, Any}}

Compare the environments of two artifacts.
Returns a dictionary of fields that differ, with (artifact1_value, artifact2_value) tuples.
"""
function compare_environments(artifact1::Dict, artifact2::Dict)::Dict{String,Tuple{Any,Any}}
    env1 = get_artifact_environment(artifact1)
    env2 = get_artifact_environment(artifact2)

    if env1 === nothing || env2 === nothing
        return Dict{String,Tuple{Any,Any}}()
    end

    differences = Dict{String,Tuple{Any,Any}}()

    for field in fieldnames(EnvironmentSnapshot)
        v1 = getfield(env1, field)
        v2 = getfield(env2, field)
        if v1 != v2
            differences[string(field)] = (v1, v2)
        end
    end

    return differences
end

export get_artifact_timestamp, get_artifact_user, get_artifact_environment
export compare_environments
