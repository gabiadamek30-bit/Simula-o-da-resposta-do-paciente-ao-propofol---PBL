# Audit trail for FDA 21 CFR Part 11 compliance
# Provides execution tracking with UUIDs, timestamps, and chain-of-custody

using UUIDs
using Printf

export AuditRecord, create_audit_record, serialize_audit_record, deserialize_audit_record
export AuditAction, AUDIT_CREATE, AUDIT_REPLAY, AUDIT_VALIDATE, AUDIT_MODIFY

"""
    AuditAction

Enumeration of audit trail action types.
"""
@enum AuditAction begin
    AUDIT_CREATE = 1    # New artifact creation
    AUDIT_REPLAY = 2    # Replay of existing artifact
    AUDIT_VALIDATE = 3  # Validation/verification
    AUDIT_MODIFY = 4    # Modification of existing artifact
end

"""
    _action_to_string(action::AuditAction) -> String

Convert AuditAction enum to string for serialization.
"""
function _action_to_string(action::AuditAction)::String
    if action == AUDIT_CREATE
        return "create"
    elseif action == AUDIT_REPLAY
        return "replay"
    elseif action == AUDIT_VALIDATE
        return "validate"
    elseif action == AUDIT_MODIFY
        return "modify"
    else
        return "unknown"
    end
end

"""
    _string_to_action(s::String) -> AuditAction

Convert string to AuditAction enum for deserialization.
"""
function _string_to_action(s::String)::AuditAction
    if s == "create"
        return AUDIT_CREATE
    elseif s == "replay"
        return AUDIT_REPLAY
    elseif s == "validate"
        return AUDIT_VALIDATE
    elseif s == "modify"
        return AUDIT_MODIFY
    else
        return AUDIT_CREATE  # Default
    end
end

"""
    AuditRecord

Complete audit record for FDA 21 CFR Part 11 compliance.
Tracks who, what, when, and where for every execution.

# Fields
- `execution_id::String`: UUID v4 uniquely identifying this execution
- `timestamp_utc::String`: ISO 8601 UTC timestamp
- `timezone_offset::String`: Local timezone offset (e.g., "-05:00")
- `system_user::String`: Operating system username
- `hostname::String`: Machine hostname
- `execution_context::String`: Context ("interactive", "ci", "batch", "jupyter")
- `action::AuditAction`: Type of action performed
- `previous_execution_id::Union{Nothing,String}`: For replays, links to original execution
- `neopkpd_version::String`: NeoPKPD version used
- `record_hash::String`: SHA-256 hash of this record (excluding record_hash itself)
"""
struct AuditRecord
    execution_id::String
    timestamp_utc::String
    timezone_offset::String
    system_user::String
    hostname::String
    execution_context::String
    action::AuditAction
    previous_execution_id::Union{Nothing,String}
    neopkpd_version::String
    record_hash::String
end

"""
    _get_system_user() -> String

Get the current system username.
"""
function _get_system_user()::String
    try
        if haskey(ENV, "USER")
            return ENV["USER"]
        elseif haskey(ENV, "USERNAME")
            return ENV["USERNAME"]
        else
            return "unknown"
        end
    catch
        return "unknown"
    end
end

"""
    _get_hostname() -> String

Get the machine hostname.
"""
function _get_hostname()::String
    try
        return gethostname()
    catch
        return "unknown"
    end
end

"""
    _get_timezone_offset() -> String

Get the local timezone offset in ISO 8601 format (e.g., "-05:00").
"""
function _get_timezone_offset()::String
    try
        # Get offset in seconds
        now_local = Dates.now()
        now_utc = Dates.now(Dates.UTC)
        offset_seconds = Dates.value(now_local - now_utc) / 1000

        # Convert to hours and minutes
        offset_hours = floor(Int, offset_seconds / 3600)
        offset_minutes = abs(floor(Int, (offset_seconds % 3600) / 60))

        # Format as ±HH:MM
        sign = offset_hours >= 0 ? "+" : "-"
        return @sprintf("%s%02d:%02d", sign, abs(offset_hours), offset_minutes)
    catch
        return "+00:00"
    end
end

"""
    _detect_execution_context() -> String

Detect the execution context (interactive, CI, batch, jupyter).
"""
function _detect_execution_context()::String
    # Check for CI environments
    if haskey(ENV, "CI") || haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI") ||
       haskey(ENV, "TRAVIS") || haskey(ENV, "CIRCLECI") || haskey(ENV, "JENKINS_URL")
        return "ci"
    end

    # Check for Jupyter
    if haskey(ENV, "JPY_PARENT_PID") || haskey(ENV, "JUPYTER_RUNTIME_DIR")
        return "jupyter"
    end

    # Check for batch/non-interactive
    if !isinteractive()
        return "batch"
    end

    return "interactive"
end

"""
    _compute_record_hash(record_data::Dict) -> String

Compute SHA-256 hash of audit record data (excluding the hash field itself).
"""
function _compute_record_hash(record_data::Dict)::String
    # Remove record_hash if present
    data_for_hash = Dict{String,Any}()
    for (k, v) in record_data
        if k != "record_hash"
            data_for_hash[k] = v
        end
    end
    return compute_content_hash(data_for_hash)
end

"""
    create_audit_record(;
        action::AuditAction = AUDIT_CREATE,
        previous_execution_id::Union{Nothing,String} = nothing
    ) -> AuditRecord

Create a new audit record for the current execution.
Automatically captures system information, timestamp, and generates UUID.

# Arguments
- `action`: The type of action being performed
- `previous_execution_id`: For replays, the ID of the original execution

# Example
```julia
# For new execution
record = create_audit_record()

# For replay
record = create_audit_record(
    action = AUDIT_REPLAY,
    previous_execution_id = "550e8400-e29b-41d4-a716-446655440000"
)
```
"""
function create_audit_record(;
    action::AuditAction=AUDIT_CREATE,
    previous_execution_id::Union{Nothing,String}=nothing
)::AuditRecord
    execution_id = string(uuid4())
    timestamp_utc = _utc_timestamp()
    timezone_offset = _get_timezone_offset()
    system_user = _get_system_user()
    hostname = _get_hostname()
    execution_context = _detect_execution_context()

    # Compute hash of record data
    record_data = Dict{String,Any}(
        "execution_id" => execution_id,
        "timestamp_utc" => timestamp_utc,
        "timezone_offset" => timezone_offset,
        "system_user" => system_user,
        "hostname" => hostname,
        "execution_context" => execution_context,
        "action" => _action_to_string(action),
        "previous_execution_id" => previous_execution_id,
        "neopkpd_version" => NEOPKPD_VERSION
    )
    record_hash = _compute_record_hash(record_data)

    return AuditRecord(
        execution_id,
        timestamp_utc,
        timezone_offset,
        system_user,
        hostname,
        execution_context,
        action,
        previous_execution_id,
        NEOPKPD_VERSION,
        record_hash
    )
end

"""
    serialize_audit_record(record::AuditRecord) -> Dict{String, Any}

Serialize an AuditRecord to a dictionary for JSON encoding.
"""
function serialize_audit_record(record::AuditRecord)::Dict{String,Any}
    return Dict{String,Any}(
        "execution_id" => record.execution_id,
        "timestamp_utc" => record.timestamp_utc,
        "timezone_offset" => record.timezone_offset,
        "system_user" => record.system_user,
        "hostname" => record.hostname,
        "execution_context" => record.execution_context,
        "action" => _action_to_string(record.action),
        "previous_execution_id" => record.previous_execution_id,
        "neopkpd_version" => record.neopkpd_version,
        "record_hash" => record.record_hash
    )
end

"""
    deserialize_audit_record(d::Dict) -> AuditRecord

Deserialize an AuditRecord from a dictionary.
"""
function deserialize_audit_record(d::Dict)::AuditRecord
    return AuditRecord(
        String(d["execution_id"]),
        String(d["timestamp_utc"]),
        String(d["timezone_offset"]),
        String(d["system_user"]),
        String(d["hostname"]),
        String(d["execution_context"]),
        _string_to_action(String(d["action"])),
        d["previous_execution_id"] === nothing ? nothing : String(d["previous_execution_id"]),
        String(d["neopkpd_version"]),
        String(d["record_hash"])
    )
end

"""
    verify_audit_record(record::AuditRecord) -> Bool

Verify the integrity of an audit record by recomputing its hash.
"""
function verify_audit_record(record::AuditRecord)::Bool
    record_data = Dict{String,Any}(
        "execution_id" => record.execution_id,
        "timestamp_utc" => record.timestamp_utc,
        "timezone_offset" => record.timezone_offset,
        "system_user" => record.system_user,
        "hostname" => record.hostname,
        "execution_context" => record.execution_context,
        "action" => _action_to_string(record.action),
        "previous_execution_id" => record.previous_execution_id,
        "neopkpd_version" => record.neopkpd_version
    )
    computed_hash = _compute_record_hash(record_data)
    return computed_hash == record.record_hash
end

export verify_audit_record
