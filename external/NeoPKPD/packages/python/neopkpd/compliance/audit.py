"""
Audit trail for FDA 21 CFR Part 11 compliance.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Optional

from ..bridge import get_julia


class AuditAction(Enum):
    """Enumeration of audit trail action types."""

    CREATE = 1
    REPLAY = 2
    VALIDATE = 3
    MODIFY = 4


@dataclass(frozen=True)
class AuditRecord:
    """
    Complete audit record for FDA 21 CFR Part 11 compliance.

    Attributes:
        execution_id: UUID v4 uniquely identifying this execution
        timestamp_utc: ISO 8601 UTC timestamp
        timezone_offset: Local timezone offset (e.g., "-05:00")
        system_user: Operating system username
        hostname: Machine hostname
        execution_context: Context ("interactive", "ci", "batch", "jupyter")
        action: Type of action performed
        previous_execution_id: For replays, links to original execution
        neopkpd_version: NeoPKPD version used
        record_hash: SHA-256 hash of this record
    """

    execution_id: str
    timestamp_utc: str
    timezone_offset: str
    system_user: str
    hostname: str
    execution_context: str
    action: AuditAction
    previous_execution_id: Optional[str]
    neopkpd_version: str
    record_hash: str


def _action_from_string(s: str) -> AuditAction:
    """Convert string to AuditAction enum."""
    mapping = {
        "create": AuditAction.CREATE,
        "replay": AuditAction.REPLAY,
        "validate": AuditAction.VALIDATE,
        "modify": AuditAction.MODIFY,
    }
    return mapping.get(s.lower(), AuditAction.CREATE)


def create_audit_record(
    action: AuditAction = AuditAction.CREATE,
    previous_execution_id: Optional[str] = None,
) -> AuditRecord:
    """
    Create a new audit record for the current execution.

    Args:
        action: The type of action being performed.
        previous_execution_id: For replays, the ID of the original execution.

    Returns:
        AuditRecord: A new audit record with UUID, timestamp, and system info.

    Example:
        >>> record = create_audit_record()
        >>> print(record.execution_id)
        550e8400-e29b-41d4-a716-446655440000
    """
    jl = get_julia()

    # Map Python enum to Julia enum
    julia_actions = {
        AuditAction.CREATE: jl.NeoPKPD.AUDIT_CREATE,
        AuditAction.REPLAY: jl.NeoPKPD.AUDIT_REPLAY,
        AuditAction.VALIDATE: jl.NeoPKPD.AUDIT_VALIDATE,
        AuditAction.MODIFY: jl.NeoPKPD.AUDIT_MODIFY,
    }

    julia_action = julia_actions[action]
    record = jl.NeoPKPD.create_audit_record(
        action=julia_action,
        previous_execution_id=previous_execution_id,
    )

    return AuditRecord(
        execution_id=str(record.execution_id),
        timestamp_utc=str(record.timestamp_utc),
        timezone_offset=str(record.timezone_offset),
        system_user=str(record.system_user),
        hostname=str(record.hostname),
        execution_context=str(record.execution_context),
        action=_action_from_string(str(record.action)),
        previous_execution_id=(
            str(record.previous_execution_id)
            if record.previous_execution_id is not None
            else None
        ),
        neopkpd_version=str(record.neopkpd_version),
        record_hash=str(record.record_hash),
    )


def verify_audit_record(record: AuditRecord) -> bool:
    """
    Verify the integrity of an audit record by recomputing its hash.

    Args:
        record: The audit record to verify.

    Returns:
        True if the record's hash is valid.

    Example:
        >>> record = create_audit_record()
        >>> assert verify_audit_record(record)
    """
    jl = get_julia()

    # Reconstruct Julia record
    julia_actions = {
        AuditAction.CREATE: jl.NeoPKPD.AUDIT_CREATE,
        AuditAction.REPLAY: jl.NeoPKPD.AUDIT_REPLAY,
        AuditAction.VALIDATE: jl.NeoPKPD.AUDIT_VALIDATE,
        AuditAction.MODIFY: jl.NeoPKPD.AUDIT_MODIFY,
    }

    julia_record = jl.NeoPKPD.AuditRecord(
        record.execution_id,
        record.timestamp_utc,
        record.timezone_offset,
        record.system_user,
        record.hostname,
        record.execution_context,
        julia_actions[record.action],
        record.previous_execution_id,
        record.neopkpd_version,
        record.record_hash,
    )

    return bool(jl.NeoPKPD.verify_audit_record(julia_record))


def get_execution_id(artifact: Dict[str, Any]) -> Optional[str]:
    """
    Extract the execution ID from an artifact.

    Args:
        artifact: The artifact dictionary.

    Returns:
        The execution ID if present, None otherwise.

    Example:
        >>> with open("simulation.json") as f:
        ...     artifact = json.load(f)
        >>> exec_id = get_execution_id(artifact)
        >>> print(exec_id)
        550e8400-e29b-41d4-a716-446655440000
    """
    if "compliance_metadata" not in artifact:
        return None

    meta = artifact["compliance_metadata"]
    if "audit_record" not in meta:
        return None

    return meta["audit_record"]["execution_id"]


def get_artifact_timestamp(artifact: Dict[str, Any]) -> Optional[str]:
    """
    Extract the creation timestamp from an artifact.

    Args:
        artifact: The artifact dictionary.

    Returns:
        The ISO 8601 UTC timestamp if present, None otherwise.

    Example:
        >>> ts = get_artifact_timestamp(artifact)
        >>> print(ts)
        2024-01-10T15:30:00.123Z
    """
    if "compliance_metadata" not in artifact:
        return None

    meta = artifact["compliance_metadata"]
    if "audit_record" not in meta:
        return None

    return meta["audit_record"]["timestamp_utc"]


def get_artifact_user(artifact: Dict[str, Any]) -> Optional[str]:
    """
    Extract the system user who created the artifact.

    Args:
        artifact: The artifact dictionary.

    Returns:
        The username if present, None otherwise.

    Example:
        >>> user = get_artifact_user(artifact)
        >>> print(user)
        analyst
    """
    if "compliance_metadata" not in artifact:
        return None

    meta = artifact["compliance_metadata"]
    if "audit_record" not in meta:
        return None

    return meta["audit_record"]["system_user"]
