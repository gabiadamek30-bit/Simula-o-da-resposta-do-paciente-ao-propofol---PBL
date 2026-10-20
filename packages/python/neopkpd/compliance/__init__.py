"""
FDA 21 CFR Part 11 Compliance Module for NeoPKPD

This module provides regulatory compliance features aligned with:
- FDA 21 CFR Part 11: Electronic Records and Signatures
- GAMP 5: Good Automated Manufacturing Practice
- ALCOA+ Principles

Features:
- Audit trails with execution IDs and timestamps
- Data integrity via SHA-256 checksums
- Environment capture for reproducibility
- Validation documentation (IQ/OQ/PQ)

Example:
    >>> from neopkpd.compliance import (
    ...     verify_artifact_integrity,
    ...     get_environment_snapshot,
    ...     get_execution_id
    ... )
    >>>
    >>> # Verify artifact integrity
    >>> with open("simulation.json") as f:
    ...     artifact = json.load(f)
    >>> result = verify_artifact_integrity(artifact)
    >>> assert result.is_valid
    >>>
    >>> # Get environment info
    >>> env = get_environment_snapshot()
    >>> print(f"Julia: {env.julia_version}")
"""

from .audit import (
    AuditRecord,
    create_audit_record,
    verify_audit_record,
    get_execution_id,
    get_artifact_timestamp,
    get_artifact_user,
)

from .integrity import (
    IntegrityMetadata,
    IntegrityVerificationResult,
    verify_artifact_integrity,
    compute_content_hash,
    has_compliance_metadata,
)

from .environment import (
    EnvironmentSnapshot,
    get_environment_snapshot,
    get_artifact_environment,
    compare_environments,
)

from .config import (
    ComplianceLevel,
    ComplianceConfig,
    get_compliance_config,
    set_compliance_config,
)

from .signatures import (
    SignatureAlgorithm,
    SignaturePurpose,
    SigningKey,
    VerificationKey,
    KeyPair,
    ArtifactSignature,
    SignatureVerificationResult,
    generate_keypair,
    sign_artifact,
    verify_signature,
    get_artifact_signatures,
    verify_all_signatures,
    export_public_key,
    import_public_key,
)

__all__ = [
    # Audit
    "AuditRecord",
    "create_audit_record",
    "verify_audit_record",
    "get_execution_id",
    "get_artifact_timestamp",
    "get_artifact_user",
    # Integrity
    "IntegrityMetadata",
    "IntegrityVerificationResult",
    "verify_artifact_integrity",
    "compute_content_hash",
    "has_compliance_metadata",
    # Environment
    "EnvironmentSnapshot",
    "get_environment_snapshot",
    "get_artifact_environment",
    "compare_environments",
    # Config
    "ComplianceLevel",
    "ComplianceConfig",
    "get_compliance_config",
    "set_compliance_config",
    # Signatures
    "SignatureAlgorithm",
    "SignaturePurpose",
    "SigningKey",
    "VerificationKey",
    "KeyPair",
    "ArtifactSignature",
    "SignatureVerificationResult",
    "generate_keypair",
    "sign_artifact",
    "verify_signature",
    "get_artifact_signatures",
    "verify_all_signatures",
    "export_public_key",
    "import_public_key",
]
