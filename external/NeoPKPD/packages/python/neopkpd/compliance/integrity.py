"""
Data integrity via SHA-256 hashing for FDA 21 CFR Part 11 compliance.
"""

import hashlib
import json
from dataclasses import dataclass
from typing import Any, Dict, Optional

from ..bridge import get_julia


@dataclass(frozen=True)
class IntegrityMetadata:
    """
    Cryptographic integrity metadata for artifact validation.

    Attributes:
        content_hash: SHA-256 hash of the entire artifact content
        hash_algorithm: Hash algorithm used (always "SHA-256")
        input_hash: SHA-256 hash of input parameters only
        output_hash: SHA-256 hash of output results only
        semantics_hash: SHA-256 hash of semantics fingerprint
        hash_computed_at: ISO 8601 UTC timestamp
    """

    content_hash: str
    hash_algorithm: str
    input_hash: str
    output_hash: str
    semantics_hash: str
    hash_computed_at: str


@dataclass(frozen=True)
class IntegrityVerificationResult:
    """
    Result of verifying artifact integrity.

    Attributes:
        is_valid: Whether the artifact passes integrity verification
        content_hash_valid: Whether the content hash matches
        input_hash_valid: Whether the input hash matches
        output_hash_valid: Whether the output hash matches
        semantics_hash_valid: Whether the semantics hash matches
        computed_content_hash: The computed content hash
        stored_content_hash: The stored content hash (if available)
        error_message: Error message if verification failed
    """

    is_valid: bool
    content_hash_valid: bool
    input_hash_valid: bool
    output_hash_valid: bool
    semantics_hash_valid: bool
    computed_content_hash: str
    stored_content_hash: Optional[str]
    error_message: Optional[str]


def compute_content_hash(content: str) -> str:
    """
    Compute SHA-256 hash of a string.

    Args:
        content: The string to hash.

    Returns:
        64-character hexadecimal SHA-256 hash.

    Example:
        >>> hash = compute_content_hash("test data")
        >>> len(hash)
        64
    """
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def has_compliance_metadata(artifact: Dict[str, Any]) -> bool:
    """
    Check if an artifact has a compliance_metadata block.

    Args:
        artifact: The artifact dictionary.

    Returns:
        True if compliance_metadata exists.

    Example:
        >>> if has_compliance_metadata(artifact):
        ...     print("This artifact has compliance tracking")
    """
    return "compliance_metadata" in artifact


def verify_artifact_integrity(artifact: Dict[str, Any]) -> IntegrityVerificationResult:
    """
    Verify the integrity of an artifact by recomputing and comparing hashes.

    For artifacts without compliance_metadata (schema < 1.1.0), returns valid=True
    with a note that verification was not performed.

    Args:
        artifact: The artifact dictionary to verify.

    Returns:
        IntegrityVerificationResult with validation status.

    Example:
        >>> with open("simulation.json") as f:
        ...     artifact = json.load(f)
        >>> result = verify_artifact_integrity(artifact)
        >>> if result.is_valid:
        ...     print("Artifact integrity verified")
        >>> else:
        ...     print(f"Integrity check failed: {result.error_message}")
    """
    jl = get_julia()

    # Convert Python dict to Julia Dict
    julia_artifact = jl.Dict(artifact)

    result = jl.NeoPKPD.verify_artifact_integrity(julia_artifact)

    return IntegrityVerificationResult(
        is_valid=bool(result.is_valid),
        content_hash_valid=bool(result.content_hash_valid),
        input_hash_valid=bool(result.input_hash_valid),
        output_hash_valid=bool(result.output_hash_valid),
        semantics_hash_valid=bool(result.semantics_hash_valid),
        computed_content_hash=str(result.computed_content_hash),
        stored_content_hash=(
            str(result.stored_content_hash)
            if result.stored_content_hash is not None
            else None
        ),
        error_message=(
            str(result.error_message)
            if result.error_message is not None
            else None
        ),
    )


def get_integrity_metadata(artifact: Dict[str, Any]) -> Optional[IntegrityMetadata]:
    """
    Extract integrity metadata from an artifact.

    Args:
        artifact: The artifact dictionary.

    Returns:
        IntegrityMetadata if present, None otherwise.

    Example:
        >>> meta = get_integrity_metadata(artifact)
        >>> if meta:
        ...     print(f"Content hash: {meta.content_hash}")
    """
    if "compliance_metadata" not in artifact:
        return None

    meta = artifact["compliance_metadata"]
    if "integrity" not in meta:
        return None

    integrity = meta["integrity"]

    return IntegrityMetadata(
        content_hash=integrity["content_hash"],
        hash_algorithm=integrity["hash_algorithm"],
        input_hash=integrity["input_hash"],
        output_hash=integrity["output_hash"],
        semantics_hash=integrity["semantics_hash"],
        hash_computed_at=integrity["hash_computed_at"],
    )
