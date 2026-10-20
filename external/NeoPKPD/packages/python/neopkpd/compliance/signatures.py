"""
Digital Signatures for FDA 21 CFR Part 11 Compliance.

Implements ECDSA-P256 and RSA-2048 signing for artifact authentication.

FDA 21 CFR 11.70(a): "Signed electronic records shall contain information
associated with the signing that clearly indicates... the meaning of the signature"
"""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

from ..bridge import get_julia


class SignatureAlgorithm(Enum):
    """Supported digital signature algorithms."""

    ECDSA_P256 = "ECDSA-P256"
    RSA_2048 = "RSA-2048"


class SignaturePurpose(Enum):
    """
    The meaning/intent of a signature per FDA 21 CFR 11.70(a).

    Attributes:
        AUTHORSHIP: "I created this artifact"
        APPROVAL: "I approve this artifact"
        REVIEW: "I have reviewed this artifact"
        VERIFICATION: "I verify this artifact is accurate"
        RESPONSIBILITY: "I am responsible for this artifact"
    """

    AUTHORSHIP = "authorship"
    APPROVAL = "approval"
    REVIEW = "review"
    VERIFICATION = "verification"
    RESPONSIBILITY = "responsibility"


@dataclass(frozen=True)
class SigningKey:
    """
    Private key for signing artifacts. Must be kept secure.

    Attributes:
        algorithm: The signature algorithm
        key_id: Unique identifier for the key (UUID)
        key_data_hex: The private key bytes as hex string
        created_at: ISO 8601 UTC timestamp of key creation
        owner: Owner/organization of the key
    """

    algorithm: SignatureAlgorithm
    key_id: str
    key_data_hex: str
    created_at: str
    owner: str


@dataclass(frozen=True)
class VerificationKey:
    """
    Public key for verifying signatures. Can be freely distributed.

    Attributes:
        algorithm: The signature algorithm
        key_id: Unique identifier matching the signing key
        key_data_hex: The public key bytes as hex string
        created_at: ISO 8601 UTC timestamp of key creation
        owner: Owner/organization of the key
    """

    algorithm: SignatureAlgorithm
    key_id: str
    key_data_hex: str
    created_at: str
    owner: str


@dataclass(frozen=True)
class KeyPair:
    """A matched pair of signing and verification keys."""

    signing_key: SigningKey
    verification_key: VerificationKey


@dataclass(frozen=True)
class ArtifactSignature:
    """
    A digital signature on an artifact with metadata.

    Attributes:
        signature_id: UUID of this signature
        key_id: ID of the key used for signing
        algorithm: Algorithm used
        signature_hex: The actual signature as hex string
        signed_hash: SHA-256 hash of the signed content
        purpose: The meaning of this signature
        signer_name: Name of the signer
        signer_title: Title/role of the signer
        signer_organization: Organization of the signer
        signed_at: ISO 8601 UTC timestamp
        comments: Optional comments about the signature
    """

    signature_id: str
    key_id: str
    algorithm: SignatureAlgorithm
    signature_hex: str
    signed_hash: str
    purpose: SignaturePurpose
    signer_name: str
    signer_title: str
    signer_organization: str
    signed_at: str
    comments: str


@dataclass(frozen=True)
class SignatureVerificationResult:
    """
    Result of verifying a signature.

    Attributes:
        is_valid: Whether the signature is valid
        signature: The signature that was verified
        error_message: Error message if verification failed
    """

    is_valid: bool
    signature: ArtifactSignature
    error_message: Optional[str]


def generate_keypair(
    algorithm: SignatureAlgorithm = SignatureAlgorithm.ECDSA_P256,
    owner: str = "",
) -> KeyPair:
    """
    Generate a new key pair for signing artifacts.

    Args:
        algorithm: The signature algorithm to use (default: ECDSA_P256)
        owner: The owner/organization of the key

    Returns:
        A KeyPair containing matched signing and verification keys.

    Example:
        >>> keypair = generate_keypair(owner="Pharma Corp")
        >>> # Store keypair.signing_key securely
        >>> # Distribute keypair.verification_key freely
    """
    jl = get_julia()

    # Map enum to Julia enum
    jl_algorithm = (
        jl.NeoPKPD.ECDSA_P256
        if algorithm == SignatureAlgorithm.ECDSA_P256
        else jl.NeoPKPD.RSA_2048
    )

    keypair = jl.NeoPKPD.generate_keypair(algorithm=jl_algorithm, owner=owner)

    signing_key = SigningKey(
        algorithm=algorithm,
        key_id=str(keypair.signing_key.key_id),
        key_data_hex=jl.bytes2hex(keypair.signing_key.key_data),
        created_at=str(keypair.signing_key.created_at),
        owner=str(keypair.signing_key.owner),
    )

    verification_key = VerificationKey(
        algorithm=algorithm,
        key_id=str(keypair.verification_key.key_id),
        key_data_hex=jl.bytes2hex(keypair.verification_key.key_data),
        created_at=str(keypair.verification_key.created_at),
        owner=str(keypair.verification_key.owner),
    )

    return KeyPair(signing_key=signing_key, verification_key=verification_key)


def sign_artifact(
    artifact: Dict[str, Any],
    signing_key: SigningKey,
    purpose: SignaturePurpose = SignaturePurpose.AUTHORSHIP,
    signer_name: str = "",
    signer_title: str = "",
    signer_organization: str = "",
    comments: str = "",
) -> ArtifactSignature:
    """
    Sign an artifact with a private key.

    Args:
        artifact: The artifact dictionary to sign
        signing_key: The private signing key
        purpose: The meaning of this signature (default: AUTHORSHIP)
        signer_name: Name of the person signing
        signer_title: Title/role of the signer
        signer_organization: Organization of the signer
        comments: Optional comments about the signature

    Returns:
        An ArtifactSignature that can be added to the artifact.

    Example:
        >>> signature = sign_artifact(
        ...     artifact,
        ...     keypair.signing_key,
        ...     purpose=SignaturePurpose.APPROVAL,
        ...     signer_name="Dr. Jane Smith",
        ...     signer_title="Clinical Pharmacologist",
        ...     signer_organization="Pharma Corp"
        ... )
    """
    jl = get_julia()

    # Map Python SigningKey to Julia SigningKey
    jl_algorithm = (
        jl.NeoPKPD.ECDSA_P256
        if signing_key.algorithm == SignatureAlgorithm.ECDSA_P256
        else jl.NeoPKPD.RSA_2048
    )

    jl_signing_key = jl.NeoPKPD.SigningKey(
        jl_algorithm,
        signing_key.key_id,
        jl.hex2bytes(signing_key.key_data_hex),
        signing_key.created_at,
        signing_key.owner,
    )

    # Map purpose
    purpose_map = {
        SignaturePurpose.AUTHORSHIP: jl.NeoPKPD.SIGNATURE_AUTHORSHIP,
        SignaturePurpose.APPROVAL: jl.NeoPKPD.SIGNATURE_APPROVAL,
        SignaturePurpose.REVIEW: jl.NeoPKPD.SIGNATURE_REVIEW,
        SignaturePurpose.VERIFICATION: jl.NeoPKPD.SIGNATURE_VERIFICATION,
        SignaturePurpose.RESPONSIBILITY: jl.NeoPKPD.SIGNATURE_RESPONSIBILITY,
    }

    julia_artifact = jl.Dict(artifact)

    sig = jl.NeoPKPD.sign_artifact(
        julia_artifact,
        jl_signing_key,
        purpose=purpose_map[purpose],
        signer_name=signer_name,
        signer_title=signer_title,
        signer_organization=signer_organization or signing_key.owner,
        comments=comments,
    )

    return ArtifactSignature(
        signature_id=str(sig.signature_id),
        key_id=str(sig.key_id),
        algorithm=signing_key.algorithm,
        signature_hex=jl.bytes2hex(sig.signature_bytes),
        signed_hash=str(sig.signed_hash),
        purpose=purpose,
        signer_name=str(sig.signer_name),
        signer_title=str(sig.signer_title),
        signer_organization=str(sig.signer_organization),
        signed_at=str(sig.signed_at),
        comments=str(sig.comments),
    )


def verify_signature(
    artifact: Dict[str, Any],
    signature: ArtifactSignature,
    verification_key: VerificationKey,
) -> bool:
    """
    Verify a signature on an artifact.

    Args:
        artifact: The artifact dictionary that was signed
        signature: The signature to verify
        verification_key: The public verification key

    Returns:
        True if the signature is valid, False otherwise.

    Example:
        >>> is_valid = verify_signature(artifact, signature, keypair.verification_key)
        >>> if is_valid:
        ...     print(f"Signature verified: {signature.signer_name}")
    """
    jl = get_julia()

    # Map Python VerificationKey to Julia VerificationKey
    jl_algorithm = (
        jl.NeoPKPD.ECDSA_P256
        if verification_key.algorithm == SignatureAlgorithm.ECDSA_P256
        else jl.NeoPKPD.RSA_2048
    )

    jl_verification_key = jl.NeoPKPD.VerificationKey(
        jl_algorithm,
        verification_key.key_id,
        jl.hex2bytes(verification_key.key_data_hex),
        verification_key.created_at,
        verification_key.owner,
    )

    # Map purpose
    purpose_map = {
        SignaturePurpose.AUTHORSHIP: jl.NeoPKPD.SIGNATURE_AUTHORSHIP,
        SignaturePurpose.APPROVAL: jl.NeoPKPD.SIGNATURE_APPROVAL,
        SignaturePurpose.REVIEW: jl.NeoPKPD.SIGNATURE_REVIEW,
        SignaturePurpose.VERIFICATION: jl.NeoPKPD.SIGNATURE_VERIFICATION,
        SignaturePurpose.RESPONSIBILITY: jl.NeoPKPD.SIGNATURE_RESPONSIBILITY,
    }

    # Map Python ArtifactSignature to Julia ArtifactSignature
    jl_signature = jl.NeoPKPD.ArtifactSignature(
        signature.signature_id,
        signature.key_id,
        jl_algorithm,
        jl.hex2bytes(signature.signature_hex),
        signature.signed_hash,
        purpose_map[signature.purpose],
        signature.signer_name,
        signature.signer_title,
        signature.signer_organization,
        signature.signed_at,
        signature.comments,
    )

    julia_artifact = jl.Dict(artifact)

    return bool(
        jl.NeoPKPD.verify_signature(
            julia_artifact, jl_signature, jl_verification_key
        )
    )


def get_artifact_signatures(artifact: Dict[str, Any]) -> List[ArtifactSignature]:
    """
    Extract all signatures from an artifact.

    Args:
        artifact: The artifact dictionary.

    Returns:
        List of ArtifactSignature objects.

    Example:
        >>> signatures = get_artifact_signatures(artifact)
        >>> for sig in signatures:
        ...     print(f"{sig.signer_name}: {sig.purpose.value}")
    """
    if "compliance_metadata" not in artifact:
        return []

    meta = artifact["compliance_metadata"]
    if "signatures" not in meta:
        return []

    sigs_data = meta["signatures"].get("signatures", [])
    signatures = []

    for sig in sigs_data:
        algorithm = (
            SignatureAlgorithm.ECDSA_P256
            if sig["algorithm"] == "ECDSA-P256"
            else SignatureAlgorithm.RSA_2048
        )

        purpose_map = {
            "authorship": SignaturePurpose.AUTHORSHIP,
            "approval": SignaturePurpose.APPROVAL,
            "review": SignaturePurpose.REVIEW,
            "verification": SignaturePurpose.VERIFICATION,
            "responsibility": SignaturePurpose.RESPONSIBILITY,
        }

        signatures.append(
            ArtifactSignature(
                signature_id=sig["signature_id"],
                key_id=sig["key_id"],
                algorithm=algorithm,
                signature_hex=sig["signature"],
                signed_hash=sig["signed_hash"],
                purpose=purpose_map.get(sig["purpose"], SignaturePurpose.AUTHORSHIP),
                signer_name=sig["signer_name"],
                signer_title=sig["signer_title"],
                signer_organization=sig["signer_organization"],
                signed_at=sig["signed_at"],
                comments=sig.get("comments", ""),
            )
        )

    return signatures


def verify_all_signatures(
    artifact: Dict[str, Any],
    trusted_keys: List[VerificationKey],
) -> List[SignatureVerificationResult]:
    """
    Verify all signatures on an artifact against trusted keys.

    Args:
        artifact: The artifact dictionary.
        trusted_keys: List of trusted verification keys.

    Returns:
        List of SignatureVerificationResult for each signature.

    Example:
        >>> results = verify_all_signatures(artifact, [trusted_key1, trusted_key2])
        >>> for result in results:
        ...     status = "VALID" if result.is_valid else "INVALID"
        ...     print(f"{result.signature.signer_name}: {status}")
    """
    signatures = get_artifact_signatures(artifact)
    results = []

    for sig in signatures:
        # Find matching key
        matching_key = None
        for key in trusted_keys:
            if key.key_id == sig.key_id:
                matching_key = key
                break

        if matching_key is None:
            results.append(
                SignatureVerificationResult(
                    is_valid=False,
                    signature=sig,
                    error_message=f"No trusted key found for key_id: {sig.key_id}",
                )
            )
            continue

        is_valid = verify_signature(artifact, sig, matching_key)
        results.append(
            SignatureVerificationResult(
                is_valid=is_valid,
                signature=sig,
                error_message=None if is_valid else "Signature verification failed",
            )
        )

    return results


def export_public_key(key: VerificationKey) -> Dict[str, Any]:
    """
    Export a verification key to a shareable format.

    Args:
        key: The verification key to export.

    Returns:
        Dictionary containing the key data.

    Example:
        >>> key_data = export_public_key(keypair.verification_key)
        >>> with open("public_key.json", "w") as f:
        ...     json.dump(key_data, f)
    """
    return {
        "key_type": "verification",
        "key_id": key.key_id,
        "algorithm": key.algorithm.value,
        "key_data": key.key_data_hex,
        "created_at": key.created_at,
        "owner": key.owner,
    }


def import_public_key(data: Dict[str, Any]) -> VerificationKey:
    """
    Import a verification key from exported format.

    Args:
        data: Dictionary containing the key data.

    Returns:
        VerificationKey object.

    Example:
        >>> with open("public_key.json") as f:
        ...     key_data = json.load(f)
        >>> key = import_public_key(key_data)
    """
    algorithm = (
        SignatureAlgorithm.ECDSA_P256
        if data["algorithm"] == "ECDSA-P256"
        else SignatureAlgorithm.RSA_2048
    )

    return VerificationKey(
        algorithm=algorithm,
        key_id=data["key_id"],
        key_data_hex=data["key_data"],
        created_at=data["created_at"],
        owner=data["owner"],
    )
