"""
Compliance configuration for FDA 21 CFR Part 11 alignment.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Optional

from ..bridge import get_julia


class ComplianceLevel(Enum):
    """Enumeration of compliance levels."""

    DISABLED = 0  # No compliance overhead (development)
    MINIMAL = 1  # Timestamps and versions only
    STANDARD = 2  # Full audit trail, hashing, environment
    STRICT = 3  # Standard plus validation reports


@dataclass(frozen=True)
class ComplianceConfig:
    """Configuration for compliance features."""

    level: ComplianceLevel
    enable_audit_trail: bool
    enable_hashing: bool
    enable_environment_capture: bool
    enable_validation_reports: bool


def get_compliance_config() -> ComplianceConfig:
    """
    Get the current global compliance configuration.

    Returns:
        ComplianceConfig: The current compliance configuration.

    Example:
        >>> config = get_compliance_config()
        >>> print(config.level)
        ComplianceLevel.STANDARD
    """
    jl = get_julia()
    config = jl.NeoPKPD.get_compliance_config()

    # Map Julia enum to Python enum
    level_value = int(jl.Int(config.level))
    level = ComplianceLevel(level_value)

    return ComplianceConfig(
        level=level,
        enable_audit_trail=bool(config.enable_audit_trail),
        enable_hashing=bool(config.enable_hashing),
        enable_environment_capture=bool(config.enable_environment_capture),
        enable_validation_reports=bool(config.enable_validation_reports),
    )


def set_compliance_config(level: ComplianceLevel) -> ComplianceConfig:
    """
    Set the global compliance configuration level.

    Args:
        level: The compliance level to set.

    Returns:
        ComplianceConfig: The new compliance configuration.

    Example:
        >>> set_compliance_config(ComplianceLevel.DISABLED)
        >>> # Fast development mode without compliance overhead
    """
    jl = get_julia()

    # Map Python enum to Julia enum
    julia_levels = {
        ComplianceLevel.DISABLED: jl.NeoPKPD.COMPLIANCE_DISABLED,
        ComplianceLevel.MINIMAL: jl.NeoPKPD.COMPLIANCE_MINIMAL,
        ComplianceLevel.STANDARD: jl.NeoPKPD.COMPLIANCE_STANDARD,
        ComplianceLevel.STRICT: jl.NeoPKPD.COMPLIANCE_STRICT,
    }

    julia_level = julia_levels[level]
    jl.NeoPKPD.set_compliance_config_b(julia_level)

    return get_compliance_config()
