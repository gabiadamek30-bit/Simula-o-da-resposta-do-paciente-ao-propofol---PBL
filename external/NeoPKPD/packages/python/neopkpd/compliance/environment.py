"""
Environment capture for FDA 21 CFR Part 11 compliance.
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional

from ..bridge import get_julia


@dataclass(frozen=True)
class EnvironmentSnapshot:
    """
    Comprehensive snapshot of the execution environment.

    Attributes:
        julia_version: Julia version (e.g., "1.10.0")
        julia_build: Julia build information
        neopkpd_version: NeoPKPD Core version
        artifact_schema_version: Artifact schema version
        event_semantics_version: Event semantics version
        solver_semantics_version: Solver semantics version
        differentialequations_version: DifferentialEquations.jl version
        scimlbase_version: SciMLBase.jl version
        os: Operating system name
        os_version: Operating system version
        cpu_model: CPU model name
        cpu_cores: Number of CPU cores
        memory_gb: Total system memory in GB
        architecture: System architecture
        word_size: Julia word size (32 or 64)
        python_version: Python version if available
        git_commit: Git commit hash if in repository
        git_branch: Git branch name if in repository
        git_dirty: Whether working directory has uncommitted changes
        capture_timestamp: ISO 8601 UTC timestamp
    """

    julia_version: str
    julia_build: str
    neopkpd_version: str
    artifact_schema_version: str
    event_semantics_version: str
    solver_semantics_version: str
    differentialequations_version: str
    scimlbase_version: str
    os: str
    os_version: str
    cpu_model: str
    cpu_cores: int
    memory_gb: float
    architecture: str
    word_size: int
    python_version: Optional[str]
    git_commit: Optional[str]
    git_branch: Optional[str]
    git_dirty: Optional[bool]
    capture_timestamp: str


def get_environment_snapshot() -> EnvironmentSnapshot:
    """
    Capture a comprehensive snapshot of the current execution environment.

    This includes Julia version, NeoPKPD version, OS information, CPU details,
    and git repository status.

    Returns:
        EnvironmentSnapshot: Complete environment information.

    Example:
        >>> env = get_environment_snapshot()
        >>> print(f"Julia: {env.julia_version}")
        Julia: 1.10.0
        >>> print(f"NeoPKPD: {env.neopkpd_version}")
        NeoPKPD: 0.1.0
    """
    jl = get_julia()
    env = jl.NeoPKPD.capture_environment()

    return EnvironmentSnapshot(
        julia_version=str(env.julia_version),
        julia_build=str(env.julia_build),
        neopkpd_version=str(env.neopkpd_version),
        artifact_schema_version=str(env.artifact_schema_version),
        event_semantics_version=str(env.event_semantics_version),
        solver_semantics_version=str(env.solver_semantics_version),
        differentialequations_version=str(env.differentialequations_version),
        scimlbase_version=str(env.scimlbase_version),
        os=str(env.os),
        os_version=str(env.os_version),
        cpu_model=str(env.cpu_model),
        cpu_cores=int(env.cpu_cores),
        memory_gb=float(env.memory_gb),
        architecture=str(env.architecture),
        word_size=int(env.word_size),
        python_version=str(env.python_version) if env.python_version is not None else None,
        git_commit=str(env.git_commit) if env.git_commit is not None else None,
        git_branch=str(env.git_branch) if env.git_branch is not None else None,
        git_dirty=bool(env.git_dirty) if env.git_dirty is not None else None,
        capture_timestamp=str(env.capture_timestamp),
    )


def get_artifact_environment(artifact: Dict[str, Any]) -> Optional[EnvironmentSnapshot]:
    """
    Extract the environment snapshot from an artifact.

    Args:
        artifact: The artifact dictionary.

    Returns:
        EnvironmentSnapshot if present, None otherwise.

    Example:
        >>> with open("simulation.json") as f:
        ...     artifact = json.load(f)
        >>> env = get_artifact_environment(artifact)
        >>> if env:
        ...     print(f"Created with Julia {env.julia_version}")
    """
    if "compliance_metadata" not in artifact:
        return None

    meta = artifact["compliance_metadata"]
    if "environment" not in meta:
        return None

    env = meta["environment"]

    return EnvironmentSnapshot(
        julia_version=env["julia_version"],
        julia_build=env["julia_build"],
        neopkpd_version=env["neopkpd_version"],
        artifact_schema_version=env["artifact_schema_version"],
        event_semantics_version=env["event_semantics_version"],
        solver_semantics_version=env["solver_semantics_version"],
        differentialequations_version=env["differentialequations_version"],
        scimlbase_version=env["scimlbase_version"],
        os=env["os"],
        os_version=env["os_version"],
        cpu_model=env["cpu_model"],
        cpu_cores=env["cpu_cores"],
        memory_gb=env["memory_gb"],
        architecture=env["architecture"],
        word_size=env["word_size"],
        python_version=env.get("python_version"),
        git_commit=env.get("git_commit"),
        git_branch=env.get("git_branch"),
        git_dirty=env.get("git_dirty"),
        capture_timestamp=env["capture_timestamp"],
    )


def compare_environments(
    artifact1: Dict[str, Any], artifact2: Dict[str, Any]
) -> Dict[str, tuple]:
    """
    Compare the environments of two artifacts.

    Args:
        artifact1: First artifact dictionary.
        artifact2: Second artifact dictionary.

    Returns:
        Dictionary of fields that differ, with (artifact1_value, artifact2_value) tuples.

    Example:
        >>> with open("old.json") as f:
        ...     old = json.load(f)
        >>> with open("new.json") as f:
        ...     new = json.load(f)
        >>> diffs = compare_environments(old, new)
        >>> for field, (v1, v2) in diffs.items():
        ...     print(f"{field}: {v1} -> {v2}")
    """
    env1 = get_artifact_environment(artifact1)
    env2 = get_artifact_environment(artifact2)

    if env1 is None or env2 is None:
        return {}

    differences = {}
    for field in env1.__dataclass_fields__:
        v1 = getattr(env1, field)
        v2 = getattr(env2, field)
        if v1 != v2:
            differences[field] = (v1, v2)

    return differences
