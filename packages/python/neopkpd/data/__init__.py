"""
NeoPKPD CDISC/SDTM Data Module

This module provides Python bindings for reading and converting
CDISC/SDTM clinical data formats (PC, EX, DM domains).
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Union
from pathlib import Path

from .._core import _require_julia


@dataclass
class SubjectData:
    """Data for a single subject."""
    subject_id: str
    times: List[float]
    observations: List[float]
    doses: List[Dict[str, float]]
    covariates: Optional[Dict[str, float]] = None


@dataclass
class ObservedData:
    """Container for observed clinical data."""
    subjects: List[SubjectData]
    n_subjects: int
    n_observations: int


def read_cdisc_pc(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Read CDISC PC (Pharmacokinetic Concentrations) domain.

    The PC domain contains concentration-time data for each subject.

    Args:
        path: Path to PC domain CSV file

    Returns:
        Dict containing parsed PC records

    Example:
        >>> pc_data = read_cdisc_pc("pc.csv")
        >>> print(f"Found {len(pc_data['records'])} concentration records")
    """
    jl = _require_julia()
    result = jl.NeoPKPD.read_cdisc_pc(str(Path(path).resolve()))
    return _convert_cdisc_domain(result)


def read_cdisc_ex(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Read CDISC EX (Exposure) domain.

    The EX domain contains dosing information for each subject.

    Args:
        path: Path to EX domain CSV file

    Returns:
        Dict containing parsed EX records

    Example:
        >>> ex_data = read_cdisc_ex("ex.csv")
        >>> print(f"Found {len(ex_data['records'])} dosing records")
    """
    jl = _require_julia()
    result = jl.NeoPKPD.read_cdisc_ex(str(Path(path).resolve()))
    return _convert_cdisc_domain(result)


def read_cdisc_dm(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Read CDISC DM (Demographics) domain.

    The DM domain contains demographic and covariate data for each subject.

    Args:
        path: Path to DM domain CSV file

    Returns:
        Dict containing parsed DM records

    Example:
        >>> dm_data = read_cdisc_dm("dm.csv")
        >>> print(f"Found {len(dm_data['records'])} subject demographics")
    """
    jl = _require_julia()
    result = jl.NeoPKPD.read_cdisc_dm(str(Path(path).resolve()))
    return _convert_cdisc_domain(result)


def cdisc_to_observed_data(
    pc_data: Dict[str, Any],
    ex_data: Optional[Dict[str, Any]] = None,
    dm_data: Optional[Dict[str, Any]] = None,
) -> ObservedData:
    """
    Convert CDISC domain data to NeoPKPD observed data format.

    Args:
        pc_data: Parsed PC domain from read_cdisc_pc
        ex_data: Optional parsed EX domain from read_cdisc_ex
        dm_data: Optional parsed DM domain from read_cdisc_dm

    Returns:
        ObservedData suitable for estimation or VPC

    Example:
        >>> pc = read_cdisc_pc("pc.csv")
        >>> ex = read_cdisc_ex("ex.csv")
        >>> observed = cdisc_to_observed_data(pc, ex)
        >>> print(f"Data for {observed.n_subjects} subjects")
    """
    jl = _require_julia()

    # Convert back to Julia format for the conversion function
    pc_jl = _convert_to_julia_domain(jl, pc_data, "PC")
    ex_jl = _convert_to_julia_domain(jl, ex_data, "EX") if ex_data else None
    dm_jl = _convert_to_julia_domain(jl, dm_data, "DM") if dm_data else None

    result = jl.NeoPKPD.cdisc_to_observed_data(pc_jl, ex_jl, dm_jl)
    return _convert_observed_data(result)


def read_cdisc_dataset(
    pc_path: Union[str, Path],
    ex_path: Optional[Union[str, Path]] = None,
    dm_path: Optional[Union[str, Path]] = None,
) -> ObservedData:
    """
    Read complete CDISC dataset and convert to NeoPKPD format.

    Convenience function that combines read_cdisc_* and cdisc_to_observed_data.

    Args:
        pc_path: Path to PC domain CSV file
        ex_path: Optional path to EX domain CSV file
        dm_path: Optional path to DM domain CSV file

    Returns:
        ObservedData ready for analysis

    Example:
        >>> observed = read_cdisc_dataset("pc.csv", "ex.csv", "dm.csv")
        >>> for subj in observed.subjects:
        ...     print(f"Subject {subj.subject_id}: {len(subj.observations)} obs")
    """
    pc = read_cdisc_pc(pc_path)
    ex = read_cdisc_ex(ex_path) if ex_path else None
    dm = read_cdisc_dm(dm_path) if dm_path else None
    return cdisc_to_observed_data(pc, ex, dm)


def validate_cdisc_dataset(
    pc_path: Union[str, Path],
    ex_path: Optional[Union[str, Path]] = None,
    dm_path: Optional[Union[str, Path]] = None,
) -> List[str]:
    """
    Validate CDISC dataset for common issues.

    Checks for:
    - Missing required columns
    - Invalid data types
    - Inconsistent subject IDs across domains
    - Missing dosing information
    - Negative concentrations

    Args:
        pc_path: Path to PC domain CSV
        ex_path: Optional path to EX domain CSV
        dm_path: Optional path to DM domain CSV

    Returns:
        List of warning/error messages (empty if valid)

    Example:
        >>> warnings = validate_cdisc_dataset("pc.csv", "ex.csv")
        >>> if warnings:
        ...     for w in warnings:
        ...         print(f"Warning: {w}")
    """
    jl = _require_julia()
    result = jl.NeoPKPD.validate_cdisc_dataset(
        str(Path(pc_path).resolve()),
        str(Path(ex_path).resolve()) if ex_path else None,
        str(Path(dm_path).resolve()) if dm_path else None,
    )
    return [str(w) for w in result]


def _convert_cdisc_domain(domain) -> Dict[str, Any]:
    """Convert Julia CDISC domain to Python dict."""
    records = []
    for rec in domain.records:
        record = {}
        for field in dir(rec):
            if not field.startswith("_"):
                val = getattr(rec, field)
                if val is not None:
                    if isinstance(val, (int, float, str, bool)):
                        record[field] = val
                    else:
                        record[field] = str(val)
        records.append(record)

    return {
        "domain": str(domain.domain) if hasattr(domain, "domain") else "UNKNOWN",
        "records": records,
        "n_records": len(records),
    }


def _convert_to_julia_domain(jl, data: Dict[str, Any], domain_type: str):
    """Convert Python dict back to Julia domain for processing."""
    # This is a simplified version - actual implementation would need
    # to properly reconstruct Julia types
    return data


def _convert_observed_data(result) -> ObservedData:
    """Convert Julia ObservedData to Python dataclass."""
    subjects = []
    for subj in result.subjects:
        doses = [{"time": float(d.time), "amount": float(d.amount)} for d in subj.doses]
        subjects.append(SubjectData(
            subject_id=str(subj.subject_id),
            times=list(subj.times),
            observations=list(subj.observations),
            doses=doses,
            covariates=None,  # Would extract if present
        ))

    return ObservedData(
        subjects=subjects,
        n_subjects=len(subjects),
        n_observations=sum(len(s.observations) for s in subjects),
    )


__all__ = [
    "read_cdisc_pc",
    "read_cdisc_ex",
    "read_cdisc_dm",
    "cdisc_to_observed_data",
    "read_cdisc_dataset",
    "validate_cdisc_dataset",
    "SubjectData",
    "ObservedData",
]
