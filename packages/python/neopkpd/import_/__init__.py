"""
NeoPKPD Model Import Module

This module provides Python bindings for importing models from
NONMEM control files (.ctl) and Monolix project files (.mlxtran).

Features:
- Parse and convert NONMEM control files (ADVAN1-4, 10, 11)
- Parse and convert Monolix project files (.mlxtran)
- Extract parameter values, bounds, and IIV structure
- Detect and report unsupported constructs
- Support for covariates and error model extraction

Example:
    >>> import neopkpd
    >>> from neopkpd.import_ import import_nonmem, import_monolix
    >>> neopkpd.init_julia()
    >>>
    >>> # Import NONMEM model
    >>> model = import_nonmem("run001.ctl", doses=[{"time": 0.0, "amount": 100.0}])
    >>> print(f"Model: {model.model_kind}")
    >>> print(f"Parameters: {model.params}")
    >>>
    >>> # Import Monolix model
    >>> model = import_monolix("project.mlxtran", doses=[{"time": 0.0, "amount": 100.0}])
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union
from pathlib import Path

from .._core import _require_julia


# ============================================================================
# Enumerations
# ============================================================================

class ErrorModelType(str, Enum):
    """Types of residual error models."""
    PROPORTIONAL = "proportional"
    ADDITIVE = "additive"
    COMBINED = "combined"
    EXPONENTIAL = "exponential"
    UNKNOWN = "unknown"


class IIVTransformation(str, Enum):
    """Transformation types for inter-individual variability."""
    EXPONENTIAL = "exponential"
    ADDITIVE = "additive"
    PROPORTIONAL = "proportional"
    NONE = "none"


class CovariateEffectType(str, Enum):
    """Types of covariate effects."""
    POWER = "power"
    LINEAR = "linear"
    EXPONENTIAL = "exponential"


class OmegaStructure(str, Enum):
    """Structure of omega matrices."""
    DIAGONAL = "diagonal"
    BLOCK = "block"
    SAME = "same"


# ============================================================================
# NONMEM Types
# ============================================================================

@dataclass
class THETASpec:
    """
    Specification for a single THETA parameter.

    Attributes:
        init: Initial estimate
        lower: Lower bound (can be -inf)
        upper: Upper bound (can be inf)
        fixed: Whether parameter is fixed
        name: Optional parameter name/label
    """
    init: float
    lower: float = float('-inf')
    upper: float = float('inf')
    fixed: bool = False
    name: str = ""


@dataclass
class OMEGABlock:
    """
    OMEGA block specification for inter-individual variability.

    Attributes:
        values: Variance/covariance values (diagonal or block)
        structure: diagonal, block, or same
        dimension: Size of the block
        fixed: Whether block is fixed
    """
    values: List[float]
    structure: OmegaStructure = OmegaStructure.DIAGONAL
    dimension: int = 0
    fixed: bool = False

    def __post_init__(self):
        if self.dimension == 0:
            if self.structure == OmegaStructure.DIAGONAL:
                self.dimension = len(self.values)
            else:
                # For block: n*(n+1)/2 = len, solve for n
                import math
                n = int((-1 + math.sqrt(1 + 8 * len(self.values))) / 2)
                self.dimension = n


@dataclass
class SIGMABlock:
    """
    SIGMA block specification for residual error.

    Attributes:
        values: Variance values
        structure: diagonal or block
        fixed: Whether fixed
    """
    values: List[float]
    structure: OmegaStructure = OmegaStructure.DIAGONAL
    fixed: bool = False


@dataclass
class SubroutineSpec:
    """
    $SUBROUTINES specification (ADVAN/TRANS).

    Attributes:
        advan: ADVAN number (1-13)
        trans: TRANS number (1-6)
        other: Other subroutines specified
    """
    advan: int
    trans: int = 1
    other: List[str] = field(default_factory=list)


@dataclass
class DataSpec:
    """
    $DATA specification.

    Attributes:
        filename: Path to data file
        ignore: Characters/conditions to ignore
        accept: Conditions to accept
    """
    filename: str
    ignore: List[str] = field(default_factory=list)
    accept: List[str] = field(default_factory=list)


@dataclass
class InputColumn:
    """
    $INPUT column specification.

    Attributes:
        name: Column name
        drop: Whether to drop this column
        alias: Alias for standard names (ID, TIME, DV, etc.)
    """
    name: str
    drop: bool = False
    alias: str = ""


@dataclass
class PKCovariateEffect:
    """
    Covariate effect extracted from $PK block.

    Represents relationships like:
    - Power: TVCL * (WT/70)**THETA(3)
    - Linear: TVCL * (1 + THETA(3)*(AGE-40))
    - Exponential: TVCL * EXP(THETA(3)*(CRCL-100))

    Attributes:
        covariate: Covariate name (WT, AGE, CRCL, etc.)
        theta_index: THETA index for the coefficient
        effect_type: power, linear, or exponential
        reference: Reference value (70 for WT, 40 for AGE, etc.)
    """
    covariate: str
    theta_index: int
    effect_type: CovariateEffectType
    reference: float = 0.0


@dataclass
class PKAssignment:
    """
    Parameter assignment extracted from $PK block.

    Represents assignments like:
    - TVCL = THETA(1) (typical value definition)
    - CL = TVCL * EXP(ETA(1)) (individual parameter with IIV)

    Attributes:
        target: Parameter name (CL, V, KA, etc.)
        tv_theta: THETA index for typical value
        eta_index: ETA index for IIV
        transformation: exponential, additive, or none
        covariate_effects: Vector of covariate effects
    """
    target: str
    tv_theta: Optional[int] = None
    eta_index: Optional[int] = None
    transformation: IIVTransformation = IIVTransformation.NONE
    covariate_effects: List[PKCovariateEffect] = field(default_factory=list)


@dataclass
class ScalingFactor:
    """
    Scaling factor from $PK block (S1 = V, S2 = V1, etc.).

    Attributes:
        compartment: Compartment number
        parameter: Parameter symbol it maps to
    """
    compartment: int
    parameter: str


@dataclass
class PKBlock:
    """
    Parsed $PK block.

    Attributes:
        tv_definitions: Map from TV name to THETA index
        assignments: Parameter assignments with ETA and covariate info
        scaling: Scaling factors (S1, S2, etc.)
        raw_code: Original lines from $PK block
        unsupported_lines: Lines containing unsupported constructs
    """
    tv_definitions: Dict[str, int] = field(default_factory=dict)
    assignments: List[PKAssignment] = field(default_factory=list)
    scaling: List[ScalingFactor] = field(default_factory=list)
    raw_code: List[str] = field(default_factory=list)
    unsupported_lines: List[str] = field(default_factory=list)


@dataclass
class ErrorBlock:
    """
    Parsed $ERROR block.

    Attributes:
        error_type: Detected error model type
        theta_indices: THETA indices used in error model
        sigma_fixed_to_1: Whether SIGMA is fixed to 1
        raw_code: Original lines
        unsupported_lines: Lines with unsupported constructs
    """
    error_type: ErrorModelType = ErrorModelType.UNKNOWN
    theta_indices: List[int] = field(default_factory=list)
    sigma_fixed_to_1: bool = False
    raw_code: List[str] = field(default_factory=list)
    unsupported_lines: List[str] = field(default_factory=list)


@dataclass
class UnsupportedConstruct:
    """
    Represents an unsupported NONMEM construct detected during parsing.

    Attributes:
        construct: Name of the unsupported construct
        location: Where it was found ($PK, $ERROR, etc.)
        line: The offending line of code
        message: Human-readable error message
    """
    construct: str
    location: str
    line: str
    message: str = ""

    def __post_init__(self):
        if not self.message:
            self.message = f"{self.construct} in {self.location} is not supported: {self.line}"


@dataclass
class NONMEMControlFile:
    """
    Complete NONMEM control file representation.

    Attributes:
        problem: Problem description from $PROBLEM
        data: Data file specification
        input_columns: Input column definitions
        subroutines: ADVAN/TRANS specification
        thetas: Vector of THETA specifications
        omegas: Vector of OMEGA blocks
        sigmas: Vector of SIGMA blocks
        pk_block: Parsed $PK block
        error_block: Parsed $ERROR block
        estimation: Estimation method settings
        tables: Table output specifications
        raw_text: Original control file text
    """
    problem: str
    data: Optional[DataSpec] = None
    input_columns: List[InputColumn] = field(default_factory=list)
    subroutines: Optional[SubroutineSpec] = None
    thetas: List[THETASpec] = field(default_factory=list)
    omegas: List[OMEGABlock] = field(default_factory=list)
    sigmas: List[SIGMABlock] = field(default_factory=list)
    pk_block: Optional[PKBlock] = None
    error_block: Optional[ErrorBlock] = None
    estimation: Dict[str, Any] = field(default_factory=dict)
    tables: List[Dict[str, Any]] = field(default_factory=list)
    raw_text: str = ""

    # Legacy compatibility properties
    @property
    def data_file(self) -> Optional[str]:
        return self.data.filename if self.data else None

    @property
    def advan(self) -> int:
        return self.subroutines.advan if self.subroutines else 0

    @property
    def trans(self) -> int:
        return self.subroutines.trans if self.subroutines else 1

    @property
    def theta_specs(self) -> List[Dict[str, Any]]:
        return [
            {"init": t.init, "lower": t.lower, "upper": t.upper, "fixed": t.fixed}
            for t in self.thetas
        ]

    @property
    def omega_block(self) -> List[List[float]]:
        return [[v] for o in self.omegas for v in o.values] if self.omegas else []

    @property
    def sigma_block(self) -> List[List[float]]:
        return [[v] for s in self.sigmas for v in s.values] if self.sigmas else []

    @property
    def pk_code(self) -> List[str]:
        return self.pk_block.raw_code if self.pk_block else []

    @property
    def error_code(self) -> List[str]:
        return self.error_block.raw_code if self.error_block else []


# ============================================================================
# Monolix Types
# ============================================================================

@dataclass
class MonolixModelType:
    """
    Monolix model type specification.

    Attributes:
        lib: Library name (e.g., "pklib")
        model: Model name (e.g., "pk_oral1cpt_kaVCl_PLASMA")
    """
    lib: str
    model: str


@dataclass
class MonolixStructuralModel:
    """
    Monolix structural model configuration.

    Attributes:
        model_type: The model library and name
        admin_type: Administration type (oral, iv, infusion)
        n_compartments: Number of compartments
        elimination: Elimination type (linear, mm, mixed)
        absorption: Absorption type (firstOrder, zeroOrder, etc.)
        has_lag: Whether model has absorption lag
        has_bioavailability: Whether model has bioavailability parameter
    """
    model_type: Optional[MonolixModelType] = None
    admin_type: str = "iv"
    n_compartments: int = 1
    elimination: str = "linear"
    absorption: str = "bolus"
    has_lag: bool = False
    has_bioavailability: bool = False


@dataclass
class MonolixParameter:
    """
    Monolix parameter definition.

    Attributes:
        name: Parameter name (e.g., "ka", "V", "Cl")
        value: Initial value
        fixed: Whether parameter is fixed
        distribution: Distribution type (logNormal, normal, logitNormal)
        omega: IIV variance (if has_iiv)
        has_iiv: Whether parameter has inter-individual variability
    """
    name: str
    value: float
    fixed: bool = False
    distribution: str = "logNormal"
    omega: float = 0.0
    has_iiv: bool = False


@dataclass
class MonolixObservation:
    """
    Monolix observation definition.

    Attributes:
        name: Observation name (e.g., "y1", "Cc")
        type: Observation type (continuous, discrete)
        error_model: Error model type
        error_params: Error model parameters
    """
    name: str
    type: str = "continuous"
    error_model: str = "combined"
    error_params: List[float] = field(default_factory=list)


@dataclass
class MonolixDataset:
    """
    Monolix dataset specification.

    Attributes:
        filename: Path to data file
        header_types: Column type mappings
        id_column: Subject ID column
        time_column: Time column
        observation_column: Observation column
        dose_column: Dose column (if present)
        rate_column: Infusion rate column (if present)
    """
    filename: str
    header_types: Dict[str, str] = field(default_factory=dict)
    id_column: str = "ID"
    time_column: str = "TIME"
    observation_column: str = "DV"
    dose_column: str = "AMT"
    rate_column: str = "RATE"


@dataclass
class MonolixProject:
    """
    Complete Monolix project representation.

    Attributes:
        description: Project description
        data: Dataset specification
        model: Structural model configuration
        parameters: Vector of parameter definitions
        observations: Vector of observation definitions
        estimation_method: Estimation algorithm used
        raw_text: Original .mlxtran file content
    """
    description: str = ""
    data: Optional[MonolixDataset] = None
    model: Optional[MonolixStructuralModel] = None
    parameters: List[MonolixParameter] = field(default_factory=list)
    observations: List[MonolixObservation] = field(default_factory=list)
    estimation_method: str = "SAEM"
    raw_text: str = ""

    # Legacy compatibility properties
    @property
    def model_type(self) -> str:
        if self.model and self.model.model_type:
            return f"{self.model.model_type.lib}:{self.model.model_type.model}"
        return "unknown"

    @property
    def structural_model(self) -> str:
        if self.model and self.model.model_type:
            return self.model.model_type.model
        return ""

    @property
    def data_file(self) -> Optional[str]:
        return self.data.filename if self.data else None


@dataclass
class UnsupportedMonolixConstruct:
    """
    Represents an unsupported Monolix construct.

    Attributes:
        construct: Name of the unsupported construct
        location: Where it was found
        line: The offending line of code
        message: Human-readable error message
    """
    construct: str
    location: str
    line: str
    message: str = ""

    def __post_init__(self):
        if not self.message:
            self.message = f"{self.construct} in {self.location} is not supported: {self.line}"


# ============================================================================
# Import Result
# ============================================================================

@dataclass
class ImportedModel:
    """
    Result from importing an external model file.

    Attributes:
        source_format: "nonmem" or "monolix"
        source_file: Path to source file
        model_kind: NeoPKPD model type name
        params: Dictionary of parameter name => value
        theta_init: Initial theta values
        theta_names: Parameter names
        omega_init: Initial omega values (list of lists)
        omega_names: Names of random effects
        sigma_type: Residual error model type
        sigma_init: Initial sigma value
        warnings: Any warnings generated during import
        metadata: Additional metadata
    """
    source_format: str
    source_file: str
    model_kind: str
    params: Dict[str, float]
    theta_init: List[float]
    theta_names: List[str]
    omega_init: List[List[float]]
    omega_names: List[str]
    sigma_type: str
    sigma_init: float
    warnings: List[str]
    metadata: Dict[str, Any]


# ============================================================================
# ADVAN/TRANS Mapping
# ============================================================================

ADVAN_TRANS_MAP: Dict[Tuple[int, int], Tuple[str, List[str]]] = {
    # ADVAN1: One-compartment IV bolus
    (1, 1): ("OneCompIVBolus", ["K", "V"]),
    (1, 2): ("OneCompIVBolus", ["CL", "V"]),

    # ADVAN2: One-compartment oral first-order
    (2, 1): ("OneCompOralFirstOrder", ["KA", "K", "V"]),
    (2, 2): ("OneCompOralFirstOrder", ["KA", "CL", "V"]),

    # ADVAN3: Two-compartment IV bolus
    (3, 1): ("TwoCompIVBolus", ["K", "K12", "K21", "V"]),
    (3, 3): ("TwoCompIVBolus", ["CL", "V", "Q", "VSS"]),
    (3, 4): ("TwoCompIVBolus", ["CL", "V1", "Q", "V2"]),

    # ADVAN4: Two-compartment oral first-order
    (4, 1): ("TwoCompOral", ["KA", "K", "K23", "K32", "V"]),
    (4, 3): ("TwoCompOral", ["KA", "CL", "V", "Q", "VSS"]),
    (4, 4): ("TwoCompOral", ["KA", "CL", "V1", "Q", "V2"]),

    # ADVAN11: Three-compartment IV bolus
    (11, 1): ("ThreeCompIVBolus", ["K", "K12", "K21", "K13", "K31", "V"]),
    (11, 4): ("ThreeCompIVBolus", ["CL", "V1", "Q2", "V2", "Q3", "V3"]),

    # ADVAN10: Michaelis-Menten elimination
    (10, 1): ("MichaelisMentenElimination", ["VM", "KM", "V"]),
}


def get_model_mapping(advan: int, trans: int) -> Optional[Tuple[str, List[str]]]:
    """
    Get the NeoPKPD model kind and expected parameters for ADVAN/TRANS.

    Args:
        advan: ADVAN number
        trans: TRANS number

    Returns:
        Tuple of (model_kind, parameter_names) or None if not supported
    """
    return ADVAN_TRANS_MAP.get((advan, trans))


# ============================================================================
# Monolix Model Mapping
# ============================================================================

MONOLIX_MODEL_MAP: Dict[str, str] = {
    # One-compartment models
    "pk_bolus1cpt_Vk_PLASMA": "OneCompIVBolus",
    "pk_bolus1cpt_VCl_PLASMA": "OneCompIVBolus",
    "pk_oral1cpt_1abs_kaTbioVk_PLASMA": "OneCompOralFirstOrder",
    "pk_oral1cpt_1abs_kaVCl_PLASMA": "OneCompOralFirstOrder",
    "pk_oral1cpt_1abs_TlagkaVCl_PLASMA": "OneCompOralFirstOrder",

    # Two-compartment models
    "pk_bolus2cpt_Vk12k21k_PLASMA": "TwoCompIVBolus",
    "pk_bolus2cpt_V1ClQ2V2_PLASMA": "TwoCompIVBolus",
    "pk_oral2cpt_1abs_kaV1ClQ2V2_PLASMA": "TwoCompOral",
    "pk_oral2cpt_1abs_TlagkaV1ClQ2V2_PLASMA": "TwoCompOral",

    # Three-compartment models
    "pk_bolus3cpt_V1ClQ2V2Q3V3_PLASMA": "ThreeCompIVBolus",

    # Michaelis-Menten models
    "pk_bolus1cpt_VVmKm_PLASMA": "MichaelisMentenElimination",
}


def get_monolix_model_mapping(model_name: str) -> Optional[str]:
    """
    Get the NeoPKPD model kind for a Monolix model name.

    Args:
        model_name: Monolix model name

    Returns:
        NeoPKPD model kind or None if not supported
    """
    # Try exact match first
    if model_name in MONOLIX_MODEL_MAP:
        return MONOLIX_MODEL_MAP[model_name]

    # Try pattern matching
    model_lower = model_name.lower()

    if "oral" in model_lower and "1cpt" in model_lower:
        return "OneCompOralFirstOrder"
    elif "bolus" in model_lower and "1cpt" in model_lower:
        return "OneCompIVBolus"
    elif "oral" in model_lower and "2cpt" in model_lower:
        return "TwoCompOral"
    elif "bolus" in model_lower and "2cpt" in model_lower:
        return "TwoCompIVBolus"
    elif "3cpt" in model_lower:
        return "ThreeCompIVBolus"
    elif "mm" in model_lower or "michaelis" in model_lower:
        return "MichaelisMentenElimination"

    return None


# ============================================================================
# Import Functions
# ============================================================================

def import_nonmem(
    path: Union[str, Path],
    doses: Optional[List[Dict[str, float]]] = None
) -> ImportedModel:
    """
    Import a model from a NONMEM control file.

    Supports ADVAN1-4, ADVAN10, ADVAN11 with common TRANS options. Extracts:
    - Model structure (ADVAN/TRANS)
    - THETA initial estimates and bounds
    - OMEGA structure and values
    - SIGMA values
    - Covariate relationships
    - Error model specification

    Args:
        path: Path to NONMEM control file (.ctl or .mod)
        doses: Optional list of dose events, e.g., [{"time": 0.0, "amount": 100.0}]

    Returns:
        ImportedModel with NeoPKPD-compatible model specification

    Supported ADVAN/TRANS combinations:
        - ADVAN1/TRANS2: One-compartment IV bolus (CL, V)
        - ADVAN2/TRANS2: One-compartment oral (Ka, CL, V)
        - ADVAN3/TRANS4: Two-compartment IV bolus (CL, V1, Q, V2)
        - ADVAN4/TRANS4: Two-compartment oral (Ka, CL, V1, Q, V2)
        - ADVAN10: Michaelis-Menten elimination
        - ADVAN11: Three-compartment model

    Example:
        >>> model = import_nonmem("run001.ctl", doses=[{"time": 0.0, "amount": 100.0}])
        >>> print(f"Model type: {model.model_kind}")
        >>> print(f"Parameters: {model.params}")
    """
    jl = _require_julia()

    # Read the control file text
    path = Path(path).resolve()
    with open(path, 'r') as f:
        ctl_text = f.read()

    # Parse the control file
    ctl = jl.NeoPKPD.parse_nonmem_control(ctl_text)

    # Create dose events if provided
    jl_doses = jl.seval("DoseEvent[]")
    if doses:
        for d in doses:
            dose = jl.NeoPKPD.DoseEvent(
                float(d.get("time", 0.0)),
                float(d.get("amount", 100.0))
            )
            jl.seval("push!")(jl_doses, dose)

    # Convert to NeoPKPD format
    result = jl.NeoPKPD.convert_nonmem_to_neopkpd(ctl, doses=jl_doses)

    return _convert_julia_import_result(result, "nonmem", str(path))


def import_monolix(
    path: Union[str, Path],
    doses: Optional[List[Dict[str, float]]] = None
) -> ImportedModel:
    """
    Import a model from a Monolix project file.

    Parses .mlxtran files and extracts:
    - Structural model definition
    - Parameter initial values
    - Random effect structure
    - Error model specification

    Args:
        path: Path to Monolix project file (.mlxtran)
        doses: Optional list of dose events

    Returns:
        ImportedModel with NeoPKPD-compatible model specification

    Example:
        >>> model = import_monolix("project.mlxtran", doses=[{"time": 0.0, "amount": 100.0}])
        >>> print(f"Model type: {model.model_kind}")
    """
    jl = _require_julia()

    # Read the project file text
    path = Path(path).resolve()
    with open(path, 'r') as f:
        mlx_text = f.read()

    # Parse the project file
    project = jl.NeoPKPD.parse_monolix_project(mlx_text)

    # Create dose events if provided
    jl_doses = jl.seval("DoseEvent[]")
    if doses:
        for d in doses:
            dose = jl.NeoPKPD.DoseEvent(
                float(d.get("time", 0.0)),
                float(d.get("amount", 100.0))
            )
            jl.seval("push!")(jl_doses, dose)

    # Convert to NeoPKPD format
    result = jl.NeoPKPD.convert_monolix_to_neopkpd(project, doses=jl_doses)

    return _convert_julia_import_result(result, "monolix", str(path))


def import_model(
    path: Union[str, Path],
    format: Optional[str] = None,
    doses: Optional[List[Dict[str, float]]] = None
) -> ImportedModel:
    """
    Import a model from NONMEM or Monolix format.

    Auto-detects format based on file extension if not specified.

    Args:
        path: Path to model file
        format: Optional format override ("nonmem" or "monolix")
        doses: Optional list of dose events

    Returns:
        ImportedModel with NeoPKPD-compatible model specification

    Example:
        >>> model = import_model("run001.ctl")
        >>> model = import_model("project.mlxtran")
    """
    path = Path(path)

    if format is None:
        if path.suffix in [".ctl", ".mod"]:
            format = "nonmem"
        elif path.suffix == ".mlxtran":
            format = "monolix"
        else:
            raise ValueError(
                f"Cannot auto-detect format for {path.suffix}. "
                "Specify format='nonmem' or format='monolix'"
            )

    if format == "nonmem":
        return import_nonmem(path, doses=doses)
    elif format == "monolix":
        return import_monolix(path, doses=doses)
    else:
        raise ValueError(f"Unknown format: {format}")


# ============================================================================
# Parse Functions
# ============================================================================

def parse_nonmem_control(path: Union[str, Path]) -> NONMEMControlFile:
    """
    Parse a NONMEM control file without converting to NeoPKPD format.

    Useful for inspecting the control file structure.

    Args:
        path: Path to NONMEM control file

    Returns:
        NONMEMControlFile with parsed sections

    Example:
        >>> ctl = parse_nonmem_control("run001.ctl")
        >>> print(f"Problem: {ctl.problem}")
        >>> print(f"ADVAN{ctl.advan} TRANS{ctl.trans}")
        >>> print(f"Number of THETAs: {len(ctl.thetas)}")
    """
    jl = _require_julia()

    path = Path(path).resolve()
    with open(path, 'r') as f:
        ctl_text = f.read()

    result = jl.NeoPKPD.parse_nonmem_control(ctl_text)
    return _convert_nonmem_control(result, ctl_text)


def parse_nonmem_control_text(text: str) -> NONMEMControlFile:
    """
    Parse NONMEM control file text directly.

    Args:
        text: Raw control file text

    Returns:
        NONMEMControlFile with parsed sections
    """
    jl = _require_julia()
    result = jl.NeoPKPD.parse_nonmem_control(text)
    return _convert_nonmem_control(result, text)


def parse_monolix_project(path: Union[str, Path]) -> MonolixProject:
    """
    Parse a Monolix project file without converting to NeoPKPD format.

    Args:
        path: Path to Monolix project file

    Returns:
        MonolixProject with parsed structure

    Example:
        >>> mlx = parse_monolix_project("project.mlxtran")
        >>> print(f"Model: {mlx.structural_model}")
        >>> print(f"Parameters: {[p.name for p in mlx.parameters]}")
    """
    jl = _require_julia()

    path = Path(path).resolve()
    with open(path, 'r') as f:
        mlx_text = f.read()

    result = jl.NeoPKPD.parse_monolix_project(mlx_text)
    return _convert_monolix_project(result, mlx_text)


def parse_monolix_project_text(text: str) -> MonolixProject:
    """
    Parse Monolix project text directly.

    Args:
        text: Raw .mlxtran file text

    Returns:
        MonolixProject with parsed structure
    """
    jl = _require_julia()
    result = jl.NeoPKPD.parse_monolix_project(text)
    return _convert_monolix_project(result, text)


# ============================================================================
# Validation Functions
# ============================================================================

def check_unsupported_constructs(ctl: NONMEMControlFile) -> List[UnsupportedConstruct]:
    """
    Check a parsed NONMEM control file for unsupported constructs.

    Args:
        ctl: Parsed NONMEMControlFile

    Returns:
        List of UnsupportedConstruct for each unsupported feature found

    Example:
        >>> ctl = parse_nonmem_control("run001.ctl")
        >>> unsupported = check_unsupported_constructs(ctl)
        >>> for u in unsupported:
        ...     print(f"Warning: {u.message}")
    """
    unsupported = []

    # Check ADVAN number
    if ctl.subroutines:
        advan = ctl.subroutines.advan
        if advan in [5, 6, 7, 8, 9, 12, 13]:
            unsupported.append(UnsupportedConstruct(
                f"ADVAN{advan}",
                "$SUBROUTINES",
                f"ADVAN{advan}",
                f"ADVAN{advan} (general linear/nonlinear ODE) is not supported. "
                "Only ADVAN1-4, 10, 11 with predefined models are supported."
            ))

    # Check $PK block
    if ctl.pk_block:
        for line in ctl.pk_block.raw_code:
            line_upper = line.upper().strip()

            if "IF" in line_upper and "(" in line_upper:
                unsupported.append(UnsupportedConstruct(
                    "IF statement", "$PK", line
                ))

            if "ALAG" in line_upper:
                unsupported.append(UnsupportedConstruct(
                    "Absorption lag time (ALAG)", "$PK", line
                ))

            if "MTIME" in line_upper:
                unsupported.append(UnsupportedConstruct(
                    "Model event time (MTIME)", "$PK", line
                ))

    # Check $ERROR block
    if ctl.error_block:
        for line in ctl.error_block.raw_code:
            line_upper = line.upper().strip()

            if "IF" in line_upper and "(" in line_upper:
                unsupported.append(UnsupportedConstruct(
                    "IF statement", "$ERROR", line
                ))

            if "CALL" in line_upper:
                unsupported.append(UnsupportedConstruct(
                    "CALL statement", "$ERROR", line
                ))

    # Check raw text for unsupported sections
    raw_upper = ctl.raw_text.upper()

    if "$DES" in raw_upper:
        unsupported.append(UnsupportedConstruct(
            "Custom differential equations ($DES)",
            "$DES",
            "",
            "Custom differential equations ($DES) are not supported. "
            "Use predefined ADVAN models."
        ))

    if "$MODEL" in raw_upper:
        unsupported.append(UnsupportedConstruct(
            "Custom compartment model ($MODEL)",
            "$MODEL",
            "",
            "Custom compartment models ($MODEL) are not supported. "
            "Use predefined ADVAN models."
        ))

    if "$MIX" in raw_upper:
        unsupported.append(UnsupportedConstruct(
            "Mixture model ($MIX)",
            "$MIX",
            "",
            "Mixture models ($MIX) are not supported."
        ))

    return unsupported


def validate_nonmem_import(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Validate a NONMEM control file for import compatibility.

    Returns validation results including supported features and warnings.

    Args:
        path: Path to NONMEM control file

    Returns:
        Dict with keys:
        - valid: Whether import is likely to succeed
        - model_kind: Detected model type (if valid)
        - parameters: Expected parameter names
        - warnings: List of warning messages
        - unsupported: List of unsupported constructs
    """
    ctl = parse_nonmem_control(path)
    unsupported = check_unsupported_constructs(ctl)

    warnings = []
    model_kind = None
    parameters = []

    if ctl.subroutines:
        mapping = get_model_mapping(ctl.subroutines.advan, ctl.subroutines.trans)
        if mapping:
            model_kind, parameters = mapping
        else:
            warnings.append(
                f"ADVAN{ctl.subroutines.advan}/TRANS{ctl.subroutines.trans} "
                "combination not directly supported"
            )

    valid = len(unsupported) == 0 and model_kind is not None

    return {
        "valid": valid,
        "model_kind": model_kind,
        "parameters": parameters,
        "warnings": warnings,
        "unsupported": [u.message for u in unsupported],
    }


def validate_monolix_import(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Validate a Monolix project file for import compatibility.

    Args:
        path: Path to Monolix project file

    Returns:
        Dict with validation results
    """
    project = parse_monolix_project(path)

    warnings = []
    model_kind = None

    if project.model and project.model.model_type:
        model_kind = get_monolix_model_mapping(project.model.model_type.model)
        if not model_kind:
            warnings.append(
                f"Model '{project.model.model_type.model}' not directly supported"
            )

    valid = model_kind is not None

    return {
        "valid": valid,
        "model_kind": model_kind,
        "parameters": [p.name for p in project.parameters],
        "warnings": warnings,
        "unsupported": [],
    }


# ============================================================================
# Conversion Helpers
# ============================================================================

def _convert_julia_import_result(result, format: str, path: str) -> ImportedModel:
    """Convert Julia NONMEMConversionResult or MonolixConversionResult to Python ImportedModel."""
    jl = _require_julia()

    # Check for conversion errors
    if hasattr(result, 'errors') and result.errors and len(result.errors) > 0:
        raise ValueError(f"Import failed: {'; '.join(str(e) for e in result.errors)}")

    model_spec = result.model_spec
    if model_spec is None:
        raise ValueError("Import failed: no model spec generated")

    # Extract model kind from the 'kind' field (ModelSpec.kind)
    # Use Julia's typeof and nameof to get the actual type name
    model_kind_obj = model_spec.kind
    try:
        # Get the Julia type name using Julia's typeof and nameof
        model_kind = str(jl.seval("x -> nameof(typeof(x))")(model_kind_obj))
    except Exception:
        # Fallback to Python type name
        model_kind = type(model_kind_obj).__name__

    # Extract params as dict
    params = {}
    model_params = model_spec.params
    try:
        # Try to get Julia field names using Julia's fieldnames function
        field_names = jl.seval("fieldnames")(jl.seval("typeof")(model_params))
        for field in field_names:
            val = getattr(model_params, str(field))
            if isinstance(val, (int, float)):
                params[str(field)] = float(val)
    except Exception:
        # Fallback: try Python attribute access
        for attr in dir(model_params):
            if not attr.startswith('_'):
                try:
                    val = getattr(model_params, attr)
                    if isinstance(val, (int, float)):
                        params[attr] = float(val)
                except Exception:
                    pass

    # Extract theta info
    theta_init = list(params.values())
    theta_names = list(params.keys())

    # Extract omega from IIV spec
    omega_init = []
    omega_names = []
    iiv_spec = result.iiv_spec
    if iiv_spec is not None:
        try:
            omegas = iiv_spec.omegas
            for k, v in omegas.items():
                omega_names.append(str(k))
                omega_init.append([float(v)])
        except Exception:
            pass

    # Extract sigma from error spec
    sigma_type = "proportional"
    sigma_init = 0.1
    error_spec = result.error_spec
    if error_spec is not None:
        try:
            sigma_type = str(error_spec.kind) if hasattr(error_spec, 'kind') else "proportional"
            sigma_init = float(error_spec.sigma) if hasattr(error_spec, 'sigma') else 0.1
        except Exception:
            pass

    # Get warnings
    warnings = [str(w) for w in result.warnings] if result.warnings else []

    return ImportedModel(
        source_format=format,
        source_file=path,
        model_kind=model_kind,
        params=params,
        theta_init=theta_init,
        theta_names=theta_names,
        omega_init=omega_init if omega_init else [[0.09]],
        omega_names=omega_names if omega_names else ["eta_1"],
        sigma_type=sigma_type,
        sigma_init=sigma_init,
        warnings=warnings,
        metadata={"parameter_mapping": dict(result.parameter_mapping) if hasattr(result, 'parameter_mapping') else {}},
    )


def _convert_nonmem_control(result, raw_text: str = "") -> NONMEMControlFile:
    """Convert Julia NONMEM control file to Python dataclass."""
    # Convert THETAs
    thetas = []
    if hasattr(result, 'thetas'):
        for spec in result.thetas:
            thetas.append(THETASpec(
                init=float(spec.init),
                lower=float(spec.lower) if hasattr(spec, 'lower') and spec.lower is not None else float('-inf'),
                upper=float(spec.upper) if hasattr(spec, 'upper') and spec.upper is not None else float('inf'),
                fixed=bool(spec.fixed) if hasattr(spec, 'fixed') else False,
                name=str(spec.name) if hasattr(spec, 'name') and spec.name else "",
            ))

    # Convert OMEGAs
    omegas = []
    if hasattr(result, 'omegas'):
        for block in result.omegas:
            omegas.append(OMEGABlock(
                values=[float(v) for v in block.values],
                structure=OmegaStructure(str(block.structure)) if hasattr(block, 'structure') else OmegaStructure.DIAGONAL,
                fixed=bool(block.fixed) if hasattr(block, 'fixed') else False,
            ))

    # Convert SIGMAs
    sigmas = []
    if hasattr(result, 'sigmas'):
        for block in result.sigmas:
            sigmas.append(SIGMABlock(
                values=[float(v) for v in block.values],
                structure=OmegaStructure(str(block.structure)) if hasattr(block, 'structure') else OmegaStructure.DIAGONAL,
                fixed=bool(block.fixed) if hasattr(block, 'fixed') else False,
            ))

    # Convert subroutines
    subroutines = None
    if hasattr(result, 'subroutines') and result.subroutines is not None:
        sub = result.subroutines
        subroutines = SubroutineSpec(
            advan=int(sub.advan),
            trans=int(sub.trans) if hasattr(sub, 'trans') else 1,
            other=[str(o) for o in sub.other] if hasattr(sub, 'other') else [],
        )

    # Convert data spec
    data = None
    if hasattr(result, 'data') and result.data is not None:
        d = result.data
        data = DataSpec(
            filename=str(d.filename),
            ignore=[str(i) for i in d.ignore] if hasattr(d, 'ignore') else [],
            accept=[str(a) for a in d.accept] if hasattr(d, 'accept') else [],
        )

    # Convert input columns
    input_columns = []
    if hasattr(result, 'input'):
        for col in result.input:
            input_columns.append(InputColumn(
                name=str(col.name),
                drop=bool(col.drop) if hasattr(col, 'drop') else False,
                alias=str(col.alias) if hasattr(col, 'alias') and col.alias else "",
            ))

    # Convert PK block
    pk_block = None
    pk_code = [str(line) for line in result.pk_code] if hasattr(result, 'pk_code') else []
    if pk_code:
        pk_block = PKBlock(raw_code=pk_code)

    # Convert error block
    error_block = None
    error_code = [str(line) for line in result.error_code] if hasattr(result, 'error_code') else []
    if error_code:
        error_block = ErrorBlock(raw_code=error_code)

    # Convert estimation
    estimation = {}
    if hasattr(result, 'estimation') and result.estimation:
        for k, v in result.estimation.items():
            estimation[str(k)] = v

    # Convert tables
    tables = []
    if hasattr(result, 'tables'):
        for t in result.tables:
            tables.append(dict(t))

    return NONMEMControlFile(
        problem=str(result.problem) if hasattr(result, 'problem') else "",
        data=data,
        input_columns=input_columns,
        subroutines=subroutines,
        thetas=thetas,
        omegas=omegas,
        sigmas=sigmas,
        pk_block=pk_block,
        error_block=error_block,
        estimation=estimation,
        tables=tables,
        raw_text=raw_text,
    )


def _convert_monolix_project(result, raw_text: str = "") -> MonolixProject:
    """Convert Julia Monolix project to Python dataclass."""
    # Convert model
    model = None
    if hasattr(result, 'model') and result.model is not None:
        m = result.model
        model_type = None
        if hasattr(m, 'model_type') and m.model_type is not None:
            mt = m.model_type
            model_type = MonolixModelType(
                lib=str(mt.lib) if hasattr(mt, 'lib') else "",
                model=str(mt.model) if hasattr(mt, 'model') else "",
            )
        model = MonolixStructuralModel(
            model_type=model_type,
            admin_type=str(m.admin_type) if hasattr(m, 'admin_type') else "iv",
            n_compartments=int(m.n_compartments) if hasattr(m, 'n_compartments') else 1,
            elimination=str(m.elimination) if hasattr(m, 'elimination') else "linear",
            absorption=str(m.absorption) if hasattr(m, 'absorption') else "bolus",
            has_lag=bool(m.has_lag) if hasattr(m, 'has_lag') else False,
            has_bioavailability=bool(m.has_bioavailability) if hasattr(m, 'has_bioavailability') else False,
        )

    # Convert data
    data = None
    if hasattr(result, 'data') and result.data is not None:
        d = result.data
        data = MonolixDataset(
            filename=str(d.filename),
            header_types={str(k): str(v) for k, v in d.header_types.items()} if hasattr(d, 'header_types') else {},
            id_column=str(d.id_column) if hasattr(d, 'id_column') else "ID",
            time_column=str(d.time_column) if hasattr(d, 'time_column') else "TIME",
            observation_column=str(d.observation_column) if hasattr(d, 'observation_column') else "DV",
            dose_column=str(d.dose_column) if hasattr(d, 'dose_column') else "AMT",
            rate_column=str(d.rate_column) if hasattr(d, 'rate_column') else "RATE",
        )

    # Convert parameters
    parameters = []
    if hasattr(result, 'parameters'):
        for p in result.parameters:
            parameters.append(MonolixParameter(
                name=str(p.name),
                value=float(p.value),
                fixed=bool(p.fixed) if hasattr(p, 'fixed') else False,
                distribution=str(p.distribution) if hasattr(p, 'distribution') else "logNormal",
                omega=float(p.omega) if hasattr(p, 'omega') else 0.0,
                has_iiv=bool(p.has_iiv) if hasattr(p, 'has_iiv') else False,
            ))

    # Convert observations
    observations = []
    if hasattr(result, 'observations'):
        for o in result.observations:
            observations.append(MonolixObservation(
                name=str(o.name),
                type=str(o.type) if hasattr(o, 'type') else "continuous",
                error_model=str(o.error_model) if hasattr(o, 'error_model') else "combined",
                error_params=[float(p) for p in o.error_params] if hasattr(o, 'error_params') else [],
            ))

    return MonolixProject(
        description=str(result.description) if hasattr(result, 'description') else "",
        data=data,
        model=model,
        parameters=parameters,
        observations=observations,
        estimation_method=str(result.estimation_method) if hasattr(result, 'estimation_method') else "SAEM",
        raw_text=raw_text,
    )


# ============================================================================
# Exports
# ============================================================================

__all__ = [
    # Enums
    "ErrorModelType",
    "IIVTransformation",
    "CovariateEffectType",
    "OmegaStructure",
    # NONMEM types
    "THETASpec",
    "OMEGABlock",
    "SIGMABlock",
    "SubroutineSpec",
    "DataSpec",
    "InputColumn",
    "PKCovariateEffect",
    "PKAssignment",
    "ScalingFactor",
    "PKBlock",
    "ErrorBlock",
    "UnsupportedConstruct",
    "NONMEMControlFile",
    # Monolix types
    "MonolixModelType",
    "MonolixStructuralModel",
    "MonolixParameter",
    "MonolixObservation",
    "MonolixDataset",
    "MonolixProject",
    "UnsupportedMonolixConstruct",
    # Import result
    "ImportedModel",
    # Mappings
    "ADVAN_TRANS_MAP",
    "MONOLIX_MODEL_MAP",
    "get_model_mapping",
    "get_monolix_model_mapping",
    # Import functions
    "import_nonmem",
    "import_monolix",
    "import_model",
    # Parse functions
    "parse_nonmem_control",
    "parse_nonmem_control_text",
    "parse_monolix_project",
    "parse_monolix_project_text",
    # Validation
    "check_unsupported_constructs",
    "validate_nonmem_import",
    "validate_monolix_import",
]
