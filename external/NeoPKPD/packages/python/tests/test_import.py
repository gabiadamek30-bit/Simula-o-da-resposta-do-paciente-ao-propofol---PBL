"""
Comprehensive tests for NeoPKPD model import module.

Tests cover:
- NONMEM control file parsing (ADVAN1-4, 10, 11)
- Monolix project file parsing (.mlxtran)
- Type definitions (dataclasses and enums)
- Model mapping functions
- Unsupported construct detection
- Import validation
"""

import pytest
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Any


# ============================================================================
# Test Data - Sample Control Files
# ============================================================================

NONMEM_ADVAN1_ONECOMP = """
$PROBLEM One-compartment IV bolus PK model
$INPUT ID TIME DV AMT MDV EVID
$DATA ../data.csv IGNORE=@

$SUBROUTINES ADVAN1 TRANS2

$PK
  TVCL = THETA(1)
  TVV  = THETA(2)

  CL = TVCL * EXP(ETA(1))
  V  = TVV  * EXP(ETA(2))

  S1 = V

$ERROR
  IPRED = F
  W = IPRED * THETA(3)
  IRES = DV - IPRED
  IWRES = IRES / W
  Y = IPRED + W * ERR(1)

$THETA
  (0, 5.0)    ; CL (L/h)
  (0, 50.0)   ; V (L)
  (0, 0.10)   ; Proportional error

$OMEGA BLOCK(2)
  0.09        ; IIV CL (CV ~30%)
  0.03 0.0625 ; IIV V (CV ~25%)

$SIGMA
  1 FIX       ; Proportional residual

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 PRINT=10 NOABORT
$COVARIANCE PRINT=E
"""

# Version with IF statement for testing unsupported construct detection
NONMEM_ADVAN1_WITH_IF = """
$PROBLEM One-compartment IV bolus PK model
$INPUT ID TIME DV AMT MDV EVID
$DATA ../data.csv IGNORE=@

$SUBROUTINES ADVAN1 TRANS2

$PK
  TVCL = THETA(1)
  TVV  = THETA(2)
  CL = TVCL * EXP(ETA(1))
  V  = TVV  * EXP(ETA(2))
  S1 = V

$ERROR
  IPRED = F
  W = IPRED * THETA(3)
  IF(W.EQ.0) W = 1
  Y = IPRED + W * ERR(1)

$THETA
  (0, 5.0)
  (0, 50.0)
  (0, 0.10)

$OMEGA
  0.09
  0.0625

$SIGMA
  1 FIX

$ESTIMATION METHOD=1
"""

NONMEM_ADVAN2_ORAL = """
$PROBLEM One-compartment oral PK model with first-order absorption
$INPUT ID TIME DV AMT MDV EVID
$DATA ../data.csv IGNORE=@

$SUBROUTINES ADVAN2 TRANS2

$PK
  TVKA = THETA(1)
  TVCL = THETA(2)
  TVV  = THETA(3)

  KA = TVKA * EXP(ETA(1))
  CL = TVCL * EXP(ETA(2))
  V  = TVV  * EXP(ETA(3))

  S2 = V

$ERROR
  IPRED = F
  W = SQRT(THETA(4)**2 + (THETA(5)*IPRED)**2)
  IRES = DV - IPRED
  IWRES = IRES / W
  Y = IPRED + W * ERR(1)

$THETA
  (0, 1.5)    ; Ka (1/h)
  (0, 5.0)    ; CL (L/h)
  (0, 50.0)   ; V (L)
  (0, 0.5)    ; Additive error (mg/L)
  (0, 0.10)   ; Proportional error

$OMEGA
  0.16        ; IIV Ka (CV ~40%)
  0.09        ; IIV CL (CV ~30%)
  0.0625      ; IIV V (CV ~25%)

$SIGMA
  1 FIX       ; Combined residual

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 PRINT=10 NOABORT
"""

NONMEM_ADVAN4_TWOCOMP = """
$PROBLEM Two-compartment IV bolus PK model
$INPUT ID TIME DV AMT MDV EVID
$DATA ../data.csv IGNORE=@

$SUBROUTINES ADVAN4 TRANS4

$PK
  TVCL = THETA(1)
  TVV1 = THETA(2)
  TVQ  = THETA(3)
  TVV2 = THETA(4)

  CL = TVCL * EXP(ETA(1))
  V1 = TVV1 * EXP(ETA(2))
  Q  = TVQ
  V2 = TVV2

  S1 = V1

$ERROR
  IPRED = F
  W = IPRED * THETA(5)
  IRES = DV - IPRED
  IWRES = IRES / W
  Y = IPRED + W * ERR(1)

$THETA
  (0, 5.0)    ; CL (L/h)
  (0, 10.0)   ; V1 (L)
  (0, 2.0)    ; Q (L/h)
  (0, 20.0)   ; V2 (L)
  (0, 0.15)   ; Proportional error

$OMEGA
  0.09        ; IIV CL (CV ~30%)
  0.0625      ; IIV V1 (CV ~25%)

$SIGMA
  1 FIX       ; Proportional residual

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 PRINT=10 NOABORT
"""

NONMEM_UNSUPPORTED = """
$PROBLEM Unsupported model with $DES
$INPUT ID TIME DV AMT EVID
$DATA data.csv

$SUBROUTINES ADVAN6 TRANS1

$MODEL
  COMP=(DEPOT,DEFDOSE)
  COMP=(CENTRAL)

$PK
  KA = THETA(1)
  CL = THETA(2)
  V = THETA(3)

$DES
  DADT(1) = -KA*A(1)
  DADT(2) = KA*A(1) - CL/V*A(2)

$ERROR
  Y = A(2)/V

$THETA
  1
  10
  100

$ESTIMATION METHOD=1
"""

MONOLIX_PK1CPT_ORAL = """
<DATAFILE>
[FILEINFO]
file = 'data.csv'
delimiter = comma
header = {ID, TIME, DV, AMT, EVID}

[CONTENT]
ID = {use=identifier}
TIME = {use=time}
DV = {use=observation, name=Cc, type=continuous}
AMT = {use=amount}
EVID = {use=eventidentifier}

<MODEL>
[INDIVIDUAL]
input = {ka_pop, V_pop, Cl_pop, omega_ka, omega_V, omega_Cl}

DEFINITION:
ka = {distribution=logNormal, typical=ka_pop, sd=omega_ka}
V = {distribution=logNormal, typical=V_pop, sd=omega_V}
Cl = {distribution=logNormal, typical=Cl_pop, sd=omega_Cl}

[LONGITUDINAL]
input = {a, b}

file = 'lib:oral1_1cpt_kaVCl.txt'

DEFINITION:
Cc = {distribution=normal, prediction=Cc, errorModel=combined2(a, b)}

<FIT>
data = Cc
model = Cc

<PARAMETER>
ka_pop = {value=1.5, method=MLE}
V_pop = {value=50, method=MLE}
Cl_pop = {value=5, method=MLE}
omega_ka = {value=0.4, method=MLE}
omega_V = {value=0.25, method=MLE}
omega_Cl = {value=0.3, method=MLE}
a = {value=0.5, method=MLE}
b = {value=0.1, method=MLE}

<MONOLIX>
[TASKS]
populationParameters()
conditionalModes(linearization = true)
"""

MONOLIX_PK2CPT_IV = """
<DATAFILE>
[FILEINFO]
file = 'data.csv'
delimiter = comma
header = {ID, TIME, DV, AMT, EVID}

[CONTENT]
ID = {use=identifier}
TIME = {use=time}
DV = {use=observation, name=Cc, type=continuous}
AMT = {use=amount}
EVID = {use=eventidentifier}

<MODEL>
[INDIVIDUAL]
input = {Cl_pop, V1_pop, Q_pop, V2_pop, omega_Cl, omega_V1}

DEFINITION:
Cl = {distribution=logNormal, typical=Cl_pop, sd=omega_Cl}
V1 = {distribution=logNormal, typical=V1_pop, sd=omega_V1}
Q = {distribution=logNormal, typical=Q_pop, no-variability}
V2 = {distribution=logNormal, typical=V2_pop, no-variability}

[LONGITUDINAL]
input = {b}

file = 'lib:bolus2_2cpt_ClV1QV2.txt'

DEFINITION:
Cc = {distribution=normal, prediction=Cc, errorModel=proportional(b)}

<FIT>
data = Cc
model = Cc

<PARAMETER>
Cl_pop = {value=5, method=MLE}
V1_pop = {value=10, method=MLE}
Q_pop = {value=2, method=MLE}
V2_pop = {value=20, method=MLE}
omega_Cl = {value=0.3, method=MLE}
omega_V1 = {value=0.25, method=MLE}
b = {value=0.15, method=MLE}

<MONOLIX>
[TASKS]
populationParameters()
"""


# ============================================================================
# Test Fixtures
# ============================================================================

# Note: Julia initialization is handled by conftest.py with proper signal
# handling. Tests use 'init' fixture from conftest.py to ensure Julia is ready.


@pytest.fixture
def temp_ctl_file():
    """Create a temporary NONMEM control file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.ctl', delete=False) as f:
        f.write(NONMEM_ADVAN1_ONECOMP)
        f.flush()
        yield f.name
    os.unlink(f.name)


@pytest.fixture
def temp_mlxtran_file():
    """Create a temporary Monolix project file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.mlxtran', delete=False) as f:
        f.write(MONOLIX_PK1CPT_ORAL)
        f.flush()
        yield f.name
    os.unlink(f.name)


# ============================================================================
# Pure Python Tests (No Julia Required)
# ============================================================================

class TestEnums:
    """Test enum definitions."""

    def test_error_model_type_values(self):
        """Test ErrorModelType enum values."""
        from neopkpd.import_ import ErrorModelType

        assert ErrorModelType.PROPORTIONAL.value == "proportional"
        assert ErrorModelType.ADDITIVE.value == "additive"
        assert ErrorModelType.COMBINED.value == "combined"
        assert ErrorModelType.EXPONENTIAL.value == "exponential"
        assert ErrorModelType.UNKNOWN.value == "unknown"

    def test_iiv_transformation_values(self):
        """Test IIVTransformation enum values."""
        from neopkpd.import_ import IIVTransformation

        assert IIVTransformation.EXPONENTIAL.value == "exponential"
        assert IIVTransformation.ADDITIVE.value == "additive"
        assert IIVTransformation.PROPORTIONAL.value == "proportional"
        assert IIVTransformation.NONE.value == "none"

    def test_covariate_effect_type_values(self):
        """Test CovariateEffectType enum values."""
        from neopkpd.import_ import CovariateEffectType

        assert CovariateEffectType.POWER.value == "power"
        assert CovariateEffectType.LINEAR.value == "linear"
        assert CovariateEffectType.EXPONENTIAL.value == "exponential"

    def test_omega_structure_values(self):
        """Test OmegaStructure enum values."""
        from neopkpd.import_ import OmegaStructure

        assert OmegaStructure.DIAGONAL.value == "diagonal"
        assert OmegaStructure.BLOCK.value == "block"
        assert OmegaStructure.SAME.value == "same"


class TestNONMEMTypes:
    """Test NONMEM dataclass types."""

    def test_theta_spec_creation(self):
        """Test THETASpec dataclass."""
        from neopkpd.import_ import THETASpec

        theta = THETASpec(init=5.0, lower=0.0, upper=100.0, fixed=False, name="CL")
        assert theta.init == 5.0
        assert theta.lower == 0.0
        assert theta.upper == 100.0
        assert theta.fixed is False
        assert theta.name == "CL"

    def test_theta_spec_defaults(self):
        """Test THETASpec default values."""
        from neopkpd.import_ import THETASpec

        theta = THETASpec(init=10.0)
        assert theta.lower == float('-inf')
        assert theta.upper == float('inf')
        assert theta.fixed is False
        assert theta.name == ""

    def test_omega_block_creation(self):
        """Test OMEGABlock dataclass."""
        from neopkpd.import_ import OMEGABlock, OmegaStructure

        # Diagonal omega
        omega_diag = OMEGABlock(values=[0.09, 0.0625])
        assert omega_diag.dimension == 2
        assert omega_diag.structure == OmegaStructure.DIAGONAL

        # Block omega (2x2 has 3 elements: lower triangular)
        omega_block = OMEGABlock(values=[0.09, 0.03, 0.0625], structure=OmegaStructure.BLOCK)
        assert omega_block.dimension == 2

    def test_sigma_block_creation(self):
        """Test SIGMABlock dataclass."""
        from neopkpd.import_ import SIGMABlock

        sigma = SIGMABlock(values=[1.0], fixed=True)
        assert sigma.values == [1.0]
        assert sigma.fixed is True

    def test_subroutine_spec(self):
        """Test SubroutineSpec dataclass."""
        from neopkpd.import_ import SubroutineSpec

        sub = SubroutineSpec(advan=1, trans=2)
        assert sub.advan == 1
        assert sub.trans == 2
        assert sub.other == []

    def test_data_spec(self):
        """Test DataSpec dataclass."""
        from neopkpd.import_ import DataSpec

        data = DataSpec(filename="data.csv", ignore=["@", "C"])
        assert data.filename == "data.csv"
        assert data.ignore == ["@", "C"]
        assert data.accept == []

    def test_input_column(self):
        """Test InputColumn dataclass."""
        from neopkpd.import_ import InputColumn

        col = InputColumn(name="AMT", drop=False, alias="DOSE")
        assert col.name == "AMT"
        assert col.drop is False
        assert col.alias == "DOSE"

    def test_pk_covariate_effect(self):
        """Test PKCovariateEffect dataclass."""
        from neopkpd.import_ import PKCovariateEffect, CovariateEffectType

        effect = PKCovariateEffect(
            covariate="WT",
            theta_index=3,
            effect_type=CovariateEffectType.POWER,
            reference=70.0
        )
        assert effect.covariate == "WT"
        assert effect.theta_index == 3
        assert effect.effect_type == CovariateEffectType.POWER
        assert effect.reference == 70.0

    def test_pk_assignment(self):
        """Test PKAssignment dataclass."""
        from neopkpd.import_ import PKAssignment, IIVTransformation

        assignment = PKAssignment(
            target="CL",
            tv_theta=1,
            eta_index=1,
            transformation=IIVTransformation.EXPONENTIAL
        )
        assert assignment.target == "CL"
        assert assignment.tv_theta == 1
        assert assignment.eta_index == 1
        assert assignment.transformation == IIVTransformation.EXPONENTIAL

    def test_scaling_factor(self):
        """Test ScalingFactor dataclass."""
        from neopkpd.import_ import ScalingFactor

        scaling = ScalingFactor(compartment=1, parameter="V")
        assert scaling.compartment == 1
        assert scaling.parameter == "V"

    def test_pk_block(self):
        """Test PKBlock dataclass."""
        from neopkpd.import_ import PKBlock

        pk = PKBlock(raw_code=["TVCL = THETA(1)", "CL = TVCL * EXP(ETA(1))"])
        assert len(pk.raw_code) == 2
        assert pk.tv_definitions == {}
        assert pk.assignments == []

    def test_error_block(self):
        """Test ErrorBlock dataclass."""
        from neopkpd.import_ import ErrorBlock, ErrorModelType

        error = ErrorBlock(error_type=ErrorModelType.PROPORTIONAL)
        assert error.error_type == ErrorModelType.PROPORTIONAL
        assert error.theta_indices == []

    def test_unsupported_construct(self):
        """Test UnsupportedConstruct dataclass."""
        from neopkpd.import_ import UnsupportedConstruct

        unsup = UnsupportedConstruct(
            construct="IF statement",
            location="$PK",
            line="IF (TIME.GT.0) CL = CL * 2"
        )
        assert unsup.construct == "IF statement"
        assert unsup.location == "$PK"
        assert "is not supported" in unsup.message

    def test_nonmem_control_file(self):
        """Test NONMEMControlFile dataclass."""
        from neopkpd.import_ import NONMEMControlFile, THETASpec, SubroutineSpec

        thetas = [THETASpec(init=5.0), THETASpec(init=50.0)]
        subroutines = SubroutineSpec(advan=1, trans=2)

        ctl = NONMEMControlFile(
            problem="Test problem",
            subroutines=subroutines,
            thetas=thetas,
            raw_text=NONMEM_ADVAN1_ONECOMP
        )

        assert ctl.problem == "Test problem"
        assert ctl.advan == 1  # Legacy property
        assert ctl.trans == 2  # Legacy property
        assert len(ctl.theta_specs) == 2  # Legacy property


class TestMonolixTypes:
    """Test Monolix dataclass types."""

    def test_monolix_model_type(self):
        """Test MonolixModelType dataclass."""
        from neopkpd.import_ import MonolixModelType

        mt = MonolixModelType(lib="pklib", model="pk_oral1cpt_kaVCl_PLASMA")
        assert mt.lib == "pklib"
        assert mt.model == "pk_oral1cpt_kaVCl_PLASMA"

    def test_monolix_structural_model(self):
        """Test MonolixStructuralModel dataclass."""
        from neopkpd.import_ import MonolixStructuralModel, MonolixModelType

        mt = MonolixModelType(lib="pklib", model="pk_oral1cpt_kaVCl_PLASMA")
        model = MonolixStructuralModel(
            model_type=mt,
            admin_type="oral",
            n_compartments=1,
            elimination="linear",
            absorption="firstOrder"
        )
        assert model.admin_type == "oral"
        assert model.n_compartments == 1
        assert model.has_lag is False

    def test_monolix_parameter(self):
        """Test MonolixParameter dataclass."""
        from neopkpd.import_ import MonolixParameter

        param = MonolixParameter(
            name="ka",
            value=1.5,
            fixed=False,
            distribution="logNormal",
            omega=0.4,
            has_iiv=True
        )
        assert param.name == "ka"
        assert param.value == 1.5
        assert param.has_iiv is True

    def test_monolix_observation(self):
        """Test MonolixObservation dataclass."""
        from neopkpd.import_ import MonolixObservation

        obs = MonolixObservation(
            name="Cc",
            type="continuous",
            error_model="combined",
            error_params=[0.5, 0.1]
        )
        assert obs.name == "Cc"
        assert obs.error_model == "combined"

    def test_monolix_dataset(self):
        """Test MonolixDataset dataclass."""
        from neopkpd.import_ import MonolixDataset

        data = MonolixDataset(filename="data.csv", id_column="SUBJ")
        assert data.filename == "data.csv"
        assert data.id_column == "SUBJ"
        assert data.time_column == "TIME"  # Default

    def test_monolix_project(self):
        """Test MonolixProject dataclass."""
        from neopkpd.import_ import MonolixProject, MonolixParameter

        params = [
            MonolixParameter(name="ka", value=1.5),
            MonolixParameter(name="V", value=50.0),
        ]
        project = MonolixProject(
            description="Test project",
            parameters=params,
            estimation_method="SAEM"
        )
        assert project.description == "Test project"
        assert len(project.parameters) == 2
        assert project.estimation_method == "SAEM"

    def test_unsupported_monolix_construct(self):
        """Test UnsupportedMonolixConstruct dataclass."""
        from neopkpd.import_ import UnsupportedMonolixConstruct

        unsup = UnsupportedMonolixConstruct(
            construct="PD model",
            location="<MODEL>",
            line="file = 'pd_turnover.txt'"
        )
        assert "PD model" in unsup.message


class TestImportedModel:
    """Test ImportedModel dataclass."""

    def test_imported_model_creation(self):
        """Test ImportedModel dataclass."""
        from neopkpd.import_ import ImportedModel

        model = ImportedModel(
            source_format="nonmem",
            source_file="/path/to/run001.ctl",
            model_kind="OneCompIVBolus",
            params={"CL": 5.0, "V": 50.0},
            theta_init=[5.0, 50.0],
            theta_names=["CL", "V"],
            omega_init=[[0.09], [0.0625]],
            omega_names=["eta_CL", "eta_V"],
            sigma_type="proportional",
            sigma_init=0.1,
            warnings=[],
            metadata={}
        )
        assert model.source_format == "nonmem"
        assert model.model_kind == "OneCompIVBolus"
        assert model.params["CL"] == 5.0


class TestModelMappings:
    """Test model mapping functions and constants."""

    def test_advan_trans_map_contents(self):
        """Test ADVAN_TRANS_MAP contains expected mappings."""
        from neopkpd.import_ import ADVAN_TRANS_MAP

        # ADVAN1
        assert (1, 1) in ADVAN_TRANS_MAP
        assert (1, 2) in ADVAN_TRANS_MAP
        assert ADVAN_TRANS_MAP[(1, 2)][0] == "OneCompIVBolus"

        # ADVAN2
        assert (2, 1) in ADVAN_TRANS_MAP
        assert (2, 2) in ADVAN_TRANS_MAP
        assert ADVAN_TRANS_MAP[(2, 2)][0] == "OneCompOralFirstOrder"

        # ADVAN3/4
        assert (3, 4) in ADVAN_TRANS_MAP
        assert (4, 4) in ADVAN_TRANS_MAP

        # ADVAN10/11
        assert (10, 1) in ADVAN_TRANS_MAP
        assert (11, 4) in ADVAN_TRANS_MAP

    def test_get_model_mapping(self):
        """Test get_model_mapping function."""
        from neopkpd.import_ import get_model_mapping

        # Supported combinations
        result = get_model_mapping(1, 2)
        assert result is not None
        model_kind, params = result
        assert model_kind == "OneCompIVBolus"
        assert "CL" in params
        assert "V" in params

        result = get_model_mapping(2, 2)
        assert result is not None
        model_kind, params = result
        assert model_kind == "OneCompOralFirstOrder"
        assert "KA" in params

        # Unsupported combination
        result = get_model_mapping(6, 1)
        assert result is None

    def test_monolix_model_map_contents(self):
        """Test MONOLIX_MODEL_MAP contains expected mappings."""
        from neopkpd.import_ import MONOLIX_MODEL_MAP

        assert "pk_bolus1cpt_VCl_PLASMA" in MONOLIX_MODEL_MAP
        assert MONOLIX_MODEL_MAP["pk_bolus1cpt_VCl_PLASMA"] == "OneCompIVBolus"

        assert "pk_oral1cpt_1abs_kaVCl_PLASMA" in MONOLIX_MODEL_MAP
        assert MONOLIX_MODEL_MAP["pk_oral1cpt_1abs_kaVCl_PLASMA"] == "OneCompOralFirstOrder"

    def test_get_monolix_model_mapping_exact(self):
        """Test get_monolix_model_mapping with exact match."""
        from neopkpd.import_ import get_monolix_model_mapping

        result = get_monolix_model_mapping("pk_bolus1cpt_VCl_PLASMA")
        assert result == "OneCompIVBolus"

        result = get_monolix_model_mapping("pk_oral2cpt_1abs_kaV1ClQ2V2_PLASMA")
        assert result == "TwoCompOral"

    def test_get_monolix_model_mapping_pattern(self):
        """Test get_monolix_model_mapping with pattern matching."""
        from neopkpd.import_ import get_monolix_model_mapping

        # Pattern matching for variations
        result = get_monolix_model_mapping("pk_oral_custom_1cpt_model")
        assert result == "OneCompOralFirstOrder"

        result = get_monolix_model_mapping("pk_bolus_2cpt_custom")
        assert result == "TwoCompIVBolus"

        result = get_monolix_model_mapping("model_3cpt_something")
        assert result == "ThreeCompIVBolus"

    def test_get_monolix_model_mapping_unknown(self):
        """Test get_monolix_model_mapping with unknown model."""
        from neopkpd.import_ import get_monolix_model_mapping

        result = get_monolix_model_mapping("completely_unknown_model")
        assert result is None


# ============================================================================
# Julia-Dependent Tests
# ============================================================================

class TestNONMEMParsing:
    """Test NONMEM control file parsing (requires Julia)."""

    def test_parse_advan1_onecomp(self, julia_initialized, temp_ctl_file):
        """Test parsing ADVAN1 one-compartment model."""
        from neopkpd.import_ import parse_nonmem_control

        ctl = parse_nonmem_control(temp_ctl_file)

        assert ctl.problem == "One-compartment IV bolus PK model"
        assert ctl.subroutines is not None
        assert ctl.subroutines.advan == 1
        assert ctl.subroutines.trans == 2
        assert len(ctl.thetas) == 3  # CL, V, prop error

    def test_parse_nonmem_control_text(self, julia_initialized):
        """Test parsing control file text directly."""
        from neopkpd.import_ import parse_nonmem_control_text

        ctl = parse_nonmem_control_text(NONMEM_ADVAN2_ORAL)

        assert "oral" in ctl.problem.lower()
        assert ctl.subroutines.advan == 2
        assert len(ctl.thetas) == 5  # Ka, CL, V, add error, prop error

    def test_parse_advan4_twocomp(self, julia_initialized):
        """Test parsing ADVAN4 two-compartment model."""
        from neopkpd.import_ import parse_nonmem_control_text

        ctl = parse_nonmem_control_text(NONMEM_ADVAN4_TWOCOMP)

        assert ctl.subroutines.advan == 4
        assert ctl.subroutines.trans == 4
        assert len(ctl.thetas) == 5  # CL, V1, Q, V2, prop error

    def test_check_unsupported_constructs_advan6(self, julia_initialized):
        """Test detection of unsupported ADVAN6."""
        from neopkpd.import_ import parse_nonmem_control_text, check_unsupported_constructs

        ctl = parse_nonmem_control_text(NONMEM_UNSUPPORTED)
        unsupported = check_unsupported_constructs(ctl)

        # Should detect ADVAN6 as unsupported
        advan_issues = [u for u in unsupported if "ADVAN6" in u.construct]
        assert len(advan_issues) > 0

    def test_check_unsupported_constructs_des(self, julia_initialized):
        """Test detection of unsupported $DES."""
        from neopkpd.import_ import parse_nonmem_control_text, check_unsupported_constructs

        ctl = parse_nonmem_control_text(NONMEM_UNSUPPORTED)
        unsupported = check_unsupported_constructs(ctl)

        # Should detect $DES as unsupported
        des_issues = [u for u in unsupported if "$DES" in u.construct]
        assert len(des_issues) > 0

    def test_validate_nonmem_import(self, julia_initialized, temp_ctl_file):
        """Test NONMEM validation function."""
        from neopkpd.import_ import validate_nonmem_import

        result = validate_nonmem_import(temp_ctl_file)

        assert result["valid"] is True
        assert result["model_kind"] == "OneCompIVBolus"
        assert "CL" in result["parameters"]
        assert "V" in result["parameters"]


class TestMonolixParsing:
    """Test Monolix project file parsing (requires Julia)."""

    def test_parse_monolix_1cpt_oral(self, julia_initialized, temp_mlxtran_file):
        """Test parsing one-compartment oral Monolix project."""
        from neopkpd.import_ import parse_monolix_project

        project = parse_monolix_project(temp_mlxtran_file)

        # Check basic parsing
        assert project.data is not None
        assert project.data.filename == "data.csv"

    def test_parse_monolix_project_text(self, julia_initialized):
        """Test parsing Monolix project text directly."""
        from neopkpd.import_ import parse_monolix_project_text

        project = parse_monolix_project_text(MONOLIX_PK2CPT_IV)

        # Check parameters were parsed
        assert len(project.parameters) > 0
        param_names = [p.name for p in project.parameters]
        # Should have population parameters
        assert any("Cl" in name or "CL" in name.upper() for name in param_names)

    def test_validate_monolix_import(self, julia_initialized, temp_mlxtran_file):
        """Test Monolix validation function."""
        from neopkpd.import_ import validate_monolix_import

        result = validate_monolix_import(temp_mlxtran_file)

        # Should have parsed parameter names
        assert len(result["parameters"]) > 0


class TestImportFunctions:
    """Test main import functions (requires Julia)."""

    def test_import_nonmem_advan1(self, julia_initialized, temp_ctl_file):
        """Test importing NONMEM ADVAN1 model."""
        from neopkpd.import_ import import_nonmem

        doses = [{"time": 0.0, "amount": 100.0}]
        model = import_nonmem(temp_ctl_file, doses=doses)

        assert model.source_format == "nonmem"
        assert model.model_kind == "OneCompIVBolus"
        assert len(model.params) > 0
        # Should have CL and V parameters
        param_keys = [k.upper() for k in model.params.keys()]
        assert any("CL" in k for k in param_keys) or any("K" in k for k in param_keys)

    def test_import_monolix_1cpt(self, julia_initialized, temp_mlxtran_file):
        """Test importing Monolix one-compartment model."""
        from neopkpd.import_ import import_monolix

        doses = [{"time": 0.0, "amount": 100.0}]
        model = import_monolix(temp_mlxtran_file, doses=doses)

        assert model.source_format == "monolix"
        assert len(model.params) > 0

    def test_import_model_auto_detect_ctl(self, julia_initialized, temp_ctl_file):
        """Test import_model auto-detects .ctl files."""
        from neopkpd.import_ import import_model

        doses = [{"time": 0.0, "amount": 100.0}]
        model = import_model(temp_ctl_file, doses=doses)

        assert model.source_format == "nonmem"

    def test_import_model_auto_detect_mlxtran(self, julia_initialized, temp_mlxtran_file):
        """Test import_model auto-detects .mlxtran files."""
        from neopkpd.import_ import import_model

        doses = [{"time": 0.0, "amount": 100.0}]
        model = import_model(temp_mlxtran_file, doses=doses)

        assert model.source_format == "monolix"

    def test_import_model_explicit_format(self, julia_initialized, temp_ctl_file):
        """Test import_model with explicit format."""
        from neopkpd.import_ import import_model

        doses = [{"time": 0.0, "amount": 100.0}]
        model = import_model(temp_ctl_file, format="nonmem", doses=doses)

        assert model.source_format == "nonmem"

    def test_import_model_invalid_format(self, julia_initialized):
        """Test import_model raises error for invalid format."""
        from neopkpd.import_ import import_model

        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
            f.write(b"some text")
            f.flush()
            try:
                with pytest.raises(ValueError, match="Cannot auto-detect"):
                    import_model(f.name)
            finally:
                os.unlink(f.name)


class TestRealExampleFiles:
    """Test with real example files from docs/examples/import."""

    @pytest.fixture
    def example_dir(self):
        """Get the examples directory."""
        return Path("/Users/shramishkafle/Desktop/neopkpd/docs/examples/import")

    def test_import_real_advan1(self, julia_initialized, example_dir):
        """Test importing real ADVAN1 example."""
        from neopkpd.import_ import import_nonmem

        ctl_path = example_dir / "nonmem" / "01_advan1_onecomp" / "run001.ctl"
        if not ctl_path.exists():
            pytest.skip(f"Example file not found: {ctl_path}")

        doses = [{"time": 0.0, "amount": 100.0}]
        model = import_nonmem(str(ctl_path), doses=doses)

        assert model.model_kind == "OneCompIVBolus"

    def test_import_real_advan2(self, julia_initialized, example_dir):
        """Test importing real ADVAN2 example."""
        from neopkpd.import_ import import_nonmem

        ctl_path = example_dir / "nonmem" / "02_advan2_oral" / "run002.ctl"
        if not ctl_path.exists():
            pytest.skip(f"Example file not found: {ctl_path}")

        doses = [{"time": 0.0, "amount": 500.0}]
        model = import_nonmem(str(ctl_path), doses=doses)

        assert model.model_kind == "OneCompOralFirstOrder"

    def test_import_real_advan4(self, julia_initialized, example_dir):
        """Test importing real ADVAN4 example."""
        from neopkpd.import_ import import_nonmem

        ctl_path = example_dir / "nonmem" / "03_advan4_twocomp" / "run003.ctl"
        if not ctl_path.exists():
            pytest.skip(f"Example file not found: {ctl_path}")

        doses = [{"time": 0.0, "amount": 100.0}]
        model = import_nonmem(str(ctl_path), doses=doses)

        assert model.model_kind == "TwoCompOral"

    def test_import_real_monolix_1cpt(self, julia_initialized, example_dir):
        """Test importing real Monolix 1-cpt example."""
        from neopkpd.import_ import import_monolix

        mlx_path = example_dir / "monolix" / "01_pk1cpt_oral" / "project.mlxtran"
        if not mlx_path.exists():
            pytest.skip(f"Example file not found: {mlx_path}")

        doses = [{"time": 0.0, "amount": 100.0}]
        model = import_monolix(str(mlx_path), doses=doses)

        assert model.source_format == "monolix"

    def test_import_real_monolix_2cpt(self, julia_initialized, example_dir):
        """Test importing real Monolix 2-cpt example."""
        from neopkpd.import_ import import_monolix

        mlx_path = example_dir / "monolix" / "02_pk2cpt_iv" / "project.mlxtran"
        if not mlx_path.exists():
            pytest.skip(f"Example file not found: {mlx_path}")

        doses = [{"time": 0.0, "amount": 100.0}]
        model = import_monolix(str(mlx_path), doses=doses)

        assert model.source_format == "monolix"


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_theta_list(self):
        """Test NONMEMControlFile with empty theta list."""
        from neopkpd.import_ import NONMEMControlFile

        ctl = NONMEMControlFile(problem="Test")
        assert ctl.theta_specs == []

    def test_legacy_properties_without_subroutines(self):
        """Test legacy properties when subroutines is None."""
        from neopkpd.import_ import NONMEMControlFile

        ctl = NONMEMControlFile(problem="Test", subroutines=None)
        assert ctl.advan == 0
        assert ctl.trans == 1

    def test_legacy_properties_without_data(self):
        """Test legacy properties when data is None."""
        from neopkpd.import_ import NONMEMControlFile

        ctl = NONMEMControlFile(problem="Test", data=None)
        assert ctl.data_file is None

    def test_monolix_project_legacy_properties(self):
        """Test MonolixProject legacy properties."""
        from neopkpd.import_ import MonolixProject

        project = MonolixProject()
        assert project.model_type == "unknown"
        assert project.structural_model == ""
        assert project.data_file is None


class TestExports:
    """Test that all expected items are exported."""

    def test_all_exports_available(self):
        """Test that __all__ exports are importable."""
        from neopkpd import import_

        expected_exports = [
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

        for name in expected_exports:
            assert hasattr(import_, name), f"Missing export: {name}"
