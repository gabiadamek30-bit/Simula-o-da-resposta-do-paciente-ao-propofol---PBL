"""
Tests for NeoPKPD NCA (Non-Compartmental Analysis) module.

These tests verify the Python bindings to the Julia NCA implementation.
Comprehensive coverage of:
- Primary exposure metrics
- AUC calculations
- Lambda-z estimation
- PK parameters
- Multiple dose metrics
- Bioequivalence analysis
- Reference-scaled BE (RSABE/ABEL)
"""

import math
import pytest
import neopkpd
from neopkpd.nca import (
    # Full workflow
    run_nca,
    NCAConfig,
    NCAResult,
    # Primary metrics
    nca_cmax,
    nca_tmax,
    nca_cmin,
    nca_clast,
    nca_tlast,
    nca_cavg,
    # Additional exposure
    nca_ctrough,
    nca_c_at_time,
    time_above_concentration,
    # AUC
    auc_0_t,
    auc_0_inf,
    auc_0_tau,
    auc_partial,
    aumc_0_t,
    aumc_0_inf,
    # Lambda-z
    estimate_lambda_z,
    nca_half_life,
    # PK parameters
    nca_mrt,
    nca_cl_f,
    nca_vz_f,
    nca_vss,
    nca_cl,
    nca_vz,
    nca_mrt_iv,
    nca_cl_ss,
    nca_vss_from_aumc,
    nca_vc,
    nca_mean_absorption_time,
    nca_c0_from_regression,
    nca_bioavailability,
    # Multiple dose
    nca_accumulation_index,
    nca_ptf,
    nca_swing,
    nca_linearity_index,
    nca_time_to_steady_state,
    nca_accumulation_predicted,
    nca_accumulation_cmax,
    nca_accumulation_cmin,
    nca_dose_normalized_auc,
    nca_dose_normalized_cmax,
    nca_time_to_steady_state_doses,
    # Bioequivalence
    bioequivalence_90ci,
    tost_analysis,
    be_conclusion,
    geometric_mean_ratio,
    geometric_mean,
    within_subject_cv,
    # Study designs
    BEStudyDesign,
    ReplicateDesign,
    RegulatoryGuidance,
    RSABEConfig,
)


# ============================================================================
# Test Data
# ============================================================================

# Standard oral PK profile
TEST_TIMES = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
TEST_CONC = [0.0, 0.82, 1.44, 1.62, 1.28, 0.94, 0.68, 0.36, 0.08]
TEST_DOSE = 100.0

# Multiple dose profile
SS_TIMES = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
SS_CONC = [2.5, 4.2, 3.8, 3.0, 2.8, 2.6, 2.5, 2.5]
SS_TAU = 12.0


# ============================================================================
# Fixtures
# ============================================================================

# Note: The 'init' fixture is now provided by conftest.py which handles
# Julia initialization with proper signal handling for CI environments.
# Tests use 'init' fixture from conftest.py to ensure Julia is ready.


# ============================================================================
# Primary Exposure Metrics Tests
# ============================================================================

class TestExposureMetrics:
    """Tests for primary exposure metrics."""

    def test_nca_cmax(self, init):
        """Test Cmax calculation."""
        cmax = nca_cmax(TEST_CONC)
        assert cmax == 1.62

    def test_nca_tmax(self, init):
        """Test Tmax calculation."""
        tmax = nca_tmax(TEST_TIMES, TEST_CONC)
        assert tmax == 2.0

    def test_nca_cmin(self, init):
        """Test Cmin calculation."""
        cmin = nca_cmin(SS_CONC)
        assert cmin == 2.5

    def test_nca_clast(self, init):
        """Test Clast calculation."""
        clast = nca_clast(TEST_TIMES, TEST_CONC)
        assert clast == 0.08

    def test_nca_tlast(self, init):
        """Test Tlast calculation."""
        tlast = nca_tlast(TEST_TIMES, TEST_CONC)
        assert tlast == 24.0

    def test_clast_with_lloq(self, init):
        """Test Clast with LLOQ."""
        clast = nca_clast(TEST_TIMES, TEST_CONC, lloq=0.1)
        assert clast == 0.36

    def test_tlast_with_lloq(self, init):
        """Test Tlast with LLOQ."""
        tlast = nca_tlast(TEST_TIMES, TEST_CONC, lloq=0.1)
        assert tlast == 12.0


# ============================================================================
# Lambda-z and AUC Tests
# ============================================================================

class TestLambdaZAndAUC:
    """Tests for lambda-z estimation and AUC calculations."""

    def test_estimate_lambda_z(self, init):
        """Test lambda-z estimation."""
        result = estimate_lambda_z(TEST_TIMES, TEST_CONC)

        assert result["lambda_z"] is not None
        assert result["lambda_z"] > 0.0
        assert result["t_half"] is not None
        assert result["t_half"] > 0.0
        assert result["r_squared"] is not None
        assert result["r_squared"] >= 0.9
        assert result["quality_flag"] in ["good", "warning", "insufficient"]

    def test_lambda_z_t_half_relationship(self, init):
        """Test that t_half = ln(2) / lambda_z."""
        result = estimate_lambda_z(TEST_TIMES, TEST_CONC)

        if result["lambda_z"] is not None and result["t_half"] is not None:
            import math
            expected_t_half = math.log(2) / result["lambda_z"]
            assert abs(result["t_half"] - expected_t_half) < 1e-10

    def test_nca_half_life(self, init):
        """Test half-life calculation from lambda_z."""
        import math
        lambda_z = 0.1
        t_half = nca_half_life(lambda_z)
        expected = math.log(2) / 0.1
        assert abs(t_half - expected) < 1e-10

    def test_auc_0_t(self, init):
        """Test AUC 0-t calculation."""
        auc = auc_0_t(TEST_TIMES, TEST_CONC)
        assert auc > 0.0

    def test_auc_methods(self, init):
        """Test different AUC calculation methods."""
        auc_lin = auc_0_t(TEST_TIMES, TEST_CONC, method="linear")
        auc_log = auc_0_t(TEST_TIMES, TEST_CONC, method="log_linear")
        auc_mixed = auc_0_t(TEST_TIMES, TEST_CONC, method="lin_log_mixed")

        # All should be positive
        assert auc_lin > 0.0
        assert auc_log > 0.0
        assert auc_mixed > 0.0

        # Results should be similar (within 20%)
        assert 0.8 < auc_lin / auc_log < 1.2


# ============================================================================
# PK Parameters Tests
# ============================================================================

class TestPKParameters:
    """Tests for secondary PK parameters."""

    def test_nca_mrt(self, init):
        """Test MRT calculation."""
        mrt = nca_mrt(100.0, 20.0)  # AUMC=100, AUC=20
        assert mrt == 5.0

    def test_nca_mrt_iv_infusion(self, init):
        """Test MRT adjustment for IV infusion."""
        mrt = nca_mrt(100.0, 20.0, route="iv_infusion", t_inf=1.0)
        assert mrt == 4.5  # 5.0 - 1.0/2

    def test_nca_cl_f(self, init):
        """Test CL/F calculation."""
        cl_f = nca_cl_f(100.0, 20.0)  # dose=100, auc=20
        assert cl_f == 5.0

    def test_nca_vz_f(self, init):
        """Test Vz/F calculation."""
        vz_f = nca_vz_f(100.0, 0.1, 20.0)  # dose=100, lambda_z=0.1, auc=20
        assert vz_f == 50.0


# ============================================================================
# Multiple Dose Metrics Tests
# ============================================================================

class TestMultipleDoseMetrics:
    """Tests for multiple dose metrics."""

    def test_nca_accumulation_index(self, init):
        """Test accumulation index calculation."""
        rac = nca_accumulation_index(25.0, 20.0)  # AUC_ss=25, AUC_sd=20
        assert rac == 1.25

    def test_nca_ptf(self, init):
        """Test PTF calculation."""
        ptf = nca_ptf(4.2, 2.5, 3.0)  # Cmax=4.2, Cmin=2.5, Cavg=3.0
        expected = 100.0 * (4.2 - 2.5) / 3.0
        assert abs(ptf - expected) < 1e-10

    def test_nca_swing(self, init):
        """Test Swing calculation."""
        swing = nca_swing(4.2, 2.5)  # Cmax=4.2, Cmin=2.5
        expected = 100.0 * (4.2 - 2.5) / 2.5
        assert abs(swing - expected) < 1e-10


# ============================================================================
# Bioequivalence Tests
# ============================================================================

class TestBioequivalence:
    """Tests for bioequivalence analysis."""

    # Paired crossover data
    test_values = [20.0, 22.0, 18.0, 25.0, 21.0, 19.0, 23.0, 20.0]
    ref_values = [19.0, 21.0, 17.0, 24.0, 20.0, 18.0, 22.0, 19.0]

    def test_geometric_mean(self, init):
        """Test geometric mean calculation."""
        gm = geometric_mean(self.test_values)
        assert gm > 0.0

        import math
        expected = math.exp(sum(math.log(v) for v in self.test_values) / len(self.test_values))
        assert abs(gm - expected) < 1e-10

    def test_geometric_mean_ratio(self, init):
        """Test GMR calculation."""
        gmr = geometric_mean_ratio(self.test_values, self.ref_values)
        assert gmr > 0.0
        assert abs(gmr - 1.05) < 0.1  # Close to 1 for similar values

    def test_within_subject_cv(self, init):
        """Test within-subject CV calculation."""
        cv = within_subject_cv(self.test_values, self.ref_values)
        assert cv > 0.0
        assert cv < 100.0

    def test_bioequivalence_90ci(self, init):
        """Test 90% CI calculation."""
        result = bioequivalence_90ci(self.test_values, self.ref_values)

        assert "gmr" in result
        assert "ci_lower" in result
        assert "ci_upper" in result
        assert "cv_intra" in result
        assert "n" in result

        assert result["gmr"] > 0.0
        assert result["ci_lower"] < result["gmr"]
        assert result["ci_upper"] > result["gmr"]
        assert result["n"] == len(self.test_values)

    def test_tost_analysis(self, init):
        """Test TOST analysis."""
        result = tost_analysis(self.test_values, self.ref_values)

        assert "conclusion" in result
        assert result["conclusion"] in ["bioequivalent", "not_bioequivalent"]

    def test_be_conclusion(self, init):
        """Test BE conclusion determination."""
        # Within limits
        assert be_conclusion(0.85, 1.20) == "bioequivalent"

        # Outside limits
        assert be_conclusion(0.75, 1.10) == "inconclusive"
        assert be_conclusion(0.90, 1.30) == "inconclusive"

        # Completely outside
        assert be_conclusion(0.70, 0.78) == "not_bioequivalent"
        assert be_conclusion(1.30, 1.50) == "not_bioequivalent"


# ============================================================================
# Full NCA Workflow Tests
# ============================================================================

class TestFullNCAWorkflow:
    """Tests for full NCA workflow."""

    def test_run_nca_single_dose(self, init):
        """Test full NCA for single dose."""
        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE)

        assert isinstance(result, NCAResult)
        assert result.cmax == 1.62
        assert result.tmax == 2.0
        assert result.clast == 0.08
        assert result.tlast == 24.0
        assert result.auc_0_t > 0.0

    def test_run_nca_multiple_dose(self, init):
        """Test full NCA for multiple dose."""
        result = run_nca(
            SS_TIMES, SS_CONC, TEST_DOSE,
            dosing_type="multiple",
            tau=SS_TAU
        )

        assert result.cmin is not None
        assert result.cavg is not None
        assert result.auc_0_tau is not None

    def test_run_nca_iv_bolus(self, init):
        """Test full NCA for IV bolus."""
        result = run_nca(
            TEST_TIMES, TEST_CONC, TEST_DOSE,
            route="iv_bolus"
        )

        assert result.cmax > 0.0
        assert result.metadata["route"] == "iv_bolus"

    def test_nca_config(self, init):
        """Test NCA with custom configuration."""
        config = NCAConfig(
            method="log_linear",
            lambda_z_min_points=4,
            lambda_z_r2_threshold=0.85,
        )

        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE, config=config)
        assert result.auc_0_t > 0.0

    def test_nca_result_to_dict(self, init):
        """Test NCAResult to dict conversion."""
        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE)
        result_dict = result.to_dict()

        assert isinstance(result_dict, dict)
        assert "cmax" in result_dict
        assert "auc_0_t" in result_dict
        assert "quality_flags" in result_dict

    def test_quality_flags(self, init):
        """Test that quality flags are properly returned."""
        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE)

        assert isinstance(result.quality_flags, list)
        assert isinstance(result.warnings, list)
        assert isinstance(result.metadata, dict)


# ============================================================================
# Edge Cases Tests
# ============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_minimum_data_points(self, init):
        """Test with minimum data points."""
        t = [0.0, 1.0, 2.0]
        c = [0.0, 1.0, 0.5]

        result = run_nca(t, c, 100.0)
        assert result.cmax == 1.0

    def test_constant_concentration(self, init):
        """Test with constant concentration."""
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [1.0, 1.0, 1.0, 1.0, 1.0]

        result = run_nca(t, c, 100.0)
        assert result.cmax == 1.0


# ============================================================================
# Additional PK Parameters Tests
# ============================================================================

class TestAdditionalPKParameters:
    """Tests for additional PK parameters."""

    def test_nca_vss(self, init):
        """Test Vss calculation."""
        vss = nca_vss(5.0, 10.0)  # CL=5, MRT=10
        assert vss == 50.0

    def test_nca_cl(self, init):
        """Test CL (IV) calculation."""
        cl = nca_cl(100.0, 20.0)
        assert cl == 5.0

    def test_nca_vz(self, init):
        """Test Vz (IV) calculation."""
        vz = nca_vz(100.0, 0.1, 20.0)
        assert vz == 50.0

    def test_nca_mrt_iv(self, init):
        """Test MRT_iv from extravascular MRT."""
        mrt_iv = nca_mrt_iv(10.0, 2.0)
        assert mrt_iv == 8.0

    def test_nca_cl_ss(self, init):
        """Test steady-state clearance."""
        cl_ss = nca_cl_ss(100.0, 25.0)
        assert cl_ss == 4.0

    def test_nca_vss_from_aumc(self, init):
        """Test Vss from moment curves."""
        vss = nca_vss_from_aumc(100.0, 20.0, 200.0)
        expected = 100.0 * 200.0 / (20.0 ** 2)
        assert abs(vss - expected) < 0.001

    def test_nca_vc(self, init):
        """Test central volume."""
        vc = nca_vc(100.0, 10.0)
        assert vc == 10.0

    def test_nca_mean_absorption_time(self, init):
        """Test MAT calculation."""
        mat = nca_mean_absorption_time(10.0, 6.0)
        assert mat == 4.0

    def test_nca_bioavailability(self, init):
        """Test bioavailability."""
        f = nca_bioavailability(18.0, 100.0, 20.0, 100.0)
        assert abs(f - 0.9) < 0.001

    def test_nca_c0_from_regression(self, init):
        """Test C0 from regression intercept."""
        c0 = nca_c0_from_regression(math.log(10.0))
        assert abs(c0 - 10.0) < 0.001


# ============================================================================
# Additional Multiple Dose Tests
# ============================================================================

class TestAdditionalMultipleDoseMetrics:
    """Tests for additional multiple dose metrics."""

    def test_nca_accumulation_predicted(self, init):
        """Test predicted accumulation."""
        rac = nca_accumulation_predicted(0.1, 12.0)
        expected = 1.0 / (1.0 - math.exp(-0.1 * 12.0))
        assert abs(rac - expected) < 0.001

    def test_nca_accumulation_cmax(self, init):
        """Test Cmax accumulation ratio."""
        rac = nca_accumulation_cmax(3.0, 2.5)
        assert rac == 1.2

    def test_nca_accumulation_cmin(self, init):
        """Test Cmin accumulation ratio."""
        rac = nca_accumulation_cmin(0.5, 0.4)
        assert rac == 1.25

    def test_nca_dose_normalized_auc(self, init):
        """Test dose-normalized AUC."""
        auc_dn = nca_dose_normalized_auc(20.0, 100.0)
        assert auc_dn == 0.2

    def test_nca_dose_normalized_cmax(self, init):
        """Test dose-normalized Cmax."""
        cmax_dn = nca_dose_normalized_cmax(2.5, 100.0)
        assert cmax_dn == 0.025

    def test_nca_time_to_steady_state(self, init):
        """Test time to steady state."""
        t_ss = nca_time_to_steady_state(0.1, fraction=0.90)
        expected = -math.log(0.10) / 0.1
        assert abs(t_ss - expected) < 0.001

    def test_nca_time_to_steady_state_doses(self, init):
        """Test number of doses to steady state."""
        n_doses = nca_time_to_steady_state_doses(0.1, 8.0, fraction=0.90)
        assert isinstance(n_doses, int)
        assert n_doses > 0

    def test_nca_linearity_index(self, init):
        """Test dose linearity assessment."""
        doses = [50.0, 100.0, 200.0]
        aucs = [10.0, 20.0, 40.0]
        result = nca_linearity_index(doses, aucs)

        assert "beta" in result
        assert "r_squared" in result
        assert "is_linear" in result
        assert abs(result["beta"] - 1.0) < 0.1


# ============================================================================
# Additional AUC Tests
# ============================================================================

class TestAdditionalAUC:
    """Tests for additional AUC functions."""

    def test_auc_0_inf(self, init):
        """Test AUC0-inf calculation."""
        lz = estimate_lambda_z(TEST_TIMES, TEST_CONC)
        if lz["lambda_z"] is not None:
            auc_inf, pct = auc_0_inf(TEST_TIMES, TEST_CONC, lz["lambda_z"], 0.08)
            assert auc_inf > auc_0_t(TEST_TIMES, TEST_CONC)
            assert 0 <= pct <= 100

    def test_auc_0_tau(self, init):
        """Test AUC over dosing interval."""
        auc = auc_0_tau(SS_TIMES, SS_CONC, SS_TAU)
        assert auc > 0

    def test_auc_partial(self, init):
        """Test partial AUC."""
        auc = auc_partial(TEST_TIMES, TEST_CONC, 2.0, 8.0)
        assert auc > 0

    def test_aumc_0_t(self, init):
        """Test AUMC0-t."""
        aumc = aumc_0_t(TEST_TIMES, TEST_CONC)
        assert aumc > 0

    def test_aumc_0_inf(self, init):
        """Test AUMC0-inf."""
        lz = estimate_lambda_z(TEST_TIMES, TEST_CONC)
        if lz["lambda_z"] is not None:
            aumc_inf = aumc_0_inf(TEST_TIMES, TEST_CONC, lz["lambda_z"], 0.08, 24.0)
            assert aumc_inf > aumc_0_t(TEST_TIMES, TEST_CONC)


# ============================================================================
# Study Design and RSABE Tests
# ============================================================================

class TestStudyDesignsAndRSABE:
    """Tests for study design enums and RSABE configuration."""

    def test_be_study_design_enum(self, init):
        """Test BE study design enum values."""
        assert BEStudyDesign.CROSSOVER_2X2 == "crossover_2x2"
        assert BEStudyDesign.CROSSOVER_2X4 == "crossover_2x4"
        assert BEStudyDesign.PARALLEL_GROUP == "parallel_group"

    def test_replicate_design_enum(self, init):
        """Test replicate design enum values."""
        assert ReplicateDesign.PARTIAL_REPLICATE_3X3 == "partial_replicate_3x3"
        assert ReplicateDesign.FULL_REPLICATE_2X4 == "full_replicate_2x4"
        assert ReplicateDesign.FULL_REPLICATE_2X3 == "full_replicate_2x3"

    def test_regulatory_guidance_enum(self, init):
        """Test regulatory guidance enum values."""
        assert RegulatoryGuidance.FDA == "fda"
        assert RegulatoryGuidance.EMA == "ema"
        assert RegulatoryGuidance.HEALTH_CANADA == "health_canada"

    def test_rsabe_config(self, init):
        """Test RSABE configuration dataclass."""
        config = RSABEConfig(
            guidance=RegulatoryGuidance.FDA,
            design=ReplicateDesign.FULL_REPLICATE_2X4,
            parameter="cmax",
            alpha=0.05
        )
        assert config.guidance == RegulatoryGuidance.FDA
        assert config.design == ReplicateDesign.FULL_REPLICATE_2X4
        assert config.parameter == "cmax"
        assert config.alpha == 0.05


# ============================================================================
# Import Tests
# ============================================================================

class TestImports:
    """Test that all NCA functions are properly importable."""

    def test_import_from_main_package(self):
        """Test imports from main neopkpd package."""
        from neopkpd import (
            run_nca, NCAConfig, NCAResult,
            nca_cmax, nca_tmax, nca_cmin,
            auc_0_t, auc_0_inf, auc_0_tau,
            estimate_lambda_z, nca_half_life,
            nca_mrt, nca_cl_f, nca_vz_f, nca_vss,
            nca_accumulation_index, nca_ptf, nca_swing,
            bioequivalence_90ci, tost_analysis, be_conclusion,
            BEStudyDesign, ReplicateDesign, RegulatoryGuidance,
        )
        assert callable(run_nca)
        assert callable(nca_cmax)

    def test_import_from_nca_module(self):
        """Test imports from neopkpd.nca module."""
        from neopkpd.nca import (
            nca_ctrough, nca_c_at_time, time_above_concentration,
            nca_cl, nca_vz, nca_mrt_iv, nca_cl_ss,
            nca_vss_from_aumc, nca_vc, nca_mean_absorption_time,
            nca_accumulation_predicted, nca_accumulation_cmax,
            nca_dose_normalized_auc, nca_dose_normalized_cmax,
            RSABEConfig, RSABEResult, ABELResult,
        )
        assert callable(nca_ctrough)
        assert callable(nca_c_at_time)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
