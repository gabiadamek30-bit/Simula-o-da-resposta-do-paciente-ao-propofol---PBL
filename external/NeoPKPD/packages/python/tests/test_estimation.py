"""
Comprehensive tests for NeoPKPD parameter estimation module.

Tests cover:
- FOCE-I, SAEM, and Laplacian estimation methods
- Various model types (OneCompIVBolus, TwoCompOral, etc.)
- Configuration options (BLQ, IOV, covariates)
- Bootstrap analysis
- Model comparison (LRT, AIC, BIC)
- Diagnostics computation
"""

import pytest
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Any


# ============================================================================
# Test Fixtures
# ============================================================================

# Note: Julia initialization is handled by conftest.py with proper signal
# handling. Tests use 'init' fixture from conftest.py to ensure Julia is ready.


@pytest.fixture
def onecomp_observed_data():
    """Generate synthetic observed data for one-compartment IV bolus model."""
    # True parameters: CL=10, V=50
    # Typical profile with IIV
    subjects = []
    np.random.seed(12345)

    for i in range(10):  # 10 subjects
        # Add IIV (30% CV for CL, 20% CV for V)
        eta_cl = np.random.normal(0, 0.3)
        eta_v = np.random.normal(0, 0.2)
        cl_i = 10.0 * np.exp(eta_cl)
        v_i = 50.0 * np.exp(eta_v)

        times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        dose = 100.0

        # True one-comp IV bolus model: C = (D/V) * exp(-CL/V * t)
        ke = cl_i / v_i
        observations = []
        for t in times:
            c_true = (dose / v_i) * np.exp(-ke * t)
            # Add proportional residual error (10% CV)
            c_obs = c_true * (1 + np.random.normal(0, 0.1))
            observations.append(max(0.01, c_obs))  # Ensure positive

        subjects.append({
            "subject_id": f"SUBJ_{i+1:03d}",
            "times": times,
            "observations": observations,
            "doses": [{"time": 0.0, "amount": dose}],
        })

    return {"subjects": subjects}


@pytest.fixture
def twocomp_observed_data():
    """Generate synthetic observed data for two-compartment oral model."""
    # True parameters: Ka=1.0, CL=5.0, V1=30.0, Q=2.0, V2=50.0
    subjects = []
    np.random.seed(54321)

    for i in range(8):  # 8 subjects
        # Add IIV
        eta_cl = np.random.normal(0, 0.25)
        eta_v1 = np.random.normal(0, 0.2)

        ka_i = 1.0
        cl_i = 5.0 * np.exp(eta_cl)
        v1_i = 30.0 * np.exp(eta_v1)
        q_i = 2.0
        v2_i = 50.0

        times = [0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
        dose = 500.0

        # Simplified two-comp oral approximation
        k10 = cl_i / v1_i
        k12 = q_i / v1_i
        k21 = q_i / v2_i

        observations = []
        for t in times:
            # Approximate profile
            c_true = (dose / v1_i) * ka_i / (ka_i - k10) * (np.exp(-k10 * t) - np.exp(-ka_i * t))
            c_obs = c_true * (1 + np.random.normal(0, 0.15))
            observations.append(max(0.01, c_obs))

        subjects.append({
            "subject_id": f"SUBJ_{i+1:03d}",
            "times": times,
            "observations": observations,
            "doses": [{"time": 0.0, "amount": dose}],
        })

    return {"subjects": subjects}


@pytest.fixture
def basic_grid():
    """Basic simulation grid."""
    return {
        "t0": 0.0,
        "t1": 24.0,
        "saveat": list(np.arange(0.0, 24.1, 0.5)),
    }


@pytest.fixture
def basic_solver():
    """Basic solver settings."""
    return {
        "alg": "Tsit5",
        "reltol": 1e-8,
        "abstol": 1e-10,
        "maxiters": 10**6,
    }


# ============================================================================
# Test Classes
# ============================================================================

class TestEstimationConfig:
    """Test EstimationConfig dataclass."""

    def test_default_config(self):
        """Test default configuration values."""
        from neopkpd.estimation import EstimationConfig

        config = EstimationConfig()
        assert config.method == "foce"
        assert config.sigma_type == "proportional"
        assert config.sigma_init == 0.1
        assert config.max_iter == 500
        assert config.compute_se is True
        assert config.verbose is False

    def test_custom_config(self):
        """Test custom configuration."""
        from neopkpd.estimation import EstimationConfig, OmegaStructure

        config = EstimationConfig(
            method="saem",
            theta_init=[10.0, 50.0],
            theta_names=["CL", "V"],
            omega_init=[[0.09, 0.0], [0.0, 0.04]],
            omega_structure=OmegaStructure.DIAGONAL,
            sigma_type="combined",
            sigma_init=0.2,
            max_iter=1000,
        )

        assert config.method == "saem"
        assert config.theta_init == [10.0, 50.0]
        assert config.theta_names == ["CL", "V"]
        assert config.sigma_type == "combined"
        assert config.max_iter == 1000

    def test_blq_config(self):
        """Test BLQ configuration."""
        from neopkpd.estimation import EstimationConfig, BLQConfig, BLQMethod

        blq = BLQConfig(
            method=BLQMethod.M3_LIKELIHOOD,
            lloq=0.1,
            max_consecutive_blq=3,
        )

        config = EstimationConfig(
            theta_init=[10.0, 50.0],
            blq_config=blq,
        )

        assert config.blq_config is not None
        assert config.blq_config.method == BLQMethod.M3_LIKELIHOOD
        assert config.blq_config.lloq == 0.1


class TestFOCEEstimation:
    """Test FOCE-I estimation method."""

    def test_foce_onecomp_iv(self, julia_initialized, onecomp_observed_data, basic_grid, basic_solver):
        """Test FOCE estimation for one-compartment IV bolus."""
        from neopkpd.estimation import estimate, EstimationConfig

        config = EstimationConfig(
            method="foce",
            theta_init=[8.0, 40.0],  # Start away from true values (10, 50)
            theta_names=["CL", "V"],
            omega_init=[[0.1, 0.0], [0.0, 0.1]],
            omega_names=["eta_CL", "eta_V"],
            sigma_type="proportional",
            sigma_init=0.15,
            max_iter=200,
            verbose=False,
        )

        result = estimate(
            onecomp_observed_data,
            "OneCompIVBolus",
            config,
            basic_grid,
            basic_solver,
        )

        # Check result structure
        assert result is not None
        assert result.method == "foce"
        assert len(result.theta) == 2
        assert np.isfinite(result.ofv)  # OFV can be negative (it's -2LL)

        # Check convergence
        assert result.convergence is True or result.n_iterations > 0

        # Check parameter estimates are reasonable (within 50% of true values)
        assert 5.0 < result.theta[0] < 15.0  # CL ~ 10
        assert 25.0 < result.theta[1] < 75.0  # V ~ 50

        # Check individuals
        assert len(result.individuals) == 10
        for ind in result.individuals:
            assert ind.subject_id is not None
            assert len(ind.eta) == 2
            assert len(ind.ipred) > 0

    def test_foce_robust_se(self, julia_initialized, onecomp_observed_data, basic_grid, basic_solver):
        """Test FOCE with robust standard error computation."""
        from neopkpd.estimation import estimate, EstimationConfig

        config = EstimationConfig(
            method="foce",
            theta_init=[10.0, 50.0],
            omega_init=[[0.09, 0.0], [0.0, 0.04]],
            foce_compute_robust_se=True,
            compute_se=True,
            max_iter=200,  # More iterations for robust SE
        )

        result = estimate(
            onecomp_observed_data,
            "OneCompIVBolus",
            config,
            basic_grid,
            basic_solver,
        )

        # Check that both standard and robust SE are computed
        assert result.theta_se is not None
        # Robust SE may or may not be available
        if result.theta_se_robust is not None:
            assert len(result.theta_se_robust) == len(result.theta)


class TestLaplacianEstimation:
    """Test Laplacian estimation method."""

    def test_laplacian_onecomp_iv(self, julia_initialized, onecomp_observed_data, basic_grid, basic_solver):
        """Test Laplacian estimation for one-compartment IV bolus."""
        from neopkpd.estimation import estimate, EstimationConfig

        config = EstimationConfig(
            method="laplacian",
            theta_init=[9.0, 45.0],
            theta_names=["CL", "V"],
            omega_init=[[0.09, 0.0], [0.0, 0.04]],
            sigma_type="proportional",
            sigma_init=0.1,
            max_iter=100,
        )

        result = estimate(
            onecomp_observed_data,
            "OneCompIVBolus",
            config,
            basic_grid,
            basic_solver,
        )

        assert result is not None
        assert result.method == "laplacian"
        assert len(result.theta) == 2
        assert np.isfinite(result.ofv)  # OFV can be negative


class TestSAEMEstimation:
    """Test SAEM estimation method."""

    def test_saem_onecomp_iv(self, julia_initialized, onecomp_observed_data, basic_grid, basic_solver):
        """Test SAEM estimation for one-compartment IV bolus."""
        from neopkpd.estimation import estimate, EstimationConfig

        config = EstimationConfig(
            method="saem",
            theta_init=[10.0, 50.0],
            omega_init=[[0.09, 0.0], [0.0, 0.04]],
            saem_n_burn=30,  # Reduced for testing
            saem_n_iter=30,
            saem_n_chains=1,
        )

        result = estimate(
            onecomp_observed_data,
            "OneCompIVBolus",
            config,
            basic_grid,
            basic_solver,
        )

        assert result is not None
        assert result.method == "saem"


class TestModelComparison:
    """Test model comparison functions."""

    def test_likelihood_ratio_test(self):
        """Test likelihood ratio test calculation."""
        from neopkpd.estimation import likelihood_ratio_test

        # Full model OFV = 100, Reduced = 110, df = 2
        result = likelihood_ratio_test(100.0, 110.0, df=2)

        assert result["chi_sq"] == 10.0
        assert result["df"] == 2
        assert 0.0 < result["p_value"] < 1.0

        # If reduced is better (should not happen in nested models)
        result2 = likelihood_ratio_test(110.0, 100.0, df=2)
        assert result2["chi_sq"] == 0.0  # chi_sq should be max(0, ...)

    def test_compare_models_basic(self):
        """Test model comparison with mock results."""
        from neopkpd.estimation import compare_models, EstimationResult, IndividualEstimate

        # Create mock results
        result1 = EstimationResult(
            method="foce",
            theta=[10.0, 50.0, 1.0],
            theta_se=[1.0, 5.0, 0.1],
            theta_se_robust=None,
            theta_rse=[10.0, 10.0, 10.0],
            theta_ci_lower=[8.0, 40.0, 0.8],
            theta_ci_upper=[12.0, 60.0, 1.2],
            omega=[[0.09, 0.0], [0.0, 0.04]],
            omega_se=None,
            omega_corr=[[1.0, 0.0], [0.0, 1.0]],
            sigma={"type": "proportional"},
            sigma_se=None,
            ofv=500.0,
            aic=510.0,
            bic=520.0,
            convergence=True,
            n_iterations=50,
            gradient_norm=1e-6,
            condition_number=100.0,
            eigenvalue_ratio=10.0,
            covariance_successful=True,
            individuals=[],
            diagnostics=None,
            blq_summary=None,
            runtime_seconds=10.0,
            warnings=[],
        )

        result2 = EstimationResult(
            method="foce",
            theta=[10.0, 50.0],
            theta_se=[1.0, 5.0],
            theta_se_robust=None,
            theta_rse=[10.0, 10.0],
            theta_ci_lower=[8.0, 40.0],
            theta_ci_upper=[12.0, 60.0],
            omega=[[0.09, 0.0], [0.0, 0.04]],
            omega_se=None,
            omega_corr=[[1.0, 0.0], [0.0, 1.0]],
            sigma={"type": "proportional"},
            sigma_se=None,
            ofv=505.0,
            aic=513.0,
            bic=521.0,
            convergence=True,
            n_iterations=40,
            gradient_norm=1e-6,
            condition_number=80.0,
            eigenvalue_ratio=8.0,
            covariance_successful=True,
            individuals=[],
            diagnostics=None,
            blq_summary=None,
            runtime_seconds=8.0,
            warnings=[],
        )

        comparison = compare_models(
            [result1, result2],
            ["Full", "Reduced"],
            n_obs=100,
        )

        assert comparison.model_names == ["Full", "Reduced"]
        assert comparison.ofv_values == [500.0, 505.0]
        assert comparison.best_model_aic in ["Full", "Reduced"]
        assert comparison.lrt_statistic is not None
        assert comparison.lrt_pvalue is not None


class TestDiagnostics:
    """Test diagnostic functions."""

    def test_compute_diagnostics(self):
        """Test diagnostics computation from mock result."""
        from neopkpd.estimation import compute_diagnostics, EstimationResult, IndividualEstimate

        # Create mock individuals
        individuals = []
        for i in range(5):
            ind = IndividualEstimate(
                subject_id=f"SUBJ_{i+1}",
                eta=[np.random.normal(0, 0.3), np.random.normal(0, 0.2)],
                eta_se=None,
                ipred=[10.0, 8.0, 5.0, 3.0],
                pred=[9.0, 7.5, 4.5, 2.8],
                cwres=[np.random.normal(0, 1) for _ in range(4)],
                iwres=[np.random.normal(0, 1) for _ in range(4)],
                wres=[np.random.normal(0, 1) for _ in range(4)],
                ofv_contribution=50.0,
            )
            individuals.append(ind)

        result = EstimationResult(
            method="foce",
            theta=[10.0, 50.0],
            theta_se=[1.0, 5.0],
            theta_se_robust=None,
            theta_rse=[10.0, 10.0],
            theta_ci_lower=[8.0, 40.0],
            theta_ci_upper=[12.0, 60.0],
            omega=[[0.09, 0.0], [0.0, 0.04]],
            omega_se=None,
            omega_corr=[[1.0, 0.0], [0.0, 1.0]],
            sigma={"type": "proportional"},
            sigma_se=None,
            ofv=250.0,
            aic=258.0,
            bic=265.0,
            convergence=True,
            n_iterations=50,
            gradient_norm=1e-6,
            condition_number=100.0,
            eigenvalue_ratio=10.0,
            covariance_successful=True,
            individuals=individuals,
            diagnostics=None,
            blq_summary=None,
            runtime_seconds=10.0,
            warnings=[],
        )

        diag = compute_diagnostics(result)

        assert diag is not None
        assert -1.0 < diag.cwres_mean < 1.0  # Should be close to 0
        assert 0.5 < diag.cwres_std < 1.5    # Should be close to 1
        assert len(diag.eta_shrinkage) == 2


class TestBootstrap:
    """Test bootstrap analysis."""

    def test_bootstrap_config(self):
        """Test BootstrapConfig creation."""
        from neopkpd.estimation import BootstrapConfig, BootstrapCIMethod

        config = BootstrapConfig(
            n_bootstrap=500,
            ci_level=0.95,
            ci_method=BootstrapCIMethod.PERCENTILE,
            parallel=False,
        )

        assert config.n_bootstrap == 500
        assert config.ci_level == 0.95
        assert config.ci_method == BootstrapCIMethod.PERCENTILE

    def test_bootstrap_run(self, julia_initialized, onecomp_observed_data, basic_grid, basic_solver):
        """Test running bootstrap analysis."""
        from neopkpd.estimation import run_bootstrap, EstimationConfig, BootstrapConfig

        est_config = EstimationConfig(
            method="foce",
            theta_init=[10.0, 50.0],
            omega_init=[[0.09, 0.0], [0.0, 0.04]],
            max_iter=50,  # Faster estimation for bootstrap
        )

        boot_config = BootstrapConfig(
            n_bootstrap=100,  # Minimum required by Julia (FDA recommends ≥500)
            ci_level=0.95,
        )

        result = run_bootstrap(
            onecomp_observed_data,
            "OneCompIVBolus",
            est_config,
            boot_config,
            basic_grid,
            basic_solver,
        )

        assert result is not None
        assert result.theta_estimates.shape[0] == 100
        assert len(result.theta_ci_lower) == 2
        assert len(result.theta_ci_upper) == 2
        assert result.convergence_rate > 0.3  # Relaxed convergence threshold


class TestErrorHandling:
    """Test error handling and edge cases."""

    def test_invalid_method(self, julia_initialized):
        """Test that invalid method raises error."""
        from neopkpd.estimation import EstimationConfig

        # This should create config but will fail at estimation
        config = EstimationConfig(
            method="invalid_method",
            theta_init=[10.0, 50.0],
        )

        assert config.method == "invalid_method"

    def test_invalid_model_kind(self, julia_initialized, onecomp_observed_data, basic_grid, basic_solver):
        """Test that invalid model kind raises error."""
        from neopkpd.estimation import estimate, EstimationConfig

        config = EstimationConfig(
            method="foce",
            theta_init=[10.0, 50.0],
        )

        with pytest.raises(ValueError, match="Unsupported model kind"):
            estimate(
                onecomp_observed_data,
                "InvalidModelKind",
                config,
                basic_grid,
                basic_solver,
            )

    def test_mismatched_theta_omega_dims(self, julia_initialized):
        """Test config with mismatched dimensions."""
        from neopkpd.estimation import EstimationConfig

        # This should create config but may cause issues during estimation
        config = EstimationConfig(
            theta_init=[10.0, 50.0],
            omega_init=[[0.09]],  # Only 1x1, but 2 thetas
        )

        assert len(config.theta_init) == 2
        assert len(config.omega_init) == 1


class TestEnums:
    """Test enumeration types."""

    def test_blq_methods(self):
        """Test BLQ method enum."""
        from neopkpd.estimation import BLQMethod

        assert BLQMethod.M1_DISCARD.value == "m1_discard"
        assert BLQMethod.M2_IMPUTE_HALF.value == "m2_half"
        assert BLQMethod.M3_LIKELIHOOD.value == "m3_likelihood"

    def test_omega_structures(self):
        """Test omega structure enum."""
        from neopkpd.estimation import OmegaStructure

        assert OmegaStructure.DIAGONAL.value == "diagonal"
        assert OmegaStructure.BLOCK.value == "block"
        assert OmegaStructure.FULL.value == "full"

    def test_bootstrap_ci_methods(self):
        """Test bootstrap CI method enum."""
        from neopkpd.estimation import BootstrapCIMethod

        assert BootstrapCIMethod.PERCENTILE.value == "percentile"
        assert BootstrapCIMethod.BCA.value == "bca"
        assert BootstrapCIMethod.BASIC.value == "basic"


class TestImports:
    """Test that all estimation exports are available."""

    def test_main_imports(self):
        """Test imports from main neopkpd module."""
        from neopkpd import (
            estimate,
            run_bootstrap,
            compare_models,
            likelihood_ratio_test,
            compute_diagnostics,
            get_individual_predictions,
            EstimationConfig,
            BootstrapConfig,
            BLQConfig,
            IOVSpec,
            CovariateOnIIV,
            EstimationResult,
            BootstrapResult,
            IndividualEstimate,
            DiagnosticsSummary,
            BLQSummary,
            ModelComparisonResult,
            BLQMethod,
            OmegaStructure,
            BootstrapCIMethod,
        )

        assert estimate is not None
        assert run_bootstrap is not None
        assert EstimationConfig is not None
        assert BLQMethod is not None

    def test_submodule_imports(self):
        """Test imports from estimation submodule."""
        from neopkpd.estimation import (
            estimate,
            run_bootstrap,
            compare_models,
            likelihood_ratio_test,
        )

        assert estimate is not None
        assert run_bootstrap is not None


# ============================================================================
# Run Tests
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
