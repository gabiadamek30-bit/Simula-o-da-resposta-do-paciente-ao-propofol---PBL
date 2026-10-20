"""
Tests for new visualization modules:
- viz/vpc.py
- viz/estimation.py
- viz/bootstrap.py
- viz/sensitivity.py
"""

import pytest
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Any, Optional


# ============================================================================
# Mock Data Classes for Testing
# ============================================================================

@dataclass
class MockVPCPercentileData:
    percentile: float
    observed: float
    simulated_median: float
    simulated_lower: float
    simulated_upper: float


@dataclass
class MockVPCBin:
    bin_id: int
    time_min: float
    time_max: float
    time_midpoint: float
    n_observed: int
    n_simulated: int
    percentiles: List[MockVPCPercentileData]


@dataclass
class MockBLQBinStats:
    bin_id: int
    n_total: int
    n_blq: int
    pct_blq_observed: float
    pct_blq_simulated_median: float
    pct_blq_simulated_lower: float
    pct_blq_simulated_upper: float


@dataclass
class MockVPCConfig:
    pi_levels: List[float] = None
    ci_level: float = 0.95
    n_simulations: int = 200
    prediction_corrected: bool = False

    def __post_init__(self):
        if self.pi_levels is None:
            self.pi_levels = [0.05, 0.50, 0.95]


@dataclass
class MockVPCResult:
    config: MockVPCConfig
    bins: List[MockVPCBin]
    n_subjects_observed: int
    n_observations_observed: int
    n_simulations: int
    strata: str = ""
    prediction_corrected: bool = False
    blq_stats: Optional[List[MockBLQBinStats]] = None


@dataclass
class MockStratifiedVPCResult:
    results: List[MockVPCResult]
    stratify_by: List[str]
    strata_names: List[str]


@dataclass
class MockDiagnosticsSummary:
    cwres_mean: float
    cwres_std: float
    iwres_mean: float
    iwres_std: float
    eta_shrinkage: List[float]
    epsilon_shrinkage: float
    condition_number: float
    eigenvalue_ratio: float


@dataclass
class MockIndividualEstimate:
    subject_id: str
    eta: List[float]
    eta_se: Optional[List[float]]
    ipred: List[float]
    pred: List[float]
    cwres: List[float]
    iwres: List[float]
    wres: List[float]
    ofv_contribution: float


@dataclass
class MockEstimationResult:
    method: str
    theta: List[float]
    theta_se: Optional[List[float]]
    theta_se_robust: Optional[List[float]]
    theta_rse: Optional[List[float]]
    theta_ci_lower: Optional[List[float]]
    theta_ci_upper: Optional[List[float]]
    omega: List[List[float]]
    omega_se: Optional[List[List[float]]]
    omega_corr: List[List[float]]
    sigma: Dict[str, Any]
    sigma_se: Optional[Dict[str, Any]]
    ofv: float
    aic: float
    bic: float
    convergence: bool
    n_iterations: int
    gradient_norm: float
    condition_number: float
    eigenvalue_ratio: float
    covariance_successful: bool
    individuals: List[MockIndividualEstimate]
    diagnostics: Optional[MockDiagnosticsSummary]
    blq_summary: Optional[Any] = None
    runtime_seconds: float = 1.0
    warnings: List[str] = None

    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []


@dataclass
class MockBootstrapResult:
    theta_estimates: np.ndarray
    theta_mean: List[float]
    theta_se: List[float]
    theta_rse: List[float]
    theta_ci_lower: List[float]
    theta_ci_upper: List[float]
    original_estimate: List[float]
    bias: List[float]
    bias_corrected: List[float]
    n_successful: int
    n_failed: int
    convergence_rate: float
    ci_level: float
    ci_method: str
    omega_mean: Optional[List[List[float]]] = None
    omega_se: Optional[List[List[float]]] = None
    omega_ci_lower: Optional[List[List[float]]] = None
    omega_ci_upper: Optional[List[List[float]]] = None
    sigma_mean: Optional[float] = None
    sigma_se: Optional[float] = None
    sigma_ci_lower: Optional[float] = None
    sigma_ci_upper: Optional[float] = None


@dataclass
class MockModelComparisonResult:
    model_names: List[str]
    ofv_values: List[float]
    aic_values: List[float]
    bic_values: List[float]
    n_params: List[int]
    best_model_aic: str
    best_model_bic: str
    lrt_statistic: Optional[float]
    lrt_pvalue: Optional[float]
    lrt_df: Optional[int]


# ============================================================================
# Helper Functions
# ============================================================================

def create_mock_vpc_result(n_bins: int = 8) -> MockVPCResult:
    """Create mock VPC result for testing."""
    config = MockVPCConfig()
    bins = []

    for i in range(n_bins):
        t_min = i * 3.0
        t_max = (i + 1) * 3.0
        t_mid = (t_min + t_max) / 2

        percentiles = []
        for level in config.pi_levels:
            base = 50.0 * np.exp(-0.1 * t_mid)
            percentiles.append(MockVPCPercentileData(
                percentile=level,
                observed=base * (0.8 + 0.4 * level),
                simulated_median=base * (0.75 + 0.5 * level),
                simulated_lower=base * (0.7 + 0.5 * level),
                simulated_upper=base * (0.8 + 0.5 * level),
            ))

        bins.append(MockVPCBin(
            bin_id=i + 1,
            time_min=t_min,
            time_max=t_max,
            time_midpoint=t_mid,
            n_observed=20,
            n_simulated=200,
            percentiles=percentiles,
        ))

    return MockVPCResult(
        config=config,
        bins=bins,
        n_subjects_observed=50,
        n_observations_observed=400,
        n_simulations=200,
    )


def create_mock_vpc_with_blq() -> MockVPCResult:
    """Create mock VPC result with BLQ stats."""
    result = create_mock_vpc_result()

    blq_stats = []
    for i, bin in enumerate(result.bins):
        blq_pct = min(0.3 + 0.1 * i, 0.8)  # Increasing BLQ over time
        blq_stats.append(MockBLQBinStats(
            bin_id=bin.bin_id,
            n_total=bin.n_observed,
            n_blq=int(bin.n_observed * blq_pct),
            pct_blq_observed=blq_pct,
            pct_blq_simulated_median=blq_pct * 0.9,
            pct_blq_simulated_lower=blq_pct * 0.7,
            pct_blq_simulated_upper=blq_pct * 1.1,
        ))

    result.blq_stats = blq_stats
    return result


def create_mock_stratified_vpc() -> MockStratifiedVPCResult:
    """Create mock stratified VPC result."""
    results = [
        create_mock_vpc_result(6),
        create_mock_vpc_result(6),
        create_mock_vpc_result(6),
    ]
    results[0].strata = "Low Weight"
    results[1].strata = "Medium Weight"
    results[2].strata = "High Weight"

    return MockStratifiedVPCResult(
        results=results,
        stratify_by=["WEIGHT_GROUP"],
        strata_names=["Low Weight", "Medium Weight", "High Weight"],
    )


def create_mock_estimation_result() -> MockEstimationResult:
    """Create mock estimation result for testing."""
    # Create individuals
    individuals = []
    for i in range(10):
        individuals.append(MockIndividualEstimate(
            subject_id=f"SUBJ{i+1}",
            eta=[np.random.normal(0, 0.3), np.random.normal(0, 0.2)],
            eta_se=[0.1, 0.08],
            ipred=list(np.exp(-0.1 * np.arange(24))),
            pred=list(np.exp(-0.1 * np.arange(24))),
            cwres=list(np.random.normal(0, 1, 24)),
            iwres=list(np.random.normal(0, 1, 24)),
            wres=list(np.random.normal(0, 1, 24)),
            ofv_contribution=np.random.uniform(10, 30),
        ))

    return MockEstimationResult(
        method="foce",
        theta=[10.0, 50.0, 1.0],
        theta_se=[1.0, 5.0, 0.1],
        theta_se_robust=[1.1, 5.5, 0.11],
        theta_rse=[10.0, 10.0, 10.0],
        theta_ci_lower=[8.0, 40.0, 0.8],
        theta_ci_upper=[12.0, 60.0, 1.2],
        omega=[[0.09, 0.01], [0.01, 0.04]],
        omega_se=[[0.02, 0.005], [0.005, 0.01]],
        omega_corr=[[1.0, 0.167], [0.167, 1.0]],
        sigma={"proportional": 0.1},
        sigma_se={"proportional": 0.01},
        ofv=500.0,
        aic=512.0,
        bic=520.0,
        convergence=True,
        n_iterations=150,
        gradient_norm=1e-5,
        condition_number=15.0,
        eigenvalue_ratio=0.01,
        covariance_successful=True,
        individuals=individuals,
        diagnostics=MockDiagnosticsSummary(
            cwres_mean=0.01,
            cwres_std=1.02,
            iwres_mean=-0.02,
            iwres_std=0.98,
            eta_shrinkage=[0.15, 0.20],
            epsilon_shrinkage=0.10,
            condition_number=15.0,
            eigenvalue_ratio=0.01,
        ),
    )


def create_mock_bootstrap_result() -> MockBootstrapResult:
    """Create mock bootstrap result for testing."""
    n_boot = 500
    n_params = 3

    # Generate bootstrap estimates
    theta_estimates = np.zeros((n_boot, n_params))
    true_values = [10.0, 50.0, 1.0]
    for i in range(n_params):
        theta_estimates[:, i] = np.random.normal(true_values[i], true_values[i] * 0.1, n_boot)

    theta_mean = list(np.mean(theta_estimates, axis=0))
    theta_se = list(np.std(theta_estimates, axis=0))

    return MockBootstrapResult(
        theta_estimates=theta_estimates,
        theta_mean=theta_mean,
        theta_se=theta_se,
        theta_rse=[se / m * 100 for m, se in zip(theta_mean, theta_se)],
        theta_ci_lower=list(np.percentile(theta_estimates, 2.5, axis=0)),
        theta_ci_upper=list(np.percentile(theta_estimates, 97.5, axis=0)),
        original_estimate=true_values,
        bias=[m - t for m, t in zip(theta_mean, true_values)],
        bias_corrected=[2 * t - m for t, m in zip(true_values, theta_mean)],
        n_successful=480,
        n_failed=20,
        convergence_rate=0.96,
        ci_level=0.95,
        ci_method="percentile",
    )


def create_mock_model_comparison() -> MockModelComparisonResult:
    """Create mock model comparison result."""
    return MockModelComparisonResult(
        model_names=["Base Model", "Full Model", "Reduced Model"],
        ofv_values=[520.0, 500.0, 515.0],
        aic_values=[530.0, 512.0, 523.0],
        bic_values=[540.0, 520.0, 530.0],
        n_params=[5, 6, 4],
        best_model_aic="Full Model",
        best_model_bic="Full Model",
        lrt_statistic=20.0,
        lrt_pvalue=0.0001,
        lrt_df=1,
    )


# ============================================================================
# VPC Visualization Tests
# ============================================================================

class TestVPCVisualization:
    """Tests for VPC visualization functions."""

    def test_plot_vpc_detailed(self):
        """Test detailed VPC plot."""
        from neopkpd.viz.vpc import plot_vpc

        vpc_result = create_mock_vpc_result()
        fig = plot_vpc(vpc_result, log_scale=False)
        assert fig is not None

    def test_plot_vpc_log_scale(self):
        """Test VPC with log scale."""
        from neopkpd.viz.vpc import plot_vpc

        vpc_result = create_mock_vpc_result()
        fig = plot_vpc(vpc_result, log_scale=True)
        assert fig is not None

    def test_plot_pcvpc(self):
        """Test prediction-corrected VPC plot."""
        from neopkpd.viz.vpc import plot_pcvpc

        vpc_result = create_mock_vpc_result()
        vpc_result.prediction_corrected = True
        fig = plot_pcvpc(vpc_result)
        assert fig is not None

    def test_plot_vpc_ci(self):
        """Test VPC CI plot."""
        from neopkpd.viz.vpc import plot_vpc_ci

        vpc_result = create_mock_vpc_result()
        fig = plot_vpc_ci(vpc_result)
        assert fig is not None

    def test_plot_stratified_vpc(self):
        """Test stratified VPC plot."""
        from neopkpd.viz.vpc import plot_stratified_vpc

        stratified_result = create_mock_stratified_vpc()
        fig = plot_stratified_vpc(stratified_result)
        assert fig is not None

    def test_plot_vpc_with_blq(self):
        """Test VPC with BLQ panel."""
        from neopkpd.viz.vpc import plot_vpc_with_blq

        vpc_result = create_mock_vpc_with_blq()
        fig = plot_vpc_with_blq(vpc_result)
        assert fig is not None


# ============================================================================
# Estimation Visualization Tests
# ============================================================================

class TestEstimationVisualization:
    """Tests for estimation visualization functions."""

    def test_plot_parameter_estimates(self):
        """Test parameter estimates forest plot."""
        from neopkpd.viz.estimation import plot_parameter_estimates

        result = create_mock_estimation_result()
        fig = plot_parameter_estimates(result)
        assert fig is not None

    def test_plot_parameter_estimates_with_omega(self):
        """Test parameter estimates with omega."""
        from neopkpd.viz.estimation import plot_parameter_estimates

        result = create_mock_estimation_result()
        fig = plot_parameter_estimates(result, include_omega=True)
        assert fig is not None

    def test_plot_omega_matrix(self):
        """Test omega matrix heatmap."""
        from neopkpd.viz.estimation import plot_omega_matrix

        result = create_mock_estimation_result()
        fig = plot_omega_matrix(result, show_correlation=True)
        assert fig is not None

    def test_plot_omega_covariance(self):
        """Test omega covariance heatmap."""
        from neopkpd.viz.estimation import plot_omega_matrix

        result = create_mock_estimation_result()
        fig = plot_omega_matrix(result, show_correlation=False)
        assert fig is not None

    def test_plot_convergence(self):
        """Test convergence plot."""
        from neopkpd.viz.estimation import plot_convergence

        result = create_mock_estimation_result()
        fig = plot_convergence(result)
        assert fig is not None

    def test_plot_shrinkage(self):
        """Test shrinkage bar chart."""
        from neopkpd.viz.estimation import plot_shrinkage

        result = create_mock_estimation_result()
        fig = plot_shrinkage(result)
        assert fig is not None

    def test_plot_shrinkage_threshold(self):
        """Test shrinkage with custom threshold."""
        from neopkpd.viz.estimation import plot_shrinkage

        result = create_mock_estimation_result()
        fig = plot_shrinkage(result, warning_threshold=25.0)
        assert fig is not None

    def test_plot_eta_distributions(self):
        """Test eta distribution histograms."""
        from neopkpd.viz.estimation import plot_eta_distributions

        result = create_mock_estimation_result()
        fig = plot_eta_distributions(result)
        assert fig is not None

    def test_plot_individual_parameters(self):
        """Test individual parameter plot."""
        from neopkpd.viz.estimation import plot_individual_parameters

        result = create_mock_estimation_result()
        fig = plot_individual_parameters(result)
        assert fig is not None

    def test_plot_ofv_comparison(self):
        """Test OFV comparison plot."""
        from neopkpd.viz.estimation import plot_ofv_comparison

        comparison = create_mock_model_comparison()
        fig = plot_ofv_comparison(comparison)
        assert fig is not None

    def test_plot_sigma_residuals(self):
        """Test sigma residuals plot."""
        from neopkpd.viz.estimation import plot_sigma_residuals

        result = create_mock_estimation_result()
        fig = plot_sigma_residuals(result)
        assert fig is not None


# ============================================================================
# Bootstrap Visualization Tests
# ============================================================================

class TestBootstrapVisualization:
    """Tests for bootstrap visualization functions."""

    def test_plot_bootstrap_distributions(self):
        """Test bootstrap distribution histograms."""
        from neopkpd.viz.bootstrap import plot_bootstrap_distributions

        result = create_mock_bootstrap_result()
        fig = plot_bootstrap_distributions(result)
        assert fig is not None

    def test_plot_bootstrap_distributions_with_ci(self):
        """Test bootstrap distributions with CI."""
        from neopkpd.viz.bootstrap import plot_bootstrap_distributions

        result = create_mock_bootstrap_result()
        fig = plot_bootstrap_distributions(result, show_ci=True, show_original=True)
        assert fig is not None

    def test_plot_bootstrap_ci(self):
        """Test bootstrap CI comparison plot."""
        from neopkpd.viz.bootstrap import plot_bootstrap_ci

        result = create_mock_bootstrap_result()
        fig = plot_bootstrap_ci(result)
        assert fig is not None

    def test_plot_bootstrap_stability(self):
        """Test bootstrap stability plot."""
        from neopkpd.viz.bootstrap import plot_bootstrap_stability

        result = create_mock_bootstrap_result()
        fig = plot_bootstrap_stability(result)
        assert fig is not None

    def test_plot_bootstrap_correlation(self):
        """Test bootstrap correlation heatmap."""
        from neopkpd.viz.bootstrap import plot_bootstrap_correlation

        result = create_mock_bootstrap_result()
        fig = plot_bootstrap_correlation(result)
        assert fig is not None


# ============================================================================
# Sensitivity Visualization Tests
# ============================================================================

class TestSensitivityVisualization:
    """Tests for sensitivity visualization functions."""

    def test_plot_tornado(self):
        """Test tornado diagram."""
        from neopkpd.viz.sensitivity import plot_tornado

        sensitivities = {
            "CL": (85.0, 115.0),
            "V": (98.0, 102.0),
            "Ka": (90.0, 110.0),
            "F": (95.0, 105.0),
        }
        fig = plot_tornado(sensitivities, baseline_value=100.0)
        assert fig is not None

    def test_plot_tornado_sorted(self):
        """Test tornado with sorting."""
        from neopkpd.viz.sensitivity import plot_tornado

        sensitivities = {
            "CL": (85.0, 115.0),
            "V": (98.0, 102.0),
            "Ka": (90.0, 110.0),
        }
        fig = plot_tornado(sensitivities, sort_by_magnitude=True)
        assert fig is not None

    def test_plot_spider(self):
        """Test spider plot."""
        from neopkpd.viz.sensitivity import plot_spider

        curves = {
            "CL": ([0.8, 0.9, 1.0, 1.1, 1.2], [125, 111, 100, 91, 83]),
            "V": ([0.8, 0.9, 1.0, 1.1, 1.2], [100, 100, 100, 100, 100]),
            "Ka": ([0.8, 0.9, 1.0, 1.1, 1.2], [95, 98, 100, 102, 105]),
        }
        fig = plot_spider(curves, normalize=True)
        assert fig is not None

    def test_plot_sensitivity_heatmap(self):
        """Test sensitivity heatmap."""
        from neopkpd.viz.sensitivity import plot_sensitivity_heatmap

        matrix = np.array([
            [0.95, 0.85, 0.10],
            [0.02, 0.30, 0.01],
            [0.15, 0.50, 0.80],
        ])
        fig = plot_sensitivity_heatmap(
            matrix,
            param_names=["CL", "V", "Ka"],
            metric_names=["AUC", "Cmax", "Tmax"]
        )
        assert fig is not None

    def test_plot_waterfall(self):
        """Test waterfall chart."""
        from neopkpd.viz.sensitivity import plot_waterfall

        contributions = {"CL": -15.0, "V": 2.0, "Ka": -5.0, "F": 3.0}
        fig = plot_waterfall(contributions, baseline_value=100.0)
        assert fig is not None


# ============================================================================
# Module Import Tests
# ============================================================================

class TestModuleImports:
    """Tests for module imports."""

    def test_vpc_module_imports(self):
        """Test VPC module imports."""
        from neopkpd.viz.vpc import (
            plot_vpc,
            plot_pcvpc,
            plot_stratified_vpc,
            plot_vpc_with_blq,
            plot_vpc_ci,
        )
        assert callable(plot_vpc)
        assert callable(plot_pcvpc)
        assert callable(plot_stratified_vpc)
        assert callable(plot_vpc_with_blq)
        assert callable(plot_vpc_ci)

    def test_estimation_module_imports(self):
        """Test estimation module imports."""
        from neopkpd.viz.estimation import (
            plot_parameter_estimates,
            plot_omega_matrix,
            plot_convergence,
            plot_shrinkage,
            plot_eta_distributions,
            plot_individual_parameters,
            plot_ofv_comparison,
            plot_correlation_matrix,
            plot_sigma_residuals,
        )
        assert callable(plot_parameter_estimates)
        assert callable(plot_omega_matrix)
        assert callable(plot_convergence)

    def test_bootstrap_module_imports(self):
        """Test bootstrap module imports."""
        from neopkpd.viz.bootstrap import (
            plot_bootstrap_distributions,
            plot_bootstrap_ci,
            plot_bootstrap_stability,
            plot_bootstrap_correlation,
        )
        assert callable(plot_bootstrap_distributions)
        assert callable(plot_bootstrap_ci)

    def test_sensitivity_module_imports(self):
        """Test sensitivity module imports."""
        from neopkpd.viz.sensitivity import (
            plot_tornado,
            plot_spider,
            plot_sensitivity_heatmap,
            plot_waterfall,
        )
        assert callable(plot_tornado)
        assert callable(plot_spider)

    def test_main_viz_module_exports(self):
        """Test main viz module exports new functions."""
        from neopkpd import viz

        # VPC functions
        assert hasattr(viz, 'plot_vpc_detailed')
        assert hasattr(viz, 'plot_pcvpc')
        assert hasattr(viz, 'plot_stratified_vpc')

        # Estimation functions
        assert hasattr(viz, 'plot_parameter_estimates')
        assert hasattr(viz, 'plot_omega_matrix')
        assert hasattr(viz, 'plot_shrinkage')

        # Bootstrap functions
        assert hasattr(viz, 'plot_bootstrap_distributions')
        assert hasattr(viz, 'plot_bootstrap_ci')

        # Sensitivity functions
        assert hasattr(viz, 'plot_tornado')
        assert hasattr(viz, 'plot_spider')
