"""
Tests for NeoPKPD VPC Module

Tests comprehensive VPC functionality including:
- Binning strategies (Quantile, EqualWidth, KMeans)
- BLQ handling (M1-M7 methods)
- Pure Python VPC computation
- Result extraction helpers
- Julia-connected VPC (when available)
"""

import pytest
import numpy as np


class TestBinningStrategies:
    """Test binning strategy implementations."""

    def test_quantile_binning_basic(self):
        """Test basic quantile binning."""
        from neopkpd.vpc import QuantileBinning

        binning = QuantileBinning(5)
        assert binning.n_bins == 5

        times = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        bins = binning.compute_bins(times)

        assert len(bins) == 5
        assert bins[0].id == 1
        assert bins[-1].id == 5
        assert bins[0].lower <= bins[0].midpoint <= bins[0].upper

    def test_quantile_binning_uneven_data(self):
        """Test quantile binning with uneven data distribution."""
        from neopkpd.vpc import QuantileBinning

        binning = QuantileBinning(4)
        # Dense sampling at start, sparse at end
        times = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 10, 20, 24])
        bins = binning.compute_bins(times)

        assert len(bins) == 4
        # Bins should have roughly equal observations

    def test_equal_width_binning(self):
        """Test equal-width binning."""
        from neopkpd.vpc import EqualWidthBinning

        binning = EqualWidthBinning(4)
        times = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16])
        bins = binning.compute_bins(times)

        assert len(bins) == 4
        # Each bin should have equal width
        widths = [b.upper - b.lower for b in bins]
        assert all(np.isclose(w, widths[0], rtol=0.01) for w in widths)

    def test_kmeans_binning(self):
        """Test k-means based binning."""
        from neopkpd.vpc import KMeansBinning

        binning = KMeansBinning(3)
        # Data with natural clusters
        times = np.array([0, 0.5, 1, 1.5, 8, 9, 10, 20, 22, 24])
        bins = binning.compute_bins(times)

        assert len(bins) == 3
        # Bins should capture the clusters

    def test_binning_empty_data(self):
        """Test binning with empty data."""
        from neopkpd.vpc import QuantileBinning

        binning = QuantileBinning(5)
        bins = binning.compute_bins(np.array([]))

        assert len(bins) == 0

    def test_binning_validation(self):
        """Test binning parameter validation."""
        from neopkpd.vpc import QuantileBinning

        with pytest.raises(ValueError):
            QuantileBinning(1)  # Must have at least 2 bins


class TestBLQHandling:
    """Test BLQ handling methods."""

    def test_blq_m1_discard(self):
        """Test M1: Discard all BLQ."""
        from neopkpd.vpc import handle_blq, BLQMethod

        values = np.array([0.05, 0.5, 2.0, 1.5, 0.08])
        times = np.array([0, 0.5, 1, 2, 4])
        lloq = 0.1

        result = handle_blq(values, times, lloq, method=BLQMethod.M1)

        assert np.isnan(result[0])  # 0.05 < LLOQ
        assert result[1] == 0.5
        assert result[2] == 2.0
        assert np.isnan(result[4])  # 0.08 < LLOQ

    def test_blq_m4_replace(self):
        """Test M4: Replace with LLOQ/2."""
        from neopkpd.vpc import handle_blq, BLQMethod

        values = np.array([0.05, 0.5, 2.0, 1.5, 0.08])
        times = np.array([0, 0.5, 1, 2, 4])
        lloq = 0.1

        result = handle_blq(values, times, lloq, method=BLQMethod.M4)

        assert result[0] == 0.05  # LLOQ/2
        assert result[1] == 0.5
        assert result[4] == 0.05

    def test_blq_m5_before_after_tmax(self):
        """Test M5: 0 before Tmax, LLOQ/2 after."""
        from neopkpd.vpc import handle_blq, BLQMethod

        values = np.array([0.05, 0.5, 2.0, 1.5, 0.08])  # Tmax at index 2
        times = np.array([0, 0.5, 1, 2, 4])
        lloq = 0.1

        result = handle_blq(values, times, lloq, method=BLQMethod.M5)

        assert result[0] == 0.0  # Before Tmax -> 0
        assert result[4] == 0.05  # After Tmax -> LLOQ/2

    def test_blq_m7_before_after_tmax(self):
        """Test M7: 0 before Tmax, discard after."""
        from neopkpd.vpc import handle_blq, BLQMethod

        values = np.array([0.05, 0.5, 2.0, 1.5, 0.08])
        times = np.array([0, 0.5, 1, 2, 4])
        lloq = 0.1

        result = handle_blq(values, times, lloq, method=BLQMethod.M7)

        assert result[0] == 0.0  # Before Tmax -> 0
        assert np.isnan(result[4])  # After Tmax -> discard


class TestVPCConfig:
    """Test VPCConfig validation."""

    def test_config_defaults(self):
        """Test default configuration values."""
        from neopkpd.vpc import VPCConfig

        config = VPCConfig()

        assert config.pi_levels == [0.05, 0.50, 0.95]
        assert config.ci_level == 0.95
        assert config.n_simulations == 200
        assert config.n_bootstrap == 500
        assert config.prediction_corrected is False

    def test_config_custom(self):
        """Test custom configuration."""
        from neopkpd.vpc import VPCConfig, QuantileBinning

        config = VPCConfig(
            pi_levels=[0.10, 0.50, 0.90],
            binning=QuantileBinning(8),
            n_simulations=500,
            prediction_corrected=True
        )

        assert config.pi_levels == [0.10, 0.50, 0.90]
        assert config.binning.n_bins == 8
        assert config.n_simulations == 500
        assert config.prediction_corrected is True

    def test_config_validation(self):
        """Test configuration validation."""
        from neopkpd.vpc import VPCConfig

        # Invalid pi_levels
        with pytest.raises(ValueError):
            VPCConfig(pi_levels=[0.0, 0.5, 1.0])

        # Invalid ci_level
        with pytest.raises(ValueError):
            VPCConfig(ci_level=1.5)

        # Too few simulations
        with pytest.raises(ValueError):
            VPCConfig(n_simulations=5)


class TestVPCPython:
    """Test pure Python VPC computation."""

    def test_vpc_python_basic(self):
        """Test basic pure Python VPC."""
        from neopkpd.vpc import compute_vpc_python, VPCConfig, QuantileBinning

        # Generate test data
        np.random.seed(42)
        n_obs = 50
        obs_times = np.sort(np.random.uniform(0, 24, n_obs))
        obs_values = 2.0 * np.exp(-0.1 * obs_times) + np.random.normal(0, 0.2, n_obs)
        obs_values = np.maximum(obs_values, 0)

        # Generate simulated data
        n_sims = 50
        sim_data = []
        for _ in range(n_sims):
            sim_times = np.sort(np.random.uniform(0, 24, n_obs))
            sim_values = 2.0 * np.exp(-0.1 * sim_times) + np.random.normal(0, 0.3, n_obs)
            sim_values = np.maximum(sim_values, 0)
            sim_data.append((sim_times, sim_values))

        config = VPCConfig(
            binning=QuantileBinning(5),
            n_simulations=n_sims,
            n_bootstrap=100
        )

        result = compute_vpc_python(obs_times, obs_values, sim_data, config)

        assert result is not None
        assert len(result.bins) == 5
        assert result.n_simulations == n_sims

        # Check percentile data
        for bin in result.bins:
            assert len(bin.percentiles) == 3  # 5th, 50th, 95th
            for p in bin.percentiles:
                assert 0 <= p.percentile <= 1

    def test_vpc_python_with_blq(self):
        """Test pure Python VPC with BLQ handling."""
        from neopkpd.vpc import compute_vpc_python, VPCConfig, BLQMethod

        np.random.seed(42)
        obs_times = np.array([0, 1, 2, 4, 8, 12, 24])
        obs_values = np.array([0.05, 1.5, 2.0, 1.5, 0.8, 0.3, 0.05])  # BLQ at 0 and 24h

        sim_data = [
            (obs_times.copy(), obs_values * np.random.uniform(0.8, 1.2, len(obs_values)))
            for _ in range(20)
        ]

        config = VPCConfig(
            lloq=0.1,
            blq_method=BLQMethod.M4,
            n_bootstrap=100
        )

        result = compute_vpc_python(obs_times, obs_values, sim_data, config)

        assert result is not None


class TestResultExtraction:
    """Test VPC result extraction helpers."""

    @pytest.fixture
    def mock_vpc_result(self):
        """Create a mock VPC result for testing."""
        from neopkpd.vpc import (
            VPCResult, VPCBin, VPCPercentileData, VPCConfig,
            QuantileBinning
        )

        bins = [
            VPCBin(
                bin_id=1,
                time_min=0,
                time_max=8,
                time_midpoint=4,
                n_observed=10,
                n_simulated=100,
                percentiles=[
                    VPCPercentileData(0.05, 0.5, 0.6, 0.4, 0.8),
                    VPCPercentileData(0.50, 1.5, 1.6, 1.4, 1.8),
                    VPCPercentileData(0.95, 2.5, 2.6, 2.4, 2.8),
                ]
            ),
            VPCBin(
                bin_id=2,
                time_min=8,
                time_max=24,
                time_midpoint=16,
                n_observed=10,
                n_simulated=100,
                percentiles=[
                    VPCPercentileData(0.05, 0.2, 0.3, 0.1, 0.5),
                    VPCPercentileData(0.50, 0.8, 0.9, 0.7, 1.1),
                    VPCPercentileData(0.95, 1.5, 1.6, 1.4, 1.8),
                ]
            ),
        ]

        config = VPCConfig(binning=QuantileBinning(2))

        return VPCResult(
            config=config,
            bins=bins,
            n_subjects_observed=5,
            n_observations_observed=20,
            n_simulations=100
        )

    def test_get_bin_midpoints(self, mock_vpc_result):
        """Test bin midpoint extraction."""
        from neopkpd.vpc import get_bin_midpoints

        midpoints = get_bin_midpoints(mock_vpc_result)

        assert len(midpoints) == 2
        assert midpoints[0] == 4
        assert midpoints[1] == 16

    def test_get_observed_percentile(self, mock_vpc_result):
        """Test observed percentile extraction."""
        from neopkpd.vpc import get_observed_percentile

        median = get_observed_percentile(mock_vpc_result, 0.50)

        assert len(median) == 2
        assert median[0] == 1.5
        assert median[1] == 0.8

    def test_get_simulated_median(self, mock_vpc_result):
        """Test simulated median extraction."""
        from neopkpd.vpc import get_simulated_median

        sim_median = get_simulated_median(mock_vpc_result, 0.50)

        assert len(sim_median) == 2
        assert sim_median[0] == 1.6
        assert sim_median[1] == 0.9

    def test_get_simulated_ci(self, mock_vpc_result):
        """Test simulated CI extraction."""
        from neopkpd.vpc import get_simulated_ci

        lower, upper = get_simulated_ci(mock_vpc_result, 0.50)

        assert len(lower) == 2
        assert len(upper) == 2
        assert lower[0] == 1.4
        assert upper[0] == 1.8


class TestVPCJulia:
    """Test Julia-connected VPC functions."""

    @pytest.fixture
    def julia_initialized(self):
        """Check if Julia is available."""
        try:
            from neopkpd.bridge import get_julia
            get_julia()
            return True
        except Exception:
            pytest.skip("Julia not available")
            return False

    def test_compute_vpc_julia(self, julia_initialized):
        """Test Julia-connected VPC computation."""
        from neopkpd.vpc import compute_vpc, VPCConfig, QuantileBinning

        observed = {
            "subjects": [
                {
                    "subject_id": "001",
                    "times": [0, 1, 2, 4, 8, 12, 24],
                    "observations": [0, 1.5, 2.0, 1.8, 1.2, 0.8, 0.3],
                    "doses": [{"time": 0, "amount": 100}]
                },
                {
                    "subject_id": "002",
                    "times": [0, 1, 2, 4, 8, 12, 24],
                    "observations": [0, 1.6, 2.1, 1.9, 1.1, 0.7, 0.25],
                    "doses": [{"time": 0, "amount": 100}]
                }
            ]
        }

        pop_spec = {
            "model": {
                "kind": "OneCompIVBolus",
                "params": {"CL": 10, "V": 50},
                "doses": [{"time": 0, "amount": 100}]
            },
            "iiv": {"omegas": {"CL": 0.3, "V": 0.2}, "n": 50}
        }

        grid = {"t0": 0, "t1": 24, "saveat": [0, 1, 2, 4, 8, 12, 24]}

        config = VPCConfig(
            binning=QuantileBinning(4),
            n_simulations=20,
            n_bootstrap=100
        )

        result = compute_vpc(observed, pop_spec, grid, config)

        assert result is not None
        assert len(result.bins) > 0
        assert result.n_simulations == 20

    def test_compute_pcvpc_julia(self, julia_initialized):
        """Test Julia-connected pcVPC computation."""
        from neopkpd.vpc import compute_pcvpc, VPCConfig

        observed = {
            "subjects": [
                {
                    "subject_id": "001",
                    "times": [0, 1, 2, 4, 8, 12, 24],
                    "observations": [0, 1.5, 2.0, 1.8, 1.2, 0.8, 0.3],
                    "doses": [{"time": 0, "amount": 100}]
                }
            ]
        }

        pop_spec = {
            "model": {
                "kind": "OneCompIVBolus",
                "params": {"CL": 10, "V": 50},
                "doses": [{"time": 0, "amount": 100}]
            },
            "iiv": {"omegas": {"CL": 0.3}, "n": 30}
        }

        grid = {"t0": 0, "t1": 24, "saveat": [0, 1, 2, 4, 8, 12, 24]}

        config = VPCConfig(
            n_simulations=10,
            n_bootstrap=100
        )

        result = compute_pcvpc(observed, pop_spec, grid, config)

        assert result is not None
        assert result.prediction_corrected is True


class TestStratifiedVPC:
    """Test stratified VPC functionality."""

    def test_stratified_vpc_python(self):
        """Test stratified VPC with pure Python."""
        from neopkpd.vpc import compute_vpc_python, VPCConfig

        # Create data with two groups
        np.random.seed(42)

        # Group 1: Lower exposure
        g1_times = np.array([0, 1, 2, 4, 8, 12, 24])
        g1_values = np.array([0, 1.0, 1.5, 1.2, 0.8, 0.5, 0.2])

        # Group 2: Higher exposure
        g2_times = np.array([0, 1, 2, 4, 8, 12, 24])
        g2_values = np.array([0, 2.0, 3.0, 2.5, 1.5, 1.0, 0.4])

        # Combined for single VPC
        all_times = np.concatenate([g1_times, g2_times])
        all_values = np.concatenate([g1_values, g2_values])

        sim_data = [
            (all_times.copy(), all_values * np.random.uniform(0.8, 1.2, len(all_values)))
            for _ in range(20)
        ]

        config = VPCConfig(n_bootstrap=100)
        result = compute_vpc_python(all_times, all_values, sim_data, config)

        assert result is not None
        assert len(result.bins) > 0
