"""
Tests for NeoPKPD Trial Core Integration

Tests the Julia-connected clinical trial functions including:
- Model-based subject exposure simulation
- Dose escalation algorithms
- Crossover analysis
- Adaptive trial simulation
"""

import pytest
import numpy as np


class TestSubjectExposure:
    """Test model-based subject exposure simulation."""

    @pytest.fixture
    def julia_initialized(self):
        """Initialize Julia."""
        try:
            from neopkpd.bridge import get_julia
            get_julia()
            return True
        except Exception:
            pytest.skip("Julia not available")
            return False

    def test_simulate_exposure_onecomp_iv(self, julia_initialized):
        """Test one-compartment IV bolus exposure simulation."""
        from neopkpd.trial import simulate_subject_exposure

        exposure = simulate_subject_exposure(
            model_kind='OneCompIVBolus',
            dose_events=[{'time': 0.0, 'amount': 100.0, 'rate': 0.0}],
            observation_times=[0, 0.5, 1, 2, 4, 8, 12, 24],
            pk_params={'CL': 10.0, 'V': 50.0}
        )

        assert exposure is not None
        assert len(exposure.times) == 8
        assert len(exposure.concentrations) == 8
        assert 'cmax' in exposure.pk_metrics
        assert exposure.pk_metrics['cmax'] > 0

    def test_simulate_exposure_twocomp_oral(self, julia_initialized):
        """Test two-compartment oral exposure simulation."""
        from neopkpd.trial import simulate_subject_exposure

        exposure = simulate_subject_exposure(
            model_kind='TwoCompOral',
            dose_events=[{'time': 0.0, 'amount': 100.0, 'rate': 0.0}],
            observation_times=[0, 0.5, 1, 2, 4, 8, 12, 24],
            pk_params={'CL': 10.0, 'V1': 50.0, 'Q': 5.0, 'V2': 100.0, 'Ka': 1.0, 'F': 1.0}
        )

        assert exposure is not None
        assert len(exposure.concentrations) > 0
        assert 'cmax' in exposure.pk_metrics
        assert 'tmax' in exposure.pk_metrics
        assert exposure.pk_metrics['tmax'] > 0  # Oral has delayed Tmax

    def test_exposure_pk_metrics(self, julia_initialized):
        """Test PK metrics calculation from exposure."""
        from neopkpd.trial import simulate_subject_exposure

        exposure = simulate_subject_exposure(
            model_kind='OneCompIVBolus',
            dose_events=[{'time': 0.0, 'amount': 100.0, 'rate': 0.0}],
            observation_times=[0, 0.5, 1, 2, 4, 8, 12, 24],
            pk_params={'CL': 10.0, 'V': 50.0}
        )

        metrics = exposure.pk_metrics
        assert 'cmax' in metrics
        assert 'tmax' in metrics
        assert 'auc_0_last' in metrics

        # Cmax should be ~2 for 100mg dose with V=50L
        assert 1.5 < metrics['cmax'] < 2.5


class TestDoseEscalation:
    """Test dose escalation simulation."""

    @pytest.fixture
    def julia_initialized(self):
        """Initialize Julia."""
        try:
            from neopkpd.bridge import get_julia
            get_julia()
            return True
        except Exception:
            pytest.skip("Julia not available")
            return False

    def test_3plus3_basic(self, julia_initialized):
        """Test basic 3+3 dose escalation."""
        from neopkpd.trial import simulate_dose_escalation_3plus3

        result = simulate_dose_escalation_3plus3(
            dose_levels=[10, 25, 50, 100, 200],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            dlt_threshold=50.0,  # Low threshold to trigger DLTs
            seed=42
        )

        assert result is not None
        assert result.design_type == "3+3"
        assert len(result.cohorts) > 0
        assert result.total_subjects > 0
        assert result.completed

    def test_3plus3_finds_mtd(self, julia_initialized):
        """Test 3+3 finds MTD at appropriate dose."""
        from neopkpd.trial import simulate_dose_escalation_3plus3

        # With V=50, Cmax = Dose/V
        # DLT threshold at 3 means MTD should be around dose=150
        result = simulate_dose_escalation_3plus3(
            dose_levels=[25, 50, 100, 150, 200, 300],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            dlt_threshold=3.0,  # Cmax > 3 causes DLT
            seed=123
        )

        # MTD should exist
        assert result.mtd_dose is not None or result.termination_reason in ['too_toxic_at_lowest', 'max_dose_reached']

    def test_mtpi_basic(self, julia_initialized):
        """Test mTPI dose escalation."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_dose_escalation_mtpi

        result = simulate_dose_escalation_mtpi(
            dose_levels=[10, 25, 50, 100, 200],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            dlt_threshold=3.0,
            target_dlt_rate=0.25,
            max_subjects=18,
            seed=42
        )

        assert result is not None
        assert result.design_type == "mTPI"
        assert len(result.cohorts) > 0


class TestCrossoverAnalysis:
    """Test crossover study analysis."""

    def test_analyze_crossover_basic(self):
        """Test basic crossover analysis."""
        from neopkpd.trial import analyze_crossover

        # Simulated 2x2 crossover data
        # AB sequence: Period 1 = Test, Period 2 = Reference
        # BA sequence: Period 1 = Reference, Period 2 = Test
        p1 = [100, 110, 95, 105, 98, 102]  # Period 1 values
        p2 = [98, 108, 92, 100, 100, 98]   # Period 2 values
        seq = ['AB', 'AB', 'AB', 'BA', 'BA', 'BA']

        result = analyze_crossover(p1, p2, seq, log_transform=True)

        assert result is not None
        assert isinstance(result.treatment_effect, float)
        assert isinstance(result.treatment_ci, tuple)
        assert len(result.treatment_ci) == 2
        assert result.within_subject_cv > 0

    def test_crossover_bioequivalence(self):
        """Test bioequivalence assessment from crossover."""
        from neopkpd.trial import analyze_crossover

        # Very similar values should be bioequivalent
        p1 = [100, 102, 98, 101, 99, 100]
        p2 = [99, 101, 97, 100, 100, 99]
        seq = ['AB', 'AB', 'AB', 'BA', 'BA', 'BA']

        result = analyze_crossover(p1, p2, seq, log_transform=True)

        assert result.be_assessment is not None
        assert 'geometric_mean_ratio' in result.be_assessment
        # GMR should be close to 1.0
        assert 0.9 < result.be_assessment['geometric_mean_ratio'] < 1.1

    def test_period_effect(self):
        """Test period effect detection."""
        from neopkpd.trial import test_period_effect

        # No period effect
        p1 = [100, 102, 98, 101, 99, 100]
        p2 = [99, 101, 97, 100, 100, 99]

        result = test_period_effect(p1, p2)

        assert 'effect' in result
        assert 'p_value' in result
        assert result['p_value'] > 0.05  # No significant period effect

    def test_sequence_effect(self):
        """Test sequence effect detection."""
        from neopkpd.trial import test_sequence_effect

        p1 = [100, 102, 98, 101, 99, 100]
        p2 = [99, 101, 97, 100, 100, 99]
        seq = ['AB', 'AB', 'AB', 'BA', 'BA', 'BA']

        result = test_sequence_effect(p1, p2, seq)

        assert 't_statistic' in result
        assert 'p_value' in result

    def test_within_subject_cv(self):
        """Test within-subject CV calculation."""
        from neopkpd.trial import compute_within_subject_cv

        v1 = [100, 105, 95, 102, 98]
        v2 = [102, 103, 97, 100, 99]

        cv = compute_within_subject_cv(v1, v2, log_scale=True)

        assert cv > 0
        assert cv < 50  # Should be reasonable CV


class TestAdaptiveTrial:
    """Test adaptive trial simulation."""

    def test_adaptive_trial_basic(self):
        """Test basic adaptive trial simulation."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_adaptive_trial

        result = simulate_adaptive_trial(
            n_per_arm=50,
            effect_size=0.5,
            sd=1.0,
            interim_fractions=[0.5],
            alpha=0.05,
            seed=42
        )

        assert result is not None
        assert len(result.interim_results) >= 1
        assert result.final_n > 0
        assert isinstance(result.stopped_early, bool)

    def test_adaptive_trial_efficacy_stopping(self):
        """Test that large effect can stop trial early."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_adaptive_trial

        # Very large effect should trigger early stopping
        result = simulate_adaptive_trial(
            n_per_arm=100,
            effect_size=2.0,  # Very large effect
            sd=1.0,
            interim_fractions=[0.5],
            alpha=0.05,
            seed=123
        )

        # Should likely stop early for efficacy
        assert result is not None
        if result.stopped_early:
            assert result.stop_reason == "efficacy"

    def test_adaptive_trial_no_effect(self):
        """Test trial with no true effect."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_adaptive_trial

        result = simulate_adaptive_trial(
            n_per_arm=50,
            effect_size=0.0,  # No effect
            sd=1.0,
            interim_fractions=[0.5],
            alpha=0.05,
            futility_threshold=0.20,
            seed=42
        )

        assert result is not None
        # With no effect, might stop for futility or complete without significance


class TestPKMetrics:
    """Test PK metrics calculation."""

    def test_calculate_pk_metrics_basic(self):
        """Test basic PK metrics calculation."""
        from neopkpd.trial import calculate_pk_metrics

        times = [0, 0.5, 1, 2, 4, 8, 12, 24]
        concs = [0, 1.5, 2.0, 1.8, 1.2, 0.6, 0.3, 0.1]

        metrics = calculate_pk_metrics(times, concs)

        assert 'cmax' in metrics
        assert 'tmax' in metrics
        assert 'auc_0_last' in metrics

        assert metrics['cmax'] == 2.0
        assert metrics['tmax'] == 1.0
        assert metrics['auc_0_last'] > 0

    def test_calculate_pk_metrics_half_life(self):
        """Test half-life calculation."""
        from neopkpd.trial import calculate_pk_metrics

        # Exponential decay with known half-life
        times = [0, 1, 2, 4, 8, 12, 24]
        # t_half = 4h means lambda_z = ln(2)/4 = 0.173
        concs = [2.0, 1.68, 1.41, 1.0, 0.5, 0.25, 0.063]

        metrics = calculate_pk_metrics(times, concs)

        assert 'cmax' in metrics
        assert metrics['cmax'] == 2.0

        # Half-life should be approximately 4 hours
        if 't_half' in metrics:
            assert 2.0 < metrics['t_half'] < 8.0


class TestIntegrationWorkflow:
    """Test complete trial simulation workflows."""

    @pytest.fixture
    def julia_initialized(self):
        """Initialize Julia."""
        try:
            from neopkpd.bridge import get_julia
            get_julia()
            return True
        except Exception:
            pytest.skip("Julia not available")
            return False

    def test_phase1_escalation_workflow(self, julia_initialized):
        """Test complete Phase I dose escalation workflow."""
        from neopkpd.trial import (
            simulate_dose_escalation_3plus3,
            calculate_pk_metrics
        )

        # Run escalation
        result = simulate_dose_escalation_3plus3(
            dose_levels=[10, 25, 50, 100, 200],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            dlt_threshold=3.0,
            seed=42,
            verbose=False
        )

        assert result.completed
        assert len(result.cohorts) > 0

        # Check cohort data
        for cohort in result.cohorts:
            assert cohort.n_subjects > 0
            assert cohort.dose_amount > 0
            assert len(cohort.pk_exposures) == cohort.n_subjects

    def test_be_study_workflow(self, julia_initialized):
        """Test bioequivalence study workflow."""
        from neopkpd.trial import (
            simulate_subject_exposure,
            analyze_crossover,
            crossover_2x2
        )

        # Simulate crossover BE study
        np.random.seed(42)
        n_subjects = 24

        # Reference formulation
        ref_exposures = []
        test_exposures = []

        for i in range(n_subjects):
            # Reference
            ref_exp = simulate_subject_exposure(
                model_kind='OneCompOral',
                dose_events=[{'time': 0.0, 'amount': 100.0, 'rate': 0.0}],
                observation_times=[0, 0.5, 1, 2, 4, 8, 12, 24],
                pk_params={'CL': 10.0, 'V': 50.0, 'Ka': 1.0, 'F': 1.0},
                seed=42 + i
            )
            ref_exposures.append(ref_exp.pk_metrics.get('cmax', 1.0) * (1 + np.random.normal(0, 0.2)))

            # Test (slightly different bioavailability)
            test_exp = simulate_subject_exposure(
                model_kind='OneCompOral',
                dose_events=[{'time': 0.0, 'amount': 100.0, 'rate': 0.0}],
                observation_times=[0, 0.5, 1, 2, 4, 8, 12, 24],
                pk_params={'CL': 10.0, 'V': 50.0, 'Ka': 1.2, 'F': 0.95},
                seed=42 + i + 1000
            )
            test_exposures.append(test_exp.pk_metrics.get('cmax', 1.0) * (1 + np.random.normal(0, 0.2)))

        # Create crossover structure
        half = n_subjects // 2
        sequences = ['AB'] * half + ['BA'] * half

        # Period 1: AB gets Test, BA gets Reference
        # Period 2: AB gets Reference, BA gets Test
        p1 = test_exposures[:half] + ref_exposures[half:]
        p2 = ref_exposures[:half] + test_exposures[half:]

        # Analyze
        analysis = analyze_crossover(p1, p2, sequences, log_transform=True)

        assert analysis.be_assessment is not None
        assert 'geometric_mean_ratio' in analysis.be_assessment
        assert analysis.within_subject_cv > 0


class TestCRMDoseEscalation:
    """Test CRM dose escalation simulation."""

    @pytest.fixture
    def julia_initialized(self):
        """Initialize Julia."""
        try:
            from neopkpd.bridge import get_julia
            get_julia()
            return True
        except Exception:
            pytest.skip("Julia not available")
            return False

    def test_crm_basic(self, julia_initialized):
        """Test basic CRM dose escalation."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_dose_escalation_crm

        result = simulate_dose_escalation_crm(
            dose_levels=[10, 25, 50, 100, 200],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            dlt_threshold=3.0,
            target_dlt_rate=0.25,
            max_subjects=15,
            seed=42
        )

        assert result is not None
        assert result.design_type == "CRM"
        assert len(result.cohorts) > 0
        assert result.total_subjects <= 15
        assert result.mtd_dose is not None

    def test_crm_with_skeleton(self, julia_initialized):
        """Test CRM with custom prior skeleton."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_dose_escalation_crm

        skeleton = [0.05, 0.10, 0.20, 0.35, 0.50]
        result = simulate_dose_escalation_crm(
            dose_levels=[10, 25, 50, 100, 200],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            dlt_threshold=3.0,
            target_dlt_rate=0.25,
            skeleton=skeleton,
            max_subjects=12,
            seed=123
        )

        assert result is not None
        assert result.completed
        assert result.mtd_level is not None


class TestModelConnectedTrial:
    """Test model-connected clinical trial simulation."""

    @pytest.fixture
    def julia_initialized(self):
        """Initialize Julia."""
        try:
            from neopkpd.bridge import get_julia
            get_julia()
            return True
        except Exception:
            pytest.skip("Julia not available")
            return False

    def test_trial_with_model_basic(self, julia_initialized):
        """Test basic model-connected trial simulation."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_trial_with_model

        result = simulate_trial_with_model(
            trial_name="Test PK Trial",
            arms=[
                {'name': 'Placebo', 'dose': 0.0, 'n_subjects': 5, 'placebo': True},
                {'name': 'Low Dose', 'dose': 50.0, 'n_subjects': 5},
                {'name': 'High Dose', 'dose': 100.0, 'n_subjects': 5},
            ],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            observation_times=[0, 1, 2, 4, 8, 12, 24],
            seed=42
        )

        assert result is not None
        assert result.trial_name == "Test PK Trial"
        assert len(result.arms) == 3
        assert 'Placebo' in result.arms
        assert 'Low Dose' in result.arms
        assert 'High Dose' in result.arms

        # Check arm results
        placebo = result.arms['Placebo']
        assert placebo.n_enrolled == 5
        assert len(placebo.subjects) == 5

        high_dose = result.arms['High Dose']
        assert high_dose.n_enrolled == 5
        # High dose should have higher exposure than placebo
        if high_dose.pk_summary.get('cmax'):
            assert high_dose.pk_summary['cmax']['mean'] > 0

    def test_trial_with_model_iiv(self, julia_initialized):
        """Test model-connected trial with inter-individual variability."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_trial_with_model

        # 2x2 omega matrix for CL and V
        omega = [[0.09, 0.0], [0.0, 0.04]]  # 30% CV for CL, 20% for V

        result = simulate_trial_with_model(
            trial_name="IIV PK Trial",
            arms=[
                {'name': 'Active', 'dose': 100.0, 'n_subjects': 10},
            ],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            observation_times=[0, 1, 2, 4, 8, 24],
            omega=omega,
            seed=42
        )

        assert result is not None
        active_arm = result.arms['Active']
        assert active_arm.n_enrolled == 10

        # With IIV, should see variability in Cmax
        cmax_values = [s.pk_metrics.get('cmax', 0) for s in active_arm.subjects]
        assert len(set(cmax_values)) > 1  # Should have different values

    def test_trial_with_dropout(self, julia_initialized):
        """Test model-connected trial with dropout."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_trial_with_model

        result = simulate_trial_with_model(
            trial_name="Dropout Trial",
            arms=[
                {'name': 'Active', 'dose': 100.0, 'n_subjects': 20},
            ],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            observation_times=[0, 1, 4, 8, 24],
            duration_days=28.0,
            dropout_rate=0.05,  # 5% daily dropout
            seed=42
        )

        assert result is not None
        active_arm = result.arms['Active']
        # With 5% daily dropout over 28 days, expect some dropouts
        # but not all dropouts (probabilistic)
        assert active_arm.n_enrolled == 20
        assert active_arm.n_completed <= active_arm.n_enrolled

    def test_trial_comparisons(self, julia_initialized):
        """Test statistical comparisons between arms."""
        pytest.importorskip("scipy")
        from neopkpd.trial import simulate_trial_with_model

        result = simulate_trial_with_model(
            trial_name="Comparison Trial",
            arms=[
                {'name': 'Placebo', 'dose': 0.0, 'n_subjects': 10, 'placebo': True},
                {'name': 'Active', 'dose': 100.0, 'n_subjects': 10},
            ],
            model_kind='OneCompIVBolus',
            pk_params={'CL': 10.0, 'V': 50.0},
            observation_times=[0, 1, 4, 8, 24],
            seed=42
        )

        assert result is not None
        # Should have comparison between Active and Placebo
        assert len(result.comparisons) > 0

        comparison_key = 'Active_vs_Placebo'
        assert comparison_key in result.comparisons

        comparison = result.comparisons[comparison_key]
        assert 'difference' in comparison
        assert 'p_value' in comparison
        # Active should have significantly higher Cmax than placebo
        assert comparison['difference'] > 0
