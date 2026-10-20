"""
Tests for NeoPKPD Trial Module.

These tests verify the trial simulation, study designs, dosing regimens,
virtual population generation, and statistical analysis functions.
"""

import pytest
import math

from neopkpd.trial import (
    # Designs
    parallel_design,
    crossover_2x2,
    crossover_3x3,
    williams_design,
    dose_escalation_3plus3,
    dose_escalation_mtpi,
    dose_escalation_crm,
    adaptive_design,
    bioequivalence_design,
    get_design_description,
    ParallelDesign,
    CrossoverDesign,
    DoseEscalationDesign,
    BioequivalenceDesign,
    AdaptiveDesign,
    # Regimens
    dosing_qd,
    dosing_bid,
    dosing_tid,
    dosing_qid,
    dosing_custom,
    titration_regimen,
    dose_event_times,
    total_regimen_duration,
    DosingRegimen,
    TitrationRegimen,
    # Population
    generate_virtual_population,
    default_demographic_spec,
    healthy_volunteer_spec,
    patient_population_spec,
    summarize_population,
    DemographicSpec,
    DiseaseSpec,
    VirtualPopulationSpec,
    VirtualSubject,
    # Simulation
    simulate_trial,
    simulate_trial_replicates,
    simulate_dropout,
    apply_compliance,
    TrialSpec,
    TrialResult,
    TreatmentArm,
    DropoutSpec,
    ComplianceSpec,
    # Analysis
    estimate_power_analytical,
    estimate_sample_size,
    alpha_spending_function,
    incremental_alpha,
    compare_arms,
    responder_analysis,
    bioequivalence_90ci,
    assess_bioequivalence,
    PowerResult,
    SampleSizeResult,
    ComparisonResult,
    ResponderResult,
)


# ============================================================================
# Study Designs Tests
# ============================================================================

class TestStudyDesigns:
    """Tests for study design creation."""

    def test_parallel_design_basic(self):
        """Test basic parallel design creation."""
        design = parallel_design(2)

        assert isinstance(design, ParallelDesign)
        assert design.n_arms == 2
        assert len(design.randomization_ratio) == 2
        assert sum(design.randomization_ratio) == pytest.approx(1.0)

    def test_parallel_design_custom_ratio(self):
        """Test parallel design with custom randomization ratio."""
        design = parallel_design(3, randomization_ratio=[0.5, 0.25, 0.25])

        assert design.n_arms == 3
        assert design.randomization_ratio == [0.5, 0.25, 0.25]

    def test_crossover_2x2(self):
        """Test 2x2 crossover design."""
        design = crossover_2x2(washout_duration=14.0)

        assert isinstance(design, CrossoverDesign)
        assert design.n_periods == 2
        assert design.n_sequences == 2
        assert design.washout_duration == 14.0
        assert design.sequence_assignments == [[1, 2], [2, 1]]

    def test_crossover_3x3(self):
        """Test 3x3 crossover design (Latin square)."""
        design = crossover_3x3()

        assert design.n_periods == 3
        assert design.n_sequences == 3
        assert len(design.sequence_assignments) == 3

    def test_williams_design(self):
        """Test Williams design for different treatment numbers."""
        design2 = williams_design(2)
        assert design2.n_periods == 2
        assert design2.n_sequences == 2

        design3 = williams_design(3)
        assert design3.n_periods == 3
        assert design3.n_sequences == 6

        design4 = williams_design(4)
        assert design4.n_periods == 4
        assert design4.n_sequences == 8

    def test_williams_design_invalid(self):
        """Test Williams design with invalid treatment count."""
        with pytest.raises(ValueError):
            williams_design(5)

    def test_dose_escalation_3plus3(self):
        """Test 3+3 dose escalation design."""
        dose_levels = [10.0, 25.0, 50.0, 100.0, 200.0]
        design = dose_escalation_3plus3(dose_levels)

        assert isinstance(design, DoseEscalationDesign)
        assert design.dose_levels == dose_levels
        assert design.starting_dose == 10.0
        assert design.escalation_rule == "3+3"
        assert design.cohort_size == 3

    def test_dose_escalation_mtpi(self):
        """Test mTPI dose escalation design."""
        dose_levels = [10.0, 25.0, 50.0, 100.0]
        design = dose_escalation_mtpi(dose_levels, target_dlt_rate=0.30)

        assert design.escalation_rule == "mTPI"
        assert design.target_dlt_rate == 0.30

    def test_dose_escalation_crm(self):
        """Test CRM dose escalation design."""
        dose_levels = [10.0, 25.0, 50.0, 100.0, 200.0]
        design = dose_escalation_crm(dose_levels)

        assert design.escalation_rule == "CRM"
        assert design.cohort_size == 1  # CRM typically uses cohort size of 1

    def test_adaptive_design(self):
        """Test adaptive design creation."""
        base = parallel_design(2)
        design = adaptive_design(base, interim_analyses=[0.5, 0.75])

        assert isinstance(design, AdaptiveDesign)
        assert design.interim_analyses == [0.5, 0.75]
        assert design.alpha_spending == "obrien_fleming"

    def test_bioequivalence_design(self):
        """Test bioequivalence design."""
        design = bioequivalence_design(regulatory_guidance="fda")

        assert isinstance(design, BioequivalenceDesign)
        assert design.bioequivalence_limits == (0.80, 1.25)
        assert "cmax" in design.parameters
        assert "auc_0_inf" in design.parameters

    def test_get_design_description(self):
        """Test design description generation."""
        parallel = parallel_design(3)
        assert "3-arm parallel" in get_design_description(parallel)

        crossover = crossover_2x2()
        assert "crossover" in get_design_description(crossover)

        escalation = dose_escalation_3plus3([10, 25, 50, 100])
        assert "escalation" in get_design_description(escalation).lower()


# ============================================================================
# Dosing Regimens Tests
# ============================================================================

class TestDosingRegimens:
    """Tests for dosing regimen creation."""

    def test_dosing_qd(self):
        """Test once-daily dosing regimen."""
        regimen = dosing_qd(100.0, 7)

        assert isinstance(regimen, DosingRegimen)
        assert regimen.dose_amount == 100.0
        assert regimen.duration_days == 7
        assert len(regimen.dose_times) == 1

    def test_dosing_qd_with_loading(self):
        """Test QD with loading dose."""
        regimen = dosing_qd(100.0, 7, loading_dose=200.0)

        assert regimen.loading_dose == 200.0

    def test_dosing_bid(self):
        """Test twice-daily dosing regimen."""
        regimen = dosing_bid(50.0, 14)

        assert regimen.dose_amount == 50.0
        assert len(regimen.dose_times) == 2

    def test_dosing_tid(self):
        """Test three-times-daily dosing regimen."""
        regimen = dosing_tid(25.0, 7)

        assert len(regimen.dose_times) == 3

    def test_dosing_qid(self):
        """Test four-times-daily dosing regimen."""
        regimen = dosing_qid(20.0, 5)

        assert len(regimen.dose_times) == 4

    def test_dosing_custom(self):
        """Test custom dosing regimen."""
        custom_times = [6.0, 22.0]
        regimen = dosing_custom(75.0, 21, custom_times)

        assert regimen.dose_times == custom_times

    def test_titration_regimen(self):
        """Test titration regimen creation."""
        regimen = titration_regimen(25.0, 100.0, 4, 7)

        assert isinstance(regimen, TitrationRegimen)
        assert len(regimen.steps) == 4
        assert regimen.steps[0].dose == 25.0
        assert regimen.steps[-1].dose == 100.0

    def test_titration_regimen_with_maintenance(self):
        """Test titration regimen with maintenance phase."""
        regimen = titration_regimen(25.0, 100.0, 4, 7, maintenance_days=28)

        assert len(regimen.steps) == 5
        assert regimen.steps[-1].duration_days == 28

    def test_titration_minimum_steps(self):
        """Test titration requires at least 2 steps."""
        with pytest.raises(ValueError):
            titration_regimen(25.0, 100.0, 1, 7)

    def test_dose_event_times_qd(self):
        """Test dose event times for QD regimen."""
        regimen = dosing_qd(100.0, 3)
        times = dose_event_times(regimen)

        assert len(times) == 3
        assert times[1] - times[0] == 24.0

    def test_dose_event_times_bid(self):
        """Test dose event times for BID regimen."""
        regimen = dosing_bid(50.0, 2)
        times = dose_event_times(regimen)

        assert len(times) == 4

    def test_total_regimen_duration(self):
        """Test total duration calculation."""
        qd = dosing_qd(100.0, 7)
        assert total_regimen_duration(qd) == 7

        titration = titration_regimen(25.0, 100.0, 4, 7)
        assert total_regimen_duration(titration) == 28


# ============================================================================
# Virtual Population Tests
# ============================================================================

class TestVirtualPopulation:
    """Tests for virtual population generation."""

    def test_generate_basic_population(self):
        """Test basic population generation."""
        pop = generate_virtual_population(100, seed=42)

        assert len(pop) == 100
        assert all(isinstance(s, VirtualSubject) for s in pop)

    def test_population_reproducibility(self):
        """Test population generation reproducibility with seed."""
        pop1 = generate_virtual_population(50, seed=123)
        pop2 = generate_virtual_population(50, seed=123)

        for s1, s2 in zip(pop1, pop2):
            assert s1.age == s2.age
            assert s1.weight == s2.weight
            assert s1.sex == s2.sex

    def test_demographic_ranges(self):
        """Test that demographics are within specified ranges."""
        spec = DemographicSpec(
            age_mean=50.0,
            age_sd=10.0,
            age_range=(40.0, 60.0),
            weight_mean=75.0,
            weight_sd=10.0,
            weight_range=(60.0, 90.0),
        )
        pop = generate_virtual_population(100, demographics=spec, seed=42)

        for s in pop:
            assert 40.0 <= s.age <= 60.0
            assert 60.0 <= s.weight <= 90.0

    def test_sex_distribution(self):
        """Test sex distribution in population."""
        spec = DemographicSpec(female_proportion=0.6)
        pop = generate_virtual_population(1000, demographics=spec, seed=42)

        n_female = sum(1 for s in pop if s.sex == "female")
        prop_female = n_female / len(pop)

        assert 0.5 < prop_female < 0.7

    def test_default_demographic_spec(self):
        """Test default demographic specification."""
        spec = default_demographic_spec()

        assert spec.age_mean == 35.0
        assert spec.female_proportion == 0.5

    def test_healthy_volunteer_spec(self):
        """Test healthy volunteer specification."""
        spec = healthy_volunteer_spec()

        assert spec.age_mean == 30.0
        assert spec.age_range[1] <= 45.0

    def test_patient_population_spec(self):
        """Test patient population specifications."""
        demo, dis = patient_population_spec("diabetes")

        assert demo.age_mean == 55.0
        assert dis.name == "diabetes"

        demo, dis = patient_population_spec("renal")
        assert "renal" in dis.name.lower()

        demo, dis = patient_population_spec("oncology")
        assert "cancer" in dis.name.lower() or "oncology" in dis.name.lower()

    def test_population_with_disease(self):
        """Test population generation with disease specification."""
        demo, dis = patient_population_spec("diabetes")
        pop = generate_virtual_population(100, demographics=demo, disease=dis, seed=42)

        assert all(s.disease_severity is not None for s in pop)
        assert all(s.baseline_biomarker is not None for s in pop)

    def test_summarize_population(self):
        """Test population summary statistics."""
        pop = generate_virtual_population(100, seed=42)
        summary = summarize_population(pop)

        assert summary["n"] == 100
        assert "age_mean" in summary
        assert "weight_mean" in summary
        assert "female_proportion" in summary
        assert "race_distribution" in summary

    def test_virtual_subject_to_dict(self):
        """Test VirtualSubject to dict conversion."""
        subject = VirtualSubject(
            id=1, age=45.0, weight=75.0,
            sex="male", race="caucasian"
        )
        d = subject.to_dict()

        assert d["id"] == 1
        assert d["age"] == 45.0
        assert d["sex"] == "male"


# ============================================================================
# Simulation Tests
# ============================================================================

class TestSimulation:
    """Tests for trial simulation."""

    def test_simulate_dropout(self):
        """Test dropout simulation."""
        spec = DropoutSpec(random_rate_per_day=0.05)
        dropout_days = simulate_dropout(100, 28.0, spec, seed=42)

        assert len(dropout_days) == 100
        n_dropouts = sum(1 for d in dropout_days if d is not None)
        assert n_dropouts > 0  # Should have some dropouts

    def test_no_dropout(self):
        """Test simulation without dropout."""
        dropout_days = simulate_dropout(100, 28.0, None)

        assert all(d is None for d in dropout_days)

    def test_apply_compliance_random(self):
        """Test random compliance application."""
        times = list(range(0, 168, 24))  # 7 days
        amounts = [100.0] * 7
        spec = ComplianceSpec(mean_compliance=0.80, pattern="random")

        actual = apply_compliance(times, amounts, spec, seed=42)

        assert len(actual) == 7
        assert any(d == 0.0 for d in actual)  # Some missed doses

    def test_apply_compliance_decay(self):
        """Test decaying compliance."""
        times = list(range(0, 168, 24))
        amounts = [100.0] * 7
        spec = ComplianceSpec(mean_compliance=0.90, pattern="decay")

        actual = apply_compliance(times, amounts, spec, seed=42)
        assert len(actual) == 7

    def test_simulate_trial_basic(self):
        """Test basic trial simulation."""
        spec = TrialSpec(
            name="Test Trial",
            design=parallel_design(2),
            arms=[
                TreatmentArm("Placebo", dosing_qd(0.0, 28), 25, placebo=True),
                TreatmentArm("Active", dosing_qd(100.0, 28), 25),
            ],
            duration_days=28,
            seed=42,
        )

        result = simulate_trial(spec)

        assert isinstance(result, TrialResult)
        assert result.trial_name == "Test Trial"
        assert len(result.arms) == 2
        assert "Placebo" in result.arms
        assert "Active" in result.arms

    def test_simulate_trial_with_dropout(self):
        """Test trial simulation with dropout."""
        spec = TrialSpec(
            name="Dropout Trial",
            design=parallel_design(2),
            arms=[
                TreatmentArm("Placebo", dosing_qd(0.0, 28), 50, placebo=True),
                TreatmentArm("Active", dosing_qd(100.0, 28), 50),
            ],
            duration_days=28,
            dropout=DropoutSpec(random_rate_per_day=0.02),
            seed=42,
        )

        result = simulate_trial(spec)

        assert result.overall_completion_rate < 1.0

    def test_simulate_trial_with_compliance(self):
        """Test trial simulation with compliance."""
        spec = TrialSpec(
            name="Compliance Trial",
            design=parallel_design(2),
            arms=[
                TreatmentArm("Active", dosing_qd(100.0, 28), 50),
            ],
            duration_days=28,
            compliance=ComplianceSpec(mean_compliance=0.85),
            seed=42,
        )

        result = simulate_trial(spec)

        assert result.overall_compliance < 1.0

    def test_simulate_trial_replicates(self):
        """Test multiple trial replicates."""
        spec = TrialSpec(
            name="Multi Trial",
            design=parallel_design(2),
            arms=[
                TreatmentArm("Placebo", dosing_qd(0.0, 14), 20, placebo=True),
                TreatmentArm("Active", dosing_qd(100.0, 14), 20),
            ],
            duration_days=14,
        )

        results = simulate_trial_replicates(spec, n_replicates=5, seed=42)

        assert len(results) == 5
        assert all(r.replicate == i + 1 for i, r in enumerate(results))


# ============================================================================
# Statistical Analysis Tests
# ============================================================================

class TestStatisticalAnalysis:
    """Tests for statistical analysis functions."""

    def test_estimate_power_analytical(self):
        """Test analytical power estimation."""
        result = estimate_power_analytical(
            n_per_arm=50,
            effect_size=0.5,
            alpha=0.05,
        )

        assert isinstance(result, PowerResult)
        assert 0.0 <= result.power <= 1.0
        assert result.alpha == 0.05
        assert result.n_per_arm == 50

    def test_power_increases_with_n(self):
        """Test that power increases with sample size."""
        power_small = estimate_power_analytical(n_per_arm=20, effect_size=0.5).power
        power_large = estimate_power_analytical(n_per_arm=100, effect_size=0.5).power

        assert power_large > power_small

    def test_power_increases_with_effect_size(self):
        """Test that power increases with effect size."""
        power_small = estimate_power_analytical(n_per_arm=50, effect_size=0.3).power
        power_large = estimate_power_analytical(n_per_arm=50, effect_size=0.8).power

        assert power_large > power_small

    def test_estimate_sample_size(self):
        """Test sample size estimation."""
        result = estimate_sample_size(
            target_power=0.80,
            effect_size=0.5,
        )

        assert isinstance(result, SampleSizeResult)
        assert result.achieved_power >= result.target_power
        assert result.total_n == result.n_per_arm * 2

    def test_sample_size_80_power(self):
        """Test sample size for 80% power with d=0.5."""
        result = estimate_sample_size(target_power=0.80, effect_size=0.5)

        # Should be around 64 per arm for d=0.5
        assert 50 <= result.n_per_arm <= 80

    def test_alpha_spending_obrien_fleming(self):
        """Test O'Brien-Fleming alpha spending function."""
        alpha_50 = alpha_spending_function(0.5, 0.05, "obrien_fleming")
        alpha_100 = alpha_spending_function(1.0, 0.05, "obrien_fleming")

        assert alpha_50 < alpha_100
        assert alpha_100 == pytest.approx(0.05, rel=0.01)

    def test_alpha_spending_pocock(self):
        """Test Pocock alpha spending function."""
        alpha_50 = alpha_spending_function(0.5, 0.05, "pocock")
        alpha_100 = alpha_spending_function(1.0, 0.05, "pocock")

        assert alpha_50 < alpha_100
        assert alpha_100 == pytest.approx(0.05, rel=0.01)

    def test_incremental_alpha(self):
        """Test incremental alpha calculation."""
        fractions = [0.5, 0.75, 1.0]
        alphas = incremental_alpha(fractions, 0.05)

        assert len(alphas) == 3
        assert sum(alphas) == pytest.approx(0.05, rel=0.01)

    def test_compare_arms(self):
        """Test arm comparison."""
        arm1_values = [10.2, 11.5, 9.8, 12.1, 10.5, 11.0, 10.8]
        arm2_values = [8.1, 7.9, 8.5, 7.2, 8.8, 7.5, 8.0]

        result = compare_arms(arm1_values, arm2_values, "Active", "Placebo")

        assert isinstance(result, ComparisonResult)
        assert result.arm1 == "Active"
        assert result.arm2 == "Placebo"
        assert result.difference > 0  # Active has higher values
        assert result.ci_lower < result.difference < result.ci_upper

    def test_compare_arms_significant(self):
        """Test that significant difference is detected."""
        arm1 = [20.0] * 30
        arm2 = [10.0] * 30

        result = compare_arms(arm1, arm2)

        assert result.significant
        assert result.p_value < 0.05

    def test_responder_analysis(self):
        """Test responder analysis."""
        values = [10.2, 15.5, 8.1, 12.3, 9.5, 14.2, 11.8, 13.0, 7.5, 16.0]
        result = responder_analysis(values, threshold=10.0, direction="greater")

        assert isinstance(result, ResponderResult)
        assert result.n_total == 10
        assert 0.0 <= result.response_rate <= 1.0
        assert result.ci_lower <= result.response_rate <= result.ci_upper

    def test_responder_analysis_less_than(self):
        """Test responder analysis with 'less' direction."""
        values = [5.0, 15.0, 3.0, 20.0, 2.0]
        result = responder_analysis(values, threshold=10.0, direction="less")

        assert result.n_responders == 3
        assert result.response_rate == 0.6


# ============================================================================
# Bioequivalence Tests
# ============================================================================

class TestBioequivalence:
    """Tests for bioequivalence analysis."""

    def test_bioequivalence_90ci(self):
        """Test 90% CI calculation for BE."""
        test = [95.2, 102.1, 98.5, 105.3, 97.8, 100.2, 99.5, 101.8]
        ref = [100.0, 98.5, 101.2, 99.8, 100.5, 97.2, 102.0, 98.0]

        lower, upper = bioequivalence_90ci(test, ref)

        assert lower < upper
        assert lower > 0

    def test_assess_bioequivalence_pass(self):
        """Test BE assessment that passes."""
        # Similar values should pass BE
        test = [100.0] * 10
        ref = [100.0] * 10

        result = assess_bioequivalence(test, ref)

        assert result["bioequivalent"]
        assert 0.80 <= result["ci_90_lower"]
        assert result["ci_90_upper"] <= 1.25

    def test_assess_bioequivalence_fail(self):
        """Test BE assessment that fails."""
        # Very different values should fail BE
        test = [50.0] * 10
        ref = [100.0] * 10

        result = assess_bioequivalence(test, ref)

        assert not result["bioequivalent"]

    def test_bioequivalence_custom_limits(self):
        """Test BE with custom limits (highly variable drugs)."""
        test = [95.0, 100.0, 105.0, 98.0, 102.0]
        ref = [100.0, 100.0, 100.0, 100.0, 100.0]

        result = assess_bioequivalence(
            test, ref,
            lower_limit=0.6984,
            upper_limit=1.4319,
        )

        assert result["be_lower_limit"] == 0.6984
        assert result["be_upper_limit"] == 1.4319


# ============================================================================
# Integration Tests
# ============================================================================

class TestIntegration:
    """Integration tests combining multiple components."""

    def test_full_trial_workflow(self):
        """Test complete trial simulation workflow."""
        # 1. Create demographics
        demo = DemographicSpec(
            age_mean=50.0,
            age_sd=12.0,
            weight_mean=75.0,
            female_proportion=0.5,
        )

        # 2. Create design
        design = parallel_design(2)

        # 3. Create regimens
        placebo_regimen = dosing_qd(0.0, 28)
        active_regimen = dosing_qd(100.0, 28)

        # 4. Create trial spec
        spec = TrialSpec(
            name="Phase 2 Study",
            design=design,
            arms=[
                TreatmentArm("Placebo", placebo_regimen, 30, placebo=True),
                TreatmentArm("Active", active_regimen, 30),
            ],
            population_spec=VirtualPopulationSpec(demographics=demo),
            duration_days=28,
            dropout=DropoutSpec(random_rate_per_day=0.005),
            compliance=ComplianceSpec(mean_compliance=0.90),
            seed=42,
        )

        # 5. Run simulation
        result = simulate_trial(spec)

        # 6. Verify results
        assert result.trial_name == "Phase 2 Study"
        assert len(result.arms) == 2
        assert result.overall_completion_rate > 0.5

        # Get arm results
        placebo_arm = result.arms["Placebo"]
        active_arm = result.arms["Active"]

        assert placebo_arm.n_enrolled == 30
        assert active_arm.n_enrolled == 30

    def test_crossover_be_workflow(self):
        """Test crossover bioequivalence workflow."""
        # Create BE design
        design = bioequivalence_design(
            n_periods=2,
            n_sequences=2,
            washout_duration=14.0,
            regulatory_guidance="fda",
        )

        # Verify design
        assert design.bioequivalence_limits == (0.80, 1.25)

        # Simulate BE data
        import random
        random.seed(42)

        test_values = [random.gauss(100, 10) for _ in range(24)]
        ref_values = [random.gauss(100, 10) for _ in range(24)]

        # Run BE analysis
        result = assess_bioequivalence(test_values, ref_values)

        assert "bioequivalent" in result
        assert "ci_90_lower" in result
        assert "ci_90_upper" in result

    def test_dose_escalation_workflow(self):
        """Test dose escalation workflow."""
        dose_levels = [10.0, 25.0, 50.0, 100.0, 200.0]
        design = dose_escalation_3plus3(dose_levels)

        # Verify design properties
        assert design.starting_dose == 10.0
        assert len(design.dose_levels) == 5
        assert design.cohort_size == 3

        # Create regimen for first dose level
        regimen = dosing_qd(design.starting_dose, 7)
        assert regimen.dose_amount == 10.0

    def test_adaptive_trial_workflow(self):
        """Test adaptive trial workflow."""
        # Create base design
        base = parallel_design(2)

        # Create adaptive wrapper
        design = adaptive_design(
            base,
            interim_analyses=[0.5],
            alpha_spending="obrien_fleming",
        )

        # Calculate alpha at interim
        alpha_interim = alpha_spending_function(0.5, 0.05, "obrien_fleming")
        alpha_final = 0.05 - alpha_interim

        assert alpha_interim < 0.05
        assert alpha_final > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
