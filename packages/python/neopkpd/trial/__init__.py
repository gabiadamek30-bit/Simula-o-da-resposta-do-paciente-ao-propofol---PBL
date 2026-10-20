"""
NeoPKPD Trial Module

Clinical trial simulation for pharmacometrics including:
- Study designs (parallel, crossover, dose-escalation, adaptive, bioequivalence)
- Dosing regimens (QD, BID, TID, titration)
- Virtual population generation
- Trial simulation with dropout, compliance, endpoints
- Power analysis and sample size estimation

Example:
    >>> import neopkpd
    >>> from neopkpd import trial
    >>>
    >>> # Create a parallel design
    >>> design = trial.parallel_design(2)
    >>>
    >>> # Create dosing regimen
    >>> regimen = trial.dosing_qd(100.0, 28)
    >>>
    >>> # Generate virtual population
    >>> pop = trial.generate_virtual_population(100)
    >>>
    >>> # Estimate sample size
    >>> result = trial.estimate_sample_size(0.80, 0.5, 1.0)
"""

from .designs import (
    parallel_design,
    crossover_2x2,
    crossover_3x3,
    williams_design,
    dose_escalation_3plus3,
    dose_escalation_mtpi,
    dose_escalation_crm,
    adaptive_design,
    adaptive_design_full,
    bioequivalence_design,
    get_design_description,
    ParallelDesign,
    CrossoverDesign,
    DoseEscalationDesign,
    BioequivalenceDesign,
    AdaptiveDesign,
    FullAdaptiveDesign,
    # Adaptive design configurations
    ResponseAdaptiveRandomization,
    SampleSizeReestimation,
    TreatmentSelection,
    BiomarkerEnrichment,
    # Enums
    RARMethod,
    SSRMethod,
    SelectionCriterion,
)

from .regimens import (
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
)

from .population import (
    generate_virtual_population,
    default_demographic_spec,
    healthy_volunteer_spec,
    patient_population_spec,
    summarize_population,
    DemographicSpec,
    DiseaseSpec,
    VirtualPopulationSpec,
    VirtualSubject,
)

from .simulation import (
    simulate_trial,
    simulate_trial_replicates,
    simulate_dropout,
    apply_compliance,
    TrialSpec,
    TrialResult,
    TreatmentArm,
    DropoutSpec,
    ComplianceSpec,
    SubjectResult,
    ArmResult,
)

from .analysis import (
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

from .core_integration import (
    # Subject exposure simulation
    simulate_subject_exposure,
    SubjectExposure,
    # Dose escalation
    simulate_dose_escalation_3plus3,
    simulate_dose_escalation_mtpi,
    simulate_dose_escalation_crm,
    CohortResult,
    EscalationResult,
    # Model-connected trial simulation
    simulate_trial_with_model,
    ModelTrialSubjectResult,
    ModelTrialArmResult,
    ModelTrialResult,
    # Crossover analysis
    analyze_crossover,
    test_period_effect,
    test_sequence_effect,
    compute_within_subject_cv,
    CrossoverAnalysis,
    # Adaptive trials
    simulate_adaptive_trial,
    InterimResult,
    AdaptiveTrialResult,
    # PK metrics
    calculate_pk_metrics,
)


__all__ = [
    # Designs
    "parallel_design",
    "crossover_2x2",
    "crossover_3x3",
    "williams_design",
    "dose_escalation_3plus3",
    "dose_escalation_mtpi",
    "dose_escalation_crm",
    "adaptive_design",
    "adaptive_design_full",
    "bioequivalence_design",
    "get_design_description",
    "ParallelDesign",
    "CrossoverDesign",
    "DoseEscalationDesign",
    "BioequivalenceDesign",
    "AdaptiveDesign",
    "FullAdaptiveDesign",
    # Adaptive design configurations
    "ResponseAdaptiveRandomization",
    "SampleSizeReestimation",
    "TreatmentSelection",
    "BiomarkerEnrichment",
    "RARMethod",
    "SSRMethod",
    "SelectionCriterion",
    # Regimens
    "dosing_qd",
    "dosing_bid",
    "dosing_tid",
    "dosing_qid",
    "dosing_custom",
    "titration_regimen",
    "dose_event_times",
    "total_regimen_duration",
    "DosingRegimen",
    "TitrationRegimen",
    # Population
    "generate_virtual_population",
    "default_demographic_spec",
    "healthy_volunteer_spec",
    "patient_population_spec",
    "summarize_population",
    "DemographicSpec",
    "DiseaseSpec",
    "VirtualPopulationSpec",
    "VirtualSubject",
    # Simulation
    "simulate_trial",
    "simulate_trial_replicates",
    "simulate_dropout",
    "apply_compliance",
    "TrialSpec",
    "TrialResult",
    "TreatmentArm",
    "DropoutSpec",
    "ComplianceSpec",
    "SubjectResult",
    "ArmResult",
    # Analysis
    "estimate_power_analytical",
    "estimate_sample_size",
    "alpha_spending_function",
    "incremental_alpha",
    "compare_arms",
    "responder_analysis",
    "bioequivalence_90ci",
    "assess_bioequivalence",
    "PowerResult",
    "SampleSizeResult",
    "ComparisonResult",
    "ResponderResult",
    # Core Integration - Subject Exposure
    "simulate_subject_exposure",
    "SubjectExposure",
    # Core Integration - Dose Escalation
    "simulate_dose_escalation_3plus3",
    "simulate_dose_escalation_mtpi",
    "simulate_dose_escalation_crm",
    "CohortResult",
    "EscalationResult",
    # Core Integration - Model-Connected Trial
    "simulate_trial_with_model",
    "ModelTrialSubjectResult",
    "ModelTrialArmResult",
    "ModelTrialResult",
    # Core Integration - Crossover Analysis
    "analyze_crossover",
    "test_period_effect",
    "test_sequence_effect",
    "compute_within_subject_cv",
    "CrossoverAnalysis",
    # Core Integration - Adaptive Trials
    "simulate_adaptive_trial",
    "InterimResult",
    "AdaptiveTrialResult",
    # Core Integration - PK Metrics
    "calculate_pk_metrics",
]
