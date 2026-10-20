# Trial Specifications
# Type definitions for clinical trial simulation

export StudyDesignKind, ParallelDesign, CrossoverDesign, DoseEscalationDesign
export AdaptiveDesign, BioequivalenceDesign
export DosingFrequency, QD, BID, TID, QID, CustomFrequency
export DosingRegimen, TitrationStep, TitrationRegimen
export DemographicSpec, DiseaseSpec, VirtualPopulationSpec, VirtualSubject
export DropoutSpec, ComplianceSpec, EnrollmentSpec
export EndpointSpec, PKEndpoint, PDEndpoint, SafetyEndpoint, CompositeEndpoint
export TreatmentArm, TrialSpec, TrialResult, ArmResult
export EscalationRule, ThreePlusThree, mTPI, CRM, BOIN

using Random
using StableRNGs

# ============================================================================
# Study Design Types
# ============================================================================

"""
Abstract type for study design kinds.
"""
abstract type StudyDesignKind end

"""
    ParallelDesign

Parallel group study design.

# Fields
- `n_arms::Int`: Number of treatment arms
- `randomization_ratio::Vector{Float64}`: Randomization ratio for each arm
- `stratification_factors::Vector{Symbol}`: Factors for stratified randomization
"""
struct ParallelDesign <: StudyDesignKind
    n_arms::Int
    randomization_ratio::Vector{Float64}
    stratification_factors::Vector{Symbol}

    function ParallelDesign(n_arms::Int;
                            randomization_ratio::Vector{Float64} = ones(n_arms) / n_arms,
                            stratification_factors::Vector{Symbol} = Symbol[])
        @assert n_arms >= 2 "Parallel design requires at least 2 arms"
        @assert length(randomization_ratio) == n_arms "Randomization ratio must match number of arms"
        @assert abs(sum(randomization_ratio) - 1.0) < 1e-10 "Randomization ratios must sum to 1"
        new(n_arms, randomization_ratio, stratification_factors)
    end
end

"""
    CrossoverDesign

Crossover study design.

# Fields
- `n_periods::Int`: Number of treatment periods
- `n_sequences::Int`: Number of treatment sequences
- `washout_duration::Float64`: Washout period duration (days)
- `sequence_assignments::Vector{Vector{Int}}`: Treatment assignments per sequence
"""
struct CrossoverDesign <: StudyDesignKind
    n_periods::Int
    n_sequences::Int
    washout_duration::Float64
    sequence_assignments::Vector{Vector{Int}}

    function CrossoverDesign(n_periods::Int, n_sequences::Int;
                              washout_duration::Float64 = 7.0,
                              sequence_assignments::Union{Nothing, Vector{Vector{Int}}} = nothing)
        if sequence_assignments === nothing
            # Default 2x2 crossover: AB, BA
            if n_periods == 2 && n_sequences == 2
                sequence_assignments = [[1, 2], [2, 1]]
            else
                # Generate Latin square for n treatments
                sequence_assignments = [[(i + j - 2) % n_periods + 1 for j in 1:n_periods] for i in 1:n_sequences]
            end
        end
        @assert length(sequence_assignments) == n_sequences
        @assert all(length(s) == n_periods for s in sequence_assignments)
        new(n_periods, n_sequences, washout_duration, sequence_assignments)
    end
end

"""
Abstract type for dose escalation rules.
"""
abstract type EscalationRule end

"""
    ThreePlusThree

Traditional 3+3 dose escalation rule.
"""
struct ThreePlusThree <: EscalationRule
    max_dlt_rate::Float64
    ThreePlusThree(; max_dlt_rate::Float64 = 0.33) = new(max_dlt_rate)
end

"""
    mTPI

Modified Toxicity Probability Interval design.
"""
struct mTPI <: EscalationRule
    target_dlt_rate::Float64
    equivalence_interval::Tuple{Float64, Float64}
    mTPI(; target_dlt_rate::Float64 = 0.25,
           equivalence_interval::Tuple{Float64, Float64} = (0.20, 0.30)) =
        new(target_dlt_rate, equivalence_interval)
end

"""
    CRM

Continual Reassessment Method.
"""
struct CRM <: EscalationRule
    target_dlt_rate::Float64
    skeleton::Vector{Float64}
    model::Symbol  # :logistic, :power
    CRM(; target_dlt_rate::Float64 = 0.25,
          skeleton::Vector{Float64} = [0.05, 0.10, 0.20, 0.30, 0.50],
          model::Symbol = :logistic) = new(target_dlt_rate, skeleton, model)
end

"""
    BOIN

Bayesian Optimal Interval Design for Phase I dose escalation.

BOIN uses pre-calculated optimal boundaries for dose escalation decisions
based on Bayesian decision theory.

# Fields
- `target_dlt_rate::Float64`: Target DLT rate (default 0.25)
- `p1::Float64`: Lower bound for acceptable DLT rate (default 0.6 * target)
- `p2::Float64`: Upper bound for acceptable DLT rate (default 1.4 * target)
- `cutoff_e::Float64`: Escalation boundary (default: computed from p1, target)
- `cutoff_d::Float64`: De-escalation boundary (default: computed from target, p2)
- `n_eliminate::Int`: Min subjects at dose before elimination check (default: 3)
- `eliminate_threshold::Float64`: Threshold for dose elimination (default: 0.95)

# Reference
Liu S, Yuan Y (2015). Bayesian Optimal Interval Designs for Phase I Clinical Trials.
Journal of the Royal Statistical Society: Series C.
"""
struct BOIN <: EscalationRule
    target_dlt_rate::Float64
    p1::Float64
    p2::Float64
    cutoff_e::Float64
    cutoff_d::Float64
    n_eliminate::Int
    eliminate_threshold::Float64

    function BOIN(; target_dlt_rate::Float64 = 0.25,
                    p1::Union{Float64, Nothing} = nothing,
                    p2::Union{Float64, Nothing} = nothing,
                    cutoff_e::Union{Float64, Nothing} = nothing,
                    cutoff_d::Union{Float64, Nothing} = nothing,
                    n_eliminate::Int = 3,
                    eliminate_threshold::Float64 = 0.95)
        # Default interval bounds
        if p1 === nothing
            p1 = 0.6 * target_dlt_rate
        end
        if p2 === nothing
            p2 = 1.4 * target_dlt_rate
        end

        @assert 0 < p1 < target_dlt_rate < p2 < 1 "Must have 0 < p1 < target < p2 < 1"

        # Compute optimal boundaries if not provided
        # λe = log((1-p1)/(1-target)) / log(target*(1-p1)/(p1*(1-target)))
        # λd = log((1-target)/(1-p2)) / log(p2*(1-target)/(target*(1-p2)))
        if cutoff_e === nothing
            cutoff_e = log((1 - p1) / (1 - target_dlt_rate)) /
                       log(target_dlt_rate * (1 - p1) / (p1 * (1 - target_dlt_rate)))
        end
        if cutoff_d === nothing
            cutoff_d = log((1 - target_dlt_rate) / (1 - p2)) /
                       log(p2 * (1 - target_dlt_rate) / (target_dlt_rate * (1 - p2)))
        end

        @assert cutoff_e < cutoff_d "Escalation cutoff must be less than de-escalation cutoff"

        new(target_dlt_rate, p1, p2, cutoff_e, cutoff_d, n_eliminate, eliminate_threshold)
    end
end

"""
    DoseEscalationDesign

Dose escalation study design for Phase I trials.

# Fields
- `starting_dose::Float64`: Starting dose level
- `dose_levels::Vector{Float64}`: Available dose levels
- `escalation_rule::EscalationRule`: Rule for dose escalation decisions
- `cohort_size::Int`: Number of subjects per cohort
- `max_subjects::Int`: Maximum total subjects
"""
struct DoseEscalationDesign <: StudyDesignKind
    starting_dose::Float64
    dose_levels::Vector{Float64}
    escalation_rule::EscalationRule
    cohort_size::Int
    max_subjects::Int

    function DoseEscalationDesign(dose_levels::Vector{Float64};
                                   starting_dose::Union{Float64, Nothing} = nothing,
                                   escalation_rule::EscalationRule = ThreePlusThree(),
                                   cohort_size::Int = 3,
                                   max_subjects::Int = 30)
        if starting_dose === nothing
            starting_dose = dose_levels[1]
        end
        @assert starting_dose in dose_levels "Starting dose must be in dose levels"
        new(starting_dose, sort(dose_levels), escalation_rule, cohort_size, max_subjects)
    end
end

"""
    AdaptiveDesign

Adaptive trial design with interim analyses.

# Fields
- `base_design::StudyDesignKind`: Underlying design
- `interim_analyses::Vector{Float64}`: Information fractions for interim analyses
- `adaptation_rules::Dict{Symbol, Any}`: Rules for adaptations
- `alpha_spending::Symbol`: Alpha spending function (:obrien_fleming, :pocock, :haybittle_peto)
"""
struct AdaptiveDesign <: StudyDesignKind
    base_design::StudyDesignKind
    interim_analyses::Vector{Float64}
    adaptation_rules::Dict{Symbol, Any}
    alpha_spending::Symbol

    function AdaptiveDesign(base_design::StudyDesignKind;
                            interim_analyses::Vector{Float64} = [0.5],
                            adaptation_rules::Dict{Symbol, Any} = Dict{Symbol, Any}(),
                            alpha_spending::Symbol = :obrien_fleming)
        @assert all(0 < ia < 1 for ia in interim_analyses) "Interim analyses must be between 0 and 1"
        new(base_design, sort(interim_analyses), adaptation_rules, alpha_spending)
    end
end

"""
    BioequivalenceDesign

Bioequivalence study design.

# Fields
- `crossover::CrossoverDesign`: Underlying crossover design
- `bioequivalence_limits::Tuple{Float64, Float64}`: BE acceptance limits (default 0.80, 1.25)
- `parameters::Vector{Symbol}`: PK parameters for BE assessment
- `regulatory_guidance::Symbol`: Regulatory guidance to follow (:fda, :ema)
"""
struct BioequivalenceDesign <: StudyDesignKind
    crossover::CrossoverDesign
    bioequivalence_limits::Tuple{Float64, Float64}
    parameters::Vector{Symbol}
    regulatory_guidance::Symbol

    function BioequivalenceDesign(;
            n_periods::Int = 2,
            n_sequences::Int = 2,
            washout_duration::Float64 = 7.0,
            bioequivalence_limits::Tuple{Float64, Float64} = (0.80, 1.25),
            parameters::Vector{Symbol} = [:cmax, :auc_0_inf],
            regulatory_guidance::Symbol = :fda)
        crossover = CrossoverDesign(n_periods, n_sequences; washout_duration = washout_duration)
        new(crossover, bioequivalence_limits, parameters, regulatory_guidance)
    end
end

# ============================================================================
# Dosing Regimen Types
# ============================================================================

"""
Abstract type for dosing frequency.
"""
abstract type DosingFrequency end

struct QD <: DosingFrequency end  # Once daily
struct BID <: DosingFrequency end  # Twice daily
struct TID <: DosingFrequency end  # Three times daily
struct QID <: DosingFrequency end  # Four times daily

"""
    CustomFrequency

Custom dosing frequency.
"""
struct CustomFrequency <: DosingFrequency
    dose_times::Vector{Float64}  # Hours after midnight
end

"""
    DosingRegimen

Dosing regimen specification.

# Fields
- `frequency::DosingFrequency`: Dosing frequency
- `dose_amount::Float64`: Dose amount per administration
- `duration_days::Int`: Duration of treatment in days
- `loading_dose::Union{Float64, Nothing}`: Optional loading dose
- `dose_times::Vector{Float64}`: Explicit dose times (hours)
"""
struct DosingRegimen
    frequency::DosingFrequency
    dose_amount::Float64
    duration_days::Int
    loading_dose::Union{Float64, Nothing}
    dose_times::Vector{Float64}  # Hours after midnight for each administration

    function DosingRegimen(frequency::DosingFrequency, dose_amount::Float64, duration_days::Int;
                           loading_dose::Union{Float64, Nothing} = nothing)
        dose_times = get_dose_times(frequency)
        new(frequency, dose_amount, duration_days, loading_dose, dose_times)
    end
end

function get_dose_times(::QD)
    return [8.0]  # 8 AM
end

function get_dose_times(::BID)
    return [8.0, 20.0]  # 8 AM, 8 PM
end

function get_dose_times(::TID)
    return [8.0, 14.0, 20.0]  # 8 AM, 2 PM, 8 PM
end

function get_dose_times(::QID)
    return [8.0, 12.0, 16.0, 20.0]  # 8 AM, 12 PM, 4 PM, 8 PM
end

function get_dose_times(f::CustomFrequency)
    return f.dose_times
end

"""
    TitrationStep

Single step in a titration regimen.
"""
struct TitrationStep
    dose::Float64
    duration_days::Int
end

"""
    TitrationRegimen

Titration dosing regimen with gradual dose changes.

# Fields
- `steps::Vector{TitrationStep}`: Titration steps
- `frequency::DosingFrequency`: Dosing frequency
- `loading_dose::Union{Float64, Nothing}`: Optional loading dose
"""
struct TitrationRegimen
    steps::Vector{TitrationStep}
    frequency::DosingFrequency
    loading_dose::Union{Float64, Nothing}

    function TitrationRegimen(steps::Vector{TitrationStep};
                               frequency::DosingFrequency = QD(),
                               loading_dose::Union{Float64, Nothing} = nothing)
        new(steps, frequency, loading_dose)
    end
end

# ============================================================================
# Virtual Population Types
# ============================================================================

"""
    DemographicSpec

Demographic distribution specification.

# Fields
- `age_mean::Float64`: Mean age
- `age_sd::Float64`: Age standard deviation
- `age_range::Tuple{Float64, Float64}`: Age range limits
- `weight_mean::Float64`: Mean weight (kg)
- `weight_sd::Float64`: Weight standard deviation
- `weight_range::Tuple{Float64, Float64}`: Weight range limits
- `female_proportion::Float64`: Proportion of females
- `race_distribution::Dict{Symbol, Float64}`: Race/ethnicity distribution
"""
struct DemographicSpec
    age_mean::Float64
    age_sd::Float64
    age_range::Tuple{Float64, Float64}
    weight_mean::Float64
    weight_sd::Float64
    weight_range::Tuple{Float64, Float64}
    female_proportion::Float64
    race_distribution::Dict{Symbol, Float64}

    function DemographicSpec(;
            age_mean::Float64 = 45.0,
            age_sd::Float64 = 15.0,
            age_range::Tuple{Float64, Float64} = (18.0, 75.0),
            weight_mean::Float64 = 75.0,
            weight_sd::Float64 = 15.0,
            weight_range::Tuple{Float64, Float64} = (40.0, 150.0),
            female_proportion::Float64 = 0.5,
            race_distribution::Dict{Symbol, Float64} = Dict(:caucasian => 0.70, :asian => 0.15, :black => 0.10, :other => 0.05))
        @assert 0 <= female_proportion <= 1 "Female proportion must be between 0 and 1"
        @assert abs(sum(values(race_distribution)) - 1.0) < 1e-10 "Race distribution must sum to 1"
        new(age_mean, age_sd, age_range, weight_mean, weight_sd, weight_range,
            female_proportion, race_distribution)
    end
end

"""
    DiseaseSpec

Disease state specification.

# Fields
- `name::Symbol`: Disease name
- `severity_distribution::Dict{Symbol, Float64}`: Distribution of severity levels
- `baseline_biomarker_mean::Float64`: Mean baseline biomarker value
- `baseline_biomarker_sd::Float64`: Baseline biomarker standard deviation
"""
struct DiseaseSpec
    name::Symbol
    severity_distribution::Dict{Symbol, Float64}
    baseline_biomarker_mean::Float64
    baseline_biomarker_sd::Float64

    function DiseaseSpec(name::Symbol;
            severity_distribution::Dict{Symbol, Float64} = Dict(:mild => 0.3, :moderate => 0.5, :severe => 0.2),
            baseline_biomarker_mean::Float64 = 100.0,
            baseline_biomarker_sd::Float64 = 25.0)
        new(name, severity_distribution, baseline_biomarker_mean, baseline_biomarker_sd)
    end
end

"""
    VirtualPopulationSpec

Specification for virtual population generation.

# Fields
- `demographics::DemographicSpec`: Demographic distribution
- `disease::Union{Nothing, DiseaseSpec}`: Disease specification
- `covariate_correlations::Dict{Tuple{Symbol, Symbol}, Float64}`: Covariate correlations
- `seed::UInt64`: Random seed for reproducibility
"""
struct VirtualPopulationSpec
    demographics::DemographicSpec
    disease::Union{Nothing, DiseaseSpec}
    covariate_correlations::Dict{Tuple{Symbol, Symbol}, Float64}
    seed::UInt64

    function VirtualPopulationSpec(;
            demographics::DemographicSpec = DemographicSpec(),
            disease::Union{Nothing, DiseaseSpec} = nothing,
            covariate_correlations::Dict{Tuple{Symbol, Symbol}, Float64} = Dict{Tuple{Symbol, Symbol}, Float64}(),
            seed::UInt64 = UInt64(12345))
        new(demographics, disease, covariate_correlations, seed)
    end
end

"""
    VirtualSubject

Virtual subject for trial simulation.

# Fields
- `id::Int`: Subject ID
- `age::Float64`: Age (years)
- `weight::Float64`: Weight (kg)
- `sex::Symbol`: Sex (:male or :female)
- `race::Symbol`: Race/ethnicity
- `disease_severity::Union{Nothing, Symbol}`: Disease severity level
- `baseline_biomarker::Union{Nothing, Float64}`: Baseline biomarker value
- `other::Dict{Symbol, Any}`: Other covariates (can store scalars, vectors, or any value type)
- `eta::Vector{Float64}`: Individual random effects (IIV parameters)
"""
struct VirtualSubject
    id::Int
    age::Float64
    weight::Float64
    sex::Symbol
    race::Symbol
    disease_severity::Union{Nothing, Symbol}
    baseline_biomarker::Union{Nothing, Float64}
    other::Dict{Symbol, Any}
    eta::Vector{Float64}
end

"""
    VirtualSubject(id, age, weight, sex, race; kwargs...)

Convenience constructor for VirtualSubject with optional fields.
"""
function VirtualSubject(
    id::Int,
    age::Float64,
    weight::Float64,
    sex::Symbol,
    race::Symbol;
    disease_severity::Union{Nothing, Symbol}=nothing,
    baseline_biomarker::Union{Nothing, Float64}=nothing,
    other::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    eta::Vector{Float64}=Float64[]
)
    VirtualSubject(id, age, weight, sex, race, disease_severity, baseline_biomarker, other, eta)
end

# ============================================================================
# Trial Event Types
# ============================================================================

"""
    DropoutSpec

Dropout specification.

# Fields
- `random_rate_per_day::Float64`: Random dropout rate per day
- `ae_threshold::Union{Float64, Nothing}`: Adverse event threshold for dropout
- `ae_dropout_prob::Float64`: Probability of dropout given AE threshold exceeded
- `non_compliance_threshold::Float64`: Compliance threshold for dropout
"""
struct DropoutSpec
    random_rate_per_day::Float64
    ae_threshold::Union{Float64, Nothing}
    ae_dropout_prob::Float64
    non_compliance_threshold::Float64

    function DropoutSpec(;
            random_rate_per_day::Float64 = 0.001,
            ae_threshold::Union{Float64, Nothing} = nothing,
            ae_dropout_prob::Float64 = 0.5,
            non_compliance_threshold::Float64 = 0.5)
        new(random_rate_per_day, ae_threshold, ae_dropout_prob, non_compliance_threshold)
    end
end

"""
    ComplianceSpec

Compliance specification.

# Fields
- `mean_compliance::Float64`: Mean compliance rate (0-1)
- `compliance_sd::Float64`: Compliance standard deviation
- `pattern::Symbol`: Compliance pattern (:random, :weekend_miss, :decay, :early_good)
- `min_compliance::Float64`: Minimum compliance rate
"""
struct ComplianceSpec
    mean_compliance::Float64
    compliance_sd::Float64
    pattern::Symbol
    min_compliance::Float64

    function ComplianceSpec(;
            mean_compliance::Float64 = 0.90,
            compliance_sd::Float64 = 0.10,
            pattern::Symbol = :random,
            min_compliance::Float64 = 0.50)
        @assert 0 < mean_compliance <= 1 "Mean compliance must be between 0 and 1"
        new(mean_compliance, compliance_sd, pattern, min_compliance)
    end
end

"""
    EnrollmentSpec

Enrollment specification.

# Fields
- `rate_per_day::Float64`: Average enrollment rate (subjects per day)
- `max_concurrent::Int`: Maximum concurrent subjects
- `enrollment_cap::Int`: Total enrollment cap
- `screening_failure_rate::Float64`: Screening failure rate
"""
struct EnrollmentSpec
    rate_per_day::Float64
    max_concurrent::Int
    enrollment_cap::Int
    screening_failure_rate::Float64

    function EnrollmentSpec(;
            rate_per_day::Float64 = 2.0,
            max_concurrent::Int = 100,
            enrollment_cap::Int = 1000,
            screening_failure_rate::Float64 = 0.20)
        new(rate_per_day, max_concurrent, enrollment_cap, screening_failure_rate)
    end
end

# ============================================================================
# Endpoint Types
# ============================================================================

"""
Abstract type for trial endpoints.
"""
abstract type EndpointSpec end

"""
    PKEndpoint

Pharmacokinetic endpoint.

# Fields
- `name::Symbol`: Endpoint name
- `metric::Symbol`: PK metric (:cmax, :auc_0_inf, :auc_0_t, :tmax, :t_half)
- `timepoints::Vector{Float64}`: Sampling timepoints
- `log_transform::Bool`: Whether to log-transform for analysis
"""
struct PKEndpoint <: EndpointSpec
    name::Symbol
    metric::Symbol
    timepoints::Vector{Float64}
    log_transform::Bool

    function PKEndpoint(name::Symbol;
            metric::Symbol = :auc_0_inf,
            timepoints::Vector{Float64} = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            log_transform::Bool = true)
        new(name, metric, timepoints, log_transform)
    end
end

"""
    PDEndpoint

Pharmacodynamic endpoint.

# Fields
- `name::Symbol`: Endpoint name
- `metric::Symbol`: PD metric (:emax, :ec50, :change_from_baseline)
- `assessment_times::Vector{Float64}`: Assessment timepoints
- `responder_threshold::Union{Float64, Nothing}`: Threshold for responder analysis
"""
struct PDEndpoint <: EndpointSpec
    name::Symbol
    metric::Symbol
    assessment_times::Vector{Float64}
    responder_threshold::Union{Float64, Nothing}

    function PDEndpoint(name::Symbol;
            metric::Symbol = :change_from_baseline,
            assessment_times::Vector{Float64} = [0.0, 7.0, 14.0, 28.0],
            responder_threshold::Union{Float64, Nothing} = nothing)
        new(name, metric, assessment_times, responder_threshold)
    end
end

"""
    SafetyEndpoint

Safety endpoint.

# Fields
- `name::Symbol`: Endpoint name
- `type::Symbol`: Type (:ae_rate, :dlt_rate, :sae_rate, :lab_change)
- `severity_levels::Vector{Symbol}`: Severity levels to track
- `threshold::Union{Float64, Nothing}`: Safety threshold
"""
struct SafetyEndpoint <: EndpointSpec
    name::Symbol
    type::Symbol
    severity_levels::Vector{Symbol}
    threshold::Union{Float64, Nothing}

    function SafetyEndpoint(name::Symbol;
            type::Symbol = :ae_rate,
            severity_levels::Vector{Symbol} = [:mild, :moderate, :severe],
            threshold::Union{Float64, Nothing} = nothing)
        new(name, type, severity_levels, threshold)
    end
end

"""
    CompositeEndpoint

Composite endpoint combining multiple endpoints.

# Fields
- `name::Symbol`: Endpoint name
- `components::Vector{EndpointSpec}`: Component endpoints
- `combination_rule::Symbol`: How to combine (:any, :all, :weighted)
- `weights::Vector{Float64}`: Weights for weighted combination
"""
struct CompositeEndpoint <: EndpointSpec
    name::Symbol
    components::Vector{EndpointSpec}
    combination_rule::Symbol
    weights::Vector{Float64}

    function CompositeEndpoint(name::Symbol, components::Vector{<:EndpointSpec};
            combination_rule::Symbol = :any,
            weights::Union{Nothing, Vector{Float64}} = nothing)
        if weights === nothing
            weights = ones(length(components)) / length(components)
        end
        new(name, components, combination_rule, weights)
    end
end

# ============================================================================
# Treatment Arm and Trial Specification
# ============================================================================

"""
    TreatmentArm

Treatment arm specification.

# Fields
- `name::String`: Arm name
- `pk_model_spec::Any`: PK model specification
- `pd_spec::Union{Nothing, Any}`: PD specification
- `regimen::Union{DosingRegimen, TitrationRegimen}`: Dosing regimen
- `n_subjects::Int`: Number of subjects
"""
struct TreatmentArm
    name::String
    pk_model_spec::Any  # ModelSpec from engine
    pd_spec::Union{Nothing, Any}
    regimen::Union{DosingRegimen, TitrationRegimen}
    n_subjects::Int

    function TreatmentArm(name::String, pk_model_spec::Any, regimen::Union{DosingRegimen, TitrationRegimen};
                          pd_spec::Union{Nothing, Any} = nothing,
                          n_subjects::Int = 50)
        new(name, pk_model_spec, pd_spec, regimen, n_subjects)
    end
end

"""
    TrialSpec

Complete trial specification.

# Fields
- `name::String`: Trial name
- `design::StudyDesignKind`: Study design
- `arms::Vector{TreatmentArm}`: Treatment arms
- `virtual_population::VirtualPopulationSpec`: Virtual population specification
- `duration_days::Float64`: Trial duration (days)
- `enrollment::EnrollmentSpec`: Enrollment specification
- `dropout::Union{Nothing, DropoutSpec}`: Dropout specification
- `compliance::Union{Nothing, ComplianceSpec}`: Compliance specification
- `pk_sampling_times::Vector{Float64}`: PK sampling times
- `endpoints::Vector{EndpointSpec}`: Trial endpoints
- `n_replicates::Int`: Number of trial replicates (for power analysis)
- `seed::UInt64`: Random seed
"""
struct TrialSpec
    name::String
    design::StudyDesignKind
    arms::Vector{TreatmentArm}
    virtual_population::VirtualPopulationSpec
    duration_days::Float64
    enrollment::EnrollmentSpec
    dropout::Union{Nothing, DropoutSpec}
    compliance::Union{Nothing, ComplianceSpec}
    pk_sampling_times::Vector{Float64}
    endpoints::Vector{EndpointSpec}
    n_replicates::Int
    seed::UInt64

    function TrialSpec(name::String, design::StudyDesignKind, arms::Vector{TreatmentArm};
            virtual_population::VirtualPopulationSpec = VirtualPopulationSpec(),
            duration_days::Float64 = 28.0,
            enrollment::EnrollmentSpec = EnrollmentSpec(),
            dropout::Union{Nothing, DropoutSpec} = nothing,
            compliance::Union{Nothing, ComplianceSpec} = nothing,
            pk_sampling_times::Vector{Float64} = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            endpoints::Vector{EndpointSpec} = EndpointSpec[],
            n_replicates::Int = 1,
            seed::UInt64 = UInt64(12345))
        new(name, design, arms, virtual_population, duration_days, enrollment,
            dropout, compliance, pk_sampling_times, endpoints, n_replicates, seed)
    end
end

# ============================================================================
# Result Types
# ============================================================================

"""
    ArmResult

Results for a single treatment arm.

# Fields
- `name::String`: Arm name
- `n_enrolled::Int`: Number enrolled
- `n_completed::Int`: Number completed
- `n_dropout::Int`: Number of dropouts
- `individual_results::Vector{Any}`: Individual simulation results
- `endpoint_values::Dict{Symbol, Vector{Float64}}`: Endpoint values
- `summary_stats::Dict{Symbol, Dict{Symbol, Float64}}`: Summary statistics
"""
struct ArmResult
    name::String
    n_enrolled::Int
    n_completed::Int
    n_dropout::Int
    individual_results::Vector{Any}
    endpoint_values::Dict{Symbol, Vector{Float64}}
    summary_stats::Dict{Symbol, Dict{Symbol, Float64}}
end

"""
    TrialResult

Complete trial simulation result.

# Fields
- `trial_name::String`: Trial name
- `design_type::String`: Design type description
- `arms::Dict{String, ArmResult}`: Results by arm
- `endpoint_analyses::Dict{Symbol, Dict{Symbol, Any}}`: Endpoint analysis results
- `power_estimates::Dict{Symbol, Float64}`: Power estimates (if replicates > 1)
- `be_results::Union{Nothing, Dict{Symbol, Any}}`: Bioequivalence results
- `n_replicates::Int`: Number of replicates run
- `seed::UInt64`: Random seed used
"""
struct TrialResult
    trial_name::String
    design_type::String
    arms::Dict{String, ArmResult}
    endpoint_analyses::Dict{Symbol, Dict{Symbol, Any}}
    power_estimates::Dict{Symbol, Float64}
    be_results::Union{Nothing, Dict{Symbol, Any}}
    n_replicates::Int
    seed::UInt64
end

# ============================================================================
# Validation Functions
# ============================================================================

export validate_trial_spec, validate_trial_spec!, TrialValidationResult

"""
    TrialValidationResult

Result of trial specification validation.

# Fields
- `valid::Bool`: Whether the specification is valid
- `errors::Vector{String}`: List of validation errors
- `warnings::Vector{String}`: List of validation warnings
"""
struct TrialValidationResult
    valid::Bool
    errors::Vector{String}
    warnings::Vector{String}
end

"""
    validate_trial_spec(spec::TrialSpec; strict::Bool=false)

Validate a trial specification for completeness and supported features.

# Arguments
- `spec`: Trial specification to validate
- `strict`: If true, treat warnings as errors

# Returns
- `TrialValidationResult` with validation status, errors, and warnings

# Notes
This function checks for:
- Missing pk_model_spec in treatment arms (required for real simulation)
- Unsupported adaptive design features
- Invalid endpoint configurations
"""
function validate_trial_spec(spec::TrialSpec; strict::Bool = false)::TrialValidationResult
    errors = String[]
    warnings = String[]

    # Check that all treatment arms have a pk_model_spec
    for arm in spec.arms
        if arm.pk_model_spec === nothing
            push!(errors, "TreatmentArm '$(arm.name)' has no pk_model_spec. " *
                          "All arms require a ModelSpec for real PK simulation.")
        end
    end

    # Validate adaptive design features (now implemented)
    if spec.design isa AdaptiveDesign
        adaptation_rules = spec.design.adaptation_rules

        # Response-adaptive randomization (RAR) - now supported
        if haskey(adaptation_rules, :response_adaptive_randomization) ||
           haskey(adaptation_rules, :rar) ||
           get(adaptation_rules, :adaptive_randomization, false) == true
            # Validate RAR settings
            rar_method = get(adaptation_rules, :rar_method, :thall_wathen)
            valid_rar_methods = [:thall_wathen, :dbcd, :bayesian]
            if !(rar_method in valid_rar_methods)
                push!(warnings, "Unknown RAR method '$(rar_method)'. " *
                                "Valid options: $(valid_rar_methods). Defaulting to Thall-Wathen.")
            end
            push!(warnings, "Response-adaptive randomization enabled. " *
                            "Ensure minimum allocation constraints are appropriate for the trial.")
        end

        # Sample size re-estimation (SSR) - now fully implemented
        if haskey(adaptation_rules, :sample_size_reestimation) ||
           haskey(adaptation_rules, :ssr)
            ssr_method = get(adaptation_rules, :ssr_method, :conditional_power)
            valid_ssr_methods = [:conditional_power, :variance_based, :promising_zone]
            if !(ssr_method in valid_ssr_methods)
                push!(warnings, "Unknown SSR method '$(ssr_method)'. " *
                                "Valid options: $(valid_ssr_methods)")
            end
            max_increase = get(adaptation_rules, :ssr_max_increase, 2.0)
            if max_increase > 3.0
                push!(warnings, "SSR max increase factor $(max_increase) is large. " *
                                "Consider regulatory implications.")
            end
        end

        # Information adaptive designs - now supported via SSR
        if haskey(adaptation_rules, :information_adaptive)
            push!(warnings, "Information-adaptive design enabled via sample size re-estimation.")
        end

        # Treatment selection (drop arms) - now implemented
        if haskey(adaptation_rules, :drop_arms) ||
           haskey(adaptation_rules, :treatment_selection)
            min_arms = get(adaptation_rules, :min_arms, 2)
            if min_arms < 2
                push!(errors, "Treatment selection must keep at least 2 arms (control + 1 treatment).")
            end
            selection_criterion = get(adaptation_rules, :selection_criterion, :posterior_probability)
            valid_criteria = [:posterior_probability, :frequentist_pvalue, :effect_size]
            if !(selection_criterion in valid_criteria)
                push!(warnings, "Unknown selection criterion '$(selection_criterion)'. " *
                                "Valid options: $(valid_criteria)")
            end
        end

        # Enrichment designs - now implemented
        if haskey(adaptation_rules, :enrichment) ||
           haskey(adaptation_rules, :biomarker_enrichment)
            biomarker = get(adaptation_rules, :biomarker_name, nothing)
            if biomarker === nothing
                push!(warnings, "Enrichment design enabled but no biomarker specified. " *
                                "Using default biomarker name ':biomarker'.")
            end
        end

        # Unimplemented adaptive features - reject with clear error
        unimplemented_features = [
            (:platform_trial, "Platform/umbrella trials"),
            (:basket_trial, "Basket trials"),
            (:seamless_phase2_3, "Seamless Phase II/III designs"),
            (:mams, "Multi-arm multi-stage (MAMS) designs"),
            (:dose_response_adaptive, "Adaptive dose-response modeling"),
            (:rolling_admission, "Rolling admission adaptive designs")
        ]

        for (feature_key, feature_name) in unimplemented_features
            if haskey(adaptation_rules, feature_key) && adaptation_rules[feature_key] == true
                push!(errors, "$feature_name not yet implemented. " *
                              "Please use standard adaptive features (RAR, SSR, treatment selection, enrichment).")
            end
        end

        # Validate alpha spending function
        valid_spending = [:obrien_fleming, :pocock, :haybittle_peto, :linear, :none]
        if !(spec.design.alpha_spending in valid_spending)
            push!(warnings, "Unknown alpha spending function '$(spec.design.alpha_spending)'. " *
                            "Defaulting to O'Brien-Fleming. Valid options: $(valid_spending)")
        end

        # Check interim analyses are reasonable
        for ia in spec.design.interim_analyses
            if ia <= 0.1
                push!(warnings, "Interim analysis at $(ia*100)% information may have unreliable estimates.")
            end
        end
    end

    # Check BioequivalenceDesign requirements
    if spec.design isa BioequivalenceDesign
        # Check for appropriate endpoints
        has_pk_endpoint = any(ep -> ep isa PKEndpoint, spec.endpoints)
        if !has_pk_endpoint
            push!(errors, "BioequivalenceDesign requires at least one PKEndpoint.")
        end

        # Check for log_transform on PK endpoints
        for ep in spec.endpoints
            if ep isa PKEndpoint && !ep.log_transform
                push!(warnings, "PKEndpoint '$(ep.name)' should use log_transform=true for BE analysis.")
            end
        end

        # Check number of arms
        if length(spec.arms) != 2
            push!(errors, "BioequivalenceDesign requires exactly 2 treatment arms " *
                          "(Reference and Test), got $(length(spec.arms)).")
        end
    end

    # Check CrossoverDesign requirements
    if spec.design isa CrossoverDesign
        if length(spec.arms) < 2
            push!(errors, "CrossoverDesign requires at least 2 treatment arms, got $(length(spec.arms)).")
        end

        # Warn if n_subjects is odd (can't balance sequences)
        total_subjects = sum(arm.n_subjects for arm in spec.arms)
        if spec.design.n_sequences > 1 && total_subjects % spec.design.n_sequences != 0
            push!(warnings, "Total subjects ($total_subjects) not evenly divisible by sequences " *
                            "($(spec.design.n_sequences)). Sequences may be unbalanced.")
        end
    end

    # Check DoseEscalationDesign requirements
    if spec.design isa DoseEscalationDesign
        if length(spec.arms) != 1
            push!(warnings, "DoseEscalationDesign typically uses a single arm definition. " *
                            "Multiple arms may cause unexpected behavior.")
        end

        # Check that dose_levels are unique
        if length(spec.design.dose_levels) != length(unique(spec.design.dose_levels))
            push!(errors, "DoseEscalationDesign dose_levels must be unique.")
        end

        # Validate CRM skeleton
        if spec.design.escalation_rule isa CRM
            skeleton = spec.design.escalation_rule.skeleton
            if length(skeleton) != length(spec.design.dose_levels)
                push!(errors, "CRM skeleton length ($(length(skeleton))) must match " *
                              "dose_levels length ($(length(spec.design.dose_levels))).")
            end
            if !issorted(skeleton)
                push!(errors, "CRM skeleton probabilities must be monotonically increasing.")
            end
        end
    end

    # Check endpoints
    if isempty(spec.endpoints)
        push!(warnings, "No endpoints defined. Trial simulation will not produce endpoint analyses.")
    end

    # Check for duplicate endpoint names
    endpoint_names = [ep.name for ep in spec.endpoints]
    if length(endpoint_names) != length(unique(endpoint_names))
        push!(errors, "Endpoint names must be unique.")
    end

    # Check trial duration
    if spec.duration_days <= 0
        push!(errors, "Trial duration must be positive, got $(spec.duration_days) days.")
    end

    # Check sampling times
    if isempty(spec.pk_sampling_times)
        push!(warnings, "No PK sampling times defined. Using default times.")
    else
        if any(t -> t < 0, spec.pk_sampling_times)
            push!(errors, "PK sampling times cannot be negative.")
        end
        if any(t -> t > spec.duration_days * 24, spec.pk_sampling_times)
            push!(warnings, "Some PK sampling times exceed trial duration.")
        end
    end

    # Determine validity
    is_valid = isempty(errors)
    if strict
        is_valid = isempty(errors) && isempty(warnings)
    end

    return TrialValidationResult(is_valid, errors, warnings)
end

"""
    validate_trial_spec!(spec::TrialSpec; strict::Bool=false)

Validate a trial specification and throw an error if invalid.

# Arguments
- `spec`: Trial specification to validate
- `strict`: If true, treat warnings as errors

# Throws
- `ArgumentError`: If validation fails
"""
function validate_trial_spec!(spec::TrialSpec; strict::Bool = false)
    result = validate_trial_spec(spec; strict = strict)

    if !result.valid
        error_msg = "Trial specification validation failed:\n"
        for err in result.errors
            error_msg *= "  - ERROR: $err\n"
        end
        if strict
            for warn in result.warnings
                error_msg *= "  - WARNING: $warn\n"
            end
        end
        throw(ArgumentError(error_msg))
    end

    # Print warnings even if valid
    if !isempty(result.warnings) && !strict
        for warn in result.warnings
            @warn "Trial validation warning: $warn"
        end
    end

    return result
end
