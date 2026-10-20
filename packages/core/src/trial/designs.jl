# Trial Designs
# Convenience functions for creating study designs

export parallel_design, crossover_2x2, crossover_3x3, williams_design
export dose_escalation_3plus3, dose_escalation_mtpi, dose_escalation_crm
export adaptive_design, bioequivalence_design
export get_design_description

"""
    parallel_design(n_arms; kwargs...)

Create a parallel group study design.

# Arguments
- `n_arms::Int`: Number of treatment arms

# Keyword Arguments
- `randomization_ratio::Vector{Float64}`: Randomization ratio (default: equal)
- `stratification_factors::Vector{Symbol}`: Stratification factors

# Returns
- `ParallelDesign`: Parallel study design

# Example
```julia
# 2-arm parallel with 1:1 randomization
design = parallel_design(2)

# 3-arm with 2:1:1 randomization
design = parallel_design(3, randomization_ratio=[0.5, 0.25, 0.25])
```
"""
function parallel_design(n_arms::Int;
        randomization_ratio::Union{Nothing, Vector{Float64}} = nothing,
        stratification_factors::Vector{Symbol} = Symbol[])
    if randomization_ratio === nothing
        randomization_ratio = ones(n_arms) / n_arms
    end
    return ParallelDesign(n_arms;
        randomization_ratio = randomization_ratio,
        stratification_factors = stratification_factors)
end


"""
    crossover_2x2(; washout_duration=7.0)

Create a 2-period, 2-sequence crossover design (AB, BA).

# Keyword Arguments
- `washout_duration::Float64`: Washout period in days (default: 7.0)

# Returns
- `CrossoverDesign`: 2x2 crossover design

# Example
```julia
design = crossover_2x2(washout_duration=14.0)
```
"""
function crossover_2x2(; washout_duration::Float64 = 7.0)
    return CrossoverDesign(2, 2;
        washout_duration = washout_duration,
        sequence_assignments = [[1, 2], [2, 1]])
end


"""
    crossover_3x3(; washout_duration=7.0)

Create a 3-period, 3-sequence crossover design (Latin square).

# Keyword Arguments
- `washout_duration::Float64`: Washout period in days (default: 7.0)

# Returns
- `CrossoverDesign`: 3x3 crossover design

# Example
```julia
design = crossover_3x3(washout_duration=21.0)
```
"""
function crossover_3x3(; washout_duration::Float64 = 7.0)
    return CrossoverDesign(3, 3;
        washout_duration = washout_duration,
        sequence_assignments = [[1, 2, 3], [2, 3, 1], [3, 1, 2]])
end


"""
    williams_design(n_treatments; washout_duration=7.0)

Create a Williams design for crossover studies.

Williams designs are balanced for first-order carryover effects.

# Arguments
- `n_treatments::Int`: Number of treatments (2, 3, or 4)

# Keyword Arguments
- `washout_duration::Float64`: Washout period in days (default: 7.0)

# Returns
- `CrossoverDesign`: Williams crossover design

# Example
```julia
# 4-treatment Williams design
design = williams_design(4, washout_duration=14.0)
```
"""
function williams_design(n_treatments::Int; washout_duration::Float64 = 7.0)
    if n_treatments == 2
        sequences = [[1, 2], [2, 1]]
    elseif n_treatments == 3
        sequences = [[1, 2, 3], [2, 3, 1], [3, 1, 2],
                     [1, 3, 2], [3, 2, 1], [2, 1, 3]]
    elseif n_treatments == 4
        sequences = [[1, 2, 4, 3], [2, 3, 1, 4], [3, 4, 2, 1], [4, 1, 3, 2],
                     [1, 4, 2, 3], [2, 1, 3, 4], [3, 2, 4, 1], [4, 3, 1, 2]]
    else
        error("Williams design only supports 2, 3, or 4 treatments")
    end

    return CrossoverDesign(n_treatments, length(sequences);
        washout_duration = washout_duration,
        sequence_assignments = sequences)
end


"""
    dose_escalation_3plus3(dose_levels; kwargs...)

Create a 3+3 dose escalation design.

# Arguments
- `dose_levels::Vector{Float64}`: Available dose levels

# Keyword Arguments
- `starting_dose::Float64`: Starting dose (default: first level)
- `max_dlt_rate::Float64`: Maximum acceptable DLT rate (default: 0.33)
- `cohort_size::Int`: Cohort size (default: 3)
- `max_subjects::Int`: Maximum total subjects (default: 30)

# Returns
- `DoseEscalationDesign`: 3+3 dose escalation design

# Example
```julia
design = dose_escalation_3plus3([10.0, 25.0, 50.0, 100.0, 200.0])
```
"""
function dose_escalation_3plus3(dose_levels::Vector{Float64};
        starting_dose::Union{Float64, Nothing} = nothing,
        max_dlt_rate::Float64 = 0.33,
        cohort_size::Int = 3,
        max_subjects::Int = 30)
    rule = ThreePlusThree(max_dlt_rate = max_dlt_rate)
    return DoseEscalationDesign(dose_levels;
        starting_dose = starting_dose,
        escalation_rule = rule,
        cohort_size = cohort_size,
        max_subjects = max_subjects)
end


"""
    dose_escalation_mtpi(dose_levels; kwargs...)

Create an mTPI (modified Toxicity Probability Interval) dose escalation design.

# Arguments
- `dose_levels::Vector{Float64}`: Available dose levels

# Keyword Arguments
- `starting_dose::Float64`: Starting dose (default: first level)
- `target_dlt_rate::Float64`: Target DLT rate (default: 0.25)
- `equivalence_interval::Tuple{Float64, Float64}`: Equivalence interval (default: (0.20, 0.30))
- `cohort_size::Int`: Cohort size (default: 3)
- `max_subjects::Int`: Maximum total subjects (default: 30)

# Returns
- `DoseEscalationDesign`: mTPI dose escalation design

# Example
```julia
design = dose_escalation_mtpi([10.0, 25.0, 50.0, 100.0, 200.0],
                               target_dlt_rate=0.30)
```
"""
function dose_escalation_mtpi(dose_levels::Vector{Float64};
        starting_dose::Union{Float64, Nothing} = nothing,
        target_dlt_rate::Float64 = 0.25,
        equivalence_interval::Tuple{Float64, Float64} = (0.20, 0.30),
        cohort_size::Int = 3,
        max_subjects::Int = 30)
    rule = mTPI(target_dlt_rate = target_dlt_rate,
                equivalence_interval = equivalence_interval)
    return DoseEscalationDesign(dose_levels;
        starting_dose = starting_dose,
        escalation_rule = rule,
        cohort_size = cohort_size,
        max_subjects = max_subjects)
end


"""
    dose_escalation_crm(dose_levels; kwargs...)

Create a CRM (Continual Reassessment Method) dose escalation design.

# Arguments
- `dose_levels::Vector{Float64}`: Available dose levels

# Keyword Arguments
- `starting_dose::Float64`: Starting dose (default: first level)
- `target_dlt_rate::Float64`: Target DLT rate (default: 0.25)
- `skeleton::Vector{Float64}`: Prior toxicity probabilities
- `model::Symbol`: Dose-toxicity model (:logistic or :power)
- `cohort_size::Int`: Cohort size (default: 1)
- `max_subjects::Int`: Maximum total subjects (default: 30)

# Returns
- `DoseEscalationDesign`: CRM dose escalation design

# Example
```julia
design = dose_escalation_crm([10.0, 25.0, 50.0, 100.0, 200.0],
                              skeleton=[0.05, 0.12, 0.25, 0.40, 0.55])
```
"""
function dose_escalation_crm(dose_levels::Vector{Float64};
        starting_dose::Union{Float64, Nothing} = nothing,
        target_dlt_rate::Float64 = 0.25,
        skeleton::Union{Nothing, Vector{Float64}} = nothing,
        model::Symbol = :logistic,
        cohort_size::Int = 1,
        max_subjects::Int = 30)
    if skeleton === nothing
        # Generate default skeleton
        n = length(dose_levels)
        skeleton = [0.05 + 0.9 * (i - 1) / (n - 1) * target_dlt_rate for i in 1:n]
    end
    rule = CRM(target_dlt_rate = target_dlt_rate,
               skeleton = skeleton,
               model = model)
    return DoseEscalationDesign(dose_levels;
        starting_dose = starting_dose,
        escalation_rule = rule,
        cohort_size = cohort_size,
        max_subjects = max_subjects)
end


"""
    adaptive_design(base_design; kwargs...)

Create an adaptive trial design.

# Arguments
- `base_design::StudyDesignKind`: Underlying study design

# Keyword Arguments
- `interim_analyses::Vector{Float64}`: Information fractions for interim analyses
- `alpha_spending::Symbol`: Alpha spending function (:obrien_fleming, :pocock, :haybittle_peto)
- `futility_boundary::Float64`: Futility boundary (conditional power threshold)
- `sample_size_reestimation::Bool`: Allow sample size re-estimation

# Returns
- `AdaptiveDesign`: Adaptive trial design

# Example
```julia
base = parallel_design(2)
design = adaptive_design(base, interim_analyses=[0.5],
                         alpha_spending=:obrien_fleming)
```
"""
function adaptive_design(base_design::StudyDesignKind;
        interim_analyses::Vector{Float64} = [0.5],
        alpha_spending::Symbol = :obrien_fleming,
        futility_boundary::Float64 = 0.10,
        sample_size_reestimation::Bool = false)

    adaptation_rules = Dict{Symbol, Any}(
        :futility_boundary => futility_boundary,
        :sample_size_reestimation => sample_size_reestimation
    )

    return AdaptiveDesign(base_design;
        interim_analyses = interim_analyses,
        adaptation_rules = adaptation_rules,
        alpha_spending = alpha_spending)
end


"""
    bioequivalence_design(; kwargs...)

Create a bioequivalence study design.

# Keyword Arguments
- `n_periods::Int`: Number of periods (default: 2)
- `n_sequences::Int`: Number of sequences (default: 2)
- `washout_duration::Float64`: Washout period in days (default: 7.0)
- `bioequivalence_limits::Tuple{Float64, Float64}`: BE limits (default: (0.80, 1.25))
- `parameters::Vector{Symbol}`: PK parameters to assess (default: [:cmax, :auc_0_inf])
- `regulatory_guidance::Symbol`: Regulatory guidance (:fda, :ema)

# Returns
- `BioequivalenceDesign`: Bioequivalence study design

# Example
```julia
# FDA BE study
design = bioequivalence_design(regulatory_guidance=:fda)

# EMA BE study with highly variable drug widened limits
design = bioequivalence_design(
    bioequivalence_limits=(0.6984, 1.4319),
    regulatory_guidance=:ema
)
```
"""
function bioequivalence_design(;
        n_periods::Int = 2,
        n_sequences::Int = 2,
        washout_duration::Float64 = 7.0,
        bioequivalence_limits::Tuple{Float64, Float64} = (0.80, 1.25),
        parameters::Vector{Symbol} = [:cmax, :auc_0_inf],
        regulatory_guidance::Symbol = :fda)

    return BioequivalenceDesign(;
        n_periods = n_periods,
        n_sequences = n_sequences,
        washout_duration = washout_duration,
        bioequivalence_limits = bioequivalence_limits,
        parameters = parameters,
        regulatory_guidance = regulatory_guidance)
end


"""
    get_design_description(design::StudyDesignKind)

Get a human-readable description of a study design.

# Arguments
- `design::StudyDesignKind`: Study design

# Returns
- `String`: Design description
"""
function get_design_description(design::ParallelDesign)
    return "$(design.n_arms)-arm parallel group design"
end

function get_design_description(design::CrossoverDesign)
    return "$(design.n_periods)-period $(design.n_sequences)-sequence crossover design"
end

function get_design_description(design::DoseEscalationDesign)
    rule_name = typeof(design.escalation_rule) |> string |> x -> replace(x, "NeoPKPD." => "")
    return "Dose escalation ($(rule_name)) with $(length(design.dose_levels)) dose levels"
end

function get_design_description(design::AdaptiveDesign)
    base_desc = get_design_description(design.base_design)
    n_interim = length(design.interim_analyses)
    return "Adaptive $(base_desc) with $(n_interim) interim analysis(es)"
end

function get_design_description(design::BioequivalenceDesign)
    return "Bioequivalence $(get_design_description(design.crossover)) ($(design.regulatory_guidance))"
end
