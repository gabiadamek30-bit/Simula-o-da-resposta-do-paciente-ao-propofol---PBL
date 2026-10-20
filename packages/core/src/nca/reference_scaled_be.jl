# Reference-Scaled Average Bioequivalence (RSABE)
# FDA/EMA compliant implementation for highly variable drugs (HVDs)
#
# References:
# - FDA Guidance: Bioequivalence Studies for Highly Variable Drugs (2021)
# - EMA Guideline: Average Bioequivalence with Expanding Limits (ABEL)
# - FDA Draft Guidance: Statistical Approaches for Reference-Scaled Average BE
# - Haidar et al. (2008): Bioequivalence Approaches for HVDs

export RSABEConfig, RSABEResult, ABELResult
export ReplicateDesign, PartialReplicate3x3, FullReplicate2x4, FullReplicate2x3
export RegulatoryGuidance, FDAGuidance, EMAGuidance, HealthCanadaGuidance
export rsabe_analysis, abel_analysis, reference_scaled_be
export compute_swr, compute_within_subject_variance_reference
export rsabe_criterion, abel_scaled_limits
export extract_replicate_data, ReplicateData
export replicate_anova

# =============================================================================
# Regulatory Guidance Types
# =============================================================================

"""
Abstract type for regulatory guidance selection.
Different agencies have different approaches to reference-scaled BE.
"""
abstract type RegulatoryGuidance end

"""
FDA guidance for highly variable drugs.

Uses RSABE criterion: (μT - μR)² ≤ θ × σ²WR
- θ = (ln(1.25)/0.25)² ≈ 0.8926
- Point estimate constraint: GMR within 80-125%
- Scaling threshold: σWR > 0.294 (CVw > 30%)
"""
struct FDAGuidance <: RegulatoryGuidance
    theta::Float64                    # Scaling factor (default: 0.8926)
    point_estimate_lower::Float64     # GMR lower bound (default: 0.80)
    point_estimate_upper::Float64     # GMR upper bound (default: 1.25)
    swr_threshold::Float64            # σWR threshold for scaling (default: 0.294)

    function FDAGuidance(;
        theta::Float64 = (log(1.25) / 0.25)^2,  # ≈ 0.8926
        point_estimate_lower::Float64 = 0.80,
        point_estimate_upper::Float64 = 1.25,
        swr_threshold::Float64 = 0.294
    )
        new(theta, point_estimate_lower, point_estimate_upper, swr_threshold)
    end
end

"""
EMA guidance using Average Bioequivalence with Expanding Limits (ABEL).

Scaled limits: [0.80, 1.25] × exp(±k × σWR)
- k = 0.760 (regulatory constant)
- Maximum expansion: 69.84% - 143.19%
- Scaling threshold: CVw > 30%
"""
struct EMAGuidance <: RegulatoryGuidance
    k_scaling::Float64                # Expansion factor (default: 0.760)
    limit_lower_min::Float64          # Minimum lower limit (default: 0.6984)
    limit_upper_max::Float64          # Maximum upper limit (default: 1.4319)
    cv_threshold::Float64             # CV% threshold for scaling (default: 30.0)

    function EMAGuidance(;
        k_scaling::Float64 = 0.760,
        limit_lower_min::Float64 = 0.6984,
        limit_upper_max::Float64 = 1.4319,
        cv_threshold::Float64 = 30.0
    )
        new(k_scaling, limit_lower_min, limit_upper_max, cv_threshold)
    end
end

"""
Health Canada guidance (similar to FDA).
"""
struct HealthCanadaGuidance <: RegulatoryGuidance
    theta::Float64
    point_estimate_lower::Float64
    point_estimate_upper::Float64
    swr_threshold::Float64

    function HealthCanadaGuidance(;
        theta::Float64 = (log(1.25) / 0.25)^2,
        point_estimate_lower::Float64 = 0.80,
        point_estimate_upper::Float64 = 1.25,
        swr_threshold::Float64 = 0.294
    )
        new(theta, point_estimate_lower, point_estimate_upper, swr_threshold)
    end
end

# =============================================================================
# Replicate Design Types
# =============================================================================

"""
Abstract type for replicate crossover designs used in RSABE studies.
"""
abstract type ReplicateDesign end

"""
3-period partial replicate design (TRR/RTR/RRT).

- 3 sequences, 3 periods
- Each subject receives reference twice, test once
- Allows estimation of within-subject reference variability
"""
struct PartialReplicate3x3 <: ReplicateDesign
    sequences::Vector{String}         # e.g., ["TRR", "RTR", "RRT"]

    function PartialReplicate3x3(;
        sequences::Vector{String} = ["TRR", "RTR", "RRT"]
    )
        new(sequences)
    end
end

"""
4-period full replicate design (TRTR/RTRT).

- 2 sequences, 4 periods
- Each subject receives both formulations twice
- Provides best precision for RSABE
"""
struct FullReplicate2x4 <: ReplicateDesign
    sequences::Vector{String}         # e.g., ["TRTR", "RTRT"]

    function FullReplicate2x4(;
        sequences::Vector{String} = ["TRTR", "RTRT"]
    )
        new(sequences)
    end
end

"""
3-period full replicate design (TRT/RTR).

- 2 sequences, 3 periods
- Similar to partial replicate but balanced
"""
struct FullReplicate2x3 <: ReplicateDesign
    sequences::Vector{String}

    function FullReplicate2x3(;
        sequences::Vector{String} = ["TRT", "RTR"]
    )
        new(sequences)
    end
end

# =============================================================================
# Replicate Data Extraction
# =============================================================================

"""
Container for extracted replicate design data.

Holds the observations organized by subject and formulation.
"""
struct ReplicateData
    subject_ids::Vector{String}       # Subject identifiers
    sequences::Vector{String}         # Sequence assignment per subject
    test_obs::Vector{Vector{Float64}} # Test observations per subject
    ref_obs::Vector{Vector{Float64}}  # Reference observations per subject
    n_subjects::Int
    n_test_per_subject::Int           # Number of test doses per subject
    n_ref_per_subject::Int            # Number of reference doses per subject
end

"""
    extract_replicate_data(data, design; subject_col, sequence_col,
                           formulation_col, period_col, value_col)

Extract and organize replicate design data for RSABE analysis.

# Arguments
- `data`: Data table with columns for subject, sequence, formulation, period, value
- `design::ReplicateDesign`: Study design specification

# Keyword Arguments
- `subject_col::Symbol`: Column name for subject ID (default: :subject)
- `sequence_col::Symbol`: Column name for sequence (default: :sequence)
- `formulation_col::Symbol`: Column name for formulation T/R (default: :formulation)
- `period_col::Symbol`: Column name for period (default: :period)
- `value_col::Symbol`: Column name for PK parameter value (default: :value)

# Returns
- `ReplicateData`: Organized replicate data structure

# Example
```julia
# For a 2x4 full replicate design
data = extract_replicate_data(pk_data, FullReplicate2x4();
    subject_col=:SUBJID, formulation_col=:TRT, value_col=:CMAX)
```
"""
function extract_replicate_data(
    subjects::Vector{String},
    sequences::Vector{String},
    formulations::Vector{String},
    periods::Vector{Int},
    values::Vector{Float64},
    design::ReplicateDesign
)
    @assert length(subjects) == length(sequences) == length(formulations) ==
            length(periods) == length(values) "All vectors must have same length"

    # Get unique subjects
    unique_subjects = unique(subjects)
    n_subjects = length(unique_subjects)

    # Organize data by subject
    test_obs = Vector{Vector{Float64}}(undef, n_subjects)
    ref_obs = Vector{Vector{Float64}}(undef, n_subjects)
    subj_sequences = Vector{String}(undef, n_subjects)

    for (i, subj) in enumerate(unique_subjects)
        # Find all observations for this subject
        mask = subjects .== subj
        subj_form = formulations[mask]
        subj_vals = values[mask]
        subj_periods = periods[mask]

        # Get sequence for this subject
        subj_sequences[i] = sequences[findfirst(mask)]

        # Sort by period
        order = sortperm(subj_periods)
        subj_form = subj_form[order]
        subj_vals = subj_vals[order]

        # Separate test and reference
        test_mask = subj_form .== "T" .|| subj_form .== "Test" .|| subj_form .== "test"
        ref_mask = subj_form .== "R" .|| subj_form .== "Reference" .|| subj_form .== "ref"

        test_obs[i] = subj_vals[test_mask]
        ref_obs[i] = subj_vals[ref_mask]
    end

    # Determine expected observations per subject based on design
    n_test = _n_test_per_subject(design)
    n_ref = _n_ref_per_subject(design)

    return ReplicateData(
        unique_subjects,
        subj_sequences,
        test_obs,
        ref_obs,
        n_subjects,
        n_test,
        n_ref
    )
end

"""
Convenience function to extract replicate data from arrays.
"""
function extract_replicate_data(
    test_values::Matrix{Float64},      # Subjects × Test observations
    ref_values::Matrix{Float64},       # Subjects × Reference observations
    design::ReplicateDesign
)
    n_subjects = size(test_values, 1)
    @assert size(ref_values, 1) == n_subjects "Test and reference must have same number of subjects"

    subjects = ["S" * string(i) for i in 1:n_subjects]
    sequences = fill("", n_subjects)

    test_obs = [test_values[i, :] for i in 1:n_subjects]
    ref_obs = [ref_values[i, :] for i in 1:n_subjects]

    n_test = size(test_values, 2)
    n_ref = size(ref_values, 2)

    return ReplicateData(subjects, sequences, test_obs, ref_obs, n_subjects, n_test, n_ref)
end

function _n_test_per_subject(design::PartialReplicate3x3)
    return 1  # T once per subject
end

function _n_test_per_subject(design::FullReplicate2x4)
    return 2  # T twice per subject
end

function _n_test_per_subject(design::FullReplicate2x3)
    return 1  # Variable
end

function _n_ref_per_subject(design::PartialReplicate3x3)
    return 2  # R twice per subject
end

function _n_ref_per_subject(design::FullReplicate2x4)
    return 2  # R twice per subject
end

function _n_ref_per_subject(design::FullReplicate2x3)
    return 2  # R twice per subject
end

# =============================================================================
# Within-Subject Reference Variance (σ²WR)
# =============================================================================

"""
    compute_within_subject_variance_reference(ref_obs) -> (swr_squared, swr, cv_wr)

Compute within-subject variance for reference formulation from replicate data.

This is the key estimate for reference-scaled bioequivalence.

# Arguments
- `ref_obs::Vector{Vector{Float64}}`: Reference observations per subject
  (each subject should have 2+ observations)

# Returns
- `NamedTuple`: (swr_squared, swr, cv_wr, df)
  - `swr_squared::Float64`: Within-subject variance on log scale (σ²WR)
  - `swr::Float64`: Within-subject standard deviation (σWR)
  - `cv_wr::Float64`: Within-subject CV for reference (%)
  - `df::Int`: Degrees of freedom

# Method
Uses ANOVA-type estimation:
σ²WR = Σᵢ Σⱼ (Yᵢⱼ - Ȳᵢ)² / Σᵢ(nᵢ - 1)

where Yᵢⱼ is log-transformed observation j for subject i.
"""
function compute_within_subject_variance_reference(
    ref_obs::Vector{Vector{Float64}}
)
    n_subjects = length(ref_obs)
    @assert n_subjects >= 2 "Need at least 2 subjects"

    # Log-transform observations
    log_ref = [log.(obs) for obs in ref_obs]

    # Compute within-subject sum of squares
    ss_within = 0.0
    df_total = 0

    for i in 1:n_subjects
        n_i = length(log_ref[i])
        if n_i >= 2
            mean_i = sum(log_ref[i]) / n_i
            ss_within += sum((log_ref[i] .- mean_i).^2)
            df_total += n_i - 1
        end
    end

    @assert df_total > 0 "Need at least one subject with 2+ reference observations"

    # Within-subject variance
    swr_squared = ss_within / df_total
    swr = sqrt(swr_squared)

    # Convert to CV (%)
    cv_wr = sqrt(exp(swr_squared) - 1.0) * 100.0

    return (swr_squared=swr_squared, swr=swr, cv_wr=cv_wr, df=df_total)
end

"""
    compute_swr(ref_obs) -> Float64

Convenience function to compute σWR (within-subject SD for reference).
"""
function compute_swr(ref_obs::Vector{Vector{Float64}})
    result = compute_within_subject_variance_reference(ref_obs)
    return result.swr
end

"""
Compute σWR from paired reference values (R1, R2 per subject).
"""
function compute_swr(ref1::Vector{Float64}, ref2::Vector{Float64})
    @assert length(ref1) == length(ref2) "R1 and R2 must have same length"

    n = length(ref1)
    @assert n >= 2 "Need at least 2 subjects"

    # Combine into per-subject observations
    ref_obs = [[ref1[i], ref2[i]] for i in 1:n]

    return compute_swr(ref_obs)
end

# =============================================================================
# FDA RSABE Criterion
# =============================================================================

"""
    rsabe_criterion(gmr, swr, guidance) -> (criterion_met, scaled_limit, criterion_value)

Evaluate FDA RSABE criterion for highly variable drugs.

# FDA Criterion
The reference-scaled criterion is:
  (μT - μR)² - θ × σ²WR ≤ 0

Equivalently, using the linearized form:
  |μT - μR| ≤ √(θ × σ²WR)

# Arguments
- `log_diff::Float64`: Mean log difference (μT - μR on log scale)
- `swr::Float64`: Within-subject reference SD (σWR)
- `guidance::FDAGuidance`: FDA regulatory parameters

# Returns
- `NamedTuple`: (criterion_met, scaled_limit, criterion_value, use_scaled)
  - `criterion_met::Bool`: Whether criterion is satisfied
  - `scaled_limit::Float64`: Scaled BE limit (ratio scale)
  - `criterion_value::Float64`: Value of RSABE criterion
  - `use_scaled::Bool`: Whether scaling was applied (σWR > threshold)
"""
function rsabe_criterion(
    log_diff::Float64,        # μT - μR on log scale
    swr::Float64,             # σWR
    guidance::FDAGuidance
)
    # Check if scaling should be applied
    use_scaled = swr > guidance.swr_threshold

    if use_scaled
        # Scaled criterion: (log_diff)² ≤ θ × σ²WR
        criterion_value = log_diff^2 - guidance.theta * swr^2
        criterion_met = criterion_value <= 0.0

        # Scaled limit on ratio scale
        scaled_limit_log = sqrt(guidance.theta) * swr
        scaled_limit_upper = exp(scaled_limit_log)
        scaled_limit_lower = exp(-scaled_limit_log)
    else
        # Standard ABE limits
        criterion_value = NaN  # Not applicable
        scaled_limit_upper = guidance.point_estimate_upper
        scaled_limit_lower = guidance.point_estimate_lower

        # Check standard criterion
        gmr = exp(log_diff)
        criterion_met = (gmr >= scaled_limit_lower) && (gmr <= scaled_limit_upper)
    end

    return (
        criterion_met = criterion_met,
        scaled_limit_lower = scaled_limit_lower,
        scaled_limit_upper = scaled_limit_upper,
        criterion_value = criterion_value,
        use_scaled = use_scaled
    )
end

"""
    rsabe_upper_bound(log_diff, swr, se_diff, se_swr, n; alpha=0.05)

Compute upper 95% confidence bound for RSABE criterion.

Uses the linearized approach recommended by FDA:
Upper bound of (μT - μR)² - θ × σ²WR

# Arguments
- `log_diff::Float64`: Point estimate of log difference
- `swr::Float64`: Point estimate of σWR
- `se_diff::Float64`: Standard error of log difference
- `se_swr::Float64`: Standard error of σWR estimate
- `n::Int`: Number of subjects

# Returns
- `Float64`: Upper 95% confidence bound for criterion
"""
function rsabe_upper_bound(
    log_diff::Float64,
    swr::Float64,
    se_diff::Float64,
    se_swr::Float64,
    n::Int,
    theta::Float64;
    alpha::Float64 = 0.05
)
    # Linearized criterion: C = (μT - μR)² - θ × σ²WR
    # Variance: Var(C) ≈ 4(μT-μR)² × Var(μT-μR) + θ² × 4σ²WR × Var(σ²WR)

    # Point estimate
    C_hat = log_diff^2 - theta * swr^2

    # Approximate variance using delta method
    # Var(log_diff²) ≈ 4 × log_diff² × se_diff²
    var_diff_sq = 4.0 * log_diff^2 * se_diff^2

    # Var(σ²WR) from se_swr (assuming chi-square distribution)
    # se_swr² ≈ σ²WR × 2/(df)
    var_swr_sq = 4.0 * swr^2 * se_swr^2

    # Total variance (assuming independence)
    var_C = var_diff_sq + theta^2 * var_swr_sq
    se_C = sqrt(max(var_C, 0.0))

    # Degrees of freedom (Satterthwaite approximation)
    df = max(n - 2, 1)

    # Upper bound
    t_crit = _t_critical(df, 2 * alpha)  # One-sided
    upper_bound = C_hat + t_crit * se_C

    return upper_bound
end

# =============================================================================
# EMA ABEL (Average Bioequivalence with Expanding Limits)
# =============================================================================

"""
    abel_scaled_limits(swr, guidance) -> (lower, upper, is_scaled)

Calculate EMA ABEL expanded bioequivalence limits.

# EMA Scaling Formula
When CV > 30% (σWR > ~0.294):
  Lower limit = 0.80 × exp(-k × σWR)
  Upper limit = 1.25 × exp(+k × σWR)

With k = 0.760 and caps at 69.84% - 143.19%

# Arguments
- `swr::Float64`: Within-subject reference SD (σWR)
- `guidance::EMAGuidance`: EMA regulatory parameters

# Returns
- `NamedTuple`: (lower, upper, is_scaled)
"""
function abel_scaled_limits(swr::Float64, guidance::EMAGuidance)
    # Convert σWR to CV for threshold check
    cv_wr = sqrt(exp(swr^2) - 1.0) * 100.0

    if cv_wr > guidance.cv_threshold
        # Apply scaling
        k = guidance.k_scaling

        # Expanded limits
        lower = 0.80 * exp(-k * swr)
        upper = 1.25 * exp(k * swr)

        # Apply caps
        lower = max(lower, guidance.limit_lower_min)
        upper = min(upper, guidance.limit_upper_max)

        is_scaled = true
    else
        # Standard limits
        lower = 0.80
        upper = 1.25
        is_scaled = false
    end

    return (lower=lower, upper=upper, is_scaled=is_scaled)
end

# =============================================================================
# RSABE Configuration
# =============================================================================

"""
Configuration for reference-scaled average bioequivalence analysis.
"""
struct RSABEConfig
    guidance::RegulatoryGuidance
    design::ReplicateDesign
    parameter::Symbol                  # :cmax, :auc_0_inf, etc.
    alpha::Float64                     # Significance level (default: 0.05)

    function RSABEConfig(;
        guidance::RegulatoryGuidance = FDAGuidance(),
        design::ReplicateDesign = FullReplicate2x4(),
        parameter::Symbol = :cmax,
        alpha::Float64 = 0.05
    )
        new(guidance, design, parameter, alpha)
    end
end

# =============================================================================
# RSABE Result Types
# =============================================================================

"""
Result of FDA Reference-Scaled Average Bioequivalence analysis.

# Fields
- `parameter::Symbol`: Parameter analyzed
- `n_subjects::Int`: Number of subjects
- `gmr::Float64`: Geometric mean ratio (point estimate)
- `log_diff::Float64`: Mean log difference (μT - μR)
- `se_diff::Float64`: Standard error of log difference
- `swr::Float64`: Within-subject reference SD (σWR)
- `swr_squared::Float64`: Within-subject reference variance (σ²WR)
- `cv_wr::Float64`: Within-subject reference CV (%)
- `use_scaled::Bool`: Whether scaled criterion was used
- `rsabe_criterion::Float64`: Value of RSABE criterion (if scaled)
- `rsabe_upper_bound::Float64`: Upper 95% CI for criterion
- `point_estimate_pass::Bool`: GMR within 80-125%
- `scaled_criterion_pass::Bool`: RSABE criterion satisfied
- `be_conclusion::Symbol`: Overall BE conclusion
- `scaled_limit_lower::Float64`: Lower BE limit used
- `scaled_limit_upper::Float64`: Upper BE limit used
"""
struct RSABEResult
    parameter::Symbol
    n_subjects::Int
    gmr::Float64
    log_diff::Float64
    se_diff::Float64
    swr::Float64
    swr_squared::Float64
    cv_wr::Float64
    use_scaled::Bool
    rsabe_criterion::Float64
    rsabe_upper_bound::Float64
    point_estimate_pass::Bool
    scaled_criterion_pass::Bool
    be_conclusion::Symbol
    scaled_limit_lower::Float64
    scaled_limit_upper::Float64
end

"""
Result of EMA Average Bioequivalence with Expanding Limits analysis.

# Fields
- `parameter::Symbol`: Parameter analyzed
- `n_subjects::Int`: Number of subjects
- `gmr::Float64`: Geometric mean ratio
- `ci_lower::Float64`: 90% CI lower bound
- `ci_upper::Float64`: 90% CI upper bound
- `swr::Float64`: Within-subject reference SD
- `cv_wr::Float64`: Within-subject reference CV (%)
- `use_scaled::Bool`: Whether expanded limits were used
- `be_limit_lower::Float64`: Lower BE limit (possibly expanded)
- `be_limit_upper::Float64`: Upper BE limit (possibly expanded)
- `be_conclusion::Symbol`: Overall BE conclusion
"""
struct ABELResult
    parameter::Symbol
    n_subjects::Int
    gmr::Float64
    ci_lower::Float64
    ci_upper::Float64
    swr::Float64
    cv_wr::Float64
    use_scaled::Bool
    be_limit_lower::Float64
    be_limit_upper::Float64
    be_conclusion::Symbol
end

# =============================================================================
# Main Analysis Functions
# =============================================================================

"""
    rsabe_analysis(data, config) -> RSABEResult

Perform FDA Reference-Scaled Average Bioequivalence analysis.

# Arguments
- `data::ReplicateData`: Replicate design data
- `config::RSABEConfig`: Analysis configuration (with FDAGuidance)

# Returns
- `RSABEResult`: Complete RSABE analysis result

# Example
```julia
config = RSABEConfig(
    guidance = FDAGuidance(),
    design = FullReplicate2x4(),
    parameter = :cmax
)

result = rsabe_analysis(data, config)

if result.be_conclusion == :bioequivalent
    println("BE demonstrated using RSABE")
end
```
"""
function rsabe_analysis(
    data::ReplicateData,
    config::RSABEConfig
)
    guidance = config.guidance
    @assert guidance isa FDAGuidance || guidance isa HealthCanadaGuidance "RSABE requires FDA or Health Canada guidance"

    n = data.n_subjects

    # Compute within-subject reference variance
    swr_result = compute_within_subject_variance_reference(data.ref_obs)
    swr = swr_result.swr
    swr_squared = swr_result.swr_squared
    cv_wr = swr_result.cv_wr

    # Compute GMR and log difference
    # For each subject, compute mean test and mean reference
    test_means = [geometric_mean(obs) for obs in data.test_obs]
    ref_means = [geometric_mean(obs) for obs in data.ref_obs]

    # Log-transform
    log_test = log.(test_means)
    log_ref = log.(ref_means)

    # Mean log difference
    log_diffs = log_test .- log_ref
    log_diff = sum(log_diffs) / n

    # Standard error of log difference
    ss_diff = sum((log_diffs .- log_diff).^2)
    mse_diff = ss_diff / (n - 1)
    se_diff = sqrt(mse_diff / n)

    # SE for σWR (approximate)
    se_swr = swr / sqrt(2 * swr_result.df)

    # GMR
    gmr = exp(log_diff)

    # Evaluate RSABE criterion
    criterion_result = rsabe_criterion(log_diff, swr, guidance)
    use_scaled = criterion_result.use_scaled

    # Point estimate constraint
    point_estimate_pass = (gmr >= guidance.point_estimate_lower) &&
                          (gmr <= guidance.point_estimate_upper)

    # Scaled criterion
    if use_scaled
        # Compute upper bound for scaled criterion
        upper_bound = rsabe_upper_bound(log_diff, swr, se_diff, se_swr, n, guidance.theta)
        scaled_criterion_pass = upper_bound <= 0.0
        rsabe_crit = criterion_result.criterion_value
    else
        # Standard ABE - use 90% CI approach
        t_crit = _t_critical(n - 1, 0.10)
        ci_lower = exp(log_diff - t_crit * se_diff)
        ci_upper = exp(log_diff + t_crit * se_diff)
        scaled_criterion_pass = (ci_lower >= 0.80) && (ci_upper <= 1.25)
        rsabe_crit = NaN
        upper_bound = NaN
    end

    # Overall conclusion
    if use_scaled
        # Scaled: need both point estimate and criterion
        be_conclusion = (point_estimate_pass && scaled_criterion_pass) ?
                        :bioequivalent : :not_bioequivalent
    else
        # Standard: just need CI within limits
        be_conclusion = scaled_criterion_pass ? :bioequivalent : :not_bioequivalent
    end

    return RSABEResult(
        config.parameter,
        n,
        gmr,
        log_diff,
        se_diff,
        swr,
        swr_squared,
        cv_wr,
        use_scaled,
        rsabe_crit,
        upper_bound,
        point_estimate_pass,
        scaled_criterion_pass,
        be_conclusion,
        criterion_result.scaled_limit_lower,
        criterion_result.scaled_limit_upper
    )
end

"""
    abel_analysis(data, config) -> ABELResult

Perform EMA Average Bioequivalence with Expanding Limits analysis.

# Arguments
- `data::ReplicateData`: Replicate design data
- `config::RSABEConfig`: Analysis configuration (with EMAGuidance)

# Returns
- `ABELResult`: Complete ABEL analysis result
"""
function abel_analysis(
    data::ReplicateData,
    config::RSABEConfig
)
    guidance = config.guidance
    @assert guidance isa EMAGuidance "ABEL requires EMA guidance"

    n = data.n_subjects

    # Compute within-subject reference variance
    swr_result = compute_within_subject_variance_reference(data.ref_obs)
    swr = swr_result.swr
    cv_wr = swr_result.cv_wr

    # Get scaled limits
    limits = abel_scaled_limits(swr, guidance)

    # Compute GMR and 90% CI
    test_means = [geometric_mean(obs) for obs in data.test_obs]
    ref_means = [geometric_mean(obs) for obs in data.ref_obs]

    log_test = log.(test_means)
    log_ref = log.(ref_means)

    log_diffs = log_test .- log_ref
    log_diff = sum(log_diffs) / n

    ss_diff = sum((log_diffs .- log_diff).^2)
    mse_diff = ss_diff / (n - 1)
    se_diff = sqrt(mse_diff / n)

    # 90% CI
    t_crit = _t_critical(n - 1, 0.10)
    ci_lower = exp(log_diff - t_crit * se_diff)
    ci_upper = exp(log_diff + t_crit * se_diff)

    gmr = exp(log_diff)

    # BE conclusion: 90% CI entirely within scaled limits
    be_conclusion = (ci_lower >= limits.lower) && (ci_upper <= limits.upper) ?
                    :bioequivalent : :not_bioequivalent

    return ABELResult(
        config.parameter,
        n,
        gmr,
        ci_lower,
        ci_upper,
        swr,
        cv_wr,
        limits.is_scaled,
        limits.lower,
        limits.upper,
        be_conclusion
    )
end

"""
    reference_scaled_be(data, config) -> Union{RSABEResult, ABELResult}

Unified entry point for reference-scaled bioequivalence analysis.
Automatically dispatches to FDA RSABE or EMA ABEL based on guidance.
"""
function reference_scaled_be(data::ReplicateData, config::RSABEConfig)
    if config.guidance isa EMAGuidance
        return abel_analysis(data, config)
    else
        return rsabe_analysis(data, config)
    end
end

"""
    reference_scaled_be(test, ref, config) -> Union{RSABEResult, ABELResult}

Convenience function for simple replicate data.

# Arguments
- `test::Matrix{Float64}`: Test observations (subjects × observations)
- `ref::Matrix{Float64}`: Reference observations (subjects × observations)
- `config::RSABEConfig`: Analysis configuration
"""
function reference_scaled_be(
    test::Matrix{Float64},
    ref::Matrix{Float64},
    config::RSABEConfig
)
    data = extract_replicate_data(test, ref, config.design)
    return reference_scaled_be(data, config)
end

# =============================================================================
# Mixed-Effects ANOVA for Replicate Designs
# =============================================================================

"""
    replicate_anova(data, design) -> NamedTuple

Perform mixed-effects ANOVA for replicate crossover designs.

Model: Y_ijk = μ + SEQ_i + SUBJ(SEQ)_ij + PER_k + FORM_l + ε_ijkl

Where:
- SEQ = sequence effect (fixed)
- SUBJ(SEQ) = subject nested within sequence (random)
- PER = period effect (fixed)
- FORM = formulation effect (fixed)

# Returns
- `NamedTuple`: ANOVA results including:
  - `mse_within`: Within-subject mean square error
  - `swr_squared`: Within-subject reference variance
  - `formulation_effect`: Estimated formulation effect
  - `se_formulation`: SE of formulation effect
"""
function replicate_anova(
    data::ReplicateData,
    design::ReplicateDesign
)
    n = data.n_subjects

    # Build design matrix and response vector
    # For simplicity, use method of moments estimation

    # Compute overall means
    all_test = vcat([obs for obs in data.test_obs]...)
    all_ref = vcat([obs for obs in data.ref_obs]...)

    log_all_test = log.(all_test)
    log_all_ref = log.(all_ref)

    mean_log_test = sum(log_all_test) / length(log_all_test)
    mean_log_ref = sum(log_all_ref) / length(log_all_ref)

    # Formulation effect
    formulation_effect = mean_log_test - mean_log_ref

    # Within-subject variance for reference
    swr_result = compute_within_subject_variance_reference(data.ref_obs)

    # Within-subject variance for test (if available)
    if data.n_test_per_subject >= 2
        swt_result = compute_within_subject_variance_reference(data.test_obs)
        swt = swt_result.swr
    else
        swt = NaN
    end

    # Total within-subject variance (pooled)
    # For reference-scaled BE, we primarily use σWR

    # SE of formulation effect (approximate)
    # SE ≈ sqrt(2 × MSE / n) for balanced design
    mse_within = swr_result.swr_squared
    se_formulation = sqrt(2.0 * mse_within / n)

    return (
        mse_within = mse_within,
        swr_squared = swr_result.swr_squared,
        swr = swr_result.swr,
        swt = swt,
        cv_wr = swr_result.cv_wr,
        formulation_effect = formulation_effect,
        se_formulation = se_formulation,
        df_within = swr_result.df
    )
end

# =============================================================================
# Utility Functions (use existing from bioequivalence.jl)
# =============================================================================

# These functions are defined in bioequivalence.jl and will be available:
# - _t_critical(df, alpha)
# - _t_cdf(t, df)
# - _t_quantile(p, df)
# - _lgamma(x)
# - geometric_mean(values)
