# NCA Specification Types
# FDA/EMA compliant Non-Compartmental Analysis

export NCAConfig, LambdaZResult, NCAResult, BioequivalenceResult, TOSTResult
export NCAMethod, LinearMethod, LogLinearMethod, LinLogMixedMethod, LinLogMixed
export BLQHandling, BLQZero, BLQMissing, BLQLLOQHalf

# =============================================================================
# NCA Method Types
# =============================================================================

"""
Abstract type for AUC calculation methods.
"""
abstract type NCAMethod end

"""
Linear trapezoidal method for all intervals.
Best for: Data with no clear log-linear terminal phase.
"""
struct LinearMethod <: NCAMethod end

"""
Log-linear trapezoidal method for all intervals.
Best for: Data where concentrations decline log-linearly.
"""
struct LogLinearMethod <: NCAMethod end

"""
Linear up / Log-linear down method (Lin-Log Mixed).
Uses linear trapezoidal during absorption (ascending) and
log-linear trapezoidal during elimination (descending).
Best for: Most PK data, FDA/EMA preferred.
"""
struct LinLogMixedMethod <: NCAMethod end

# Alias for convenience
const LinLogMixed = LinLogMixedMethod

# =============================================================================
# BLQ Handling Types
# =============================================================================

"""
Abstract type for below limit of quantification (BLQ) handling.
"""
abstract type BLQHandling end

"""
Treat BLQ values as zero.
"""
struct BLQZero <: BLQHandling end

"""
Treat BLQ values as missing (exclude from calculations).
"""
struct BLQMissing <: BLQHandling end

"""
Treat BLQ values as LLOQ/2.
"""
struct BLQLLOQHalf <: BLQHandling end

# =============================================================================
# NCA Configuration
# =============================================================================

"""
    NCAConfig

Configuration for NCA analysis following FDA/EMA guidance.

# Fields
- `method::NCAMethod`: AUC calculation method (default: LinLogMixedMethod)
- `lambda_z_min_points::Int`: Minimum points for lambda_z regression (default: 3)
- `lambda_z_r2_threshold::Float64`: Minimum R² for acceptable lambda_z (default: 0.9)
- `lambda_z_span_half_lives::Float64`: Minimum span in half-lives (default: 2.0)
- `extrapolation_max_pct::Float64`: Maximum allowed AUC extrapolation % (default: 20.0)
- `significant_digits::Int`: Significant digits for regulatory reporting (default: 3)
- `blq_handling::BLQHandling`: How to handle BLQ values (default: BLQZero)
- `lloq::Union{Float64,Nothing}`: Lower limit of quantification (default: nothing)
"""
struct NCAConfig
    method::NCAMethod
    lambda_z_min_points::Int
    lambda_z_r2_threshold::Float64
    lambda_z_span_half_lives::Float64
    extrapolation_max_pct::Float64
    significant_digits::Int
    blq_handling::BLQHandling
    lloq::Union{Float64,Nothing}
end

function NCAConfig(;
    method::NCAMethod = LinLogMixedMethod(),
    lambda_z_min_points::Int = 3,
    lambda_z_r2_threshold::Float64 = 0.9,
    lambda_z_span_half_lives::Float64 = 2.0,
    extrapolation_max_pct::Float64 = 20.0,
    significant_digits::Int = 3,
    blq_handling::BLQHandling = BLQZero(),
    lloq::Union{Float64,Nothing} = nothing
)
    @assert lambda_z_min_points >= 2 "lambda_z_min_points must be at least 2"
    @assert 0.0 < lambda_z_r2_threshold <= 1.0 "lambda_z_r2_threshold must be in (0, 1]"
    @assert lambda_z_span_half_lives > 0.0 "lambda_z_span_half_lives must be positive"
    @assert 0.0 < extrapolation_max_pct <= 100.0 "extrapolation_max_pct must be in (0, 100]"
    @assert significant_digits >= 1 "significant_digits must be at least 1"

    return NCAConfig(
        method,
        lambda_z_min_points,
        lambda_z_r2_threshold,
        lambda_z_span_half_lives,
        extrapolation_max_pct,
        significant_digits,
        blq_handling,
        lloq
    )
end

# =============================================================================
# Lambda-z Result
# =============================================================================

"""
    LambdaZResult

Result of terminal elimination rate constant (λz) estimation.

# Fields
- `lambda_z::Union{Float64,Nothing}`: Terminal elimination rate constant (1/time)
- `t_half::Union{Float64,Nothing}`: Terminal half-life = ln(2)/λz
- `r_squared::Union{Float64,Nothing}`: R² of log-linear regression
- `adjusted_r_squared::Union{Float64,Nothing}`: Adjusted R²
- `intercept::Union{Float64,Nothing}`: Y-intercept of log-linear regression
- `n_points::Int`: Number of points used in regression
- `start_time::Float64`: Start time of terminal phase
- `end_time::Float64`: End time of terminal phase
- `points_used::Vector{Int}`: Indices of points used in regression
- `quality_flag::Symbol`: Quality assessment (:good, :warning, :insufficient)
- `warnings::Vector{String}`: Quality warnings
"""
struct LambdaZResult
    lambda_z::Union{Float64,Nothing}
    t_half::Union{Float64,Nothing}
    r_squared::Union{Float64,Nothing}
    adjusted_r_squared::Union{Float64,Nothing}
    intercept::Union{Float64,Nothing}
    n_points::Int
    start_time::Float64
    end_time::Float64
    points_used::Vector{Int}
    quality_flag::Symbol
    warnings::Vector{String}
end

function LambdaZResult(;
    lambda_z::Union{Float64,Nothing} = nothing,
    t_half::Union{Float64,Nothing} = nothing,
    r_squared::Union{Float64,Nothing} = nothing,
    adjusted_r_squared::Union{Float64,Nothing} = nothing,
    intercept::Union{Float64,Nothing} = nothing,
    n_points::Int = 0,
    start_time::Float64 = 0.0,
    end_time::Float64 = 0.0,
    points_used::Vector{Int} = Int[],
    quality_flag::Symbol = :insufficient,
    warnings::Vector{String} = String[]
)
    return LambdaZResult(
        lambda_z, t_half, r_squared, adjusted_r_squared, intercept,
        n_points, start_time, end_time, points_used, quality_flag, warnings
    )
end

# =============================================================================
# NCA Result
# =============================================================================

"""
    NCAResult

Complete NCA analysis results following FDA/EMA guidance.

# Primary Exposure Metrics
- `cmax::Float64`: Maximum observed concentration
- `tmax::Float64`: Time of Cmax
- `cmin::Union{Float64,Nothing}`: Minimum concentration (multiple dose)
- `clast::Float64`: Last measurable concentration
- `tlast::Float64`: Time of last measurable concentration
- `cavg::Union{Float64,Nothing}`: Average concentration over dosing interval

# AUC Metrics
- `auc_0_t::Float64`: AUC from 0 to last measurable concentration
- `auc_0_inf::Union{Float64,Nothing}`: AUC extrapolated to infinity
- `auc_extra_pct::Union{Float64,Nothing}`: Percent of AUC extrapolated
- `auc_0_tau::Union{Float64,Nothing}`: AUC over dosing interval (steady state)
- `aumc_0_t::Float64`: Area under first moment curve (0 to tlast)
- `aumc_0_inf::Union{Float64,Nothing}`: AUMC extrapolated to infinity

# Terminal Phase
- `lambda_z_result::LambdaZResult`: Full lambda_z estimation result

# PK Parameters
- `t_half::Union{Float64,Nothing}`: Terminal half-life
- `mrt::Union{Float64,Nothing}`: Mean residence time
- `cl_f::Union{Float64,Nothing}`: Apparent clearance (extravascular)
- `vz_f::Union{Float64,Nothing}`: Apparent volume of distribution
- `vss::Union{Float64,Nothing}`: Volume of distribution at steady state (IV)

# Multiple Dose Metrics
- `accumulation_index::Union{Float64,Nothing}`: AUC ratio at steady state
- `ptf::Union{Float64,Nothing}`: Peak-trough fluctuation (%)
- `swing::Union{Float64,Nothing}`: Swing (%)

# Dose-Normalized Metrics
- `cmax_dn::Union{Float64,Nothing}`: Dose-normalized Cmax
- `auc_dn::Union{Float64,Nothing}`: Dose-normalized AUC

# Quality
- `quality_flags::Vector{Symbol}`: Quality assessment flags
- `warnings::Vector{String}`: Analysis warnings
- `metadata::Dict{String,Any}`: Additional metadata
"""
struct NCAResult
    # Primary exposure metrics
    cmax::Float64
    tmax::Float64
    cmin::Union{Float64,Nothing}
    clast::Float64
    tlast::Float64
    cavg::Union{Float64,Nothing}

    # AUC metrics
    auc_0_t::Float64
    auc_0_inf::Union{Float64,Nothing}
    auc_extra_pct::Union{Float64,Nothing}
    auc_0_tau::Union{Float64,Nothing}
    aumc_0_t::Float64
    aumc_0_inf::Union{Float64,Nothing}

    # Terminal phase
    lambda_z_result::LambdaZResult

    # PK parameters
    t_half::Union{Float64,Nothing}
    mrt::Union{Float64,Nothing}
    cl_f::Union{Float64,Nothing}
    vz_f::Union{Float64,Nothing}
    vss::Union{Float64,Nothing}

    # Multiple dose metrics
    accumulation_index::Union{Float64,Nothing}
    ptf::Union{Float64,Nothing}
    swing::Union{Float64,Nothing}

    # Dose-normalized metrics
    cmax_dn::Union{Float64,Nothing}
    auc_dn::Union{Float64,Nothing}

    # Quality
    quality_flags::Vector{Symbol}
    warnings::Vector{String}
    metadata::Dict{String,Any}
end

# =============================================================================
# Bioequivalence Results
# =============================================================================

"""
    BioequivalenceResult

Result of bioequivalence analysis.

# Fields
- `parameter::Symbol`: Parameter analyzed (:cmax, :auc, etc.)
- `n_test::Int`: Number of test subjects
- `n_reference::Int`: Number of reference subjects
- `gmr::Float64`: Geometric mean ratio (test/reference)
- `ci_lower::Float64`: Lower bound of 90% CI
- `ci_upper::Float64`: Upper bound of 90% CI
- `cv_intra::Float64`: Intra-subject CV (for crossover)
- `be_conclusion::Symbol`: Conclusion (:bioequivalent, :not_bioequivalent, :inconclusive)
- `be_limits::Tuple{Float64,Float64}`: BE acceptance limits used
"""
struct BioequivalenceResult
    parameter::Symbol
    n_test::Int
    n_reference::Int
    gmr::Float64
    ci_lower::Float64
    ci_upper::Float64
    cv_intra::Float64
    be_conclusion::Symbol
    be_limits::Tuple{Float64,Float64}
end

"""
    TOSTResult

Result of Two One-Sided Tests (TOST) procedure for bioequivalence.

# Fields
- `parameter::Symbol`: Parameter analyzed
- `t_lower::Float64`: T-statistic for lower bound test
- `t_upper::Float64`: T-statistic for upper bound test
- `p_lower::Float64`: P-value for lower bound test
- `p_upper::Float64`: P-value for upper bound test
- `reject_lower::Bool`: Reject H0 for lower bound
- `reject_upper::Bool`: Reject H0 for upper bound
- `be_conclusion::Symbol`: Overall conclusion
"""
struct TOSTResult
    parameter::Symbol
    t_lower::Float64
    t_upper::Float64
    p_lower::Float64
    p_upper::Float64
    reject_lower::Bool
    reject_upper::Bool
    be_conclusion::Symbol
end
