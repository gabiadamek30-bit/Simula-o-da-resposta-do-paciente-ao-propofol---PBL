# VPC Types - Visual Predictive Check
# Types and structures for VPC analysis

export VPCConfig, VPCBin, VPCResult, VPCPercentileData
export BinningStrategy, QuantileBinning, EqualWidthBinning, KMeansBinning

# ============================================================================
# Binning Strategies
# ============================================================================

"""
Abstract type for binning strategies used in VPC.
"""
abstract type BinningStrategy end

"""
Quantile-based binning - bins have equal number of observations.

Fields:
- n_bins: Number of bins to create
"""
struct QuantileBinning <: BinningStrategy
    n_bins::Int

    function QuantileBinning(n_bins::Int=10)
        @assert n_bins >= 2 "Must have at least 2 bins"
        new(n_bins)
    end
end

"""
Equal-width binning - bins have equal time ranges.

Fields:
- n_bins: Number of bins to create
"""
struct EqualWidthBinning <: BinningStrategy
    n_bins::Int

    function EqualWidthBinning(n_bins::Int=10)
        @assert n_bins >= 2 "Must have at least 2 bins"
        new(n_bins)
    end
end

"""
K-means based binning - bins are formed by clustering time points.

Fields:
- n_bins: Number of bins/clusters to create
- max_iter: Maximum iterations for k-means
"""
struct KMeansBinning <: BinningStrategy
    n_bins::Int
    max_iter::Int

    function KMeansBinning(n_bins::Int=10; max_iter::Int=100)
        @assert n_bins >= 2 "Must have at least 2 bins"
        new(n_bins, max_iter)
    end
end

# ============================================================================
# VPC Configuration
# ============================================================================

"""
Configuration for Visual Predictive Check analysis.

Fields:
- pi_levels: Prediction interval levels (e.g., [0.05, 0.50, 0.95] for 5th, 50th, 95th percentiles)
- ci_level: Confidence interval level for percentile uncertainty (e.g., 0.95)
- binning: Binning strategy for grouping time points
- prediction_corrected: If true, compute prediction-corrected VPC (pcVPC)
- stratify_by: Covariate names to stratify by
- lloq: Lower limit of quantitation (for handling BLQ)
- n_simulations: Number of simulations to run
- n_bootstrap: Number of bootstrap samples for CI calculation
- seed: Random seed for reproducibility
"""
struct VPCConfig
    pi_levels::Vector{Float64}
    ci_level::Float64
    binning::BinningStrategy
    prediction_corrected::Bool
    stratify_by::Vector{Symbol}
    lloq::Union{Nothing,Float64}
    n_simulations::Int
    n_bootstrap::Int
    seed::UInt64

    function VPCConfig(;
        pi_levels::Vector{Float64}=[0.05, 0.50, 0.95],
        ci_level::Float64=0.95,
        binning::BinningStrategy=QuantileBinning(10),
        prediction_corrected::Bool=false,
        stratify_by::Vector{Symbol}=Symbol[],
        lloq::Union{Nothing,Float64}=nothing,
        n_simulations::Int=200,
        n_bootstrap::Int=500,
        seed::UInt64=UInt64(12345)
    )
        @assert all(0.0 < p < 1.0 for p in pi_levels) "PI levels must be between 0 and 1"
        @assert 0.0 < ci_level < 1.0 "CI level must be between 0 and 1"
        @assert n_simulations >= 10 "Need at least 10 simulations"
        @assert n_bootstrap >= 100 "Need at least 100 bootstrap samples"
        new(sort(pi_levels), ci_level, binning, prediction_corrected, stratify_by,
            lloq, n_simulations, n_bootstrap, seed)
    end
end

# ============================================================================
# VPC Bin Data
# ============================================================================

"""
Percentile data for a single bin in VPC.

Fields:
- percentile: The percentile level (e.g., 0.05, 0.50, 0.95)
- observed: Observed percentile value
- simulated_median: Median of simulated percentiles
- simulated_lower: Lower CI bound of simulated percentiles
- simulated_upper: Upper CI bound of simulated percentiles
"""
struct VPCPercentileData
    percentile::Float64
    observed::Float64
    simulated_median::Float64
    simulated_lower::Float64
    simulated_upper::Float64
end

"""
Data for a single time bin in VPC analysis.

Fields:
- bin_id: Bin identifier (1-indexed)
- time_min: Minimum time in bin
- time_max: Maximum time in bin
- time_midpoint: Midpoint time for plotting
- n_observed: Number of observed data points in bin
- n_simulated: Number of simulated data points per simulation
- percentiles: Vector of VPCPercentileData for each PI level
"""
struct VPCBin
    bin_id::Int
    time_min::Float64
    time_max::Float64
    time_midpoint::Float64
    n_observed::Int
    n_simulated::Int
    percentiles::Vector{VPCPercentileData}
end

# ============================================================================
# VPC Result
# ============================================================================

"""
Result of Visual Predictive Check analysis.

Fields:
- config: The VPCConfig used
- bins: Vector of VPCBin with computed statistics
- n_subjects_observed: Number of subjects in observed data
- n_observations_observed: Total observations in observed data
- n_simulations: Number of simulations performed
- strata: Strata label (if stratified)
- simulation_seed: Seed used for simulations
"""
struct VPCResult
    config::VPCConfig
    bins::Vector{VPCBin}
    n_subjects_observed::Int
    n_observations_observed::Int
    n_simulations::Int
    strata::String
    simulation_seed::UInt64
end

"""
Get the time midpoints for all bins (useful for plotting).
"""
function bin_midpoints(result::VPCResult)::Vector{Float64}
    return [b.time_midpoint for b in result.bins]
end

"""
Get observed percentile values for a specific PI level.
"""
function observed_percentile(result::VPCResult, level::Float64)::Vector{Float64}
    values = Float64[]
    for bin in result.bins
        for p in bin.percentiles
            if isapprox(p.percentile, level; atol=1e-6)
                push!(values, p.observed)
                break
            end
        end
    end
    return values
end

"""
Get simulated median percentile values for a specific PI level.
"""
function simulated_median(result::VPCResult, level::Float64)::Vector{Float64}
    values = Float64[]
    for bin in result.bins
        for p in bin.percentiles
            if isapprox(p.percentile, level; atol=1e-6)
                push!(values, p.simulated_median)
                break
            end
        end
    end
    return values
end

"""
Get simulated CI lower bounds for a specific PI level.
"""
function simulated_lower(result::VPCResult, level::Float64)::Vector{Float64}
    values = Float64[]
    for bin in result.bins
        for p in bin.percentiles
            if isapprox(p.percentile, level; atol=1e-6)
                push!(values, p.simulated_lower)
                break
            end
        end
    end
    return values
end

"""
Get simulated CI upper bounds for a specific PI level.
"""
function simulated_upper(result::VPCResult, level::Float64)::Vector{Float64}
    values = Float64[]
    for bin in result.bins
        for p in bin.percentiles
            if isapprox(p.percentile, level; atol=1e-6)
                push!(values, p.simulated_upper)
                break
            end
        end
    end
    return values
end

export bin_midpoints, observed_percentile, simulated_median, simulated_lower, simulated_upper
