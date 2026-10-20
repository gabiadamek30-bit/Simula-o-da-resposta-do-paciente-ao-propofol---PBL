# Binning Module for VPC
# Implements various binning strategies for time-concentration data

export compute_bins, BinDefinition

"""
Definition of a single bin.

Fields:
- id: Bin identifier (1-indexed)
- lower: Lower time bound (inclusive)
- upper: Upper time bound (exclusive for all except last bin)
- midpoint: Midpoint for plotting
"""
struct BinDefinition
    id::Int
    lower::Float64
    upper::Float64
    midpoint::Float64
end

"""
Compute bin definitions using the specified strategy.

Arguments:
- times: Vector of time points
- strategy: Binning strategy (QuantileBinning, EqualWidthBinning, or KMeansBinning)

Returns:
- Vector{BinDefinition}
"""
function compute_bins(times::Vector{Float64}, strategy::BinningStrategy)::Vector{BinDefinition}
    if isempty(times)
        return BinDefinition[]
    end

    sorted_times = sort(times)

    if strategy isa QuantileBinning
        return _quantile_bins(sorted_times, strategy.n_bins)
    elseif strategy isa EqualWidthBinning
        return _equal_width_bins(sorted_times, strategy.n_bins)
    elseif strategy isa KMeansBinning
        return _kmeans_bins(sorted_times, strategy.n_bins, strategy.max_iter)
    else
        error("Unknown binning strategy: $(typeof(strategy))")
    end
end

"""
Quantile-based binning - each bin contains approximately equal number of observations.
"""
function _quantile_bins(sorted_times::Vector{Float64}, n_bins::Int)::Vector{BinDefinition}
    n = length(sorted_times)
    if n < n_bins
        n_bins = n
    end

    bins = BinDefinition[]

    for i in 1:n_bins
        # Calculate quantile boundaries
        lower_q = (i - 1) / n_bins
        upper_q = i / n_bins

        lower_idx = max(1, Int(floor(lower_q * n)) + 1)
        upper_idx = min(n, Int(ceil(upper_q * n)))

        lower = sorted_times[lower_idx]
        upper = sorted_times[upper_idx]

        # Handle edge cases
        if i == 1
            lower = sorted_times[1] - eps(sorted_times[1])
        end
        if i == n_bins
            upper = sorted_times[end] + eps(sorted_times[end])
        end

        midpoint = (lower + upper) / 2
        push!(bins, BinDefinition(i, lower, upper, midpoint))
    end

    return bins
end

"""
Equal-width binning - each bin covers an equal time range.
"""
function _equal_width_bins(sorted_times::Vector{Float64}, n_bins::Int)::Vector{BinDefinition}
    t_min = sorted_times[1]
    t_max = sorted_times[end]

    # Add small buffer to ensure all points are included
    range = t_max - t_min
    if range == 0.0
        range = 1.0  # Handle constant time
    end
    t_min -= range * 0.001
    t_max += range * 0.001

    bin_width = (t_max - t_min) / n_bins

    bins = BinDefinition[]
    for i in 1:n_bins
        lower = t_min + (i - 1) * bin_width
        upper = t_min + i * bin_width
        midpoint = (lower + upper) / 2
        push!(bins, BinDefinition(i, lower, upper, midpoint))
    end

    return bins
end

"""
K-means based binning using Lloyd's algorithm.
"""
function _kmeans_bins(sorted_times::Vector{Float64}, n_bins::Int, max_iter::Int)::Vector{BinDefinition}
    n = length(sorted_times)
    if n <= n_bins
        return _equal_width_bins(sorted_times, min(n, n_bins))
    end

    # Initialize centroids evenly spaced
    t_min, t_max = sorted_times[1], sorted_times[end]
    centroids = [t_min + (t_max - t_min) * (i - 0.5) / n_bins for i in 1:n_bins]

    assignments = zeros(Int, n)

    for _ in 1:max_iter
        # Assignment step
        new_assignments = zeros(Int, n)
        for (i, t) in enumerate(sorted_times)
            min_dist = Inf
            min_cluster = 1
            for (j, c) in enumerate(centroids)
                dist = abs(t - c)
                if dist < min_dist
                    min_dist = dist
                    min_cluster = j
                end
            end
            new_assignments[i] = min_cluster
        end

        # Check for convergence
        if new_assignments == assignments
            break
        end
        assignments = new_assignments

        # Update step
        for j in 1:n_bins
            cluster_points = [sorted_times[i] for i in 1:n if assignments[i] == j]
            if !isempty(cluster_points)
                centroids[j] = sum(cluster_points) / length(cluster_points)
            end
        end
    end

    # Sort centroids and create bins
    sort!(centroids)

    bins = BinDefinition[]
    for i in 1:n_bins
        if i == 1
            lower = sorted_times[1] - eps(sorted_times[1])
        else
            lower = (centroids[i-1] + centroids[i]) / 2
        end

        if i == n_bins
            upper = sorted_times[end] + eps(sorted_times[end])
        else
            upper = (centroids[i] + centroids[i+1]) / 2
        end

        midpoint = centroids[i]
        push!(bins, BinDefinition(i, lower, upper, midpoint))
    end

    return bins
end

"""
Assign data points to bins.

Arguments:
- times: Vector of time points
- values: Vector of corresponding values
- bins: Vector of BinDefinition

Returns:
- Vector of (bin_id, Vector{Float64}) tuples for each bin
"""
function assign_to_bins(
    times::Vector{Float64},
    values::Vector{Float64},
    bins::Vector{BinDefinition}
)::Vector{Tuple{Int,Vector{Float64}}}
    @assert length(times) == length(values)

    binned_values = [Float64[] for _ in bins]

    for (t, v) in zip(times, values)
        for bin in bins
            if bin.lower <= t < bin.upper || (bin.id == length(bins) && t <= bin.upper)
                push!(binned_values[bin.id], v)
                break
            end
        end
    end

    return [(bin.id, binned_values[bin.id]) for bin in bins]
end

"""
Compute percentile for a vector of values.

Arguments:
- values: Vector of values
- p: Percentile level (0-1)

Returns:
- Percentile value, or NaN if empty
"""
function compute_percentile(values::Vector{Float64}, p::Float64)::Float64
    if isempty(values)
        return NaN
    end
    if length(values) == 1
        return values[1]
    end

    sorted = sort(values)
    n = length(sorted)

    # Linear interpolation method (Type 7 in R)
    h = (n - 1) * p + 1
    lo = Int(floor(h))
    hi = Int(ceil(h))

    if lo == hi
        return sorted[lo]
    end

    return sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo])
end

"""
Bootstrap confidence interval for a percentile.

Arguments:
- values: Vector of values
- p: Percentile level (0-1)
- ci_level: Confidence level (e.g., 0.95)
- n_bootstrap: Number of bootstrap samples
- rng: Random number generator

Returns:
- (lower, median, upper) tuple
"""
function bootstrap_percentile_ci(
    values::Vector{Float64},
    p::Float64,
    ci_level::Float64,
    n_bootstrap::Int,
    rng
)::Tuple{Float64,Float64,Float64}
    if isempty(values) || length(values) < 2
        val = isempty(values) ? NaN : values[1]
        return (val, val, val)
    end

    n = length(values)
    bootstrap_percentiles = Float64[]

    for _ in 1:n_bootstrap
        # Sample with replacement
        sample_indices = rand(rng, 1:n, n)
        sample = [values[i] for i in sample_indices]
        push!(bootstrap_percentiles, compute_percentile(sample, p))
    end

    # Compute CI
    alpha = 1 - ci_level
    lower = compute_percentile(bootstrap_percentiles, alpha / 2)
    median = compute_percentile(bootstrap_percentiles, 0.5)
    upper = compute_percentile(bootstrap_percentiles, 1 - alpha / 2)

    return (lower, median, upper)
end

export assign_to_bins, compute_percentile, bootstrap_percentile_ci
