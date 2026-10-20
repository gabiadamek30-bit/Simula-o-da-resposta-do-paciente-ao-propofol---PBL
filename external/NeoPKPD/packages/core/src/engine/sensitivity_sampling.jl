# Sampling algorithms for Global Sensitivity Analysis
# Implements Saltelli sampling for Sobol' and trajectory sampling for Morris

using Random
using StableRNGs

export generate_saltelli_samples, generate_saltelli_samples_second_order
export generate_morris_trajectories
export transform_samples, generate_sobol_sequence

"""
    generate_sobol_sequence(n::Int, d::Int, rng::AbstractRNG)

Generate quasi-random samples using Sobol' low-discrepancy sequence approximation.

For true Sobol' sequences, consider using Sobol.jl package.
This implementation uses scrambled uniform sampling as a reasonable approximation.

# Arguments
- `n::Int`: Number of samples
- `d::Int`: Number of dimensions (parameters)
- `rng::AbstractRNG`: Random number generator

# Returns
- `Matrix{Float64}`: n×d matrix with values in [0, 1]
"""
function generate_sobol_sequence(n::Int, d::Int, rng::AbstractRNG)::Matrix{Float64}
    # Use Latin Hypercube Sampling for better space-filling than pure random
    # This is a reasonable approximation when Sobol.jl is not available
    samples = Matrix{Float64}(undef, n, d)

    for j in 1:d
        # Create stratified samples
        perm = randperm(rng, n)
        for i in 1:n
            # Stratified uniform sampling within each stratum
            samples[i, j] = (perm[i] - rand(rng)) / n
        end
    end

    return samples
end

"""
    generate_saltelli_samples(bounds::ParameterBounds, N::Int, rng::AbstractRNG)

Generate sample matrices for Sobol' sensitivity analysis using Saltelli scheme.

The Saltelli scheme generates samples efficiently for computing first-order (Si)
and total-order (STi) Sobol' indices.

# Arguments
- `bounds::ParameterBounds`: Parameter bounds specification
- `N::Int`: Base sample size
- `rng::AbstractRNG`: Random number generator for reproducibility

# Returns
Named tuple with:
- `A::Matrix{Float64}`: N×d matrix A
- `B::Matrix{Float64}`: N×d matrix B
- `AB::Vector{Matrix{Float64}}`: d matrices where AB[i] = A with column i from B
- `BA::Vector{Matrix{Float64}}`: d matrices where BA[i] = B with column i from A (for STi)

Total model evaluations required: N*(d+2) for Si and STi

# References
- Saltelli et al. (2010) "Variance based sensitivity analysis of model output"
- Sobol' (2001) "Global sensitivity indices for nonlinear mathematical models"
"""
function generate_saltelli_samples(
    bounds::ParameterBounds,
    N::Int,
    rng::AbstractRNG
)::NamedTuple{(:A, :B, :AB, :BA), Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Matrix{Float64}}, Vector{Matrix{Float64}}}}
    d = length(bounds)

    # Generate two independent sample matrices in [0,1]^d
    A_unit = generate_sobol_sequence(N, d, rng)
    B_unit = generate_sobol_sequence(N, d, rng)

    # Transform to parameter space
    A = transform_samples(A_unit, bounds)
    B = transform_samples(B_unit, bounds)

    # Generate AB matrices: A with column i replaced by B's column i
    # Used for first-order indices
    AB = Vector{Matrix{Float64}}(undef, d)
    for i in 1:d
        AB[i] = copy(A)
        AB[i][:, i] = B[:, i]
    end

    # Generate BA matrices: B with column i replaced by A's column i
    # Used for total-order indices (Jansen estimator)
    BA = Vector{Matrix{Float64}}(undef, d)
    for i in 1:d
        BA[i] = copy(B)
        BA[i][:, i] = A[:, i]
    end

    return (A=A, B=B, AB=AB, BA=BA)
end

"""
    generate_saltelli_samples_second_order(bounds::ParameterBounds, N::Int, rng::AbstractRNG)

Generate sample matrices for second-order Sobol' indices using extended Saltelli scheme.

In addition to the standard matrices (A, B, AB), this generates ABij matrices where
columns i and j are both taken from B. These are required for computing second-order
interaction indices Sij.

# Arguments
- `bounds::ParameterBounds`: Parameter bounds specification
- `N::Int`: Base sample size
- `rng::AbstractRNG`: Random number generator for reproducibility

# Returns
Named tuple with:
- `A::Matrix{Float64}`: N×d matrix A
- `B::Matrix{Float64}`: N×d matrix B
- `AB::Vector{Matrix{Float64}}`: d matrices where AB[i] = A with column i from B
- `BA::Vector{Matrix{Float64}}`: d matrices where BA[i] = B with column i from A
- `ABij::Dict{Tuple{Int,Int},Matrix{Float64}}`: (d*(d-1)/2) matrices where ABij[(i,j)] = A with columns i,j from B

Total model evaluations required: N*(d+2) + N*(d*(d-1)/2) for Si, STi, and Sij

# References
- Saltelli (2002) "Making best use of model evaluations to compute sensitivity indices"
- Homma & Saltelli (1996) "Importance measures in global sensitivity analysis"
"""
function generate_saltelli_samples_second_order(
    bounds::ParameterBounds,
    N::Int,
    rng::AbstractRNG
)::NamedTuple{(:A, :B, :AB, :BA, :ABij), Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Matrix{Float64}}, Vector{Matrix{Float64}}, Dict{Tuple{Int,Int},Matrix{Float64}}}}
    d = length(bounds)

    # Generate two independent sample matrices in [0,1]^d
    A_unit = generate_sobol_sequence(N, d, rng)
    B_unit = generate_sobol_sequence(N, d, rng)

    # Transform to parameter space
    A = transform_samples(A_unit, bounds)
    B = transform_samples(B_unit, bounds)

    # Generate AB matrices: A with column i replaced by B's column i
    AB = Vector{Matrix{Float64}}(undef, d)
    for i in 1:d
        AB[i] = copy(A)
        AB[i][:, i] = B[:, i]
    end

    # Generate BA matrices: B with column i replaced by A's column i
    BA = Vector{Matrix{Float64}}(undef, d)
    for i in 1:d
        BA[i] = copy(B)
        BA[i][:, i] = A[:, i]
    end

    # Generate ABij matrices for second-order indices
    # ABij[(i,j)] = A with columns i AND j replaced by B's columns i and j
    ABij = Dict{Tuple{Int,Int},Matrix{Float64}}()
    for i in 1:(d-1)
        for j in (i+1):d
            ABij_mat = copy(A)
            ABij_mat[:, i] = B[:, i]
            ABij_mat[:, j] = B[:, j]
            ABij[(i, j)] = ABij_mat
        end
    end

    return (A=A, B=B, AB=AB, BA=BA, ABij=ABij)
end

"""
    transform_samples(samples::Matrix{Float64}, bounds::ParameterBounds)

Transform samples from [0,1]^d unit hypercube to parameter space.

# Arguments
- `samples::Matrix{Float64}`: N×d matrix with values in [0, 1]
- `bounds::ParameterBounds`: Parameter bounds specification

# Returns
- `Matrix{Float64}`: N×d matrix with values in [lower, upper] for each parameter
"""
function transform_samples(samples::Matrix{Float64}, bounds::ParameterBounds)::Matrix{Float64}
    N, d = size(samples)
    @assert d == length(bounds) "Sample dimensions must match number of parameters"

    transformed = Matrix{Float64}(undef, N, d)
    for j in 1:d
        range = bounds.upper[j] - bounds.lower[j]
        for i in 1:N
            transformed[i, j] = bounds.lower[j] + samples[i, j] * range
        end
    end

    return transformed
end

"""
    generate_morris_trajectories(bounds::ParameterBounds, r::Int, p::Int, delta::Float64, rng::AbstractRNG)

Generate optimized Morris trajectories for elementary effects computation.

Each trajectory consists of d+1 points where d is the number of parameters.
Starting from a random base point, each step changes exactly one parameter by ±delta.

# Arguments
- `bounds::ParameterBounds`: Parameter bounds specification
- `r::Int`: Number of trajectories
- `p::Int`: Number of grid levels
- `delta::Float64`: Step size in normalized [0,1] space
- `rng::AbstractRNG`: Random number generator

# Returns
- `Vector{Matrix{Float64}}`: r trajectories, each (d+1)×d matrix
  - Row 1 is the base point
  - Rows 2 to d+1 are points after perturbing each parameter once

# Algorithm
For each trajectory:
1. Generate random base point on the p-level grid
2. Generate random permutation of parameters
3. For each parameter in permutation order:
   - Move by +delta or -delta (randomly chosen, respecting bounds)
   - Record the new point

Total evaluations: r × (d + 1)

# References
- Morris (1991) "Factorial Sampling Plans for Preliminary Computational Experiments"
- Campolongo et al. (2007) "An effective screening design for sensitivity analysis"
"""
function generate_morris_trajectories(
    bounds::ParameterBounds,
    r::Int,
    p::Int,
    delta::Float64,
    rng::AbstractRNG
)::Vector{Matrix{Float64}}
    d = length(bounds)
    trajectories = Vector{Matrix{Float64}}(undef, r)

    # Grid values in [0,1] space
    grid_values = collect(0:1/(p-1):1)

    for traj_idx in 1:r
        # Initialize trajectory matrix: (d+1) points × d parameters
        trajectory = Matrix{Float64}(undef, d + 1, d)

        # Generate random base point on the grid (in [0,1] space)
        # Ensure we can move in either direction
        base_unit = Vector{Float64}(undef, d)
        for j in 1:d
            # Choose from grid values that allow movement by delta
            valid_levels = filter(v -> (v + delta <= 1.0) || (v - delta >= 0.0), grid_values)
            if isempty(valid_levels)
                valid_levels = grid_values
            end
            base_unit[j] = rand(rng, valid_levels)
        end

        # Random permutation of parameters (order of perturbation)
        perm = randperm(rng, d)

        # Transform base point to parameter space
        current_unit = copy(base_unit)
        trajectory[1, :] = transform_point(current_unit, bounds)

        # Generate trajectory by perturbing one parameter at a time
        for step in 1:d
            param_idx = perm[step]

            # Determine direction of perturbation
            can_increase = current_unit[param_idx] + delta <= 1.0
            can_decrease = current_unit[param_idx] - delta >= 0.0

            if can_increase && can_decrease
                # Choose randomly
                direction = rand(rng, Bool) ? 1.0 : -1.0
            elseif can_increase
                direction = 1.0
            elseif can_decrease
                direction = -1.0
            else
                # Fallback: use smaller step
                direction = rand(rng, Bool) ? 1.0 : -1.0
            end

            # Apply perturbation
            current_unit[param_idx] = clamp(current_unit[param_idx] + direction * delta, 0.0, 1.0)

            # Transform to parameter space and store
            trajectory[step + 1, :] = transform_point(current_unit, bounds)
        end

        trajectories[traj_idx] = trajectory
    end

    return trajectories
end

"""
Transform a single point from [0,1]^d to parameter space.
"""
function transform_point(point_unit::Vector{Float64}, bounds::ParameterBounds)::Vector{Float64}
    d = length(bounds)
    point = Vector{Float64}(undef, d)
    for j in 1:d
        range = bounds.upper[j] - bounds.lower[j]
        point[j] = bounds.lower[j] + point_unit[j] * range
    end
    return point
end

"""
    get_trajectory_parameter_order(trajectory::Matrix{Float64}, bounds::ParameterBounds, delta::Float64)

Determine which parameter changed at each step of a Morris trajectory.

# Returns
- `Vector{Int}`: Parameter indices in order of perturbation
"""
function get_trajectory_parameter_order(
    trajectory::Matrix{Float64},
    bounds::ParameterBounds,
    delta::Float64
)::Vector{Int}
    d = size(trajectory, 2)
    n_steps = size(trajectory, 1) - 1

    order = Vector{Int}(undef, n_steps)

    for step in 1:n_steps
        prev_point = trajectory[step, :]
        curr_point = trajectory[step + 1, :]

        # Find which parameter changed
        for j in 1:d
            diff = abs(curr_point[j] - prev_point[j])
            range = bounds.upper[j] - bounds.lower[j]
            normalized_diff = diff / range

            # Check if this parameter changed (approximately by delta)
            if normalized_diff > delta * 0.5  # Allow some tolerance
                order[step] = j
                break
            end
        end
    end

    return order
end
