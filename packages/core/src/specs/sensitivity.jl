export PerturbationKind, RelativePerturbation, AbsolutePerturbation, Perturbation
export PerturbationPlan, SensitivityMetric

# Global Sensitivity Analysis exports
export GlobalSensitivityMethod, SobolMethod, MorrisMethod
export ParameterBounds, GlobalSensitivitySpec

abstract type PerturbationKind end

"""
Relative perturbation of a parameter: new = base * (1 + delta)
delta is a fraction, e.g. 0.1 means +10 percent
"""
struct RelativePerturbation <: PerturbationKind end

"""
Absolute perturbation of a parameter: new = base + delta
"""
struct AbsolutePerturbation <: PerturbationKind end

struct Perturbation{K<:PerturbationKind}
    kind::K
    param::Symbol
    delta::Float64
end

"""
A named set of perturbations.
"""
struct PerturbationPlan
    name::String
    perturbations::Vector{Perturbation}
end

"""
Simple sensitivity metrics for time series comparison.
"""
struct SensitivityMetric
    max_abs_delta::Float64
    max_rel_delta::Float64
    l2_norm_delta::Float64
end

# ============================================================================
# Global Sensitivity Analysis Types
# ============================================================================

"""
Abstract type for global sensitivity analysis methods.
"""
abstract type GlobalSensitivityMethod end

"""
    SobolMethod(; base_sample_size=1024, compute_second_order=false, bootstrap_samples=1000, bootstrap_ci_level=0.95)

Sobol' sensitivity analysis using Saltelli sampling scheme.

Computes variance-based sensitivity indices:
- First-order indices (Si): Main effect of each parameter
- Total-order indices (STi): Including all interactions

# Arguments
- `base_sample_size::Int`: N in Saltelli scheme. Total evaluations = N*(d+2) where d is number of parameters
- `compute_second_order::Bool`: Whether to compute second-order interaction indices Sij
- `bootstrap_samples::Int`: Number of bootstrap samples for confidence intervals (0 = no CI)
- `bootstrap_ci_level::Float64`: Confidence level for intervals (default 0.95)

# Example
```julia
method = SobolMethod(base_sample_size=1024, bootstrap_samples=1000)
# For 3 parameters: 1024 * (3+2) = 5120 model evaluations
```
"""
struct SobolMethod <: GlobalSensitivityMethod
    base_sample_size::Int
    compute_second_order::Bool
    bootstrap_samples::Int
    bootstrap_ci_level::Float64

    function SobolMethod(;
        base_sample_size::Int=1024,
        compute_second_order::Bool=false,
        bootstrap_samples::Int=1000,
        bootstrap_ci_level::Float64=0.95
    )
        base_sample_size >= 64 || error("base_sample_size must be >= 64, got $base_sample_size")
        bootstrap_samples >= 0 || error("bootstrap_samples must be >= 0")
        0.0 < bootstrap_ci_level < 1.0 || error("bootstrap_ci_level must be in (0, 1), got $bootstrap_ci_level")
        new(base_sample_size, compute_second_order, bootstrap_samples, bootstrap_ci_level)
    end
end

"""
    MorrisMethod(; n_trajectories=10, n_levels=4, delta=nothing)

Morris Elementary Effects (screening) method.

A computationally efficient screening method to identify important parameters
before running more expensive variance-based analysis.

# Arguments
- `n_trajectories::Int`: Number of trajectories r. Total evaluations = r*(d+1)
- `n_levels::Int`: Number of grid levels p (typically 4)
- `delta::Float64`: Step size. Default is p/(2*(p-1)) which ensures even sampling

# Outputs
- μ (mu): Mean elementary effect
- μ* (mu_star): Mean absolute elementary effect (importance measure)
- σ (sigma): Standard deviation (indicates interactions/nonlinearity)

# Example
```julia
method = MorrisMethod(n_trajectories=20)
# For 3 parameters: 20 * (3+1) = 80 model evaluations
```
"""
struct MorrisMethod <: GlobalSensitivityMethod
    n_trajectories::Int
    n_levels::Int
    delta::Float64

    function MorrisMethod(;
        n_trajectories::Int=10,
        n_levels::Int=4,
        delta::Union{Nothing,Float64}=nothing
    )
        n_trajectories >= 4 || error("n_trajectories must be >= 4, got $n_trajectories")
        n_levels >= 2 || error("n_levels must be >= 2, got $n_levels")
        d = delta === nothing ? n_levels / (2 * (n_levels - 1)) : delta
        d > 0.0 || error("delta must be positive, got $d")
        new(n_trajectories, n_levels, d)
    end
end

"""
    ParameterBounds(params, lower, upper)
    ParameterBounds(bounds::Dict{Symbol,Tuple{Float64,Float64}})

Specification of parameter bounds for global sensitivity analysis.

Each parameter has a lower and upper bound defining its variation range.
Parameters are varied uniformly within these bounds.

# Examples
```julia
# Vector constructor
bounds = ParameterBounds([:CL, :V], [0.5, 5.0], [2.0, 20.0])

# Dict constructor
bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
```
"""
struct ParameterBounds
    params::Vector{Symbol}
    lower::Vector{Float64}
    upper::Vector{Float64}

    function ParameterBounds(params::Vector{Symbol}, lower::Vector{Float64}, upper::Vector{Float64})
        length(params) == length(lower) == length(upper) ||
            error("params, lower, upper must have same length")
        all(lower .< upper) ||
            error("All lower bounds must be less than upper bounds")
        new(params, lower, upper)
    end
end

# Convenience constructor from Dict
function ParameterBounds(bounds::Dict{Symbol,Tuple{Float64,Float64}})
    params = collect(keys(bounds))
    lower = [bounds[p][1] for p in params]
    upper = [bounds[p][2] for p in params]
    ParameterBounds(params, lower, upper)
end

# Get number of parameters
Base.length(b::ParameterBounds) = length(b.params)

"""
    GlobalSensitivitySpec(method, bounds; observation=:conc, seed=12345)

Complete specification for global sensitivity analysis.

# Arguments
- `method::GlobalSensitivityMethod`: Analysis method (SobolMethod or MorrisMethod)
- `bounds::ParameterBounds`: Parameter bounds for variation
- `observation::Symbol`: Which model output to analyze (default :conc)
- `seed::UInt64`: Random seed for reproducibility

# Example
```julia
bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
spec = GlobalSensitivitySpec(SobolMethod(), bounds; observation=:conc)
```
"""
struct GlobalSensitivitySpec{M<:GlobalSensitivityMethod}
    method::M
    bounds::ParameterBounds
    observation::Symbol
    seed::UInt64

    function GlobalSensitivitySpec(
        method::M,
        bounds::ParameterBounds;
        observation::Symbol=:conc,
        seed::UInt64=UInt64(12345)
    ) where {M<:GlobalSensitivityMethod}
        length(bounds) >= 1 || error("At least one parameter required")
        new{M}(method, bounds, observation, seed)
    end
end
