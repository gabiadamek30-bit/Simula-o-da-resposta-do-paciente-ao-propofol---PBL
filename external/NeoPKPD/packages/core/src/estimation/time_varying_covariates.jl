# Time-Varying Covariates in Estimation
# Professional-grade implementation following FDA/EMA guidance
#
# Supports:
# - Time-varying effects on fixed effects (theta)
# - Multiple interpolation strategies (step, linear, LOCF, NOCB, spline)
# - Covariate centering and scaling
# - Integration with FOCE-I and SAEM estimation
# - Automatic covariate extraction from observation times
#
# References:
# - FDA Guidance: Population PK Analysis (2022)
# - Mould & Upton (2013): Basic Concepts in Population Modeling
# - Ette & Williams (2007): Pharmacometrics - The Science of Quantitative Pharmacology

using LinearAlgebra

export TimeVaryingCovariateEffect, CovariateInterpolation
export StepInterpolation, LinearInterpolation, LOCFInterpolation, NOCBInterpolation, SplineInterpolation
export CovariateTransform, NoTransform, LogTransform, PowerTransform, NormalizedTransform
export TimeVaryingCovariateSpec, TimeVaryingCovariateConfig
export get_covariate_at_time, get_covariates_at_times
export apply_time_varying_covariates, extract_time_varying_covariates
export CovariateTimeSeries, SubjectTimeVaryingCovariates
export compute_predictions_with_tv_covariates, apply_time_varying_covariates_vector

# =============================================================================
# Interpolation Strategy Types
# =============================================================================

"""
Abstract type for covariate interpolation strategies.
"""
abstract type CovariateInterpolation end

"""
Step interpolation (Last Observation Carried Forward within interval).
Value at time t is the most recent observation before or at t.
Most common for lab values, weight, etc.
"""
struct StepInterpolation <: CovariateInterpolation end

"""
Linear interpolation between observations.
Good for continuous physiological measurements.
"""
struct LinearInterpolation <: CovariateInterpolation end

"""
Last Observation Carried Forward (explicit LOCF).
Like step, but emphasizes the clinical convention.
"""
struct LOCFInterpolation <: CovariateInterpolation end

"""
Next Observation Carried Backward.
Value at time t is the next observation at or after t.
"""
struct NOCBInterpolation <: CovariateInterpolation end

"""
Cubic spline interpolation.
Smooth interpolation for dense covariate measurements.
"""
struct SplineInterpolation <: CovariateInterpolation end

# =============================================================================
# Covariate Transform Types
# =============================================================================

"""
Abstract type for covariate transformations.
"""
abstract type CovariateTransform end

"""
No transformation - use raw covariate value.
"""
struct NoTransform <: CovariateTransform end

"""
Log transformation: effect = θ * log(COV/REF).
Common for weight, creatinine clearance.
"""
struct LogTransform <: CovariateTransform
    reference::Float64

    LogTransform(reference::Float64=1.0) = new(reference)
end

"""
Power transformation: effect = θ * (COV/REF)^power.
Common for allometric scaling (power = 0.75 for CL, 1.0 for V).
"""
struct PowerTransform <: CovariateTransform
    reference::Float64
    power::Float64

    PowerTransform(reference::Float64, power::Float64=1.0) = new(reference, power)
end

"""
Normalized transformation: effect = θ * (COV - REF)/SCALE.
Centers and scales the covariate.
"""
struct NormalizedTransform <: CovariateTransform
    reference::Float64
    scale::Float64

    NormalizedTransform(reference::Float64, scale::Float64=1.0) = new(reference, scale)
end

# =============================================================================
# Covariate Time Series
# =============================================================================

"""
Time series data for a single covariate.

# Fields
- `name`: Covariate name (e.g., :WT, :CRCL, :ALB)
- `times`: Measurement times (sorted)
- `values`: Covariate values at each time
- `interpolation`: Interpolation strategy
- `units`: Units string (optional, for documentation)

# Example
```julia
# Weight measured at baseline and day 14
wt_series = CovariateTimeSeries(
    :WT,
    [0.0, 336.0],  # hours
    [70.0, 68.5],  # kg
    LinearInterpolation()
)
```
"""
struct CovariateTimeSeries
    name::Symbol
    times::Vector{Float64}
    values::Vector{Float64}
    interpolation::CovariateInterpolation
    units::String

    function CovariateTimeSeries(
        name::Symbol,
        times::Vector{Float64},
        values::Vector{Float64},
        interpolation::CovariateInterpolation=StepInterpolation();
        units::String=""
    )
        @assert length(times) == length(values) "times and values must have same length"
        @assert length(times) >= 1 "Need at least one covariate measurement"
        @assert issorted(times) "times must be sorted"
        new(name, times, values, interpolation, units)
    end
end

"""
Collection of time-varying covariates for a single subject.

# Fields
- `subject_id`: Subject identifier
- `covariates`: Dict mapping covariate names to time series
- `baseline_covariates`: Static/baseline covariates (single values)
"""
struct SubjectTimeVaryingCovariates
    subject_id::String
    covariates::Dict{Symbol, CovariateTimeSeries}
    baseline_covariates::Dict{Symbol, Float64}

    function SubjectTimeVaryingCovariates(
        subject_id::String;
        covariates::Dict{Symbol, CovariateTimeSeries}=Dict{Symbol, CovariateTimeSeries}(),
        baseline_covariates::Dict{Symbol, Float64}=Dict{Symbol, Float64}()
    )
        new(subject_id, covariates, baseline_covariates)
    end
end

# =============================================================================
# Time-Varying Covariate Effect Specification
# =============================================================================

"""
Specification for a time-varying covariate effect on a parameter.

# Fields
- `parameter`: Name of the parameter affected (e.g., :CL, :V)
- `parameter_index`: Index in theta vector (1-based)
- `covariate`: Name of the covariate (e.g., :WT, :CRCL)
- `transform`: Transformation to apply
- `interpolation`: Default interpolation (can be overridden by series)
- `theta_index`: Index in theta vector for effect magnitude (optional)
- `fixed_effect`: Fixed effect magnitude if not estimated (optional)

# Example
```julia
# Weight effect on clearance (allometric scaling)
effect = TimeVaryingCovariateEffect(
    parameter = :CL,
    parameter_index = 1,
    covariate = :WT,
    transform = PowerTransform(70.0, 0.75)  # (WT/70)^0.75
)

# Time-varying creatinine clearance effect
effect = TimeVaryingCovariateEffect(
    parameter = :CL,
    parameter_index = 1,
    covariate = :CRCL,
    transform = NormalizedTransform(100.0, 50.0),  # (CRCL - 100)/50
    theta_index = 3  # Effect magnitude estimated as theta[3]
)
```
"""
struct TimeVaryingCovariateEffect
    parameter::Symbol
    parameter_index::Int
    covariate::Symbol
    transform::CovariateTransform
    interpolation::CovariateInterpolation
    theta_index::Union{Nothing, Int}
    fixed_effect::Union{Nothing, Float64}

    function TimeVaryingCovariateEffect(;
        parameter::Symbol,
        parameter_index::Int,
        covariate::Symbol,
        transform::CovariateTransform=NoTransform(),
        interpolation::CovariateInterpolation=StepInterpolation(),
        theta_index::Union{Nothing, Int}=nothing,
        fixed_effect::Union{Nothing, Float64}=nothing
    )
        @assert parameter_index >= 1 "parameter_index must be positive"
        if theta_index !== nothing
            @assert theta_index >= 1 "theta_index must be positive"
        end
        if theta_index === nothing && fixed_effect === nothing
            # Default: multiplicative effect with magnitude 1
            fixed_effect = 1.0
        end
        new(parameter, parameter_index, covariate, transform, interpolation,
            theta_index, fixed_effect)
    end
end

"""
Configuration for time-varying covariates in estimation.

# Fields
- `effects`: Vector of TimeVaryingCovariateEffect specifications
- `center_at_baseline`: Use baseline value as reference (default: true)
- `handle_missing`: How to handle missing covariate values (:error, :locf, :interpolate)
- `extrapolate`: Allow extrapolation beyond measurement times (default: true)
- `verbose`: Print covariate processing info (default: false)
"""
struct TimeVaryingCovariateConfig
    effects::Vector{TimeVaryingCovariateEffect}
    center_at_baseline::Bool
    handle_missing::Symbol
    extrapolate::Bool
    verbose::Bool

    function TimeVaryingCovariateConfig(;
        effects::Vector{TimeVaryingCovariateEffect}=TimeVaryingCovariateEffect[],
        center_at_baseline::Bool=true,
        handle_missing::Symbol=:locf,
        extrapolate::Bool=true,
        verbose::Bool=false
    )
        @assert handle_missing in [:error, :locf, :interpolate] "Invalid handle_missing option"
        new(effects, center_at_baseline, handle_missing, extrapolate, verbose)
    end
end

# =============================================================================
# Interpolation Functions
# =============================================================================

"""
    get_covariate_at_time(series::CovariateTimeSeries, t::Float64) -> Float64

Get interpolated covariate value at time t.
"""
function get_covariate_at_time(series::CovariateTimeSeries, t::Float64)::Float64
    return interpolate_covariate(series.times, series.values, t, series.interpolation)
end

"""
    interpolate_covariate(times, values, t, interp) -> Float64

Core interpolation function.
"""
function interpolate_covariate(
    times::Vector{Float64},
    values::Vector{Float64},
    t::Float64,
    interp::StepInterpolation
)::Float64
    # Step/LOCF: use value from most recent measurement
    if t <= times[1]
        return values[1]
    end
    if t >= times[end]
        return values[end]
    end

    # Find rightmost time <= t
    idx = searchsortedlast(times, t)
    return values[idx]
end

function interpolate_covariate(
    times::Vector{Float64},
    values::Vector{Float64},
    t::Float64,
    interp::LOCFInterpolation
)::Float64
    # Same as step interpolation
    return interpolate_covariate(times, values, t, StepInterpolation())
end

function interpolate_covariate(
    times::Vector{Float64},
    values::Vector{Float64},
    t::Float64,
    interp::LinearInterpolation
)::Float64
    if t <= times[1]
        return values[1]
    end
    if t >= times[end]
        return values[end]
    end

    # Find interval containing t
    hi = searchsortedfirst(times, t)
    lo = hi - 1

    # Linear interpolation
    t0, t1 = times[lo], times[hi]
    v0, v1 = values[lo], values[hi]
    w = (t - t0) / (t1 - t0)

    return (1.0 - w) * v0 + w * v1
end

function interpolate_covariate(
    times::Vector{Float64},
    values::Vector{Float64},
    t::Float64,
    interp::NOCBInterpolation
)::Float64
    # Next Observation Carried Backward
    if t <= times[1]
        return values[1]
    end
    if t >= times[end]
        return values[end]
    end

    # Find leftmost time >= t
    idx = searchsortedfirst(times, t)
    return values[idx]
end

function interpolate_covariate(
    times::Vector{Float64},
    values::Vector{Float64},
    t::Float64,
    interp::SplineInterpolation
)::Float64
    if length(times) < 3
        # Fall back to linear for insufficient points
        return interpolate_covariate(times, values, t, LinearInterpolation())
    end

    # Clamp to range
    t_clamped = clamp(t, times[1], times[end])

    # Create cubic spline interpolation
    try
        itp = cubic_spline_interpolation(times, values)
        return itp(t_clamped)
    catch
        # Fall back to linear on failure
        return interpolate_covariate(times, values, t, LinearInterpolation())
    end
end

"""
    CubicSplineInterpolator

Natural cubic spline interpolator struct.
Stores spline coefficients for efficient repeated evaluation.
"""
struct CubicSplineInterpolator
    x::Vector{Float64}
    y::Vector{Float64}
    a::Vector{Float64}  # y values (coefficients)
    b::Vector{Float64}  # first derivative coefficients
    c::Vector{Float64}  # second derivative coefficients / 2
    d::Vector{Float64}  # third derivative coefficients / 6
end

function (itp::CubicSplineInterpolator)(t::Float64)
    return evaluate_cubic_spline(itp, t)
end

"""
    cubic_spline_interpolation(x::Vector{Float64}, y::Vector{Float64}) -> CubicSplineInterpolator

Create a natural cubic spline interpolator.
Natural boundary conditions: second derivative = 0 at endpoints.

Uses tridiagonal system solver for O(n) construction.
"""
function cubic_spline_interpolation(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)

    if n < 2
        error("Need at least 2 points for spline interpolation")
    end

    if n == 2
        # Degenerate case - linear interpolation
        slope = (y[2] - y[1]) / (x[2] - x[1])
        return CubicSplineInterpolator(x, y, copy(y), [slope, slope], [0.0, 0.0], [0.0, 0.0])
    end

    # Compute interval widths
    h = diff(x)

    # Set up tridiagonal system for second derivatives (c)
    # Natural boundary: c[1] = c[n] = 0

    # Interior equations from continuity of second derivative
    # h[i-1] * c[i-1] + 2*(h[i-1] + h[i]) * c[i] + h[i] * c[i+1] = 3 * (dy[i]/h[i] - dy[i-1]/h[i-1])

    # For natural spline, c[1] = c[n] = 0, so we solve for c[2:n-1]

    m = n - 2  # Number of interior points

    if m == 0
        # Only 2 points - linear interpolation already handled above
        error("Unexpected state in cubic spline")
    end

    # Build tridiagonal system
    # Lower diagonal, main diagonal, upper diagonal
    lower = zeros(m - 1)
    main = zeros(m)
    upper = zeros(m - 1)
    rhs = zeros(m)

    for i in 1:m
        idx = i + 1  # actual index in x, y (2 to n-1)

        # Main diagonal coefficient
        main[i] = 2.0 * (h[idx - 1] + h[idx])

        # RHS
        dy_curr = (y[idx + 1] - y[idx]) / h[idx]
        dy_prev = (y[idx] - y[idx - 1]) / h[idx - 1]
        rhs[i] = 3.0 * (dy_curr - dy_prev)

        # Off-diagonals
        if i > 1
            lower[i - 1] = h[idx - 1]
        end
        if i < m
            upper[i] = h[idx]
        end
    end

    # Solve tridiagonal system using Thomas algorithm
    c_interior = solve_tridiagonal(lower, main, upper, rhs)

    # Build full c vector with boundary conditions
    c = zeros(n)
    c[2:n-1] = c_interior
    # c[1] = c[n] = 0 (natural boundary conditions)

    # Compute b and d coefficients
    a = copy(y)
    b = zeros(n)
    d = zeros(n)

    for i in 1:n-1
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i])
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0
    end

    # Last segment coefficients (for extrapolation)
    b[n] = b[n - 1] + 2.0 * c[n - 1] * h[n - 1] + 3.0 * d[n - 1] * h[n - 1]^2
    d[n] = d[n - 1]

    return CubicSplineInterpolator(x, y, a, b, c, d)
end

"""
    solve_tridiagonal(lower, main, upper, rhs) -> Vector{Float64}

Solve tridiagonal system using Thomas algorithm (O(n)).
"""
function solve_tridiagonal(
    lower::Vector{Float64},
    main::Vector{Float64},
    upper::Vector{Float64},
    rhs::Vector{Float64}
)
    n = length(main)

    if n == 0
        return Float64[]
    end

    if n == 1
        return [rhs[1] / main[1]]
    end

    # Forward elimination
    c_prime = zeros(n)
    d_prime = zeros(n)

    c_prime[1] = upper[1] / main[1]
    d_prime[1] = rhs[1] / main[1]

    for i in 2:n
        denom = main[i] - lower[i - 1] * c_prime[i - 1]
        if i < n
            c_prime[i] = upper[i] / denom
        end
        d_prime[i] = (rhs[i] - lower[i - 1] * d_prime[i - 1]) / denom
    end

    # Back substitution
    x = zeros(n)
    x[n] = d_prime[n]

    for i in n-1:-1:1
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    end

    return x
end

"""
    evaluate_cubic_spline(itp::CubicSplineInterpolator, t::Float64) -> Float64

Evaluate cubic spline at point t.
"""
function evaluate_cubic_spline(itp::CubicSplineInterpolator, t::Float64)
    x = itp.x
    n = length(x)

    # Handle extrapolation
    if t <= x[1]
        # Extrapolate using first segment
        dx = t - x[1]
        return itp.a[1] + itp.b[1] * dx + itp.c[1] * dx^2 + itp.d[1] * dx^3
    end

    if t >= x[n]
        # Extrapolate using last segment
        dx = t - x[n - 1]
        return itp.a[n - 1] + itp.b[n - 1] * dx + itp.c[n - 1] * dx^2 + itp.d[n - 1] * dx^3
    end

    # Find interval containing t
    i = searchsortedlast(x, t)
    if i == n
        i = n - 1
    end

    # Evaluate polynomial on interval [x[i], x[i+1]]
    dx = t - x[i]
    return itp.a[i] + itp.b[i] * dx + itp.c[i] * dx^2 + itp.d[i] * dx^3
end

"""
    get_covariates_at_times(subject_cov::SubjectTimeVaryingCovariates, cov_name::Symbol, times::Vector{Float64}) -> Vector{Float64}

Get covariate values at multiple times.
"""
function get_covariates_at_times(
    subject_cov::SubjectTimeVaryingCovariates,
    cov_name::Symbol,
    times::Vector{Float64}
)::Vector{Float64}
    if haskey(subject_cov.covariates, cov_name)
        series = subject_cov.covariates[cov_name]
        return [get_covariate_at_time(series, t) for t in times]
    elseif haskey(subject_cov.baseline_covariates, cov_name)
        # Static covariate - same value at all times
        return fill(subject_cov.baseline_covariates[cov_name], length(times))
    else
        error("Covariate $cov_name not found for subject $(subject_cov.subject_id)")
    end
end

# =============================================================================
# Transform Application
# =============================================================================

"""
    apply_transform(value, transform) -> Float64

Apply covariate transformation.
"""
function apply_transform(value::Float64, transform::NoTransform)::Float64
    return value
end

function apply_transform(value::Float64, transform::LogTransform)::Float64
    return log(value / transform.reference)
end

function apply_transform(value::Float64, transform::PowerTransform)::Float64
    return (value / transform.reference) ^ transform.power
end

function apply_transform(value::Float64, transform::NormalizedTransform)::Float64
    return (value - transform.reference) / transform.scale
end

# =============================================================================
# Time-Varying Covariate Application
# =============================================================================

"""
    apply_time_varying_covariates(theta, t, effects, subject_cov) -> Vector{Float64}

Apply time-varying covariate effects to theta at time t.

Returns modified theta vector with covariate effects applied.

# Arguments
- `theta`: Population parameters
- `t`: Time point
- `effects`: Vector of TimeVaryingCovariateEffect
- `subject_cov`: Subject's time-varying covariates

# Returns
Modified theta vector with covariate effects applied.
"""
function apply_time_varying_covariates(
    theta::Vector{Float64},
    t::Float64,
    effects::Vector{TimeVaryingCovariateEffect},
    subject_cov::SubjectTimeVaryingCovariates
)::Vector{Float64}
    if isempty(effects)
        return theta
    end

    theta_modified = copy(theta)

    for effect in effects
        # Get covariate value at time t
        cov_value = get_covariate_value(subject_cov, effect.covariate, t)

        # Apply transformation
        transformed_value = apply_transform(cov_value, effect.transform)

        # Get effect magnitude
        effect_magnitude = if effect.theta_index !== nothing
            theta[effect.theta_index]
        else
            effect.fixed_effect
        end

        # Apply effect to parameter
        # For power transform: theta_i * transformed_value
        # For other transforms: theta_i * (1 + effect_magnitude * transformed_value)
        if effect.transform isa PowerTransform
            theta_modified[effect.parameter_index] *= transformed_value
        else
            theta_modified[effect.parameter_index] *= (1.0 + effect_magnitude * transformed_value)
        end
    end

    return theta_modified
end

"""
Helper to get covariate value from subject covariates.
"""
function get_covariate_value(
    subject_cov::SubjectTimeVaryingCovariates,
    cov_name::Symbol,
    t::Float64
)::Float64
    if haskey(subject_cov.covariates, cov_name)
        return get_covariate_at_time(subject_cov.covariates[cov_name], t)
    elseif haskey(subject_cov.baseline_covariates, cov_name)
        return subject_cov.baseline_covariates[cov_name]
    else
        error("Covariate $cov_name not found for subject $(subject_cov.subject_id)")
    end
end

"""
    apply_time_varying_covariates_vector(theta, times, effects, subject_cov) -> Vector{Vector{Float64}}

Apply time-varying covariate effects at multiple time points.

Returns vector of modified theta vectors, one for each time point.
"""
function apply_time_varying_covariates_vector(
    theta::Vector{Float64},
    times::Vector{Float64},
    effects::Vector{TimeVaryingCovariateEffect},
    subject_cov::SubjectTimeVaryingCovariates
)::Vector{Vector{Float64}}
    return [apply_time_varying_covariates(theta, t, effects, subject_cov) for t in times]
end

# =============================================================================
# Extraction from SubjectData
# =============================================================================

"""
    extract_time_varying_covariates(subject::SubjectData) -> SubjectTimeVaryingCovariates

Extract time-varying covariates from SubjectData structure.

Handles both scalar (baseline) and time-varying covariates.
"""
function extract_time_varying_covariates(subject::SubjectData)::SubjectTimeVaryingCovariates
    tv_covariates = Dict{Symbol, CovariateTimeSeries}()
    baseline_covariates = Dict{Symbol, Float64}()

    for (name, value) in subject.covariates
        if value isa Number
            # Scalar covariate (baseline)
            baseline_covariates[name] = Float64(value)
        elseif value isa Vector
            # Time-varying covariate stored as vector at observation times
            times = subject.times
            values = Float64.(value)
            if length(values) == length(times)
                tv_covariates[name] = CovariateTimeSeries(
                    name, times, values, StepInterpolation()
                )
            elseif length(values) == 1
                # Single value - treat as baseline
                baseline_covariates[name] = values[1]
            end
        elseif value isa Tuple && length(value) == 2
            # (times, values) tuple
            times, values = value
            tv_covariates[name] = CovariateTimeSeries(
                name, Float64.(times), Float64.(values), StepInterpolation()
            )
        elseif value isa CovariateTimeSeries
            tv_covariates[name] = value
        elseif value isa TimeCovariateSeries
            # Convert from existing TimeCovariateSeries type
            interp = value.kind isa StepTimeCovariate ? StepInterpolation() : LinearInterpolation()
            tv_covariates[name] = CovariateTimeSeries(
                name, value.times, value.values, interp
            )
        end
    end

    return SubjectTimeVaryingCovariates(
        subject.subject_id;
        covariates=tv_covariates,
        baseline_covariates=baseline_covariates
    )
end

"""
    extract_all_time_varying_covariates(observed::ObservedData) -> Vector{SubjectTimeVaryingCovariates}

Extract time-varying covariates for all subjects.
"""
function extract_all_time_varying_covariates(
    observed::ObservedData
)::Vector{SubjectTimeVaryingCovariates}
    return [extract_time_varying_covariates(subj) for subj in observed.subjects]
end

# =============================================================================
# Prediction with Time-Varying Covariates
# =============================================================================

"""
    compute_predictions_with_tv_covariates(theta, eta, times, doses, model_spec, effects, subject_cov) -> Vector{Float64}

Compute predictions with time-varying covariate effects.

This is the main function for integrating time-varying covariates into model predictions.
"""
function compute_predictions_with_tv_covariates(
    theta::Vector{Float64},
    eta::Vector{Float64},
    times::Vector{Float64},
    doses::Vector{DoseEvent},
    model_spec::ModelSpec,
    effects::Vector{TimeVaryingCovariateEffect},
    subject_cov::SubjectTimeVaryingCovariates
)::Vector{Float64}
    n_obs = length(times)
    pred = Vector{Float64}(undef, n_obs)

    # Get total dose
    total_dose = sum(d.amount for d in doses)

    for (i, t) in enumerate(times)
        # Apply time-varying covariate effects to theta at this time
        theta_t = apply_time_varying_covariates(theta, t, effects, subject_cov)

        # Compute prediction with modified theta
        pred[i] = compute_single_prediction(theta_t, eta, t, total_dose, model_spec)
    end

    return pred
end

"""
Compute single prediction at time t with given parameters.
"""
function compute_single_prediction(
    theta::Vector{Float64},
    eta::Vector{Float64},
    t::Float64,
    dose::Float64,
    model_spec::ModelSpec
)::Float64
    if model_spec.kind isa OneCompIVBolus
        CL = theta[1] * exp(eta[1])
        V = theta[2] * exp(eta[2])
        return one_comp_iv_bolus_analytic(t, dose, CL, V)

    elseif model_spec.kind isa OneCompOralFirstOrder
        CL = theta[1] * exp(eta[1])
        V = theta[2] * exp(eta[2])
        ka = length(theta) > 2 ? theta[3] * (length(eta) > 2 ? exp(eta[3]) : 1.0) : 1.0
        return one_comp_oral_analytic(t, dose, CL, V, ka)

    elseif model_spec.kind isa TwoCompIVBolus
        CL = theta[1] * exp(eta[1])
        V1 = theta[2] * exp(eta[2])
        Q = length(theta) > 2 ? theta[3] : theta[1] * 0.5
        V2 = length(theta) > 3 ? theta[4] : theta[2] * 2.0
        return two_comp_iv_bolus_analytic(t, dose, CL, V1, Q, V2)

    else
        # Fallback - return first theta as concentration proxy
        return theta[1] * exp(-theta[1] / theta[2] * t) * dose / theta[2]
    end
end

# Analytic solutions (duplicated for standalone use)
function one_comp_iv_bolus_analytic(t::Float64, dose::Float64, CL::Float64, V::Float64)::Float64
    k = CL / V
    return (dose / V) * exp(-k * t)
end

function one_comp_oral_analytic(t::Float64, dose::Float64, CL::Float64, V::Float64, ka::Float64)::Float64
    ke = CL / V
    if abs(ka - ke) < 1e-10
        # Handle ka ≈ ke case
        return (dose * ka / V) * t * exp(-ke * t)
    end
    return (dose * ka / (V * (ka - ke))) * (exp(-ke * t) - exp(-ka * t))
end

function two_comp_iv_bolus_analytic(t::Float64, dose::Float64, CL::Float64, V1::Float64, Q::Float64, V2::Float64)::Float64
    k10 = CL / V1
    k12 = Q / V1
    k21 = Q / V2

    # Eigenvalues
    sum_k = k10 + k12 + k21
    prod_k = k10 * k21
    discriminant = sum_k^2 - 4 * prod_k

    if discriminant < 0
        discriminant = 0.0
    end

    sqrt_disc = sqrt(discriminant)
    alpha = (sum_k + sqrt_disc) / 2
    beta = (sum_k - sqrt_disc) / 2

    # Macro-constants
    A = dose / V1 * (alpha - k21) / (alpha - beta)
    B = dose / V1 * (k21 - beta) / (alpha - beta)

    return A * exp(-alpha * t) + B * exp(-beta * t)
end

# =============================================================================
# Convenience Constructors
# =============================================================================

"""
    weight_effect_on_cl(; reference=70.0, power=0.75, theta_index=nothing) -> TimeVaryingCovariateEffect

Create allometric weight effect on clearance.
"""
function weight_effect_on_cl(;
    reference::Float64=70.0,
    power::Float64=0.75,
    theta_index::Union{Nothing, Int}=nothing
)::TimeVaryingCovariateEffect
    return TimeVaryingCovariateEffect(
        parameter=:CL,
        parameter_index=1,
        covariate=:WT,
        transform=PowerTransform(reference, power),
        interpolation=StepInterpolation(),
        theta_index=theta_index,
        fixed_effect=1.0
    )
end

"""
    weight_effect_on_v(; reference=70.0, power=1.0, theta_index=nothing) -> TimeVaryingCovariateEffect

Create allometric weight effect on volume.
"""
function weight_effect_on_v(;
    reference::Float64=70.0,
    power::Float64=1.0,
    theta_index::Union{Nothing, Int}=nothing
)::TimeVaryingCovariateEffect
    return TimeVaryingCovariateEffect(
        parameter=:V,
        parameter_index=2,
        covariate=:WT,
        transform=PowerTransform(reference, power),
        interpolation=StepInterpolation(),
        theta_index=theta_index,
        fixed_effect=1.0
    )
end

"""
    crcl_effect_on_cl(; reference=100.0, scale=50.0, theta_index=3) -> TimeVaryingCovariateEffect

Create creatinine clearance effect on drug clearance.
"""
function crcl_effect_on_cl(;
    reference::Float64=100.0,
    scale::Float64=50.0,
    theta_index::Int=3
)::TimeVaryingCovariateEffect
    return TimeVaryingCovariateEffect(
        parameter=:CL,
        parameter_index=1,
        covariate=:CRCL,
        transform=NormalizedTransform(reference, scale),
        interpolation=LinearInterpolation(),
        theta_index=theta_index,
        fixed_effect=nothing
    )
end

"""
    albumin_effect_on_cl(; reference=4.0, theta_index=nothing) -> TimeVaryingCovariateEffect

Create albumin effect on clearance (for highly protein-bound drugs).
"""
function albumin_effect_on_cl(;
    reference::Float64=4.0,
    theta_index::Union{Nothing, Int}=nothing
)::TimeVaryingCovariateEffect
    return TimeVaryingCovariateEffect(
        parameter=:CL,
        parameter_index=1,
        covariate=:ALB,
        transform=PowerTransform(reference, 1.0),
        interpolation=LinearInterpolation(),
        theta_index=theta_index,
        fixed_effect=1.0
    )
end

export weight_effect_on_cl, weight_effect_on_v, crcl_effect_on_cl, albumin_effect_on_cl

# =============================================================================
# Validation and Diagnostics
# =============================================================================

"""
    validate_covariate_coverage(subject_cov::SubjectTimeVaryingCovariates, obs_times::Vector{Float64}, effects::Vector{TimeVaryingCovariateEffect}) -> Vector{String}

Check that all required covariates are available for all observation times.
Returns vector of warning messages.
"""
function validate_covariate_coverage(
    subject_cov::SubjectTimeVaryingCovariates,
    obs_times::Vector{Float64},
    effects::Vector{TimeVaryingCovariateEffect}
)::Vector{String}
    warnings = String[]

    for effect in effects
        cov_name = effect.covariate

        if haskey(subject_cov.covariates, cov_name)
            series = subject_cov.covariates[cov_name]

            # Check time coverage
            min_cov_time = minimum(series.times)
            max_cov_time = maximum(series.times)
            min_obs_time = minimum(obs_times)
            max_obs_time = maximum(obs_times)

            if min_obs_time < min_cov_time
                push!(warnings,
                    "Subject $(subject_cov.subject_id): Covariate $cov_name not available before t=$(min_cov_time), " *
                    "but observations start at t=$(min_obs_time). Will extrapolate."
                )
            end

            if max_obs_time > max_cov_time
                push!(warnings,
                    "Subject $(subject_cov.subject_id): Covariate $cov_name not available after t=$(max_cov_time), " *
                    "but observations extend to t=$(max_obs_time). Will extrapolate."
                )
            end

        elseif !haskey(subject_cov.baseline_covariates, cov_name)
            push!(warnings,
                "Subject $(subject_cov.subject_id): Covariate $cov_name not found."
            )
        end
    end

    return warnings
end

"""
    summarize_time_varying_covariates(subject_cov::SubjectTimeVaryingCovariates) -> String

Generate summary of time-varying covariates for a subject.
"""
function summarize_time_varying_covariates(
    subject_cov::SubjectTimeVaryingCovariates
)::String
    lines = String["Time-Varying Covariates for $(subject_cov.subject_id):"]

    if !isempty(subject_cov.baseline_covariates)
        push!(lines, "  Baseline covariates:")
        for (name, value) in subject_cov.baseline_covariates
            push!(lines, "    $name = $value")
        end
    end

    if !isempty(subject_cov.covariates)
        push!(lines, "  Time-varying covariates:")
        for (name, series) in subject_cov.covariates
            n_points = length(series.times)
            t_range = "($(minimum(series.times)) - $(maximum(series.times)))"
            v_range = "[$(minimum(series.values)) - $(maximum(series.values))]"
            interp = typeof(series.interpolation)
            push!(lines, "    $name: $n_points points, t=$t_range, v=$v_range, interp=$interp")
        end
    end

    return join(lines, "\n")
end

export validate_covariate_coverage, summarize_time_varying_covariates

# =============================================================================
# Pretty Printing
# =============================================================================

function Base.show(io::IO, effect::TimeVaryingCovariateEffect)
    println(io, "TimeVaryingCovariateEffect")
    println(io, "  Parameter: $(effect.parameter) (index $(effect.parameter_index))")
    println(io, "  Covariate: $(effect.covariate)")
    println(io, "  Transform: $(typeof(effect.transform))")
    if effect.theta_index !== nothing
        println(io, "  Effect θ index: $(effect.theta_index)")
    else
        println(io, "  Fixed effect: $(effect.fixed_effect)")
    end
end

function Base.show(io::IO, series::CovariateTimeSeries)
    n = length(series.times)
    println(io, "CovariateTimeSeries(:$(series.name))")
    println(io, "  $n measurements from t=$(series.times[1]) to t=$(series.times[end])")
    println(io, "  Values: $(round.(extrema(series.values), digits=2))")
    println(io, "  Interpolation: $(typeof(series.interpolation))")
end

function Base.show(io::IO, config::TimeVaryingCovariateConfig)
    println(io, "TimeVaryingCovariateConfig")
    println(io, "  Effects: $(length(config.effects))")
    for effect in config.effects
        println(io, "    - $(effect.covariate) → $(effect.parameter)")
    end
    println(io, "  Center at baseline: $(config.center_at_baseline)")
    println(io, "  Handle missing: $(config.handle_missing)")
end
