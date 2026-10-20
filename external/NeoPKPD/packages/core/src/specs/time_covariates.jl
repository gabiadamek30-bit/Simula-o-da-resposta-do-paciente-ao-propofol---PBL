export TimeCovariateKind, StepTimeCovariate, LinearTimeCovariate
export TimeCovariateSeries, TimeVaryingCovariates

abstract type TimeCovariateKind end

struct StepTimeCovariate <: TimeCovariateKind end
struct LinearTimeCovariate <: TimeCovariateKind end

"""
A time series for one covariate.

times:
- sorted, unique time points

values:
- same length as times

kind:
- StepTimeCovariate: value is held constant until next time
- LinearTimeCovariate: linear interpolation between knots
"""
struct TimeCovariateSeries{K<:TimeCovariateKind}
    kind::K
    times::Vector{Float64}
    values::Vector{Float64}
end

"""
Holds time-varying covariates for one individual.
"""
struct TimeVaryingCovariates
    series::Dict{Symbol,Any}  # Symbol => TimeCovariateSeries
end
