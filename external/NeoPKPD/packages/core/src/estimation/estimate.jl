# Main Estimation Entry Point
# Dispatches to appropriate method implementation

using Optim
using LineSearches
using StableRNGs

export estimate

"""
    estimate(observed, model_spec, config; grid, solver, parallel_config) -> EstimationResult

Estimate population pharmacokinetic/pharmacodynamic parameters using NLME methods.

This is the main entry point for parameter estimation. The function dispatches
to the appropriate algorithm based on the method specified in the config.

# Arguments
- `observed::ObservedData`: Observed data with subject concentrations/effects
- `model_spec::ModelSpec`: Model specification with structural model
- `config::EstimationConfig`: Configuration including method, initial values, bounds

# Keyword Arguments
- `grid::SimGrid`: Simulation time grid
- `solver::SolverSpec`: ODE solver specification
- `parallel_config::ParallelConfig`: Parallel execution configuration (optional)
  - Use `ParallelConfig(ThreadedBackend())` for multi-threaded execution
  - Use `ParallelConfig(SerialBackend())` for serial execution (default)
  - Set `seed` in parallel_config for reproducible parallel execution

# Returns
- `EstimationResult`: Estimated parameters, standard errors, and diagnostics

# Examples

```julia
# Set up estimation config
config = EstimationConfig(
    FOCEIMethod();
    theta_init = [10.0, 50.0],  # CL, V
    theta_names = [:CL, :V],
    omega_init = diagm([0.09, 0.04]),  # 30% CV for CL, 20% CV for V
    omega_names = [:eta_CL, :eta_V],
    sigma_init = ResidualErrorSpec(ProportionalError(), ProportionalErrorParams(0.1), :conc, UInt64(1))
)

# Run estimation (serial)
result = estimate(observed, model_spec, config; grid=grid, solver=solver)

# Run estimation (parallel with all available threads)
parallel = ParallelConfig(ThreadedBackend())
result = estimate(observed, model_spec, config; grid=grid, solver=solver, parallel_config=parallel)

# Check results
println("OFV: ", result.ofv)
println("CL: ", result.theta[1], " ± ", result.theta_se[1])
```
"""
function estimate(
    observed::ObservedData,
    model_spec::ModelSpec,
    config::EstimationConfig;
    grid::SimGrid,
    solver::SolverSpec,
    parallel_config::ParallelConfig=ParallelConfig(SerialBackend())
)::EstimationResult
    start_time = time()
    rng = StableRNG(config.seed)

    if config.verbose
        println("Starting estimation with $(typeof(config.method))")
        println("  Subjects: $(n_subjects(observed))")
        println("  Observations: $(n_observations(observed))")
        println("  Parameters: theta=$(length(config.theta_init)), omega=$(size(config.omega_init, 1))")
        if is_parallel(parallel_config)
            println("  Parallel: $(n_workers(parallel_config)) workers ($(typeof(parallel_config.backend)))")
        else
            println("  Parallel: disabled (serial execution)")
        end
    end

    # Dispatch to appropriate method with parallel config
    result = _estimate_impl(observed, model_spec, config, grid, solver, rng, parallel_config)

    elapsed = time() - start_time
    if config.verbose
        println("Estimation completed in $(round(elapsed, digits=1)) seconds")
        println("  Convergence: $(result.convergence)")
        println("  OFV: $(round(result.ofv, digits=3))")
    end

    return result
end

"""
Internal dispatch to method-specific implementation.
"""
function _estimate_impl(
    observed::ObservedData,
    model_spec::ModelSpec,
    config::EstimationConfig{FOCEIMethod},
    grid::SimGrid,
    solver::SolverSpec,
    rng,
    parallel_config::ParallelConfig
)::EstimationResult
    return foce_estimate(observed, model_spec, config, grid, solver, rng; parallel_config=parallel_config)
end

function _estimate_impl(
    observed::ObservedData,
    model_spec::ModelSpec,
    config::EstimationConfig{SAEMMethod},
    grid::SimGrid,
    solver::SolverSpec,
    rng,
    parallel_config::ParallelConfig
)::EstimationResult
    return saem_estimate(observed, model_spec, config, grid, solver, rng; parallel_config=parallel_config)
end

function _estimate_impl(
    observed::ObservedData,
    model_spec::ModelSpec,
    config::EstimationConfig{LaplacianMethod},
    grid::SimGrid,
    solver::SolverSpec,
    rng,
    parallel_config::ParallelConfig
)::EstimationResult
    return laplacian_estimate(observed, model_spec, config, grid, solver, rng; parallel_config=parallel_config)
end

# ------------------------------------------------------------------
# Shared Utility Functions
# ------------------------------------------------------------------

"""
Extract subject data for estimation.
Returns vector of (subject_id, times, observations, doses) tuples.
"""
function extract_subject_data(observed::ObservedData)
    return [(s.subject_id, s.times, s.observations, s.doses) for s in observed.subjects]
end

"""
Extract subject data with covariates for estimation with covariate effects.
Returns vector of (subject_id, times, observations, doses, covariates) tuples.
"""
function extract_subject_data_with_covariates(observed::ObservedData)
    return [(s.subject_id, s.times, s.observations, s.doses, s.covariates) for s in observed.subjects]
end

"""
Extract occasion indices from subject covariates.
Looks for :OCC or :OCCASION key in covariates.
Returns vector of occasion indices (1-indexed), or vector of 1s if no occasion data.
"""
function extract_occasion_indices(covariates::Dict{Symbol,Any}, n_obs::Int)::Vector{Int}
    # Try common occasion variable names
    for key in [:OCC, :OCCASION, :occ, :occasion]
        if haskey(covariates, key)
            occ_data = covariates[key]
            if occ_data isa Vector
                return Int.(occ_data)
            elseif occ_data isa Number
                return fill(Int(occ_data), n_obs)
            end
        end
    end
    # Default: all observations on occasion 1
    return ones(Int, n_obs)
end

"""
Create individual predictions given theta, eta, and model.
"""
function compute_individual_predictions(
    theta::Vector{Float64},
    eta::Vector{Float64},
    subject_times::Vector{Float64},
    subject_doses::Vector{DoseEvent},
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec
)::Vector{Float64}
    # Apply eta to get individual parameters
    individual_params = apply_eta_to_params(theta, eta, model_spec)

    # Create individual model spec
    ind_model_spec = ModelSpec(
        model_spec.kind,
        model_spec.name,
        individual_params,
        subject_doses
    )

    # Simulate
    result = simulate(ind_model_spec, grid, solver)

    # Interpolate to observation times
    ipred = Vector{Float64}(undef, length(subject_times))
    for (i, t) in enumerate(subject_times)
        idx = searchsortedfirst(result.t, t)
        if idx > length(result.t)
            idx = length(result.t)
        elseif idx > 1 && abs(result.t[idx] - t) > abs(result.t[idx-1] - t)
            idx = idx - 1
        end
        ipred[i] = result.observations[:conc][idx]
    end

    return ipred
end

"""
Apply eta (random effects) to population parameters to get individual parameters.
"""
function apply_eta_to_params(
    theta::Vector{Float64},
    eta::Vector{Float64},
    model_spec::ModelSpec
)
    # This function needs to map theta + eta to individual parameters
    # The mapping depends on the model type and parameterization

    # For now, assume log-normal IIV: individual_param = theta * exp(eta)
    # This is the standard NONMEM approach

    # The number of etas should match the number of IIV parameters
    n_eta = length(eta)

    # Create individual parameters by applying exponential transformation
    # theta[1:n_eta] are the population means that get individual variability
    individual_theta = copy(theta)
    for i in 1:min(n_eta, length(theta))
        individual_theta[i] = theta[i] * exp(eta[i])
    end

    # Convert theta vector to model params type
    return theta_to_params(individual_theta, model_spec)
end

"""
Convert theta vector to appropriate model parameters struct.
"""
function theta_to_params(theta::Vector{Float64}, model_spec::ModelSpec{OneCompIVBolus})
    return OneCompIVBolusParams(theta[1], theta[2])  # CL, V
end

function theta_to_params(theta::Vector{Float64}, model_spec::ModelSpec{OneCompOralFirstOrder})
    return OneCompOralFirstOrderParams(theta[1], theta[2], theta[3])  # Ka, CL, V
end

function theta_to_params(theta::Vector{Float64}, model_spec::ModelSpec{TwoCompIVBolus})
    return TwoCompIVBolusParams(theta[1], theta[2], theta[3], theta[4])  # CL, V1, Q, V2
end

function theta_to_params(theta::Vector{Float64}, model_spec::ModelSpec{TwoCompOral})
    return TwoCompOralParams(theta[1], theta[2], theta[3], theta[4], theta[5])  # Ka, CL, V1, Q, V2
end

function theta_to_params(theta::Vector{Float64}, model_spec::ModelSpec{ThreeCompIVBolus})
    return ThreeCompIVBolusParams(theta[1], theta[2], theta[3], theta[4], theta[5], theta[6])
end

# Generic fallback - uses existing params type
function theta_to_params(theta::Vector{Float64}, model_spec::ModelSpec)
    # For other models, try to construct from theta
    ParamType = typeof(model_spec.params)
    return ParamType(theta...)
end

export extract_subject_data, extract_subject_data_with_covariates, extract_occasion_indices
export compute_individual_predictions, apply_eta_to_params, theta_to_params

"""
Compute -2 log-likelihood contribution for a single observation.
"""
function observation_log_likelihood(
    y::Float64,           # Observed value
    f::Float64,           # Predicted value (IPRED)
    spec::ResidualErrorSpec
)::Float64
    # Compute residual variance
    var_res = residual_variance(f, spec)

    # Handle special cases
    if var_res <= 0.0 || !isfinite(var_res)
        return Inf
    end

    # Gaussian log-likelihood: -0.5 * (log(2π) + log(var) + (y-f)²/var)
    residual = y - f
    ll = -0.5 * (log(2π) + log(var_res) + residual^2 / var_res)

    return -2.0 * ll  # Return -2LL contribution
end

export observation_log_likelihood

"""
Compute AIC from OFV and number of parameters.
"""
function compute_aic(ofv::Float64, n_params::Int)::Float64
    return ofv + 2.0 * n_params
end

"""
Compute BIC from OFV, number of parameters, and number of observations.
"""
function compute_bic(ofv::Float64, n_params::Int, n_obs::Int)::Float64
    return ofv + log(n_obs) * n_params
end

export compute_aic, compute_bic
