# Real Subject Exposure Simulation
# Connects trial simulation to the actual PK engine

using StableRNGs

export simulate_subject_exposure, SubjectExposure, apply_individual_variability

"""
    SubjectExposure

Results of exposure simulation for a single subject.

# Fields
- `subject_id`: Subject identifier
- `times`: Observation time points
- `concentrations`: Plasma concentrations at each time point
- `pd_response`: PD response values (if applicable)
- `pk_metrics`: Derived PK metrics (Cmax, Tmax, AUC)
"""
struct SubjectExposure
    subject_id::Union{Int, String}
    times::Vector{Float64}
    concentrations::Vector{Float64}
    pd_response::Union{Nothing, Vector{Float64}}
    pk_metrics::Dict{Symbol, Float64}
end

"""
    apply_individual_variability(params, eta::Vector{Float64}, model_kind::ModelKind)

Apply inter-individual variability (IIV) to population parameters.

# Arguments
- `params`: Population-level PK parameters
- `eta`: Vector of random effects for this subject
- `model_kind`: The model type

# Returns
- Individual parameters with IIV applied (exponential scaling)
"""
function apply_individual_variability(params::OneCompIVBolusParams, eta::Vector{Float64}, ::OneCompIVBolus)
    # Exponential IIV: param_i = param_pop * exp(eta)
    CL_i = params.CL * exp(length(eta) >= 1 ? eta[1] : 0.0)
    V_i = params.V * exp(length(eta) >= 2 ? eta[2] : 0.0)
    return OneCompIVBolusParams(CL_i, V_i)
end

function apply_individual_variability(params::OneCompOralFirstOrderParams, eta::Vector{Float64}, ::OneCompOralFirstOrder)
    Ka_i = params.Ka * exp(length(eta) >= 1 ? eta[1] : 0.0)
    CL_i = params.CL * exp(length(eta) >= 2 ? eta[2] : 0.0)
    V_i = params.V * exp(length(eta) >= 3 ? eta[3] : 0.0)
    return OneCompOralFirstOrderParams(Ka_i, CL_i, V_i)
end

function apply_individual_variability(params::TwoCompIVBolusParams, eta::Vector{Float64}, ::TwoCompIVBolus)
    CL_i = params.CL * exp(length(eta) >= 1 ? eta[1] : 0.0)
    V1_i = params.V1 * exp(length(eta) >= 2 ? eta[2] : 0.0)
    Q_i = params.Q * exp(length(eta) >= 3 ? eta[3] : 0.0)
    V2_i = params.V2 * exp(length(eta) >= 4 ? eta[4] : 0.0)
    return TwoCompIVBolusParams(CL_i, V1_i, Q_i, V2_i)
end

function apply_individual_variability(params::TwoCompOralParams, eta::Vector{Float64}, ::TwoCompOral)
    Ka_i = params.Ka * exp(length(eta) >= 1 ? eta[1] : 0.0)
    CL_i = params.CL * exp(length(eta) >= 2 ? eta[2] : 0.0)
    V1_i = params.V1 * exp(length(eta) >= 3 ? eta[3] : 0.0)
    Q_i = params.Q * exp(length(eta) >= 4 ? eta[4] : 0.0)
    V2_i = params.V2 * exp(length(eta) >= 5 ? eta[5] : 0.0)
    return TwoCompOralParams(Ka_i, CL_i, V1_i, Q_i, V2_i)
end

# Generic fallback for other model types
function apply_individual_variability(params, eta::Vector{Float64}, model_kind::ModelKind)
    # Return params unchanged if no specific method exists
    return params
end

"""
    calculate_pk_metrics(times::Vector{Float64}, concentrations::Vector{Float64})

Calculate standard PK metrics from time-concentration data.

# Returns
- Dict with :cmax, :tmax, :auc_0_t, :auc_0_inf, :t_half
"""
function calculate_pk_metrics(times::Vector{Float64}, concentrations::Vector{Float64})
    metrics = Dict{Symbol, Float64}()

    if isempty(concentrations) || all(isnan, concentrations)
        return metrics
    end

    # Filter valid values
    valid_idx = findall(c -> !isnan(c) && c >= 0, concentrations)
    if isempty(valid_idx)
        return metrics
    end

    valid_conc = concentrations[valid_idx]
    valid_times = times[valid_idx]

    # Cmax and Tmax
    cmax_idx = argmax(valid_conc)
    metrics[:cmax] = valid_conc[cmax_idx]
    metrics[:tmax] = valid_times[cmax_idx]

    # AUC by trapezoidal rule
    auc_0_t = 0.0
    for i in 1:(length(valid_times)-1)
        dt = valid_times[i+1] - valid_times[i]
        auc_0_t += dt * (valid_conc[i] + valid_conc[i+1]) / 2
    end
    metrics[:auc_0_t] = auc_0_t

    # Estimate terminal phase and extrapolate AUC
    if length(valid_times) >= 3 && valid_conc[end] > 0 && valid_conc[end-1] > valid_conc[end]
        lambda_z = log(valid_conc[end-1] / valid_conc[end]) / (valid_times[end] - valid_times[end-1])
        if lambda_z > 0
            auc_extra = valid_conc[end] / lambda_z
            metrics[:auc_0_inf] = auc_0_t + auc_extra
            metrics[:t_half] = log(2) / lambda_z
            metrics[:lambda_z] = lambda_z
        else
            metrics[:auc_0_inf] = auc_0_t
        end
    else
        metrics[:auc_0_inf] = auc_0_t
    end

    return metrics
end

"""
    simulate_subject_exposure(model_spec, subject, dose_events, observation_times;
                               grid=nothing, solver=nothing, include_iiv=true)

Simulate PK/PD exposure for a single subject using the actual model engine.

# Arguments
- `model_spec`: ModelSpec containing model type and population parameters
- `subject`: VirtualSubject with covariates and random effects (eta)
- `dose_events`: Vector of DoseEvent for this subject
- `observation_times`: Times at which to observe concentrations

# Keyword Arguments
- `grid`: Optional SimGrid (created automatically if not provided)
- `solver`: Optional SolverSpec (defaults to Tsit5)
- `include_iiv`: Whether to apply inter-individual variability (default true)

# Returns
- `SubjectExposure` with times, concentrations, and PK metrics

# Example
```julia
model_spec = ModelSpec(OneCompIVBolus(), "pk", OneCompIVBolusParams(10.0, 50.0), DoseEvent[])
subject = VirtualSubject(1, 45.0, 70.0, :male, :caucasian)  # Uses convenience constructor
doses = [DoseEvent(0.0, 100.0)]
exposure = simulate_subject_exposure(model_spec, subject, doses, [0.0, 1.0, 2.0, 4.0, 8.0])
```
"""
function simulate_subject_exposure(
    model_spec::ModelSpec,
    subject::VirtualSubject,
    dose_events::Vector{DoseEvent},
    observation_times::Vector{Float64};
    grid::Union{Nothing, SimGrid} = nothing,
    solver::Union{Nothing, SolverSpec} = nothing,
    include_iiv::Bool = true
)::SubjectExposure

    # Get individual parameters (apply IIV if applicable)
    individual_params = if include_iiv && hasfield(typeof(subject), :other) &&
                          haskey(subject.other, :eta) && !isempty(subject.other[:eta])
        eta = subject.other[:eta]
        apply_individual_variability(model_spec.params, eta, model_spec.kind)
    else
        model_spec.params
    end

    # Create individual model spec with this subject's doses
    ind_model_spec = ModelSpec(
        model_spec.kind,
        model_spec.name,
        individual_params,
        dose_events
    )

    # Create grid if not provided
    if grid === nothing
        t_start = isempty(observation_times) ? 0.0 : minimum(observation_times)
        t_end = isempty(observation_times) ? 24.0 : maximum(observation_times)
        saveat = isempty(observation_times) ? collect(0.0:0.5:t_end) : sort(unique(observation_times))
        grid = SimGrid(t_start, t_end, saveat)
    end

    # Create solver if not provided
    if solver === nothing
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)
    end

    # Run the actual simulation using the PK engine
    sim_result = simulate(ind_model_spec, grid, solver)

    # Extract concentrations
    concentrations = if haskey(sim_result.observations, :conc)
        sim_result.observations[:conc]
    else
        # Find any concentration-like key
        conc_key = nothing
        for k in keys(sim_result.observations)
            if contains(string(k), "conc") || contains(string(k), "Concentration")
                conc_key = k
                break
            end
        end
        conc_key !== nothing ? sim_result.observations[conc_key] : Float64[]
    end

    # Extract PD response if available
    pd_response = nothing
    for k in keys(sim_result.observations)
        if contains(string(k), "effect") || contains(string(k), "response")
            pd_response = sim_result.observations[k]
            break
        end
    end

    # Calculate PK metrics
    pk_metrics = calculate_pk_metrics(sim_result.t, concentrations)

    # Get subject ID
    subject_id = subject.id

    return SubjectExposure(
        subject_id,
        sim_result.t,
        concentrations,
        pd_response,
        pk_metrics
    )
end

"""
    simulate_subject_exposure(model_spec, dose_events, observation_times; grid, solver)

Simplified version without VirtualSubject (for population-typical simulation).

# Arguments
- `model_spec`: ModelSpec containing model type and parameters
- `dose_events`: Vector of DoseEvent
- `observation_times`: Times at which to observe concentrations

# Returns
- `SubjectExposure` with times, concentrations, and PK metrics
"""
function simulate_subject_exposure(
    model_spec::ModelSpec,
    dose_events::Vector{DoseEvent},
    observation_times::Vector{Float64};
    grid::Union{Nothing, SimGrid} = nothing,
    solver::Union{Nothing, SolverSpec} = nothing
)::SubjectExposure

    # Create model spec with doses
    ind_model_spec = ModelSpec(
        model_spec.kind,
        model_spec.name,
        model_spec.params,
        dose_events
    )

    # Create grid if not provided
    if grid === nothing
        t_start = isempty(observation_times) ? 0.0 : minimum(observation_times)
        t_end = isempty(observation_times) ? 24.0 : maximum(observation_times)
        saveat = isempty(observation_times) ? collect(0.0:0.5:t_end) : sort(unique(observation_times))
        grid = SimGrid(t_start, t_end, saveat)
    end

    # Create solver if not provided
    if solver === nothing
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)
    end

    # Run simulation
    sim_result = simulate(ind_model_spec, grid, solver)

    # Extract concentrations
    concentrations = if haskey(sim_result.observations, :conc)
        sim_result.observations[:conc]
    else
        conc_key = nothing
        for k in keys(sim_result.observations)
            if contains(string(k), "conc") || contains(string(k), "Concentration")
                conc_key = k
                break
            end
        end
        conc_key !== nothing ? sim_result.observations[conc_key] : Float64[]
    end

    # Calculate PK metrics
    pk_metrics = calculate_pk_metrics(sim_result.t, concentrations)

    return SubjectExposure(
        0,  # No subject ID
        sim_result.t,
        concentrations,
        nothing,
        pk_metrics
    )
end

"""
    generate_eta(rng::AbstractRNG, omega::Matrix{Float64})

Generate random effects (eta) from the omega variance-covariance matrix.

# Arguments
- `rng`: Random number generator
- `omega`: Variance-covariance matrix for IIV

# Returns
- Vector of random effects drawn from N(0, omega)
"""
function generate_eta(rng::AbstractRNG, omega::Matrix{Float64})
    n = size(omega, 1)
    if n == 0
        return Float64[]
    end

    # Cholesky decomposition for sampling
    L = cholesky(omega).L
    z = randn(rng, n)
    return L * z
end
