# Trial Execution Engine
# Main trial simulation functionality

export simulate_trial, simulate_trial_replicate, run_bioequivalence_analysis

using StableRNGs

"""
    simulate_trial(trial_spec::TrialSpec; grid=nothing, solver=nothing)

Simulate a clinical trial.

# Arguments
- `trial_spec::TrialSpec`: Trial specification

# Keyword Arguments
- `grid`: ODE solver grid (uses default if not specified)
- `solver`: ODE solver (uses default if not specified)

# Returns
- `TrialResult`: Complete trial result

# Example
```julia
spec = TrialSpec("Phase 1 Study", parallel_design(2), arms)
result = simulate_trial(spec)
```
"""
function simulate_trial(trial_spec::TrialSpec; grid = nothing, solver = nothing)
    rng = StableRNG(trial_spec.seed)

    # Generate virtual population
    total_subjects = sum(arm.n_subjects for arm in trial_spec.arms)
    population = generate_virtual_population(trial_spec.virtual_population, total_subjects)

    # Simulate enrollment
    enrollments = simulate_enrollment(trial_spec.enrollment, trial_spec.arms, population; rng = rng)

    # Apply compliance
    compliance_values = if trial_spec.compliance !== nothing
        apply_compliance(trial_spec.compliance, total_subjects, trial_spec.duration_days; rng = rng)
    else
        ones(total_subjects)
    end

    # Simulate adverse events (placeholder rate)
    aes = generate_adverse_events(total_subjects, trial_spec.duration_days, 0.05; rng = rng)

    # Simulate dropout
    dropouts = if trial_spec.dropout !== nothing
        simulate_dropout(trial_spec.dropout, trial_spec.duration_days, total_subjects;
                         adverse_events = aes, compliance_values = compliance_values, rng = rng)
    else
        DropoutEvent[]
    end

    # Group enrollments by arm
    arm_enrollments = Dict{String, Vector{EnrollmentEvent}}()
    for arm in trial_spec.arms
        arm_enrollments[arm.name] = filter(e -> e.arm_assignment == arm.name, enrollments)
    end

    # Simulate each arm
    arm_results = Dict{String, ArmResult}()

    for arm in trial_spec.arms
        arm_result = simulate_arm(arm, arm_enrollments[arm.name],
                                   dropouts, trial_spec; rng = rng,
                                   grid = grid, solver = solver)
        arm_results[arm.name] = arm_result
    end

    # Analyze endpoints - pass design for paired detection (crossover)
    endpoint_analyses = Dict{Symbol, Dict{Symbol, Any}}()
    for endpoint in trial_spec.endpoints
        endpoint_analyses[endpoint.name] = analyze_trial_endpoint(arm_results, endpoint, trial_spec.design)
    end

    # Power estimates (if replicates > 1)
    power_estimates = Dict{Symbol, Float64}()

    # BE results
    be_results = if trial_spec.design isa BioequivalenceDesign
        run_bioequivalence_analysis(arm_results, trial_spec.design)
    else
        nothing
    end

    return TrialResult(
        trial_spec.name,
        get_design_description(trial_spec.design),
        arm_results,
        endpoint_analyses,
        power_estimates,
        be_results,
        trial_spec.n_replicates,
        trial_spec.seed
    )
end


"""
    simulate_arm(arm::TreatmentArm, enrollments::Vector{EnrollmentEvent},
                 dropouts::Vector{DropoutEvent}, trial_spec::TrialSpec; kwargs...)

Simulate a single treatment arm.
"""
function simulate_arm(arm::TreatmentArm,
                       enrollments::Vector{EnrollmentEvent},
                       dropouts::Vector{DropoutEvent},
                       trial_spec::TrialSpec;
                       rng::AbstractRNG = StableRNG(12345),
                       grid = nothing, solver = nothing)

    individual_results = []
    endpoint_values = Dict{Symbol, Vector{Float64}}()

    n_enrolled = length(enrollments)
    n_completed = 0
    n_dropout = 0

    # Initialize endpoint vectors
    for endpoint in trial_spec.endpoints
        endpoint_values[endpoint.name] = Float64[]
    end

    for enrollment in enrollments
        subject_id = enrollment.subject_id
        covariates = enrollment.covariates

        # Check dropout
        survival_time, completed = calculate_survival_time(dropouts, subject_id,
                                                            trial_spec.duration_days)
        if completed
            n_completed += 1
        else
            n_dropout += 1
        end

        # Simulate PK/PD for this subject
        result = simulate_individual(arm, covariates, survival_time,
                                      trial_spec.pk_sampling_times;
                                      rng = rng, grid = grid, solver = solver)

        push!(individual_results, result)

        # Calculate endpoints for this subject
        for endpoint in trial_spec.endpoints
            value = calculate_individual_endpoint(result, endpoint, covariates)
            push!(endpoint_values[endpoint.name], value)
        end
    end

    # Calculate summary statistics
    summary_stats = Dict{Symbol, Dict{Symbol, Float64}}()
    for endpoint in trial_spec.endpoints
        values = endpoint_values[endpoint.name]
        valid = filter(!isnan, values)
        if !isempty(valid)
            summary_stats[endpoint.name] = Dict{Symbol, Float64}(
                :mean => sum(valid) / length(valid),
                :sd => length(valid) > 1 ? sqrt(sum((valid .- sum(valid)/length(valid)).^2) / (length(valid) - 1)) : 0.0,
                :n => Float64(length(valid))
            )
        end
    end

    return ArmResult(arm.name, n_enrolled, n_completed, n_dropout,
                     individual_results, endpoint_values, summary_stats)
end


"""
    simulate_individual(arm::TreatmentArm, covariates::VirtualSubject,
                        observation_duration::Float64, sampling_times::Vector{Float64}; kwargs...)

Simulate an individual subject using the real PK engine.
"""
function simulate_individual(arm::TreatmentArm,
                              covariates::VirtualSubject,
                              observation_duration::Float64,
                              sampling_times::Vector{Float64};
                              rng::AbstractRNG = StableRNG(12345),
                              grid = nothing, solver = nothing)

    # Get dose events
    dose_times = dose_event_times(arm.regimen)
    dose_amounts = generate_doses(arm.regimen; rng = rng)

    # Filter to observation duration and convert to DoseEvent objects
    dose_events = DoseEvent[]
    for (t, d) in zip(dose_times, dose_amounts)
        if t <= observation_duration && d > 0
            push!(dose_events, DoseEvent(t, d))
        end
    end

    # Create time grid for simulation
    max_time = min(observation_duration, maximum(sampling_times))
    t_grid = sort(unique(vcat(sampling_times, [d.time for d in dose_events])))
    t_grid = filter(t -> t <= max_time, t_grid)

    if isempty(t_grid)
        t_grid = [0.0]
    end

    # Check if we have a valid model spec to use the real PK engine
    if arm.pk_model_spec !== nothing && arm.pk_model_spec isa ModelSpec
        # Use the real PK simulation engine
        result = simulate_with_real_engine(arm.pk_model_spec, covariates, dose_events, t_grid, rng; grid=grid, solver=solver)
    else
        # Fallback to placeholder simulation when no model spec provided
        # This maintains backward compatibility
        valid_doses = [(d.time, d.amount) for d in dose_events]
        result = create_placeholder_simulation(t_grid, arm, covariates, valid_doses, rng)
    end

    return result
end

"""
    simulate_with_real_engine(model_spec, subject, dose_events, observation_times, rng; grid, solver)

Simulate using the real PK engine and return results in the expected Dict format.
"""
function simulate_with_real_engine(
    model_spec::ModelSpec,
    subject::VirtualSubject,
    dose_events::Vector{DoseEvent},
    observation_times::Vector{Float64},
    rng::AbstractRNG;
    grid = nothing,
    solver = nothing
)
    # Generate random effects for IIV
    n_params = get_n_params(model_spec.kind)
    eta = 0.3 .* randn(rng, n_params)  # 30% CV for IIV

    # Store eta in subject's eta field for IIV application
    subject_with_eta = VirtualSubject(
        subject.id,
        subject.age,
        subject.weight,
        subject.sex,
        subject.race,
        subject.disease_severity,
        subject.baseline_biomarker,
        subject.other,
        eta  # Use dedicated eta field
    )

    # Call the real simulation
    exposure = simulate_subject_exposure(
        model_spec,
        subject_with_eta,
        dose_events,
        observation_times;
        grid = grid,
        solver = solver,
        include_iiv = true
    )

    # Convert SubjectExposure to Dict format for compatibility with existing endpoint code
    return Dict{String, Any}(
        "t" => exposure.times,
        "observations" => Dict{String, Vector{Float64}}(
            "conc" => exposure.concentrations,
            "effect" => exposure.pd_response !== nothing ? exposure.pd_response : zeros(length(exposure.times))
        ),
        "covariates" => Dict{Symbol, Any}(
            :weight => subject.weight,
            :age => subject.age,
            :sex => subject.sex
        ),
        "pk_metrics" => exposure.pk_metrics
    )
end

"""
Get the number of parameters for a model type (for IIV generation).
"""
function get_n_params(::OneCompIVBolus)
    return 2  # CL, V
end

function get_n_params(::OneCompOralFirstOrder)
    return 3  # Ka, CL, V
end

function get_n_params(::TwoCompIVBolus)
    return 4  # CL, V1, Q, V2
end

function get_n_params(::TwoCompOral)
    return 5  # Ka, CL, V1, Q, V2
end

function get_n_params(::ModelKind)
    return 3  # Default fallback
end

"""
Create a placeholder simulation result for testing.
Used when no ModelSpec is provided (backward compatibility).
"""
function create_placeholder_simulation(t_grid::Vector{Float64},
                                        arm::TreatmentArm,
                                        covariates::VirtualSubject,
                                        doses::Vector{Tuple{Float64, Float64}},
                                        rng::AbstractRNG)

    n_points = length(t_grid)

    # Simple exponential decay PK model
    conc = zeros(n_points)
    effect = zeros(n_points)

    # Individual variability factors
    cl_factor = exp(0.3 * randn(rng))  # 30% CV
    v_factor = exp(0.25 * randn(rng))  # 25% CV

    # Weight-based adjustment
    weight_ref = 70.0
    weight_factor = (covariates.weight / weight_ref)^0.75

    cl = 10.0 * cl_factor * weight_factor  # L/h
    v = 100.0 * v_factor * weight_factor  # L
    ka = 1.5  # h^-1

    ke = cl / v

    for (i, t) in enumerate(t_grid)
        c_total = 0.0
        for (dose_time, dose_amount) in doses
            if t >= dose_time && dose_amount > 0
                dt = t - dose_time
                # One-compartment absorption
                c = (dose_amount / v) * (ka / (ka - ke)) *
                    (exp(-ke * dt) - exp(-ka * dt))
                c_total += max(0.0, c)
            end
        end
        conc[i] = c_total

        # Simple Emax PD model
        emax = 100.0
        ec50 = 5.0
        e0 = 20.0 + 10.0 * randn(rng)  # Baseline with variability
        effect[i] = e0 + emax * conc[i] / (ec50 + conc[i])
    end

    return Dict{String, Any}(
        "t" => t_grid,
        "observations" => Dict{String, Vector{Float64}}(
            "conc" => conc,
            "effect" => effect
        ),
        "covariates" => Dict{Symbol, Any}(
            :weight => covariates.weight,
            :age => covariates.age,
            :sex => covariates.sex
        )
    )
end


"""
    calculate_individual_endpoint(result, endpoint::EndpointSpec,
                                   covariates::VirtualSubject)

Calculate endpoint value for an individual.
"""
function calculate_individual_endpoint(result, endpoint::PKEndpoint,
                                         covariates::VirtualSubject)
    return extract_pk_endpoint(result, endpoint)
end

function calculate_individual_endpoint(result, endpoint::PDEndpoint,
                                         covariates::VirtualSubject)
    baseline = covariates.baseline_biomarker
    return extract_pd_endpoint(result, endpoint, baseline)
end

function calculate_individual_endpoint(result, endpoint::SafetyEndpoint,
                                         covariates::VirtualSubject)
    return NaN  # Safety endpoints handled at population level
end

function calculate_individual_endpoint(result, endpoint::CompositeEndpoint,
                                         covariates::VirtualSubject)
    # Calculate each component
    values = [calculate_individual_endpoint(result, comp, covariates)
              for comp in endpoint.components]

    if endpoint.combination_rule == :any
        # Any component meets criteria
        return any(!isnan(v) && v > 0 for v in values) ? 1.0 : 0.0
    elseif endpoint.combination_rule == :all
        # All components meet criteria
        return all(!isnan(v) && v > 0 for v in values) ? 1.0 : 0.0
    elseif endpoint.combination_rule == :weighted
        # Weighted average
        valid = [(v, w) for (v, w) in zip(values, endpoint.weights) if !isnan(v)]
        if isempty(valid)
            return NaN
        end
        return sum(v * w for (v, w) in valid) / sum(w for (_, w) in valid)
    end

    return NaN
end


"""
    analyze_trial_endpoint(arm_results::Dict{String, ArmResult}, endpoint::EndpointSpec, design=nothing)

Analyze an endpoint across all arms, with automatic paired analysis detection for crossover designs.

# Arguments
- `arm_results`: Results by arm name
- `endpoint`: Endpoint specification
- `design`: Study design (optional). If CrossoverDesign, uses paired analysis.

# Returns
- Dict with summary statistics and comparisons (paired if crossover)
"""
function analyze_trial_endpoint(arm_results::Dict{String, ArmResult},
                                  endpoint::EndpointSpec,
                                  design::Union{Nothing, StudyDesignKind} = nothing)

    result = Dict{Symbol, Any}()

    # Auto-detect if this is a crossover design requiring paired analysis
    is_crossover = design isa CrossoverDesign
    is_be_design = design isa BioequivalenceDesign

    result[:design_type] = if is_crossover
        :crossover
    elseif is_be_design
        :bioequivalence
    else
        :parallel
    end

    arm_names = collect(keys(arm_results))

    # Collect values by arm
    arm_values = Dict{String, Vector{Float64}}()
    for (name, arm_result) in arm_results
        if haskey(arm_result.endpoint_values, endpoint.name)
            arm_values[name] = arm_result.endpoint_values[endpoint.name]
        else
            arm_values[name] = Float64[]
        end
    end

    # Summary by arm
    result[:by_arm] = Dict{String, Dict{Symbol, Any}}()
    for (name, values) in arm_values
        result[:by_arm][name] = analyze_endpoint(values, endpoint)
    end

    # For crossover designs, use specialized crossover analysis
    if (is_crossover || is_be_design) && length(arm_names) == 2
        # Use crossover-specific analysis with paired comparisons
        crossover_design = is_crossover ? design : design.crossover
        result[:crossover_analysis] = analyze_crossover_endpoint(arm_results, endpoint, crossover_design)
        result[:is_paired] = true

        # Also include standard comparison format for backward compatibility
        name1, name2 = arm_names[1], arm_names[2]
        key = "$(name1)_vs_$(name2)"
        result[:comparisons] = Dict{String, Dict{Symbol, Any}}()

        test_type = endpoint isa PKEndpoint && endpoint.log_transform ? :ratio : :ttest

        result[:comparisons][key] = compare_arms(
            arm_values[name1], arm_values[name2];
            test = test_type,
            paired = true  # Use paired analysis for crossover
        )
    else
        # Standard parallel group analysis
        result[:is_paired] = false

        if length(arm_names) >= 2
            result[:comparisons] = Dict{String, Dict{Symbol, Any}}()

            for i in 1:(length(arm_names)-1)
                for j in (i+1):length(arm_names)
                    name1, name2 = arm_names[i], arm_names[j]
                    key = "$(name1)_vs_$(name2)"

                    test_type = endpoint isa PKEndpoint && endpoint.log_transform ?
                        :ratio : :ttest

                    result[:comparisons][key] = compare_arms(
                        arm_values[name1], arm_values[name2];
                        test = test_type,
                        paired = false
                    )
                end
            end
        end
    end

    return result
end

"""
    run_bioequivalence_analysis(arm_results::Dict{String, ArmResult},
                                 design::BioequivalenceDesign)

Run bioequivalence analysis.
"""
function run_bioequivalence_analysis(arm_results::Dict{String, ArmResult},
                                       design::BioequivalenceDesign)

    arm_names = collect(keys(arm_results))
    if length(arm_names) != 2
        return Dict{Symbol, Any}(:error => "BE requires exactly 2 arms")
    end

    test_arm = arm_names[1]
    ref_arm = arm_names[2]

    result = Dict{Symbol, Any}()

    for param in design.parameters
        param_key = param

        test_values = Float64[]
        ref_values = Float64[]

        # Collect values
        if haskey(arm_results[test_arm].endpoint_values, param_key)
            test_values = filter(!isnan, arm_results[test_arm].endpoint_values[param_key])
        end
        if haskey(arm_results[ref_arm].endpoint_values, param_key)
            ref_values = filter(!isnan, arm_results[ref_arm].endpoint_values[param_key])
        end

        if isempty(test_values) || isempty(ref_values)
            result[param] = Dict{Symbol, Any}(:error => "Insufficient data")
            continue
        end

        # Calculate geometric mean ratio and 90% CI
        pos_test = filter(v -> v > 0, test_values)
        pos_ref = filter(v -> v > 0, ref_values)

        if isempty(pos_test) || isempty(pos_ref)
            result[param] = Dict{Symbol, Any}(:error => "No positive values")
            continue
        end

        log_test = log.(pos_test)
        log_ref = log.(pos_ref)

        n_test, n_ref = length(log_test), length(log_ref)
        mean_log_test = sum(log_test) / n_test
        mean_log_ref = sum(log_ref) / n_ref

        var_log_test = n_test > 1 ? sum((log_test .- mean_log_test).^2) / (n_test - 1) : 0.0
        var_log_ref = n_ref > 1 ? sum((log_ref .- mean_log_ref).^2) / (n_ref - 1) : 0.0

        # Pooled variance
        pooled_var = ((n_test - 1) * var_log_test + (n_ref - 1) * var_log_ref) /
                     (n_test + n_ref - 2)
        se = sqrt(pooled_var * (1/n_test + 1/n_ref))

        # GMR and 90% CI
        diff_log = mean_log_test - mean_log_ref
        t_90 = 1.645  # Approximate for 90% CI

        gmr = exp(diff_log)
        ci_lower = exp(diff_log - t_90 * se)
        ci_upper = exp(diff_log + t_90 * se)

        # Check BE criteria
        be_lower, be_upper = design.bioequivalence_limits
        is_be = ci_lower >= be_lower && ci_upper <= be_upper

        result[param] = Dict{Symbol, Any}(
            :gmr => gmr,
            :ci_90_lower => ci_lower,
            :ci_90_upper => ci_upper,
            :bioequivalent => is_be,
            :n_test => n_test,
            :n_ref => n_ref
        )
    end

    # Overall BE conclusion
    all_params_be = all(haskey(result[p], :bioequivalent) && result[p][:bioequivalent]
                        for p in design.parameters if haskey(result, p))
    result[:overall_bioequivalent] = all_params_be

    return result
end


"""
    simulate_trial_replicate(trial_spec::TrialSpec, replicate_id::Int; kwargs...)

Simulate one replicate of a trial for power analysis.
"""
function simulate_trial_replicate(trial_spec::TrialSpec, replicate_id::Int;
                                   grid = nothing, solver = nothing)
    # Modify seed for this replicate
    modified_spec = TrialSpec(
        trial_spec.name,
        trial_spec.design,
        trial_spec.arms;
        virtual_population = VirtualPopulationSpec(
            demographics = trial_spec.virtual_population.demographics,
            disease = trial_spec.virtual_population.disease,
            covariate_correlations = trial_spec.virtual_population.covariate_correlations,
            seed = trial_spec.seed + UInt64(replicate_id * 1000)
        ),
        duration_days = trial_spec.duration_days,
        enrollment = trial_spec.enrollment,
        dropout = trial_spec.dropout,
        compliance = trial_spec.compliance,
        pk_sampling_times = trial_spec.pk_sampling_times,
        endpoints = trial_spec.endpoints,
        n_replicates = 1,
        seed = trial_spec.seed + UInt64(replicate_id)
    )

    return simulate_trial(modified_spec; grid = grid, solver = solver)
end


"""
    run_power_simulation(trial_spec::TrialSpec; kwargs...)

Run power simulation with multiple replicates.
"""
function run_power_simulation(trial_spec::TrialSpec;
                               alpha::Float64 = 0.05,
                               grid = nothing, solver = nothing)

    if trial_spec.n_replicates <= 1
        return simulate_trial(trial_spec; grid = grid, solver = solver)
    end

    n_sig = Dict{Symbol, Int}()
    for endpoint in trial_spec.endpoints
        n_sig[endpoint.name] = 0
    end

    for rep in 1:trial_spec.n_replicates
        result = simulate_trial_replicate(trial_spec, rep; grid = grid, solver = solver)

        # Check significance for each endpoint
        for (endpoint_name, analysis) in result.endpoint_analyses
            if haskey(analysis, :comparisons)
                for (_, comp) in analysis[:comparisons]
                    if haskey(comp, :ci_lower) && haskey(comp, :ci_upper)
                        # Significant if CI doesn't include 0
                        if comp[:ci_lower] > 0 || comp[:ci_upper] < 0
                            n_sig[endpoint_name] += 1
                            break
                        end
                    end
                end
            end
        end
    end

    # Calculate power
    power_estimates = Dict{Symbol, Float64}()
    for (name, n) in n_sig
        power_estimates[name] = n / trial_spec.n_replicates
    end

    # Return final result with power estimates
    final_result = simulate_trial(trial_spec; grid = grid, solver = solver)

    return TrialResult(
        final_result.trial_name,
        final_result.design_type,
        final_result.arms,
        final_result.endpoint_analyses,
        power_estimates,
        final_result.be_results,
        trial_spec.n_replicates,
        trial_spec.seed
    )
end
