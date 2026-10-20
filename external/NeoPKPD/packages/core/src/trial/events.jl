# Trial Events
# Dropout, compliance, and enrollment event handling

export simulate_dropout, simulate_enrollment, apply_compliance
export calculate_survival_time, DropoutEvent, EnrollmentEvent

using StableRNGs
using Random

"""
    DropoutEvent

Record of a dropout event.

# Fields
- `subject_id::Int`: Subject ID
- `time_days::Float64`: Time of dropout
- `reason::Symbol`: Dropout reason (:random, :ae, :non_compliance, :lost_to_followup)
"""
struct DropoutEvent
    subject_id::Int
    time_days::Float64
    reason::Symbol
end


"""
    EnrollmentEvent

Record of an enrollment event.

# Fields
- `subject_id::Int`: Subject ID
- `enrollment_day::Float64`: Day of enrollment
- `arm_assignment::String`: Treatment arm assignment
- `covariates::VirtualSubject`: Subject covariates
"""
struct EnrollmentEvent
    subject_id::Int
    enrollment_day::Float64
    arm_assignment::String
    covariates::VirtualSubject
end


"""
    simulate_dropout(spec::DropoutSpec, duration_days::Float64, n_subjects::Int; kwargs...)

Simulate dropout events for a cohort of subjects.

# Arguments
- `spec::DropoutSpec`: Dropout specification
- `duration_days::Float64`: Trial duration in days
- `n_subjects::Int`: Number of subjects

# Keyword Arguments
- `adverse_events::Vector{Tuple{Int, Float64, Float64}}`: (subject_id, time, severity)
- `compliance_values::Vector{Float64}`: Individual compliance rates
- `rng::AbstractRNG`: Random number generator

# Returns
- `Vector{DropoutEvent}`: Dropout events

# Example
```julia
spec = DropoutSpec(random_rate_per_day=0.002)
dropouts = simulate_dropout(spec, 28.0, 100)
```
"""
function simulate_dropout(spec::DropoutSpec, duration_days::Float64, n_subjects::Int;
                           adverse_events::Vector{Tuple{Int, Float64, Float64}} = Tuple{Int, Float64, Float64}[],
                           compliance_values::Union{Nothing, Vector{Float64}} = nothing,
                           rng::Union{Nothing, AbstractRNG} = nothing)
    if rng === nothing
        rng = StableRNG(12345)
    end

    dropouts = DropoutEvent[]

    for i in 1:n_subjects
        dropout_time = nothing
        dropout_reason = :none

        # Random dropout (exponential survival)
        if spec.random_rate_per_day > 0
            random_survival = -log(rand(rng)) / spec.random_rate_per_day
            if random_survival < duration_days
                if dropout_time === nothing || random_survival < dropout_time
                    dropout_time = random_survival
                    dropout_reason = :random
                end
            end
        end

        # AE-related dropout
        if spec.ae_threshold !== nothing
            subject_aes = filter(ae -> ae[1] == i, adverse_events)
            for (_, ae_time, ae_severity) in subject_aes
                if ae_severity >= spec.ae_threshold
                    if rand(rng) < spec.ae_dropout_prob
                        if dropout_time === nothing || ae_time < dropout_time
                            dropout_time = ae_time
                            dropout_reason = :ae
                        end
                    end
                end
            end
        end

        # Non-compliance dropout
        if compliance_values !== nothing && i <= length(compliance_values)
            if compliance_values[i] < spec.non_compliance_threshold
                # Higher chance of early dropout for non-compliant subjects
                nc_time = rand(rng) * duration_days * 0.5
                if dropout_time === nothing || nc_time < dropout_time
                    dropout_time = nc_time
                    dropout_reason = :non_compliance
                end
            end
        end

        if dropout_time !== nothing
            push!(dropouts, DropoutEvent(i, dropout_time, dropout_reason))
        end
    end

    return dropouts
end


"""
    calculate_survival_time(dropout_events, subject_id, duration_days)

Calculate survival time (time in trial) for a subject.

# Arguments
- `dropout_events::Vector{DropoutEvent}`: Dropout events
- `subject_id::Int`: Subject ID
- `duration_days::Float64`: Trial duration

# Returns
- `Tuple{Float64, Bool}`: (survival_time, completed)
"""
function calculate_survival_time(dropout_events::Vector{DropoutEvent},
                                  subject_id::Int, duration_days::Float64)
    for event in dropout_events
        if event.subject_id == subject_id
            return (event.time_days, false)
        end
    end
    return (duration_days, true)
end


"""
    simulate_enrollment(spec::EnrollmentSpec, arms::Vector{TreatmentArm},
                        population::Vector{VirtualSubject}; kwargs...)

Simulate enrollment of subjects into a trial.

# Arguments
- `spec::EnrollmentSpec`: Enrollment specification
- `arms::Vector{TreatmentArm}`: Treatment arms
- `population::Vector{VirtualSubject}`: Virtual population

# Keyword Arguments
- `rng::AbstractRNG`: Random number generator

# Returns
- `Vector{EnrollmentEvent}`: Enrollment events

# Example
```julia
spec = EnrollmentSpec(rate_per_day=5.0)
enrollments = simulate_enrollment(spec, arms, population)
```
"""
function simulate_enrollment(spec::EnrollmentSpec,
                              arms::Vector{TreatmentArm},
                              population::Vector{VirtualSubject};
                              rng::Union{Nothing, AbstractRNG} = nothing)
    if rng === nothing
        rng = StableRNG(12345)
    end

    # Calculate total subjects needed
    total_needed = sum(arm.n_subjects for arm in arms)
    total_available = min(length(population), spec.enrollment_cap)

    # Account for screening failures
    total_to_screen = ceil(Int, total_needed / (1 - spec.screening_failure_rate))
    total_to_screen = min(total_to_screen, total_available)

    enrollments = EnrollmentEvent[]
    current_day = 0.0
    enrolled_count = 0
    arm_counts = Dict(arm.name => 0 for arm in arms)

    # Randomization probabilities based on target enrollment
    arm_probs = [arm.n_subjects for arm in arms]
    arm_probs = arm_probs / sum(arm_probs)

    pop_idx = 1
    while enrolled_count < total_needed && pop_idx <= total_to_screen

        # Time to next enrollment (Poisson process)
        inter_arrival = -log(rand(rng)) / spec.rate_per_day
        current_day += inter_arrival

        # Check concurrent enrollment limit
        concurrent = count(e -> current_day - e.enrollment_day < 1.0, enrollments)
        if concurrent >= spec.max_concurrent
            continue
        end

        # Screening
        if rand(rng) < spec.screening_failure_rate
            pop_idx += 1
            continue
        end

        # Find arm that needs subjects
        available_arms = filter(arm -> arm_counts[arm.name] < arm.n_subjects, arms)
        if isempty(available_arms)
            break
        end

        # Randomize to arm
        avail_probs = [arm.n_subjects - arm_counts[arm.name] for arm in available_arms]
        avail_probs = avail_probs / sum(avail_probs)

        u = rand(rng)
        cumsum = 0.0
        selected_arm = available_arms[1]
        for (arm, prob) in zip(available_arms, avail_probs)
            cumsum += prob
            if u <= cumsum
                selected_arm = arm
                break
            end
        end

        # Enroll subject
        covariates = population[pop_idx]
        push!(enrollments, EnrollmentEvent(
            pop_idx,
            current_day,
            selected_arm.name,
            covariates
        ))

        arm_counts[selected_arm.name] += 1
        enrolled_count += 1
        pop_idx += 1
    end

    return enrollments
end


"""
    apply_compliance(spec::ComplianceSpec, n_subjects::Int, duration_days::Float64; kwargs...)

Generate individual compliance values.

# Arguments
- `spec::ComplianceSpec`: Compliance specification
- `n_subjects::Int`: Number of subjects
- `duration_days::Float64`: Trial duration

# Keyword Arguments
- `rng::AbstractRNG`: Random number generator

# Returns
- `Vector{Float64}`: Compliance rates for each subject

# Example
```julia
spec = ComplianceSpec(mean_compliance=0.85, compliance_sd=0.10)
compliance = apply_compliance(spec, 100, 28.0)
```
"""
function apply_compliance(spec::ComplianceSpec, n_subjects::Int, duration_days::Float64;
                           rng::Union{Nothing, AbstractRNG} = nothing)
    if rng === nothing
        rng = StableRNG(12345)
    end

    compliance = zeros(n_subjects)

    for i in 1:n_subjects
        base_compliance = spec.mean_compliance + spec.compliance_sd * randn(rng)
        base_compliance = clamp(base_compliance, spec.min_compliance, 1.0)

        # Apply pattern effects
        if spec.pattern == :random
            compliance[i] = base_compliance
        elseif spec.pattern == :decay
            # Compliance decays - use average over trial
            decay_factor = 0.7 + 0.3 * rand(rng)  # 70-100% of initial
            compliance[i] = base_compliance * decay_factor
        elseif spec.pattern == :early_good
            # Better early, worse later - depends on individual
            if rand(rng) < 0.3  # 30% have early good pattern
                compliance[i] = base_compliance * 1.1
                compliance[i] = min(compliance[i], 1.0)
            else
                compliance[i] = base_compliance
            end
        elseif spec.pattern == :weekend_miss
            # Approximately 2/7 doses potentially missed on weekends
            weekend_effect = 2/7 * (1 - base_compliance)
            compliance[i] = base_compliance - weekend_effect
            compliance[i] = max(compliance[i], spec.min_compliance)
        else
            compliance[i] = base_compliance
        end
    end

    return compliance
end


"""
    generate_adverse_events(n_subjects::Int, duration_days::Float64,
                            base_rate::Float64; kwargs...)

Generate adverse events for trial subjects.

# Arguments
- `n_subjects::Int`: Number of subjects
- `duration_days::Float64`: Trial duration
- `base_rate::Float64`: Base AE rate per subject per day

# Keyword Arguments
- `severity_probs::Vector{Float64}`: Probabilities for mild, moderate, severe
- `rng::AbstractRNG`: Random number generator

# Returns
- `Vector{Tuple{Int, Float64, Float64}}`: (subject_id, time, severity)
"""
function generate_adverse_events(n_subjects::Int, duration_days::Float64,
                                  base_rate::Float64;
                                  severity_probs::Vector{Float64} = [0.6, 0.3, 0.1],
                                  rng::Union{Nothing, AbstractRNG} = nothing)
    if rng === nothing
        rng = StableRNG(12345)
    end

    aes = Tuple{Int, Float64, Float64}[]

    for i in 1:n_subjects
        # Number of AEs for this subject (Poisson)
        expected_aes = base_rate * duration_days
        n_aes = rand(rng) < 1 - exp(-expected_aes) ? 1 : 0  # Simplified

        if expected_aes > 1
            # More sophisticated Poisson sampling
            n_aes = 0
            for _ in 1:100
                if rand(rng) < expected_aes / 100
                    n_aes += 1
                end
            end
        end

        for _ in 1:n_aes
            ae_time = rand(rng) * duration_days

            # Severity (1=mild, 2=moderate, 3=severe)
            u = rand(rng)
            if u < severity_probs[1]
                severity = 1.0
            elseif u < severity_probs[1] + severity_probs[2]
                severity = 2.0
            else
                severity = 3.0
            end

            push!(aes, (i, ae_time, severity))
        end
    end

    return aes
end


"""
    stratified_randomization(population::Vector{VirtualSubject},
                              arms::Vector{TreatmentArm},
                              factors::Vector{Symbol}; rng=nothing)

Perform stratified randomization.

# Arguments
- `population::Vector{VirtualSubject}`: Population to randomize
- `arms::Vector{TreatmentArm}`: Treatment arms
- `factors::Vector{Symbol}`: Stratification factors (:age_group, :sex, :race)

# Keyword Arguments
- `rng::AbstractRNG`: Random number generator

# Returns
- `Dict{Int, String}`: Subject ID to arm assignment mapping
"""
function stratified_randomization(population::Vector{VirtualSubject},
                                   arms::Vector{TreatmentArm},
                                   factors::Vector{Symbol};
                                   rng::Union{Nothing, AbstractRNG} = nothing)
    if rng === nothing
        rng = StableRNG(12345)
    end

    n = length(population)
    arm_names = [arm.name for arm in arms]
    n_arms = length(arms)

    assignments = Dict{Int, String}()

    # Create strata
    function get_stratum(ind::VirtualSubject)
        parts = String[]
        for f in factors
            if f == :age_group
                age_group = ind.age < 40 ? "young" : (ind.age < 60 ? "middle" : "elderly")
                push!(parts, age_group)
            elseif f == :sex
                push!(parts, string(ind.sex))
            elseif f == :race
                push!(parts, string(ind.race))
            elseif f == :disease_severity && ind.disease_severity !== nothing
                push!(parts, string(ind.disease_severity))
            end
        end
        return join(parts, "_")
    end

    # Block randomization within strata
    strata_counts = Dict{String, Dict{String, Int}}()

    for ind in population
        stratum = get_stratum(ind)

        if !haskey(strata_counts, stratum)
            strata_counts[stratum] = Dict(name => 0 for name in arm_names)
        end

        # Find arm with lowest count in this stratum
        min_count = minimum(values(strata_counts[stratum]))
        available_arms = [name for name in arm_names
                          if strata_counts[stratum][name] == min_count]

        # Random selection among tied arms
        selected_arm = available_arms[rand(rng, 1:length(available_arms))]

        assignments[ind.id] = selected_arm
        strata_counts[stratum][selected_arm] += 1
    end

    return assignments
end
