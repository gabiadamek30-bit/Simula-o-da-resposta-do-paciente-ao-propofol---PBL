# Dosing Regimens
# Functions for creating dosing regimens

export dosing_qd, dosing_bid, dosing_tid, dosing_qid, dosing_custom
export titration_regimen, dose_event_times, generate_doses
export total_regimen_duration


"""
    dosing_qd(dose, duration_days; loading_dose=nothing)

Create a once-daily dosing regimen.

# Arguments
- `dose::Float64`: Daily dose amount
- `duration_days::Int`: Treatment duration in days

# Keyword Arguments
- `loading_dose::Float64`: Optional loading dose on day 1

# Returns
- `DosingRegimen`: QD dosing regimen

# Example
```julia
regimen = dosing_qd(100.0, 28)
regimen = dosing_qd(100.0, 28, loading_dose=200.0)
```
"""
function dosing_qd(dose::Float64, duration_days::Int;
                   loading_dose::Union{Float64, Nothing} = nothing)
    return DosingRegimen(QD(), dose, duration_days; loading_dose = loading_dose)
end


"""
    dosing_bid(dose, duration_days; loading_dose=nothing)

Create a twice-daily dosing regimen.

# Arguments
- `dose::Float64`: Dose amount per administration
- `duration_days::Int`: Treatment duration in days

# Keyword Arguments
- `loading_dose::Float64`: Optional loading dose for first administration

# Returns
- `DosingRegimen`: BID dosing regimen

# Example
```julia
regimen = dosing_bid(50.0, 14)  # 50mg twice daily for 14 days
```
"""
function dosing_bid(dose::Float64, duration_days::Int;
                    loading_dose::Union{Float64, Nothing} = nothing)
    return DosingRegimen(BID(), dose, duration_days; loading_dose = loading_dose)
end


"""
    dosing_tid(dose, duration_days; loading_dose=nothing)

Create a three-times-daily dosing regimen.

# Arguments
- `dose::Float64`: Dose amount per administration
- `duration_days::Int`: Treatment duration in days

# Keyword Arguments
- `loading_dose::Float64`: Optional loading dose

# Returns
- `DosingRegimen`: TID dosing regimen

# Example
```julia
regimen = dosing_tid(25.0, 7)  # 25mg three times daily for 7 days
```
"""
function dosing_tid(dose::Float64, duration_days::Int;
                    loading_dose::Union{Float64, Nothing} = nothing)
    return DosingRegimen(TID(), dose, duration_days; loading_dose = loading_dose)
end


"""
    dosing_qid(dose, duration_days; loading_dose=nothing)

Create a four-times-daily dosing regimen.

# Arguments
- `dose::Float64`: Dose amount per administration
- `duration_days::Int`: Treatment duration in days

# Keyword Arguments
- `loading_dose::Float64`: Optional loading dose

# Returns
- `DosingRegimen`: QID dosing regimen

# Example
```julia
regimen = dosing_qid(25.0, 5)  # 25mg four times daily for 5 days
```
"""
function dosing_qid(dose::Float64, duration_days::Int;
                    loading_dose::Union{Float64, Nothing} = nothing)
    return DosingRegimen(QID(), dose, duration_days; loading_dose = loading_dose)
end


"""
    dosing_custom(dose, duration_days, dose_times; loading_dose=nothing)

Create a custom dosing regimen with specified dose times.

# Arguments
- `dose::Float64`: Dose amount per administration
- `duration_days::Int`: Treatment duration in days
- `dose_times::Vector{Float64}`: Times of administration (hours after midnight)

# Keyword Arguments
- `loading_dose::Float64`: Optional loading dose

# Returns
- `DosingRegimen`: Custom dosing regimen

# Example
```julia
# Dose at 6 AM and 10 PM
regimen = dosing_custom(75.0, 21, [6.0, 22.0])
```
"""
function dosing_custom(dose::Float64, duration_days::Int, dose_times::Vector{Float64};
                       loading_dose::Union{Float64, Nothing} = nothing)
    freq = CustomFrequency(dose_times)
    return DosingRegimen(freq, dose, duration_days; loading_dose = loading_dose)
end


"""
    titration_regimen(start_dose, target_dose, n_steps, days_per_step; kwargs...)

Create a titration dosing regimen.

# Arguments
- `start_dose::Float64`: Starting dose
- `target_dose::Float64`: Target maintenance dose
- `n_steps::Int`: Number of titration steps
- `days_per_step::Int`: Days at each dose level

# Keyword Arguments
- `frequency::DosingFrequency`: Dosing frequency (default: QD)
- `loading_dose::Float64`: Optional loading dose
- `maintenance_days::Int`: Days at maintenance dose (default: 0)

# Returns
- `TitrationRegimen`: Titration regimen

# Example
```julia
# Titrate from 25mg to 100mg in 4 steps, 7 days each
regimen = titration_regimen(25.0, 100.0, 4, 7)

# With 28 days maintenance at target dose
regimen = titration_regimen(25.0, 100.0, 4, 7, maintenance_days=28)
```
"""
function titration_regimen(start_dose::Float64, target_dose::Float64,
                           n_steps::Int, days_per_step::Int;
                           frequency::DosingFrequency = QD(),
                           loading_dose::Union{Float64, Nothing} = nothing,
                           maintenance_days::Int = 0)

    # Calculate dose increment
    dose_increment = (target_dose - start_dose) / (n_steps - 1)

    # Create titration steps
    steps = TitrationStep[]
    for i in 1:n_steps
        dose = start_dose + (i - 1) * dose_increment
        push!(steps, TitrationStep(dose, days_per_step))
    end

    # Add maintenance phase if specified
    if maintenance_days > 0
        push!(steps, TitrationStep(target_dose, maintenance_days))
    end

    return TitrationRegimen(steps; frequency = frequency, loading_dose = loading_dose)
end


"""
    dose_event_times(regimen::DosingRegimen)

Calculate all dose event times for a regimen.

# Arguments
- `regimen::DosingRegimen`: Dosing regimen

# Returns
- `Vector{Float64}`: Dose times in hours from start

# Example
```julia
regimen = dosing_bid(50.0, 3)
times = dose_event_times(regimen)  # [8, 20, 32, 44, 56, 68] hours
```
"""
function dose_event_times(regimen::DosingRegimen)
    times = Float64[]
    dose_times = regimen.dose_times

    for day in 0:(regimen.duration_days - 1)
        for t in dose_times
            push!(times, day * 24.0 + t)
        end
    end

    return times
end

function dose_event_times(regimen::TitrationRegimen)
    times = Float64[]
    dose_times = get_dose_times(regimen.frequency)

    current_time = 0.0
    for step in regimen.steps
        for day in 0:(step.duration_days - 1)
            for t in dose_times
                push!(times, current_time + day * 24.0 + t)
            end
        end
        current_time += step.duration_days * 24.0
    end

    return times
end


"""
    generate_doses(regimen::DosingRegimen; compliance=nothing, rng=nothing)

Generate dose amounts for each dose event, accounting for compliance.

# Arguments
- `regimen::DosingRegimen`: Dosing regimen

# Keyword Arguments
- `compliance::ComplianceSpec`: Compliance specification (optional)
- `rng::AbstractRNG`: Random number generator

# Returns
- `Vector{Float64}`: Dose amounts (0 for missed doses)

# Example
```julia
regimen = dosing_qd(100.0, 7)
doses = generate_doses(regimen)  # [100, 100, 100, 100, 100, 100, 100]

# With compliance
compliance = ComplianceSpec(mean_compliance=0.8)
doses = generate_doses(regimen, compliance=compliance, rng=StableRNG(42))
```
"""
function generate_doses(regimen::DosingRegimen;
                        compliance::Union{Nothing, ComplianceSpec} = nothing,
                        rng::Union{Nothing, AbstractRNG} = nothing)
    if rng === nothing
        rng = StableRNG(12345)
    end

    times = dose_event_times(regimen)
    n_doses = length(times)
    doses = fill(regimen.dose_amount, n_doses)

    # Apply loading dose
    if regimen.loading_dose !== nothing
        doses[1] = regimen.loading_dose
    end

    # Apply compliance
    if compliance !== nothing
        for i in 1:n_doses
            if compliance.pattern == :random
                # Random compliance - each dose has probability of being taken
                if rand(rng) > compliance.mean_compliance
                    doses[i] = 0.0
                end
            elseif compliance.pattern == :weekend_miss
                # Higher miss rate on weekends (days 6 and 7)
                day = floor(Int, times[i] / 24.0) % 7 + 1
                miss_prob = day >= 6 ? (1 - compliance.mean_compliance) * 2 : (1 - compliance.mean_compliance)
                if rand(rng) < miss_prob
                    doses[i] = 0.0
                end
            elseif compliance.pattern == :decay
                # Compliance decays over time
                relative_time = i / n_doses
                current_compliance = compliance.mean_compliance * (1 - 0.3 * relative_time)
                if rand(rng) > current_compliance
                    doses[i] = 0.0
                end
            elseif compliance.pattern == :early_good
                # Better compliance early, worse later
                relative_time = i / n_doses
                current_compliance = relative_time < 0.5 ?
                    compliance.mean_compliance + 0.1 :
                    compliance.mean_compliance - 0.1
                current_compliance = clamp(current_compliance, 0.0, 1.0)
                if rand(rng) > current_compliance
                    doses[i] = 0.0
                end
            end
        end
    end

    return doses
end

function generate_doses(regimen::TitrationRegimen;
                        compliance::Union{Nothing, ComplianceSpec} = nothing,
                        rng::Union{Nothing, AbstractRNG} = nothing)
    if rng === nothing
        rng = StableRNG(12345)
    end

    times = dose_event_times(regimen)
    n_doses = length(times)
    doses = Float64[]

    dose_times = get_dose_times(regimen.frequency)
    doses_per_day = length(dose_times)

    idx = 1
    for step in regimen.steps
        for _ in 1:step.duration_days
            for _ in 1:doses_per_day
                push!(doses, step.dose)
                idx += 1
            end
        end
    end

    # Apply loading dose
    if regimen.loading_dose !== nothing && length(doses) > 0
        doses[1] = regimen.loading_dose
    end

    # Apply compliance (same logic as above)
    if compliance !== nothing
        for i in eachindex(doses)
            if compliance.pattern == :random
                if rand(rng) > compliance.mean_compliance
                    doses[i] = 0.0
                end
            end
        end
    end

    return doses
end


"""
    total_regimen_duration(regimen)

Calculate total duration of a dosing regimen in days.

# Arguments
- `regimen`: Dosing regimen (DosingRegimen or TitrationRegimen)

# Returns
- `Int`: Total duration in days
"""
function total_regimen_duration(regimen::DosingRegimen)
    return regimen.duration_days
end

function total_regimen_duration(regimen::TitrationRegimen)
    return sum(step.duration_days for step in regimen.steps)
end


"""
    get_dose_at_time(regimen, time_hours)

Get the dose amount at a specific time.

# Arguments
- `regimen`: Dosing regimen
- `time_hours::Float64`: Time in hours from start

# Returns
- `Float64`: Dose amount (0 if no dose at this time)
"""
function get_dose_at_time(regimen::DosingRegimen, time_hours::Float64;
                          tolerance::Float64 = 0.1)
    times = dose_event_times(regimen)

    for (i, t) in enumerate(times)
        if abs(t - time_hours) < tolerance
            if i == 1 && regimen.loading_dose !== nothing
                return regimen.loading_dose
            end
            return regimen.dose_amount
        end
    end

    return 0.0
end
