# Infusion handling for zero-order drug input
# Provides infrastructure for computing time-varying infusion rates

export InfusionSchedule, compute_infusion_rate_at_time, build_infusion_schedule
export normalize_doses_with_infusion, get_infusion_events, has_infusions
export make_infusion_rate_function, get_infusion_tstops, build_bolus_callback

"""
InfusionSchedule represents the complete infusion profile for a simulation.

Contains preprocessed data for efficient infusion rate lookup during ODE solving.
The schedule stores start times, end times, and rates for all infusion doses.
"""
struct InfusionSchedule
    # Sorted vectors of infusion events
    start_times::Vector{Float64}
    end_times::Vector{Float64}
    rates::Vector{Float64}

    # For bolus doses (handled separately via callbacks)
    bolus_times::Vector{Float64}
    bolus_amounts::Vector{Float64}

    # Initial amount at t0 (sum of bolus doses at t0)
    initial_amount::Float64
end

"""
Build an InfusionSchedule from a vector of DoseEvents.

Separates bolus doses (duration=0) from infusions (duration>0),
and preprocesses infusion data for efficient lookup.

Returns:
- InfusionSchedule with all dose information organized for simulation
"""
function build_infusion_schedule(
    doses::Vector{DoseEvent},
    t0::Float64,
    t1::Float64
)::InfusionSchedule
    # Separate bolus and infusion doses
    bolus_times = Float64[]
    bolus_amounts = Float64[]
    infusion_starts = Float64[]
    infusion_ends = Float64[]
    infusion_rates = Float64[]

    initial_amount = 0.0

    for dose in doses
        if is_bolus(dose)
            # Bolus dose handling
            if dose.time == t0
                initial_amount += dose.amount
            elseif dose.time > t0 && dose.time <= t1
                push!(bolus_times, dose.time)
                push!(bolus_amounts, dose.amount)
            end
        else
            # Infusion dose handling
            start_t = max(dose.time, t0)
            end_t = min(dose_end_time(dose), t1)

            if start_t < end_t
                # Infusion overlaps with simulation window
                rate = infusion_rate(dose)

                # Handle partial infusions at t0
                if dose.time < t0
                    # Infusion started before t0, compute amount already delivered
                    # Don't add to initial_amount since infusion continues
                end

                push!(infusion_starts, start_t)
                push!(infusion_ends, end_t)
                push!(infusion_rates, rate)
            elseif dose.time < t0 && dose_end_time(dose) <= t0
                # Infusion completed before t0, add full amount as initial
                initial_amount += dose.amount
            end
        end
    end

    # Aggregate multiple bolus doses at same time
    if !isempty(bolus_times)
        unique_times = unique(bolus_times)
        aggregated_amounts = [sum(bolus_amounts[bolus_times .== t]) for t in unique_times]
        bolus_times = unique_times
        bolus_amounts = aggregated_amounts
    end

    # Sort infusions by start time
    if !isempty(infusion_starts)
        perm = sortperm(infusion_starts)
        infusion_starts = infusion_starts[perm]
        infusion_ends = infusion_ends[perm]
        infusion_rates = infusion_rates[perm]
    end

    return InfusionSchedule(
        infusion_starts,
        infusion_ends,
        infusion_rates,
        bolus_times,
        bolus_amounts,
        initial_amount
    )
end

"""
Compute the total infusion rate at time t from all active infusions.

Returns the sum of all infusion rates that are active at time t.
An infusion is active if start_time <= t < end_time.

Note: Accepts any type to support ForwardDiff autodiff (Dual numbers).
Comparisons with Dual numbers use the primal value automatically.
"""
function compute_infusion_rate_at_time(schedule::InfusionSchedule, t)::Float64
    total_rate = 0.0

    for i in eachindex(schedule.start_times)
        start_t = schedule.start_times[i]
        end_t = schedule.end_times[i]

        # Direct comparison works with Dual numbers (compares primal values)
        if start_t <= t && t < end_t
            total_rate += schedule.rates[i]
        end
    end

    return total_rate
end

"""
Check if any infusions are present in the dose schedule.
"""
has_infusions(schedule::InfusionSchedule) = !isempty(schedule.start_times)

"""
Get all event times (starts and ends) for infusions.
Used for ensuring ODE solver steps at infusion boundaries.
"""
function get_infusion_events(schedule::InfusionSchedule)::Vector{Float64}
    events = Float64[]
    append!(events, schedule.start_times)
    append!(events, schedule.end_times)
    return sort(unique(events))
end

"""
Normalize doses for simulation, handling both bolus and infusion doses.

This is an extended version of normalize_doses_for_sim that properly
handles the duration field of DoseEvent.

Returns:
- InfusionSchedule containing all preprocessed dose information
"""
function normalize_doses_with_infusion(
    doses::Vector{DoseEvent},
    t0::Float64,
    t1::Float64
)::InfusionSchedule
    return build_infusion_schedule(doses, t0, t1)
end

"""
Create a closure for computing infusion rate that can be passed to ODE.

Returns a function R(t) that gives the total infusion rate at time t.
"""
function make_infusion_rate_function(schedule::InfusionSchedule)
    return t -> compute_infusion_rate_at_time(schedule, t)
end

"""
Build callbacks for infusion events (start/stop) to ensure solver accuracy.

Uses tstops to force the solver to step exactly at infusion boundaries.
This ensures accurate handling of rate discontinuities.
"""
function get_infusion_tstops(schedule::InfusionSchedule, t0::Float64, t1::Float64)::Vector{Float64}
    tstops = Float64[]

    for i in eachindex(schedule.start_times)
        start_t = schedule.start_times[i]
        end_t = schedule.end_times[i]

        if start_t > t0 && start_t < t1
            push!(tstops, start_t)
        end
        if end_t > t0 && end_t < t1
            push!(tstops, end_t)
        end
    end

    return sort(unique(tstops))
end

"""
Create a PresetTimeCallback for bolus doses only.

This is used in combination with infusion rate in the ODE.
Bolus doses are still handled via callbacks (instantaneous addition),
while infusions are handled via the rate term in the ODE.
"""
function build_bolus_callback(
    schedule::InfusionSchedule,
    target_index::Int
)
    if isempty(schedule.bolus_times)
        return nothing
    end

    bolus_times = schedule.bolus_times
    bolus_amounts = schedule.bolus_amounts

    function affect!(integrator)
        idx = findfirst(==(integrator.t), bolus_times)
        if idx !== nothing
            integrator.u[target_index] += bolus_amounts[idx]
        end
    end

    return PresetTimeCallback(bolus_times, affect!)
end
