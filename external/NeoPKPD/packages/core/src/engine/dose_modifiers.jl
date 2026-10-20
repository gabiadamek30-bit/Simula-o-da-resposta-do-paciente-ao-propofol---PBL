# Dose Modifiers: ALAG (Absorption Lag Time) and Bioavailability (F)
# Industry-standard implementation following NONMEM conventions
#
# ALAG: Absorption lag time - delays the onset of absorption by a specified time
#   - ALAG1 delays doses entering compartment 1 (typically depot/gut)
#   - ALAG2 delays doses entering compartment 2, etc.
#   - Implemented by shifting the effective dose time: t_effective = t_nominal + ALAG
#
# F: Bioavailability fraction - fraction of dose that reaches systemic circulation
#   - F1 applies to doses entering compartment 1
#   - F2 applies to doses entering compartment 2, etc.
#   - Implemented by multiplying dose amount: amount_effective = amount_nominal * F
#
# Usage patterns:
#   1. Per-dose modifiers: Specify ALAG/F for each DoseEvent
#   2. Model-level modifiers: Default ALAG/F applied to all doses targeting a compartment
#   3. Population simulation: ALAG/F can have IIV (inter-individual variability)

export DoseModifiers, apply_dose_modifiers, validate_dose_modifiers
export ModifiedDoseEvent, apply_alag_to_dose, apply_bioavailability_to_dose
export DoseModifierSpec, CompartmentDoseModifiers
export process_doses_with_modifiers

"""
    DoseModifiers

Per-compartment dose modification parameters.

# Fields
- `alag::Float64`: Absorption lag time (time units, >= 0)
- `bioavailability::Float64`: Bioavailability fraction (0 < F <= 1, or > 1 for supersaturation)

# NONMEM Equivalence
- `alag` corresponds to ALAG1, ALAG2, etc.
- `bioavailability` corresponds to F1, F2, etc.

# Examples
```julia
# Standard oral absorption with 0.5 hour lag and 80% bioavailability
modifiers = DoseModifiers(alag=0.5, bioavailability=0.8)

# IV bolus with no lag and 100% bioavailability (default)
modifiers = DoseModifiers()
```
"""
struct DoseModifiers
    alag::Float64           # Absorption lag time (time units)
    bioavailability::Float64  # Bioavailability fraction (typically 0-1)

    function DoseModifiers(; alag::Float64=0.0, bioavailability::Float64=1.0)
        new(alag, bioavailability)
    end
end

# Convenience constructor with positional arguments
DoseModifiers(alag::Float64, bioavailability::Float64) = DoseModifiers(alag=alag, bioavailability=bioavailability)

# Default modifiers (no lag, 100% bioavailability)
const DEFAULT_DOSE_MODIFIERS = DoseModifiers()

"""
    validate_dose_modifiers(modifiers::DoseModifiers)

Validate dose modifiers for physiological plausibility.

# Validation Rules
- ALAG must be non-negative
- Bioavailability must be positive (can exceed 1.0 in special cases like supersaturation)
- Warning if bioavailability > 1.0 (unusual but valid in some scenarios)
"""
function validate_dose_modifiers(modifiers::DoseModifiers)
    if modifiers.alag < 0.0
        error("ALAG must be non-negative, got $(modifiers.alag)")
    end

    if modifiers.bioavailability <= 0.0
        error("Bioavailability must be positive, got $(modifiers.bioavailability)")
    end

    if modifiers.bioavailability > 1.5
        @warn "Bioavailability $(modifiers.bioavailability) exceeds 150%. This is unusual - verify this is intended."
    end

    return nothing
end

"""
    CompartmentDoseModifiers

Collection of dose modifiers indexed by target compartment.

# Fields
- `modifiers::Dict{Int, DoseModifiers}`: Compartment index -> modifiers

# NONMEM Equivalence
```
; In NONMEM \$PK block:
ALAG1 = THETA(5)
F1 = THETA(6)
ALAG2 = THETA(7)
F2 = 1.0

; Equivalent in NeoPKPD:
CompartmentDoseModifiers(Dict(
    1 => DoseModifiers(alag=theta5, bioavailability=theta6),
    2 => DoseModifiers(alag=theta7, bioavailability=1.0)
))
```
"""
struct CompartmentDoseModifiers
    modifiers::Dict{Int, DoseModifiers}

    function CompartmentDoseModifiers(modifiers::Dict{Int, DoseModifiers})
        # Validate all modifiers
        for (comp, mod) in modifiers
            if comp < 1
                error("Compartment index must be >= 1, got $comp")
            end
            validate_dose_modifiers(mod)
        end
        new(modifiers)
    end
end

# Convenience constructors - empty dict
function CompartmentDoseModifiers()
    return CompartmentDoseModifiers(Dict{Int, DoseModifiers}())
end

# Convenience constructor from pairs
function CompartmentDoseModifiers(pairs::Pair{Int, DoseModifiers}...)
    return CompartmentDoseModifiers(Dict(pairs...))
end

"""
Get dose modifiers for a specific compartment, returning defaults if not specified.
"""
function get_modifiers_for_compartment(cdm::CompartmentDoseModifiers, compartment::Int)::DoseModifiers
    return get(cdm.modifiers, compartment, DEFAULT_DOSE_MODIFIERS)
end

export get_modifiers_for_compartment

"""
    ModifiedDoseEvent

A dose event with applied modifiers (ALAG and F).

# Fields
- `original::DoseEvent`: The original unmodified dose
- `effective_time::Float64`: Time after applying ALAG (original.time + alag)
- `effective_amount::Float64`: Amount after applying F (original.amount * bioavailability)
- `effective_duration::Float64`: Duration (unchanged, F only affects amount)
- `target_compartment::Int`: Which compartment receives this dose
- `modifiers::DoseModifiers`: The modifiers that were applied

This struct is used internally during simulation to track how doses are modified.
"""
struct ModifiedDoseEvent
    original::DoseEvent
    effective_time::Float64
    effective_amount::Float64
    effective_duration::Float64
    target_compartment::Int
    modifiers::DoseModifiers
end

"""
Check if modified dose event is a bolus.
"""
is_bolus(dose::ModifiedDoseEvent) = dose.effective_duration == 0.0

"""
Check if modified dose event is an infusion.
"""
is_infusion(dose::ModifiedDoseEvent) = dose.effective_duration > 0.0

"""
Get infusion rate for modified dose event.
"""
function infusion_rate(dose::ModifiedDoseEvent)
    if is_bolus(dose)
        return Inf
    else
        return dose.effective_amount / dose.effective_duration
    end
end

"""
Get end time of modified dose event.
"""
dose_end_time(dose::ModifiedDoseEvent) = dose.effective_time + dose.effective_duration

"""
    apply_alag_to_dose(dose::DoseEvent, alag::Float64)

Apply absorption lag time to a dose event.

# Arguments
- `dose`: Original dose event
- `alag`: Absorption lag time (time units)

# Returns
- New dose time = original time + alag
"""
function apply_alag_to_dose(dose::DoseEvent, alag::Float64)::Float64
    return dose.time + alag
end

"""
    apply_bioavailability_to_dose(dose::DoseEvent, bioavailability::Float64)

Apply bioavailability fraction to a dose event.

# Arguments
- `dose`: Original dose event
- `bioavailability`: Bioavailability fraction (F)

# Returns
- Effective amount = original amount * F
"""
function apply_bioavailability_to_dose(dose::DoseEvent, bioavailability::Float64)::Float64
    return dose.amount * bioavailability
end

"""
    apply_dose_modifiers(dose::DoseEvent, modifiers::DoseModifiers, target_compartment::Int)

Apply both ALAG and F to a dose event, returning a ModifiedDoseEvent.

# Arguments
- `dose`: Original dose event
- `modifiers`: ALAG and F values to apply
- `target_compartment`: Which compartment this dose enters

# Returns
- `ModifiedDoseEvent` with effective time and amount
"""
function apply_dose_modifiers(
    dose::DoseEvent,
    modifiers::DoseModifiers,
    target_compartment::Int
)::ModifiedDoseEvent
    effective_time = apply_alag_to_dose(dose, modifiers.alag)
    effective_amount = apply_bioavailability_to_dose(dose, modifiers.bioavailability)

    return ModifiedDoseEvent(
        dose,
        effective_time,
        effective_amount,
        dose.duration,  # Duration is not affected by modifiers
        target_compartment,
        modifiers
    )
end

"""
    process_doses_with_modifiers(doses, compartment_modifiers, default_compartment, t0, t1)

Process a vector of dose events, applying ALAG and bioavailability modifiers.

# Arguments
- `doses::Vector{DoseEvent}`: Original dose events
- `compartment_modifiers::CompartmentDoseModifiers`: Per-compartment modifiers
- `default_compartment::Int`: Default target compartment if not specified per-dose
- `t0::Float64`: Simulation start time
- `t1::Float64`: Simulation end time

# Returns
- `modified_doses::Vector{ModifiedDoseEvent}`: Doses with applied modifiers, filtered to [t0, t1]
- `initial_amount::Float64`: Amount to add to initial condition (doses with effective_time == t0)
- `callback_doses::Vector{ModifiedDoseEvent}`: Doses requiring callback (effective_time > t0)
"""
function process_doses_with_modifiers(
    doses::Vector{DoseEvent},
    compartment_modifiers::CompartmentDoseModifiers,
    default_compartment::Int,
    t0::Float64,
    t1::Float64
)
    modified_doses = ModifiedDoseEvent[]
    initial_amount = 0.0
    callback_doses = ModifiedDoseEvent[]

    for dose in doses
        # Get modifiers for target compartment
        modifiers = get_modifiers_for_compartment(compartment_modifiers, default_compartment)

        # Apply modifiers
        modified = apply_dose_modifiers(dose, modifiers, default_compartment)

        # Check if modified dose falls within simulation window
        if modified.effective_time < t0
            # Dose time is before simulation start
            if is_bolus(modified)
                # Bolus before t0 - add to initial amount
                initial_amount += modified.effective_amount
            else
                # Infusion that started before t0
                end_time = dose_end_time(modified)
                if end_time <= t0
                    # Infusion completed before t0 - add full amount to initial
                    initial_amount += modified.effective_amount
                else
                    # Infusion overlaps t0 - will be handled by infusion schedule
                    push!(modified_doses, modified)
                end
            end
        elseif modified.effective_time == t0
            # Dose at exactly t0
            if is_bolus(modified)
                initial_amount += modified.effective_amount
            else
                push!(modified_doses, modified)
            end
        elseif modified.effective_time <= t1
            # Dose within (t0, t1]
            push!(callback_doses, modified)
            push!(modified_doses, modified)
        end
        # Doses after t1 are ignored
    end

    return modified_doses, initial_amount, callback_doses
end

"""
    DoseModifierSpec

Specification for dose modifiers to be used with a model.

This is the high-level specification that users provide, which gets converted
to CompartmentDoseModifiers during simulation setup.

# Fields
- `compartment_modifiers::Dict{Int, DoseModifiers}`: Per-compartment ALAG/F
- `default_alag::Float64`: Default ALAG for compartments not in dict (typically 0.0)
- `default_bioavailability::Float64`: Default F for compartments not in dict (typically 1.0)

# Example
```julia
# Model with absorption lag for oral dosing (compartment 1 = gut)
dose_spec = DoseModifierSpec(
    compartment_modifiers = Dict(
        1 => DoseModifiers(alag=0.5, bioavailability=0.8)
    )
)
```
"""
struct DoseModifierSpec
    compartment_modifiers::Dict{Int, DoseModifiers}
    default_alag::Float64
    default_bioavailability::Float64

    function DoseModifierSpec(;
        compartment_modifiers::Dict{Int, DoseModifiers} = Dict{Int, DoseModifiers}(),
        default_alag::Float64 = 0.0,
        default_bioavailability::Float64 = 1.0
    )
        # Validate
        if default_alag < 0.0
            error("default_alag must be non-negative")
        end
        if default_bioavailability <= 0.0
            error("default_bioavailability must be positive")
        end

        for (comp, mod) in compartment_modifiers
            validate_dose_modifiers(mod)
        end

        new(compartment_modifiers, default_alag, default_bioavailability)
    end
end

"""
Convert DoseModifierSpec to CompartmentDoseModifiers.
"""
function to_compartment_modifiers(spec::DoseModifierSpec)::CompartmentDoseModifiers
    return CompartmentDoseModifiers(spec.compartment_modifiers)
end

export DoseModifierSpec, to_compartment_modifiers

# ============================================================================
# Modified Infusion Schedule with ALAG/F Support
# ============================================================================

"""
    ModifiedInfusionSchedule

Extension of InfusionSchedule that incorporates ALAG and bioavailability.

Same structure as InfusionSchedule but created from ModifiedDoseEvents.
"""
struct ModifiedInfusionSchedule
    # Sorted vectors of infusion events
    start_times::Vector{Float64}
    end_times::Vector{Float64}
    rates::Vector{Float64}
    target_compartments::Vector{Int}

    # For bolus doses (handled via callbacks)
    bolus_times::Vector{Float64}
    bolus_amounts::Vector{Float64}
    bolus_compartments::Vector{Int}

    # Initial amount at t0 (sum of modified bolus doses at t0)
    initial_amount::Float64
end

export ModifiedInfusionSchedule

"""
    build_modified_infusion_schedule(modified_doses, initial_amount, t0, t1)

Build a ModifiedInfusionSchedule from pre-processed modified doses.
"""
function build_modified_infusion_schedule(
    modified_doses::Vector{ModifiedDoseEvent},
    initial_amount::Float64,
    t0::Float64,
    t1::Float64
)::ModifiedInfusionSchedule
    bolus_times = Float64[]
    bolus_amounts = Float64[]
    bolus_compartments = Int[]

    infusion_starts = Float64[]
    infusion_ends = Float64[]
    infusion_rates = Float64[]
    infusion_compartments = Int[]

    for dose in modified_doses
        if is_bolus(dose)
            if dose.effective_time > t0 && dose.effective_time <= t1
                push!(bolus_times, dose.effective_time)
                push!(bolus_amounts, dose.effective_amount)
                push!(bolus_compartments, dose.target_compartment)
            end
        else
            # Infusion
            start_t = max(dose.effective_time, t0)
            end_t = min(dose_end_time(dose), t1)

            if start_t < end_t
                push!(infusion_starts, start_t)
                push!(infusion_ends, end_t)
                push!(infusion_rates, infusion_rate(dose))
                push!(infusion_compartments, dose.target_compartment)
            end
        end
    end

    # Aggregate bolus doses at same time and compartment
    if !isempty(bolus_times)
        aggregated = Dict{Tuple{Float64, Int}, Float64}()
        for (t, amt, comp) in zip(bolus_times, bolus_amounts, bolus_compartments)
            key = (t, comp)
            aggregated[key] = get(aggregated, key, 0.0) + amt
        end

        bolus_times = Float64[]
        bolus_amounts = Float64[]
        bolus_compartments = Int[]

        for ((t, comp), amt) in sort(collect(aggregated), by=x->x[1][1])
            push!(bolus_times, t)
            push!(bolus_amounts, amt)
            push!(bolus_compartments, comp)
        end
    end

    # Sort infusions by start time
    if !isempty(infusion_starts)
        perm = sortperm(infusion_starts)
        infusion_starts = infusion_starts[perm]
        infusion_ends = infusion_ends[perm]
        infusion_rates = infusion_rates[perm]
        infusion_compartments = infusion_compartments[perm]
    end

    return ModifiedInfusionSchedule(
        infusion_starts,
        infusion_ends,
        infusion_rates,
        infusion_compartments,
        bolus_times,
        bolus_amounts,
        bolus_compartments,
        initial_amount
    )
end

export build_modified_infusion_schedule

"""
Compute total infusion rate at time t for a specific compartment.
"""
function compute_modified_infusion_rate_at_time(
    schedule::ModifiedInfusionSchedule,
    t,
    compartment::Int
)::Float64
    total_rate = 0.0

    for i in eachindex(schedule.start_times)
        if schedule.target_compartments[i] == compartment
            start_t = schedule.start_times[i]
            end_t = schedule.end_times[i]

            if start_t <= t && t < end_t
                total_rate += schedule.rates[i]
            end
        end
    end

    return total_rate
end

export compute_modified_infusion_rate_at_time

"""
Check if schedule has any infusions.
"""
has_infusions(schedule::ModifiedInfusionSchedule) = !isempty(schedule.start_times)

"""
Get all time stops for solver accuracy (infusion starts/ends).
"""
function get_modified_infusion_tstops(
    schedule::ModifiedInfusionSchedule,
    t0::Float64,
    t1::Float64
)::Vector{Float64}
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

    # Also add bolus times for ALAG-delayed doses
    for t in schedule.bolus_times
        if t > t0 && t < t1
            push!(tstops, t)
        end
    end

    return sort(unique(tstops))
end

export get_modified_infusion_tstops

"""
Build callback for bolus doses with multi-compartment support.
"""
function build_modified_bolus_callback(schedule::ModifiedInfusionSchedule)
    if isempty(schedule.bolus_times)
        return nothing
    end

    bolus_times = schedule.bolus_times
    bolus_amounts = schedule.bolus_amounts
    bolus_compartments = schedule.bolus_compartments

    function affect!(integrator)
        # Find all doses at this time
        for i in eachindex(bolus_times)
            if bolus_times[i] == integrator.t
                comp = bolus_compartments[i]
                integrator.u[comp] += bolus_amounts[i]
            end
        end
    end

    unique_times = unique(bolus_times)
    return PresetTimeCallback(unique_times, affect!)
end

export build_modified_bolus_callback

# ============================================================================
# Normalized Dose Processing with ALAG/F
# ============================================================================

"""
    normalize_doses_with_modifiers(doses, modifiers, target_compartment, t0, t1)

Normalize doses for simulation, applying ALAG and bioavailability.

This is the main entry point for processing doses with modifiers.

# Arguments
- `doses::Vector{DoseEvent}`: Original dose events
- `modifiers::DoseModifiers`: ALAG and F to apply
- `target_compartment::Int`: Which compartment receives these doses
- `t0::Float64`: Simulation start time
- `t1::Float64`: Simulation end time

# Returns
- `initial_amount::Float64`: Amount for initial condition at t0
- `dose_times::Vector{Float64}`: Times of callback doses (after ALAG)
- `dose_amounts::Vector{Float64}`: Amounts of callback doses (after F)
"""
function normalize_doses_with_modifiers(
    doses::Vector{DoseEvent},
    modifiers::DoseModifiers,
    target_compartment::Int,
    t0::Float64,
    t1::Float64
)
    validate_dose_modifiers(modifiers)

    initial_amount = 0.0
    dose_times_dict = Dict{Float64, Float64}()  # time -> summed amount

    for dose in doses
        # Apply modifiers
        effective_time = apply_alag_to_dose(dose, modifiers.alag)
        effective_amount = apply_bioavailability_to_dose(dose, modifiers.bioavailability)

        if is_bolus(dose)
            if effective_time < t0
                # Dose completed before t0, add to initial
                initial_amount += effective_amount
            elseif effective_time == t0
                # Dose at exactly t0
                initial_amount += effective_amount
            elseif effective_time <= t1
                # Dose within (t0, t1]
                dose_times_dict[effective_time] = get(dose_times_dict, effective_time, 0.0) + effective_amount
            end
            # Doses after t1 are ignored
        else
            # Infusion - handle more complex timing
            # For now, treat similarly to bolus but consider duration
            effective_end = effective_time + dose.duration

            if effective_end <= t0
                # Infusion completed before t0
                initial_amount += effective_amount
            elseif effective_time >= t1
                # Infusion starts after simulation ends - ignore
            else
                # Infusion overlaps with simulation - add to schedule
                # Note: This is simplified; full infusion handling would need
                # the infusion schedule approach
                dose_times_dict[effective_time] = get(dose_times_dict, effective_time, 0.0) + effective_amount
            end
        end
    end

    # Convert dict to sorted vectors
    dose_times = sort(collect(keys(dose_times_dict)))
    dose_amounts = [dose_times_dict[t] for t in dose_times]

    return initial_amount, dose_times, dose_amounts
end

export normalize_doses_with_modifiers
