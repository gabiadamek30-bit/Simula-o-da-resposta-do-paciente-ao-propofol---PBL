# CDISC to NeoPKPD Converter
# Convert CDISC datasets to NeoPKPD format for simulation and analysis

export cdisc_to_observed, CDISCConversionConfig, CDISCConversionResult

"""
Configuration for CDISC to NeoPKPD conversion.

Fields:
- time_reference: How to calculate time (`:first_dose`, `:reference_time`, `:study_day`)
- time_units: Output time units (`:hours`, `:days`, `:minutes`)
- include_blq: Whether to include below LOQ observations
- blq_handling: How to handle BLQ (`:as_is`, `:half_lloq`, `:zero`, `:exclude`)
- analyte_filter: Filter to specific analyte/test code (empty = all)
- covariate_mapping: Mapping of DM fields to covariate names
"""
struct CDISCConversionConfig
    time_reference::Symbol
    time_units::Symbol
    include_blq::Bool
    blq_handling::Symbol
    analyte_filter::String
    covariate_mapping::Dict{String,Symbol}

    function CDISCConversionConfig(;
        time_reference::Symbol=:first_dose,
        time_units::Symbol=:hours,
        include_blq::Bool=true,
        blq_handling::Symbol=:half_lloq,
        analyte_filter::String="",
        covariate_mapping::Dict{String,Symbol}=Dict{String,Symbol}(
            "AGE" => :age,
            "SEX" => :sex,
            "RACE" => :race
        )
    )
        @assert time_reference in [:first_dose, :reference_time, :study_day]
        @assert time_units in [:hours, :days, :minutes]
        @assert blq_handling in [:as_is, :half_lloq, :zero, :exclude]
        new(time_reference, time_units, include_blq, blq_handling, analyte_filter, covariate_mapping)
    end
end

"""
Result of CDISC to NeoPKPD conversion.

Fields:
- observed_data: The converted ObservedData
- warnings: Any warnings generated during conversion
- errors: Any errors encountered
- n_subjects: Number of subjects converted
- n_observations: Total observations converted
- n_excluded: Number of observations excluded
"""
struct CDISCConversionResult
    observed_data::Union{Nothing,ObservedData}
    warnings::Vector{String}
    errors::Vector{String}
    n_subjects::Int
    n_observations::Int
    n_excluded::Int
end

"""
Convert CDISC dataset to NeoPKPD ObservedData format.

Arguments:
- dataset: CDISCDataset to convert
- config: Conversion configuration (optional)

Returns:
- CDISCConversionResult
"""
function cdisc_to_observed(
    dataset::CDISCDataset;
    config::CDISCConversionConfig=CDISCConversionConfig()
)::CDISCConversionResult
    warnings = String[]
    errors = String[]
    n_excluded = 0

    if isempty(dataset.pc)
        push!(errors, "No PC (concentration) records found in dataset")
        return CDISCConversionResult(nothing, warnings, errors, 0, 0, 0)
    end

    # Group PC records by subject
    pc_by_subject = Dict{String,Vector{PCRecord}}()
    for rec in dataset.pc
        # Filter by analyte if specified
        if !isempty(config.analyte_filter) && rec.pctestcd != config.analyte_filter
            continue
        end

        if !haskey(pc_by_subject, rec.usubjid)
            pc_by_subject[rec.usubjid] = PCRecord[]
        end
        push!(pc_by_subject[rec.usubjid], rec)
    end

    # Group EX records by subject
    ex_by_subject = Dict{String,Vector{EXRecord}}()
    for rec in dataset.ex
        if !haskey(ex_by_subject, rec.usubjid)
            ex_by_subject[rec.usubjid] = EXRecord[]
        end
        push!(ex_by_subject[rec.usubjid], rec)
    end

    # Create DM lookup
    dm_by_subject = Dict{String,DMRecord}()
    for rec in dataset.dm
        dm_by_subject[rec.usubjid] = rec
    end

    # Convert each subject
    subjects = SubjectData[]
    analyte = ""
    units = ""

    for (usubjid, pc_records) in pc_by_subject
        # Sort by time
        sort!(pc_records, by=r -> something(r.pctptnum, 0.0))

        # Get reference time
        ref_time = _get_reference_time(pc_records, ex_by_subject, dm_by_subject, usubjid, config)

        # Extract times and observations
        times = Float64[]
        observations = Float64[]
        blq_flags = Bool[]
        lloq = 0.0

        for rec in pc_records
            # Get time relative to reference
            t = _get_observation_time(rec, ref_time, config)
            if t === nothing
                n_excluded += 1
                continue
            end

            # Get observation value
            conc = rec.pcstresn
            is_blq = rec.pcstat == "BLQ" || rec.pcstat == "<LLOQ" || conc === nothing

            if is_blq
                if !config.include_blq && config.blq_handling == :exclude
                    n_excluded += 1
                    continue
                end

                # Handle BLQ based on config
                if rec.pclloq !== nothing
                    lloq = rec.pclloq
                end

                if config.blq_handling == :half_lloq && lloq > 0
                    conc = lloq / 2.0
                elseif config.blq_handling == :zero
                    conc = 0.0
                elseif config.blq_handling == :as_is && conc === nothing
                    conc = 0.0
                end
            end

            if conc === nothing
                n_excluded += 1
                continue
            end

            push!(times, t)
            push!(observations, conc)
            push!(blq_flags, is_blq)

            # Capture analyte and units from first record
            if isempty(analyte)
                analyte = rec.pctestcd
            end
            if isempty(units)
                units = rec.pcstresu
            end
        end

        if isempty(times)
            push!(warnings, "Subject $usubjid has no valid observations after filtering")
            continue
        end

        # Convert EX records to DoseEvents
        doses = _extract_doses(ex_by_subject, usubjid, ref_time, config)

        if isempty(doses)
            push!(warnings, "Subject $usubjid has no dosing records")
        end

        # Extract covariates from DM
        covariates = _extract_covariates(dm_by_subject, usubjid, config)

        # Create SubjectData
        subject = SubjectData(
            usubjid,
            times,
            observations,
            doses;
            covariates=covariates,
            lloq=lloq,
            blq_flags=blq_flags
        )
        push!(subjects, subject)
    end

    if isempty(subjects)
        push!(errors, "No subjects with valid data after conversion")
        return CDISCConversionResult(nothing, warnings, errors, 0, 0, n_excluded)
    end

    # Create ObservedData
    time_units_str = config.time_units == :hours ? "h" : config.time_units == :days ? "d" : "min"
    observed = ObservedData(
        subjects;
        study_id=dataset.study_id,
        analyte=analyte,
        units=units,
        time_units=time_units_str
    )

    n_obs = sum(length(s.observations) for s in subjects)

    return CDISCConversionResult(observed, warnings, errors, length(subjects), n_obs, n_excluded)
end

"""
Get reference time for a subject based on config.
"""
function _get_reference_time(
    pc_records::Vector{PCRecord},
    ex_by_subject::Dict{String,Vector{EXRecord}},
    dm_by_subject::Dict{String,DMRecord},
    usubjid::String,
    config::CDISCConversionConfig
)::Float64
    if config.time_reference == :first_dose
        # Use first dose time as reference
        if haskey(ex_by_subject, usubjid) && !isempty(ex_by_subject[usubjid])
            ex_records = ex_by_subject[usubjid]
            sort!(ex_records, by=r -> r.exstdtc)
            first_dose_time = parse_iso8601_datetime(ex_records[1].exstdtc)
            if first_dose_time !== nothing
                return first_dose_time
            end
        end
        # Fall back to first PC record reference time
        for rec in pc_records
            ref_time = parse_iso8601_datetime(rec.pcrftdtc)
            if ref_time !== nothing
                return ref_time
            end
        end
    elseif config.time_reference == :reference_time
        # Use RFSTDTC from DM
        if haskey(dm_by_subject, usubjid)
            ref_time = parse_iso8601_datetime(dm_by_subject[usubjid].rfstdtc)
            if ref_time !== nothing
                return ref_time
            end
        end
    end

    # Default: use study day 1 as reference
    return 0.0
end

"""
Get observation time relative to reference.
"""
function _get_observation_time(
    rec::PCRecord,
    ref_time::Float64,
    config::CDISCConversionConfig
)::Union{Nothing,Float64}
    # Try PCTPTNUM first (planned time point number in hours)
    if rec.pctptnum !== nothing
        t = rec.pctptnum
        # Convert to target units
        if config.time_units == :days
            t = t / 24.0
        elseif config.time_units == :minutes
            t = t * 60.0
        end
        return t
    end

    # Try PCELTM (elapsed time from reference)
    if !isempty(rec.pceltm)
        t = parse_elapsed_time(rec.pceltm)
        if config.time_units == :days
            t = t / 24.0
        elseif config.time_units == :minutes
            t = t * 60.0
        end
        return t
    end

    # Try PCSTDTC - ref_time
    obs_time = parse_iso8601_datetime(rec.pcstdtc)
    if obs_time !== nothing
        t = obs_time - ref_time  # Hours
        if config.time_units == :days
            t = t / 24.0
        elseif config.time_units == :minutes
            t = t * 60.0
        end
        return t
    end

    # Try PCDY (study day)
    if rec.pcdy !== nothing
        t = Float64(rec.pcdy - 1) * 24.0  # Day 1 = hour 0
        if config.time_units == :days
            t = t / 24.0
        elseif config.time_units == :minutes
            t = t * 60.0
        end
        return t
    end

    return nothing
end

"""
Extract DoseEvents from EX records for a subject.
"""
function _extract_doses(
    ex_by_subject::Dict{String,Vector{EXRecord}},
    usubjid::String,
    ref_time::Float64,
    config::CDISCConversionConfig
)::Vector{DoseEvent}
    doses = DoseEvent[]

    if !haskey(ex_by_subject, usubjid)
        return doses
    end

    for rec in ex_by_subject[usubjid]
        # Get dose time relative to reference
        dose_time = parse_iso8601_datetime(rec.exstdtc)
        if dose_time === nothing && rec.exdy !== nothing
            dose_time = Float64(rec.exdy - 1) * 24.0
        end

        if dose_time === nothing
            continue
        end

        t = dose_time - ref_time  # Hours
        if config.time_units == :days
            t = t / 24.0
        elseif config.time_units == :minutes
            t = t * 60.0
        end

        # Get dose amount
        amount = rec.exdose

        # Get duration for infusions
        duration = 0.0
        if !isempty(rec.exdur)
            duration = parse_iso8601_duration(rec.exdur)
            if config.time_units == :days
                duration = duration / 24.0
            elseif config.time_units == :minutes
                duration = duration * 60.0
            end
        end

        push!(doses, DoseEvent(t, amount, duration))
    end

    # Sort by time
    sort!(doses, by=d -> d.time)

    return doses
end

"""
Extract covariates from DM record.
"""
function _extract_covariates(
    dm_by_subject::Dict{String,DMRecord},
    usubjid::String,
    config::CDISCConversionConfig
)::Dict{Symbol,Any}
    covariates = Dict{Symbol,Any}()

    if !haskey(dm_by_subject, usubjid)
        return covariates
    end

    dm = dm_by_subject[usubjid]

    # Map standard covariates
    if dm.age !== nothing
        covariates[:age] = dm.age
    end

    if !isempty(dm.sex)
        covariates[:sex] = dm.sex
    end

    if !isempty(dm.race)
        covariates[:race] = dm.race
    end

    if !isempty(dm.ethnic)
        covariates[:ethnic] = dm.ethnic
    end

    if !isempty(dm.armcd)
        covariates[:arm] = dm.armcd
    end

    if !isempty(dm.country)
        covariates[:country] = dm.country
    end

    return covariates
end
