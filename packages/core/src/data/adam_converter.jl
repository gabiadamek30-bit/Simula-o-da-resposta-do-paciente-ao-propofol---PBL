# ADaM to NeoPKPD Format Converter
# Convert ADaM datasets to NeoPKPD analysis-ready format

export adam_to_observed_data, adam_to_subject_data
export AdamConversionOptions, AdamBLQMethod

# ============================================================================
# Conversion Options
# ============================================================================

"""
BLQ (Below Limit of Quantitation) handling methods for ADaM conversion.
Note: This is separate from the VPC BLQMethod enum used for VPC analysis.
"""
@enum AdamBLQMethod begin
    ADAM_BLQ_MISSING      # Treat as missing
    ADAM_BLQ_ZERO         # Replace with 0
    ADAM_BLQ_LLOQ_HALF    # Replace with LLOQ/2
    ADAM_BLQ_LLOQ         # Replace with LLOQ
    ADAM_BLQ_KEEP         # Keep original value (if any)
end

"""
    AdamConversionOptions

Options for converting ADaM datasets to NeoPKPD format.

Fields:
- time_variable: Which time variable to use (:afrlt, :arrlt, :areltm, :nfrlt, :nrrlt, :nreltm)
- blq_method: How to handle BLQ observations (AdamBLQMethod enum)
- include_blq: Whether to include BLQ observations in output
- analysis_flag: Which analysis flag records to include (e.g., "ANL01FL")
- filter_paramcd: Filter by PARAMCD (empty string = all)
- time_units: Output time units ("h", "d", "min")
- conc_units: Output concentration units
"""
struct AdamConversionOptions
    time_variable::Symbol
    blq_method::AdamBLQMethod
    include_blq::Bool
    analysis_flag::String
    filter_paramcd::String
    time_units::String
    conc_units::String

    function AdamConversionOptions(;
        time_variable::Symbol=:afrlt,
        blq_method::AdamBLQMethod=ADAM_BLQ_LLOQ_HALF,
        include_blq::Bool=true,
        analysis_flag::String="Y",
        filter_paramcd::String="",
        time_units::String="h",
        conc_units::String=""
    )
        new(time_variable, blq_method, include_blq, analysis_flag,
            filter_paramcd, time_units, conc_units)
    end
end

export AdamConversionOptions

# ============================================================================
# Helper Functions
# ============================================================================

"""
Get time value from ADPC record based on options.
"""
function _get_adpc_time(record::ADPCRecord, time_var::Symbol)::Union{Nothing, Float64}
    if time_var == :afrlt
        return record.afrlt
    elseif time_var == :arrlt
        return record.arrlt
    elseif time_var == :areltm
        return record.areltm
    elseif time_var == :nfrlt
        return record.nfrlt
    elseif time_var == :nrrlt
        return record.nrrlt
    elseif time_var == :nreltm
        return record.nreltm
    elseif time_var == :atptn
        return record.atptn
    else
        return record.afrlt  # Default
    end
end

"""
Handle BLQ observation value.
"""
function _handle_blq(aval::Union{Nothing, Float64}, lloq::Union{Nothing, Float64},
                      method::AdamBLQMethod)::Union{Nothing, Float64}
    if aval !== nothing
        return aval
    end

    lloq_val = something(lloq, 0.0)

    if method == ADAM_BLQ_MISSING
        return nothing
    elseif method == ADAM_BLQ_ZERO
        return 0.0
    elseif method == ADAM_BLQ_LLOQ_HALF
        return lloq_val / 2.0
    elseif method == ADAM_BLQ_LLOQ
        return lloq_val
    elseif method == ADAM_BLQ_KEEP
        return aval
    else
        return nothing
    end
end

"""
Extract covariates from ADSL record.
"""
function _extract_covariates(adsl::ADSLRecord)::Dict{Symbol, Any}
    covs = Dict{Symbol, Any}()

    # Demographics
    adsl.age !== nothing && (covs[:age] = adsl.age)
    !isempty(adsl.sex) && (covs[:sex] = adsl.sex)
    adsl.sexn !== nothing && (covs[:sexn] = adsl.sexn)
    !isempty(adsl.race) && (covs[:race] = adsl.race)
    !isempty(adsl.ethnic) && (covs[:ethnic] = adsl.ethnic)

    # Baseline covariates
    adsl.bweight !== nothing && (covs[:bweight] = adsl.bweight)
    adsl.bheight !== nothing && (covs[:bheight] = adsl.bheight)
    adsl.bbmi !== nothing && (covs[:bbmi] = adsl.bbmi)
    adsl.bcrcl !== nothing && (covs[:bcrcl] = adsl.bcrcl)
    adsl.begfr !== nothing && (covs[:begfr] = adsl.begfr)
    adsl.bast !== nothing && (covs[:bast] = adsl.bast)
    adsl.balt !== nothing && (covs[:balt] = adsl.balt)
    adsl.bbili !== nothing && (covs[:bbili] = adsl.bbili)
    adsl.balb !== nothing && (covs[:balb] = adsl.balb)

    # Treatment
    !isempty(adsl.trt01a) && (covs[:treatment] = adsl.trt01a)
    !isempty(adsl.arm) && (covs[:arm] = adsl.arm)

    return covs
end

"""
Convert ADEX records to DoseEvent vector for a subject.
"""
function _convert_doses(adex_records::Vector{ADEXRecord})::Vector{DoseEvent}
    doses = DoseEvent[]

    for rec in adex_records
        # Skip if no dose value
        rec.aval === nothing && continue

        # Determine dose time (assuming ADEX has time info)
        # In ADaM, timing often comes from ASTDT/AENDT
        # For simplicity, we'll use sequence order and assume first dose at t=0
        dose_time = 0.0

        # Duration (for infusions)
        duration = 0.0
        if rec.adurn !== nothing && rec.adurn > 0
            # Convert days to hours if needed
            if uppercase(rec.aduru) == "DAYS"
                duration = rec.adurn * 24.0
            elseif uppercase(rec.aduru) == "HOURS" || uppercase(rec.aduru) == "H"
                duration = rec.adurn
            else
                duration = rec.adurn  # Assume hours
            end
        end

        push!(doses, DoseEvent(dose_time, rec.aval, duration))
    end

    return doses
end

# ============================================================================
# Main Conversion Functions
# ============================================================================

"""
    adam_to_subject_data(adpc_records, adex_records, adsl_record, options) -> SubjectData

Convert ADaM records for a single subject to SubjectData.

Arguments:
- adpc_records: ADPC records for this subject
- adex_records: ADEX records for this subject
- adsl_record: ADSL record for this subject (optional)
- options: Conversion options

Returns:
- SubjectData object
"""
function adam_to_subject_data(
    adpc_records::Vector{ADPCRecord},
    adex_records::Vector{ADEXRecord},
    adsl_record::Union{Nothing, ADSLRecord},
    options::AdamConversionOptions
)::Union{Nothing, SubjectData}
    isempty(adpc_records) && return nothing

    subject_id = adpc_records[1].usubjid

    # Filter ADPC records based on options
    filtered_records = filter(adpc_records) do rec
        # Filter by analysis flag
        if !isempty(options.analysis_flag)
            if uppercase(rec.anl01fl) != uppercase(options.analysis_flag)
                return false
            end
        end

        # Filter by PARAMCD
        if !isempty(options.filter_paramcd)
            if uppercase(rec.paramcd) != uppercase(options.filter_paramcd)
                return false
            end
        end

        # Filter BLQ if not including
        if !options.include_blq && uppercase(rec.blqfl) == "Y"
            return false
        end

        return true
    end

    isempty(filtered_records) && return nothing

    # Sort by time
    sorted_records = sort(filtered_records, by=r -> something(_get_adpc_time(r, options.time_variable), 0.0))

    # Extract times and observations
    times = Float64[]
    observations = Float64[]
    blq_flags = Bool[]
    lloq_val = 0.0

    for rec in sorted_records
        time = _get_adpc_time(rec, options.time_variable)
        time === nothing && continue

        # Determine if BLQ
        is_blq = uppercase(rec.blqfl) == "Y"

        # Get observation value
        obs_val = if is_blq
            _handle_blq(rec.aval, rec.lloq, options.blq_method)
        else
            rec.aval
        end

        obs_val === nothing && continue

        push!(times, time)
        push!(observations, obs_val)
        push!(blq_flags, is_blq)

        # Track LLOQ
        if rec.lloq !== nothing
            lloq_val = rec.lloq
        end
    end

    isempty(times) && return nothing

    # Convert doses
    doses = _convert_doses(adex_records)

    # If no doses from ADEX, try to infer from ADPC
    if isempty(doses) && !isempty(adpc_records)
        first_dose = adpc_records[1].dose
        if first_dose !== nothing && first_dose > 0
            push!(doses, DoseEvent(0.0, first_dose))
        end
    end

    # Extract covariates from ADSL if available
    covariates = if adsl_record !== nothing
        _extract_covariates(adsl_record)
    else
        Dict{Symbol, Any}()
    end

    return SubjectData(
        subject_id,
        times,
        observations,
        doses;
        covariates=covariates,
        lloq=lloq_val,
        blq_flags=blq_flags
    )
end

"""
    adam_to_observed_data(dataset; options) -> ObservedData

Convert an ADaM dataset to ObservedData format for population analysis.

Arguments:
- dataset: ADaMDataset or ADaMPKDataset
- options: AdamConversionOptions (optional)

Returns:
- ObservedData object suitable for population modeling and VPC analysis
"""
function adam_to_observed_data(
    dataset::Union{ADaMDataset, ADaMPKDataset};
    options::AdamConversionOptions=AdamConversionOptions()
)::ObservedData
    subjects = SubjectData[]

    # Get unique subject IDs from ADPC
    subject_ids = unique([r.usubjid for r in dataset.adpc])

    # Build a lookup for ADSL records
    adsl_lookup = Dict{String, ADSLRecord}()
    for rec in dataset.adsl
        adsl_lookup[rec.usubjid] = rec
    end

    for subj_id in subject_ids
        # Get all records for this subject
        adpc_records = filter(r -> r.usubjid == subj_id, dataset.adpc)
        adex_records = filter(r -> r.usubjid == subj_id, dataset.adex)
        adsl_record = get(adsl_lookup, subj_id, nothing)

        # Convert to SubjectData
        subj_data = adam_to_subject_data(adpc_records, adex_records, adsl_record, options)

        if subj_data !== nothing
            push!(subjects, subj_data)
        end
    end

    # Determine units from data
    conc_units = if !isempty(options.conc_units)
        options.conc_units
    elseif !isempty(dataset.adpc)
        dataset.adpc[1].avalu
    else
        "ng/mL"
    end

    # Get analyte name from PARAMCD
    analyte = if !isempty(dataset.adpc)
        dataset.adpc[1].paramcd
    else
        ""
    end

    return ObservedData(
        subjects;
        study_id=dataset.study_id,
        analyte=analyte,
        units=conc_units,
        time_units=options.time_units
    )
end

# ============================================================================
# Bidirectional Conversion: NeoPKPD to ADaM
# ============================================================================

"""
    observed_data_to_adpc(data::ObservedData; paramcd, param, study_id) -> Vector{ADPCRecord}

Convert ObservedData to ADPC records for export.

Arguments:
- data: ObservedData object
- paramcd: Parameter code (default: "CONC")
- param: Parameter description (default: "Concentration")
- study_id: Study identifier (default from data)

Returns:
- Vector of ADPCRecord objects
"""
function observed_data_to_adpc(
    data::ObservedData;
    paramcd::String="CONC",
    param::String="Concentration",
    study_id::String=""
)::Vector{ADPCRecord}
    records = ADPCRecord[]
    study_id = isempty(study_id) ? data.study_id : study_id
    seq = 0

    for subj in data.subjects
        for (i, (t, obs)) in enumerate(zip(subj.times, subj.observations))
            seq += 1
            is_blq = i <= length(subj.blq_flags) ? subj.blq_flags[i] : false

            record = ADPCRecord(
                studyid = study_id,
                usubjid = subj.subject_id,
                pcseq = seq,
                paramcd = paramcd,
                param = param,
                aval = obs,
                avalu = data.units,
                avalc = is_blq ? "<LLOQ" : string(obs),
                afrlt = t,
                nfrlt = t,
                anl01fl = "Y",
                blqfl = is_blq ? "Y" : "N",
                lloq = subj.lloq > 0 ? subj.lloq : nothing
            )

            push!(records, record)
        end
    end

    return records
end

export observed_data_to_adpc

"""
    nca_results_to_adpp(nca_results; study_id, usubjid, dose) -> Vector{ADPPRecord}

Convert NCA results to ADPP records for export.

Arguments:
- nca_results: NCAResult object from NeoPKPD NCA module
- study_id: Study identifier
- usubjid: Subject identifier
- dose: Dose amount (optional)

Returns:
- Vector of ADPPRecord objects
"""
function nca_results_to_adpp(
    nca_results;
    study_id::String="",
    usubjid::String="",
    dose::Union{Nothing, Float64}=nothing
)::Vector{ADPPRecord}
    records = ADPPRecord[]
    seq = 0

    # Map NCA result fields to PARAMCD
    param_mapping = [
        (:auc_last, "AUCLST", "AUC to Last Measurable Conc"),
        (:auc_inf, "AUCINF", "AUC Extrapolated to Infinity"),
        (:cmax, "CMAX", "Maximum Concentration"),
        (:tmax, "TMAX", "Time to Cmax"),
        (:half_life, "THALF", "Terminal Half-Life"),
        (:lambda_z, "LAMZ", "Terminal Rate Constant"),
        (:cl_f, "CLF", "Apparent Clearance"),
        (:vz_f, "VZF", "Apparent Volume of Distribution"),
        (:mrt, "MRT", "Mean Residence Time"),
        (:auc_extrap_pct, "AUCPEXT", "Percent AUC Extrapolated")
    ]

    for (field, paramcd, param) in param_mapping
        if hasproperty(nca_results, field)
            val = getproperty(nca_results, field)
            if val !== nothing
                seq += 1
                record = ADPPRecord(
                    studyid = study_id,
                    usubjid = usubjid,
                    ppseq = seq,
                    paramcd = paramcd,
                    param = param,
                    parcat1 = "NCA",
                    aval = val,
                    anl01fl = "Y",
                    dose = dose
                )
                push!(records, record)
            end
        end
    end

    return records
end

export nca_results_to_adpp

# ============================================================================
# Validation Functions
# ============================================================================

"""
    validate_adam_dataset(dataset) -> (valid, messages)

Validate an ADaM dataset for common issues.

Returns:
- valid: Boolean indicating if dataset is valid
- messages: Vector of validation messages/warnings
"""
function validate_adam_dataset(dataset::Union{ADaMDataset, ADaMPKDataset})::Tuple{Bool, Vector{String}}
    messages = String[]
    valid = true

    # Check for empty datasets
    if isempty(dataset.adpc)
        push!(messages, "WARNING: No ADPC records found")
    end

    if isempty(dataset.adsl)
        push!(messages, "WARNING: No ADSL records found")
    end

    # Check for subjects in ADPC without ADSL
    adpc_subjects = unique([r.usubjid for r in dataset.adpc])
    adsl_subjects = Set([r.usubjid for r in dataset.adsl])

    for subj in adpc_subjects
        if !(subj in adsl_subjects)
            push!(messages, "WARNING: Subject $subj in ADPC but not in ADSL")
        end
    end

    # Check for missing time values
    missing_times = 0
    for rec in dataset.adpc
        if rec.afrlt === nothing && rec.arrlt === nothing && rec.areltm === nothing
            missing_times += 1
        end
    end
    if missing_times > 0
        push!(messages, "WARNING: $missing_times ADPC records missing time values")
    end

    # Check for missing concentrations
    missing_conc = count(r -> r.aval === nothing && uppercase(r.blqfl) != "Y", dataset.adpc)
    if missing_conc > 0
        push!(messages, "WARNING: $missing_conc non-BLQ ADPC records missing AVAL")
        valid = false
    end

    return valid, messages
end

export validate_adam_dataset
