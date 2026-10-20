# ADaM CSV Reader
# Read ADaM datasets from CSV files

export read_adsl, read_adpc, read_adpp, read_adex
export read_adam_dataset, read_adam_pk_dataset

# ============================================================================
# CSV Parsing Utilities
# ============================================================================

"""
Parse a CSV line handling quoted fields.
"""
function _parse_csv_line(line::String)::Vector{String}
    fields = String[]
    current = ""
    in_quotes = false

    for char in line
        if char == '"'
            in_quotes = !in_quotes
        elseif char == ',' && !in_quotes
            push!(fields, strip(current))
            current = ""
        else
            current *= char
        end
    end
    push!(fields, strip(current))

    return fields
end

"""
Parse a value that might be numeric or missing.
"""
function _parse_numeric(val::String)::Union{Nothing, Float64}
    val = strip(val)
    if isempty(val) || val == "." || lowercase(val) == "na" || lowercase(val) == "nan"
        return nothing
    end
    try
        return parse(Float64, val)
    catch
        return nothing
    end
end

"""
Parse an integer value that might be missing.
"""
function _parse_int(val::String)::Union{Nothing, Int}
    val = strip(val)
    if isempty(val) || val == "." || lowercase(val) == "na"
        return nothing
    end
    try
        return parse(Int, val)
    catch
        # Try parsing as float first, then convert
        try
            return round(Int, parse(Float64, val))
        catch
            return nothing
        end
    end
end

"""
Get column index by name (case-insensitive).
"""
function _get_col_idx(headers::Vector{String}, name::String)::Int
    name_upper = uppercase(name)
    for (i, h) in enumerate(headers)
        if uppercase(strip(h)) == name_upper
            return i
        end
    end
    return 0  # Not found
end

"""
Get value from row by column name.
"""
function _get_val(row::Vector{String}, headers::Vector{String}, name::String)::String
    idx = _get_col_idx(headers, name)
    if idx == 0 || idx > length(row)
        return ""
    end
    return strip(row[idx])
end

# ============================================================================
# ADSL Reader
# ============================================================================

"""
    read_adsl(filepath::String) -> Vector{ADSLRecord}

Read ADSL (Subject Level Analysis Dataset) from a CSV file.

Arguments:
- filepath: Path to the ADSL CSV file

Returns:
- Vector of ADSLRecord objects
"""
function read_adsl(filepath::String)::Vector{ADSLRecord}
    records = ADSLRecord[]

    lines = readlines(filepath)
    isempty(lines) && return records

    headers = _parse_csv_line(lines[1])

    for i in 2:length(lines)
        line = strip(lines[i])
        isempty(line) && continue

        row = _parse_csv_line(line)

        record = ADSLRecord(
            studyid = _get_val(row, headers, "STUDYID"),
            usubjid = _get_val(row, headers, "USUBJID"),
            subjid = _get_val(row, headers, "SUBJID"),
            siteid = _get_val(row, headers, "SITEID"),
            age = _parse_numeric(_get_val(row, headers, "AGE")),
            ageu = _get_val(row, headers, "AGEU"),
            agegr1 = _get_val(row, headers, "AGEGR1"),
            sex = _get_val(row, headers, "SEX"),
            sexn = _parse_int(_get_val(row, headers, "SEXN")),
            race = _get_val(row, headers, "RACE"),
            racen = _parse_int(_get_val(row, headers, "RACEN")),
            ethnic = _get_val(row, headers, "ETHNIC"),
            saffl = _get_val(row, headers, "SAFFL"),
            pkfl = _get_val(row, headers, "PKFL"),
            ittfl = _get_val(row, headers, "ITTFL"),
            pkpopfl = _get_val(row, headers, "PKPOPFL"),
            trt01p = _get_val(row, headers, "TRT01P"),
            trt01a = _get_val(row, headers, "TRT01A"),
            trt01pn = _parse_int(_get_val(row, headers, "TRT01PN")),
            trt01an = _parse_int(_get_val(row, headers, "TRT01AN")),
            trtsdt = _parse_int(_get_val(row, headers, "TRTSDT")),
            trtedt = _parse_int(_get_val(row, headers, "TRTEDT")),
            bweight = _parse_numeric(_get_val(row, headers, "BWEIGHT")),
            bheight = _parse_numeric(_get_val(row, headers, "BHEIGHT")),
            bbmi = _parse_numeric(_get_val(row, headers, "BBMI")),
            bcrcl = _parse_numeric(_get_val(row, headers, "BCRCL")),
            begfr = _parse_numeric(_get_val(row, headers, "BEGFR")),
            bast = _parse_numeric(_get_val(row, headers, "BAST")),
            balt = _parse_numeric(_get_val(row, headers, "BALT")),
            bbili = _parse_numeric(_get_val(row, headers, "BBILI")),
            balb = _parse_numeric(_get_val(row, headers, "BALB")),
            complfl = _get_val(row, headers, "COMPLFL"),
            dcsreas = _get_val(row, headers, "DCSREAS")
        )

        push!(records, record)
    end

    return records
end

# ============================================================================
# ADPC Reader
# ============================================================================

"""
    read_adpc(filepath::String) -> Vector{ADPCRecord}

Read ADPC (PK Concentration Analysis Dataset) from a CSV file.

Arguments:
- filepath: Path to the ADPC CSV file

Returns:
- Vector of ADPCRecord objects
"""
function read_adpc(filepath::String)::Vector{ADPCRecord}
    records = ADPCRecord[]

    lines = readlines(filepath)
    isempty(lines) && return records

    headers = _parse_csv_line(lines[1])

    for i in 2:length(lines)
        line = strip(lines[i])
        isempty(line) && continue

        row = _parse_csv_line(line)

        record = ADPCRecord(
            studyid = _get_val(row, headers, "STUDYID"),
            usubjid = _get_val(row, headers, "USUBJID"),
            pcseq = something(_parse_int(_get_val(row, headers, "PCSEQ")), i-1),
            paramcd = _get_val(row, headers, "PARAMCD"),
            param = _get_val(row, headers, "PARAM"),
            paramn = _parse_int(_get_val(row, headers, "PARAMN")),
            aval = _parse_numeric(_get_val(row, headers, "AVAL")),
            avalu = _get_val(row, headers, "AVALU"),
            avalc = _get_val(row, headers, "AVALC"),
            base = _parse_numeric(_get_val(row, headers, "BASE")),
            chg = _parse_numeric(_get_val(row, headers, "CHG")),
            pchg = _parse_numeric(_get_val(row, headers, "PCHG")),
            adt = _parse_int(_get_val(row, headers, "ADT")),
            atm = _parse_numeric(_get_val(row, headers, "ATM")),
            adtm = _parse_numeric(_get_val(row, headers, "ADTM")),
            atpt = _get_val(row, headers, "ATPT"),
            atptn = _parse_numeric(_get_val(row, headers, "ATPTN")),
            areltm = _parse_numeric(_get_val(row, headers, "ARELTM")),
            nreltm = _parse_numeric(_get_val(row, headers, "NRELTM")),
            afrlt = _parse_numeric(_get_val(row, headers, "AFRLT")),
            nfrlt = _parse_numeric(_get_val(row, headers, "NFRLT")),
            arrlt = _parse_numeric(_get_val(row, headers, "ARRLT")),
            nrrlt = _parse_numeric(_get_val(row, headers, "NRRLT")),
            dosession = _parse_int(_get_val(row, headers, "DOSESSION")),
            doseno = _parse_int(_get_val(row, headers, "DOSENO")),
            dose = _parse_numeric(_get_val(row, headers, "DOSE")),
            doseu = _get_val(row, headers, "DOSEU"),
            trtp = _get_val(row, headers, "TRTP"),
            trta = _get_val(row, headers, "TRTA"),
            anl01fl = _get_val(row, headers, "ANL01FL"),
            anl02fl = _get_val(row, headers, "ANL02FL"),
            blqfl = _get_val(row, headers, "BLQFL"),
            ablfl = _get_val(row, headers, "ABLFL"),
            avisit = _get_val(row, headers, "AVISIT"),
            avisitn = _parse_int(_get_val(row, headers, "AVISITN")),
            lloq = _parse_numeric(_get_val(row, headers, "LLOQ")),
            blq_handling = _get_val(row, headers, "BLQHANDL"),
            spec = _get_val(row, headers, "SPEC")
        )

        push!(records, record)
    end

    return records
end

# ============================================================================
# ADPP Reader
# ============================================================================

"""
    read_adpp(filepath::String) -> Vector{ADPPRecord}

Read ADPP (PK Parameters Analysis Dataset) from a CSV file.

Arguments:
- filepath: Path to the ADPP CSV file

Returns:
- Vector of ADPPRecord objects
"""
function read_adpp(filepath::String)::Vector{ADPPRecord}
    records = ADPPRecord[]

    lines = readlines(filepath)
    isempty(lines) && return records

    headers = _parse_csv_line(lines[1])

    for i in 2:length(lines)
        line = strip(lines[i])
        isempty(line) && continue

        row = _parse_csv_line(line)

        record = ADPPRecord(
            studyid = _get_val(row, headers, "STUDYID"),
            usubjid = _get_val(row, headers, "USUBJID"),
            ppseq = something(_parse_int(_get_val(row, headers, "PPSEQ")), i-1),
            paramcd = _get_val(row, headers, "PARAMCD"),
            param = _get_val(row, headers, "PARAM"),
            paramn = _parse_int(_get_val(row, headers, "PARAMN")),
            parcat1 = _get_val(row, headers, "PARCAT1"),
            parcat2 = _get_val(row, headers, "PARCAT2"),
            aval = _parse_numeric(_get_val(row, headers, "AVAL")),
            avalu = _get_val(row, headers, "AVALU"),
            avalc = _get_val(row, headers, "AVALC"),
            base = _parse_numeric(_get_val(row, headers, "BASE")),
            chg = _parse_numeric(_get_val(row, headers, "CHG")),
            trtp = _get_val(row, headers, "TRTP"),
            trta = _get_val(row, headers, "TRTA"),
            trtpn = _parse_int(_get_val(row, headers, "TRTPN")),
            trtan = _parse_int(_get_val(row, headers, "TRTAN")),
            aperiod = _parse_int(_get_val(row, headers, "APERIOD")),
            aperiodc = _get_val(row, headers, "APERIODC"),
            aseq = _parse_int(_get_val(row, headers, "ASEQ")),
            avisit = _get_val(row, headers, "AVISIT"),
            avisitn = _parse_int(_get_val(row, headers, "AVISITN")),
            anl01fl = _get_val(row, headers, "ANL01FL"),
            anl02fl = _get_val(row, headers, "ANL02FL"),
            ablfl = _get_val(row, headers, "ABLFL"),
            doseno = _parse_int(_get_val(row, headers, "DOSENO")),
            dose = _parse_numeric(_get_val(row, headers, "DOSE")),
            doseu = _get_val(row, headers, "DOSEU"),
            method = _get_val(row, headers, "METHOD")
        )

        push!(records, record)
    end

    return records
end

# ============================================================================
# ADEX Reader
# ============================================================================

"""
    read_adex(filepath::String) -> Vector{ADEXRecord}

Read ADEX (Exposure Analysis Dataset) from a CSV file.

Arguments:
- filepath: Path to the ADEX CSV file

Returns:
- Vector of ADEXRecord objects
"""
function read_adex(filepath::String)::Vector{ADEXRecord}
    records = ADEXRecord[]

    lines = readlines(filepath)
    isempty(lines) && return records

    headers = _parse_csv_line(lines[1])

    for i in 2:length(lines)
        line = strip(lines[i])
        isempty(line) && continue

        row = _parse_csv_line(line)

        record = ADEXRecord(
            studyid = _get_val(row, headers, "STUDYID"),
            usubjid = _get_val(row, headers, "USUBJID"),
            exseq = something(_parse_int(_get_val(row, headers, "EXSEQ")), i-1),
            paramcd = _get_val(row, headers, "PARAMCD"),
            param = _get_val(row, headers, "PARAM"),
            paramn = _parse_int(_get_val(row, headers, "PARAMN")),
            parcat1 = _get_val(row, headers, "PARCAT1"),
            aval = _parse_numeric(_get_val(row, headers, "AVAL")),
            avalu = _get_val(row, headers, "AVALU"),
            avalc = _get_val(row, headers, "AVALC"),
            astdt = _parse_int(_get_val(row, headers, "ASTDT")),
            astdtc = _get_val(row, headers, "ASTDTC"),
            aendt = _parse_int(_get_val(row, headers, "AENDT")),
            aendtc = _get_val(row, headers, "AENDTC"),
            adurn = _parse_numeric(_get_val(row, headers, "ADURN")),
            aduru = _get_val(row, headers, "ADURU"),
            trtp = _get_val(row, headers, "TRTP"),
            trta = _get_val(row, headers, "TRTA"),
            exdosfrm = _get_val(row, headers, "EXDOSFRM"),
            exroute = _get_val(row, headers, "EXROUTE"),
            anl01fl = _get_val(row, headers, "ANL01FL"),
            avisit = _get_val(row, headers, "AVISIT"),
            avisitn = _parse_int(_get_val(row, headers, "AVISITN"))
        )

        push!(records, record)
    end

    return records
end

# ============================================================================
# Combined Readers
# ============================================================================

"""
    read_adam_dataset(; adsl_path, adpc_path, adpp_path, adex_path, study_id) -> ADaMDataset

Read a complete ADaM dataset from CSV files.

Arguments:
- adsl_path: Path to ADSL CSV file (optional)
- adpc_path: Path to ADPC CSV file (optional)
- adpp_path: Path to ADPP CSV file (optional)
- adex_path: Path to ADEX CSV file (optional)
- study_id: Study identifier (default: "")

Returns:
- ADaMDataset object
"""
function read_adam_dataset(;
    adsl_path::Union{Nothing, String}=nothing,
    adpc_path::Union{Nothing, String}=nothing,
    adpp_path::Union{Nothing, String}=nothing,
    adex_path::Union{Nothing, String}=nothing,
    study_id::String=""
)::ADaMDataset
    adsl = adsl_path !== nothing ? read_adsl(adsl_path) : ADSLRecord[]
    adpc = adpc_path !== nothing ? read_adpc(adpc_path) : ADPCRecord[]
    adpp = adpp_path !== nothing ? read_adpp(adpp_path) : ADPPRecord[]
    adex = adex_path !== nothing ? read_adex(adex_path) : ADEXRecord[]

    # Infer study_id from data if not provided
    if isempty(study_id)
        if !isempty(adsl)
            study_id = adsl[1].studyid
        elseif !isempty(adpc)
            study_id = adpc[1].studyid
        end
    end

    return ADaMDataset(
        adsl=adsl,
        adpc=adpc,
        adpp=adpp,
        adex=adex,
        study_id=study_id
    )
end

"""
    read_adam_pk_dataset(; adsl_path, adpc_path, adex_path, study_id) -> ADaMPKDataset

Read an ADaM PK-focused dataset from CSV files.

Arguments:
- adsl_path: Path to ADSL CSV file (optional)
- adpc_path: Path to ADPC CSV file (optional)
- adex_path: Path to ADEX CSV file (optional)
- study_id: Study identifier (default: "")

Returns:
- ADaMPKDataset object
"""
function read_adam_pk_dataset(;
    adsl_path::Union{Nothing, String}=nothing,
    adpc_path::Union{Nothing, String}=nothing,
    adex_path::Union{Nothing, String}=nothing,
    study_id::String=""
)::ADaMPKDataset
    adsl = adsl_path !== nothing ? read_adsl(adsl_path) : ADSLRecord[]
    adpc = adpc_path !== nothing ? read_adpc(adpc_path) : ADPCRecord[]
    adex = adex_path !== nothing ? read_adex(adex_path) : ADEXRecord[]

    # Infer study_id from data if not provided
    if isempty(study_id)
        if !isempty(adsl)
            study_id = adsl[1].studyid
        elseif !isempty(adpc)
            study_id = adpc[1].studyid
        end
    end

    return ADaMPKDataset(
        adsl=adsl,
        adpc=adpc,
        adex=adex,
        study_id=study_id
    )
end
