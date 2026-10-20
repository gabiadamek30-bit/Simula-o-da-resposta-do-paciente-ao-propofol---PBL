# CDISC Data Reader
# Read CSV files containing CDISC/SDTM domain data

export read_cdisc_csv, read_pc_csv, read_ex_csv, read_dm_csv, read_pp_csv
export parse_iso8601_datetime, parse_iso8601_duration

"""
Parse ISO 8601 datetime string to DateTime or nothing.
Supports formats: YYYY-MM-DDTHH:MM:SS, YYYY-MM-DD, etc.
"""
function parse_iso8601_datetime(s::AbstractString)::Union{Nothing,Float64}
    s = strip(s)
    if isempty(s) || s == "." || uppercase(s) == "NA"
        return nothing
    end

    # Try full datetime format: YYYY-MM-DDTHH:MM:SS
    m = match(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})", s)
    if m !== nothing
        year, month, day = parse(Int, m.captures[1]), parse(Int, m.captures[2]), parse(Int, m.captures[3])
        hour, minute, second = parse(Int, m.captures[4]), parse(Int, m.captures[5]), parse(Int, m.captures[6])
        # Return hours since midnight of day 0
        days = (year - 2000) * 365 + (month - 1) * 30 + day  # Approximate
        return days * 24.0 + hour + minute / 60.0 + second / 3600.0
    end

    # Try date only format: YYYY-MM-DD
    m = match(r"^(\d{4})-(\d{2})-(\d{2})$", s)
    if m !== nothing
        year, month, day = parse(Int, m.captures[1]), parse(Int, m.captures[2]), parse(Int, m.captures[3])
        days = (year - 2000) * 365 + (month - 1) * 30 + day
        return days * 24.0
    end

    return nothing
end

"""
Parse ISO 8601 duration string to hours.
Supports formats: PT1H30M, P1D, PT30M, etc.
"""
function parse_iso8601_duration(s::AbstractString)::Float64
    s = strip(s)
    if isempty(s) || s == "." || uppercase(s) == "NA"
        return 0.0
    end

    hours = 0.0

    # Days
    m = match(r"P(\d+)D", s)
    if m !== nothing
        hours += parse(Float64, m.captures[1]) * 24.0
    end

    # Hours
    m = match(r"T.*?(\d+(?:\.\d+)?)H", s)
    if m !== nothing
        hours += parse(Float64, m.captures[1])
    end

    # Minutes
    m = match(r"T.*?(\d+(?:\.\d+)?)M", s)
    if m !== nothing
        hours += parse(Float64, m.captures[1]) / 60.0
    end

    # Seconds
    m = match(r"T.*?(\d+(?:\.\d+)?)S", s)
    if m !== nothing
        hours += parse(Float64, m.captures[1]) / 3600.0
    end

    return hours
end

"""
Parse elapsed time string (PCELTM format) to hours.
Supports: PT1H, PT30M, -PT1H, PT1H30M, etc.
"""
function parse_elapsed_time(s::AbstractString)::Float64
    s = strip(s)
    if isempty(s) || s == "." || uppercase(s) == "NA"
        return 0.0
    end

    sign = 1.0
    if startswith(s, "-")
        sign = -1.0
        s = s[2:end]
    end

    return sign * parse_iso8601_duration(s)
end

"""
Normalize column name to uppercase without special characters.
"""
function normalize_column_name(name::AbstractString)::String
    return uppercase(replace(strip(name), r"[^A-Za-z0-9]" => ""))
end

"""
Safely get a value from a row dict, returning default if not present or empty.
"""
function get_value(row::Dict, key::String, default::T)::T where T
    if !haskey(row, key)
        return default
    end
    val = row[key]
    if val === nothing || val === missing || (val isa AbstractString && (isempty(strip(val)) || val == "." || uppercase(val) == "NA"))
        return default
    end
    if T == Float64 && val isa AbstractString
        return parse(Float64, val)
    elseif T == Int && val isa AbstractString
        return parse(Int, val)
    elseif T == String && !(val isa String)
        return string(val)
    end
    return convert(T, val)
end

function get_value_or_nothing(row::Dict, key::String, ::Type{T})::Union{Nothing,T} where T
    if !haskey(row, key)
        return nothing
    end
    val = row[key]
    if val === nothing || val === missing || (val isa AbstractString && (isempty(strip(val)) || val == "." || uppercase(val) == "NA"))
        return nothing
    end
    if T == Float64 && val isa AbstractString
        return parse(Float64, val)
    elseif T == Int && val isa AbstractString
        return parse(Int, val)
    end
    return convert(T, val)
end

"""
Read a CSV file and return rows as Vector of Dicts with normalized column names.
"""
function read_csv_rows(filepath::AbstractString)::Vector{Dict{String,Any}}
    lines = readlines(filepath)
    if isempty(lines)
        return Dict{String,Any}[]
    end

    # Parse header
    header_line = lines[1]
    # Handle quoted fields
    headers = String[]
    in_quote = false
    current = ""
    for c in header_line
        if c == '"'
            in_quote = !in_quote
        elseif c == ',' && !in_quote
            push!(headers, normalize_column_name(current))
            current = ""
        else
            current *= c
        end
    end
    push!(headers, normalize_column_name(current))

    # Parse data rows
    rows = Dict{String,Any}[]
    for i in 2:length(lines)
        line = lines[i]
        if isempty(strip(line))
            continue
        end

        # Parse fields
        fields = String[]
        in_quote = false
        current = ""
        for c in line
            if c == '"'
                in_quote = !in_quote
            elseif c == ',' && !in_quote
                push!(fields, strip(current))
                current = ""
            else
                current *= c
            end
        end
        push!(fields, strip(current))

        # Create row dict
        row = Dict{String,Any}()
        for (j, header) in enumerate(headers)
            if j <= length(fields)
                row[header] = fields[j]
            else
                row[header] = ""
            end
        end
        push!(rows, row)
    end

    return rows
end

"""
Read PC domain from CSV file.

Arguments:
- filepath: Path to PC domain CSV file

Returns:
- Vector{PCRecord}
"""
function read_pc_csv(filepath::AbstractString)::Vector{PCRecord}
    rows = read_csv_rows(filepath)
    records = PCRecord[]

    for row in rows
        record = PCRecord(
            studyid = get_value(row, "STUDYID", ""),
            usubjid = get_value(row, "USUBJID", ""),
            pcseq = get_value(row, "PCSEQ", 0),
            pctestcd = get_value(row, "PCTESTCD", ""),
            pctest = get_value(row, "PCTEST", ""),
            pcorres = get_value(row, "PCORRES", ""),
            pcorresu = get_value(row, "PCORRESU", ""),
            pcstresn = get_value_or_nothing(row, "PCSTRESN", Float64),
            pcstresu = get_value(row, "PCSTRESU", ""),
            pcstat = get_value(row, "PCSTAT", ""),
            pclloq = get_value_or_nothing(row, "PCLLOQ", Float64),
            pcspec = get_value(row, "PCSPEC", "PLASMA"),
            pcdy = get_value_or_nothing(row, "PCDY", Int),
            pctpt = get_value(row, "PCTPT", ""),
            pctptnum = get_value_or_nothing(row, "PCTPTNUM", Float64),
            pceltm = get_value(row, "PCELTM", ""),
            pcrftdtc = get_value(row, "PCRFTDTC", ""),
            pcstdtc = get_value(row, "PCSTDTC", ""),
            pcblfl = get_value(row, "PCBLFL", ""),
            pcmethod = get_value(row, "PCMETHOD", "")
        )
        push!(records, record)
    end

    return records
end

"""
Read EX domain from CSV file.

Arguments:
- filepath: Path to EX domain CSV file

Returns:
- Vector{EXRecord}
"""
function read_ex_csv(filepath::AbstractString)::Vector{EXRecord}
    rows = read_csv_rows(filepath)
    records = EXRecord[]

    for row in rows
        record = EXRecord(
            studyid = get_value(row, "STUDYID", ""),
            usubjid = get_value(row, "USUBJID", ""),
            exseq = get_value(row, "EXSEQ", 0),
            extrt = get_value(row, "EXTRT", ""),
            excat = get_value(row, "EXCAT", ""),
            exdose = get_value(row, "EXDOSE", 0.0),
            exdosu = get_value(row, "EXDOSU", "mg"),
            exdosfrm = get_value(row, "EXDOSFRM", ""),
            exdosfrq = get_value(row, "EXDOSFRQ", ""),
            exroute = get_value(row, "EXROUTE", ""),
            exstdtc = get_value(row, "EXSTDTC", ""),
            exendtc = get_value(row, "EXENDTC", ""),
            exdy = get_value_or_nothing(row, "EXDY", Int),
            exendy = get_value_or_nothing(row, "EXENDY", Int),
            exdur = get_value(row, "EXDUR", "")
        )
        push!(records, record)
    end

    return records
end

"""
Read DM domain from CSV file.

Arguments:
- filepath: Path to DM domain CSV file

Returns:
- Vector{DMRecord}
"""
function read_dm_csv(filepath::AbstractString)::Vector{DMRecord}
    rows = read_csv_rows(filepath)
    records = DMRecord[]

    for row in rows
        record = DMRecord(
            studyid = get_value(row, "STUDYID", ""),
            usubjid = get_value(row, "USUBJID", ""),
            subjid = get_value(row, "SUBJID", ""),
            rfstdtc = get_value(row, "RFSTDTC", ""),
            rfendtc = get_value(row, "RFENDTC", ""),
            siteid = get_value(row, "SITEID", ""),
            brthdtc = get_value(row, "BRTHDTC", ""),
            age = get_value_or_nothing(row, "AGE", Float64),
            ageu = get_value(row, "AGEU", "YEARS"),
            sex = get_value(row, "SEX", ""),
            race = get_value(row, "RACE", ""),
            ethnic = get_value(row, "ETHNIC", ""),
            armcd = get_value(row, "ARMCD", ""),
            arm = get_value(row, "ARM", ""),
            country = get_value(row, "COUNTRY", ""),
            dmdtc = get_value(row, "DMDTC", ""),
            dmdy = get_value_or_nothing(row, "DMDY", Int)
        )
        push!(records, record)
    end

    return records
end

"""
Read PP domain from CSV file.

Arguments:
- filepath: Path to PP domain CSV file

Returns:
- Vector{PPRecord}
"""
function read_pp_csv(filepath::AbstractString)::Vector{PPRecord}
    rows = read_csv_rows(filepath)
    records = PPRecord[]

    for row in rows
        record = PPRecord(
            studyid = get_value(row, "STUDYID", ""),
            usubjid = get_value(row, "USUBJID", ""),
            ppseq = get_value(row, "PPSEQ", 0),
            ppgrpid = get_value(row, "PPGRPID", ""),
            pptestcd = get_value(row, "PPTESTCD", ""),
            pptest = get_value(row, "PPTEST", ""),
            ppcat = get_value(row, "PPCAT", ""),
            pporres = get_value(row, "PPORRES", ""),
            pporresu = get_value(row, "PPORRESU", ""),
            ppstresn = get_value_or_nothing(row, "PPSTRESN", Float64),
            ppstresu = get_value(row, "PPSTRESU", ""),
            ppstat = get_value(row, "PPSTAT", ""),
            ppreasnd = get_value(row, "PPREASND", ""),
            ppspec = get_value(row, "PPSPEC", "PLASMA"),
            ppdtc = get_value(row, "PPDTC", "")
        )
        push!(records, record)
    end

    return records
end

"""
Read complete CDISC dataset from CSV files.

Arguments:
- pc_path: Path to PC domain CSV (optional)
- ex_path: Path to EX domain CSV (optional)
- dm_path: Path to DM domain CSV (optional)
- pp_path: Path to PP domain CSV (optional)

Returns:
- CDISCDataset
"""
function read_cdisc_csv(;
    pc_path::Union{Nothing,AbstractString}=nothing,
    ex_path::Union{Nothing,AbstractString}=nothing,
    dm_path::Union{Nothing,AbstractString}=nothing,
    pp_path::Union{Nothing,AbstractString}=nothing
)::CDISCDataset
    pc = pc_path !== nothing ? read_pc_csv(pc_path) : PCRecord[]
    ex = ex_path !== nothing ? read_ex_csv(ex_path) : EXRecord[]
    dm = dm_path !== nothing ? read_dm_csv(dm_path) : DMRecord[]
    pp = pp_path !== nothing ? read_pp_csv(pp_path) : PPRecord[]

    # Determine study ID from any available domain
    study_id = ""
    if !isempty(pc)
        study_id = pc[1].studyid
    elseif !isempty(ex)
        study_id = ex[1].studyid
    elseif !isempty(dm)
        study_id = dm[1].studyid
    end

    return CDISCDataset(pc=pc, ex=ex, dm=dm, pp=pp, study_id=study_id)
end
