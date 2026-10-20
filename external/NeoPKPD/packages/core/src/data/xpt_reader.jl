# XPT (SAS Transport) File Reader
# Pure Julia implementation for CDISC XPT format

export read_xpt, read_cdisc_xpt, XPTVariable, XPTDataset

"""
Variable metadata from XPT file.
"""
struct XPTVariable
    name::String
    ntype::Int  # 1=numeric, 2=character
    length::Int
    label::String
    format::String
    informat::String
end

"""
Dataset read from XPT file.
"""
struct XPTDataset
    name::String
    label::String
    variables::Vector{XPTVariable}
    data::Vector{Dict{String,Any}}
    nobs::Int
end

# XPT format constants
const XPT_HEADER_SIZE = 80
const XPT_NAMESTR_SIZE = 140
const IBM_MAGIC = UInt8[0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48]  # "ABCDEFGH" in EBCDIC

"""
Convert IBM 370 floating point to IEEE 754 double.
XPT uses IBM mainframe floating point format.
"""
function ibm_to_ieee(bytes::Vector{UInt8})::Float64
    if length(bytes) != 8
        return 0.0
    end

    # Check for zero
    if all(b -> b == 0x00, bytes)
        return 0.0
    end

    # Check for missing value (all 0x2E which is '.' in EBCDIC)
    if all(b -> b == 0x2E, bytes)
        return NaN
    end

    # IBM 370 floating point format:
    # Bit 0: sign (0=positive, 1=negative)
    # Bits 1-7: exponent (excess 64, base 16)
    # Bits 8-63: mantissa (normalized fraction)

    sign = (bytes[1] & 0x80) != 0
    exponent = Int(bytes[1] & 0x7F) - 64  # Remove excess 64

    # Build mantissa from remaining bytes
    mantissa = 0.0
    for i in 2:8
        mantissa = mantissa * 256.0 + Float64(bytes[i])
    end

    # Normalize: mantissa is in range [0, 256^7)
    # So we divide by 256^7 to get a fraction
    mantissa = mantissa / (256.0^7)

    # IBM uses base 16 for exponent
    value = mantissa * (16.0^exponent)

    return sign ? -value : value
end

"""
Convert EBCDIC to ASCII.
"""
function ebcdic_to_ascii(bytes::Vector{UInt8})::String
    # EBCDIC to ASCII mapping table for common characters
    ebcdic_map = Dict{UInt8,Char}(
        0x40 => ' ',
        0x4B => '.', 0x4C => '<', 0x4D => '(', 0x4E => '+', 0x4F => '|',
        0x50 => '&',
        0x5A => '!', 0x5B => '$', 0x5C => '*', 0x5D => ')', 0x5E => ';',
        0x5F => '^',
        0x60 => '-', 0x61 => '/',
        0x6B => ',', 0x6C => '%', 0x6D => '_', 0x6E => '>', 0x6F => '?',
        0x7A => ':', 0x7B => '#', 0x7C => '@', 0x7D => '\'', 0x7E => '=', 0x7F => '"',
        0x81 => 'a', 0x82 => 'b', 0x83 => 'c', 0x84 => 'd', 0x85 => 'e',
        0x86 => 'f', 0x87 => 'g', 0x88 => 'h', 0x89 => 'i',
        0x91 => 'j', 0x92 => 'k', 0x93 => 'l', 0x94 => 'm', 0x95 => 'n',
        0x96 => 'o', 0x97 => 'p', 0x98 => 'q', 0x99 => 'r',
        0xA1 => '~', 0xA2 => 's', 0xA3 => 't', 0xA4 => 'u', 0xA5 => 'v',
        0xA6 => 'w', 0xA7 => 'x', 0xA8 => 'y', 0xA9 => 'z',
        0xC0 => '{', 0xC1 => 'A', 0xC2 => 'B', 0xC3 => 'C', 0xC4 => 'D',
        0xC5 => 'E', 0xC6 => 'F', 0xC7 => 'G', 0xC8 => 'H', 0xC9 => 'I',
        0xD0 => '}', 0xD1 => 'J', 0xD2 => 'K', 0xD3 => 'L', 0xD4 => 'M',
        0xD5 => 'N', 0xD6 => 'O', 0xD7 => 'P', 0xD8 => 'Q', 0xD9 => 'R',
        0xE0 => '\\', 0xE2 => 'S', 0xE3 => 'T', 0xE4 => 'U', 0xE5 => 'V',
        0xE6 => 'W', 0xE7 => 'X', 0xE8 => 'Y', 0xE9 => 'Z',
        0xF0 => '0', 0xF1 => '1', 0xF2 => '2', 0xF3 => '3', 0xF4 => '4',
        0xF5 => '5', 0xF6 => '6', 0xF7 => '7', 0xF8 => '8', 0xF9 => '9',
    )

    result = Char[]
    for b in bytes
        if b == 0x00
            break
        end
        push!(result, get(ebcdic_map, b, ' '))
    end
    return String(result)
end

"""
Read a fixed-length string from bytes, trimming trailing spaces.
"""
function read_string(bytes::Vector{UInt8}, start::Int, len::Int)::String
    if start + len - 1 > length(bytes)
        return ""
    end
    s = String(bytes[start:start+len-1])
    return strip(s)
end

"""
Read a big-endian 16-bit integer.
"""
function read_int16_be(bytes::Vector{UInt8}, pos::Int)::Int
    return Int(bytes[pos]) * 256 + Int(bytes[pos+1])
end

"""
Read a big-endian 32-bit integer.
"""
function read_int32_be(bytes::Vector{UInt8}, pos::Int)::Int
    return (Int(bytes[pos]) << 24) + (Int(bytes[pos+1]) << 16) +
           (Int(bytes[pos+2]) << 8) + Int(bytes[pos+3])
end

"""
Parse the variable descriptor (namestr) record.
"""
function parse_namestr(bytes::Vector{UInt8})::XPTVariable
    # Namestr record format (140 bytes):
    # 0-1: ntype (short)
    # 2-3: nhfun (short)
    # 4-5: nlng (short) - length of variable
    # 6-7: nvar0 (short) - variable number
    # 8-15: nname (8 bytes) - variable name
    # 16-55: nlabel (40 bytes) - variable label
    # 56-63: nform (8 bytes) - format name
    # 64-65: nfl (short) - format field length
    # 66-67: nfd (short) - format decimal places
    # 68-69: nfj (short) - format justification
    # 70-71: nfill (short) - unused
    # 72-79: niform (8 bytes) - informat name
    # 80-81: nifl (short) - informat field length
    # 82-83: nifd (short) - informat decimal places
    # 84-91: npos (8 bytes) - position in observation
    # ... remaining bytes unused

    ntype = read_int16_be(bytes, 1)
    nlng = read_int16_be(bytes, 5)
    nname = read_string(bytes, 9, 8)
    nlabel = read_string(bytes, 17, 40)
    nform = read_string(bytes, 57, 8)
    niform = read_string(bytes, 73, 8)

    return XPTVariable(
        name = uppercase(nname),
        ntype = ntype,
        length = nlng,
        label = nlabel,
        format = nform,
        informat = niform
    )
end

"""
Read XPT file and return dataset.

Arguments:
- filepath: Path to XPT file

Returns:
- XPTDataset containing variables and data
"""
function read_xpt(filepath::AbstractString)::XPTDataset
    bytes = read(filepath)
    pos = 1

    # Read library header (first 80 bytes)
    # Should start with "HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!"
    library_header = read_string(bytes, pos, 80)
    pos += 80

    # Verify this is an XPT file
    if !startswith(library_header, "HEADER RECORD")
        error("Invalid XPT file: missing library header")
    end

    # Skip to real header (second 80-byte record)
    # Contains SAS version info
    pos += 80

    # Third record: created/modified dates
    pos += 80

    # Fourth record: member header record marker
    member_marker = read_string(bytes, pos, 80)
    pos += 80

    if !contains(member_marker, "MEMBER")
        error("Invalid XPT file: missing member header marker at position $(pos-80)")
    end

    # Fifth record: dataset header record marker
    dataset_marker = read_string(bytes, pos, 80)
    pos += 80

    # Sixth record: member header info
    # Contains dataset name (8 bytes at position 8)
    member_info = bytes[pos:pos+79]
    dataset_name = read_string(member_info, 9, 8)
    pos += 80

    # Seventh record: dataset label and type
    label_record = bytes[pos:pos+79]
    dataset_label = read_string(label_record, 1, 40)
    pos += 80

    # Namestr header record
    namestr_marker = read_string(bytes, pos, 80)
    pos += 80

    if !contains(namestr_marker, "NAMESTR")
        error("Invalid XPT file: missing namestr header marker")
    end

    # Parse number of variables from the namestr marker
    # Format: "HEADER RECORD*******NAMESTR HEADER RECORD!!!!!!!000000XXXX"
    # where XXXX is the number of variables (padded with spaces, right-justified)
    m = match(r"(\d+)\s*$", namestr_marker)
    nvars = m !== nothing ? parse(Int, m.captures[1]) : 0

    if nvars == 0
        # Try alternative format
        m = match(r"!{7}0*(\d+)", namestr_marker)
        nvars = m !== nothing ? parse(Int, m.captures[1]) : 0
    end

    # Read variable descriptors (namestr records)
    # Each is 140 bytes
    variables = XPTVariable[]
    for i in 1:nvars
        namestr_bytes = bytes[pos:pos+139]
        push!(variables, parse_namestr(namestr_bytes))
        pos += 140
    end

    # Pad to 80-byte boundary
    if (nvars * 140) % 80 != 0
        padding = 80 - ((nvars * 140) % 80)
        pos += padding
    end

    # Observation header record
    obs_marker = read_string(bytes, pos, 80)
    pos += 80

    if !contains(obs_marker, "OBS")
        # Some files don't have observation marker, backtrack
        pos -= 80
    end

    # Calculate observation length
    obs_length = sum(v -> v.ntype == 1 ? 8 : v.length, variables)

    # Round up to 80-byte boundary for record length
    # XPT packs multiple observations per record

    # Read observations
    data = Dict{String,Any}[]

    while pos + obs_length <= length(bytes) + 1
        row = Dict{String,Any}()
        var_pos = pos

        for var in variables
            if var.ntype == 1  # Numeric
                if var_pos + 7 <= length(bytes)
                    num_bytes = bytes[var_pos:var_pos+7]
                    value = ibm_to_ieee(num_bytes)
                    row[var.name] = isnan(value) ? nothing : value
                else
                    row[var.name] = nothing
                end
                var_pos += 8
            else  # Character
                if var_pos + var.length - 1 <= length(bytes)
                    char_bytes = bytes[var_pos:var_pos+var.length-1]
                    row[var.name] = strip(String(char_bytes))
                else
                    row[var.name] = ""
                end
                var_pos += var.length
            end
        end

        # Check if this is a padding row (all empty)
        all_empty = all(v -> v === nothing || (v isa String && isempty(v)), values(row))

        if !all_empty
            push!(data, row)
        end

        pos += obs_length
    end

    return XPTDataset(
        name = dataset_name,
        label = dataset_label,
        variables = variables,
        data = data,
        nobs = length(data)
    )
end

"""
Read CDISC domain from XPT file.

Arguments:
- filepath: Path to XPT file containing CDISC domain

Returns:
- Vector of rows as Dict{String,Any}
"""
function read_cdisc_xpt(filepath::AbstractString)::Vector{Dict{String,Any}}
    dataset = read_xpt(filepath)
    return dataset.data
end

"""
Read PC domain from XPT file and convert to PCRecord vector.

Arguments:
- filepath: Path to PC XPT file

Returns:
- Vector{PCRecord}
"""
function read_pc_xpt(filepath::AbstractString)::Vector{PCRecord}
    rows = read_cdisc_xpt(filepath)
    records = PCRecord[]

    for row in rows
        # Normalize column names to uppercase
        norm_row = Dict{String,Any}(uppercase(string(k)) => v for (k, v) in row)

        record = PCRecord(
            studyid = get_value(norm_row, "STUDYID", ""),
            usubjid = get_value(norm_row, "USUBJID", ""),
            pcseq = get_value(norm_row, "PCSEQ", 0),
            pctestcd = get_value(norm_row, "PCTESTCD", ""),
            pctest = get_value(norm_row, "PCTEST", ""),
            pcorres = get_value(norm_row, "PCORRES", ""),
            pcorresu = get_value(norm_row, "PCORRESU", ""),
            pcstresn = get_value_or_nothing(norm_row, "PCSTRESN", Float64),
            pcstresu = get_value(norm_row, "PCSTRESU", ""),
            pcstat = get_value(norm_row, "PCSTAT", ""),
            pclloq = get_value_or_nothing(norm_row, "PCLLOQ", Float64),
            pcspec = get_value(norm_row, "PCSPEC", "PLASMA"),
            pcdy = get_value_or_nothing(norm_row, "PCDY", Int),
            pctpt = get_value(norm_row, "PCTPT", ""),
            pctptnum = get_value_or_nothing(norm_row, "PCTPTNUM", Float64),
            pceltm = get_value(norm_row, "PCELTM", ""),
            pcrftdtc = get_value(norm_row, "PCRFTDTC", ""),
            pcstdtc = get_value(norm_row, "PCSTDTC", ""),
            pcblfl = get_value(norm_row, "PCBLFL", ""),
            pcmethod = get_value(norm_row, "PCMETHOD", "")
        )
        push!(records, record)
    end

    return records
end

"""
Read EX domain from XPT file and convert to EXRecord vector.

Arguments:
- filepath: Path to EX XPT file

Returns:
- Vector{EXRecord}
"""
function read_ex_xpt(filepath::AbstractString)::Vector{EXRecord}
    rows = read_cdisc_xpt(filepath)
    records = EXRecord[]

    for row in rows
        norm_row = Dict{String,Any}(uppercase(string(k)) => v for (k, v) in row)

        record = EXRecord(
            studyid = get_value(norm_row, "STUDYID", ""),
            usubjid = get_value(norm_row, "USUBJID", ""),
            exseq = get_value(norm_row, "EXSEQ", 0),
            extrt = get_value(norm_row, "EXTRT", ""),
            excat = get_value(norm_row, "EXCAT", ""),
            exdose = get_value(norm_row, "EXDOSE", 0.0),
            exdosu = get_value(norm_row, "EXDOSU", "mg"),
            exdosfrm = get_value(norm_row, "EXDOSFRM", ""),
            exdosfrq = get_value(norm_row, "EXDOSFRQ", ""),
            exroute = get_value(norm_row, "EXROUTE", ""),
            exstdtc = get_value(norm_row, "EXSTDTC", ""),
            exendtc = get_value(norm_row, "EXENDTC", ""),
            exdy = get_value_or_nothing(norm_row, "EXDY", Int),
            exendy = get_value_or_nothing(norm_row, "EXENDY", Int),
            exdur = get_value(norm_row, "EXDUR", "")
        )
        push!(records, record)
    end

    return records
end

"""
Read DM domain from XPT file and convert to DMRecord vector.

Arguments:
- filepath: Path to DM XPT file

Returns:
- Vector{DMRecord}
"""
function read_dm_xpt(filepath::AbstractString)::Vector{DMRecord}
    rows = read_cdisc_xpt(filepath)
    records = DMRecord[]

    for row in rows
        norm_row = Dict{String,Any}(uppercase(string(k)) => v for (k, v) in row)

        record = DMRecord(
            studyid = get_value(norm_row, "STUDYID", ""),
            usubjid = get_value(norm_row, "USUBJID", ""),
            subjid = get_value(norm_row, "SUBJID", ""),
            rfstdtc = get_value(norm_row, "RFSTDTC", ""),
            rfendtc = get_value(norm_row, "RFENDTC", ""),
            siteid = get_value(norm_row, "SITEID", ""),
            brthdtc = get_value(norm_row, "BRTHDTC", ""),
            age = get_value_or_nothing(norm_row, "AGE", Float64),
            ageu = get_value(norm_row, "AGEU", "YEARS"),
            sex = get_value(norm_row, "SEX", ""),
            race = get_value(norm_row, "RACE", ""),
            ethnic = get_value(norm_row, "ETHNIC", ""),
            armcd = get_value(norm_row, "ARMCD", ""),
            arm = get_value(norm_row, "ARM", ""),
            country = get_value(norm_row, "COUNTRY", ""),
            dmdtc = get_value(norm_row, "DMDTC", ""),
            dmdy = get_value_or_nothing(norm_row, "DMDY", Int)
        )
        push!(records, record)
    end

    return records
end

"""
Read complete CDISC dataset from XPT files.

Arguments:
- pc_path: Path to PC domain XPT file (optional)
- ex_path: Path to EX domain XPT file (optional)
- dm_path: Path to DM domain XPT file (optional)
- pp_path: Path to PP domain XPT file (optional)

Returns:
- CDISCDataset
"""
function read_cdisc_from_xpt(;
    pc_path::Union{Nothing,AbstractString}=nothing,
    ex_path::Union{Nothing,AbstractString}=nothing,
    dm_path::Union{Nothing,AbstractString}=nothing,
    pp_path::Union{Nothing,AbstractString}=nothing
)::CDISCDataset
    pc = pc_path !== nothing ? read_pc_xpt(pc_path) : PCRecord[]
    ex = ex_path !== nothing ? read_ex_xpt(ex_path) : EXRecord[]
    dm = dm_path !== nothing ? read_dm_xpt(dm_path) : DMRecord[]
    pp = PPRecord[]  # PP XPT reader would follow same pattern

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
