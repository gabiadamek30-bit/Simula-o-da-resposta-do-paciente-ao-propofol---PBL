# NONMEM Control File Parser
# Parses .ctl files into NONMEMControlFile structures

export parse_nonmem_control, read_nonmem_control
export parse_pk_block, parse_error_block, check_unsupported_constructs

"""
Read and parse a NONMEM control file from disk.

Arguments:
- filepath: Path to .ctl file

Returns:
- NONMEMControlFile structure
"""
function read_nonmem_control(filepath::AbstractString)::NONMEMControlFile
    text = read(filepath, String)
    return parse_nonmem_control(text)
end

"""
Parse NONMEM control file text into a structured representation.

Arguments:
- text: Raw control file text

Returns:
- NONMEMControlFile structure
"""
function parse_nonmem_control(text::String)::NONMEMControlFile
    # Remove comments (lines starting with ; or text after ;)
    lines = String[]
    for line in split(text, '\n')
        # Remove inline comments
        comment_idx = findfirst(';', line)
        if comment_idx !== nothing
            line = line[1:comment_idx-1]
        end
        push!(lines, rstrip(line))
    end

    # Find all $SECTION blocks
    sections = _extract_sections(lines)

    # Parse each section
    problem = _parse_problem(get(sections, "PROBLEM", String[]))
    data = _parse_data(get(sections, "DATA", String[]))
    input = _parse_input(get(sections, "INPUT", String[]))
    subroutines = _parse_subroutines(get(sections, "SUBROUTINES", String[]))
    thetas = _parse_theta(get(sections, "THETA", String[]))
    omegas = _parse_omega(get(sections, "OMEGA", String[]))
    sigmas = _parse_sigma(get(sections, "SIGMA", String[]))
    pk_code = get(sections, "PK", String[])
    error_code = get(sections, "ERROR", String[])
    estimation = _parse_estimation(get(sections, "ESTIMATION", String[]))
    tables = _parse_tables(get(sections, "TABLE", String[]))

    return NONMEMControlFile(
        problem,
        data,
        input,
        subroutines,
        thetas,
        omegas,
        sigmas,
        pk_code,
        error_code,
        estimation,
        tables,
        text
    )
end

"""
Extract sections from control file lines.
Returns Dict mapping section name to lines within that section.
"""
function _extract_sections(lines::Vector{String})::Dict{String,Vector{String}}
    sections = Dict{String,Vector{String}}()
    current_section = nothing
    current_lines = String[]

    for line in lines
        stripped = strip(line)
        if isempty(stripped)
            continue
        end

        # Check if this is a section header
        if startswith(stripped, '$')
            # Save previous section
            if current_section !== nothing
                sections[current_section] = current_lines
            end

            # Parse new section name
            parts = split(stripped)
            section_name = uppercase(replace(parts[1], '$' => ""))

            # Handle $SUB as alias for $SUBROUTINES
            if section_name == "SUB"
                section_name = "SUBROUTINES"
            end

            current_section = section_name
            # Include rest of line if present
            if length(parts) > 1
                current_lines = [join(parts[2:end], " ")]
            else
                current_lines = String[]
            end
        elseif current_section !== nothing
            push!(current_lines, stripped)
        end
    end

    # Save last section
    if current_section !== nothing
        sections[current_section] = current_lines
    end

    return sections
end

"""
Parse \$PROBLEM section.
"""
function _parse_problem(lines::Vector{String})::String
    return join(lines, " ")
end

"""
Parse \$DATA section.
"""
function _parse_data(lines::Vector{String})::Union{Nothing,DataSpec}
    if isempty(lines)
        return nothing
    end

    text = join(lines, " ")
    parts = split(text)

    if isempty(parts)
        return nothing
    end

    filename = String(parts[1])
    ignore = String[]
    accept = String[]

    for i in 2:length(parts)
        part = uppercase(String(parts[i]))
        if startswith(part, "IGNORE=")
            push!(ignore, replace(part, "IGNORE=" => ""))
        elseif startswith(part, "ACCEPT=")
            push!(accept, replace(part, "ACCEPT=" => ""))
        end
    end

    return DataSpec(filename, ignore, accept)
end

"""
Parse \$INPUT section.
"""
function _parse_input(lines::Vector{String})::Vector{InputColumn}
    if isempty(lines)
        return InputColumn[]
    end

    text = join(lines, " ")
    # Handle continuation and multiple spaces
    text = replace(text, r"\s+" => " ")

    columns = InputColumn[]
    for part in split(text)
        part = strip(part)
        if isempty(part)
            continue
        end

        # Check for DROP
        if uppercase(part) == "DROP" || uppercase(part) == "SKIP"
            continue
        end

        # Check for NAME=DROP or NAME=ALIAS format
        if contains(part, "=")
            name, value = split(part, "=", limit=2)
            name = String(strip(name))
            value = String(strip(uppercase(value)))

            if value == "DROP" || value == "SKIP"
                push!(columns, InputColumn(name, true))
            else
                push!(columns, InputColumn(name, false, value))
            end
        else
            push!(columns, InputColumn(String(part), false))
        end
    end

    return columns
end

"""
Parse \$SUBROUTINES section.
"""
function _parse_subroutines(lines::Vector{String})::Union{Nothing,SubroutineSpec}
    if isempty(lines)
        return nothing
    end

    text = uppercase(join(lines, " "))
    text = replace(text, r"\s+" => " ")

    advan = 0
    trans = 1
    other = String[]

    for part in split(text)
        part = strip(part)
        if startswith(part, "ADVAN")
            advan_str = replace(part, "ADVAN" => "")
            advan = parse(Int, advan_str)
        elseif startswith(part, "TRANS")
            trans_str = replace(part, "TRANS" => "")
            trans = parse(Int, trans_str)
        elseif !isempty(part)
            push!(other, part)
        end
    end

    if advan == 0
        return nothing
    end

    return SubroutineSpec(advan, trans, other)
end

"""
Parse \$THETA section.
"""
function _parse_theta(lines::Vector{String})::Vector{THETASpec}
    if isempty(lines)
        return THETASpec[]
    end

    text = join(lines, " ")
    thetas = THETASpec[]

    # Handle parenthesized specifications
    # Formats: (LOWER, INIT, UPPER), (INIT), INIT, (INIT) FIX
    i = 1
    while i <= length(text)
        # Skip whitespace
        while i <= length(text) && isspace(text[i])
            i += 1
        end
        if i > length(text)
            break
        end

        if text[i] == '('
            # Find matching closing paren
            j = findnext(')', text, i)
            if j === nothing
                break
            end

            inner = text[i+1:j-1]
            parts = [strip(p) for p in split(inner, ',')]

            # Check for FIX after closing paren
            rest = j < length(text) ? text[j+1:end] : ""
            fixed = occursin(r"^\s*FIX"i, rest)

            if length(parts) == 1
                init = parse(Float64, parts[1])
                push!(thetas, THETASpec(-Inf, init, Inf, fixed))
            elseif length(parts) == 2
                # Format: (lower, init) - no upper bound
                lower = parts[1] == "" || lowercase(parts[1]) == "-inf" ? -Inf : parse(Float64, parts[1])
                init = parse(Float64, parts[2])
                push!(thetas, THETASpec(lower, init, Inf, fixed))
            elseif length(parts) == 3
                lower = parts[1] == "" || lowercase(parts[1]) == "-inf" ? -Inf : parse(Float64, parts[1])
                init = parse(Float64, parts[2])
                upper = parts[3] == "" || lowercase(parts[3]) == "inf" ? Inf : parse(Float64, parts[3])
                push!(thetas, THETASpec(lower, init, upper, fixed))
            end

            # Move past FIX if present
            if fixed
                i = j + 1
                while i <= length(text) && (isspace(text[i]) || uppercase(text[i]) in ['F','I','X'])
                    i += 1
                end
            else
                i = j + 1
            end
        else
            # Simple number format
            j = i
            while j <= length(text) && (isdigit(text[j]) || text[j] in ['.', '-', '+', 'E', 'e'])
                j += 1
            end
            if j > i
                num_str = text[i:j-1]
                try
                    init = parse(Float64, num_str)
                    push!(thetas, THETASpec(init))
                catch
                end
            end
            i = j
        end
    end

    return thetas
end

"""
Parse \$OMEGA section.
"""
function _parse_omega(lines::Vector{String})::Vector{OMEGABlock}
    if isempty(lines)
        return OMEGABlock[]
    end

    text = join(lines, " ")
    text = replace(text, r"\s+" => " ")

    blocks = OMEGABlock[]
    structure = :diagonal
    fixed = false

    # Check for BLOCK or DIAGONAL
    if occursin(r"BLOCK"i, text)
        structure = :block
        # Extract block size if specified
        m = match(r"BLOCK\s*\((\d+)\)"i, text)
        if m !== nothing
            # block_size = parse(Int, m.captures[1])
        end
        text = replace(text, r"BLOCK\s*\(\d+\)"i => "")
    end

    if occursin(r"DIAGONAL"i, text)
        structure = :diagonal
        text = replace(text, r"DIAGONAL\s*\(\d+\)"i => "")
        text = replace(text, r"DIAGONAL"i => "")
    end

    if occursin(r"\bFIX\b"i, text)
        fixed = true
        text = replace(text, r"\bFIX\b"i => "")
    end

    # Extract numeric values
    values = Float64[]
    for m in eachmatch(r"[\d.eE+-]+", text)
        try
            push!(values, parse(Float64, m.match))
        catch
        end
    end

    if !isempty(values)
        push!(blocks, OMEGABlock(values, structure, fixed))
    end

    return blocks
end

"""
Parse \$SIGMA section.
"""
function _parse_sigma(lines::Vector{String})::Vector{SIGMABlock}
    if isempty(lines)
        return SIGMABlock[]
    end

    text = join(lines, " ")
    text = replace(text, r"\s+" => " ")

    blocks = SIGMABlock[]
    structure = :diagonal
    fixed = false

    if occursin(r"BLOCK"i, text)
        structure = :block
        text = replace(text, r"BLOCK\s*\(\d+\)"i => "")
    end

    if occursin(r"\bFIX\b"i, text)
        fixed = true
        text = replace(text, r"\bFIX\b"i => "")
    end

    # Extract numeric values
    values = Float64[]
    for m in eachmatch(r"[\d.eE+-]+", text)
        try
            push!(values, parse(Float64, m.match))
        catch
        end
    end

    if !isempty(values)
        push!(blocks, SIGMABlock(values, structure, fixed))
    end

    return blocks
end

"""
Parse \$ESTIMATION section.
"""
function _parse_estimation(lines::Vector{String})::Dict{String,Any}
    if isempty(lines)
        return Dict{String,Any}()
    end

    text = uppercase(join(lines, " "))
    result = Dict{String,Any}()

    # Common estimation options
    if occursin("METHOD=0", text) || occursin("METHOD=ZERO", text)
        result["method"] = "FO"
    elseif occursin("METHOD=1", text) || occursin("METHOD=COND", text)
        result["method"] = "FOCE"
    elseif occursin("LAPLACIAN", text) || occursin("LAPLACE", text)
        result["method"] = "LAPLACIAN"
    elseif occursin("SAEM", text)
        result["method"] = "SAEM"
    end

    if occursin("INTERACTION", text) || occursin("INTER", text)
        result["interaction"] = true
    end

    # Parse MAXEVAL
    m = match(r"MAXEVAL\s*=\s*(\d+)", text)
    if m !== nothing
        result["maxeval"] = parse(Int, m.captures[1])
    end

    return result
end

"""
Parse \$TABLE section.
"""
function _parse_tables(lines::Vector{String})::Vector{Dict{String,Any}}
    if isempty(lines)
        return Dict{String,Any}[]
    end

    # For now, just store raw table specifications
    tables = Dict{String,Any}[]
    push!(tables, Dict{String,Any}("raw" => join(lines, " ")))

    return tables
end

# ============================================================================
# $PK Block Parser
# ============================================================================

"""
Parse \$PK block to extract parameter assignments, ETA relationships, and covariates.

Returns a PKBlock containing:
- tv_definitions: Mapping from TV names to THETA indices
- assignments: Parameter assignments with ETA and covariate info
- scaling: Scaling factors (S1, S2, etc.)
- unsupported_lines: Lines that couldn't be parsed
"""
function parse_pk_block(lines::Vector{String})::PKBlock
    if isempty(lines)
        return PKBlock()
    end

    tv_definitions = Dict{Symbol,Int}()
    assignments = PKAssignment[]
    scaling = ScalingFactor[]
    unsupported_lines = String[]

    # Track intermediate values for two-pass processing
    tv_to_theta = Dict{Symbol,Int}()  # TVCL => 1
    tv_covariates = Dict{Symbol,Vector{PKCovariateEffect}}()  # TVCL => [covariates]

    for line in lines
        line_upper = uppercase(strip(line))
        if isempty(line_upper)
            continue
        end

        # Check for unsupported constructs first
        if _is_unsupported_pk_line(line_upper)
            push!(unsupported_lines, line)
            continue
        end

        # Try to parse the line
        parsed = false

        # Pattern 1: TV definition with covariates
        # TVCL = THETA(1) * (WT/70)**THETA(3)
        m = match(r"^(TV\w+)\s*=\s*THETA\((\d+)\)(.*)$"i, line_upper)
        if m !== nothing
            tv_name = Symbol(m.captures[1])
            theta_idx = parse(Int, m.captures[2])
            rest = m.captures[3]

            tv_definitions[tv_name] = theta_idx
            tv_to_theta[tv_name] = theta_idx

            # Check for covariate effects in the rest
            covariates = _extract_covariates_from_expression(rest)
            if !isempty(covariates)
                tv_covariates[tv_name] = covariates
            end
            parsed = true
        end

        # Pattern 2: Parameter with exponential ETA
        # CL = TVCL * EXP(ETA(1))
        if !parsed
            m = match(r"^(\w+)\s*=\s*(TV\w+)\s*\*\s*EXP\s*\(\s*ETA\s*\(\s*(\d+)\s*\)\s*\)"i, line_upper)
            if m !== nothing
                param = Symbol(m.captures[1])
                tv_name = Symbol(m.captures[2])
                eta_idx = parse(Int, m.captures[3])

                # Get THETA index from TV definition
                theta_idx = get(tv_to_theta, tv_name, nothing)
                covs = get(tv_covariates, tv_name, PKCovariateEffect[])

                push!(assignments, PKAssignment(param, theta_idx, eta_idx, :exponential, covs))
                parsed = true
            end
        end

        # Pattern 3: Parameter with additive ETA
        # CL = TVCL + ETA(1)
        if !parsed
            m = match(r"^(\w+)\s*=\s*(TV\w+)\s*\+\s*ETA\s*\(\s*(\d+)\s*\)"i, line_upper)
            if m !== nothing
                param = Symbol(m.captures[1])
                tv_name = Symbol(m.captures[2])
                eta_idx = parse(Int, m.captures[3])

                theta_idx = get(tv_to_theta, tv_name, nothing)
                covs = get(tv_covariates, tv_name, PKCovariateEffect[])

                push!(assignments, PKAssignment(param, theta_idx, eta_idx, :additive, covs))
                parsed = true
            end
        end

        # Pattern 4: Direct parameter from THETA with ETA
        # CL = THETA(1) * EXP(ETA(1))
        if !parsed
            m = match(r"^(\w+)\s*=\s*THETA\s*\(\s*(\d+)\s*\)\s*\*\s*EXP\s*\(\s*ETA\s*\(\s*(\d+)\s*\)\s*\)"i, line_upper)
            if m !== nothing
                param = Symbol(m.captures[1])
                theta_idx = parse(Int, m.captures[2])
                eta_idx = parse(Int, m.captures[3])

                push!(assignments, PKAssignment(param, theta_idx, eta_idx, :exponential, PKCovariateEffect[]))
                parsed = true
            end
        end

        # Pattern 5: Scaling factor
        # S1 = V, S2 = V1
        if !parsed
            m = match(r"^S(\d+)\s*=\s*(\w+)$"i, line_upper)
            if m !== nothing
                comp = parse(Int, m.captures[1])
                param = Symbol(m.captures[2])
                push!(scaling, ScalingFactor(comp, param))
                parsed = true
            end
        end

        # Pattern 6: Simple TV = THETA assignment (already handled in pattern 1)

        # If not parsed and not a comment/empty, log as potentially important
        if !parsed && !startswith(line_upper, ";") && !isempty(line_upper)
            # Some lines are OK to skip (like intermediate calculations)
            # but track ones that look like assignments
            if contains(line_upper, "=")
                # Could be an intermediate variable - that's OK
            end
        end
    end

    return PKBlock(tv_definitions, assignments, scaling, lines, unsupported_lines)
end

"""
Extract covariate effects from an expression like `* (WT/70)**THETA(3)`.
"""
function _extract_covariates_from_expression(expr::AbstractString)::Vector{PKCovariateEffect}
    effects = PKCovariateEffect[]

    # Power covariate: * (COV/REF)**THETA(n)
    for m in eachmatch(r"\*\s*\(\s*(\w+)\s*/\s*([\d.]+)\s*\)\s*\*\*\s*THETA\s*\(\s*(\d+)\s*\)"i, expr)
        cov = Symbol(m.captures[1])
        ref = parse(Float64, m.captures[2])
        theta_idx = parse(Int, m.captures[3])
        push!(effects, PKCovariateEffect(cov, theta_idx, :power, ref))
    end

    # Linear covariate: * (1 + THETA(n)*(COV-REF))
    for m in eachmatch(r"\*\s*\(\s*1\s*\+\s*THETA\s*\(\s*(\d+)\s*\)\s*\*\s*\(\s*(\w+)\s*-\s*([\d.]+)\s*\)\s*\)"i, expr)
        theta_idx = parse(Int, m.captures[1])
        cov = Symbol(m.captures[2])
        ref = parse(Float64, m.captures[3])
        push!(effects, PKCovariateEffect(cov, theta_idx, :linear, ref))
    end

    # Exponential covariate: * EXP(THETA(n)*(COV-REF))
    for m in eachmatch(r"\*\s*EXP\s*\(\s*THETA\s*\(\s*(\d+)\s*\)\s*\*\s*\(\s*(\w+)\s*-\s*([\d.]+)\s*\)\s*\)"i, expr)
        theta_idx = parse(Int, m.captures[1])
        cov = Symbol(m.captures[2])
        ref = parse(Float64, m.captures[3])
        push!(effects, PKCovariateEffect(cov, theta_idx, :exponential, ref))
    end

    return effects
end

"""
Check if a \$PK line contains unsupported constructs.
"""
function _is_unsupported_pk_line(line::AbstractString)::Bool
    # IF statements (but allow safe patterns)
    if occursin(r"\bIF\s*\("i, line)
        if !_is_safe_if_pattern(line)
            return true
        end
    end

    # ALAG (absorption lag time)
    if occursin(r"\bALAG\d"i, line)
        return true
    end

    # F1, F2, etc. (bioavailability fractions)
    if occursin(r"\bF\d\s*="i, line)
        return true
    end

    # MTIME (model event time)
    if occursin(r"\bMTIME"i, line)
        return true
    end

    # R1, R2 (infusion rates)
    if occursin(r"\bR\d\s*="i, line)
        return true
    end

    # D1, D2 (infusion durations)
    if occursin(r"\bD\d\s*="i, line)
        return true
    end

    return false
end

# ============================================================================
# $ERROR Block Parser
# ============================================================================

"""
Parse \$ERROR block to extract error model type and THETA indices.

Returns an ErrorBlock containing:
- error_type: :proportional, :additive, :combined, :exponential, or :unknown
- theta_indices: THETA indices used in the error model
- sigma_fixed_to_1: Whether SIGMA appears to be fixed to 1
"""
function parse_error_block(lines::Vector{String})::ErrorBlock
    if isempty(lines)
        return ErrorBlock()
    end

    text = uppercase(join(lines, " "))
    theta_indices = Int[]
    unsupported_lines = String[]
    error_type = :unknown

    # Check for unsupported constructs
    for line in lines
        line_upper = uppercase(strip(line))
        if _is_unsupported_error_line(line_upper)
            push!(unsupported_lines, line)
        end
    end

    # Extract all THETA references in error block
    for m in eachmatch(r"THETA\s*\(\s*(\d+)\s*\)", text)
        theta_idx = parse(Int, m.captures[1])
        if !(theta_idx in theta_indices)
            push!(theta_indices, theta_idx)
        end
    end

    # Detect error model type based on W definition patterns

    # Combined error with SQRT: W = SQRT(THETA(n)**2 + (THETA(m)*IPRED)**2)
    if occursin(r"W\s*=\s*SQRT\s*\(.*THETA.*\*\*\s*2.*THETA.*IPRED"i, text) ||
       occursin(r"W\s*=\s*SQRT\s*\(.*THETA.*\+.*THETA.*IPRED"i, text)
        error_type = :combined

    # Combined error additive form: W = THETA(n) + THETA(m)*IPRED
    elseif occursin(r"W\s*=\s*THETA\s*\(\d+\)\s*\+\s*THETA\s*\(\d+\)\s*\*\s*(?:IPRED|F)"i, text)
        error_type = :combined

    # Proportional error: W = IPRED * THETA(n) or W = F * THETA(n)
    elseif occursin(r"W\s*=\s*(?:IPRED|F)\s*\*\s*THETA\s*\(\d+\)"i, text)
        error_type = :proportional

    # Additive error: W = THETA(n) (not followed by *)
    elseif occursin(r"W\s*=\s*THETA\s*\(\d+\)\s*$"i, text) ||
           occursin(r"W\s*=\s*THETA\s*\(\d+\)\s+[^*]"i, text)
        error_type = :additive

    # Exponential error: Y = F * EXP(ERR(1))
    elseif occursin(r"Y\s*=\s*(?:IPRED|F)\s*\*\s*EXP\s*\(.*ERR"i, text)
        error_type = :exponential

    # Simple case where error is purely from SIGMA (no THETA in W)
    elseif isempty(theta_indices) && (occursin(r"Y\s*=.*ERR"i, text) || occursin(r"W\s*="i, text))
        # Check what W is set to
        if occursin(r"W\s*=\s*(?:IPRED|F)"i, text)
            error_type = :proportional
        elseif occursin(r"W\s*=\s*1\b"i, text) || occursin(r"W\s*=\s*\d"i, text)
            error_type = :additive
        end
    end

    # Check if SIGMA is fixed to 1 (common pattern when error is parameterized via THETA)
    sigma_fixed_to_1 = !isempty(theta_indices) && error_type != :unknown

    return ErrorBlock(error_type, theta_indices, sigma_fixed_to_1, lines, unsupported_lines)
end

"""
Check if an IF statement is a safe/common pattern that can be ignored.

Safe patterns include:
- IF(W.EQ.0) W = 1 (prevent division by zero)
- IF(W.LE.0) W = 1
- IF(F.EQ.0) F = 1
- IF(IPRED.EQ.0) ... (protect IPRED)
"""
function _is_safe_if_pattern(line::AbstractString)::Bool
    line_upper = uppercase(strip(line))

    # Pattern: IF(X.EQ.0) X = 1 or IF(X.LE.0) X = 1 (division by zero protection)
    if occursin(r"IF\s*\(\s*\w+\s*\.(EQ|LE|LT)\s*\.?\s*0\s*\)\s*\w+\s*=\s*\d"i, line_upper)
        return true
    end

    # Pattern: IF(W.EQ.0) W = 1 specifically
    if occursin(r"IF\s*\(\s*W\s*\.(EQ|LE)\s*\.?\s*0\s*\)"i, line_upper)
        return true
    end

    # Pattern: IF(F.EQ.0) F = ... (bioavailability protection)
    if occursin(r"IF\s*\(\s*F\s*\.(EQ|LE)\s*\.?\s*0\s*\)"i, line_upper)
        return true
    end

    # Pattern: IF(IPRED.EQ.0) or IF(IPRED.LE.0) (prediction protection)
    if occursin(r"IF\s*\(\s*IPRED\s*\.(EQ|LE)\s*\.?\s*0\s*\)"i, line_upper)
        return true
    end

    return false
end

"""
Check if a \$ERROR line contains unsupported constructs.
"""
function _is_unsupported_error_line(line::AbstractString)::Bool
    # IF statements (but allow safe patterns)
    if occursin(r"\bIF\s*\("i, line)
        if !_is_safe_if_pattern(line)
            return true
        end
    end

    # CALL statements
    if occursin(r"\bCALL\b"i, line)
        return true
    end

    return false
end

# ============================================================================
# Unsupported Construct Detection
# ============================================================================

"""
Check a NONMEM control file for unsupported constructs.

Returns a vector of UnsupportedConstruct for each unsupported feature found.
"""
function check_unsupported_constructs(ctl::NONMEMControlFile)::Vector{UnsupportedConstruct}
    unsupported = UnsupportedConstruct[]

    # Check ADVAN number
    if ctl.subroutines !== nothing
        advan = ctl.subroutines.advan
        if advan in [5, 6, 7, 8, 9, 12, 13]
            push!(unsupported, UnsupportedConstruct(
                "ADVAN$advan",
                "\$SUBROUTINES",
                "ADVAN$advan",
                "ADVAN$advan (general linear/nonlinear ODE) is not supported. Only ADVAN1-4, 10, 11 with predefined models are supported."
            ))
        end
    end

    # Check $PK block
    for line in ctl.pk_code
        line_upper = uppercase(strip(line))

        # IF statements (but allow safe patterns)
        if occursin(r"\bIF\s*\("i, line_upper) && !_is_safe_if_pattern(line)
            push!(unsupported, UnsupportedConstruct("IF statement", "\$PK", line))
        end

        if occursin(r"\bALAG\d"i, line_upper)
            push!(unsupported, UnsupportedConstruct("Absorption lag time (ALAG)", "\$PK", line))
        end

        if occursin(r"\bF\d\s*="i, line_upper)
            push!(unsupported, UnsupportedConstruct("Bioavailability fraction (Fn)", "\$PK", line))
        end

        if occursin(r"\bMTIME"i, line_upper)
            push!(unsupported, UnsupportedConstruct("Model event time (MTIME)", "\$PK", line))
        end

        if occursin(r"\bR\d\s*="i, line_upper)
            push!(unsupported, UnsupportedConstruct("Infusion rate (Rn)", "\$PK", line))
        end

        if occursin(r"\bD\d\s*="i, line_upper)
            push!(unsupported, UnsupportedConstruct("Infusion duration (Dn)", "\$PK", line))
        end
    end

    # Check $ERROR block
    for line in ctl.error_code
        line_upper = uppercase(strip(line))

        # IF statements (but allow safe patterns)
        if occursin(r"\bIF\s*\("i, line_upper) && !_is_safe_if_pattern(line)
            push!(unsupported, UnsupportedConstruct("IF statement", "\$ERROR", line))
        end

        if occursin(r"\bCALL\b"i, line_upper)
            push!(unsupported, UnsupportedConstruct("CALL statement", "\$ERROR", line))
        end
    end

    # Check for $DES, $MODEL, $MIX sections (detected by checking raw text)
    raw_upper = uppercase(ctl.raw_text)

    if occursin(r"\$DES\b", raw_upper)
        push!(unsupported, UnsupportedConstruct(
            "Custom differential equations (\$DES)",
            "\$DES",
            "",
            "Custom differential equations (\$DES) are not supported. Use predefined ADVAN models."
        ))
    end

    if occursin(r"\$MODEL\b", raw_upper)
        push!(unsupported, UnsupportedConstruct(
            "Custom compartment model (\$MODEL)",
            "\$MODEL",
            "",
            "Custom compartment models (\$MODEL) are not supported. Use predefined ADVAN models."
        ))
    end

    if occursin(r"\$MIX\b", raw_upper)
        push!(unsupported, UnsupportedConstruct(
            "Mixture model (\$MIX)",
            "\$MIX",
            "",
            "Mixture models (\$MIX) are not supported."
        ))
    end

    # Check for SAME in OMEGA
    for line in get(Dict("OMEGA" => String[]), "OMEGA", ctl.raw_text)
        if occursin(r"\bSAME\b"i, line)
            push!(unsupported, UnsupportedConstruct("SAME omega structure", "\$OMEGA", line))
            break
        end
    end
    # Also check in raw text for SAME
    if occursin(r"\$OMEGA[^\$]*\bSAME\b"i, raw_upper)
        # Only add if not already added
        if !any(u -> u.construct == "SAME omega structure" for u in unsupported)
            push!(unsupported, UnsupportedConstruct(
                "SAME omega structure",
                "\$OMEGA",
                "",
                "SAME omega structure is not supported."
            ))
        end
    end

    return unsupported
end
