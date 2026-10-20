# Monolix Project File Parser
# Parses .mlxtran files into MonolixProject structures
# Supports both legacy [SECTION] format and modern <SECTION> format

export parse_monolix_project, read_monolix_project

"""
Read and parse a Monolix project file from disk.

Arguments:
- filepath: Path to .mlxtran file

Returns:
- MonolixProject structure
"""
function read_monolix_project(filepath::AbstractString)::MonolixProject
    text = read(filepath, String)
    return parse_monolix_project(text)
end

"""
Parse Monolix project file text into a structured representation.

Supports both formats:
- Legacy: [SECTION_NAME]
- Modern (Monolix 2019+): <SECTION_NAME> with [SUBSECTION] inside

Arguments:
- text: Raw .mlxtran file text

Returns:
- MonolixProject structure
"""
function parse_monolix_project(text::String)::MonolixProject
    # Detect format and extract blocks
    if occursin(r"<\w+>", text)
        # Modern format with angle brackets
        blocks = _extract_modern_monolix_blocks(text)
    else
        # Legacy format with square brackets only
        blocks = _extract_legacy_monolix_blocks(text)
    end

    # Parse each section based on available blocks
    description = _parse_description(blocks)
    data = _parse_datafile(blocks)
    model = _parse_model(blocks)
    parameters = _parse_parameters(blocks)

    # Extract raw parameter values for observation parsing
    param_values = _extract_parameter_values(blocks)
    observations = _parse_observations(blocks, param_values)
    estimation = _parse_estimation(blocks)

    return MonolixProject(
        description,
        data,
        model,
        parameters,
        observations,
        estimation,
        text
    )
end

# ============================================================================
# Modern Format Parsing (<SECTION> with [SUBSECTION])
# ============================================================================

"""
Extract blocks from modern Monolix format using <SECTION> tags.
Returns nested Dict with section -> subsection -> content
"""
function _extract_modern_monolix_blocks(text::String)::Dict{String,Any}
    blocks = Dict{String,Any}()
    current_section = nothing
    current_subsection = nothing
    section_content = Dict{String,Vector{String}}()
    subsection_lines = String[]

    for line in split(text, '\n')
        stripped = strip(line)

        # Check for section start <SECTION>
        m = match(r"^<(\w+)>$", stripped)
        if m !== nothing
            # Save previous section
            if current_section !== nothing
                if current_subsection !== nothing
                    section_content[current_subsection] = subsection_lines
                end
                blocks[current_section] = section_content
            end
            current_section = uppercase(String(m.captures[1]))
            section_content = Dict{String,Vector{String}}()
            current_subsection = nothing
            subsection_lines = String[]
            continue
        end

        # Check for subsection start [SUBSECTION]
        m = match(r"^\[(\w+)\]$", stripped)
        if m !== nothing && current_section !== nothing
            # Save previous subsection
            if current_subsection !== nothing
                section_content[current_subsection] = subsection_lines
            end
            current_subsection = uppercase(String(m.captures[1]))
            subsection_lines = String[]
            continue
        end

        # Add line to current content
        if current_section !== nothing
            if current_subsection !== nothing
                push!(subsection_lines, stripped)
            else
                # Lines before any subsection go to "_ROOT"
                if !haskey(section_content, "_ROOT")
                    section_content["_ROOT"] = String[]
                end
                push!(section_content["_ROOT"], stripped)
            end
        end
    end

    # Save last section
    if current_section !== nothing
        if current_subsection !== nothing
            section_content[current_subsection] = subsection_lines
        end
        blocks[current_section] = section_content
    end

    return blocks
end

"""
Extract blocks from legacy Monolix format using [SECTION] tags only.
"""
function _extract_legacy_monolix_blocks(text::String)::Dict{String,Any}
    blocks = Dict{String,Any}()
    current_block = nothing
    current_content = String[]

    for line in split(text, '\n')
        stripped = strip(line)

        # Check for block start [SECTION]
        m = match(r"^\[([A-Z_]+)\]$", stripped)
        if m !== nothing
            # Save previous block
            if current_block !== nothing
                blocks[current_block] = Dict("_ROOT" => current_content)
            end
            current_block = String(m.captures[1])
            current_content = String[]
        elseif current_block !== nothing
            push!(current_content, stripped)
        end
    end

    # Save last block
    if current_block !== nothing
        blocks[current_block] = Dict("_ROOT" => current_content)
    end

    return blocks
end

# ============================================================================
# Section Parsers
# ============================================================================

"""
Parse description from blocks.
"""
function _parse_description(blocks::Dict{String,Any})::String
    if haskey(blocks, "DESCRIPTION")
        content = blocks["DESCRIPTION"]
        if content isa Dict && haskey(content, "_ROOT")
            # Filter out empty lines and join with space
            lines = filter(!isempty, content["_ROOT"])
            return strip(join(lines, " "))
        end
    end
    return ""
end

"""
Parse <DATAFILE> section.
"""
function _parse_datafile(blocks::Dict{String,Any})::Union{Nothing,MonolixDataset}
    if !haskey(blocks, "DATAFILE")
        return nothing
    end

    content = blocks["DATAFILE"]
    filename = ""
    header_types = Dict{String,String}()
    id_col = "ID"
    time_col = "TIME"
    obs_col = "DV"
    dose_col = "AMT"
    rate_col = "RATE"

    # Parse [FILEINFO] subsection
    if content isa Dict && haskey(content, "FILEINFO")
        for line in content["FILEINFO"]
            if isempty(line) || startswith(line, ";")
                continue
            end

            # Parse file = 'path'
            m = match(r"file\s*=\s*'([^']+)'", line)
            if m !== nothing
                filename = String(m.captures[1])
                continue
            end

            # Parse header = {...}
            m = match(r"header\s*=\s*\{([^}]+)\}", line)
            if m !== nothing
                for col in split(m.captures[1], ',')
                    col = strip(String(col))
                    if !isempty(col)
                        header_types[col] = col
                    end
                end
            end
        end
    end

    # Parse [CONTENT] subsection for column uses
    if content isa Dict && haskey(content, "CONTENT")
        for line in content["CONTENT"]
            if isempty(line) || startswith(line, ";")
                continue
            end

            # Parse COLUMN = {use=type, ...}
            m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
            if m !== nothing
                col_name = String(m.captures[1])
                props = String(m.captures[2])

                # Extract use type
                use_match = match(r"use\s*=\s*(\w+)", props)
                if use_match !== nothing
                    use_type = lowercase(String(use_match.captures[1]))
                    header_types[col_name] = use_type

                    if use_type == "identifier"
                        id_col = col_name
                    elseif use_type == "time"
                        time_col = col_name
                    elseif use_type == "observation"
                        obs_col = col_name
                    elseif use_type == "amount"
                        dose_col = col_name
                    elseif use_type == "rate"
                        rate_col = col_name
                    end
                end
            end
        end
    end

    # Legacy format fallback
    if content isa Dict && haskey(content, "_ROOT") && isempty(filename)
        for line in content["_ROOT"]
            m = match(r"file\s*=\s*'([^']+)'", line)
            if m !== nothing
                filename = String(m.captures[1])
            end
        end
    end

    if isempty(filename)
        return nothing
    end

    return MonolixDataset(
        filename;
        header_types=header_types,
        id_column=id_col,
        time_column=time_col,
        observation_column=obs_col,
        dose_column=dose_col,
        rate_column=rate_col
    )
end

"""
Parse <MODEL> section.
"""
function _parse_model(blocks::Dict{String,Any})::Union{Nothing,MonolixStructuralModel}
    model_type = nothing
    admin_type = "iv"
    n_compartments = 1
    elimination = "linear"
    absorption = "bolus"
    has_lag = false
    has_bio = false

    # Parse from MODEL block
    if haskey(blocks, "MODEL")
        content = blocks["MODEL"]

        # Parse [LONGITUDINAL] for model file
        if content isa Dict && haskey(content, "LONGITUDINAL")
            for line in content["LONGITUDINAL"]
                if isempty(line) || startswith(line, ";")
                    continue
                end

                # Parse file = 'lib:model_name.txt'
                m = match(r"file\s*=\s*'lib:([^']+)'", line)
                if m !== nothing
                    model_name = String(m.captures[1])
                    # Remove .txt extension if present
                    model_name = replace(model_name, ".txt" => "")
                    model_type = MonolixModelType("pklib", model_name)

                    # Infer properties from model name
                    _infer_model_properties_from_name!(
                        model_name, admin_type, n_compartments,
                        elimination, absorption, has_lag, has_bio
                    )

                    # Update local variables based on inference
                    model_lower = lowercase(model_name)
                    if occursin("oral", model_lower)
                        admin_type = "oral"
                        absorption = "firstOrder"
                    end
                    if occursin("1cpt", model_lower) || occursin("1_1cpt", model_lower)
                        n_compartments = 1
                    elseif occursin("2cpt", model_lower) || occursin("2_2cpt", model_lower)
                        n_compartments = 2
                    elseif occursin("3cpt", model_lower)
                        n_compartments = 3
                    end
                    if occursin("tlag", model_lower)
                        has_lag = true
                    end
                    if occursin("bio", model_lower)
                        has_bio = true
                    end
                    if occursin("mm", model_lower) || occursin("michaelis", model_lower)
                        elimination = "mm"
                    end
                end
            end
        end
    end

    # Legacy format: STRUCTURAL_MODEL block
    if haskey(blocks, "STRUCTURAL_MODEL")
        content = blocks["STRUCTURAL_MODEL"]
        lines = content isa Dict && haskey(content, "_ROOT") ? content["_ROOT"] : String[]

        for line in lines
            if isempty(line)
                continue
            end

            # Parse lib = pklib:model_name
            m = match(r"lib\s*=\s*(\w+):(\S+)", line)
            if m !== nothing
                model_type = MonolixModelType(String(m.captures[1]), String(m.captures[2]))
                model_name = lowercase(String(m.captures[2]))

                if occursin("oral", model_name)
                    admin_type = "oral"
                    absorption = "firstOrder"
                end
                if occursin("1cpt", model_name)
                    n_compartments = 1
                elseif occursin("2cpt", model_name)
                    n_compartments = 2
                elseif occursin("3cpt", model_name)
                    n_compartments = 3
                end
                if occursin("tlag", model_name)
                    has_lag = true
                end
                if occursin("bio", model_name)
                    has_bio = true
                end
                if occursin("mm", model_name)
                    elimination = "mm"
                end
            end
        end
    end

    return MonolixStructuralModel(
        model_type,
        admin_type,
        n_compartments,
        elimination,
        absorption,
        has_lag,
        has_bio
    )
end

"""
Helper to infer model properties from model name.
"""
function _infer_model_properties_from_name!(
    model_name::String,
    admin_type::String,
    n_compartments::Int,
    elimination::String,
    absorption::String,
    has_lag::Bool,
    has_bio::Bool
)
    model_lower = lowercase(model_name)
    # Properties are updated via the calling function since Julia doesn't have pass-by-reference for primitives
end

"""
Parse <PARAMETER> section and <MODEL>[INDIVIDUAL] for IIV.
"""
function _parse_parameters(blocks::Dict{String,Any})::Vector{MonolixParameter}
    parameters = MonolixParameter[]
    param_values = Dict{String,Float64}()
    param_fixed = Dict{String,Bool}()
    omega_values = Dict{String,Float64}()
    distributions = Dict{String,String}()
    typical_mapping = Dict{String,String}()  # Maps param name to typical value param

    # Parse [INDIVIDUAL] block for IIV definitions
    if haskey(blocks, "MODEL")
        content = blocks["MODEL"]
        if content isa Dict && haskey(content, "INDIVIDUAL")
            in_definition = false
            for line in content["INDIVIDUAL"]
                if isempty(line) || startswith(line, ";")
                    continue
                end

                # Track DEFINITION: section
                if occursin("DEFINITION:", uppercase(line))
                    in_definition = true
                    continue
                end

                if in_definition
                    # Parse: param = {distribution=logNormal, typical=param_pop, sd=omega_param}
                    m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
                    if m !== nothing
                        param_name = String(m.captures[1])
                        props = String(m.captures[2])

                        # Extract distribution
                        dist_match = match(r"distribution\s*=\s*(\w+)", props)
                        if dist_match !== nothing
                            distributions[param_name] = String(dist_match.captures[1])
                        end

                        # Extract typical value reference
                        typical_match = match(r"typical\s*=\s*(\w+)", props)
                        if typical_match !== nothing
                            typical_mapping[param_name] = String(typical_match.captures[1])
                        end

                        # Extract sd (omega)
                        sd_match = match(r"sd\s*=\s*(\w+)", props)
                        if sd_match !== nothing
                            omega_ref = String(sd_match.captures[1])
                            # Store reference, will resolve later
                            omega_values[param_name] = -1.0  # Placeholder, will be resolved
                            # Store the omega parameter name
                            typical_mapping[param_name * "_omega"] = omega_ref
                        end

                        # Check for no-variability
                        if occursin("no-variability", lowercase(props))
                            omega_values[param_name] = 0.0
                        end
                    end
                end
            end
        end
    end

    # Parse <PARAMETER> section
    if haskey(blocks, "PARAMETER")
        content = blocks["PARAMETER"]
        lines = String[]
        if content isa Dict
            if haskey(content, "_ROOT")
                lines = content["_ROOT"]
            end
        end

        for line in lines
            if isempty(line) || startswith(line, ";")
                continue
            end

            # Parse: param = {value=..., method=...}
            m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
            if m !== nothing
                param_name = String(m.captures[1])
                props = String(m.captures[2])

                # Extract value
                val_match = match(r"value\s*=\s*([0-9.eE+-]+)", props)
                if val_match !== nothing
                    param_values[param_name] = parse(Float64, val_match.captures[1])
                end

                # Check if fixed
                param_fixed[param_name] = occursin("method=FIXED", uppercase(props))
                continue
            end

            # Parse simple: param = value
            m = match(r"(\w+)\s*=\s*([0-9.eE+-]+)$", line)
            if m !== nothing
                param_name = String(m.captures[1])
                param_values[param_name] = parse(Float64, m.captures[2])
            end
        end
    end

    # Resolve omega values from parameter values
    for (param_name, omega_ref) in typical_mapping
        if endswith(param_name, "_omega")
            base_param = replace(param_name, "_omega" => "")
            if haskey(param_values, omega_ref)
                omega_values[base_param] = param_values[omega_ref]
            end
        end
    end

    # Build parameter list - process individual parameters (with IIV)
    processed = Set{String}()
    for (param_name, typical_ref) in typical_mapping
        if endswith(param_name, "_omega")
            continue
        end

        # Get typical value from referenced parameter
        value = get(param_values, typical_ref, 1.0)
        omega = get(omega_values, param_name, 0.0)
        dist = get(distributions, param_name, "logNormal")
        has_iiv = omega > 0.0

        push!(parameters, MonolixParameter(
            param_name, value;
            distribution=dist,
            omega=omega,
            has_iiv=has_iiv
        ))
        push!(processed, param_name)
        push!(processed, typical_ref)
        if haskey(typical_mapping, param_name * "_omega")
            push!(processed, typical_mapping[param_name * "_omega"])
        end
    end

    # Add remaining parameters (error params, etc.)
    for (param_name, value) in param_values
        if param_name in processed
            continue
        end
        # Skip omega parameters (they're processed above)
        if startswith(param_name, "omega_")
            continue
        end

        fixed = get(param_fixed, param_name, false)
        push!(parameters, MonolixParameter(
            param_name, value;
            fixed=fixed,
            distribution="normal",
            omega=0.0,
            has_iiv=false
        ))
    end

    # Legacy format: Parse POPULATION block for IIV
    if haskey(blocks, "POPULATION")
        content = blocks["POPULATION"]
        lines = content isa Dict && haskey(content, "_ROOT") ? content["_ROOT"] : String[]

        for line in lines
            if isempty(line) || startswith(line, ";")
                continue
            end

            # Parse: param = {value=..., omega=..., distribution=...}
            m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
            if m !== nothing
                param_name = String(m.captures[1])
                props = String(m.captures[2])

                # Extract omega
                omega_match = match(r"omega\s*=\s*([0-9.eE+-]+)", props)
                if omega_match !== nothing
                    omega_val = parse(Float64, omega_match.captures[1])

                    # Update existing parameter or create new one
                    idx = findfirst(p -> p.name == param_name, parameters)
                    if idx !== nothing
                        # Replace with updated parameter including IIV
                        old_p = parameters[idx]
                        parameters[idx] = MonolixParameter(
                            old_p.name, old_p.value;
                            fixed=old_p.fixed,
                            distribution=get(distributions, param_name, "logNormal"),
                            omega=omega_val,
                            has_iiv=omega_val > 0.0
                        )
                    end
                end

                # Extract distribution
                dist_match = match(r"distribution\s*=\s*(\w+)", props)
                if dist_match !== nothing
                    distributions[param_name] = String(dist_match.captures[1])
                end
            end
        end
    end

    return parameters
end

"""
Extract raw parameter values from blocks for error param resolution.
"""
function _extract_parameter_values(blocks::Dict{String,Any})::Dict{String,Float64}
    param_values = Dict{String,Float64}()

    # Parse <PARAMETER> section
    if haskey(blocks, "PARAMETER")
        content = blocks["PARAMETER"]
        lines = String[]
        if content isa Dict
            if haskey(content, "_ROOT")
                lines = content["_ROOT"]
            end
        end

        for line in lines
            if isempty(line) || startswith(line, ";")
                continue
            end

            # Parse: param = {value=..., method=...}
            m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
            if m !== nothing
                param_name = String(m.captures[1])
                props = String(m.captures[2])

                # Extract value
                val_match = match(r"value\s*=\s*([0-9.eE+-]+)", props)
                if val_match !== nothing
                    param_values[param_name] = parse(Float64, val_match.captures[1])
                end
                continue
            end

            # Parse simple: param = value
            m = match(r"(\w+)\s*=\s*([0-9.eE+-]+)$", line)
            if m !== nothing
                param_name = String(m.captures[1])
                param_values[param_name] = parse(Float64, m.captures[2])
            end
        end
    end

    return param_values
end

"""
Parse observation model from <MODEL>[LONGITUDINAL].
Parameter values dict is used to resolve error parameter names to values.
"""
function _parse_observations(blocks::Dict{String,Any}, param_values::Dict{String,Float64}=Dict{String,Float64}())::Vector{MonolixObservation}
    observations = MonolixObservation[]

    # Parse from MODEL[LONGITUDINAL]
    if haskey(blocks, "MODEL")
        content = blocks["MODEL"]
        if content isa Dict && haskey(content, "LONGITUDINAL")
            in_definition = false
            for line in content["LONGITUDINAL"]
                if isempty(line) || startswith(line, ";")
                    continue
                end

                if occursin("DEFINITION:", uppercase(line))
                    in_definition = true
                    continue
                end

                if in_definition
                    # Parse: Cc = {distribution=normal, prediction=Cc, errorModel=combined2(a, b)}
                    m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
                    if m !== nothing
                        obs_name = String(m.captures[1])
                        props = String(m.captures[2])

                        obs_type = "continuous"
                        error_model = "combined"
                        error_params = Float64[]

                        # Extract type
                        type_match = match(r"type\s*=\s*(\w+)", props)
                        if type_match !== nothing
                            obs_type = lowercase(String(type_match.captures[1]))
                        end

                        # Extract error model: errorModel=combined2(a, b)
                        error_match = match(r"errorModel\s*=\s*(\w+)\s*\(([^)]*)\)", props)
                        if error_match !== nothing
                            error_str = lowercase(String(error_match.captures[1]))
                            param_str = String(error_match.captures[2])

                            # Keep original error model name for tests
                            if occursin("constant", error_str)
                                error_model = "constant"
                            elseif occursin("additive", error_str)
                                error_model = "additive"
                            elseif occursin("proportional", error_str)
                                error_model = "proportional"
                            elseif occursin("combined", error_str)
                                error_model = "combined"
                            elseif occursin("exponential", error_str)
                                error_model = "exponential"
                            end

                            # Extract parameter names and resolve values
                            for param_name in split(param_str, ',')
                                param_name = strip(String(param_name))
                                if !isempty(param_name) && haskey(param_values, param_name)
                                    push!(error_params, param_values[param_name])
                                end
                            end
                        end

                        push!(observations, MonolixObservation(
                            obs_name;
                            type=obs_type,
                            error_model=error_model,
                            error_params=error_params
                        ))
                    end
                end
            end
        end
    end

    # Legacy format: OBSERVATION_MODEL
    if isempty(observations) && haskey(blocks, "OBSERVATION_MODEL")
        content = blocks["OBSERVATION_MODEL"]
        lines = content isa Dict && haskey(content, "_ROOT") ? content["_ROOT"] : String[]

        for line in lines
            if isempty(line)
                continue
            end

            m = match(r"(\w+)\s*=\s*\{([^}]+)\}", line)
            if m !== nothing
                obs_name = String(m.captures[1])
                props = String(m.captures[2])

                obs_type = "continuous"
                error_model = "combined"

                type_match = match(r"type\s*=\s*(\w+)", props)
                if type_match !== nothing
                    obs_type = lowercase(String(type_match.captures[1]))
                end

                error_match = match(r"error\s*=\s*(\w+)", props)
                if error_match !== nothing
                    error_str = lowercase(String(error_match.captures[1]))
                    if occursin("constant", error_str)
                        error_model = "constant"
                    elseif occursin("additive", error_str)
                        error_model = "additive"
                    elseif occursin("proportional", error_str)
                        error_model = "proportional"
                    elseif occursin("combined", error_str)
                        error_model = "combined"
                    elseif occursin("exponential", error_str)
                        error_model = "exponential"
                    end
                end

                push!(observations, MonolixObservation(
                    obs_name;
                    type=obs_type,
                    error_model=error_model,
                    error_params=Float64[]
                ))
            end
        end
    end

    return observations
end

"""
Parse estimation method from <FIT> or <MONOLIX> sections.
"""
function _parse_estimation(blocks::Dict{String,Any})::String
    # Check MONOLIX[TASKS]
    if haskey(blocks, "MONOLIX")
        content = blocks["MONOLIX"]
        if content isa Dict && haskey(content, "TASKS")
            text = join(content["TASKS"], " ")
            if occursin("populationParameters", text)
                return "SAEM"
            end
        end
    end

    # Check FIT block
    if haskey(blocks, "FIT")
        content = blocks["FIT"]
        lines = content isa Dict && haskey(content, "_ROOT") ? content["_ROOT"] : String[]
        text = uppercase(join(lines, " "))

        if occursin("SAEM", text)
            return "SAEM"
        elseif occursin("FOCE", text)
            return "FOCE"
        elseif occursin("LAPLACE", text)
            return "LAPLACIAN"
        end
    end

    return "SAEM"  # Default Monolix method
end
