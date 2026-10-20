# Monolix to NeoPKPD Converter
# Converts parsed Monolix projects to NeoPKPD model specifications

export convert_monolix_to_neopkpd, MonolixConversionResult, check_unsupported_monolix_constructs

"""
Result of converting a Monolix project file to NeoPKPD format.

Fields:
- model_spec: The converted ModelSpec (or nothing if conversion failed)
- iiv_spec: Inter-individual variability specification
- error_spec: Residual error specification
- warnings: Any warnings generated during conversion
- errors: Any errors that prevented conversion
- parameter_mapping: Mapping from Monolix parameters to NeoPKPD parameters
- unsupported: List of unsupported constructs detected
"""
struct MonolixConversionResult
    model_spec::Union{Nothing,ModelSpec}
    iiv_spec::Union{Nothing,IIVSpec}
    error_spec::Union{Nothing,ResidualErrorSpec}
    warnings::Vector{String}
    errors::Vector{String}
    parameter_mapping::Dict{String,Symbol}
    unsupported::Vector{UnsupportedMonolixConstruct}
end

# Backward compatible constructor without unsupported field
function MonolixConversionResult(
    model_spec, iiv_spec, error_spec, warnings, errors, parameter_mapping
)
    MonolixConversionResult(model_spec, iiv_spec, error_spec, warnings, errors, parameter_mapping, UnsupportedMonolixConstruct[])
end

"""
    check_unsupported_monolix_constructs(mlx::MonolixProject) -> Vector{UnsupportedMonolixConstruct}

Check a Monolix project for unsupported constructs.

Returns a vector of UnsupportedMonolixConstruct describing each unsupported feature found.
"""
function check_unsupported_monolix_constructs(mlx::MonolixProject)::Vector{UnsupportedMonolixConstruct}
    unsupported = UnsupportedMonolixConstruct[]

    # Check model type for unsupported features
    if mlx.model !== nothing && mlx.model.model_type !== nothing
        model_name = mlx.model.model_type.model

        # Check against unsupported feature patterns
        for (pattern, name, message) in UNSUPPORTED_MONOLIX_FEATURES
            if occursin(pattern, model_name)
                push!(unsupported, UnsupportedMonolixConstruct(
                    name,
                    "<MODEL>",
                    "lib:$(model_name)",
                    message
                ))
            end
        end

        # Check if model is in supported list
        is_supported = any(p -> occursin(p, model_name), SUPPORTED_MONOLIX_PATTERNS)
        if !is_supported && isempty(unsupported)
            # Check if we can at least identify it
            mapping = get_monolix_model_mapping(model_name)
            if mapping === nothing
                push!(unsupported, UnsupportedMonolixConstruct(
                    "Unknown model",
                    "<MODEL>",
                    "lib:$(model_name)",
                    "Model '$model_name' is not recognized. Supported: 1/2/3-compartment PK with bolus, oral, or infusion administration"
                ))
            end
        end
    end

    # Check for lag time (Tlag) - supported with warning
    if mlx.model !== nothing && mlx.model.has_lag
        push!(unsupported, UnsupportedMonolixConstruct(
            "Lag time (Tlag)",
            "<MODEL>",
            "Tlag parameter detected",
            "Lag time (Tlag) is not yet implemented - model will be imported without lag time"
        ))
    end

    # Check for bioavailability fraction
    if mlx.model !== nothing && mlx.model.has_bioavailability
        push!(unsupported, UnsupportedMonolixConstruct(
            "Bioavailability (F)",
            "<MODEL>",
            "Bioavailability parameter detected",
            "Bioavailability fraction is not yet implemented - model will assume F=1"
        ))
    end

    # Check observations for unsupported types
    for obs in mlx.observations
        if obs.type != "continuous"
            push!(unsupported, UnsupportedMonolixConstruct(
                "Non-continuous observation",
                "<MODEL>[LONGITUDINAL]",
                "type=$(obs.type)",
                "Only continuous observations are supported. Found: $(obs.type)"
            ))
        end
    end

    # Check raw text for advanced features
    raw_text = mlx.raw_text

    # Check for covariates (complex covariate models)
    if occursin(r"<COVARIATE>"i, raw_text) || occursin(r"\[COVARIATE\]"i, raw_text)
        push!(unsupported, UnsupportedMonolixConstruct(
            "Covariate model",
            "<COVARIATE>",
            "Covariate block detected",
            "Complex covariate models are not yet imported automatically"
        ))
    end

    # Check for custom ODE/PK equations
    if occursin(r"EQUATION:"i, raw_text) || occursin(r"ode\s*="i, raw_text)
        push!(unsupported, UnsupportedMonolixConstruct(
            "Custom ODE",
            "<MODEL>",
            "EQUATION block detected",
            "Custom differential equations are not supported - use library models"
        ))
    end

    # Check for inter-occasion variability (IOV)
    if occursin(r"occasion"i, raw_text) || occursin(r"\biov\b"i, raw_text)
        push!(unsupported, UnsupportedMonolixConstruct(
            "Inter-occasion variability",
            "<MODEL>",
            "IOV detected",
            "Inter-occasion variability is detected but not automatically imported"
        ))
    end

    # Check for correlation structures
    if occursin(r"correlation"i, raw_text) && occursin(r"block"i, raw_text)
        push!(unsupported, UnsupportedMonolixConstruct(
            "Omega correlation",
            "<MODEL>[INDIVIDUAL]",
            "Correlation block detected",
            "Correlated random effects (block omega) are not yet supported"
        ))
    end

    return unsupported
end

"""
Convert a Monolix project to NeoPKPD format.

Arguments:
- mlx: Parsed MonolixProject
- doses: Vector of DoseEvent (required, as Monolix data is external)
- name: Model name (optional, defaults to description)
- strict: If true, errors on unsupported constructs; if false, adds warnings (default: false)

Returns:
- MonolixConversionResult containing converted specs and diagnostics
"""
function convert_monolix_to_neopkpd(
    mlx::MonolixProject;
    doses::Vector{DoseEvent}=DoseEvent[],
    name::String="",
    strict::Bool=false
)::MonolixConversionResult
    warnings = String[]
    errors = String[]
    parameter_mapping = Dict{String,Symbol}()

    # Check for unsupported constructs
    unsupported = check_unsupported_monolix_constructs(mlx)

    # Handle unsupported constructs based on strict mode
    for u in unsupported
        if strict
            # In strict mode, blocking unsupported features cause errors
            if u.construct in ["Unknown model", "Custom ODE", "Non-continuous observation",
                              "Mixture model", "PD model", "Markov model"]
                push!(errors, u.message)
            else
                push!(warnings, u.message)
            end
        else
            push!(warnings, u.message)
        end
    end

    # Use description as name if not provided
    if isempty(name)
        name = isempty(mlx.description) ? "imported_monolix_model" : mlx.description
    end

    # Check for structural model
    if mlx.model === nothing
        push!(errors, "No structural model found in project")
        return MonolixConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping, unsupported)
    end

    # Determine NeoPKPD model kind
    model_kind_sym = nothing
    if mlx.model.model_type !== nothing
        model_kind_sym = get_monolix_model_mapping(mlx.model.model_type.model)
    end

    if model_kind_sym === nothing
        # Try to infer from model properties
        model_kind_sym = _infer_model_kind(mlx.model)
    end

    if model_kind_sym === nothing
        push!(errors, "Could not determine model type from Monolix project")
        return MonolixConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping, unsupported)
    end

    # Build parameter mapping from Monolix parameters
    param_map = Dict{Symbol,Float64}()
    for p in mlx.parameters
        name_lower = lowercase(p.name)
        name_upper = uppercase(p.name)

        # Map common parameter names
        if name_lower in ["ka", "kabs", "k_a"]
            param_map[:Ka] = p.value
            parameter_mapping[p.name] = :Ka
        elseif name_lower in ["cl", "clearance"]
            param_map[:CL] = p.value
            parameter_mapping[p.name] = :CL
        elseif name_lower in ["v", "v1", "vc", "vd", "volume"]
            if model_kind_sym in [:TwoCompIVBolus, :TwoCompOral, :ThreeCompIVBolus]
                param_map[:V1] = p.value
                parameter_mapping[p.name] = :V1
            else
                param_map[:V] = p.value
                parameter_mapping[p.name] = :V
            end
        elseif name_lower in ["v2", "vp", "vp1"]
            param_map[:V2] = p.value
            parameter_mapping[p.name] = :V2
        elseif name_lower in ["v3", "vp2"]
            param_map[:V3] = p.value
            parameter_mapping[p.name] = :V3
        elseif name_lower in ["q", "q1", "cld", "cld1"]
            param_map[:Q] = p.value
            parameter_mapping[p.name] = :Q
        elseif name_lower in ["q2", "cld2"]
            param_map[:Q2] = p.value
            parameter_mapping[p.name] = :Q2
        elseif name_lower in ["q3", "cld3"]
            param_map[:Q3] = p.value
            parameter_mapping[p.name] = :Q3
        elseif name_lower in ["vm", "vmax"]
            param_map[:Vmax] = p.value
            parameter_mapping[p.name] = :Vmax
        elseif name_lower in ["km", "k_m"]
            param_map[:Km] = p.value
            parameter_mapping[p.name] = :Km
        elseif name_lower == "k"
            param_map[:K] = p.value
            parameter_mapping[p.name] = :K
        end
    end

    # Create model spec
    model_spec = _create_monolix_model_spec(model_kind_sym, param_map, name, doses, warnings)

    if model_spec === nothing
        push!(errors, "Failed to create model spec for $(model_kind_sym)")
        return MonolixConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping, unsupported)
    end

    # Convert IIV to IIVSpec
    iiv_spec = _convert_monolix_iiv(mlx.parameters, warnings)

    # Convert error model to ResidualErrorSpec
    error_spec = _convert_monolix_error(mlx.observations, mlx.parameters, warnings)

    return MonolixConversionResult(model_spec, iiv_spec, error_spec, warnings, errors, parameter_mapping, unsupported)
end

"""
Infer NeoPKPD model kind from Monolix structural model properties.
"""
function _infer_model_kind(model::MonolixStructuralModel)::Union{Nothing,Symbol}
    is_oral = model.admin_type == "oral" || model.absorption == "firstorder"
    is_mm = model.elimination == "mm"

    if is_mm
        return :MichaelisMentenElimination
    elseif model.n_compartments == 1
        return is_oral ? :OneCompOralFirstOrder : :OneCompIVBolus
    elseif model.n_compartments == 2
        return is_oral ? :TwoCompOral : :TwoCompIVBolus
    elseif model.n_compartments == 3
        return :ThreeCompIVBolus
    end

    return nothing
end

"""
Create a ModelSpec from Monolix model kind and parameters.
"""
function _create_monolix_model_spec(
    kind_sym::Symbol,
    param_map::Dict{Symbol,Float64},
    name::String,
    doses::Vector{DoseEvent},
    warnings::Vector{String}
)::Union{Nothing,ModelSpec}

    if kind_sym == :OneCompIVBolus
        CL = get(param_map, :CL, get(param_map, :K, 0.0) * get(param_map, :V, 1.0))
        V = get(param_map, :V, 1.0)
        params = OneCompIVBolusParams(CL, V)
        return ModelSpec(OneCompIVBolus(), name, params, doses)

    elseif kind_sym == :OneCompOralFirstOrder
        Ka = get(param_map, :Ka, 1.0)
        CL = get(param_map, :CL, get(param_map, :K, 0.0) * get(param_map, :V, 1.0))
        V = get(param_map, :V, 1.0)
        params = OneCompOralFirstOrderParams(Ka, CL, V)
        return ModelSpec(OneCompOralFirstOrder(), name, params, doses)

    elseif kind_sym == :TwoCompIVBolus
        CL = get(param_map, :CL, 1.0)
        V1 = get(param_map, :V1, get(param_map, :V, 1.0))
        Q = get(param_map, :Q, 1.0)
        V2 = get(param_map, :V2, 1.0)
        params = TwoCompIVBolusParams(CL, V1, Q, V2)
        return ModelSpec(TwoCompIVBolus(), name, params, doses)

    elseif kind_sym == :TwoCompOral
        Ka = get(param_map, :Ka, 1.0)
        CL = get(param_map, :CL, 1.0)
        V1 = get(param_map, :V1, get(param_map, :V, 1.0))
        Q = get(param_map, :Q, 1.0)
        V2 = get(param_map, :V2, 1.0)
        params = TwoCompOralParams(Ka, CL, V1, Q, V2)
        return ModelSpec(TwoCompOral(), name, params, doses)

    elseif kind_sym == :ThreeCompIVBolus
        CL = get(param_map, :CL, 1.0)
        V1 = get(param_map, :V1, 1.0)
        Q2 = get(param_map, :Q2, get(param_map, :Q, 1.0))
        V2 = get(param_map, :V2, 1.0)
        Q3 = get(param_map, :Q3, 1.0)
        V3 = get(param_map, :V3, 1.0)
        params = ThreeCompIVBolusParams(CL, V1, Q2, V2, Q3, V3)
        return ModelSpec(ThreeCompIVBolus(), name, params, doses)

    elseif kind_sym == :MichaelisMentenElimination
        Vmax = get(param_map, :Vmax, 1.0)
        Km = get(param_map, :Km, 1.0)
        V = get(param_map, :V, 1.0)
        params = MichaelisMentenEliminationParams(Vmax, Km, V)
        return ModelSpec(MichaelisMentenElimination(), name, params, doses)
    end

    push!(warnings, "Model kind $(kind_sym) not yet supported for conversion")
    return nothing
end

"""
Convert Monolix IIV parameters to IIVSpec.
"""
function _convert_monolix_iiv(
    parameters::Vector{MonolixParameter},
    warnings::Vector{String}
)::Union{Nothing,IIVSpec}
    omegas_dict = Dict{Symbol,Float64}()

    for p in parameters
        if p.has_iiv && p.omega > 0.0
            # Map parameter name to NeoPKPD symbol
            name_lower = lowercase(p.name)
            neopkpd_sym = if name_lower in ["ka", "kabs", "k_a"]
                :Ka
            elseif name_lower in ["cl", "clearance"]
                :CL
            elseif name_lower in ["v", "v1", "vc", "vd", "volume"]
                :V
            elseif name_lower in ["v2", "vp"]
                :V2
            elseif name_lower in ["q", "cld"]
                :Q
            elseif name_lower in ["vm", "vmax"]
                :Vmax
            elseif name_lower in ["km"]
                :Km
            else
                Symbol(uppercase(p.name))
            end

            # Monolix stores omega as SD for log-normal distribution
            # If distribution is logNormal, omega is on log-scale
            if p.distribution == "logNormal"
                omegas_dict[neopkpd_sym] = p.omega  # Already SD
            else
                omegas_dict[neopkpd_sym] = p.omega
                push!(warnings, "Non-logNormal distribution $(p.distribution) for $(p.name) - may need manual adjustment")
            end
        end
    end

    if isempty(omegas_dict)
        return nothing
    end

    return IIVSpec(LogNormalIIV(), omegas_dict, UInt64(12345), 1)
end

"""
Convert Monolix observation error model to ResidualErrorSpec.
"""
function _convert_monolix_error(
    observations::Vector{MonolixObservation},
    parameters::Vector{MonolixParameter},
    warnings::Vector{String}
)::Union{Nothing,ResidualErrorSpec}
    if isempty(observations)
        return nothing
    end

    # Use first observation's error model
    obs = observations[1]
    error_model = obs.error_model

    # Use error params from observation if available, otherwise find from parameter list
    error_param_a = 0.1  # Default additive
    error_param_b = 0.1  # Default proportional

    # First try to use the extracted error_params from the observation
    if !isempty(obs.error_params)
        if error_model in ["constant", "additive"]
            error_param_a = obs.error_params[1]
        elseif error_model == "proportional"
            error_param_b = obs.error_params[1]
        elseif error_model == "combined"
            if length(obs.error_params) >= 2
                error_param_a = obs.error_params[1]
                error_param_b = obs.error_params[2]
            elseif length(obs.error_params) >= 1
                error_param_a = obs.error_params[1]
            end
        elseif error_model == "exponential"
            error_param_b = obs.error_params[1]
        end
    else
        # Fallback: Try to find error parameters from the parameter list
        for p in parameters
            name_lower = lowercase(p.name)
            if name_lower in ["a", "a1", "sigma_a", "sigma_add"]
                error_param_a = p.value
            elseif name_lower in ["b", "b1", "sigma_b", "sigma_prop"]
                error_param_b = p.value
            elseif name_lower in ["sigma", "sigma1"]
                # Could be either depending on model
                if error_model in ["additive", "constant"]
                    error_param_a = p.value
                else
                    error_param_b = p.value
                end
            end
        end
    end

    # Map error model names (constant is the same as additive)
    if error_model in ["additive", "constant"]
        return ResidualErrorSpec(
            AdditiveError(),
            AdditiveErrorParams(error_param_a),
            :conc,
            UInt64(12345)
        )
    elseif error_model == "proportional"
        return ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(error_param_b),
            :conc,
            UInt64(12345)
        )
    elseif error_model == "combined"
        return ResidualErrorSpec(
            CombinedError(),
            CombinedErrorParams(error_param_a, error_param_b),
            :conc,
            UInt64(12345)
        )
    elseif error_model == "exponential"
        return ResidualErrorSpec(
            ExponentialError(),
            ExponentialErrorParams(error_param_b),
            :conc,
            UInt64(12345)
        )
    end

    push!(warnings, "Unknown error model $(error_model) - using combined error with default parameters")
    return ResidualErrorSpec(
        CombinedError(),
        CombinedErrorParams(0.1, 0.1),
        :conc,
        UInt64(12345)
    )
end
