# NONMEM to NeoPKPD Converter
# Converts parsed NONMEM control files to NeoPKPD model specifications

export convert_nonmem_to_neopkpd, NONMEMConversionResult
export validate_nonmem_conversion, ValidationResult

"""
Result of converting a NONMEM control file to NeoPKPD format.

Fields:
- model_spec: The converted ModelSpec (or nothing if conversion failed)
- iiv_spec: Inter-individual variability specification
- error_spec: Residual error specification
- warnings: Any warnings generated during conversion
- errors: Any errors that prevented conversion
- parameter_mapping: Mapping from NONMEM THETAs to NeoPKPD parameters
- pk_block: Parsed \$PK block information
- error_block: Parsed \$ERROR block information
- covariate_effects: Extracted covariate effects (for informational purposes)
"""
struct NONMEMConversionResult
    model_spec::Union{Nothing,ModelSpec}
    iiv_spec::Union{Nothing,IIVSpec}
    error_spec::Union{Nothing,ResidualErrorSpec}
    warnings::Vector{String}
    errors::Vector{String}
    parameter_mapping::Dict{Int,Symbol}
    pk_block::Union{Nothing,PKBlock}
    error_block::Union{Nothing,ErrorBlock}
    covariate_effects::Vector{PKCovariateEffect}

    # Backward compatible constructor
    function NONMEMConversionResult(
        model_spec, iiv_spec, error_spec, warnings, errors, parameter_mapping,
        pk_block=nothing, error_block=nothing, covariate_effects=PKCovariateEffect[]
    )
        new(model_spec, iiv_spec, error_spec, warnings, errors, parameter_mapping,
            pk_block, error_block, covariate_effects)
    end
end

"""
Result of validation checks on NONMEM conversion.

Fields:
- is_valid: Whether conversion can proceed
- errors: Fatal errors that prevent conversion
- warnings: Non-fatal warnings
"""
struct ValidationResult
    is_valid::Bool
    errors::Vector{String}
    warnings::Vector{String}
end

"""
Convert a NONMEM control file to NeoPKPD format.

Arguments:
- ctl: Parsed NONMEMControlFile
- doses: Vector of DoseEvent (required, as NONMEM data is external)
- name: Model name (optional, defaults to problem description)
- strict: If true, fail on any unsupported construct (default: true)

Returns:
- NONMEMConversionResult containing converted specs and diagnostics
"""
function convert_nonmem_to_neopkpd(
    ctl::NONMEMControlFile;
    doses::Vector{DoseEvent}=DoseEvent[],
    name::String="",
    strict::Bool=true
)::NONMEMConversionResult
    warnings = String[]
    errors = String[]
    parameter_mapping = Dict{Int,Symbol}()

    # Use problem description as name if not provided
    if isempty(name)
        name = isempty(ctl.problem) ? "imported_model" : ctl.problem
    end

    # =========================================================================
    # Step 1: Check for unsupported constructs FIRST
    # =========================================================================
    unsupported = check_unsupported_constructs(ctl)
    if !isempty(unsupported)
        for u in unsupported
            if strict
                push!(errors, u.message)
            else
                push!(warnings, "WARNING: " * u.message)
            end
        end

        # If strict mode and there are errors, return early
        if strict && !isempty(errors)
            return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping)
        end
    end

    # =========================================================================
    # Step 2: Parse $PK and $ERROR blocks
    # =========================================================================
    pk_block = parse_pk_block(ctl.pk_code)
    error_block = parse_error_block(ctl.error_code)

    # Report unsupported lines from parsing
    for line in pk_block.unsupported_lines
        push!(warnings, "Unsupported construct in \$PK: $line")
    end
    for line in error_block.unsupported_lines
        push!(warnings, "Unsupported construct in \$ERROR: $line")
    end

    # =========================================================================
    # Step 3: Validate the control file
    # =========================================================================
    validation = validate_nonmem_conversion(ctl, pk_block, error_block)
    append!(errors, validation.errors)
    append!(warnings, validation.warnings)

    if !validation.is_valid
        return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping,
                                       pk_block, error_block, PKCovariateEffect[])
    end

    # =========================================================================
    # Step 4: Check for subroutines and get model mapping
    # =========================================================================
    if ctl.subroutines === nothing
        push!(errors, "No \$SUBROUTINES found - cannot determine model type")
        return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping,
                                       pk_block, error_block, PKCovariateEffect[])
    end

    mapping = get_model_mapping(ctl.subroutines.advan, ctl.subroutines.trans)
    if mapping === nothing
        push!(errors, "Unsupported ADVAN$(ctl.subroutines.advan) TRANS$(ctl.subroutines.trans) combination")
        return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping,
                                       pk_block, error_block, PKCovariateEffect[])
    end

    model_kind_sym, param_symbols = mapping

    # =========================================================================
    # Step 5: Build parameter mapping from parsed $PK block
    # =========================================================================
    n_params = length(param_symbols)

    # Try to get parameter values from $PK assignments first
    theta_values = Float64[]
    for sym in param_symbols
        # Look for this parameter in PK assignments
        found = false
        for assignment in pk_block.assignments
            if assignment.target == sym && assignment.tv_theta !== nothing
                if assignment.tv_theta <= length(ctl.thetas)
                    push!(theta_values, ctl.thetas[assignment.tv_theta].init)
                    parameter_mapping[assignment.tv_theta] = sym
                    found = true
                    break
                end
            end
        end

        # Fallback: use positional THETA mapping
        if !found
            idx = length(theta_values) + 1
            if idx <= length(ctl.thetas)
                push!(theta_values, ctl.thetas[idx].init)
                parameter_mapping[idx] = sym
            else
                push!(errors, "Not enough THETAs for parameter $sym (need at least $n_params)")
                return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping,
                                               pk_block, error_block, PKCovariateEffect[])
            end
        end
    end

    # =========================================================================
    # Step 6: Create model spec
    # =========================================================================
    model_spec = _create_model_spec(model_kind_sym, theta_values, param_symbols, name, doses, warnings)

    if model_spec === nothing
        push!(errors, "Failed to create model spec for $(model_kind_sym)")
        return NONMEMConversionResult(nothing, nothing, nothing, warnings, errors, parameter_mapping,
                                       pk_block, error_block, PKCovariateEffect[])
    end

    # =========================================================================
    # Step 7: Convert IIV using parsed $PK block
    # =========================================================================
    iiv_spec = _convert_omega_to_iiv_from_pk(ctl.omegas, pk_block, param_symbols, warnings)

    # =========================================================================
    # Step 8: Convert error model using parsed $ERROR block
    # =========================================================================
    error_spec = _convert_error_from_blocks(error_block, ctl.thetas, ctl.sigmas, warnings)

    # =========================================================================
    # Step 9: Extract covariate effects
    # =========================================================================
    all_covariates = PKCovariateEffect[]
    for assignment in pk_block.assignments
        append!(all_covariates, assignment.covariate_effects)
    end

    if !isempty(all_covariates)
        cov_names = unique([c.covariate for c in all_covariates])
        push!(warnings, "Covariate effects detected on $(join(cov_names, ", ")). " *
                        "These are stored in covariate_effects but not yet applied during simulation.")
    end

    # =========================================================================
    # Step 10: Additional warnings
    # =========================================================================
    # Warn about extra THETAs
    max_used_theta = maximum(keys(parameter_mapping), init=0)
    for assignment in pk_block.assignments
        if assignment.tv_theta !== nothing
            max_used_theta = max(max_used_theta, assignment.tv_theta)
        end
        for cov in assignment.covariate_effects
            max_used_theta = max(max_used_theta, cov.theta_index)
        end
    end
    for idx in error_block.theta_indices
        max_used_theta = max(max_used_theta, idx)
    end

    # Detect unused THETAs
    unused_thetas = Int[]
    for i in 1:length(ctl.thetas)
        if i in keys(parameter_mapping) || i in error_block.theta_indices
            continue
        end
        # Check if THETA(i) is used in any assignment
        used_in_assignment = false
        for a in pk_block.assignments
            if a.tv_theta == i
                used_in_assignment = true
                break
            end
            for c in a.covariate_effects
                if c.theta_index == i
                    used_in_assignment = true
                    break
                end
            end
            if used_in_assignment
                break
            end
        end
        if !used_in_assignment
            push!(unused_thetas, i)
        end
    end
    if !isempty(unused_thetas)
        push!(warnings, "THETAs not detected in \$PK or \$ERROR: $(join(unused_thetas, ", "))")
    end

    return NONMEMConversionResult(model_spec, iiv_spec, error_spec, warnings, errors, parameter_mapping,
                                   pk_block, error_block, all_covariates)
end

"""
Validate NONMEM conversion for consistency.

Checks:
- THETA indices referenced don't exceed defined THETAs
- ETA indices referenced don't exceed OMEGA dimensions
- Scaling factors reference valid compartments
- Parameter assignments are complete
"""
function validate_nonmem_conversion(
    ctl::NONMEMControlFile,
    pk_block::PKBlock,
    error_block::ErrorBlock
)::ValidationResult
    errors = String[]
    warnings = String[]

    n_thetas = length(ctl.thetas)

    # Calculate OMEGA dimension
    omega_dim = 0
    for block in ctl.omegas
        omega_dim += block.dimension
    end

    # 1. Check THETA indices in $PK block
    for assignment in pk_block.assignments
        if assignment.tv_theta !== nothing && assignment.tv_theta > n_thetas
            push!(errors, "THETA($(assignment.tv_theta)) referenced in \$PK but only $n_thetas THETAs defined")
        end
        for cov in assignment.covariate_effects
            if cov.theta_index > n_thetas
                push!(errors, "THETA($(cov.theta_index)) referenced for covariate $(cov.covariate) but only $n_thetas THETAs defined")
            end
        end
    end

    # 2. Check ETA indices
    for assignment in pk_block.assignments
        if assignment.eta_index !== nothing && assignment.eta_index > omega_dim
            push!(errors, "ETA($(assignment.eta_index)) referenced but OMEGA only has $omega_dim dimensions")
        end
    end

    # 3. Check THETA indices in $ERROR block
    for theta_idx in error_block.theta_indices
        if theta_idx > n_thetas
            push!(errors, "THETA($theta_idx) referenced in \$ERROR but only $n_thetas THETAs defined")
        end
    end

    # 4. Check scaling factors
    if ctl.subroutines !== nothing
        n_compartments = _get_n_compartments(ctl.subroutines.advan)
        for sf in pk_block.scaling
            if sf.compartment > n_compartments
                push!(errors, "S$(sf.compartment) references compartment $( sf.compartment) but ADVAN$(ctl.subroutines.advan) only has $n_compartments compartments")
            end
        end
    end

    # 5. Warn about error model type if unknown
    if error_block.error_type == :unknown && !isempty(ctl.error_code)
        push!(warnings, "Could not determine error model type from \$ERROR block. Using SIGMA-based inference.")
    end

    # 6. Check for TV definitions without corresponding parameter assignments
    for (tv_name, theta_idx) in pk_block.tv_definitions
        # Extract parameter name from TV name (e.g., TVCL -> CL)
        param_name = Symbol(replace(String(tv_name), r"^TV"i => ""))
        has_assignment = any(a.target == param_name for a in pk_block.assignments)
        if !has_assignment
            push!(warnings, "TV definition $tv_name (THETA($theta_idx)) found but no corresponding individual parameter assignment (e.g., $param_name = $tv_name * EXP(ETA(n)))")
        end
    end

    return ValidationResult(isempty(errors), errors, warnings)
end

"""
Get number of NONMEM compartments for an ADVAN number.

NONMEM compartment numbers:
- ADVAN1: 1 (central only for IV bolus)
- ADVAN2: 2 (depot + central for oral 1-comp) - S2 scales central
- ADVAN3: 2 (central + peripheral for IV 2-comp)
- ADVAN4: 3 (depot + central + peripheral for oral 2-comp) - S2 scales central
- ADVAN10: 1 (central with MM elimination)
- ADVAN11: 3 (central + 2 peripherals for IV 3-comp)
"""
function _get_n_compartments(advan::Int)::Int
    if advan == 1
        return 1  # IV bolus 1-comp: central only
    elseif advan == 2
        return 2  # Oral 1-comp: depot + central
    elseif advan == 3
        return 2  # IV 2-comp: central + peripheral
    elseif advan == 4
        return 3  # Oral 2-comp: depot + central + peripheral
    elseif advan == 10
        return 1  # MM elimination: central only
    elseif advan == 11
        return 3  # IV 3-comp: central + 2 peripherals
    elseif advan == 12
        return 4  # Oral 3-comp: depot + central + 2 peripherals
    else
        return 10  # Default for general ADVANs
    end
end

"""
Create a ModelSpec from model kind and parameters.
"""
function _create_model_spec(
    kind_sym::Symbol,
    values::Vector{Float64},
    param_symbols::Vector{Symbol},
    name::String,
    doses::Vector{DoseEvent},
    warnings::Vector{String}
)::Union{Nothing,ModelSpec}

    # Map NONMEM parameter names to NeoPKPD
    param_map = Dict{Symbol,Float64}()
    for (i, sym) in enumerate(param_symbols)
        if i <= length(values)
            param_map[sym] = values[i]
        end
    end

    if kind_sym == :OneCompIVBolus
        CL = get(param_map, :CL, get(param_map, :K, 0.0) * get(param_map, :V, 1.0))
        V = get(param_map, :V, 1.0)
        params = OneCompIVBolusParams(CL, V)
        return ModelSpec(OneCompIVBolus(), name, params, doses)

    elseif kind_sym == :OneCompOralFirstOrder
        Ka = get(param_map, :KA, get(param_map, :Ka, 1.0))
        CL = get(param_map, :CL, get(param_map, :K, 0.0) * get(param_map, :V, 1.0))
        V = get(param_map, :V, 1.0)
        params = OneCompOralFirstOrderParams(Ka, CL, V)
        return ModelSpec(OneCompOralFirstOrder(), name, params, doses)

    elseif kind_sym == :TwoCompIVBolus
        CL = get(param_map, :CL, 1.0)
        V1 = get(param_map, :V1, get(param_map, :V, 1.0))
        Q = get(param_map, :Q, 1.0)
        V2 = get(param_map, :V2, get(param_map, :VSS, V1) - V1)
        params = TwoCompIVBolusParams(CL, V1, Q, V2)
        return ModelSpec(TwoCompIVBolus(), name, params, doses)

    elseif kind_sym == :TwoCompOral
        Ka = get(param_map, :KA, get(param_map, :Ka, 1.0))
        CL = get(param_map, :CL, 1.0)
        V1 = get(param_map, :V1, get(param_map, :V, 1.0))
        Q = get(param_map, :Q, 1.0)
        V2 = get(param_map, :V2, get(param_map, :VSS, V1) - V1)
        params = TwoCompOralParams(Ka, CL, V1, Q, V2)
        return ModelSpec(TwoCompOral(), name, params, doses)

    elseif kind_sym == :ThreeCompIVBolus
        CL = get(param_map, :CL, 1.0)
        V1 = get(param_map, :V1, 1.0)
        Q2 = get(param_map, :Q2, 1.0)
        V2 = get(param_map, :V2, 1.0)
        Q3 = get(param_map, :Q3, 1.0)
        V3 = get(param_map, :V3, 1.0)
        params = ThreeCompIVBolusParams(CL, V1, Q2, V2, Q3, V3)
        return ModelSpec(ThreeCompIVBolus(), name, params, doses)

    elseif kind_sym == :MichaelisMentenElimination
        Vmax = get(param_map, :VM, get(param_map, :Vmax, 1.0))
        Km = get(param_map, :KM, get(param_map, :Km, 1.0))
        V = get(param_map, :V, 1.0)
        params = MichaelisMentenEliminationParams(Vmax, Km, V)
        return ModelSpec(MichaelisMentenElimination(), name, params, doses)
    end

    push!(warnings, "Model kind $(kind_sym) not yet supported for conversion")
    return nothing
end

"""
Convert OMEGA blocks to IIV specification using parsed \$PK block.

Uses ETA assignments from PKBlock to map omegas to the correct parameters.
"""
function _convert_omega_to_iiv_from_pk(
    omegas::Vector{OMEGABlock},
    pk_block::PKBlock,
    param_symbols::Vector{Symbol},
    warnings::Vector{String}
)::Union{Nothing,IIVSpec}
    if isempty(omegas)
        return nothing
    end

    # Extract diagonal variances from OMEGA blocks
    omega_variances = Float64[]
    for block in omegas
        if block.structure == :diagonal
            append!(omega_variances, block.values)
        else
            # For block structure, extract diagonal elements
            # Block is stored as lower triangular: [ω11, ω21, ω22, ω31, ω32, ω33, ...]
            dim = block.dimension
            idx = 1
            for i in 1:dim
                # Diagonal element is at position: sum(1:i) = i*(i+1)/2
                diag_idx = div(i * (i + 1), 2)
                if diag_idx <= length(block.values)
                    push!(omega_variances, block.values[diag_idx])
                end
            end
        end
    end

    if isempty(omega_variances)
        return nothing
    end

    # Build ETA to parameter mapping from PKBlock
    eta_to_param = Dict{Int,Symbol}()
    for assignment in pk_block.assignments
        if assignment.eta_index !== nothing
            eta_to_param[assignment.eta_index] = assignment.target
        end
    end

    # Map omegas to parameters
    omegas_dict = Dict{Symbol,Float64}()

    for (eta_idx, omega_var) in enumerate(omega_variances)
        param = get(eta_to_param, eta_idx, nothing)

        if param !== nothing
            # Convert variance to SD
            omegas_dict[param] = sqrt(omega_var)
        elseif eta_idx <= length(param_symbols)
            # Fallback: use positional mapping
            omegas_dict[param_symbols[eta_idx]] = sqrt(omega_var)
        end
    end

    if isempty(omegas_dict)
        return nothing
    end

    return IIVSpec(LogNormalIIV(), omegas_dict, UInt64(12345), 1)
end

"""
Convert error model using parsed \$ERROR block and THETAs.
"""
function _convert_error_from_blocks(
    error_block::ErrorBlock,
    thetas::Vector{THETASpec},
    sigmas::Vector{SIGMABlock},
    warnings::Vector{String}
)::Union{Nothing,ResidualErrorSpec}

    # If error model is parameterized via THETA (sigma fixed to 1)
    if error_block.sigma_fixed_to_1 && !isempty(error_block.theta_indices)
        if error_block.error_type == :proportional
            if !isempty(error_block.theta_indices)
                theta_idx = error_block.theta_indices[1]
                if theta_idx <= length(thetas)
                    sigma = thetas[theta_idx].init
                    return ResidualErrorSpec(
                        ProportionalError(),
                        ProportionalErrorParams(sigma),
                        :conc,
                        UInt64(12345)
                    )
                end
            end
        elseif error_block.error_type == :additive
            if !isempty(error_block.theta_indices)
                theta_idx = error_block.theta_indices[1]
                if theta_idx <= length(thetas)
                    sigma = thetas[theta_idx].init
                    return ResidualErrorSpec(
                        AdditiveError(),
                        AdditiveErrorParams(sigma),
                        :conc,
                        UInt64(12345)
                    )
                end
            end
        elseif error_block.error_type == :combined
            if length(error_block.theta_indices) >= 2
                add_idx = error_block.theta_indices[1]
                prop_idx = error_block.theta_indices[2]
                if add_idx <= length(thetas) && prop_idx <= length(thetas)
                    sigma_add = thetas[add_idx].init
                    sigma_prop = thetas[prop_idx].init
                    return ResidualErrorSpec(
                        CombinedError(),
                        CombinedErrorParams(sigma_add, sigma_prop),
                        :conc,
                        UInt64(12345)
                    )
                end
            end
        end
    end

    # Fallback: use SIGMA values directly
    return _convert_sigma_to_error(sigmas, error_block.error_type, warnings)
end

"""
Convert SIGMA to residual error specification.
"""
function _convert_sigma_to_error(
    sigmas::Vector{SIGMABlock},
    error_type::Symbol,
    warnings::Vector{String}
)::Union{Nothing,ResidualErrorSpec}
    if isempty(sigmas)
        return nothing
    end

    sigma_values = sigmas[1].values
    if isempty(sigma_values)
        return nothing
    end

    # Use detected error type to guide interpretation
    if error_type == :additive || (error_type == :unknown && length(sigma_values) == 1)
        # Check if sigma is close to 1 (likely fixed to 1 with THETA parameterization)
        if abs(sigma_values[1] - 1.0) < 1e-6
            push!(warnings, "SIGMA = 1 detected. Error may be parameterized via THETA in \$ERROR block.")
        end
        sigma = sqrt(sigma_values[1])
        return ResidualErrorSpec(
            AdditiveError(),
            AdditiveErrorParams(sigma),
            :conc,
            UInt64(12345)
        )
    elseif error_type == :proportional
        sigma = sqrt(sigma_values[1])
        return ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(sigma),
            :conc,
            UInt64(12345)
        )
    elseif error_type == :combined && length(sigma_values) >= 2
        sigma_add = sqrt(sigma_values[1])
        sigma_prop = sqrt(sigma_values[2])
        return ResidualErrorSpec(
            CombinedError(),
            CombinedErrorParams(sigma_add, sigma_prop),
            :conc,
            UInt64(12345)
        )
    else
        # Default: treat as proportional
        sigma = sqrt(sigma_values[1])
        return ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(sigma),
            :conc,
            UInt64(12345)
        )
    end
end
