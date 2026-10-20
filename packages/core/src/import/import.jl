# Model Import Module
# Import NONMEM and Monolix models to NeoPKPD format
#
# This module provides parsers for external pharmacometrics tools:
# - NONMEM: Parse .ctl control files
# - Monolix: Parse .mlxtran project files

# NONMEM import
include("nonmem_types.jl")
include("nonmem_parser.jl")
include("nonmem_converter.jl")

# Monolix import
include("monolix_types.jl")
include("monolix_parser.jl")
include("monolix_converter.jl")

# ============================================================================
# Convenience Functions
# ============================================================================

export import_nonmem, import_monolix, ImportResult

"""
Result of importing an external model file.

Fields:
- model_type: The NeoPKPD model type symbol (e.g., :OneCompOralFirstOrder)
- parameters: Dictionary of parameter name => value
- spec: Full ModelSpec ready for simulation
- iiv: Inter-individual variability specification (if present)
- error: Residual error specification (if present)
- source: Original source format ("NONMEM" or "Monolix")
- warnings: Any warnings generated during import
"""
struct ImportResult
    model_type::Symbol
    parameters::Dict{Symbol,Float64}
    spec::ModelSpec
    iiv::Union{Nothing,IIVSpec}
    error::Union{Nothing,ResidualErrorSpec}
    source::String
    warnings::Vector{String}
end

"""
    import_nonmem(filepath::AbstractString; doses::Vector{DoseEvent}=DoseEvent[], name::String="") -> ImportResult

Import a NONMEM control file and convert to NeoPKPD format.

# Arguments
- `filepath`: Path to the NONMEM .ctl control file
- `doses`: Dose events for simulation (required since NONMEM data is external)
- `name`: Model name (optional, defaults to problem description)

# Returns
- `ImportResult` containing the converted model specification

# Example
```julia
result = import_nonmem("run001.ctl", doses=[DoseEvent(0.0, 100.0)])
println(result.model_type)  # :OneCompOralFirstOrder
println(result.parameters)  # Dict(:Ka => 1.5, :CL => 10.0, :V => 50.0)
sim = simulate(result.spec, SimGrid(0.0, 24.0, 0.5), SolverSpec())
```
"""
function import_nonmem(
    filepath::AbstractString;
    doses::Vector{DoseEvent}=DoseEvent[],
    name::String=""
)::ImportResult
    # Read and parse the control file
    ctl = read_nonmem_control(filepath)

    # Convert to NeoPKPD
    conversion = convert_nonmem_to_neopkpd(ctl; doses=doses, name=name)

    # Check for errors
    if !isempty(conversion.errors)
        error("NONMEM import failed: $(join(conversion.errors, "; "))")
    end

    if conversion.model_spec === nothing
        error("NONMEM import failed: no model spec generated")
    end

    # Extract model type and parameters
    model_spec = conversion.model_spec
    model_type = Symbol(typeof(model_spec.kind))

    # Build parameters dictionary from the params struct
    params = model_spec.params
    param_dict = Dict{Symbol,Float64}()
    for field in fieldnames(typeof(params))
        param_dict[field] = getfield(params, field)
    end

    return ImportResult(
        model_type,
        param_dict,
        model_spec,
        conversion.iiv_spec,
        conversion.error_spec,
        "NONMEM",
        conversion.warnings
    )
end

"""
    import_monolix(filepath::AbstractString; doses::Vector{DoseEvent}=DoseEvent[], name::String="") -> ImportResult

Import a Monolix project file and convert to NeoPKPD format.

# Arguments
- `filepath`: Path to the Monolix .mlxtran project file
- `doses`: Dose events for simulation (required since Monolix data is external)
- `name`: Model name (optional, defaults to project description)

# Returns
- `ImportResult` containing the converted model specification

# Example
```julia
result = import_monolix("project.mlxtran", doses=[DoseEvent(0.0, 100.0)])
println(result.model_type)  # :OneCompOralFirstOrder
sim = simulate(result.spec, SimGrid(0.0, 24.0, 0.5), SolverSpec())
```
"""
function import_monolix(
    filepath::AbstractString;
    doses::Vector{DoseEvent}=DoseEvent[],
    name::String=""
)::ImportResult
    # Read and parse the project file
    project = read_monolix_project(filepath)

    # Convert to NeoPKPD
    conversion = convert_monolix_to_neopkpd(project; doses=doses, name=name)

    # Check for errors
    if !isempty(conversion.errors)
        error("Monolix import failed: $(join(conversion.errors, "; "))")
    end

    if conversion.model_spec === nothing
        error("Monolix import failed: no model spec generated")
    end

    # Extract model type and parameters
    model_spec = conversion.model_spec
    model_type = Symbol(typeof(model_spec.kind))

    # Build parameters dictionary from the params struct
    params = model_spec.params
    param_dict = Dict{Symbol,Float64}()
    for field in fieldnames(typeof(params))
        param_dict[field] = getfield(params, field)
    end

    return ImportResult(
        model_type,
        param_dict,
        model_spec,
        conversion.iiv_spec,
        conversion.error_spec,
        "Monolix",
        conversion.warnings
    )
end
