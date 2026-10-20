# JSON Schema Validation for FDA 21 CFR Part 11 Compliance
# Provides proactive validation of artifacts against formal schemas
#
# FDA 21 CFR 11.10(h): "...procedures... validation of data input"

using JSON

export SchemaValidationResult, SchemaValidationError
export validate_artifact_schema, get_schema_definition
export SCHEMA_REGISTRY

"""
    SchemaValidationError

A single validation error in an artifact.

# Fields
- `path::String`: JSON path to the invalid field (e.g., "model_spec.params.CL")
- `error_type::String`: Type of error (missing_field, invalid_type, out_of_range, etc.)
- `message::String`: Human-readable error message
- `expected::Any`: Expected value/type
- `actual::Any`: Actual value found
"""
struct SchemaValidationError
    path::String
    error_type::String
    message::String
    expected::Any
    actual::Any
end

"""
    SchemaValidationResult

Result of validating an artifact against its schema.

# Fields
- `is_valid::Bool`: Whether the artifact passes validation
- `schema_version::String`: Schema version validated against
- `artifact_type::String`: Type of artifact (pk_single, population, gsa, etc.)
- `errors::Vector{SchemaValidationError}`: All validation errors
- `warnings::Vector{String}`: Non-fatal validation warnings
- `validated_at::String`: ISO 8601 UTC timestamp
"""
struct SchemaValidationResult
    is_valid::Bool
    schema_version::String
    artifact_type::String
    errors::Vector{SchemaValidationError}
    warnings::Vector{String}
    validated_at::String
end

# =============================================================================
# Schema Definitions
# =============================================================================

"""
Schema registry containing validation rules for each artifact type.
"""
const SCHEMA_REGISTRY = Dict{String,Dict{String,Any}}(
    # -------------------------------------------------------------------------
    # Common fields required in all artifacts
    # -------------------------------------------------------------------------
    "common" => Dict{String,Any}(
        "required_fields" => [
            "artifact_schema_version",
            "semantics_fingerprint"
        ],
        "field_types" => Dict(
            "artifact_schema_version" => String,
            "execution_mode" => String,
            "semantics_fingerprint" => Dict
        )
    ),

    # -------------------------------------------------------------------------
    # PK Single Subject Artifact
    # -------------------------------------------------------------------------
    "pk_single" => Dict{String,Any}(
        "artifact_type" => "pk_single",
        "required_fields" => [
            "model_spec",
            "grid",
            "solver",
            "result"
        ],
        "model_spec" => Dict{String,Any}(
            "required_fields" => ["kind", "name", "params", "doses"],
            "field_types" => Dict(
                "kind" => String,
                "name" => String,
                "params" => Dict,
                "doses" => Vector
            ),
            "valid_kinds" => [
                "OneCompIVBolus",
                "OneCompOralFirstOrder",
                "TwoCompIVBolus",
                "TwoCompOral",
                "ThreeCompIVBolus",
                "TransitAbsorption",
                "MichaelisMentenElimination"
            ]
        ),
        "grid" => Dict{String,Any}(
            "required_fields" => ["t0", "t1", "saveat"],
            "field_types" => Dict(
                "t0" => Number,
                "t1" => Number,
                "saveat" => Vector
            ),
            "constraints" => Dict(
                "t0" => (min=0.0, max=Inf),
                "t1" => (min=0.0, max=Inf)
            )
        ),
        "solver" => Dict{String,Any}(
            "required_fields" => ["alg", "reltol", "abstol", "maxiters"],
            "field_types" => Dict(
                "alg" => String,
                "reltol" => Number,
                "abstol" => Number,
                "maxiters" => Integer
            ),
            "valid_algorithms" => [
                "Tsit5", "Vern7", "Vern9", "Rosenbrock23", "Rodas5P"
            ],
            "constraints" => Dict(
                "reltol" => (min=1e-14, max=1.0),
                "abstol" => (min=1e-14, max=1.0),
                "maxiters" => (min=1, max=10000000)
            )
        ),
        "result" => Dict{String,Any}(
            "required_fields" => ["t", "states", "observations", "metadata"],
            "field_types" => Dict(
                "t" => Vector,
                "states" => Dict,
                "observations" => Dict,
                "metadata" => Dict
            )
        )
    ),

    # -------------------------------------------------------------------------
    # Population Artifact
    # -------------------------------------------------------------------------
    "population" => Dict{String,Any}(
        "artifact_type" => "population",
        "required_fields" => [
            "population_spec",
            "grid",
            "solver",
            "result"
        ],
        "population_spec" => Dict{String,Any}(
            "required_fields" => ["model_kind", "name", "params", "doses", "n_subjects"],
            "field_types" => Dict(
                "n_subjects" => Integer
            ),
            "constraints" => Dict(
                "n_subjects" => (min=1, max=100000)
            )
        )
    ),

    # -------------------------------------------------------------------------
    # Sensitivity (Local) Artifact
    # -------------------------------------------------------------------------
    "sensitivity_single" => Dict{String,Any}(
        "artifact_type" => "sensitivity_single",
        "required_fields" => [
            "model_spec",
            "grid",
            "solver",
            "plan",
            "observation",
            "base_series",
            "pert_series",
            "metrics"
        ]
    ),

    # -------------------------------------------------------------------------
    # Global Sensitivity Analysis - Sobol'
    # -------------------------------------------------------------------------
    "sobol_sensitivity" => Dict{String,Any}(
        "artifact_type" => "sobol_sensitivity",
        "required_fields" => [
            "indices",
            "n_evaluations",
            "convergence_metric",
            "metadata"
        ],
        "indices" => Dict{String,Any}(
            "required_per_param" => ["Si", "STi"]
        )
    ),

    # -------------------------------------------------------------------------
    # Global Sensitivity Analysis - Morris
    # -------------------------------------------------------------------------
    "morris_sensitivity" => Dict{String,Any}(
        "artifact_type" => "morris_sensitivity",
        "required_fields" => [
            "indices",
            "elementary_effects",
            "n_evaluations",
            "metadata"
        ],
        "indices" => Dict{String,Any}(
            "required_per_param" => ["mu", "mu_star", "sigma"]
        )
    ),

    # -------------------------------------------------------------------------
    # Compliance Metadata
    # -------------------------------------------------------------------------
    "compliance_metadata" => Dict{String,Any}(
        "optional" => true,
        "audit_record" => Dict{String,Any}(
            "required_fields" => [
                "execution_id",
                "timestamp_utc",
                "system_user",
                "hostname",
                "action",
                "neopkpd_version",
                "record_hash"
            ]
        ),
        "integrity" => Dict{String,Any}(
            "required_fields" => [
                "content_hash",
                "hash_algorithm",
                "input_hash",
                "output_hash",
                "hash_computed_at"
            ]
        ),
        "environment" => Dict{String,Any}(
            "required_fields" => [
                "julia_version",
                "neopkpd_version",
                "os",
                "capture_timestamp"
            ]
        )
    )
)

# =============================================================================
# Validation Functions
# =============================================================================

"""
    _validate_field_type(value::Any, expected_type::Type) -> Bool

Check if a value matches the expected type.
"""
function _validate_field_type(value::Any, expected_type::Type)::Bool
    if expected_type == Number
        return value isa Number
    elseif expected_type == Integer
        return value isa Integer || (value isa Number && isinteger(value))
    elseif expected_type == String
        return value isa AbstractString
    elseif expected_type == Vector
        return value isa AbstractVector
    elseif expected_type == Dict
        return value isa AbstractDict
    elseif expected_type == Bool
        return value isa Bool
    else
        return value isa expected_type
    end
end

"""
    _validate_constraint(value::Number, constraint::NamedTuple) -> Bool

Check if a numeric value satisfies min/max constraints.
"""
function _validate_constraint(value::Number, constraint::NamedTuple)::Bool
    if haskey(constraint, :min) && value < constraint.min
        return false
    end
    if haskey(constraint, :max) && value > constraint.max
        return false
    end
    return true
end

"""
    _validate_section(
        data::Dict,
        schema::Dict,
        path::String,
        errors::Vector{SchemaValidationError}
    )

Recursively validate a section of the artifact.
"""
function _validate_section(
    data::Dict,
    schema::Dict,
    path::String,
    errors::Vector{SchemaValidationError}
)
    # Check required fields
    if haskey(schema, "required_fields")
        for field in schema["required_fields"]
            if !haskey(data, field)
                push!(errors, SchemaValidationError(
                    isempty(path) ? field : "$path.$field",
                    "missing_field",
                    "Required field '$field' is missing",
                    "present",
                    "missing"
                ))
            end
        end
    end

    # Check field types
    if haskey(schema, "field_types")
        for (field, expected_type) in schema["field_types"]
            if haskey(data, field)
                value = data[field]
                if !_validate_field_type(value, expected_type)
                    push!(errors, SchemaValidationError(
                        isempty(path) ? field : "$path.$field",
                        "invalid_type",
                        "Field '$field' has wrong type",
                        string(expected_type),
                        string(typeof(value))
                    ))
                end
            end
        end
    end

    # Check valid values
    if haskey(schema, "valid_kinds") && haskey(data, "kind")
        if !(data["kind"] in schema["valid_kinds"])
            push!(errors, SchemaValidationError(
                isempty(path) ? "kind" : "$path.kind",
                "invalid_value",
                "Invalid model kind '$(data["kind"])'",
                schema["valid_kinds"],
                data["kind"]
            ))
        end
    end

    if haskey(schema, "valid_algorithms") && haskey(data, "alg")
        if !(data["alg"] in schema["valid_algorithms"])
            push!(errors, SchemaValidationError(
                isempty(path) ? "alg" : "$path.alg",
                "invalid_value",
                "Invalid solver algorithm '$(data["alg"])'",
                schema["valid_algorithms"],
                data["alg"]
            ))
        end
    end

    # Check constraints
    if haskey(schema, "constraints")
        for (field, constraint) in schema["constraints"]
            if haskey(data, field) && data[field] isa Number
                if !_validate_constraint(data[field], constraint)
                    push!(errors, SchemaValidationError(
                        isempty(path) ? field : "$path.$field",
                        "out_of_range",
                        "Field '$field' value $(data[field]) is out of valid range",
                        constraint,
                        data[field]
                    ))
                end
            end
        end
    end
end

# Overload for JSON.Object to convert to Dict
function _validate_section(
    data::JSON.Object,
    schema::Dict,
    path::String,
    errors::Vector{SchemaValidationError}
)
    _validate_section(Dict{String,Any}(data), schema, path, errors)
end

"""
    _detect_artifact_type(artifact::Dict) -> String

Detect the type of artifact from its structure.
"""
function _detect_artifact_type(artifact::Dict)::String
    if haskey(artifact, "artifact_type")
        return String(artifact["artifact_type"])
    end

    # Infer from structure
    if haskey(artifact, "indices") && haskey(artifact, "elementary_effects")
        return "morris_sensitivity"
    elseif haskey(artifact, "indices") && haskey(artifact, "convergence_metric")
        return "sobol_sensitivity"
    elseif haskey(artifact, "population_spec")
        return "population"
    elseif haskey(artifact, "plan") && haskey(artifact, "observation")
        return "sensitivity_single"
    elseif haskey(artifact, "model_spec") && haskey(artifact, "result")
        return "pk_single"
    else
        return "unknown"
    end
end

"""
    validate_artifact_schema(artifact::Dict) -> SchemaValidationResult

Validate an artifact against its schema.

# Arguments
- `artifact`: The artifact dictionary to validate

# Returns
SchemaValidationResult with validation status and any errors.

# Example
```julia
result = validate_artifact_schema(artifact)
if !result.is_valid
    for error in result.errors
        println("Error at \$(error.path): \$(error.message)")
    end
end
```
"""
function validate_artifact_schema(artifact::Dict)::SchemaValidationResult
    errors = SchemaValidationError[]
    warnings = String[]

    # Get schema version
    schema_version = get(artifact, "artifact_schema_version", "unknown")

    # Detect artifact type
    artifact_type = _detect_artifact_type(artifact)

    # Validate common fields
    common_schema = SCHEMA_REGISTRY["common"]
    _validate_section(artifact, common_schema, "", errors)

    # Validate type-specific fields
    if haskey(SCHEMA_REGISTRY, artifact_type)
        type_schema = SCHEMA_REGISTRY[artifact_type]

        # Validate top-level
        _validate_section(artifact, type_schema, "", errors)

        # Validate nested sections
        for section_name in ["model_spec", "grid", "solver", "result", "population_spec"]
            if haskey(type_schema, section_name) && haskey(artifact, section_name)
                _validate_section(
                    artifact[section_name],
                    type_schema[section_name],
                    section_name,
                    errors
                )
            end
        end
    else
        push!(warnings, "Unknown artifact type: $artifact_type")
    end

    # Validate compliance_metadata if present
    if haskey(artifact, "compliance_metadata")
        meta = artifact["compliance_metadata"]
        compliance_schema = SCHEMA_REGISTRY["compliance_metadata"]

        for section in ["audit_record", "integrity", "environment"]
            if haskey(meta, section) && haskey(compliance_schema, section)
                _validate_section(
                    meta[section],
                    compliance_schema[section],
                    "compliance_metadata.$section",
                    errors
                )
            end
        end
    end

    is_valid = isempty(errors)

    return SchemaValidationResult(
        is_valid,
        string(schema_version),
        artifact_type,
        errors,
        warnings,
        _utc_timestamp()
    )
end

"""
    get_schema_definition(artifact_type::String) -> Union{Nothing, Dict}

Get the schema definition for an artifact type.
"""
function get_schema_definition(artifact_type::String)::Union{Nothing,Dict}
    return get(SCHEMA_REGISTRY, artifact_type, nothing)
end

"""
    serialize_validation_result(result::SchemaValidationResult) -> Dict{String, Any}

Serialize a SchemaValidationResult to a dictionary.
"""
function serialize_validation_result(result::SchemaValidationResult)::Dict{String,Any}
    errors_serialized = [
        Dict{String,Any}(
            "path" => e.path,
            "error_type" => e.error_type,
            "message" => e.message,
            "expected" => string(e.expected),
            "actual" => string(e.actual)
        )
        for e in result.errors
    ]

    return Dict{String,Any}(
        "is_valid" => result.is_valid,
        "schema_version" => result.schema_version,
        "artifact_type" => result.artifact_type,
        "errors" => errors_serialized,
        "warnings" => result.warnings,
        "validated_at" => result.validated_at
    )
end

export serialize_validation_result
