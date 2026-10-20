# Data integrity via SHA-256 hashing
# FDA 21 CFR Part 11 requires assurance that electronic records are complete, accurate, and unaltered

using SHA
using JSON

export IntegrityMetadata, compute_artifact_hash, compute_content_hash
export verify_artifact_integrity, IntegrityVerificationResult

"""
    IntegrityMetadata

Cryptographic integrity metadata for artifact validation.
Uses SHA-256 hashing per industry best practices.

# Fields
- `content_hash::String`: SHA-256 hash of the entire artifact content
- `hash_algorithm::String`: Hash algorithm used (always "SHA-256")
- `input_hash::String`: SHA-256 hash of input parameters only
- `output_hash::String`: SHA-256 hash of output results only
- `semantics_hash::String`: SHA-256 hash of semantics fingerprint
- `hash_computed_at::String`: ISO 8601 UTC timestamp when hashes were computed
"""
struct IntegrityMetadata
    content_hash::String
    hash_algorithm::String
    input_hash::String
    output_hash::String
    semantics_hash::String
    hash_computed_at::String
end

"""
    compute_content_hash(content::String) -> String

Compute SHA-256 hash of a string and return as hex string.
"""
function compute_content_hash(content::String)::String
    hash_bytes = sha256(content)
    return bytes2hex(hash_bytes)
end

"""
    compute_content_hash(d::Dict) -> String

Compute SHA-256 hash of a dictionary by serializing to canonical JSON.
Uses sorted keys for deterministic hashing.
"""
function compute_content_hash(d::Dict)::String
    # Sort keys for deterministic serialization
    canonical_json = _canonical_json(d)
    return compute_content_hash(canonical_json)
end

"""
    compute_content_hash(obj::JSON.Object) -> String

Compute SHA-256 hash of a JSON.Object by converting to Dict first.
"""
function compute_content_hash(obj::JSON.Object)::String
    return compute_content_hash(Dict{String,Any}(obj))
end

"""
    _canonical_json(d::Dict) -> String

Serialize a dictionary to canonical JSON with sorted keys.
This ensures deterministic hashing regardless of insertion order.
"""
function _canonical_json(d::Dict)::String
    io = IOBuffer()
    _write_canonical_json(io, d)
    return String(take!(io))
end

function _canonical_json(obj::JSON.Object)::String
    return _canonical_json(Dict{String,Any}(obj))
end

function _write_canonical_json(io::IO, d::Dict)
    write(io, "{")
    sorted_keys = sort(collect(keys(d)))
    for (i, k) in enumerate(sorted_keys)
        if i > 1
            write(io, ",")
        end
        _write_canonical_json(io, string(k))
        write(io, ":")
        _write_canonical_json(io, d[k])
    end
    write(io, "}")
end

function _write_canonical_json(io::IO, obj::JSON.Object)
    _write_canonical_json(io, Dict{String,Any}(obj))
end

function _write_canonical_json(io::IO, v::Vector)
    write(io, "[")
    for (i, elem) in enumerate(v)
        if i > 1
            write(io, ",")
        end
        _write_canonical_json(io, elem)
    end
    write(io, "]")
end

function _write_canonical_json(io::IO, s::AbstractString)
    write(io, "\"")
    for c in s
        if c == '"'
            write(io, "\\\"")
        elseif c == '\\'
            write(io, "\\\\")
        elseif c == '\n'
            write(io, "\\n")
        elseif c == '\r'
            write(io, "\\r")
        elseif c == '\t'
            write(io, "\\t")
        else
            write(io, c)
        end
    end
    write(io, "\"")
end

function _write_canonical_json(io::IO, n::Number)
    if isnan(n)
        write(io, "null")
    elseif isinf(n)
        write(io, n > 0 ? "1e308" : "-1e308")
    else
        write(io, string(n))
    end
end

function _write_canonical_json(io::IO, b::Bool)
    write(io, b ? "true" : "false")
end

function _write_canonical_json(io::IO, ::Nothing)
    write(io, "null")
end

function _write_canonical_json(io::IO, s::Symbol)
    _write_canonical_json(io, string(s))
end

function _write_canonical_json(io::IO, t::Tuple)
    write(io, "[")
    for (i, elem) in enumerate(t)
        if i > 1
            write(io, ",")
        end
        _write_canonical_json(io, elem)
    end
    write(io, "]")
end

# Catch-all for other types - convert to string
function _write_canonical_json(io::IO, x)
    _write_canonical_json(io, string(x))
end

"""
    compute_artifact_hash(artifact::Dict; exclude_keys::Vector{String}=String[]) -> String

Compute SHA-256 hash of an artifact dictionary, optionally excluding certain keys.
This is used to hash the entire artifact for integrity verification.
"""
function compute_artifact_hash(artifact::Dict; exclude_keys::Vector{String}=String[])::String
    # Create copy without excluded keys
    filtered = Dict{String,Any}()
    for (k, v) in artifact
        if !(string(k) in exclude_keys)
            filtered[string(k)] = v
        end
    end
    return compute_content_hash(filtered)
end

"""
    compute_input_hash(artifact::Dict) -> String

Compute SHA-256 hash of input-related fields only.
"""
function compute_input_hash(artifact::Dict)::String
    input_keys = ["model_spec", "population_spec", "grid", "solver", "plan", "bounds", "method"]
    inputs = Dict{String,Any}()
    for k in input_keys
        if haskey(artifact, k)
            inputs[k] = artifact[k]
        end
    end
    return compute_content_hash(inputs)
end

"""
    compute_output_hash(artifact::Dict) -> String

Compute SHA-256 hash of output-related fields only.
"""
function compute_output_hash(artifact::Dict)::String
    output_keys = ["result", "base_series", "pert_series", "metrics", "indices", "elementary_effects"]
    outputs = Dict{String,Any}()
    for k in output_keys
        if haskey(artifact, k)
            outputs[k] = artifact[k]
        end
    end
    return compute_content_hash(outputs)
end

"""
    compute_semantics_hash(artifact::Dict) -> String

Compute SHA-256 hash of semantics fingerprint.
"""
function compute_semantics_hash(artifact::Dict)::String
    if haskey(artifact, "semantics_fingerprint")
        return compute_content_hash(artifact["semantics_fingerprint"])
    else
        return compute_content_hash("")
    end
end

"""
    compute_integrity_metadata(artifact::Dict) -> IntegrityMetadata

Compute complete integrity metadata for an artifact.
"""
function compute_integrity_metadata(artifact::Dict)::IntegrityMetadata
    # Exclude compliance_metadata from content hash to avoid circular dependency
    content_hash = compute_artifact_hash(artifact; exclude_keys=["compliance_metadata"])
    input_hash = compute_input_hash(artifact)
    output_hash = compute_output_hash(artifact)
    semantics_hash = compute_semantics_hash(artifact)

    return IntegrityMetadata(
        content_hash,
        "SHA-256",
        input_hash,
        output_hash,
        semantics_hash,
        _utc_timestamp()
    )
end

"""
    serialize_integrity(meta::IntegrityMetadata) -> Dict{String, Any}

Serialize IntegrityMetadata to a dictionary.
"""
function serialize_integrity(meta::IntegrityMetadata)::Dict{String,Any}
    return Dict{String,Any}(
        "content_hash" => meta.content_hash,
        "hash_algorithm" => meta.hash_algorithm,
        "input_hash" => meta.input_hash,
        "output_hash" => meta.output_hash,
        "semantics_hash" => meta.semantics_hash,
        "hash_computed_at" => meta.hash_computed_at
    )
end

"""
    deserialize_integrity(d::Dict) -> IntegrityMetadata

Deserialize IntegrityMetadata from a dictionary.
"""
function deserialize_integrity(d::Dict)::IntegrityMetadata
    return IntegrityMetadata(
        String(d["content_hash"]),
        String(d["hash_algorithm"]),
        String(d["input_hash"]),
        String(d["output_hash"]),
        String(d["semantics_hash"]),
        String(d["hash_computed_at"])
    )
end

"""
    IntegrityVerificationResult

Result of verifying artifact integrity.

# Fields
- `is_valid::Bool`: Whether the artifact passes integrity verification
- `content_hash_valid::Bool`: Whether the content hash matches
- `input_hash_valid::Bool`: Whether the input hash matches
- `output_hash_valid::Bool`: Whether the output hash matches
- `semantics_hash_valid::Bool`: Whether the semantics hash matches
- `computed_content_hash::String`: The computed content hash
- `stored_content_hash::Union{Nothing,String}`: The stored content hash (if available)
- `error_message::Union{Nothing,String}`: Error message if verification failed
"""
struct IntegrityVerificationResult
    is_valid::Bool
    content_hash_valid::Bool
    input_hash_valid::Bool
    output_hash_valid::Bool
    semantics_hash_valid::Bool
    computed_content_hash::String
    stored_content_hash::Union{Nothing,String}
    error_message::Union{Nothing,String}
end

"""
    verify_artifact_integrity(artifact::Dict) -> IntegrityVerificationResult

Verify the integrity of an artifact by recomputing and comparing hashes.

For artifacts without compliance_metadata (schema < 1.1.0), returns valid=true
with a note that verification was not performed.
"""
function verify_artifact_integrity(artifact::Dict)::IntegrityVerificationResult
    # Check if artifact has compliance metadata
    if !haskey(artifact, "compliance_metadata")
        # Pre-1.1.0 artifact - no integrity metadata to verify
        return IntegrityVerificationResult(
            true, true, true, true, true,
            "", nothing,
            "No compliance_metadata found - artifact predates schema 1.1.0"
        )
    end

    compliance = artifact["compliance_metadata"]
    if !haskey(compliance, "integrity")
        return IntegrityVerificationResult(
            true, true, true, true, true,
            "", nothing,
            "No integrity metadata found in compliance_metadata"
        )
    end

    integrity = compliance["integrity"]

    # Compute current hashes
    computed_content = compute_artifact_hash(artifact; exclude_keys=["compliance_metadata"])
    computed_input = compute_input_hash(artifact)
    computed_output = compute_output_hash(artifact)
    computed_semantics = compute_semantics_hash(artifact)

    # Compare with stored hashes
    stored_content = String(integrity["content_hash"])
    stored_input = String(integrity["input_hash"])
    stored_output = String(integrity["output_hash"])
    stored_semantics = String(integrity["semantics_hash"])

    content_valid = computed_content == stored_content
    input_valid = computed_input == stored_input
    output_valid = computed_output == stored_output
    semantics_valid = computed_semantics == stored_semantics

    is_valid = content_valid && input_valid && output_valid && semantics_valid

    error_msg = nothing
    if !is_valid
        parts = String[]
        if !content_valid
            push!(parts, "content hash mismatch")
        end
        if !input_valid
            push!(parts, "input hash mismatch")
        end
        if !output_valid
            push!(parts, "output hash mismatch")
        end
        if !semantics_valid
            push!(parts, "semantics hash mismatch")
        end
        error_msg = "Integrity verification failed: " * join(parts, ", ")
    end

    return IntegrityVerificationResult(
        is_valid,
        content_valid,
        input_valid,
        output_valid,
        semantics_valid,
        computed_content,
        stored_content,
        error_msg
    )
end

export compute_integrity_metadata, serialize_integrity, deserialize_integrity
