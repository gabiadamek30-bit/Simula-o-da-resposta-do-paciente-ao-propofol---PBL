export semantics_fingerprint

"""
Return a stable fingerprint representing numerical semantics.
"""
function semantics_fingerprint()
    return Dict(
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
        "artifact_schema_version" => ARTIFACT_SCHEMA_VERSION,
    )
end
