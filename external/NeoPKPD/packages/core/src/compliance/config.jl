# Compliance configuration for FDA 21 CFR Part 11 alignment
# Provides configurable compliance levels for different deployment scenarios

export ComplianceLevel, ComplianceConfig
export COMPLIANCE_DISABLED, COMPLIANCE_MINIMAL, COMPLIANCE_STANDARD, COMPLIANCE_STRICT
export get_compliance_config, set_compliance_config!, with_compliance

"""
    ComplianceLevel

Enumeration of compliance levels for FDA 21 CFR Part 11 alignment.

- `COMPLIANCE_DISABLED`: No compliance overhead (development/testing)
- `COMPLIANCE_MINIMAL`: Timestamps and version information only
- `COMPLIANCE_STANDARD`: Full audit trail, hashing, and environment capture
- `COMPLIANCE_STRICT`: Standard plus validation report generation
"""
@enum ComplianceLevel begin
    COMPLIANCE_DISABLED = 0
    COMPLIANCE_MINIMAL = 1
    COMPLIANCE_STANDARD = 2
    COMPLIANCE_STRICT = 3
end

"""
    ComplianceConfig

Configuration for compliance features.

# Fields
- `level::ComplianceLevel`: Overall compliance level
- `enable_audit_trail::Bool`: Whether to generate audit records with execution IDs
- `enable_hashing::Bool`: Whether to compute SHA-256 integrity hashes
- `enable_environment_capture::Bool`: Whether to capture environment snapshots
- `enable_validation_reports::Bool`: Whether to generate validation documentation
"""
struct ComplianceConfig
    level::ComplianceLevel
    enable_audit_trail::Bool
    enable_hashing::Bool
    enable_environment_capture::Bool
    enable_validation_reports::Bool
end

"""
    ComplianceConfig(level::ComplianceLevel)

Create a ComplianceConfig with settings appropriate for the given level.
"""
function ComplianceConfig(level::ComplianceLevel)
    if level == COMPLIANCE_DISABLED
        return ComplianceConfig(level, false, false, false, false)
    elseif level == COMPLIANCE_MINIMAL
        return ComplianceConfig(level, false, false, true, false)
    elseif level == COMPLIANCE_STANDARD
        return ComplianceConfig(level, true, true, true, false)
    else  # COMPLIANCE_STRICT
        return ComplianceConfig(level, true, true, true, true)
    end
end

"""
    ComplianceConfig()

Create a ComplianceConfig with default COMPLIANCE_STANDARD level.
"""
ComplianceConfig() = ComplianceConfig(COMPLIANCE_STANDARD)

# Global compliance configuration (thread-local for safety)
const _COMPLIANCE_CONFIG = Ref{ComplianceConfig}(ComplianceConfig(COMPLIANCE_STANDARD))

"""
    get_compliance_config() -> ComplianceConfig

Get the current global compliance configuration.
"""
function get_compliance_config()::ComplianceConfig
    return _COMPLIANCE_CONFIG[]
end

"""
    set_compliance_config!(config::ComplianceConfig)

Set the global compliance configuration.
"""
function set_compliance_config!(config::ComplianceConfig)
    _COMPLIANCE_CONFIG[] = config
    return config
end

"""
    set_compliance_config!(level::ComplianceLevel)

Set the global compliance configuration to the specified level.
"""
function set_compliance_config!(level::ComplianceLevel)
    return set_compliance_config!(ComplianceConfig(level))
end

"""
    with_compliance(f, level::ComplianceLevel)
    with_compliance(f, config::ComplianceConfig)

Execute function `f` with a temporary compliance configuration.
The original configuration is restored after `f` completes.

# Example
```julia
with_compliance(COMPLIANCE_DISABLED) do
    # Fast simulation without compliance overhead
    result = simulate(spec, grid, solver)
end
```
"""
function with_compliance(f::Function, config::ComplianceConfig)
    old_config = get_compliance_config()
    try
        set_compliance_config!(config)
        return f()
    finally
        set_compliance_config!(old_config)
    end
end

function with_compliance(f::Function, level::ComplianceLevel)
    return with_compliance(f, ComplianceConfig(level))
end

# Helper functions for checking compliance features
"""
    is_audit_trail_enabled() -> Bool

Check if audit trail generation is enabled in current configuration.
"""
is_audit_trail_enabled() = get_compliance_config().enable_audit_trail

"""
    is_hashing_enabled() -> Bool

Check if SHA-256 integrity hashing is enabled in current configuration.
"""
is_hashing_enabled() = get_compliance_config().enable_hashing

"""
    is_environment_capture_enabled() -> Bool

Check if environment snapshot capture is enabled in current configuration.
"""
is_environment_capture_enabled() = get_compliance_config().enable_environment_capture

"""
    is_validation_reports_enabled() -> Bool

Check if validation report generation is enabled in current configuration.
"""
is_validation_reports_enabled() = get_compliance_config().enable_validation_reports

export is_audit_trail_enabled, is_hashing_enabled
export is_environment_capture_enabled, is_validation_reports_enabled
