# Environment capture for FDA 21 CFR Part 11 compliance
# Captures comprehensive system and package information for reproducibility

export EnvironmentSnapshot, capture_environment, serialize_environment

"""
    EnvironmentSnapshot

Comprehensive snapshot of the execution environment for reproducibility tracking.
Aligned with ALCOA+ principles for complete and enduring records.

# Fields
- `julia_version::String`: Julia version (e.g., "1.10.0")
- `julia_build::String`: Julia build information
- `neopkpd_version::String`: NeoPKPD Core version
- `artifact_schema_version::String`: Artifact schema version
- `event_semantics_version::String`: Event semantics version
- `solver_semantics_version::String`: Solver semantics version
- `differentialequations_version::String`: DifferentialEquations.jl version
- `scimlbase_version::String`: SciMLBase.jl version
- `os::String`: Operating system name
- `os_version::String`: Operating system version
- `cpu_model::String`: CPU model name
- `cpu_cores::Int`: Number of CPU cores
- `memory_gb::Float64`: Total system memory in GB
- `architecture::String`: System architecture (x86_64, aarch64, etc.)
- `word_size::Int`: Julia word size (32 or 64)
- `python_version::Union{Nothing,String}`: Python version if available
- `git_commit::Union{Nothing,String}`: Git commit hash if in repository
- `git_branch::Union{Nothing,String}`: Git branch name if in repository
- `git_dirty::Union{Nothing,Bool}`: Whether working directory has uncommitted changes
- `capture_timestamp::String`: ISO 8601 UTC timestamp of capture
"""
struct EnvironmentSnapshot
    julia_version::String
    julia_build::String
    neopkpd_version::String
    artifact_schema_version::String
    event_semantics_version::String
    solver_semantics_version::String
    differentialequations_version::String
    scimlbase_version::String
    os::String
    os_version::String
    cpu_model::String
    cpu_cores::Int
    memory_gb::Float64
    architecture::String
    word_size::Int
    python_version::Union{Nothing,String}
    git_commit::Union{Nothing,String}
    git_branch::Union{Nothing,String}
    git_dirty::Union{Nothing,Bool}
    capture_timestamp::String
end

"""
    _get_package_version(pkg_name::String) -> String

Get the version of an installed package. Returns "unknown" if not found.
"""
function _get_package_version(pkg_name::String)::String
    try
        # Try to get version from Pkg
        deps = Pkg.dependencies()
        for (uuid, pkg_info) in deps
            if pkg_info.name == pkg_name
                return string(pkg_info.version)
            end
        end
        return "unknown"
    catch
        return "unknown"
    end
end

"""
    _get_cpu_model() -> String

Get the CPU model name. Platform-specific implementation.
"""
function _get_cpu_model()::String
    try
        if Sys.isapple()
            result = read(`sysctl -n machdep.cpu.brand_string`, String)
            return strip(result)
        elseif Sys.islinux()
            for line in eachline("/proc/cpuinfo")
                if startswith(line, "model name")
                    parts = split(line, ":")
                    if length(parts) >= 2
                        return strip(parts[2])
                    end
                end
            end
            return "unknown"
        elseif Sys.iswindows()
            result = read(`wmic cpu get name`, String)
            lines = split(result, "\n")
            if length(lines) >= 2
                return strip(lines[2])
            end
            return "unknown"
        else
            return "unknown"
        end
    catch
        return "unknown"
    end
end

"""
    _get_memory_gb() -> Float64

Get total system memory in gigabytes.
"""
function _get_memory_gb()::Float64
    try
        return Sys.total_memory() / (1024^3)
    catch
        return 0.0
    end
end

"""
    _get_os_version() -> String

Get the operating system version string.
"""
function _get_os_version()::String
    try
        if Sys.isapple()
            result = read(`sw_vers -productVersion`, String)
            return strip(result)
        elseif Sys.islinux()
            if isfile("/etc/os-release")
                for line in eachline("/etc/os-release")
                    if startswith(line, "VERSION=")
                        return replace(strip(split(line, "=")[2]), "\"" => "")
                    end
                end
            end
            return string(Sys.KERNEL)
        elseif Sys.iswindows()
            return string(Sys.windows_version())
        else
            return string(Sys.KERNEL)
        end
    catch
        return "unknown"
    end
end

"""
    _get_git_info() -> NamedTuple{(:commit, :branch, :dirty), Tuple{Union{Nothing,String}, Union{Nothing,String}, Union{Nothing,Bool}}}

Get git repository information if available.
"""
function _get_git_info()
    commit = nothing
    branch = nothing
    dirty = nothing

    try
        # Check if we're in a git repository
        run(pipeline(`git rev-parse --git-dir`, devnull, devnull))

        # Get commit hash
        commit = strip(read(`git rev-parse HEAD`, String))

        # Get branch name
        branch_result = read(`git rev-parse --abbrev-ref HEAD`, String)
        branch = strip(branch_result)

        # Check if dirty
        status = read(`git status --porcelain`, String)
        dirty = !isempty(strip(status))
    catch
        # Not in a git repository or git not available
    end

    return (commit=commit, branch=branch, dirty=dirty)
end

"""
    _get_python_version() -> Union{Nothing, String}

Get Python version if available through PyCall or system Python.
"""
function _get_python_version()::Union{Nothing,String}
    try
        result = read(`python3 --version`, String)
        parts = split(strip(result))
        if length(parts) >= 2
            return parts[2]
        end
        return strip(result)
    catch
        try
            result = read(`python --version`, String)
            parts = split(strip(result))
            if length(parts) >= 2
                return parts[2]
            end
            return strip(result)
        catch
            return nothing
        end
    end
end

"""
    _utc_timestamp() -> String

Get current UTC timestamp in ISO 8601 format with milliseconds.
"""
function _utc_timestamp()::String
    dt = Dates.now(Dates.UTC)
    return Dates.format(dt, "yyyy-mm-ddTHH:MM:SS.sss") * "Z"
end

"""
    capture_environment() -> EnvironmentSnapshot

Capture a comprehensive snapshot of the current execution environment.
This function gathers all relevant system, package, and repository information
needed for FDA 21 CFR Part 11 compliance and reproducibility tracking.

# Example
```julia
env = capture_environment()
println("Julia version: ", env.julia_version)
println("NeoPKPD version: ", env.neopkpd_version)
println("Git commit: ", env.git_commit)
```
"""
function capture_environment()::EnvironmentSnapshot
    git_info = _get_git_info()

    return EnvironmentSnapshot(
        string(VERSION),
        string(VERSION.build),
        NEOPKPD_VERSION,
        ARTIFACT_SCHEMA_VERSION,
        string(EVENT_SEMANTICS_VERSION),
        string(SOLVER_SEMANTICS_VERSION),
        _get_package_version("DifferentialEquations"),
        _get_package_version("SciMLBase"),
        string(Sys.KERNEL),
        _get_os_version(),
        _get_cpu_model(),
        Sys.CPU_THREADS,
        _get_memory_gb(),
        string(Sys.ARCH),
        Sys.WORD_SIZE,
        _get_python_version(),
        git_info.commit,
        git_info.branch,
        git_info.dirty,
        _utc_timestamp()
    )
end

"""
    serialize_environment(env::EnvironmentSnapshot) -> Dict{String, Any}

Serialize an EnvironmentSnapshot to a dictionary for JSON encoding.
"""
function serialize_environment(env::EnvironmentSnapshot)::Dict{String,Any}
    return Dict{String,Any}(
        "julia_version" => env.julia_version,
        "julia_build" => env.julia_build,
        "neopkpd_version" => env.neopkpd_version,
        "artifact_schema_version" => env.artifact_schema_version,
        "event_semantics_version" => env.event_semantics_version,
        "solver_semantics_version" => env.solver_semantics_version,
        "differentialequations_version" => env.differentialequations_version,
        "scimlbase_version" => env.scimlbase_version,
        "os" => env.os,
        "os_version" => env.os_version,
        "cpu_model" => env.cpu_model,
        "cpu_cores" => env.cpu_cores,
        "memory_gb" => env.memory_gb,
        "architecture" => env.architecture,
        "word_size" => env.word_size,
        "python_version" => env.python_version,
        "git_commit" => env.git_commit,
        "git_branch" => env.git_branch,
        "git_dirty" => env.git_dirty,
        "capture_timestamp" => env.capture_timestamp
    )
end

"""
    deserialize_environment(d::Dict) -> EnvironmentSnapshot

Deserialize an EnvironmentSnapshot from a dictionary.
"""
function deserialize_environment(d::Dict)::EnvironmentSnapshot
    return EnvironmentSnapshot(
        String(d["julia_version"]),
        String(d["julia_build"]),
        String(d["neopkpd_version"]),
        String(d["artifact_schema_version"]),
        String(d["event_semantics_version"]),
        String(d["solver_semantics_version"]),
        String(d["differentialequations_version"]),
        String(d["scimlbase_version"]),
        String(d["os"]),
        String(d["os_version"]),
        String(d["cpu_model"]),
        Int(d["cpu_cores"]),
        Float64(d["memory_gb"]),
        String(d["architecture"]),
        Int(d["word_size"]),
        d["python_version"] === nothing ? nothing : String(d["python_version"]),
        d["git_commit"] === nothing ? nothing : String(d["git_commit"]),
        d["git_branch"] === nothing ? nothing : String(d["git_branch"]),
        d["git_dirty"] === nothing ? nothing : Bool(d["git_dirty"]),
        String(d["capture_timestamp"])
    )
end

export deserialize_environment
