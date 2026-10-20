# NONMEM Control File Types
# Structures for representing parsed NONMEM .ctl files

export NONMEMControlFile, THETASpec, OMEGABlock, SIGMABlock
export SubroutineSpec, DataSpec, InputColumn
export PKCovariateEffect, PKAssignment, ScalingFactor, PKBlock
export ErrorBlock, UnsupportedConstruct

"""
Specification for a single THETA parameter.

Fields:
- lower: Lower bound (can be -Inf)
- init: Initial estimate
- upper: Upper bound (can be Inf)
- fixed: Whether parameter is fixed
- name: Optional parameter name/label
"""
struct THETASpec
    lower::Float64
    init::Float64
    upper::Float64
    fixed::Bool
    name::String

    function THETASpec(lower::Float64, init::Float64, upper::Float64, fixed::Bool=false, name::String="")
        if !fixed && lower > init
            error("THETA lower bound ($lower) > init ($init)")
        end
        if !fixed && init > upper
            error("THETA init ($init) > upper bound ($upper)")
        end
        new(lower, init, upper, fixed, name)
    end
end

# Convenience constructor for simple (INIT) format
THETASpec(init::Float64) = THETASpec(-Inf, init, Inf, false, "")

"""
OMEGA block specification.

Represents variance-covariance matrix for inter-individual variability.

Fields:
- values: Variance/covariance values (diagonal or block)
- structure: :diagonal, :block, or :same
- dimension: Size of the block
- fixed: Whether block is fixed
"""
struct OMEGABlock
    values::Vector{Float64}
    structure::Symbol
    dimension::Int
    fixed::Bool

    function OMEGABlock(values::Vector{Float64}, structure::Symbol=:diagonal, fixed::Bool=false)
        dim = if structure == :diagonal
            length(values)
        else
            # For block, n*(n+1)/2 = length, solve for n
            n = Int((-1 + sqrt(1 + 8*length(values))) / 2)
            n
        end
        new(values, structure, dim, fixed)
    end
end

"""
SIGMA block specification.

Represents residual error variance.

Fields:
- values: Variance values
- structure: :diagonal or :block
- fixed: Whether fixed
"""
struct SIGMABlock
    values::Vector{Float64}
    structure::Symbol
    fixed::Bool

    function SIGMABlock(values::Vector{Float64}, structure::Symbol=:diagonal, fixed::Bool=false)
        new(values, structure, fixed)
    end
end

"""
\$SUBROUTINES specification (ADVAN/TRANS).

Fields:
- advan: ADVAN number (1-13)
- trans: TRANS number (1-6)
- other: Other subroutines specified
"""
struct SubroutineSpec
    advan::Int
    trans::Int
    other::Vector{String}

    function SubroutineSpec(advan::Int, trans::Int=1, other::Vector{String}=String[])
        advan in 1:13 || error("ADVAN must be 1-13, got $advan")
        trans in 1:6 || error("TRANS must be 1-6, got $trans")
        new(advan, trans, other)
    end
end

"""
\$DATA specification.

Fields:
- filename: Path to data file
- ignore: Characters/conditions to ignore
- accept: Conditions to accept
"""
struct DataSpec
    filename::String
    ignore::Vector{String}
    accept::Vector{String}

    DataSpec(filename::String) = new(filename, String[], String[])
    DataSpec(filename::String, ignore::Vector{String}, accept::Vector{String}) = new(filename, ignore, accept)
end

"""
\$INPUT column specification.

Fields:
- name: Column name
- drop: Whether to drop this column
- alias: Alias for standard names (ID, TIME, DV, etc.)
"""
struct InputColumn
    name::String
    drop::Bool
    alias::String

    InputColumn(name::String) = new(name, false, "")
    InputColumn(name::String, drop::Bool) = new(name, drop, "")
    InputColumn(name::String, drop::Bool, alias::String) = new(name, drop, alias)
end

# ============================================================================
# $PK Block Parsing Types
# ============================================================================

"""
Covariate effect extracted from \$PK block.

Represents covariate relationships like:
- Power: `TVCL * (WT/70)**THETA(3)`
- Linear: `TVCL * (1 + THETA(3)*(AGE-40))`
- Exponential: `TVCL * EXP(THETA(3)*(CRCL-100))`

Fields:
- covariate: Covariate name (WT, AGE, CRCL, etc.)
- theta_index: THETA index for the coefficient
- effect_type: :power, :linear, or :exponential
- reference: Reference value (70 for WT, 40 for AGE, etc.)
"""
struct PKCovariateEffect
    covariate::Symbol
    theta_index::Int
    effect_type::Symbol
    reference::Float64

    function PKCovariateEffect(covariate::Symbol, theta_index::Int, effect_type::Symbol, reference::Float64=0.0)
        effect_type in (:power, :linear, :exponential) || error("Invalid effect_type: $effect_type")
        theta_index > 0 || error("THETA index must be positive")
        new(covariate, theta_index, effect_type, reference)
    end
end

"""
Parameter assignment extracted from \$PK block.

Represents assignments like:
- `TVCL = THETA(1)` (typical value definition)
- `CL = TVCL * EXP(ETA(1))` (individual parameter with IIV)
- `CL = TVCL + ETA(1)` (additive random effect)

Fields:
- target: Parameter name (CL, V, KA, etc.)
- tv_theta: THETA index for typical value (nothing if not directly from THETA)
- eta_index: ETA index for IIV (nothing if no IIV)
- transformation: :exponential, :additive, or :none
- covariate_effects: Vector of covariate effects on this parameter
"""
struct PKAssignment
    target::Symbol
    tv_theta::Union{Int,Nothing}
    eta_index::Union{Int,Nothing}
    transformation::Symbol
    covariate_effects::Vector{PKCovariateEffect}

    function PKAssignment(
        target::Symbol,
        tv_theta::Union{Int,Nothing}=nothing,
        eta_index::Union{Int,Nothing}=nothing,
        transformation::Symbol=:none,
        covariate_effects::Vector{PKCovariateEffect}=PKCovariateEffect[]
    )
        transformation in (:exponential, :additive, :proportional, :none) || error("Invalid transformation: $transformation")
        new(target, tv_theta, eta_index, transformation, covariate_effects)
    end
end

"""
Scaling factor from \$PK block (S1 = V, S2 = V1, etc.).

Fields:
- compartment: Compartment number (1, 2, etc.)
- parameter: Parameter symbol it maps to (V, V1, etc.)
"""
struct ScalingFactor
    compartment::Int
    parameter::Symbol

    function ScalingFactor(compartment::Int, parameter::Symbol)
        compartment > 0 || error("Compartment must be positive")
        new(compartment, parameter)
    end
end

"""
Parsed \$PK block.

Fields:
- tv_definitions: Map from TV name to THETA index (e.g., :TVCL => 1)
- assignments: Parameter assignments with ETA and covariate info
- scaling: Scaling factors (S1, S2, etc.)
- raw_code: Original lines from \$PK block
- unsupported_lines: Lines containing unsupported constructs
"""
struct PKBlock
    tv_definitions::Dict{Symbol,Int}
    assignments::Vector{PKAssignment}
    scaling::Vector{ScalingFactor}
    raw_code::Vector{String}
    unsupported_lines::Vector{String}

    PKBlock() = new(Dict{Symbol,Int}(), PKAssignment[], ScalingFactor[], String[], String[])
    PKBlock(tv_defs, assignments, scaling, raw, unsupported) = new(tv_defs, assignments, scaling, raw, unsupported)
end

# ============================================================================
# $ERROR Block Parsing Types
# ============================================================================

"""
Parsed \$ERROR block.

Fields:
- error_type: Detected error model type (:proportional, :additive, :combined, :exponential, :unknown)
- theta_indices: THETA indices used in error model
- sigma_fixed_to_1: Whether SIGMA is fixed to 1 (common pattern)
- raw_code: Original lines from \$ERROR block
- unsupported_lines: Lines containing unsupported constructs
"""
struct ErrorBlock
    error_type::Symbol
    theta_indices::Vector{Int}
    sigma_fixed_to_1::Bool
    raw_code::Vector{String}
    unsupported_lines::Vector{String}

    ErrorBlock() = new(:unknown, Int[], false, String[], String[])
    ErrorBlock(error_type, theta_indices, sigma_fixed, raw, unsupported) = new(error_type, theta_indices, sigma_fixed, raw, unsupported)
end

# ============================================================================
# Unsupported Construct Detection
# ============================================================================

"""
Represents an unsupported NONMEM construct detected during parsing.

Fields:
- construct: Name of the unsupported construct (e.g., "IF statement", "ALAG")
- location: Where it was found ("\$PK", "\$ERROR", etc.)
- line: The offending line of code
- message: Human-readable error message
"""
struct UnsupportedConstruct
    construct::String
    location::String
    line::String
    message::String

    # 3-argument constructor with auto-generated message
    UnsupportedConstruct(construct, location, line) = new(
        construct, location, line,
        "$construct in $location is not supported: $line"
    )

    # 4-argument constructor with custom message
    UnsupportedConstruct(construct, location, line, message) = new(
        construct, location, line, message
    )
end

"""
Complete NONMEM control file representation.

This structure captures the essential elements needed to convert
a NONMEM model to NeoPKPD format.

Fields:
- problem: Problem description from \$PROBLEM
- data: Data file specification
- input: Input column definitions
- subroutines: ADVAN/TRANS specification
- thetas: Vector of THETA specifications
- omegas: Vector of OMEGA blocks
- sigmas: Vector of SIGMA blocks
- pk_code: Lines from \$PK block
- error_code: Lines from \$ERROR block
- estimation: Estimation method settings
- tables: Table output specifications
- raw_text: Original control file text
"""
struct NONMEMControlFile
    problem::String
    data::Union{Nothing,DataSpec}
    input::Vector{InputColumn}
    subroutines::Union{Nothing,SubroutineSpec}
    thetas::Vector{THETASpec}
    omegas::Vector{OMEGABlock}
    sigmas::Vector{SIGMABlock}
    pk_code::Vector{String}
    error_code::Vector{String}
    estimation::Dict{String,Any}
    tables::Vector{Dict{String,Any}}
    raw_text::String
end

# ADVAN/TRANS to Model Kind mapping
const ADVAN_TRANS_MAP = Dict{Tuple{Int,Int},Tuple{Symbol,Vector{Symbol}}}(
    # ADVAN1: One-compartment IV bolus
    (1, 1) => (:OneCompIVBolus, [:K, :V]),      # K parameterization
    (1, 2) => (:OneCompIVBolus, [:CL, :V]),     # CL/V parameterization

    # ADVAN2: One-compartment oral first-order
    (2, 1) => (:OneCompOralFirstOrder, [:KA, :K, :V]),
    (2, 2) => (:OneCompOralFirstOrder, [:KA, :CL, :V]),

    # ADVAN3: Two-compartment IV bolus
    (3, 1) => (:TwoCompIVBolus, [:K, :K12, :K21, :V]),     # Micro-constants
    (3, 3) => (:TwoCompIVBolus, [:CL, :V, :Q, :VSS]),      # CL, V, Q, Vss
    (3, 4) => (:TwoCompIVBolus, [:CL, :V1, :Q, :V2]),      # CL, V1, Q, V2

    # ADVAN4: Two-compartment oral first-order
    (4, 1) => (:TwoCompOral, [:KA, :K, :K23, :K32, :V]),
    (4, 3) => (:TwoCompOral, [:KA, :CL, :V, :Q, :VSS]),
    (4, 4) => (:TwoCompOral, [:KA, :CL, :V1, :Q, :V2]),

    # ADVAN11: Three-compartment IV bolus
    (11, 1) => (:ThreeCompIVBolus, [:K, :K12, :K21, :K13, :K31, :V]),
    (11, 4) => (:ThreeCompIVBolus, [:CL, :V1, :Q2, :V2, :Q3, :V3]),

    # ADVAN10: Michaelis-Menten elimination
    (10, 1) => (:MichaelisMentenElimination, [:VM, :KM, :V]),
)

"""
Get the NeoPKPD model kind and expected parameters for an ADVAN/TRANS combination.

Returns (model_kind_symbol, parameter_symbols) or nothing if not supported.
"""
function get_model_mapping(advan::Int, trans::Int)
    key = (advan, trans)
    if haskey(ADVAN_TRANS_MAP, key)
        return ADVAN_TRANS_MAP[key]
    end
    return nothing
end

export ADVAN_TRANS_MAP, get_model_mapping
