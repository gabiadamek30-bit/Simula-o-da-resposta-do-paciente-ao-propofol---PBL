# Monolix Project File Types
# Structures for representing parsed .mlxtran files

export MonolixProject, MonolixDataset, MonolixParameter, MonolixObservation
export MonolixModelType, MonolixStructuralModel, UnsupportedMonolixConstruct

"""
Monolix model type specification.

Fields:
- lib: Library name (e.g., "pklib")
- model: Model name (e.g., "pk_oral1cpt_1abs_TlagkaTbioVmClQ2V2_PLASMA")
"""
struct MonolixModelType
    lib::String
    model::String
end

"""
Monolix structural model configuration.

Fields:
- model_type: The model library and name
- admin_type: Administration type (e.g., "oral", "iv", "infusion")
- n_compartments: Number of compartments
- elimination: Elimination type ("linear", "mm", "mixed")
- absorption: Absorption type ("firstOrder", "zeroOrder", etc.)
"""
struct MonolixStructuralModel
    model_type::Union{Nothing,MonolixModelType}
    admin_type::String
    n_compartments::Int
    elimination::String
    absorption::String
    has_lag::Bool
    has_bioavailability::Bool
end

"""
Monolix parameter definition.

Fields:
- name: Parameter name (e.g., "ka", "V", "Cl")
- value: Initial value
- fixed: Whether parameter is fixed
- distribution: Distribution type ("logNormal", "normal", "logitNormal")
- omega: IIV variance (if has_iiv)
- has_iiv: Whether parameter has inter-individual variability
"""
struct MonolixParameter
    name::String
    value::Float64
    fixed::Bool
    distribution::String
    omega::Float64
    has_iiv::Bool

    function MonolixParameter(
        name::String,
        value::Float64;
        fixed::Bool=false,
        distribution::String="logNormal",
        omega::Float64=0.0,
        has_iiv::Bool=false
    )
        new(name, value, fixed, distribution, omega, has_iiv)
    end
end

"""
Monolix observation definition.

Fields:
- name: Observation name (e.g., "y1", "Cc")
- type: Observation type ("continuous", "discrete")
- error_model: Error model type ("constant", "proportional", "combined", "exponential")
- error_params: Error model parameters
"""
struct MonolixObservation
    name::String
    type::String
    error_model::String
    error_params::Vector{Float64}

    function MonolixObservation(
        name::String;
        type::String="continuous",
        error_model::String="combined",
        error_params::Vector{Float64}=Float64[]
    )
        new(name, type, error_model, error_params)
    end
end

"""
Monolix dataset specification.

Fields:
- filename: Path to data file
- header_types: Column type mappings
- id_column: Subject ID column
- time_column: Time column
- observation_column: Observation column
- dose_column: Dose column (if present)
- rate_column: Infusion rate column (if present)
"""
struct MonolixDataset
    filename::String
    header_types::Dict{String,String}
    id_column::String
    time_column::String
    observation_column::String
    dose_column::String
    rate_column::String

    function MonolixDataset(
        filename::String;
        header_types::Dict{String,String}=Dict{String,String}(),
        id_column::String="ID",
        time_column::String="TIME",
        observation_column::String="DV",
        dose_column::String="AMT",
        rate_column::String="RATE"
    )
        new(filename, header_types, id_column, time_column, observation_column, dose_column, rate_column)
    end
end

"""
Complete Monolix project representation.

This structure captures the essential elements needed to convert
a Monolix project to NeoPKPD format.

Fields:
- description: Project description
- data: Dataset specification
- model: Structural model configuration
- parameters: Vector of parameter definitions
- observations: Vector of observation definitions
- estimation_method: Estimation algorithm used
- raw_text: Original .mlxtran file content
"""
struct MonolixProject
    description::String
    data::Union{Nothing,MonolixDataset}
    model::Union{Nothing,MonolixStructuralModel}
    parameters::Vector{MonolixParameter}
    observations::Vector{MonolixObservation}
    estimation_method::String
    raw_text::String
end

# Mapping from Monolix model names to NeoPKPD model kinds
const MONOLIX_MODEL_MAP = Dict{String,Symbol}(
    # One-compartment models
    "pk_bolus1cpt_Vk_PLASMA" => :OneCompIVBolus,
    "pk_bolus1cpt_VCl_PLASMA" => :OneCompIVBolus,
    "pk_oral1cpt_1abs_kaTbioVk_PLASMA" => :OneCompOralFirstOrder,
    "pk_oral1cpt_1abs_kaVCl_PLASMA" => :OneCompOralFirstOrder,
    "pk_oral1cpt_1abs_TlagkaVCl_PLASMA" => :OneCompOralFirstOrder,

    # Two-compartment models
    "pk_bolus2cpt_Vk12k21k_PLASMA" => :TwoCompIVBolus,
    "pk_bolus2cpt_V1ClQ2V2_PLASMA" => :TwoCompIVBolus,
    "pk_oral2cpt_1abs_kaV1ClQ2V2_PLASMA" => :TwoCompOral,
    "pk_oral2cpt_1abs_TlagkaV1ClQ2V2_PLASMA" => :TwoCompOral,

    # Three-compartment models
    "pk_bolus3cpt_V1ClQ2V2Q3V3_PLASMA" => :ThreeCompIVBolus,

    # Michaelis-Menten models
    "pk_bolus1cpt_VVmKm_PLASMA" => :MichaelisMentenElimination,
)

"""
Get the NeoPKPD model kind for a Monolix model name.

Returns the corresponding Symbol or nothing if not supported.
"""
function get_monolix_model_mapping(model_name::String)
    # Try exact match first
    if haskey(MONOLIX_MODEL_MAP, model_name)
        return MONOLIX_MODEL_MAP[model_name]
    end

    # Try pattern matching for common variations
    model_lower = lowercase(model_name)

    if occursin("oral", model_lower) && occursin("1cpt", model_lower)
        return :OneCompOralFirstOrder
    elseif occursin("bolus", model_lower) && occursin("1cpt", model_lower)
        return :OneCompIVBolus
    elseif occursin("oral", model_lower) && occursin("2cpt", model_lower)
        return :TwoCompOral
    elseif occursin("bolus", model_lower) && occursin("2cpt", model_lower)
        return :TwoCompIVBolus
    elseif occursin("3cpt", model_lower)
        return :ThreeCompIVBolus
    elseif occursin("mm", model_lower) || occursin("michaelis", model_lower)
        return :MichaelisMentenElimination
    end

    return nothing
end

# ============================================================================
# Unsupported Construct Detection
# ============================================================================

"""
Represents an unsupported Monolix construct detected during parsing.

Fields:
- construct: Name of the unsupported construct (e.g., "Custom ODE", "Mixture model")
- location: Where it was found ("<MODEL>", "<COVARIATE>", etc.)
- line: The offending line of code
- message: Human-readable error message
"""
struct UnsupportedMonolixConstruct
    construct::String
    location::String
    line::String
    message::String

    # 3-argument constructor with auto-generated message
    UnsupportedMonolixConstruct(construct, location, line) = new(
        construct, location, line,
        "$construct in $location is not supported: $line"
    )

    # 4-argument constructor with custom message
    UnsupportedMonolixConstruct(construct, location, line, message) = new(
        construct, location, line, message
    )
end

# Supported Monolix model patterns
const SUPPORTED_MONOLIX_PATTERNS = [
    r"oral.*1cpt"i,
    r"oral.*2cpt"i,
    r"bolus.*1cpt"i,
    r"bolus.*2cpt"i,
    r"bolus.*3cpt"i,
    r"infusion.*1cpt"i,
    r"infusion.*2cpt"i,
    r"mm"i,
    r"michaelis"i,
]

# Unsupported Monolix features
const UNSUPPORTED_MONOLIX_FEATURES = [
    (r"PD"i, "PD model", "Pharmacodynamic models are not yet supported"),
    (r"turnover"i, "Turnover model", "Turnover models require custom ODE support"),
    (r"transit"i, "Transit compartment", "Transit compartment models are not yet supported"),
    (r"mixture"i, "Mixture model", "Mixture models are not supported"),
    (r"markov"i, "Markov model", "Markov chain models are not supported"),
    (r"tte"i, "Time-to-event", "Time-to-event models are not supported"),
    (r"count"i, "Count data", "Count data models are not supported"),
    (r"categorical"i, "Categorical", "Categorical response models are not supported"),
    (r"ordered"i, "Ordered categorical", "Ordered categorical models are not supported"),
]

export MONOLIX_MODEL_MAP, get_monolix_model_mapping, SUPPORTED_MONOLIX_PATTERNS, UNSUPPORTED_MONOLIX_FEATURES
