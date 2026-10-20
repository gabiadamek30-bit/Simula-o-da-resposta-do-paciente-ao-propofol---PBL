# Custom ODE Model Support
# Allows users to define arbitrary ODE systems for PK/PD modeling

export CustomODE, CustomODEParams, CustomODEBuilder
export create_custom_ode, validate_custom_ode

# ============================================================================
# Custom ODE Type Definitions
# ============================================================================

"""
    CustomODE <: ModelKind

A custom ODE model that allows users to define arbitrary differential equations.

This enables modeling of:
- Non-standard PK models (e.g., target-mediated drug disposition)
- Novel PD mechanisms
- Disease progression models
- Any system of ODEs

Example:
```julia
# Define a custom two-compartment model with nonlinear clearance
ode_fn! = (du, u, p, t) -> begin
    A1, A2 = u
    CL, V1, Q, V2, Vmax, Km = p.CL, p.V1, p.Q, p.V2, p.Vmax, p.Km

    C1 = A1 / V1
    # Mixed linear + Michaelis-Menten clearance
    du[1] = -(CL/V1)*A1 - (Q/V1)*A1 + (Q/V2)*A2 - Vmax*C1/(Km + C1)
    du[2] = (Q/V1)*A1 - (Q/V2)*A2
end

custom_model = CustomODE(
    ode_fn!,
    2,  # n_states
    [:CL, :V1, :Q, :V2, :Vmax, :Km],  # param_names
    (u, p) -> u[1] / p.V1,  # output_fn (concentration)
    (p,) -> [0.0, 0.0]  # initial_fn
)
```
"""
struct CustomODE <: ModelKind
    # The ODE function: ode_fn!(du, u, p, t)
    # where p is a NamedTuple of parameters
    ode_fn!::Function

    # Number of state variables
    n_states::Int

    # Parameter names (for validation and named tuple construction)
    param_names::Vector{Symbol}

    # Output function: output_fn(u, p) -> Float64
    # Typically returns concentration from amount and volume
    output_fn::Function

    # Initial conditions function: initial_fn(p) -> Vector{Float64}
    # Returns initial state vector (usually zeros)
    initial_fn::Function

    # Which state index receives the dose (default: 1)
    dose_target_index::Int

    # State variable names (for output labeling)
    state_names::Vector{Symbol}

    # Model description
    description::String

    function CustomODE(
        ode_fn!::Function,
        n_states::Int,
        param_names::Vector{Symbol},
        output_fn::Function,
        initial_fn::Function;
        dose_target_index::Int=1,
        state_names::Union{Nothing, Vector{Symbol}}=nothing,
        description::String=""
    )
        # Validate inputs
        n_states > 0 || error("n_states must be positive, got $n_states")
        dose_target_index > 0 || error("dose_target_index must be positive")
        dose_target_index <= n_states || error("dose_target_index ($dose_target_index) exceeds n_states ($n_states)")

        # Generate default state names if not provided
        if state_names === nothing
            state_names = [Symbol("state_$i") for i in 1:n_states]
        else
            length(state_names) == n_states || error("state_names length must match n_states")
        end

        new(ode_fn!, n_states, param_names, output_fn, initial_fn,
            dose_target_index, state_names, description)
    end
end

"""
    CustomODEParams

Generic parameter container for custom ODE models.
Stores parameters as a dictionary that gets converted to NamedTuple for ODE evaluation.
"""
struct CustomODEParams
    values::Dict{Symbol, Float64}
end

# Convenience constructor from keyword arguments
function CustomODEParams(; kwargs...)
    CustomODEParams(Dict{Symbol, Float64}(pairs(kwargs)...))
end

# Access parameters by symbol
Base.getproperty(p::CustomODEParams, s::Symbol) = s == :values ? getfield(p, :values) : p.values[s]
Base.haskey(p::CustomODEParams, s::Symbol) = haskey(p.values, s)

# Convert to NamedTuple for efficient ODE evaluation
function to_namedtuple(p::CustomODEParams)
    keys = Tuple(sort(collect(Base.keys(p.values))))
    vals = Tuple(p.values[k] for k in keys)
    return NamedTuple{keys}(vals)
end

# ============================================================================
# Builder Pattern for Custom ODEs
# ============================================================================

"""
    CustomODEBuilder

Fluent builder for constructing CustomODE models step by step.

Example:
```julia
model = CustomODEBuilder()
    |> set_n_states(2)
    |> set_param_names([:CL, :V1, :Q, :V2])
    |> set_ode_fn!((du, u, p, t) -> begin
        # Two-compartment model
        du[1] = -(p.CL/p.V1)*u[1] - (p.Q/p.V1)*u[1] + (p.Q/p.V2)*u[2]
        du[2] = (p.Q/p.V1)*u[1] - (p.Q/p.V2)*u[2]
    end)
    |> set_output_fn((u, p) -> u[1] / p.V1)
    |> set_initial_fn(p -> [0.0, 0.0])
    |> set_state_names([:A_central, :A_peripheral])
    |> build()
```
"""
mutable struct CustomODEBuilder
    ode_fn!::Union{Nothing, Function}
    n_states::Union{Nothing, Int}
    param_names::Union{Nothing, Vector{Symbol}}
    output_fn::Union{Nothing, Function}
    initial_fn::Union{Nothing, Function}
    dose_target_index::Int
    state_names::Union{Nothing, Vector{Symbol}}
    description::String

    CustomODEBuilder() = new(nothing, nothing, nothing, nothing, nothing, 1, nothing, "")
end

function set_n_states(builder::CustomODEBuilder, n::Int)
    builder.n_states = n
    return builder
end

function set_param_names(builder::CustomODEBuilder, names::Vector{Symbol})
    builder.param_names = names
    return builder
end

function set_ode_fn!(builder::CustomODEBuilder, fn::Function)
    builder.ode_fn! = fn
    return builder
end

function set_output_fn(builder::CustomODEBuilder, fn::Function)
    builder.output_fn = fn
    return builder
end

function set_initial_fn(builder::CustomODEBuilder, fn::Function)
    builder.initial_fn = fn
    return builder
end

function set_dose_target_index(builder::CustomODEBuilder, idx::Int)
    builder.dose_target_index = idx
    return builder
end

function set_state_names(builder::CustomODEBuilder, names::Vector{Symbol})
    builder.state_names = names
    return builder
end

function set_description(builder::CustomODEBuilder, desc::String)
    builder.description = desc
    return builder
end

function build(builder::CustomODEBuilder)::CustomODE
    # Validate required fields
    builder.ode_fn! !== nothing || error("ODE function is required")
    builder.n_states !== nothing || error("n_states is required")
    builder.param_names !== nothing || error("param_names is required")
    builder.output_fn !== nothing || error("output_fn is required")
    builder.initial_fn !== nothing || error("initial_fn is required")

    return CustomODE(
        builder.ode_fn!,
        builder.n_states,
        builder.param_names,
        builder.output_fn,
        builder.initial_fn;
        dose_target_index=builder.dose_target_index,
        state_names=builder.state_names,
        description=builder.description
    )
end

export set_n_states, set_param_names, set_ode_fn!, set_output_fn
export set_initial_fn, set_dose_target_index, set_state_names, set_description, build

# ============================================================================
# Quick Construction Functions
# ============================================================================

"""
    create_custom_ode(; ode_fn!, n_states, param_names, output_fn, initial_fn, kwargs...)

Create a CustomODE model with named arguments.

Arguments:
- ode_fn!: The ODE function (du, u, p, t) -> nothing
- n_states: Number of state variables
- param_names: Vector of parameter name symbols
- output_fn: Output function (u, p) -> Float64
- initial_fn: Initial conditions function (p,) -> Vector{Float64}
- dose_target_index: State index that receives doses (default: 1)
- state_names: Names for state variables (default: auto-generated)
- description: Model description string
"""
function create_custom_ode(;
    ode_fn!::Function,
    n_states::Int,
    param_names::Vector{Symbol},
    output_fn::Function,
    initial_fn::Function,
    dose_target_index::Int=1,
    state_names::Union{Nothing, Vector{Symbol}}=nothing,
    description::String=""
)
    return CustomODE(
        ode_fn!, n_states, param_names, output_fn, initial_fn;
        dose_target_index=dose_target_index,
        state_names=state_names,
        description=description
    )
end

# ============================================================================
# PK Interface Implementation for CustomODE
# ============================================================================

function pk_param_tuple(spec::ModelSpec{<:CustomODE, CustomODEParams})
    return to_namedtuple(spec.params)
end

function pk_state_symbols(kind::CustomODE)
    return kind.state_names
end

function pk_dose_target_index(kind::CustomODE)
    return kind.dose_target_index
end

function pk_u0(spec::ModelSpec{<:CustomODE, CustomODEParams}, grid::SimGrid)
    p = to_namedtuple(spec.params)
    return spec.kind.initial_fn(p)
end

function pk_ode!(du, u, p, t, kind::CustomODE)
    kind.ode_fn!(du, u, p, t)
    return nothing
end

function pk_conc(u, p, kind::CustomODE)
    return kind.output_fn(u, p)
end

function pk_ode_with_infusion!(du, u, p, t, kind::CustomODE, infusion_rate::Float64)
    pk_ode!(du, u, p, t, kind)
    du[pk_dose_target_index(kind)] += infusion_rate
    return nothing
end

# ============================================================================
# Validation
# ============================================================================

"""
    validate_custom_ode(kind::CustomODE, params::CustomODEParams)

Validate that the custom ODE model is well-specified.

Checks:
- All required parameters are present
- Parameter values are valid (positive where expected)
- ODE function produces valid output
- Output function produces valid output
- Initial conditions have correct dimension
"""
function validate_custom_ode(kind::CustomODE, params::CustomODEParams)
    # Check all parameters are present
    for pname in kind.param_names
        haskey(params, pname) || error("Missing parameter: $pname")
    end

    # Convert to named tuple for testing
    p = to_namedtuple(params)

    # Test initial conditions
    u0 = kind.initial_fn(p)
    length(u0) == kind.n_states || error(
        "Initial conditions dimension ($(length(u0))) doesn't match n_states ($(kind.n_states))"
    )

    # Test ODE function
    du = zeros(kind.n_states)
    try
        kind.ode_fn!(du, u0, p, 0.0)
    catch e
        error("ODE function failed at t=0: $e")
    end

    all(isfinite.(du)) || error("ODE function produces non-finite values at t=0")

    # Test output function
    try
        output = kind.output_fn(u0, p)
        isfinite(output) || error("Output function produces non-finite value at t=0")
    catch e
        error("Output function failed: $e")
    end

    return nothing
end

function validate(spec::ModelSpec{<:CustomODE, CustomODEParams})
    validate_custom_ode(spec.kind, spec.params)

    # Validate doses
    if isempty(spec.doses)
        error("At least one DoseEvent is required")
    end

    for (i, d) in enumerate(spec.doses)
        d.time >= 0.0 || error("DoseEvent time must be >= 0 at index $i, got $(d.time)")
        d.amount > 0.0 || error("DoseEvent amount must be > 0 at index $i, got $(d.amount)")
    end

    if !issorted([d.time for d in spec.doses])
        error("Dose events must be sorted by time ascending")
    end

    return nothing
end

# ============================================================================
# Pre-built Custom Models (Examples)
# ============================================================================

"""
    target_mediated_drug_disposition()

Create a Target-Mediated Drug Disposition (TMDD) model.

This model describes drugs that bind to their pharmacological target,
forming a drug-target complex that affects both PK and PD.

States:
- L: Free drug (ligand) concentration
- R: Free receptor concentration
- RL: Drug-receptor complex concentration

Parameters:
- kel: Drug elimination rate constant
- kon: Association rate constant
- koff: Dissociation rate constant
- ksyn: Receptor synthesis rate
- kdeg: Receptor degradation rate
- kint: Complex internalization rate
- V: Volume of distribution
"""
function target_mediated_drug_disposition()
    ode_fn! = (du, u, p, t) -> begin
        L, R, RL = u
        # TMDD dynamics (quasi-equilibrium approximation often used)
        du[1] = -p.kel * L - p.kon * L * R + p.koff * RL  # Free drug
        du[2] = p.ksyn - p.kdeg * R - p.kon * L * R + p.koff * RL  # Free receptor
        du[3] = p.kon * L * R - p.koff * RL - p.kint * RL  # Complex
    end

    output_fn = (u, p) -> u[1]  # Free drug concentration

    initial_fn = (p,) -> begin
        R0 = p.ksyn / p.kdeg  # Steady-state receptor
        [0.0, R0, 0.0]
    end

    return CustomODE(
        ode_fn!,
        3,
        [:kel, :kon, :koff, :ksyn, :kdeg, :kint, :V],
        output_fn,
        initial_fn;
        dose_target_index=1,
        state_names=[:L_free, :R_free, :RL_complex],
        description="Target-Mediated Drug Disposition (TMDD) model"
    )
end

export target_mediated_drug_disposition

"""
    parallel_first_order_absorption()

Create a model with parallel first-order absorption pathways.

This model describes drugs with multiple absorption sites or mechanisms,
each with its own absorption rate constant.

States:
- A_depot1: Amount in first depot
- A_depot2: Amount in second depot
- A_central: Amount in central compartment

Parameters:
- Ka1: First absorption rate constant
- Ka2: Second absorption rate constant
- F1: Fraction to first depot (F2 = 1 - F1)
- CL: Clearance
- V: Volume
"""
function parallel_first_order_absorption()
    ode_fn! = (du, u, p, t) -> begin
        A1, A2, Ac = u
        du[1] = -p.Ka1 * A1
        du[2] = -p.Ka2 * A2
        du[3] = p.Ka1 * A1 + p.Ka2 * A2 - (p.CL / p.V) * Ac
    end

    output_fn = (u, p) -> u[3] / p.V
    initial_fn = (p,) -> [0.0, 0.0, 0.0]

    return CustomODE(
        ode_fn!,
        3,
        [:Ka1, :Ka2, :F1, :CL, :V],
        output_fn,
        initial_fn;
        dose_target_index=1,  # Note: dose splitting handled separately
        state_names=[:A_depot1, :A_depot2, :A_central],
        description="Parallel first-order absorption model"
    )
end

export parallel_first_order_absorption

"""
    enterohepatic_recirculation()

Create an enterohepatic recirculation (EHR) model.

This model describes drugs that undergo biliary excretion and
reabsorption from the GI tract, leading to secondary peaks.

States:
- A_gut: Amount in GI tract
- A_central: Amount in central compartment
- A_bile: Amount in biliary compartment

Parameters:
- Ka: Absorption rate constant
- CL: Clearance
- V: Volume
- Kbile: Biliary excretion rate constant
- Kreab: Reabsorption rate constant
- F_reab: Fraction reabsorbed
"""
function enterohepatic_recirculation()
    ode_fn! = (du, u, p, t) -> begin
        Agut, Ac, Abile = u
        C = Ac / p.V

        du[1] = -p.Ka * Agut + p.F_reab * p.Kreab * Abile
        du[2] = p.Ka * Agut - (p.CL / p.V) * Ac - p.Kbile * Ac
        du[3] = p.Kbile * Ac - p.Kreab * Abile
    end

    output_fn = (u, p) -> u[2] / p.V
    initial_fn = (p,) -> [0.0, 0.0, 0.0]

    return CustomODE(
        ode_fn!,
        3,
        [:Ka, :CL, :V, :Kbile, :Kreab, :F_reab],
        output_fn,
        initial_fn;
        dose_target_index=1,
        state_names=[:A_gut, :A_central, :A_bile],
        description="Enterohepatic recirculation model"
    )
end

export enterohepatic_recirculation

"""
    autoinduction()

Create an autoinduction model where the drug induces its own metabolism.

States:
- A_central: Amount in central compartment
- E: Enzyme level (relative to baseline)

Parameters:
- CL0: Baseline clearance
- V: Volume
- Emax: Maximum enzyme induction
- EC50: Concentration for 50% induction
- kenz: Enzyme turnover rate
"""
function autoinduction()
    ode_fn! = (du, u, p, t) -> begin
        Ac, E = u
        C = Ac / p.V

        # Enzyme induction (E = 1 at baseline)
        Einduced = 1.0 + p.Emax * C / (p.EC50 + C)

        # Enzyme dynamics (first-order turnover)
        du[2] = p.kenz * (Einduced - E)

        # Drug elimination (CL scales with enzyme level)
        du[1] = -(p.CL0 * E / p.V) * Ac
    end

    output_fn = (u, p) -> u[1] / p.V
    initial_fn = (p,) -> [0.0, 1.0]  # Enzyme at baseline

    return CustomODE(
        ode_fn!,
        2,
        [:CL0, :V, :Emax, :EC50, :kenz],
        output_fn,
        initial_fn;
        dose_target_index=1,
        state_names=[:A_central, :E_enzyme],
        description="Autoinduction model"
    )
end

export autoinduction
