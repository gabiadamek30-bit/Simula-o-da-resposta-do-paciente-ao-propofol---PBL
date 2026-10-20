# =============================================================================
# TMDD Population Simulation
# =============================================================================
#
# Population-level simulation for TMDD models with IIV support.
# Industry-standard implementation compatible with NONMEM/Monolix.
# =============================================================================

using StableRNGs
using Statistics

export TMDDPopulationSpec, TMDDPopulationResult
export simulate_tmdd_population

# =============================================================================
# Population Specification
# =============================================================================

"""
    TMDDPopulationSpec{K,P}

Population specification for TMDD models.

Supports inter-individual variability (IIV) on TMDD parameters following
exponential random effects model (log-normal distribution).

# Fields
- `base_spec`: Base TMDDSpec with typical values
- `iiv`: IIVSpec with omega values for each parameter
- `n_subjects`: Number of subjects to simulate

# Example
```julia
base = TMDDSpec(
    TwoCptTMDD(QSS, IVBolus),
    "mAb PopPK",
    TwoCptTMDDParams(
        CL=0.22, V1=3.28, V2=2.66, Q=0.54,
        KSS=0.087, kint=0.037, ksyn=0.001, kdeg=0.012, R0=0.083
    ),
    [DoseEvent(0.0, 200.0)]
)

# Typical omegas for mAb PK
iiv = IIVSpec(
    Dict(:CL => 0.3, :V1 => 0.25, :V2 => 0.3, :Q => 0.3),
    100,  # n subjects
    12345  # seed
)

pop_spec = TMDDPopulationSpec(base, iiv, 100)
```
"""
struct TMDDPopulationSpec{K<:TMDDModelKind,P}
    base_spec::TMDDSpec{K,P}
    iiv::Union{Nothing,IIVSpec}
    n_subjects::Int
end

"""
Create a population spec from a base TMDD spec.
"""
function TMDDPopulationSpec(
    spec::TMDDSpec{K,P},
    iiv::IIVSpec
) where {K<:TMDDModelKind,P}
    return TMDDPopulationSpec(spec, iiv, iiv.n)
end

"""
Create population spec without IIV (for simulation at typical values).
"""
function TMDDPopulationSpec(
    spec::TMDDSpec{K,P},
    n_subjects::Int
) where {K<:TMDDModelKind,P}
    return TMDDPopulationSpec(spec, nothing, n_subjects)
end

# =============================================================================
# Population Result
# =============================================================================

"""
    TMDDPopulationResult

Result container for TMDD population simulation.

# Fields
- `individual_results`: Vector of TMDDSimResult for each subject
- `summary`: Summary statistics by observation (mean, SD, percentiles)
- `parameters`: Individual parameter values
- `metadata`: Simulation metadata
"""
struct TMDDPopulationResult
    individual_results::Vector{TMDDSimResult}
    summary::Dict{Symbol,Dict{String,Vector{Float64}}}
    parameters::Dict{Symbol,Vector{Float64}}
    metadata::Dict{String,Any}
end

# =============================================================================
# Population Simulation
# =============================================================================

"""
    simulate_tmdd_population(pop_spec, grid, solver) -> TMDDPopulationResult

Simulate TMDD model for a population of subjects with IIV.

# Arguments
- `pop_spec::TMDDPopulationSpec`: Population specification
- `grid::SimGrid`: Simulation time grid
- `solver::SolverSpec`: ODE solver configuration

# Returns
- `TMDDPopulationResult`: Individual and summary results

# Example
```julia
grid = SimGrid(0.0, 672.0, collect(0.0:1.0:672.0))
solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

result = simulate_tmdd_population(pop_spec, grid, solver)

# Access summary statistics
mean_conc = result.summary[:conc]["mean"]
p95_conc = result.summary[:conc]["p95"]
```
"""
function simulate_tmdd_population(
    pop_spec::TMDDPopulationSpec{K,P},
    grid::SimGrid,
    solver::SolverSpec
) where {K<:TMDDModelKind,P}

    base_spec = pop_spec.base_spec
    n_subjects = pop_spec.n_subjects

    # Generate individual parameters
    individual_params = _generate_individual_params(pop_spec)

    # Simulate each individual
    individual_results = Vector{TMDDSimResult}(undef, n_subjects)

    for i in 1:n_subjects
        # Create individual spec with perturbed parameters
        ind_spec = TMDDSpec(
            base_spec.kind,
            "$(base_spec.name)_ID$(i)",
            individual_params[i],
            base_spec.doses,
            base_spec.target_units,
            base_spec.drug_units
        )

        # Simulate
        individual_results[i] = solve_tmdd(ind_spec, grid, solver)
    end

    # Compute summary statistics
    summary = _compute_population_summary(individual_results)

    # Collect parameter values
    param_values = _collect_parameters(individual_params, P)

    # Metadata
    metadata = Dict{String,Any}(
        "model_type" => string(K),
        "approximation" => string(_get_approximation(base_spec.kind)),
        "route" => string(_get_route(base_spec.kind)),
        "n_subjects" => n_subjects,
        "has_iiv" => pop_spec.iiv !== nothing,
        "solver" => string(solver.alg)
    )

    return TMDDPopulationResult(individual_results, summary, param_values, metadata)
end

# =============================================================================
# Individual Parameter Generation
# =============================================================================

"""
Generate individual parameters with IIV for OneCptTMDD.
"""
function _generate_individual_params(
    pop_spec::TMDDPopulationSpec{OneCptTMDD,OneCptTMDDParams}
)
    base_params = pop_spec.base_spec.params
    n = pop_spec.n_subjects
    iiv = pop_spec.iiv

    if iiv === nothing
        return [base_params for _ in 1:n]
    end

    rng = StableRNG(iiv.seed)
    omegas = iiv.omegas

    params = Vector{OneCptTMDDParams}(undef, n)

    for i in 1:n
        CL = base_params.CL * exp(get(omegas, :CL, 0.0) * randn(rng))
        V = base_params.V * exp(get(omegas, :V, 0.0) * randn(rng))
        KSS = base_params.KSS * exp(get(omegas, :KSS, 0.0) * randn(rng))
        kint = base_params.kint * exp(get(omegas, :kint, 0.0) * randn(rng))
        ksyn = base_params.ksyn * exp(get(omegas, :ksyn, 0.0) * randn(rng))
        kdeg = base_params.kdeg * exp(get(omegas, :kdeg, 0.0) * randn(rng))
        R0 = base_params.R0 * exp(get(omegas, :R0, 0.0) * randn(rng))
        ka = base_params.ka * exp(get(omegas, :ka, 0.0) * randn(rng))
        F = base_params.F  # Usually no IIV on F directly

        params[i] = OneCptTMDDParams(CL, V, KSS, kint, ksyn, kdeg, R0, ka, F)
    end

    return params
end

"""
Generate individual parameters with IIV for TwoCptTMDD.
"""
function _generate_individual_params(
    pop_spec::TMDDPopulationSpec{TwoCptTMDD,TwoCptTMDDParams}
)
    base_params = pop_spec.base_spec.params
    n = pop_spec.n_subjects
    iiv = pop_spec.iiv

    if iiv === nothing
        return [base_params for _ in 1:n]
    end

    rng = StableRNG(iiv.seed)
    omegas = iiv.omegas

    params = Vector{TwoCptTMDDParams}(undef, n)

    for i in 1:n
        CL = base_params.CL * exp(get(omegas, :CL, 0.0) * randn(rng))
        V1 = base_params.V1 * exp(get(omegas, :V1, 0.0) * randn(rng))
        V2 = base_params.V2 * exp(get(omegas, :V2, 0.0) * randn(rng))
        Q = base_params.Q * exp(get(omegas, :Q, 0.0) * randn(rng))
        KSS = base_params.KSS * exp(get(omegas, :KSS, 0.0) * randn(rng))
        kint = base_params.kint * exp(get(omegas, :kint, 0.0) * randn(rng))
        ksyn = base_params.ksyn * exp(get(omegas, :ksyn, 0.0) * randn(rng))
        kdeg = base_params.kdeg * exp(get(omegas, :kdeg, 0.0) * randn(rng))
        R0 = base_params.R0 * exp(get(omegas, :R0, 0.0) * randn(rng))
        ka = base_params.ka * exp(get(omegas, :ka, 0.0) * randn(rng))
        F = base_params.F  # Usually no IIV on F
        Tlag = base_params.Tlag

        params[i] = TwoCptTMDDParams(CL, V1, V2, Q, KSS, kint, ksyn, kdeg, R0, ka, F, Tlag)
    end

    return params
end

"""
Generate individual parameters with IIV for TwoCptTMDDFcRn.
"""
function _generate_individual_params(
    pop_spec::TMDDPopulationSpec{TwoCptTMDDFcRn,TwoCptTMDDFcRnParams}
)
    base_params = pop_spec.base_spec.params
    n = pop_spec.n_subjects
    iiv = pop_spec.iiv

    if iiv === nothing
        return [base_params for _ in 1:n]
    end

    rng = StableRNG(iiv.seed)
    omegas = iiv.omegas

    params = Vector{TwoCptTMDDFcRnParams}(undef, n)

    for i in 1:n
        V1 = base_params.V1 * exp(get(omegas, :V1, 0.0) * randn(rng))
        V2 = base_params.V2 * exp(get(omegas, :V2, 0.0) * randn(rng))
        Q = base_params.Q * exp(get(omegas, :Q, 0.0) * randn(rng))
        CLup = base_params.CLup * exp(get(omegas, :CLup, 0.0) * randn(rng))
        FR = base_params.FR  # Usually fixed
        KSS = base_params.KSS * exp(get(omegas, :KSS, 0.0) * randn(rng))
        kint = base_params.kint * exp(get(omegas, :kint, 0.0) * randn(rng))
        ksyn = base_params.ksyn * exp(get(omegas, :ksyn, 0.0) * randn(rng))
        kdeg = base_params.kdeg * exp(get(omegas, :kdeg, 0.0) * randn(rng))
        R0 = base_params.R0 * exp(get(omegas, :R0, 0.0) * randn(rng))
        ka = base_params.ka * exp(get(omegas, :ka, 0.0) * randn(rng))
        F = base_params.F
        Tlag = base_params.Tlag

        params[i] = TwoCptTMDDFcRnParams(V1, V2, Q, CLup, FR, KSS, kint, ksyn, kdeg, R0, ka, F, Tlag)
    end

    return params
end

"""
Generate individual parameters with IIV for TwoCptTMDDADA.
"""
function _generate_individual_params(
    pop_spec::TMDDPopulationSpec{TwoCptTMDDADA,TwoCptTMDDADAParams}
)
    base_params = pop_spec.base_spec.params
    n = pop_spec.n_subjects
    iiv = pop_spec.iiv

    if iiv === nothing
        return [base_params for _ in 1:n]
    end

    rng = StableRNG(iiv.seed)
    omegas = iiv.omegas

    params = Vector{TwoCptTMDDADAParams}(undef, n)

    for i in 1:n
        CL = base_params.CL * exp(get(omegas, :CL, 0.0) * randn(rng))
        V1 = base_params.V1 * exp(get(omegas, :V1, 0.0) * randn(rng))
        V2 = base_params.V2 * exp(get(omegas, :V2, 0.0) * randn(rng))
        Q = base_params.Q * exp(get(omegas, :Q, 0.0) * randn(rng))
        KSS = base_params.KSS * exp(get(omegas, :KSS, 0.0) * randn(rng))
        kint = base_params.kint * exp(get(omegas, :kint, 0.0) * randn(rng))
        ksyn = base_params.ksyn * exp(get(omegas, :ksyn, 0.0) * randn(rng))
        kdeg = base_params.kdeg * exp(get(omegas, :kdeg, 0.0) * randn(rng))
        R0 = base_params.R0 * exp(get(omegas, :R0, 0.0) * randn(rng))
        kADA_prod = base_params.kADA_prod * exp(get(omegas, :kADA_prod, 0.0) * randn(rng))
        kADA_deg = base_params.kADA_deg
        kon_ADA = base_params.kon_ADA
        koff_ADA = base_params.koff_ADA
        CL_complex = base_params.CL_complex
        T_onset = base_params.T_onset * exp(get(omegas, :T_onset, 0.0) * randn(rng))
        ka = base_params.ka * exp(get(omegas, :ka, 0.0) * randn(rng))
        F = base_params.F

        params[i] = TwoCptTMDDADAParams(
            CL, V1, V2, Q, KSS, kint, ksyn, kdeg, R0,
            kADA_prod, kADA_deg, kon_ADA, koff_ADA, CL_complex, T_onset,
            ka, F
        )
    end

    return params
end

# Generic fallback for other model types
function _generate_individual_params(pop_spec::TMDDPopulationSpec)
    base_params = pop_spec.base_spec.params
    n = pop_spec.n_subjects
    # No IIV - return copies of base params
    return [base_params for _ in 1:n]
end

# =============================================================================
# Summary Statistics
# =============================================================================

"""
Compute summary statistics across population.
"""
function _compute_population_summary(results::Vector{TMDDSimResult})
    if isempty(results)
        return Dict{Symbol,Dict{String,Vector{Float64}}}()
    end

    t = results[1].t
    n_times = length(t)
    n_subjects = length(results)

    summary = Dict{Symbol,Dict{String,Vector{Float64}}}()

    # Get all observation keys from first result
    obs_keys = keys(results[1].observations)

    for key in obs_keys
        # Collect values across subjects
        values_matrix = zeros(n_subjects, n_times)
        for (i, res) in enumerate(results)
            if haskey(res.observations, key)
                obs_vals = res.observations[key]
                # Handle potential length mismatch
                len = min(length(obs_vals), n_times)
                values_matrix[i, 1:len] = obs_vals[1:len]
            end
        end

        # Compute statistics
        mean_values = vec(mean(values_matrix, dims=1))
        std_values = vec(std(values_matrix, dims=1))
        median_values = vec(median(values_matrix, dims=1))

        # Percentiles
        p5 = [quantile(values_matrix[:, j], 0.05) for j in 1:n_times]
        p10 = [quantile(values_matrix[:, j], 0.10) for j in 1:n_times]
        p25 = [quantile(values_matrix[:, j], 0.25) for j in 1:n_times]
        p75 = [quantile(values_matrix[:, j], 0.75) for j in 1:n_times]
        p90 = [quantile(values_matrix[:, j], 0.90) for j in 1:n_times]
        p95 = [quantile(values_matrix[:, j], 0.95) for j in 1:n_times]

        # Geometric mean and CV (for log-normally distributed data like conc)
        log_values = log.(max.(values_matrix, 1e-20))
        geom_mean = exp.(vec(mean(log_values, dims=1)))
        geom_cv = sqrt.(exp.(vec(var(log_values, dims=1))) .- 1.0) .* 100  # CV%

        summary[key] = Dict{String,Vector{Float64}}(
            "mean" => mean_values,
            "std" => std_values,
            "median" => median_values,
            "geom_mean" => geom_mean,
            "geom_cv" => geom_cv,
            "p5" => p5,
            "p10" => p10,
            "p25" => p25,
            "p75" => p75,
            "p90" => p90,
            "p95" => p95
        )
    end

    return summary
end

# =============================================================================
# Parameter Collection
# =============================================================================

"""
Collect parameter values from individual params for analysis.
"""
function _collect_parameters(params::Vector{OneCptTMDDParams}, ::Type{OneCptTMDDParams})
    return Dict{Symbol,Vector{Float64}}(
        :CL => [p.CL for p in params],
        :V => [p.V for p in params],
        :KSS => [p.KSS for p in params],
        :kint => [p.kint for p in params],
        :ksyn => [p.ksyn for p in params],
        :kdeg => [p.kdeg for p in params],
        :R0 => [p.R0 for p in params],
        :ka => [p.ka for p in params],
        :F => [p.F for p in params]
    )
end

function _collect_parameters(params::Vector{TwoCptTMDDParams}, ::Type{TwoCptTMDDParams})
    return Dict{Symbol,Vector{Float64}}(
        :CL => [p.CL for p in params],
        :V1 => [p.V1 for p in params],
        :V2 => [p.V2 for p in params],
        :Q => [p.Q for p in params],
        :KSS => [p.KSS for p in params],
        :kint => [p.kint for p in params],
        :ksyn => [p.ksyn for p in params],
        :kdeg => [p.kdeg for p in params],
        :R0 => [p.R0 for p in params],
        :ka => [p.ka for p in params],
        :F => [p.F for p in params],
        :Tlag => [p.Tlag for p in params]
    )
end

function _collect_parameters(params::Vector{TwoCptTMDDFcRnParams}, ::Type{TwoCptTMDDFcRnParams})
    return Dict{Symbol,Vector{Float64}}(
        :V1 => [p.V1 for p in params],
        :V2 => [p.V2 for p in params],
        :Q => [p.Q for p in params],
        :CLup => [p.CLup for p in params],
        :FR => [p.FR for p in params],
        :KSS => [p.KSS for p in params],
        :kint => [p.kint for p in params],
        :ksyn => [p.ksyn for p in params],
        :kdeg => [p.kdeg for p in params],
        :R0 => [p.R0 for p in params]
    )
end

function _collect_parameters(params::Vector{TwoCptTMDDADAParams}, ::Type{TwoCptTMDDADAParams})
    return Dict{Symbol,Vector{Float64}}(
        :CL => [p.CL for p in params],
        :V1 => [p.V1 for p in params],
        :V2 => [p.V2 for p in params],
        :Q => [p.Q for p in params],
        :KSS => [p.KSS for p in params],
        :kint => [p.kint for p in params],
        :ksyn => [p.ksyn for p in params],
        :kdeg => [p.kdeg for p in params],
        :R0 => [p.R0 for p in params],
        :kADA_prod => [p.kADA_prod for p in params],
        :T_onset => [p.T_onset for p in params]
    )
end

# Generic fallback
function _collect_parameters(params::Vector, ::Type)
    return Dict{Symbol,Vector{Float64}}()
end
