# NCA Module
# Non-Compartmental Analysis following FDA/EMA guidance

export run_nca, run_population_nca, summarize_population_nca
export nca_from_simresult, round_nca_result, round_nca_to_regulatory
export run_sparse_nca, SparseNCAResult, SparseNCAConfig
export validate_c0_extrapolation, C0ValidationResult
export RoundingRule, SIGNIFICANT_FIGURES, DECIMAL_PLACES, PHARMACOPEIAL

# Include NCA components
include("specs.jl")
include("lambda_z.jl")
include("auc.jl")
include("exposure_metrics.jl")
include("pk_parameters.jl")
include("multiple_dose.jl")
include("bioequivalence.jl")
include("reference_scaled_be.jl")

# =============================================================================
# Main NCA Entry Point
# =============================================================================

"""
    run_nca(t, c, dose; config=NCAConfig(), dosing_type=:single, tau=nothing, route=:extravascular)

Perform complete Non-Compartmental Analysis on concentration-time data.

# Arguments
- `t::Vector{Float64}`: Time points (sorted, ascending)
- `c::Vector{Float64}`: Concentration values
- `dose::Float64`: Administered dose

# Keyword Arguments
- `config::NCAConfig`: NCA configuration (default: standard FDA/EMA settings)
- `dosing_type::Symbol`: `:single`, `:multiple`, or `:steady_state` (default: :single)
- `tau::Union{Float64,Nothing}`: Dosing interval for multiple dose (required if dosing_type != :single)
- `route::Symbol`: Administration route - `:extravascular`, `:iv_bolus`, or `:iv_infusion`
- `t_inf::Float64`: Infusion duration (required for :iv_infusion)

# Returns
- `NCAResult`: Complete NCA results

# Example
```julia
t = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
c = [0.0, 1.2, 2.0, 1.8, 1.2, 0.6, 0.3, 0.075]
dose = 100.0

result = run_nca(t, c, dose)
println("Cmax: \$(result.cmax)")
println("AUC0-inf: \$(result.auc_0_inf)")
println("t1/2: \$(result.t_half)")
```
"""
function run_nca(
    t::Vector{Float64},
    c::Vector{Float64},
    dose::Float64;
    config::NCAConfig = NCAConfig(),
    dosing_type::Symbol = :single,
    tau::Union{Float64,Nothing} = nothing,
    route::Symbol = :extravascular,
    t_inf::Float64 = 0.0
)
    # Validate inputs
    n = length(t)
    @assert n == length(c) "Time and concentration vectors must have same length"
    @assert n >= 3 "Need at least 3 data points for NCA"
    @assert dose > 0.0 "Dose must be positive"
    @assert issorted(t) "Time points must be sorted ascending"

    if dosing_type in [:multiple, :steady_state]
        @assert tau !== nothing "Dosing interval (tau) required for multiple dose analysis"
        @assert tau > 0.0 "Dosing interval must be positive"
    end

    warnings = String[]
    quality_flags = Symbol[]

    # Apply BLQ handling
    c_clean = _apply_blq_handling(c, config)

    # ==========================================================================
    # Primary Exposure Metrics
    # ==========================================================================

    cmax, tmax_idx = find_cmax(t, c_clean)
    tmax = t[tmax_idx]

    clast, tlast, tlast_idx = find_clast(t, c_clean; lloq=config.lloq !== nothing ? config.lloq : 0.0)

    cmin = nothing
    cavg = nothing

    if dosing_type in [:multiple, :steady_state] && tau !== nothing
        cmin = nca_cmin(c_clean[findall(ti -> 0.0 <= ti <= tau, t)])
        cavg = nca_cavg(t, c_clean, tau, config)
    end

    # ==========================================================================
    # AUC Calculations
    # ==========================================================================

    auc_0t = auc_0_t(t, c_clean, config)
    aumc_0t = aumc_0_t(t, c_clean, config)

    # ==========================================================================
    # Lambda-z Estimation
    # ==========================================================================

    lambda_z_result = estimate_lambda_z(t, c_clean, config; tmax_idx=tmax_idx)

    # Extract lambda_z values
    lambda_z = lambda_z_result.lambda_z
    t_half = lambda_z_result.t_half

    if lambda_z_result.quality_flag == :insufficient
        push!(quality_flags, :lambda_z_insufficient)
        append!(warnings, lambda_z_result.warnings)
    elseif lambda_z_result.quality_flag == :warning
        push!(quality_flags, :lambda_z_warning)
        append!(warnings, lambda_z_result.warnings)
    end

    # ==========================================================================
    # AUC Extrapolation (if lambda_z available)
    # ==========================================================================

    auc_inf = nothing
    auc_extra_pct = nothing
    aumc_inf = nothing

    if lambda_z !== nothing && clast > 0.0
        auc_inf, auc_extra_pct = auc_0_inf(t, c_clean, lambda_z, clast, config)
        aumc_inf = aumc_0_inf(t, c_clean, lambda_z, clast, tlast, config)

        if auc_extra_pct > config.extrapolation_max_pct
            push!(quality_flags, :high_extrapolation)
            push!(warnings, "AUC extrapolation $(round(auc_extra_pct, digits=1))% exceeds $(config.extrapolation_max_pct)% threshold")
        end
    end

    # ==========================================================================
    # Multiple Dose / Steady State AUC
    # ==========================================================================

    auc_tau = nothing
    if tau !== nothing
        auc_tau = auc_0_tau(t, c_clean, tau, config)
    end

    # ==========================================================================
    # PK Parameters
    # ==========================================================================

    mrt = nothing
    cl_f = nothing
    vz_f = nothing
    vss = nothing

    if auc_inf !== nothing && aumc_inf !== nothing
        mrt = nca_mrt(aumc_inf, auc_inf; route=route, t_inf=t_inf)

        cl_f = nca_cl_f(dose, auc_inf)

        if lambda_z !== nothing
            vz_f = nca_vz_f(dose, lambda_z, auc_inf)
        end

        if route == :iv_bolus && mrt !== nothing && cl_f !== nothing
            vss = nca_vss(cl_f, mrt)
        end
    end

    # For IV bolus, CL and Vz (not CL/F and Vz/F)
    # The naming convention changes but calculation is same
    # Users should interpret based on route

    # ==========================================================================
    # Multiple Dose Metrics
    # ==========================================================================

    accumulation_index = nothing
    ptf = nothing
    swing = nothing

    if dosing_type in [:multiple, :steady_state] && tau !== nothing && cmin !== nothing
        if auc_inf !== nothing && auc_tau !== nothing
            accumulation_index = nca_accumulation_index(auc_tau, auc_inf)
        end

        if cavg !== nothing && cmin > 0.0
            ptf = nca_ptf(cmax, cmin, cavg)
            swing = nca_swing(cmax, cmin)
        end
    end

    # ==========================================================================
    # Dose-Normalized Metrics
    # ==========================================================================

    cmax_dn = nca_dose_normalized_cmax(cmax, dose)
    auc_dn = auc_inf !== nothing ? nca_dose_normalized_auc(auc_inf, dose) : nothing

    # ==========================================================================
    # Build Result
    # ==========================================================================

    metadata = Dict{String,Any}(
        "dose" => dose,
        "route" => String(route),
        "dosing_type" => String(dosing_type),
        "n_points" => n,
        "config_method" => string(typeof(config.method))
    )

    if tau !== nothing
        metadata["tau"] = tau
    end

    return NCAResult(
        # Primary exposure
        cmax, tmax, cmin, clast, tlast, cavg,
        # AUC
        auc_0t, auc_inf, auc_extra_pct, auc_tau, aumc_0t, aumc_inf,
        # Terminal phase
        lambda_z_result,
        # PK parameters
        t_half, mrt, cl_f, vz_f, vss,
        # Multiple dose
        accumulation_index, ptf, swing,
        # Dose-normalized
        cmax_dn, auc_dn,
        # Quality
        quality_flags, warnings, metadata
    )
end

# =============================================================================
# Population NCA
# =============================================================================

"""
    run_population_nca(pop_result, dose; config=NCAConfig(), observation=:conc)

Perform NCA analysis on each individual in a population simulation result.

# Arguments
- `pop_result::PopulationResult`: Population simulation result
- `dose::Float64`: Administered dose

# Keyword Arguments
- `config::NCAConfig`: NCA configuration
- `observation::Symbol`: Observation to analyze (default: :conc)
- `dosing_type::Symbol`: Dosing type (default: :single)
- `tau::Union{Float64,Nothing}`: Dosing interval for multiple dose
- `route::Symbol`: Administration route (default: :extravascular)

# Returns
- `Vector{NCAResult}`: NCA results for each individual
"""
function run_population_nca(
    pop_result::PopulationResult,
    dose::Float64;
    config::NCAConfig = NCAConfig(),
    observation::Symbol = :conc,
    dosing_type::Symbol = :single,
    tau::Union{Float64,Nothing} = nothing,
    route::Symbol = :extravascular
)
    n_individuals = length(pop_result.individuals)
    results = Vector{NCAResult}(undef, n_individuals)

    for i in 1:n_individuals
        ind = pop_result.individuals[i]

        t = ind.t
        c = ind.observations[observation]

        results[i] = run_nca(
            t, c, dose;
            config=config,
            dosing_type=dosing_type,
            tau=tau,
            route=route
        )
    end

    return results
end

"""
    summarize_population_nca(nca_results; parameters=[:cmax, :auc_0_inf, :t_half])

Summarize NCA results across a population.

# Arguments
- `nca_results::Vector{NCAResult}`: NCA results from population
- `parameters::Vector{Symbol}`: Parameters to summarize

# Returns
- `Dict{Symbol, NamedTuple}`: Summary statistics for each parameter
"""
function summarize_population_nca(
    nca_results::Vector{NCAResult};
    parameters::Vector{Symbol} = [:cmax, :auc_0_inf, :t_half, :cl_f]
)
    n = length(nca_results)
    summaries = Dict{Symbol, NamedTuple}()

    for param in parameters
        values = Float64[]

        for res in nca_results
            val = getfield(res, param)
            if val !== nothing
                push!(values, val)
            end
        end

        if !isempty(values)
            n_valid = length(values)
            mean_val = sum(values) / n_valid
            sorted_vals = sort(values)
            median_val = n_valid % 2 == 0 ?
                (sorted_vals[n_valid÷2] + sorted_vals[n_valid÷2+1]) / 2 :
                sorted_vals[(n_valid+1)÷2]
            sd_val = sqrt(sum((values .- mean_val).^2) / (n_valid - 1))
            cv_val = sd_val / mean_val * 100.0
            gm_val = exp(sum(log.(values)) / n_valid)

            summaries[param] = (
                n = n_valid,
                mean = mean_val,
                median = median_val,
                sd = sd_val,
                cv_pct = cv_val,
                geometric_mean = gm_val,
                min = minimum(values),
                max = maximum(values)
            )
        end
    end

    return summaries
end

# =============================================================================
# Helper Functions
# =============================================================================

"""
    _apply_blq_handling(c, config)

Apply BLQ handling to concentration data.
"""
function _apply_blq_handling(c::Vector{Float64}, config::NCAConfig)
    if config.lloq === nothing
        return c
    end

    c_clean = copy(c)
    lloq = config.lloq

    for i in eachindex(c_clean)
        if c_clean[i] < lloq
            if config.blq_handling isa BLQZero
                c_clean[i] = 0.0
            elseif config.blq_handling isa BLQLLOQHalf
                c_clean[i] = lloq / 2.0
            elseif config.blq_handling isa BLQMissing
                c_clean[i] = NaN  # Will be handled during calculations
            end
        end
    end

    return c_clean
end

"""
    nca_from_simresult(result, dose; config=NCAConfig(), route=:extravascular)

Convenience function to run NCA directly on a SimResult.

# Arguments
- `result::SimResult`: Simulation result
- `dose::Float64`: Administered dose

# Returns
- `NCAResult`: NCA analysis results
"""
function nca_from_simresult(
    result::SimResult,
    dose::Float64;
    config::NCAConfig = NCAConfig(),
    route::Symbol = :extravascular,
    observation::Symbol = :conc
)
    t = result.t
    c = result.observations[observation]

    return run_nca(t, c, dose; config=config, route=route)
end

# =============================================================================
# Regulatory Rounding
# =============================================================================

"""
    round_nca_result(value, significant_digits)

Round NCA result to specified significant digits for regulatory reporting.
"""
function round_nca_result(value::Float64, significant_digits::Int)
    if value == 0.0
        return 0.0
    end

    d = ceil(Int, log10(abs(value)))
    factor = 10.0^(significant_digits - d)

    return round(value * factor) / factor
end

# =============================================================================
# Regulatory-Style Rounding for All Outputs
# =============================================================================

"""
    RoundingRule

Regulatory rounding rules for NCA outputs.

- `SIGNIFICANT_FIGURES`: Standard significant figures rounding
- `DECIMAL_PLACES`: Fixed decimal places
- `PHARMACOPEIAL`: USP/EP style (3 sig figs, no trailing zeros)
"""
@enum RoundingRule begin
    SIGNIFICANT_FIGURES
    DECIMAL_PLACES
    PHARMACOPEIAL
end

"""
    round_nca_to_regulatory(value, rule, precision; parameter=nothing)

Apply regulatory rounding rules to NCA values.

# Arguments
- `value::Float64`: Value to round
- `rule::RoundingRule`: Rounding rule to apply
- `precision::Int`: Number of significant figures or decimal places

# Keyword Arguments
- `parameter::Union{Symbol,Nothing}`: Parameter type for smart rounding

# Returns
- `Float64`: Rounded value
"""
function round_nca_to_regulatory(
    value::Float64,
    rule::RoundingRule = SIGNIFICANT_FIGURES,
    precision::Int = 3;
    parameter::Union{Symbol,Nothing} = nothing
)
    if isnan(value) || isinf(value)
        return value
    end

    if value == 0.0
        return 0.0
    end

    # Parameter-specific rounding (FDA/EMA guidance)
    if parameter !== nothing
        precision = _get_parameter_precision(parameter)
    end

    if rule == SIGNIFICANT_FIGURES
        return round_nca_result(value, precision)
    elseif rule == DECIMAL_PLACES
        return round(value, digits=precision)
    elseif rule == PHARMACOPEIAL
        # USP/EP: 3 significant figures, half-up rounding
        return _pharmacopeial_round(value, precision)
    end

    return value
end

"""
Get recommended precision for specific NCA parameters per FDA/EMA guidance.
"""
function _get_parameter_precision(param::Symbol)::Int
    # High precision parameters (4 sig figs)
    if param in [:lambda_z, :t_half, :r_squared]
        return 4
    # Standard precision parameters (3 sig figs)
    elseif param in [:cmax, :auc_0_t, :auc_0_inf, :cl_f, :vz_f, :vss, :mrt]
        return 3
    # Time parameters (2 decimal places typically, but use sig figs)
    elseif param in [:tmax, :tlast]
        return 3
    # Percentage parameters
    elseif param in [:auc_extra_pct, :cv_pct, :ptf, :swing]
        return 3
    else
        return 3  # Default
    end
end

"""
Pharmacopeial rounding (USP/EP style).
Uses round-half-up (banker's rounding alternative).
"""
function _pharmacopeial_round(value::Float64, sig_figs::Int)
    if value == 0.0
        return 0.0
    end

    d = ceil(Int, log10(abs(value)))
    factor = 10.0^(sig_figs - d)

    # Use half-up rounding (not Julia's default half-even)
    scaled = value * factor
    if scaled >= 0
        rounded = floor(scaled + 0.5)
    else
        rounded = ceil(scaled - 0.5)
    end

    return rounded / factor
end

"""
    round_nca_result_all(result; config=NCAConfig())

Apply regulatory rounding to all numeric fields in an NCAResult.

# Arguments
- `result::NCAResult`: NCA result to round

# Returns
- `NCAResult`: New result with rounded values
"""
function round_nca_result_all(result::NCAResult; config::NCAConfig = NCAConfig())
    sig_dig = config.significant_digits

    return NCAResult(
        # Primary exposure
        round_nca_result(result.cmax, sig_dig),
        round_nca_result(result.tmax, sig_dig),
        result.cmin === nothing ? nothing : round_nca_result(result.cmin, sig_dig),
        round_nca_result(result.clast, sig_dig),
        round_nca_result(result.tlast, sig_dig),
        result.cavg === nothing ? nothing : round_nca_result(result.cavg, sig_dig),
        # AUC
        round_nca_result(result.auc_0_t, sig_dig),
        result.auc_0_inf === nothing ? nothing : round_nca_result(result.auc_0_inf, sig_dig),
        result.auc_extra_pct === nothing ? nothing : round_nca_result(result.auc_extra_pct, sig_dig),
        result.auc_0_tau === nothing ? nothing : round_nca_result(result.auc_0_tau, sig_dig),
        round_nca_result(result.aumc_0_t, sig_dig),
        result.aumc_0_inf === nothing ? nothing : round_nca_result(result.aumc_0_inf, sig_dig),
        # Terminal phase (keep as-is, has its own rounding)
        result.lambda_z_result,
        # PK parameters
        result.t_half === nothing ? nothing : round_nca_result(result.t_half, sig_dig),
        result.mrt === nothing ? nothing : round_nca_result(result.mrt, sig_dig),
        result.cl_f === nothing ? nothing : round_nca_result(result.cl_f, sig_dig),
        result.vz_f === nothing ? nothing : round_nca_result(result.vz_f, sig_dig),
        result.vss === nothing ? nothing : round_nca_result(result.vss, sig_dig),
        # Multiple dose
        result.accumulation_index === nothing ? nothing : round_nca_result(result.accumulation_index, sig_dig),
        result.ptf === nothing ? nothing : round_nca_result(result.ptf, sig_dig),
        result.swing === nothing ? nothing : round_nca_result(result.swing, sig_dig),
        # Dose-normalized
        result.cmax_dn === nothing ? nothing : round_nca_result(result.cmax_dn, sig_dig),
        result.auc_dn === nothing ? nothing : round_nca_result(result.auc_dn, sig_dig),
        # Quality
        result.quality_flags,
        result.warnings,
        result.metadata
    )
end

export round_nca_result_all

# =============================================================================
# C0 Back-Extrapolation Validation
# =============================================================================

"""
    C0ValidationResult

Result of C0 back-extrapolation validation.

# Fields
- `c0::Float64`: Back-extrapolated C0 value
- `is_valid::Bool`: Whether C0 passes validation rules
- `validation_flags::Vector{Symbol}`: Flags for validation issues
- `warnings::Vector{String}`: Warning messages
- `dose_c0_ratio::Union{Float64,Nothing}`: Dose/C0 ratio (should equal Vc)
- `expected_vc_range::Tuple{Float64,Float64}`: Expected Vc range for drug
"""
struct C0ValidationResult
    c0::Float64
    is_valid::Bool
    validation_flags::Vector{Symbol}
    warnings::Vector{String}
    dose_c0_ratio::Union{Float64,Nothing}
    expected_vc_range::Tuple{Float64,Float64}
end

"""
    validate_c0_extrapolation(c0, dose; expected_vc_range=(0.05, 1.0), max_extrapolation_time=nothing, first_sample_time=nothing)

Validate C0 back-extrapolation against physiological constraints.

# Validation Rules
1. C0 must be positive
2. Dose/C0 ratio should be physiologically plausible (Vc in L/kg range)
3. Extrapolation time should not exceed one half-life
4. C0 should not exceed 2x the first measured concentration (for IV bolus)

# Arguments
- `c0::Float64`: Back-extrapolated C0
- `dose::Float64`: Administered dose

# Keyword Arguments
- `expected_vc_range::Tuple{Float64,Float64}`: Expected Vc range in L/kg (default: 0.05-1.0)
- `body_weight::Float64`: Body weight in kg (default: 70)
- `max_extrapolation_time::Union{Float64,Nothing}`: Max acceptable extrapolation time
- `first_sample_time::Union{Float64,Nothing}`: Time of first sample
- `first_sample_conc::Union{Float64,Nothing}`: Concentration at first sample
- `t_half::Union{Float64,Nothing}`: Terminal half-life for extrapolation check

# Returns
- `C0ValidationResult`: Validation results
"""
function validate_c0_extrapolation(
    c0::Float64,
    dose::Float64;
    expected_vc_range::Tuple{Float64,Float64} = (0.05, 1.0),
    body_weight::Float64 = 70.0,
    max_extrapolation_time::Union{Float64,Nothing} = nothing,
    first_sample_time::Union{Float64,Nothing} = nothing,
    first_sample_conc::Union{Float64,Nothing} = nothing,
    t_half::Union{Float64,Nothing} = nothing
)
    validation_flags = Symbol[]
    warnings = String[]
    is_valid = true

    # Rule 1: C0 must be positive
    if c0 <= 0.0
        push!(validation_flags, :non_positive_c0)
        push!(warnings, "C0 must be positive, got $(c0)")
        is_valid = false
        return C0ValidationResult(c0, false, validation_flags, warnings, nothing, expected_vc_range)
    end

    # Calculate Vc = Dose / C0
    vc = dose / c0
    vc_per_kg = vc / body_weight

    # Rule 2: Check Vc is physiologically plausible
    vc_lower, vc_upper = expected_vc_range
    if vc_per_kg < vc_lower
        push!(validation_flags, :vc_too_low)
        push!(warnings, "Calculated Vc ($(round(vc_per_kg, digits=3)) L/kg) below expected range $(vc_lower)-$(vc_upper) L/kg")
        is_valid = false
    elseif vc_per_kg > vc_upper
        push!(validation_flags, :vc_too_high)
        push!(warnings, "Calculated Vc ($(round(vc_per_kg, digits=3)) L/kg) above expected range $(vc_lower)-$(vc_upper) L/kg")
        is_valid = false
    end

    # Rule 3: Extrapolation time should not exceed one half-life
    if max_extrapolation_time !== nothing && t_half !== nothing
        if max_extrapolation_time > t_half
            push!(validation_flags, :long_extrapolation)
            push!(warnings, "Extrapolation time ($(round(max_extrapolation_time, digits=2))) exceeds half-life ($(round(t_half, digits=2)))")
            # This is a warning, not necessarily invalid
        end
    end

    # Rule 4: C0 should not exceed 2x first measured concentration (for rapid distribution)
    if first_sample_conc !== nothing && first_sample_time !== nothing && first_sample_time > 0.0
        if c0 > 2.0 * first_sample_conc
            push!(validation_flags, :c0_too_high)
            push!(warnings, "C0 ($(round(c0, digits=3))) exceeds 2x first sample concentration ($(round(2.0 * first_sample_conc, digits=3)))")
            # This is a warning for possible multicompartment behavior
        end
    end

    return C0ValidationResult(
        c0,
        is_valid,
        validation_flags,
        warnings,
        vc,
        expected_vc_range
    )
end

# =============================================================================
# Sparse NCA for Population PK
# =============================================================================

"""
    SparseNCAConfig

Configuration for sparse NCA analysis (pooled data approach).

# Fields
- `pooling_method::Symbol`: `:naive_pooled`, `:destructive_sampling`, or `:linear_mixed`
- `time_bins::Union{Vector{Float64},Nothing}`: Time bins for grouping samples
- `time_tolerance::Float64`: Tolerance for assigning samples to nominal times
- `min_samples_per_bin::Int`: Minimum samples required per time bin
- `base_config::NCAConfig`: Underlying NCA configuration
"""
struct SparseNCAConfig
    pooling_method::Symbol
    time_bins::Union{Vector{Float64},Nothing}
    time_tolerance::Float64
    min_samples_per_bin::Int
    base_config::NCAConfig

    function SparseNCAConfig(;
        pooling_method::Symbol = :naive_pooled,
        time_bins::Union{Vector{Float64},Nothing} = nothing,
        time_tolerance::Float64 = 0.1,
        min_samples_per_bin::Int = 3,
        base_config::NCAConfig = NCAConfig()
    )
        @assert pooling_method in [:naive_pooled, :destructive_sampling, :linear_mixed] "Invalid pooling method"
        @assert time_tolerance >= 0.0 "Time tolerance must be non-negative"
        @assert min_samples_per_bin >= 1 "Minimum samples per bin must be at least 1"

        new(pooling_method, time_bins, time_tolerance, min_samples_per_bin, base_config)
    end
end

"""
    SparseNCAResult

Result of sparse NCA analysis.

# Fields
- `nca_result::NCAResult`: Pooled NCA results
- `n_subjects::Int`: Number of subjects contributing data
- `n_samples_per_time::Dict{Float64,Int}`: Number of samples at each nominal time
- `mean_profile::NamedTuple`: Mean concentration profile (time, mean_conc, sd, n)
- `individual_auc_contributions::Vector{Float64}`: Per-subject AUC contributions (if calculable)
- `pooling_method::Symbol`: Method used for pooling
- `quality_flags::Vector{Symbol}`: Quality warnings
"""
struct SparseNCAResult
    nca_result::NCAResult
    n_subjects::Int
    n_samples_per_time::Dict{Float64,Int}
    mean_profile::NamedTuple
    individual_auc_contributions::Vector{Float64}
    pooling_method::Symbol
    quality_flags::Vector{Symbol}
end

"""
    run_sparse_nca(subject_data, dose; config=SparseNCAConfig())

Perform sparse/pooled NCA for population PK analysis.

This method pools sparse concentration data from multiple subjects to
estimate population-level PK parameters, as is common in pediatric
studies or destructive sampling (e.g., mouse PK).

# Arguments
- `subject_data::Vector{NamedTuple}`: Vector of (subject_id, times, concentrations) tuples
- `dose::Float64`: Administered dose

# Keyword Arguments
- `config::SparseNCAConfig`: Sparse NCA configuration

# Returns
- `SparseNCAResult`: Pooled NCA results

# Example
```julia
# Destructive sampling from multiple animals
data = [
    (subject_id="Mouse1", times=[0.5], concs=[10.2]),
    (subject_id="Mouse2", times=[0.5], concs=[9.8]),
    (subject_id="Mouse3", times=[1.0], concs=[8.1]),
    (subject_id="Mouse4", times=[1.0], concs=[7.5]),
    (subject_id="Mouse5", times=[2.0], concs=[4.2]),
    (subject_id="Mouse6", times=[2.0], concs=[4.0]),
]

result = run_sparse_nca(data, 100.0)
```
"""
function run_sparse_nca(
    subject_data::Vector{<:NamedTuple},
    dose::Float64;
    config::SparseNCAConfig = SparseNCAConfig(),
    route::Symbol = :extravascular
)
    @assert !isempty(subject_data) "Subject data cannot be empty"
    @assert dose > 0.0 "Dose must be positive"

    n_subjects = length(subject_data)
    quality_flags = Symbol[]

    # Pool all data
    all_times = Float64[]
    all_concs = Float64[]
    subject_ids = String[]

    for subj in subject_data
        for (t, c) in zip(subj.times, subj.concs)
            push!(all_times, t)
            push!(all_concs, c)
            push!(subject_ids, string(subj.subject_id))
        end
    end

    # Bin data to nominal times if specified
    if config.time_bins !== nothing
        t_binned, c_binned, n_per_bin = _bin_sparse_data(
            all_times, all_concs, config.time_bins, config.time_tolerance
        )
    else
        # Auto-detect nominal times
        t_binned, c_binned, n_per_bin = _auto_bin_sparse_data(
            all_times, all_concs, config.time_tolerance, config.min_samples_per_bin
        )
    end

    # Check minimum samples per time point
    for (t, n) in n_per_bin
        if n < config.min_samples_per_bin
            push!(quality_flags, :insufficient_samples)
        end
    end

    # Calculate mean concentration profile
    mean_profile = _calculate_mean_profile(t_binned, c_binned, n_per_bin)

    # Run NCA on mean profile - sort by time
    t_mean = sort(collect(keys(mean_profile)))
    c_mean = [mean_profile[t].mean for t in t_mean]

    # Ensure time 0 is included if not present
    if t_mean[1] > 0.0
        pushfirst!(t_mean, 0.0)
        pushfirst!(c_mean, 0.0)  # Assume C=0 at t=0 for extravascular
    end

    nca_result = run_nca(
        t_mean, c_mean, dose;
        config=config.base_config,
        route=route
    )

    # Calculate individual AUC contributions (naive approach)
    individual_aucs = Float64[]
    for subj in subject_data
        if length(subj.times) >= 2
            # Simple linear AUC for available points
            t_subj = subj.times
            c_subj = subj.concs

            # Sort by time
            perm = sortperm(t_subj)
            t_sorted = t_subj[perm]
            c_sorted = c_subj[perm]

            auc_subj = 0.0
            for i in 2:length(t_sorted)
                auc_subj += 0.5 * (t_sorted[i] - t_sorted[i-1]) * (c_sorted[i-1] + c_sorted[i])
            end
            push!(individual_aucs, auc_subj)
        end
    end

    return SparseNCAResult(
        nca_result,
        n_subjects,
        n_per_bin,
        (times=t_mean, mean_conc=c_mean, profile=mean_profile),
        individual_aucs,
        config.pooling_method,
        quality_flags
    )
end

"""
Bin sparse data to nominal time points.
"""
function _bin_sparse_data(
    times::Vector{Float64},
    concs::Vector{Float64},
    nominal_times::Vector{Float64},
    tolerance::Float64
)
    binned_times = Dict{Float64,Vector{Float64}}()
    binned_concs = Dict{Float64,Vector{Float64}}()

    for nom_t in nominal_times
        binned_times[nom_t] = Float64[]
        binned_concs[nom_t] = Float64[]
    end

    for (t, c) in zip(times, concs)
        # Find closest nominal time within tolerance
        best_nom = nothing
        best_dist = Inf

        for nom_t in nominal_times
            dist = abs(t - nom_t)
            if dist <= tolerance * max(nom_t, 1.0) && dist < best_dist
                best_nom = nom_t
                best_dist = dist
            end
        end

        if best_nom !== nothing
            push!(binned_times[best_nom], t)
            push!(binned_concs[best_nom], c)
        end
    end

    n_per_bin = Dict{Float64,Int}()
    for nom_t in nominal_times
        n_per_bin[nom_t] = length(binned_concs[nom_t])
    end

    return binned_times, binned_concs, n_per_bin
end

"""
Auto-detect nominal times from sparse data.
"""
function _auto_bin_sparse_data(
    times::Vector{Float64},
    concs::Vector{Float64},
    tolerance::Float64,
    min_samples::Int
)
    # Sort unique times
    unique_times = sort(unique(times))

    # Group nearby times
    nominal_times = Float64[]
    current_group = [unique_times[1]]

    for t in unique_times[2:end]
        if t - current_group[end] <= tolerance * max(current_group[end], 1.0)
            push!(current_group, t)
        else
            # Finish current group
            push!(nominal_times, sum(current_group) / length(current_group))
            current_group = [t]
        end
    end
    push!(nominal_times, sum(current_group) / length(current_group))

    return _bin_sparse_data(times, concs, nominal_times, tolerance)
end

"""
Calculate mean profile from binned data.
"""
function _calculate_mean_profile(
    binned_times::Dict{Float64,Vector{Float64}},
    binned_concs::Dict{Float64,Vector{Float64}},
    n_per_bin::Dict{Float64,Int}
)
    profile = Dict{Float64,NamedTuple}()

    for (nom_t, concs) in binned_concs
        if !isempty(concs)
            n = length(concs)
            mean_c = sum(concs) / n
            sd_c = n > 1 ? sqrt(sum((concs .- mean_c).^2) / (n - 1)) : 0.0
            se_c = sd_c / sqrt(n)

            profile[nom_t] = (
                mean = mean_c,
                sd = sd_c,
                se = se_c,
                n = n,
                min = minimum(concs),
                max = maximum(concs)
            )
        end
    end

    return profile
end
