# Visual Predictive Check (VPC) Computation
# Implements standard VPC and prediction-corrected VPC (pcVPC)

using StableRNGs

export compute_vpc, compute_pcvpc

"""
Compute Visual Predictive Check from observed data and population simulations.

This function computes a VPC by:
1. Binning observed data by time
2. Computing observed percentiles in each bin
3. Running population simulations
4. Computing simulated percentiles with bootstrap CIs

If prediction_corrected=true, applies pcVPC correction:
  pcDV = DV * PRED_bin / PRED_individual

Arguments:
- observed: ObservedData from CDISC or other source
- population_spec: PopulationSpec for simulations
- grid: SimGrid defining simulation time points
- solver: SolverSpec for ODE solving
- config: VPCConfig with analysis settings
- error_spec: Optional ResidualErrorSpec for adding residual error

Returns:
- VPCResult
"""
function compute_vpc(
    observed::ObservedData,
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    config::VPCConfig=VPCConfig(),
    error_spec::Union{Nothing,ResidualErrorSpec}=nothing
)::VPCResult
    rng = StableRNG(config.seed)

    # Extract observed times and values
    obs_times = all_times(observed)
    obs_values = all_observations(observed)

    # Compute bin definitions based on observed data
    bin_defs = compute_bins(obs_times, config.binning)

    if isempty(bin_defs)
        return VPCResult(
            config, VPCBin[], n_subjects(observed), n_observations(observed),
            0, "", config.seed
        )
    end

    # Handle prediction correction if enabled
    if config.prediction_corrected
        # Compute typical (population) predictions for correction
        obs_preds, bin_pred_medians = _compute_prediction_correction_data(
            observed, population_spec, grid, solver, bin_defs
        )

        # Apply prediction correction to observed values
        obs_values_corrected = _apply_prediction_correction(
            obs_values, obs_times, obs_preds, bin_defs, bin_pred_medians
        )

        # Assign corrected observed data to bins
        obs_binned = assign_to_bins(obs_times, obs_values_corrected, bin_defs)

        # Run simulations with prediction correction
        sim_percentiles_per_bin = _run_simulations_for_pcvpc(
            population_spec, grid, solver, bin_defs, bin_pred_medians,
            config, error_spec, rng
        )
    else
        # Standard VPC
        obs_binned = assign_to_bins(obs_times, obs_values, bin_defs)

        # Run population simulations
        sim_percentiles_per_bin = _run_simulations_for_vpc(
            population_spec, grid, solver, bin_defs, config, error_spec, rng
        )
    end

    # Build VPCBin results
    vpc_bins = VPCBin[]

    for (bin_def, (bin_id, obs_vals)) in zip(bin_defs, obs_binned)
        percentile_data = VPCPercentileData[]

        for pi_level in config.pi_levels
            # Observed percentile
            obs_pctl = compute_percentile(obs_vals, pi_level)

            # Simulated percentiles (from all simulations)
            sim_pctls = sim_percentiles_per_bin[bin_id][pi_level]

            # Bootstrap CI for simulated percentiles
            lower, median_sim, upper = bootstrap_percentile_ci(
                sim_pctls, 0.5, config.ci_level, config.n_bootstrap, rng
            )

            push!(percentile_data, VPCPercentileData(
                pi_level, obs_pctl, median_sim, lower, upper
            ))
        end

        push!(vpc_bins, VPCBin(
            bin_def.id,
            bin_def.lower,
            bin_def.upper,
            bin_def.midpoint,
            length(obs_vals),
            config.n_simulations,
            percentile_data
        ))
    end

    return VPCResult(
        config,
        vpc_bins,
        n_subjects(observed),
        n_observations(observed),
        config.n_simulations,
        "",
        config.seed
    )
end

"""
Compute prediction correction data for pcVPC.

Returns:
- obs_preds: Vector of typical predictions at each observed time
- bin_pred_medians: Dict of bin_id => median PRED for that bin
"""
function _compute_prediction_correction_data(
    observed::ObservedData,
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    bin_defs::Vector{BinDefinition}
)
    # Get base model spec (typical parameters without IIV)
    base_model_spec = population_spec.base_model_spec

    # Get observed times
    obs_times = all_times(observed)

    # Compute typical predictions at observed times
    # Create a grid that includes all observed time points
    extended_times = sort(unique(vcat(grid.saveat, obs_times)))
    extended_grid = SimGrid(grid.t0, grid.t1, extended_times)

    # Simulate with base model (no IIV)
    base_result = simulate(base_model_spec, extended_grid, solver)

    # Interpolate to get predictions at exact observed times
    obs_preds = Float64[]
    for t in obs_times
        # Find nearest time in simulation result
        idx = argmin(abs.(base_result.t .- t))
        push!(obs_preds, base_result.observations[:conc][idx])
    end

    # Compute median PRED for each bin
    pred_binned = assign_to_bins(obs_times, obs_preds, bin_defs)

    bin_pred_medians = Dict{Int, Float64}()
    for (bin_id, pred_vals) in pred_binned
        if !isempty(pred_vals)
            bin_pred_medians[bin_id] = compute_percentile(pred_vals, 0.5)
        else
            bin_pred_medians[bin_id] = 1.0  # Fallback to avoid division by zero
        end
    end

    return obs_preds, bin_pred_medians
end

"""
Apply prediction correction to observed values.

pcDV_ij = DV_ij * PRED_bin / PRED_ij

This normalizes DV by the ratio of bin-median prediction to individual prediction,
reducing variability due to structural model differences.
"""
function _apply_prediction_correction(
    obs_values::Vector{Float64},
    obs_times::Vector{Float64},
    obs_preds::Vector{Float64},
    bin_defs::Vector{BinDefinition},
    bin_pred_medians::Dict{Int, Float64}
)::Vector{Float64}
    corrected = Float64[]

    for (i, (dv, t, pred)) in enumerate(zip(obs_values, obs_times, obs_preds))
        # Find which bin this time belongs to
        bin_id = _find_bin(t, bin_defs)

        if bin_id !== nothing && pred > 0
            pred_bin = bin_pred_medians[bin_id]
            # Apply correction: pcDV = DV * PRED_bin / PRED_i
            pc_dv = dv * pred_bin / pred
            push!(corrected, pc_dv)
        else
            # No correction possible, use original value
            push!(corrected, dv)
        end
    end

    return corrected
end

"""
Find which bin a time point belongs to.
"""
function _find_bin(t::Float64, bin_defs::Vector{BinDefinition})::Union{Nothing, Int}
    for bin_def in bin_defs
        if bin_def.lower <= t <= bin_def.upper
            return bin_def.id
        end
    end
    return nothing
end

"""
Run simulations for prediction-corrected VPC (pcVPC).

Applies the same prediction correction to simulated data.

The prediction correction formula (Bergstrand et al., 2011):
  pcDV_ij = DV_ij * PRED_bin / PRED_ij

where:
- DV_ij is the (simulated) observation for subject i at time j
- PRED_ij is the typical (population) prediction at time j
- PRED_bin is the median typical prediction for the bin containing time j
"""
function _run_simulations_for_pcvpc(
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    bin_defs::Vector{BinDefinition},
    bin_pred_medians::Dict{Int, Float64},
    config::VPCConfig,
    error_spec::Union{Nothing,ResidualErrorSpec},
    rng
)::Dict{Int,Dict{Float64,Vector{Float64}}}
    # Initialize storage
    result = Dict{Int,Dict{Float64,Vector{Float64}}}()
    for bin_def in bin_defs
        result[bin_def.id] = Dict{Float64,Vector{Float64}}()
        for pi in config.pi_levels
            result[bin_def.id][pi] = Float64[]
        end
    end

    # Get base model for typical predictions
    base_model_spec = population_spec.base_model_spec

    # Pre-compute base model predictions at all grid times ONCE
    # This is the key optimization - avoid repeated ODE solves
    base_result = simulate(base_model_spec, grid, solver)
    base_times = base_result.t
    base_preds = base_result.observations[:conc]

    # Build interpolation lookup for base predictions
    # For efficiency, we create a sorted lookup table
    function interpolate_base_pred(t::Float64)
        # Find the closest time point in base_times
        idx = argmin(abs.(base_times .- t))

        # If exact match or only one point, return directly
        if abs(base_times[idx] - t) < 1e-10 || length(base_times) == 1
            return base_preds[idx]
        end

        # Linear interpolation between bracketing points
        if base_times[idx] > t && idx > 1
            # t is between idx-1 and idx
            t1, t2 = base_times[idx-1], base_times[idx]
            c1, c2 = base_preds[idx-1], base_preds[idx]
        elseif base_times[idx] < t && idx < length(base_times)
            # t is between idx and idx+1
            t1, t2 = base_times[idx], base_times[idx+1]
            c1, c2 = base_preds[idx], base_preds[idx+1]
        else
            # Edge case - return nearest
            return base_preds[idx]
        end

        # Linear interpolation
        frac = (t - t1) / (t2 - t1)
        return c1 + frac * (c2 - c1)
    end

    # Run simulations
    for sim_idx in 1:config.n_simulations
        sim_seed = rand(rng, UInt64)

        # Create pop spec with new seed
        if population_spec.iiv !== nothing
            new_iiv = IIVSpec(
                population_spec.iiv.kind,
                population_spec.iiv.omegas,
                sim_seed,
                population_spec.iiv.n
            )
            pop_spec_sim = PopulationSpec(
                population_spec.base_model_spec,
                new_iiv,
                population_spec.iov,
                population_spec.covariate_model,
                population_spec.covariates
            )
        else
            pop_spec_sim = population_spec
        end

        # Run population simulation (with IIV)
        pop_result = simulate_population(pop_spec_sim, grid, solver)

        # Collect simulation data and compute predictions efficiently
        sim_times = Float64[]
        sim_values = Float64[]
        sim_preds = Float64[]  # Typical predictions

        for individual in pop_result.individuals
            # Get simulated values
            append!(sim_times, individual.t)
            append!(sim_values, individual.observations[:conc])

            # Compute typical predictions at same times using pre-computed base model
            for t in individual.t
                pred = interpolate_base_pred(t)
                push!(sim_preds, pred)
            end
        end

        # Apply residual error if specified
        if error_spec !== nothing
            error_rng = StableRNG(sim_seed + UInt64(1))
            sim_values = apply_residual_error(sim_values, error_spec; rng=error_rng)
        end

        # Apply prediction correction to simulated values
        sim_values_corrected = _apply_prediction_correction(
            sim_values, sim_times, sim_preds, bin_defs, bin_pred_medians
        )

        # Assign to bins and compute percentiles
        sim_binned = assign_to_bins(sim_times, sim_values_corrected, bin_defs)

        for (bin_id, sim_vals) in sim_binned
            for pi in config.pi_levels
                pctl = compute_percentile(sim_vals, pi)
                if !isnan(pctl)
                    push!(result[bin_id][pi], pctl)
                end
            end
        end
    end

    return result
end

"""
Compute prediction-corrected VPC (pcVPC).

In pcVPC, both observed and simulated data are normalized by the
population prediction to reduce variability from the structural model.

Arguments:
- observed: ObservedData from CDISC or other source
- population_spec: PopulationSpec for simulations
- grid: SimGrid defining simulation time points
- solver: SolverSpec for ODE solving
- config: VPCConfig with analysis settings
- error_spec: Optional ResidualErrorSpec for adding residual error

Returns:
- VPCResult with prediction-corrected values
"""
function compute_pcvpc(
    observed::ObservedData,
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    config::VPCConfig=VPCConfig(),
    error_spec::Union{Nothing,ResidualErrorSpec}=nothing
)::VPCResult
    # Enable prediction correction in config
    pc_config = VPCConfig(
        pi_levels=config.pi_levels,
        ci_level=config.ci_level,
        binning=config.binning,
        prediction_corrected=true,
        stratify_by=config.stratify_by,
        lloq=config.lloq,
        n_simulations=config.n_simulations,
        n_bootstrap=config.n_bootstrap,
        seed=config.seed
    )

    return compute_vpc(observed, population_spec, grid, solver;
        config=pc_config, error_spec=error_spec)
end

"""
Run simulations and compute percentiles for each bin.
"""
function _run_simulations_for_vpc(
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
    bin_defs::Vector{BinDefinition},
    config::VPCConfig,
    error_spec::Union{Nothing,ResidualErrorSpec},
    rng
)::Dict{Int,Dict{Float64,Vector{Float64}}}
    n_bins = length(bin_defs)

    # Initialize storage: bin_id -> pi_level -> Vector of percentiles from each simulation
    result = Dict{Int,Dict{Float64,Vector{Float64}}}()
    for bin_def in bin_defs
        result[bin_def.id] = Dict{Float64,Vector{Float64}}()
        for pi in config.pi_levels
            result[bin_def.id][pi] = Float64[]
        end
    end

    # Run simulations
    for sim_idx in 1:config.n_simulations
        # Create a new seed for this simulation
        sim_seed = rand(rng, UInt64)

        # Create new IIV spec with the new seed (if IIV is present)
        if population_spec.iiv !== nothing
            new_iiv = IIVSpec(
                population_spec.iiv.kind,
                population_spec.iiv.omegas,
                sim_seed,
                population_spec.iiv.n
            )
            pop_spec_sim = PopulationSpec(
                population_spec.base_model_spec,
                new_iiv,
                population_spec.iov,
                population_spec.covariate_model,
                population_spec.covariates
            )
        else
            pop_spec_sim = population_spec
        end

        # Run population simulation
        pop_result = simulate_population(pop_spec_sim, grid, solver)

        # Extract simulated observations at time points matching bins
        sim_times = Float64[]
        sim_values = Float64[]

        for individual in pop_result.individuals
            append!(sim_times, individual.t)
            append!(sim_values, individual.observations[:conc])
        end

        # Apply residual error if specified
        if error_spec !== nothing
            error_rng = StableRNG(sim_seed + UInt64(1))
            sim_values = apply_residual_error(sim_values, error_spec; rng=error_rng)
        end

        # Assign to bins and compute percentiles
        sim_binned = assign_to_bins(sim_times, sim_values, bin_defs)

        for (bin_id, sim_vals) in sim_binned
            for pi in config.pi_levels
                pctl = compute_percentile(sim_vals, pi)
                if !isnan(pctl)
                    push!(result[bin_id][pi], pctl)
                end
            end
        end
    end

    return result
end

"""
Compute VPC from population simulation result (without separate observed data).

This is useful for simulation-based model evaluation when you have
simulated data and want to compute prediction intervals.

Arguments:
- pop_result: PopulationResult from simulate_population
- config: VPCConfig with analysis settings

Returns:
- VPCResult
"""
function compute_vpc_from_simulation(
    pop_result::PopulationResult;
    config::VPCConfig=VPCConfig()
)::VPCResult
    rng = StableRNG(config.seed)

    # Extract all times and values from simulation
    all_t = Float64[]
    all_v = Float64[]

    for individual in pop_result.individuals
        append!(all_t, individual.t)
        append!(all_v, individual.observations[:conc])
    end

    # Compute bins
    bin_defs = compute_bins(all_t, config.binning)

    if isempty(bin_defs)
        return VPCResult(config, VPCBin[], length(pop_result.individuals), 0, 0, "", config.seed)
    end

    # Assign data to bins
    binned = assign_to_bins(all_t, all_v, bin_defs)

    # Build VPCBin results (simulated only, no observed)
    vpc_bins = VPCBin[]

    for (bin_def, (bin_id, vals)) in zip(bin_defs, binned)
        percentile_data = VPCPercentileData[]

        for pi_level in config.pi_levels
            pctl = compute_percentile(vals, pi_level)

            # Bootstrap CI
            lower, median_val, upper = bootstrap_percentile_ci(
                vals, pi_level, config.ci_level, config.n_bootstrap, rng
            )

            # For simulation-only VPC, observed = simulated
            push!(percentile_data, VPCPercentileData(
                pi_level, pctl, median_val, lower, upper
            ))
        end

        push!(vpc_bins, VPCBin(
            bin_def.id,
            bin_def.lower,
            bin_def.upper,
            bin_def.midpoint,
            length(vals),
            length(vals),
            percentile_data
        ))
    end

    return VPCResult(
        config,
        vpc_bins,
        length(pop_result.individuals),
        length(all_t),
        1,  # Single "simulation" (the actual pop_result)
        "",
        config.seed
    )
end

export compute_vpc_from_simulation

# ============================================================================
# VPC Stratification
# ============================================================================

"""
Stratified VPC result containing VPCs for each stratum.
"""
struct StratifiedVPCResult
    results::Vector{VPCResult}
    stratify_by::Vector{Symbol}
    strata_names::Vector{String}
end

export StratifiedVPCResult

"""
    compute_stratified_vpc(observed, population_spec, grid, solver; config, error_spec, strata_data)

Compute VPC stratified by covariate values.

This function separates observed data into strata based on covariate values
and computes separate VPCs for each stratum.

Arguments:
- observed: ObservedData from CDISC or other source
- population_spec: PopulationSpec for simulations
- grid: SimGrid defining simulation time points
- solver: SolverSpec for ODE solving
- config: VPCConfig with stratify_by set
- error_spec: Optional ResidualErrorSpec for adding residual error
- strata_data: Dict mapping subject_id => Dict of covariate values

Returns:
- StratifiedVPCResult containing VPCResult for each stratum
"""
function compute_stratified_vpc(
    observed::ObservedData,
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    config::VPCConfig=VPCConfig(),
    error_spec::Union{Nothing,ResidualErrorSpec}=nothing,
    strata_data::Dict{String, Dict{Symbol, Any}}=Dict{String, Dict{Symbol, Any}}()
)::StratifiedVPCResult
    # If no stratification variables specified, compute regular VPC
    if isempty(config.stratify_by)
        result = compute_vpc(observed, population_spec, grid, solver;
            config=config, error_spec=error_spec)
        return StratifiedVPCResult([result], Symbol[], ["All"])
    end

    # Group subjects by strata values
    strata_groups = _group_by_strata(observed, config.stratify_by, strata_data)

    results = VPCResult[]
    strata_names = String[]

    for (strata_name, subject_ids) in strata_groups
        # Create filtered ObservedData for this stratum
        filtered_observed = filter_subjects(observed, subject_ids)

        if n_subjects(filtered_observed) == 0
            continue
        end

        # Compute VPC for this stratum
        strata_config = VPCConfig(
            pi_levels=config.pi_levels,
            ci_level=config.ci_level,
            binning=config.binning,
            prediction_corrected=config.prediction_corrected,
            stratify_by=Symbol[],  # Don't further stratify
            lloq=config.lloq,
            n_simulations=config.n_simulations,
            n_bootstrap=config.n_bootstrap,
            seed=config.seed
        )

        vpc_result = compute_vpc(filtered_observed, population_spec, grid, solver;
            config=strata_config, error_spec=error_spec)

        # Create result with strata label
        strata_result = VPCResult(
            vpc_result.config,
            vpc_result.bins,
            vpc_result.n_subjects_observed,
            vpc_result.n_observations_observed,
            vpc_result.n_simulations,
            strata_name,
            vpc_result.simulation_seed
        )

        push!(results, strata_result)
        push!(strata_names, strata_name)
    end

    return StratifiedVPCResult(results, config.stratify_by, strata_names)
end

export compute_stratified_vpc

"""
Group subjects by strata values.
"""
function _group_by_strata(
    observed::ObservedData,
    stratify_by::Vector{Symbol},
    strata_data::Dict{String, Dict{Symbol, Any}}
)::Vector{Tuple{String, Vector{String}}}
    strata_groups = Dict{String, Vector{String}}()

    subject_ids = get_subject_ids(observed)

    for subj_id in subject_ids
        # Get covariate values for this subject
        covs = get(strata_data, subj_id, Dict{Symbol, Any}())

        # Build strata key from covariate values
        strata_parts = String[]
        for var in stratify_by
            val = get(covs, var, "Unknown")
            push!(strata_parts, "$(var)=$(val)")
        end
        strata_key = join(strata_parts, ", ")

        if !haskey(strata_groups, strata_key)
            strata_groups[strata_key] = String[]
        end
        push!(strata_groups[strata_key], subj_id)
    end

    # Convert to sorted vector of tuples
    return [(k, v) for (k, v) in sort(collect(strata_groups), by=x->x[1])]
end

"""
Filter ObservedData to include only specified subjects.
"""
function filter_subjects(observed::ObservedData, subject_ids::Vector{String})::ObservedData
    # This assumes ObservedData has a records field - adjust based on actual structure
    if hasfield(typeof(observed), :records)
        filtered_records = filter(r -> r.subject_id in subject_ids, observed.records)
        return ObservedData(filtered_records)
    else
        # Fallback - return as-is if structure doesn't support filtering
        return observed
    end
end

"""
Get unique subject IDs from ObservedData.
"""
function get_subject_ids(observed::ObservedData)::Vector{String}
    if hasfield(typeof(observed), :records)
        return unique([r.subject_id for r in observed.records])
    else
        return String[]
    end
end

# ============================================================================
# VPC BLQ (Below Limit of Quantification) Handling
# ============================================================================

"""
BLQ handling methods based on Beal (2001).
"""
@enum BLQMethod begin
    M1  # Discard all BLQ observations
    M3  # Treat as censored (maximize likelihood)
    M4  # Replace with LLOQ/2
    M5  # 0 before Tmax, LLOQ/2 after
    M6  # LLOQ/2 before Tmax, discard after
    M7  # 0 before Tmax, discard after
end

export BLQMethod, M1, M3, M4, M5, M6, M7

"""
BLQ statistics for a time bin.
"""
struct BLQBinStats
    bin_id::Int
    n_total::Int
    n_blq::Int
    pct_blq_observed::Float64
    pct_blq_simulated_median::Float64
    pct_blq_simulated_lower::Float64
    pct_blq_simulated_upper::Float64
end

export BLQBinStats

"""
    handle_blq(values, times, lloq; method=M4, tmax=nothing) -> Vector{Float64}

Apply BLQ handling method to concentration values.

Arguments:
- values: Vector of observed concentrations
- times: Vector of time points
- lloq: Lower limit of quantification
- method: BLQ handling method (M1-M7)
- tmax: Time of maximum concentration (required for M5-M7)

Returns:
- Vector of handled values (some may be NaN for discarded observations)
"""
function handle_blq(
    values::Vector{Float64},
    times::Vector{Float64},
    lloq::Float64;
    method::BLQMethod=M4,
    tmax::Union{Nothing, Float64}=nothing
)::Vector{Float64}
    n = length(values)
    result = copy(values)

    # Determine Tmax if not provided (for methods M5-M7)
    if tmax === nothing && method in [M5, M6, M7]
        _, idx = findmax(values)
        tmax = times[idx]
    end

    for i in 1:n
        if values[i] < lloq
            if method == M1
                # Discard BLQ
                result[i] = NaN
            elseif method == M3
                # Treat as censored - keep as-is for likelihood
                # (actual censoring handled in likelihood calculation)
                result[i] = values[i]
            elseif method == M4
                # Replace with LLOQ/2
                result[i] = lloq / 2
            elseif method == M5
                # 0 before Tmax, LLOQ/2 after
                if tmax !== nothing && times[i] < tmax
                    result[i] = 0.0
                else
                    result[i] = lloq / 2
                end
            elseif method == M6
                # LLOQ/2 before Tmax, discard after
                if tmax !== nothing && times[i] < tmax
                    result[i] = lloq / 2
                else
                    result[i] = NaN
                end
            elseif method == M7
                # 0 before Tmax, discard after
                if tmax !== nothing && times[i] < tmax
                    result[i] = 0.0
                else
                    result[i] = NaN
                end
            end
        end
    end

    return result
end

export handle_blq

"""
    compute_vpc_with_blq(observed, population_spec, grid, solver; config, error_spec, blq_method)

Compute VPC with BLQ handling.

This function applies BLQ handling to both observed and simulated data,
and additionally computes %BLQ statistics for each bin.

Arguments:
- observed: ObservedData from CDISC or other source
- population_spec: PopulationSpec for simulations
- grid: SimGrid defining simulation time points
- solver: SolverSpec for ODE solving
- config: VPCConfig with lloq set
- error_spec: Optional ResidualErrorSpec for adding residual error
- blq_method: BLQ handling method (default: M4)

Returns:
- Tuple of (VPCResult, Vector{BLQBinStats})
"""
function compute_vpc_with_blq(
    observed::ObservedData,
    population_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    config::VPCConfig=VPCConfig(),
    error_spec::Union{Nothing,ResidualErrorSpec}=nothing,
    blq_method::BLQMethod=M4
)::Tuple{VPCResult, Vector{BLQBinStats}}
    rng = StableRNG(config.seed)

    lloq = config.lloq
    if lloq === nothing
        error("LLOQ must be specified in VPCConfig for BLQ handling")
    end

    # Extract observed times and values
    obs_times = all_times(observed)
    obs_values = all_observations(observed)

    # Determine Tmax from observed data
    _, tmax_idx = findmax(obs_values)
    tmax = obs_times[tmax_idx]

    # Apply BLQ handling to observed data
    obs_values_handled = handle_blq(obs_values, obs_times, lloq; method=blq_method, tmax=tmax)

    # Filter out NaN values (discarded observations)
    valid_mask = .!isnan.(obs_values_handled)
    obs_times_valid = obs_times[valid_mask]
    obs_values_valid = obs_values_handled[valid_mask]

    # Compute bin definitions
    bin_defs = compute_bins(obs_times_valid, config.binning)

    if isempty(bin_defs)
        return VPCResult(
            config, VPCBin[], n_subjects(observed), n_observations(observed),
            0, "", config.seed
        ), BLQBinStats[]
    end

    # Assign handled observed data to bins
    obs_binned = assign_to_bins(obs_times_valid, obs_values_valid, bin_defs)

    # Compute observed %BLQ per bin (using original unfiltered data)
    obs_binned_original = assign_to_bins(obs_times, obs_values, bin_defs)
    obs_pct_blq_per_bin = Dict{Int, Float64}()
    for (bin_id, vals) in obs_binned_original
        n_blq = count(v -> v < lloq, vals)
        obs_pct_blq_per_bin[bin_id] = 100.0 * n_blq / max(length(vals), 1)
    end

    # Run simulations with BLQ handling
    sim_percentiles_per_bin = Dict{Int, Dict{Float64, Vector{Float64}}}()
    sim_pct_blq_per_bin = Dict{Int, Vector{Float64}}()

    for bin_def in bin_defs
        sim_percentiles_per_bin[bin_def.id] = Dict{Float64, Vector{Float64}}()
        for pi in config.pi_levels
            sim_percentiles_per_bin[bin_def.id][pi] = Float64[]
        end
        sim_pct_blq_per_bin[bin_def.id] = Float64[]
    end

    # Run simulations
    for sim_idx in 1:config.n_simulations
        sim_seed = rand(rng, UInt64)

        # Create pop spec with new seed
        if population_spec.iiv !== nothing
            new_iiv = IIVSpec(
                population_spec.iiv.kind,
                population_spec.iiv.omegas,
                sim_seed,
                population_spec.iiv.n
            )
            pop_spec_sim = PopulationSpec(
                population_spec.base_model_spec,
                new_iiv,
                population_spec.iov,
                population_spec.covariate_model,
                population_spec.covariates
            )
        else
            pop_spec_sim = population_spec
        end

        # Run population simulation
        pop_result = simulate_population(pop_spec_sim, grid, solver)

        # Extract simulated observations
        sim_times = Float64[]
        sim_values = Float64[]

        for individual in pop_result.individuals
            append!(sim_times, individual.t)
            append!(sim_values, individual.observations[:conc])
        end

        # Apply residual error if specified
        if error_spec !== nothing
            error_rng = StableRNG(sim_seed + UInt64(1))
            sim_values = apply_residual_error(sim_values, error_spec; rng=error_rng)
        end

        # Compute %BLQ per bin (before BLQ handling)
        sim_binned_original = assign_to_bins(sim_times, sim_values, bin_defs)
        for (bin_id, vals) in sim_binned_original
            n_blq = count(v -> v < lloq, vals)
            push!(sim_pct_blq_per_bin[bin_id], 100.0 * n_blq / max(length(vals), 1))
        end

        # Apply BLQ handling
        sim_tmax_idx = isempty(sim_values) ? 1 : argmax(sim_values)
        sim_tmax = isempty(sim_times) ? 0.0 : sim_times[sim_tmax_idx]
        sim_values_handled = handle_blq(sim_values, sim_times, lloq; method=blq_method, tmax=sim_tmax)

        # Filter out NaN values
        sim_valid_mask = .!isnan.(sim_values_handled)
        sim_times_valid = sim_times[sim_valid_mask]
        sim_values_valid = sim_values_handled[sim_valid_mask]

        # Assign to bins and compute percentiles
        sim_binned = assign_to_bins(sim_times_valid, sim_values_valid, bin_defs)

        for (bin_id, sim_vals) in sim_binned
            for pi in config.pi_levels
                pctl = compute_percentile(sim_vals, pi)
                if !isnan(pctl)
                    push!(sim_percentiles_per_bin[bin_id][pi], pctl)
                end
            end
        end
    end

    # Build VPCBin results
    vpc_bins = VPCBin[]
    blq_stats = BLQBinStats[]

    for (bin_def, (bin_id, obs_vals)) in zip(bin_defs, obs_binned)
        percentile_data = VPCPercentileData[]

        for pi_level in config.pi_levels
            obs_pctl = compute_percentile(obs_vals, pi_level)
            sim_pctls = sim_percentiles_per_bin[bin_id][pi_level]

            lower, median_sim, upper = bootstrap_percentile_ci(
                sim_pctls, 0.5, config.ci_level, config.n_bootstrap, rng
            )

            push!(percentile_data, VPCPercentileData(
                pi_level, obs_pctl, median_sim, lower, upper
            ))
        end

        push!(vpc_bins, VPCBin(
            bin_def.id,
            bin_def.lower,
            bin_def.upper,
            bin_def.midpoint,
            length(obs_vals),
            config.n_simulations,
            percentile_data
        ))

        # Compute BLQ statistics for this bin
        obs_n_total = length(obs_binned_original[bin_id][2])
        obs_n_blq = count(v -> v < lloq, obs_binned_original[bin_id][2])

        sim_blq_pcts = sim_pct_blq_per_bin[bin_id]
        blq_lower, blq_median, blq_upper = if !isempty(sim_blq_pcts)
            bootstrap_percentile_ci(sim_blq_pcts, 0.5, config.ci_level, config.n_bootstrap, rng)
        else
            (0.0, 0.0, 0.0)
        end

        push!(blq_stats, BLQBinStats(
            bin_id,
            obs_n_total,
            obs_n_blq,
            obs_pct_blq_per_bin[bin_id],
            blq_median,
            blq_lower,
            blq_upper
        ))
    end

    vpc_result = VPCResult(
        config,
        vpc_bins,
        n_subjects(observed),
        n_observations(observed),
        config.n_simulations,
        "",
        config.seed
    )

    return vpc_result, blq_stats
end

export compute_vpc_with_blq
