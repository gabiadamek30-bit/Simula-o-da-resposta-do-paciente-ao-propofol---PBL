# =============================================================================
# TMDD Analysis Utilities
# =============================================================================
#
# Industry-standard pharmacokinetic analysis functions for TMDD models.
# Includes target occupancy, steady-state calculations, and model selection.
#
# References:
# - Gibiansky L, et al. J Pharmacokinet Pharmacodyn. 2008
# - Mager DE, Jusko WJ. J Pharmacokinet Pharmacodyn. 2001
# =============================================================================

using Statistics

export calculate_target_occupancy, calculate_tmdd_steady_state
export calculate_tmdd_half_lives, identify_tmdd_regime
export tmdd_exposure_metrics, compare_tmdd_approximations

# =============================================================================
# Target Occupancy Analysis
# =============================================================================

"""
    calculate_target_occupancy(C, KD) -> Float64

Calculate target occupancy at a given drug concentration.

Target Occupancy = C / (K + C)

where K is KD (equilibrium) or KSS (quasi-steady-state).

# Arguments
- `C`: Drug concentration
- `KD`: Equilibrium dissociation constant or KSS

# Returns
- Fraction of target occupied (0 to 1)

# Example
```julia
# 50% occupancy at KD
calculate_target_occupancy(1.0, 1.0)  # 0.5

# 90% occupancy at 9×KD
calculate_target_occupancy(9.0, 1.0)  # 0.9
```
"""
function calculate_target_occupancy(C::Float64, KD::Float64)::Float64
    if C <= 0.0 || KD <= 0.0
        return 0.0
    end
    return C / (KD + C)
end

"""
Calculate EC50 for target occupancy.
Returns K (by definition, 50% occupancy occurs at C = K).
"""
target_occupancy_EC50(KD::Float64) = KD

"""
Calculate EC90 for target occupancy (90% occupancy).
"""
target_occupancy_EC90(KD::Float64) = 9.0 * KD

"""
Calculate EC99 for target occupancy (99% occupancy).
"""
target_occupancy_EC99(KD::Float64) = 99.0 * KD

# =============================================================================
# Steady-State Calculations
# =============================================================================

"""
    calculate_tmdd_steady_state(params) -> NamedTuple

Calculate steady-state conditions for TMDD model.

Returns:
- R_ss: Free target at steady state (no drug)
- Rtot_ss: Total target at steady state (no drug)
- t_half_target: Target half-life
"""
function calculate_tmdd_steady_state(params::OneCptTMDDParams)
    # At steady state without drug
    R_ss = params.ksyn / params.kdeg
    Rtot_ss = R_ss
    t_half_target = log(2) / params.kdeg

    return (
        R_ss = R_ss,
        Rtot_ss = Rtot_ss,
        t_half_target = t_half_target,
        KSS = params.KSS
    )
end

function calculate_tmdd_steady_state(params::TwoCptTMDDParams)
    R_ss = params.ksyn / params.kdeg
    Rtot_ss = R_ss
    t_half_target = log(2) / params.kdeg
    t_half_linear = log(2) * params.V1 / params.CL

    # Derived PK parameters
    kel = params.CL / params.V1
    k12 = params.Q / params.V1
    k21 = params.Q / params.V2

    # Two-compartment eigenvalues
    a = kel + k12 + k21
    b = kel * k21
    discriminant = sqrt(max(0.0, a^2 - 4b))
    alpha = (a + discriminant) / 2
    beta = (a - discriminant) / 2

    t_half_alpha = log(2) / alpha
    t_half_beta = log(2) / beta

    return (
        R_ss = R_ss,
        Rtot_ss = Rtot_ss,
        t_half_target = t_half_target,
        t_half_linear = t_half_linear,
        t_half_alpha = t_half_alpha,
        t_half_beta = t_half_beta,
        KSS = params.KSS,
        Vss = params.V1 + params.V2
    )
end

function calculate_tmdd_steady_state(params::TwoCptTMDDFcRnParams)
    R_ss = params.ksyn / params.kdeg
    t_half_target = log(2) / params.kdeg

    # FcRn-mediated half-life depends on FR
    # Higher FR → longer half-life
    CL_eff = params.CLup * (1.0 - params.FR)
    t_half_fcrn = log(2) * params.V1 / CL_eff

    return (
        R_ss = R_ss,
        t_half_target = t_half_target,
        t_half_fcrn = t_half_fcrn,
        fraction_recycled = params.FR,
        KSS = params.KSS
    )
end

# =============================================================================
# Half-Life Calculations
# =============================================================================

"""
    calculate_tmdd_half_lives(params) -> NamedTuple

Calculate effective half-lives for TMDD system.

At high concentrations: t1/2 ≈ ln(2) / kel (linear elimination dominates)
At low concentrations: t1/2 much shorter (target-mediated elimination dominates)

Returns multiple half-life estimates for different concentration regimes.
"""
function calculate_tmdd_half_lives(params::OneCptTMDDParams)
    kel = params.CL / params.V

    t_half_linear = log(2) / kel
    t_half_target = log(2) / params.kdeg
    t_half_complex = log(2) / params.kint

    return (
        t_half_linear = t_half_linear,
        t_half_target = t_half_target,
        t_half_complex = t_half_complex,
        CL_linear = params.CL
    )
end

function calculate_tmdd_half_lives(params::TwoCptTMDDParams)
    # Two-compartment: distribution and terminal phases
    V1 = params.V1
    V2 = params.V2
    Q = params.Q
    CL = params.CL

    kel = CL / V1
    k12 = Q / V1
    k21 = Q / V2

    # Eigenvalues of the system (simplified, ignoring target)
    a = kel + k12 + k21
    b = kel * k21
    discriminant = sqrt(max(0.0, a^2 - 4b))

    alpha = (a + discriminant) / 2
    beta = (a - discriminant) / 2

    t_half_alpha = log(2) / alpha  # Distribution phase
    t_half_beta = log(2) / beta    # Terminal phase

    return (
        t_half_alpha = t_half_alpha,
        t_half_beta = t_half_beta,
        t_half_target = log(2) / params.kdeg,
        t_half_complex = log(2) / params.kint,
        CL_linear = CL
    )
end

# =============================================================================
# TMDD Regime Identification
# =============================================================================

"""
    identify_tmdd_regime(C, KSS, R_total) -> Symbol

Identify the TMDD kinetic regime based on drug concentration.

Returns:
- :linear - High concentration, linear kinetics dominate
- :mixed - Intermediate, both linear and TMDD
- :tmdd - Low concentration, TMDD kinetics dominate
"""
function identify_tmdd_regime(C::Float64, KSS::Float64, R_total::Float64)
    if C <= 0.0 || R_total <= 0.0
        return :tmdd
    end

    # Target occupancy
    occupancy = C / (KSS + C)

    # Ratio of free drug to target
    ratio = C / R_total

    if ratio > 10.0 || occupancy > 0.95
        return :linear  # Linear kinetics dominate
    elseif ratio < 0.1 || occupancy < 0.5
        return :tmdd    # TMDD kinetics dominate
    else
        return :mixed   # Mixed kinetics
    end
end

"""
Identify regime over time from simulation result.
"""
function identify_tmdd_regime(result::TMDDSimResult, KSS::Float64)
    regimes = Symbol[]

    if !haskey(result.observations, :conc) || !haskey(result.observations, :R_total)
        return regimes
    end

    conc = result.observations[:conc]
    R_total = result.observations[:R_total]

    for i in 1:length(result.t)
        push!(regimes, identify_tmdd_regime(conc[i], KSS, R_total[i]))
    end

    return regimes
end

# =============================================================================
# Exposure Metrics
# =============================================================================

"""
    tmdd_exposure_metrics(result::TMDDSimResult) -> NamedTuple

Calculate PK exposure metrics from TMDD simulation.

Returns:
- AUC: Area under concentration curve
- AUClast: AUC to last observation
- Cmax: Maximum concentration
- Tmax: Time of maximum concentration
- Ctrough: Concentration at end of simulation
- AUC_occupancy: Area under target occupancy curve
- mean_occupancy: Mean target occupancy
"""
function tmdd_exposure_metrics(result::TMDDSimResult)
    t = result.t
    conc = result.observations[:conc]

    # AUC using trapezoidal rule
    AUC = 0.0
    for i in 2:length(t)
        dt = t[i] - t[i-1]
        AUC += 0.5 * (conc[i] + conc[i-1]) * dt
    end

    # Cmax and Tmax
    Cmax_idx = argmax(conc)
    Cmax = conc[Cmax_idx]
    Tmax = t[Cmax_idx]

    # Ctrough (last concentration)
    Ctrough = conc[end]

    # C at 24h (if available)
    C24 = NaN
    idx_24 = findfirst(x -> x >= 24.0, t)
    if idx_24 !== nothing
        C24 = conc[idx_24]
    end

    # Target occupancy AUC (if available)
    AUC_occupancy = 0.0
    if haskey(result.observations, :target_occupancy)
        occ = result.observations[:target_occupancy]
        for i in 2:length(t)
            dt = t[i] - t[i-1]
            AUC_occupancy += 0.5 * (occ[i] + occ[i-1]) * dt
        end
    end

    # Mean target occupancy
    mean_occupancy = haskey(result.observations, :target_occupancy) ?
        mean(result.observations[:target_occupancy]) : NaN

    # Max target occupancy
    max_occupancy = haskey(result.observations, :target_occupancy) ?
        maximum(result.observations[:target_occupancy]) : NaN

    # Total drug (bound + free) AUC
    AUC_total = 0.0
    if haskey(result.observations, :conc_total)
        conc_total = result.observations[:conc_total]
        for i in 2:length(t)
            dt = t[i] - t[i-1]
            AUC_total += 0.5 * (conc_total[i] + conc_total[i-1]) * dt
        end
    else
        AUC_total = AUC
    end

    return (
        AUC = AUC,
        AUC_total = AUC_total,
        Cmax = Cmax,
        Tmax = Tmax,
        Ctrough = Ctrough,
        C24 = C24,
        AUC_occupancy = AUC_occupancy,
        mean_occupancy = mean_occupancy,
        max_occupancy = max_occupancy
    )
end

# =============================================================================
# Model Comparison/Selection
# =============================================================================

"""
    compare_tmdd_approximations(full_result, approx_result) -> NamedTuple

Compare full TMDD model to an approximation (QSS, QE, MM).

Returns metrics for assessing approximation validity.
"""
function compare_tmdd_approximations(
    full_result::TMDDSimResult,
    approx_result::TMDDSimResult
)
    # Ensure same time points
    t_full = full_result.t
    t_approx = approx_result.t

    if length(t_full) != length(t_approx) || t_full != t_approx
        @warn "Time grids differ between models - comparing common subset"
    end

    # Compare concentrations
    conc_full = full_result.observations[:conc]
    conc_approx = approx_result.observations[:conc]

    n = min(length(conc_full), length(conc_approx))

    # Relative error
    rel_errors = Float64[]
    for i in 1:n
        if conc_full[i] > 1e-10
            push!(rel_errors, abs(conc_full[i] - conc_approx[i]) / conc_full[i])
        end
    end

    max_rel_error = isempty(rel_errors) ? 0.0 : maximum(rel_errors)
    mean_rel_error = isempty(rel_errors) ? 0.0 : mean(rel_errors)

    # Compare AUC
    metrics_full = tmdd_exposure_metrics(full_result)
    metrics_approx = tmdd_exposure_metrics(approx_result)

    AUC_rel_error = abs(metrics_full.AUC - metrics_approx.AUC) / max(metrics_full.AUC, 1e-10)
    Cmax_rel_error = abs(metrics_full.Cmax - metrics_approx.Cmax) / max(metrics_full.Cmax, 1e-10)

    # Compare target occupancy if available
    occupancy_rel_error = NaN
    if haskey(full_result.observations, :target_occupancy) &&
       haskey(approx_result.observations, :target_occupancy)
        occ_full = full_result.observations[:target_occupancy]
        occ_approx = approx_result.observations[:target_occupancy]

        occ_errors = Float64[]
        for i in 1:min(length(occ_full), length(occ_approx))
            if occ_full[i] > 0.01
                push!(occ_errors, abs(occ_full[i] - occ_approx[i]) / occ_full[i])
            end
        end
        occupancy_rel_error = isempty(occ_errors) ? 0.0 : mean(occ_errors)
    end

    # Recommendation
    if max_rel_error < 0.05 && AUC_rel_error < 0.05
        recommendation = :excellent
    elseif max_rel_error < 0.15 && AUC_rel_error < 0.10
        recommendation = :acceptable
    elseif max_rel_error < 0.25
        recommendation = :marginal
    else
        recommendation = :use_full_model
    end

    return (
        max_rel_error = max_rel_error,
        mean_rel_error = mean_rel_error,
        AUC_rel_error = AUC_rel_error,
        Cmax_rel_error = Cmax_rel_error,
        occupancy_rel_error = occupancy_rel_error,
        recommendation = recommendation
    )
end

# =============================================================================
# Dose Selection Guidance
# =============================================================================

"""
    target_trough_dose(params, target_occupancy, dosing_interval) -> Float64

Estimate dose needed to achieve target trough occupancy.

Simplified calculation based on steady-state assumptions.
"""
function target_trough_dose(
    params::TwoCptTMDDParams,
    target_occupancy::Float64,
    dosing_interval::Float64
)
    # Target concentration for desired occupancy
    # occupancy = C / (KSS + C)
    # C = KSS * occupancy / (1 - occupancy)
    C_target = params.KSS * target_occupancy / (1.0 - target_occupancy)

    # Rough estimate of kel
    kel = params.CL / params.V1

    # C_trough = C_max * exp(-kel * tau)
    # We want C_trough = C_target
    # C_max = C_target * exp(kel * tau)
    C_max_needed = C_target * exp(kel * dosing_interval)

    # Dose ≈ Cmax * V (very simplified, ignoring TMDD during infusion)
    dose = C_max_needed * params.V1

    return dose
end

"""
Calculate loading dose to achieve rapid target saturation.
"""
function loading_dose(params::TwoCptTMDDParams, target_occupancy::Float64=0.90)
    # Concentration for target occupancy
    C_target = params.KSS * target_occupancy / (1.0 - target_occupancy)

    # Need to saturate target
    R0_conc = params.R0 / params.V1

    # Loading dose to achieve C_target and saturate R0
    Vss = params.V1 + params.V2
    dose = C_target * Vss + params.R0  # Extra to saturate target

    return dose
end

# =============================================================================
# Time Above Target Metrics
# =============================================================================

"""
    time_above_occupancy(result::TMDDSimResult, threshold) -> Float64

Calculate time (hours) above a target occupancy threshold.
"""
function time_above_occupancy(result::TMDDSimResult, threshold::Float64)
    if !haskey(result.observations, :target_occupancy)
        return NaN
    end

    t = result.t
    occ = result.observations[:target_occupancy]

    time_above = 0.0
    for i in 2:length(t)
        if occ[i] >= threshold && occ[i-1] >= threshold
            time_above += t[i] - t[i-1]
        elseif occ[i] >= threshold || occ[i-1] >= threshold
            # Linear interpolation for crossing point
            if occ[i] >= threshold
                # Crossing upward
                frac = (threshold - occ[i-1]) / (occ[i] - occ[i-1])
                time_above += (1 - frac) * (t[i] - t[i-1])
            else
                # Crossing downward
                frac = (threshold - occ[i-1]) / (occ[i] - occ[i-1])
                time_above += frac * (t[i] - t[i-1])
            end
        end
    end

    return time_above
end

"""
    time_above_concentration(result::TMDDSimResult, threshold) -> Float64

Calculate time (hours) above a concentration threshold.
"""
function time_above_concentration(result::TMDDSimResult, threshold::Float64)
    t = result.t
    conc = result.observations[:conc]

    time_above = 0.0
    for i in 2:length(t)
        if conc[i] >= threshold && conc[i-1] >= threshold
            time_above += t[i] - t[i-1]
        elseif conc[i] >= threshold || conc[i-1] >= threshold
            # Linear interpolation
            if conc[i] >= threshold
                frac = (threshold - conc[i-1]) / (conc[i] - conc[i-1])
                time_above += (1 - frac) * (t[i] - t[i-1])
            else
                frac = (threshold - conc[i-1]) / (conc[i] - conc[i-1])
                time_above += frac * (t[i] - t[i-1])
            end
        end
    end

    return time_above
end

# =============================================================================
# Dose Escalation Analysis
# =============================================================================

"""
    dose_proportionality(results::Vector{TMDDSimResult}, doses::Vector{Float64})

Assess dose proportionality in TMDD system.

TMDD typically shows more-than-dose-proportional increase in AUC at low doses
(due to target saturation) and dose-proportional at high doses.
"""
function dose_proportionality(results::Vector{TMDDSimResult}, doses::Vector{Float64})
    if length(results) != length(doses)
        error("Number of results must match number of doses")
    end

    n = length(doses)
    aucs = [tmdd_exposure_metrics(r).AUC for r in results]
    cmax_values = [tmdd_exposure_metrics(r).Cmax for r in results]

    # Dose-normalized AUC
    auc_dn = aucs ./ doses

    # Power model fit: AUC = a * Dose^b
    # If b ≈ 1, dose proportional
    log_doses = log.(doses)
    log_aucs = log.(aucs)

    # Simple linear regression for slope (b)
    mean_x = mean(log_doses)
    mean_y = mean(log_aucs)
    b_auc = sum((log_doses .- mean_x) .* (log_aucs .- mean_y)) /
            sum((log_doses .- mean_x).^2)

    # For Cmax
    log_cmax = log.(cmax_values)
    mean_y_cmax = mean(log_cmax)
    b_cmax = sum((log_doses .- mean_x) .* (log_cmax .- mean_y_cmax)) /
             sum((log_doses .- mean_x).^2)

    # Interpretation
    if b_auc > 1.2
        auc_interpretation = :more_than_proportional
    elseif b_auc < 0.8
        auc_interpretation = :less_than_proportional
    else
        auc_interpretation = :approximately_proportional
    end

    return (
        doses = doses,
        aucs = aucs,
        auc_dn = auc_dn,
        cmax_values = cmax_values,
        power_auc = b_auc,
        power_cmax = b_cmax,
        auc_interpretation = auc_interpretation
    )
end

export target_occupancy_EC50, target_occupancy_EC90, target_occupancy_EC99
export target_trough_dose, loading_dose
export time_above_occupancy, time_above_concentration
export dose_proportionality
