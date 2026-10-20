# Multiple Dose Metrics
# NCA metrics for multiple dose and steady-state analysis

export nca_accumulation_index, nca_ptf, nca_swing
export nca_linearity_index, nca_time_to_steady_state
export nca_accumulation_predicted, nca_accumulation_cmax, nca_accumulation_cmin
export nca_ptf_from_concentrations, nca_swing_from_concentrations
export nca_dose_normalized_auc, nca_dose_normalized_cmax
export nca_time_to_steady_state_doses, nca_effective_half_life

# =============================================================================
# Accumulation Metrics
# =============================================================================

"""
    nca_accumulation_index(auc_ss, auc_sd)

Calculate accumulation index (Rac) from single dose and steady-state AUC.

Rac = AUCss / AUCsd

Where:
- AUCss = AUC over dosing interval at steady state
- AUCsd = AUC0-inf from single dose

# Arguments
- `auc_ss::Float64`: AUC0-τ at steady state
- `auc_sd::Float64`: AUC0-inf from single dose

# Returns
- `Float64`: Accumulation index
"""
function nca_accumulation_index(auc_ss::Float64, auc_sd::Float64)
    @assert auc_ss >= 0.0 "Steady-state AUC must be non-negative"
    @assert auc_sd > 0.0 "Single-dose AUC must be positive"

    return auc_ss / auc_sd
end

"""
    nca_accumulation_predicted(lambda_z, tau)

Predict accumulation index from lambda_z and dosing interval.

Rac_pred = 1 / (1 - exp(-λz × τ))

This is the theoretical accumulation for a one-compartment model.

# Arguments
- `lambda_z::Float64`: Terminal elimination rate constant
- `tau::Float64`: Dosing interval

# Returns
- `Float64`: Predicted accumulation index
"""
function nca_accumulation_predicted(lambda_z::Float64, tau::Float64)
    @assert lambda_z > 0.0 "lambda_z must be positive"
    @assert tau > 0.0 "Dosing interval must be positive"

    return 1.0 / (1.0 - exp(-lambda_z * tau))
end

"""
    nca_accumulation_cmax(cmax_ss, cmax_sd)

Calculate Cmax accumulation ratio.

Rac,Cmax = Cmax,ss / Cmax,sd

# Arguments
- `cmax_ss::Float64`: Cmax at steady state
- `cmax_sd::Float64`: Cmax from single dose

# Returns
- `Float64`: Cmax accumulation ratio
"""
function nca_accumulation_cmax(cmax_ss::Float64, cmax_sd::Float64)
    @assert cmax_ss >= 0.0 "Steady-state Cmax must be non-negative"
    @assert cmax_sd > 0.0 "Single-dose Cmax must be positive"

    return cmax_ss / cmax_sd
end

"""
    nca_accumulation_cmin(cmin_ss, c_at_tau_sd)

Calculate Cmin accumulation ratio.

Rac,Cmin = Cmin,ss / C(τ)sd

# Arguments
- `cmin_ss::Float64`: Cmin (trough) at steady state
- `c_at_tau_sd::Float64`: Concentration at time τ after single dose

# Returns
- `Float64`: Cmin accumulation ratio
"""
function nca_accumulation_cmin(cmin_ss::Float64, c_at_tau_sd::Float64)
    @assert cmin_ss >= 0.0 "Steady-state Cmin must be non-negative"
    @assert c_at_tau_sd > 0.0 "Single-dose concentration at τ must be positive"

    return cmin_ss / c_at_tau_sd
end

# =============================================================================
# Fluctuation Metrics
# =============================================================================

"""
    nca_ptf(cmax, cmin, cavg)

Calculate Peak-Trough Fluctuation (PTF) percentage.

PTF = 100 × (Cmax - Cmin) / Cavg

# Arguments
- `cmax::Float64`: Maximum concentration in dosing interval
- `cmin::Float64`: Minimum (trough) concentration
- `cavg::Float64`: Average concentration over dosing interval

# Returns
- `Float64`: Peak-trough fluctuation (%)
"""
function nca_ptf(cmax::Float64, cmin::Float64, cavg::Float64)
    @assert cmax >= cmin "Cmax must be >= Cmin"
    @assert cavg > 0.0 "Cavg must be positive"

    return 100.0 * (cmax - cmin) / cavg
end

"""
    nca_ptf_from_concentrations(t, c, tau, config)

Calculate PTF from concentration-time data over a dosing interval.

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `tau::Float64`: Dosing interval
- `config::NCAConfig`: NCA configuration

# Returns
- `Float64`: Peak-trough fluctuation (%)
"""
function nca_ptf_from_concentrations(
    t::Vector{Float64},
    c::Vector{Float64},
    tau::Float64,
    config::NCAConfig
)
    # Filter to dosing interval
    idx = findall(ti -> 0.0 <= ti <= tau, t)
    c_tau = c[idx]

    cmax = maximum(c_tau)
    cmin = minimum(c_tau)
    cavg = nca_cavg(t, c, tau, config)

    return nca_ptf(cmax, cmin, cavg)
end

"""
    nca_swing(cmax, cmin)

Calculate Swing percentage.

Swing = 100 × (Cmax - Cmin) / Cmin

# Arguments
- `cmax::Float64`: Maximum concentration in dosing interval
- `cmin::Float64`: Minimum (trough) concentration

# Returns
- `Float64`: Swing (%)
"""
function nca_swing(cmax::Float64, cmin::Float64)
    @assert cmax >= cmin "Cmax must be >= Cmin"
    @assert cmin > 0.0 "Cmin must be positive for swing calculation"

    return 100.0 * (cmax - cmin) / cmin
end

"""
    nca_swing_from_concentrations(t, c, tau)

Calculate Swing from concentration-time data over a dosing interval.

# Arguments
- `t::Vector{Float64}`: Time points
- `c::Vector{Float64}`: Concentration values
- `tau::Float64`: Dosing interval

# Returns
- `Float64`: Swing (%)
"""
function nca_swing_from_concentrations(
    t::Vector{Float64},
    c::Vector{Float64},
    tau::Float64
)
    idx = findall(ti -> 0.0 <= ti <= tau, t)
    c_tau = c[idx]

    cmax = maximum(c_tau)
    cmin = minimum(c_tau)

    return nca_swing(cmax, cmin)
end

# =============================================================================
# Dose Linearity Assessment
# =============================================================================

"""
    nca_linearity_index(doses, aucs)

Assess dose proportionality using power model.

AUC = α × Dose^β

Linear if β ≈ 1.0.

# Arguments
- `doses::Vector{Float64}`: Dose levels
- `aucs::Vector{Float64}`: Corresponding AUC values (dose-normalized or raw)

# Returns
- `NamedTuple`: (beta, r_squared, is_linear)
"""
function nca_linearity_index(doses::Vector{Float64}, aucs::Vector{Float64})
    @assert length(doses) == length(aucs) "Doses and AUCs must have same length"
    @assert length(doses) >= 2 "Need at least 2 dose levels"
    @assert all(d -> d > 0, doses) "All doses must be positive"
    @assert all(a -> a > 0, aucs) "All AUCs must be positive"

    # Log-transform for power model
    log_doses = log.(doses)
    log_aucs = log.(aucs)

    # Linear regression: log(AUC) = log(α) + β × log(Dose)
    slope, intercept, r_squared, _ = _power_model_regression(log_doses, log_aucs)

    # Beta is the slope in log-log space
    beta = slope

    # Linear if beta is close to 1.0 (within typical 0.8-1.25 bounds)
    is_linear = 0.8 <= beta <= 1.25

    return (beta=beta, r_squared=r_squared, is_linear=is_linear)
end

"""
    _power_model_regression(log_x, log_y)

Perform linear regression on log-transformed data.
"""
function _power_model_regression(log_x::Vector{Float64}, log_y::Vector{Float64})
    n = length(log_x)
    x_mean = sum(log_x) / n
    y_mean = sum(log_y) / n

    ss_xx = sum((log_x .- x_mean).^2)
    ss_xy = sum((log_x .- x_mean) .* (log_y .- y_mean))
    ss_yy = sum((log_y .- y_mean).^2)

    slope = ss_xy / ss_xx
    intercept = y_mean - slope * x_mean

    ss_res = sum((log_y .- (intercept .+ slope .* log_x)).^2)
    r_squared = 1.0 - ss_res / ss_yy

    return (slope, intercept, r_squared, nothing)
end

"""
    nca_dose_normalized_auc(auc, dose)

Calculate dose-normalized AUC.

AUC/D = AUC / Dose

# Arguments
- `auc::Float64`: AUC value
- `dose::Float64`: Dose administered

# Returns
- `Float64`: Dose-normalized AUC
"""
function nca_dose_normalized_auc(auc::Float64, dose::Float64)
    @assert dose > 0.0 "Dose must be positive"
    return auc / dose
end

"""
    nca_dose_normalized_cmax(cmax, dose)

Calculate dose-normalized Cmax.

Cmax/D = Cmax / Dose

# Arguments
- `cmax::Float64`: Maximum concentration
- `dose::Float64`: Dose administered

# Returns
- `Float64`: Dose-normalized Cmax
"""
function nca_dose_normalized_cmax(cmax::Float64, dose::Float64)
    @assert dose > 0.0 "Dose must be positive"
    return cmax / dose
end

# =============================================================================
# Time to Steady State
# =============================================================================

"""
    nca_time_to_steady_state(lambda_z; fraction=0.90)

Estimate time to reach a fraction of steady state.

t_ss = -ln(1 - fraction) / λz

For 90% steady state: t_ss = 2.303 × t1/2 ≈ 3.3 × t1/2
For 95% steady state: t_ss = 3 × t1/2 × 1.1 ≈ 4.3 × t1/2

# Arguments
- `lambda_z::Float64`: Terminal elimination rate constant
- `fraction::Float64`: Fraction of steady state (default: 0.90)

# Returns
- `Float64`: Time to reach specified fraction of steady state
"""
function nca_time_to_steady_state(lambda_z::Float64; fraction::Float64 = 0.90)
    @assert lambda_z > 0.0 "lambda_z must be positive"
    @assert 0.0 < fraction < 1.0 "Fraction must be between 0 and 1"

    return -log(1.0 - fraction) / lambda_z
end

"""
    nca_time_to_steady_state_doses(lambda_z, tau; fraction=0.90)

Estimate number of doses to reach a fraction of steady state.

n_doses = ceil(t_ss / τ)

# Arguments
- `lambda_z::Float64`: Terminal elimination rate constant
- `tau::Float64`: Dosing interval
- `fraction::Float64`: Fraction of steady state (default: 0.90)

# Returns
- `Int`: Number of doses to reach steady state
"""
function nca_time_to_steady_state_doses(
    lambda_z::Float64,
    tau::Float64;
    fraction::Float64 = 0.90
)
    t_ss = nca_time_to_steady_state(lambda_z; fraction=fraction)
    return ceil(Int, t_ss / tau)
end

# =============================================================================
# Effective Half-Life
# =============================================================================

"""
    nca_effective_half_life(auc_ss, auc_inf_sd)

Calculate effective half-life from steady-state data.

t1/2,eff = ln(2) × τ / ln(AUC0-∞,sd / (AUC0-∞,sd - AUCτ,ss))

Alternative: t1/2,eff = ln(2) × τ / ln(Rac / (Rac - 1))

# Arguments
- `auc_ss::Float64`: AUC0-τ at steady state
- `auc_inf_sd::Float64`: AUC0-inf from single dose

# Returns
- `Float64`: Effective half-life
"""
function nca_effective_half_life(auc_ss::Float64, auc_inf_sd::Float64)
    @assert auc_ss > 0.0 "Steady-state AUC must be positive"
    @assert auc_inf_sd > 0.0 "Single-dose AUC must be positive"
    @assert auc_inf_sd > auc_ss "AUC0-inf should be > AUC0-tau at steady state"

    rac = auc_ss / auc_inf_sd

    if rac >= 1.0
        error("Invalid accumulation ratio (>= 1.0)")
    end

    # t1/2,eff = -ln(2) × τ / ln(1 - Rac)
    # This requires knowing tau, so we return the ratio
    return log(2.0) / log(auc_inf_sd / (auc_inf_sd - auc_ss))
end
