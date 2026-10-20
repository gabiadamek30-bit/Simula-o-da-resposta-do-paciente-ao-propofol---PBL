# PK Parameters
# Secondary PK parameters derived from NCA metrics

export nca_half_life, nca_mrt, nca_cl_f, nca_vz_f, nca_vss
export nca_cl, nca_vz, nca_mrt_iv, nca_cl_ss
export nca_vss_from_aumc, nca_vc
export nca_mean_absorption_time, nca_bioavailability
export nca_c0_backextrap, nca_c0_from_regression

# =============================================================================
# Terminal Half-Life
# =============================================================================

"""
    nca_half_life(lambda_z)

Calculate terminal half-life from lambda_z.

t1/2 = ln(2) / λz

# Arguments
- `lambda_z::Float64`: Terminal elimination rate constant

# Returns
- `Float64`: Terminal half-life
"""
function nca_half_life(lambda_z::Float64)
    @assert lambda_z > 0.0 "lambda_z must be positive"
    return log(2.0) / lambda_z
end

# =============================================================================
# Mean Residence Time
# =============================================================================

"""
    nca_mrt(aumc_0_inf, auc_0_inf; route=:extravascular)

Calculate Mean Residence Time (MRT).

For extravascular administration:
    MRT = AUMC0-∞ / AUC0-∞

For IV bolus:
    MRTiv = AUMC0-∞ / AUC0-∞

For IV infusion:
    MRTinf = AUMC0-∞ / AUC0-∞ - Tinf/2

# Arguments
- `aumc_0_inf::Float64`: AUMC from 0 to infinity
- `auc_0_inf::Float64`: AUC from 0 to infinity
- `route::Symbol`: Administration route (:extravascular, :iv_bolus, :iv_infusion)
- `t_inf::Float64`: Infusion duration (only for :iv_infusion)

# Returns
- `Float64`: Mean residence time
"""
function nca_mrt(
    aumc_0_inf::Float64,
    auc_0_inf::Float64;
    route::Symbol = :extravascular,
    t_inf::Float64 = 0.0
)
    @assert auc_0_inf > 0.0 "AUC0-inf must be positive"
    @assert aumc_0_inf >= 0.0 "AUMC0-inf must be non-negative"

    mrt = aumc_0_inf / auc_0_inf

    if route == :iv_infusion && t_inf > 0.0
        mrt = mrt - t_inf / 2.0
    end

    return mrt
end

"""
    nca_mrt_iv(mrt_extravascular, absorption_time)

Estimate MRT for IV administration from extravascular MRT.

MRTiv = MRText - MAT

where MAT (Mean Absorption Time) is estimated from absorption_time.

# Arguments
- `mrt_extravascular::Float64`: MRT from extravascular administration
- `absorption_time::Float64`: Mean absorption time estimate

# Returns
- `Float64`: Estimated IV MRT
"""
function nca_mrt_iv(mrt_extravascular::Float64, absorption_time::Float64)
    return mrt_extravascular - absorption_time
end

# =============================================================================
# Clearance
# =============================================================================

"""
    nca_cl_f(dose, auc_0_inf)

Calculate apparent clearance (CL/F) for extravascular administration.

CL/F = Dose / AUC0-∞

# Arguments
- `dose::Float64`: Administered dose
- `auc_0_inf::Float64`: AUC from 0 to infinity

# Returns
- `Float64`: Apparent clearance
"""
function nca_cl_f(dose::Float64, auc_0_inf::Float64)
    @assert dose > 0.0 "Dose must be positive"
    @assert auc_0_inf > 0.0 "AUC0-inf must be positive"

    return dose / auc_0_inf
end

"""
    nca_cl(dose, auc_0_inf)

Calculate systemic clearance (CL) for IV administration.

CL = Dose / AUC0-∞

# Arguments
- `dose::Float64`: Administered dose
- `auc_0_inf::Float64`: AUC from 0 to infinity

# Returns
- `Float64`: Systemic clearance
"""
function nca_cl(dose::Float64, auc_0_inf::Float64)
    return nca_cl_f(dose, auc_0_inf)
end

"""
    nca_cl_ss(dose, auc_0_tau)

Calculate clearance at steady state.

CLss = Dose / AUC0-τ

# Arguments
- `dose::Float64`: Administered dose per interval
- `auc_0_tau::Float64`: AUC over dosing interval

# Returns
- `Float64`: Steady-state clearance
"""
function nca_cl_ss(dose::Float64, auc_0_tau::Float64)
    @assert dose > 0.0 "Dose must be positive"
    @assert auc_0_tau > 0.0 "AUC0-tau must be positive"

    return dose / auc_0_tau
end

# =============================================================================
# Volume of Distribution
# =============================================================================

"""
    nca_vz_f(dose, lambda_z, auc_0_inf)

Calculate apparent volume of distribution (Vz/F) for extravascular administration.

Vz/F = Dose / (λz × AUC0-∞)

Also known as terminal volume of distribution.

# Arguments
- `dose::Float64`: Administered dose
- `lambda_z::Float64`: Terminal elimination rate constant
- `auc_0_inf::Float64`: AUC from 0 to infinity

# Returns
- `Float64`: Apparent terminal volume of distribution
"""
function nca_vz_f(dose::Float64, lambda_z::Float64, auc_0_inf::Float64)
    @assert dose > 0.0 "Dose must be positive"
    @assert lambda_z > 0.0 "lambda_z must be positive"
    @assert auc_0_inf > 0.0 "AUC0-inf must be positive"

    return dose / (lambda_z * auc_0_inf)
end

"""
    nca_vz(dose, lambda_z, auc_0_inf)

Calculate terminal volume of distribution (Vz) for IV administration.

Vz = Dose / (λz × AUC0-∞)

# Arguments
- `dose::Float64`: Administered dose
- `lambda_z::Float64`: Terminal elimination rate constant
- `auc_0_inf::Float64`: AUC from 0 to infinity

# Returns
- `Float64`: Terminal volume of distribution
"""
function nca_vz(dose::Float64, lambda_z::Float64, auc_0_inf::Float64)
    return nca_vz_f(dose, lambda_z, auc_0_inf)
end

"""
    nca_vss(cl, mrt)

Calculate volume of distribution at steady state (Vss).

Vss = CL × MRT

# Arguments
- `cl::Float64`: Clearance (CL or CL/F)
- `mrt::Float64`: Mean residence time

# Returns
- `Float64`: Volume at steady state
"""
function nca_vss(cl::Float64, mrt::Float64)
    @assert cl > 0.0 "Clearance must be positive"
    @assert mrt > 0.0 "MRT must be positive"

    return cl * mrt
end

"""
    nca_vss_from_aumc(dose, auc_0_inf, aumc_0_inf)

Calculate Vss directly from moment curves.

Vss = Dose × AUMC0-∞ / (AUC0-∞)²

# Arguments
- `dose::Float64`: Administered dose
- `auc_0_inf::Float64`: AUC from 0 to infinity
- `aumc_0_inf::Float64`: AUMC from 0 to infinity

# Returns
- `Float64`: Volume at steady state
"""
function nca_vss_from_aumc(dose::Float64, auc_0_inf::Float64, aumc_0_inf::Float64)
    @assert dose > 0.0 "Dose must be positive"
    @assert auc_0_inf > 0.0 "AUC0-inf must be positive"
    @assert aumc_0_inf >= 0.0 "AUMC0-inf must be non-negative"

    return dose * aumc_0_inf / (auc_0_inf^2)
end

"""
    nca_vc(dose, c0)

Calculate central volume of distribution from back-extrapolated C0.

Vc = Dose / C0

For IV bolus administration.

# Arguments
- `dose::Float64`: Administered dose
- `c0::Float64`: Back-extrapolated concentration at time 0

# Returns
- `Float64`: Central volume of distribution
"""
function nca_vc(dose::Float64, c0::Float64)
    @assert dose > 0.0 "Dose must be positive"
    @assert c0 > 0.0 "C0 must be positive"

    return dose / c0
end

# =============================================================================
# Absorption Parameters
# =============================================================================

"""
    nca_mean_absorption_time(mrt_po, mrt_iv)

Calculate Mean Absorption Time (MAT).

MAT = MRTpo - MRTiv

# Arguments
- `mrt_po::Float64`: MRT from oral/extravascular administration
- `mrt_iv::Float64`: MRT from IV administration

# Returns
- `Float64`: Mean absorption time
"""
function nca_mean_absorption_time(mrt_po::Float64, mrt_iv::Float64)
    @assert mrt_po >= mrt_iv "MRT(po) should be >= MRT(iv)"
    return mrt_po - mrt_iv
end

"""
    nca_bioavailability(auc_test, dose_test, auc_reference, dose_reference)

Calculate relative bioavailability (F).

F = (AUCtest / Dosetest) / (AUCref / Doseref)

# Arguments
- `auc_test::Float64`: AUC of test formulation
- `dose_test::Float64`: Dose of test formulation
- `auc_reference::Float64`: AUC of reference formulation
- `dose_reference::Float64`: Dose of reference formulation

# Returns
- `Float64`: Relative bioavailability
"""
function nca_bioavailability(
    auc_test::Float64,
    dose_test::Float64,
    auc_reference::Float64,
    dose_reference::Float64
)
    @assert dose_test > 0.0 "Test dose must be positive"
    @assert dose_reference > 0.0 "Reference dose must be positive"
    @assert auc_test >= 0.0 "Test AUC must be non-negative"
    @assert auc_reference > 0.0 "Reference AUC must be positive"

    return (auc_test / dose_test) / (auc_reference / dose_reference)
end

# =============================================================================
# Back-Extrapolation
# =============================================================================

"""
    nca_c0_backextrap(t, c, lambda_z)

Back-extrapolate C0 for IV bolus administration.

C0 = Clast × exp(λz × tlast)

Uses the lambda_z regression parameters.

# Arguments
- `t::Vector{Float64}`: Time points (from terminal phase regression)
- `c::Vector{Float64}`: Concentration values
- `lambda_z::Float64`: Terminal elimination rate constant

# Returns
- `Float64`: Back-extrapolated C0
"""
function nca_c0_backextrap(t::Vector{Float64}, c::Vector{Float64}, lambda_z::Float64)
    @assert lambda_z > 0.0 "lambda_z must be positive"
    @assert !isempty(t) "Time vector cannot be empty"

    # Use first point to back-extrapolate
    return c[1] * exp(lambda_z * t[1])
end

"""
    nca_c0_from_regression(intercept)

Get C0 from lambda_z regression intercept.

For log-linear regression: ln(C) = intercept - λz × t
C0 = exp(intercept)

# Arguments
- `intercept::Float64`: Y-intercept from log-linear regression

# Returns
- `Float64`: Back-extrapolated C0
"""
function nca_c0_from_regression(intercept::Float64)
    return exp(intercept)
end
