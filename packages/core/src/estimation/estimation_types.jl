# Estimation Types for NLME Parameter Estimation
# Defines methods, configurations, and results

using LinearAlgebra

export EstimationMethod, FOCEIMethod, SAEMMethod, LaplacianMethod
export EstimationConfig, IndividualEstimate, EstimationResult
export OmegaStructure, DiagonalOmega, BlockOmega
export CovariateOnIIV, EstimationIOVSpec, OccasionData
export EstimationBLQMethod, BLQ_M1_DISCARD, BLQ_M2_IMPUTE_HALF, BLQ_M2_IMPUTE_ZERO, BLQ_M3_LIKELIHOOD
export BLQConfig, BLQSummary

# ------------------------------------------------------------------
# Estimation Methods
# ------------------------------------------------------------------

"""
Abstract type for estimation methods.
"""
abstract type EstimationMethod end

"""
First-Order Conditional Estimation with Interaction (FOCE-I).

This is the most commonly used method in NONMEM. It linearizes the model
around the conditional mode of the random effects (eta).

Fields:
- max_inner_iter: Maximum iterations for finding eta mode (default: 50)
- inner_tol: Convergence tolerance for inner optimization (default: 1e-6)
- centered: Use centered FO approximation at η=0 instead of η=η* (default: false)
  - false (FOCE-I): Linearize at conditional mode η* (includes interaction)
  - true (FOCE): Linearize at population mean η=0 (no interaction)
- compute_robust_se: Compute sandwich estimator for robust SEs (default: true)
"""
struct FOCEIMethod <: EstimationMethod
    max_inner_iter::Int
    inner_tol::Float64
    centered::Bool
    compute_robust_se::Bool

    function FOCEIMethod(;
        max_inner_iter::Int=50,
        inner_tol::Float64=1e-6,
        centered::Bool=false,
        compute_robust_se::Bool=true
    )
        new(max_inner_iter, inner_tol, centered, compute_robust_se)
    end
end

"""
Stochastic Approximation Expectation Maximization (SAEM).

SAEM uses MCMC sampling to approximate the E-step and stochastic
approximation for the M-step. Often more robust than FOCE-I for
complex models.

Fields:
- n_burn: Number of burn-in iterations (default: 300)
- n_iter: Number of main iterations after burn-in (default: 200)
- n_chains: Number of MCMC chains per subject (default: 3)
- n_mcmc_steps: Number of MCMC steps per E-step iteration (default: 100)
- step_size_schedule: Step size decay schedule (:harmonic or :constant)
- adapt_proposal: Whether to adapt MCMC proposal variance (default: true)
- target_acceptance: Target MCMC acceptance rate for adaptation (default: 0.234)
- adaptation_interval: How often to adapt proposal (default: 50 iterations)
- track_diagnostics: Whether to compute convergence diagnostics (default: true)
- use_all_chains: Whether to use all chains or just first chain (default: true)
- parallel_chains: Run MCMC chains in parallel using threads (default: true)
  This provides ~Nx speedup where N = min(n_chains, n_threads).
  Set to false for debugging or when running subjects in parallel already.
"""
struct SAEMMethod <: EstimationMethod
    n_burn::Int
    n_iter::Int
    n_chains::Int
    n_mcmc_steps::Int
    step_size_schedule::Symbol
    adapt_proposal::Bool
    target_acceptance::Float64
    adaptation_interval::Int
    track_diagnostics::Bool
    use_all_chains::Bool
    parallel_chains::Bool

    function SAEMMethod(;
        n_burn::Int=300,
        n_iter::Int=200,
        n_chains::Int=3,
        n_mcmc_steps::Int=100,
        step_size_schedule::Symbol=:harmonic,
        adapt_proposal::Bool=true,
        target_acceptance::Float64=0.234,
        adaptation_interval::Int=50,
        track_diagnostics::Bool=true,
        use_all_chains::Bool=true,
        parallel_chains::Bool=true
    )
        @assert step_size_schedule in (:harmonic, :constant) "Invalid step size schedule"
        @assert n_mcmc_steps >= 10 "n_mcmc_steps must be at least 10"
        @assert 0.0 < target_acceptance < 1.0 "target_acceptance must be in (0, 1)"
        @assert adaptation_interval >= 1 "adaptation_interval must be positive"
        new(n_burn, n_iter, n_chains, n_mcmc_steps, step_size_schedule,
            adapt_proposal, target_acceptance, adaptation_interval,
            track_diagnostics, use_all_chains, parallel_chains)
    end
end

"""
Laplacian Estimation Method.

The simplest approximation method. Uses the mode of the posterior
for random effects and the Laplace approximation for the marginal
likelihood. Less accurate than FOCE-I but computationally faster.
"""
struct LaplacianMethod <: EstimationMethod
    max_inner_iter::Int
    inner_tol::Float64

    function LaplacianMethod(;
        max_inner_iter::Int=50,
        inner_tol::Float64=1e-6
    )
        new(max_inner_iter, inner_tol)
    end
end

# ------------------------------------------------------------------
# Omega Structure
# ------------------------------------------------------------------

"""
Abstract type for omega matrix structure.
"""
abstract type OmegaStructure end

"""
Diagonal omega matrix - no correlations between random effects.
"""
struct DiagonalOmega <: OmegaStructure end

"""
Block diagonal omega matrix - allows correlations within blocks.

Fields:
- block_sizes: Vector of block sizes (e.g., [2, 1] for a 2x2 block and 1x1 block)
"""
struct BlockOmega <: OmegaStructure
    block_sizes::Vector{Int}
end

"""
Full omega matrix - all random effects correlated.
"""
struct FullOmega <: OmegaStructure end

export FullOmega

# ------------------------------------------------------------------
# Covariate Effects on IIV (Inter-Individual Variability)
# ------------------------------------------------------------------

"""
Specification for covariate effects on inter-individual variability.

This allows variance parameters to depend on covariates:
    Var(η_i) = ω² * exp(θ_cov * (COV_i - COV_ref))

Example: Weight effect on CL variability
    CovariateOnIIV(:CL, :WT, 70.0)  # Weight centered on 70 kg

Fields:
- eta_name: Name of the random effect (e.g., :CL, :V)
- covariate_name: Name of the covariate (e.g., :WT, :CRCL)
- reference_value: Reference/centering value for the covariate
- effect_type: :exponential (default) or :linear
"""
struct CovariateOnIIV
    eta_name::Symbol
    covariate_name::Symbol
    reference_value::Float64
    effect_type::Symbol  # :exponential or :linear

    function CovariateOnIIV(
        eta_name::Symbol,
        covariate_name::Symbol,
        reference_value::Float64=0.0;
        effect_type::Symbol=:exponential
    )
        @assert effect_type in [:exponential, :linear] "effect_type must be :exponential or :linear"
        new(eta_name, covariate_name, reference_value, effect_type)
    end
end

# ------------------------------------------------------------------
# Inter-Occasion Variability (IOV)
# ------------------------------------------------------------------

"""
Specification for inter-occasion variability in estimation.

IOV allows random effects to vary across occasions within a subject:
    η_total = η_IIV + η_IOV_occasion

Example: Drug absorption varies between study visits
    EstimationIOVSpec(:ka, [:V1, :V2, :V3])  # IOV on ka across 3 visits

Fields:
- eta_name: Name of the random effect with IOV
- occasion_names: Names/labels for each occasion
- omega_iov: Variance of IOV (separate from IIV variance)
"""
struct EstimationIOVSpec
    eta_name::Symbol
    occasion_names::Vector{Symbol}
    omega_iov::Float64

    function EstimationIOVSpec(
        eta_name::Symbol,
        occasion_names::Vector{Symbol};
        omega_iov::Float64=0.04  # Default 20% CV
    )
        @assert length(occasion_names) >= 2 "Need at least 2 occasions for IOV"
        @assert omega_iov >= 0.0 "omega_iov must be non-negative"
        new(eta_name, occasion_names, omega_iov)
    end
end

"""
Subject occasion data for IOV models.

Maps observations to occasions for IOV computation.

Fields:
- subject_id: Subject identifier
- occasion_indices: Vector mapping each observation to an occasion (1-indexed)
"""
struct OccasionData
    subject_id::String
    occasion_indices::Vector{Int}

    function OccasionData(subject_id::String, occasion_indices::Vector{Int})
        @assert all(occasion_indices .>= 1) "Occasion indices must be >= 1"
        new(subject_id, occasion_indices)
    end
end

# ------------------------------------------------------------------
# BLQ/Censoring Handling
# ------------------------------------------------------------------

"""
Method for handling Below Limit of Quantification (BLQ) observations in estimation.

Following FDA/EMA guidelines and NONMEM conventions:
- M1: Discard all BLQ observations (simplest, may bias)
- M2: Impute BLQ with LLOQ/2 or 0 (legacy, simple)
- M3: Censored likelihood P(Y < LLOQ) (gold standard)

Note: This is distinct from BLQMethod in vpc.jl which is used for VPC analysis.
"""
@enum EstimationBLQMethod begin
    BLQ_M1_DISCARD      # Discard all BLQ observations
    BLQ_M2_IMPUTE_HALF  # Replace BLQ with LLOQ/2
    BLQ_M2_IMPUTE_ZERO  # Replace BLQ with 0
    BLQ_M3_LIKELIHOOD   # Censored likelihood (gold standard)
end

"""
Configuration for BLQ/censoring handling in estimation.

# Fields
- `method`: BLQ handling method (M1, M2, or M3)
- `lloq`: Global lower limit of quantification (can be overridden per subject)
- `max_consecutive_blq`: Maximum consecutive BLQ before warning (default: 5)
- `report_blq_summary`: Include BLQ summary in results (default: true)

# Example
```julia
# Use M3 censored likelihood (recommended)
blq_config = BLQConfig(BLQ_M3_LIKELIHOOD, 0.5)

# Use M1 to discard BLQ observations
blq_config = BLQConfig(BLQ_M1_DISCARD, 0.5)
```
"""
struct BLQConfig
    method::EstimationBLQMethod
    lloq::Float64
    max_consecutive_blq::Int
    report_blq_summary::Bool

    function BLQConfig(
        method::EstimationBLQMethod=BLQ_M3_LIKELIHOOD,
        lloq::Float64=0.0;
        max_consecutive_blq::Int=5,
        report_blq_summary::Bool=true
    )
        @assert lloq >= 0.0 "LLOQ must be non-negative"
        @assert max_consecutive_blq >= 1 "max_consecutive_blq must be positive"
        new(method, lloq, max_consecutive_blq, report_blq_summary)
    end
end

# Convenience constructors
BLQConfig(lloq::Float64) = BLQConfig(BLQ_M3_LIKELIHOOD, lloq)

"""
Summary of BLQ observations in estimation.

# Fields
- `total_observations`: Total number of observations
- `blq_observations`: Number of BLQ observations
- `blq_percentage`: Percentage of BLQ observations
- `blq_by_subject`: Number of BLQ observations per subject
- `method_used`: BLQ method that was used
- `warnings`: Any BLQ-related warnings
"""
struct BLQSummary
    total_observations::Int
    blq_observations::Int
    blq_percentage::Float64
    blq_by_subject::Dict{String,Int}
    method_used::EstimationBLQMethod
    warnings::Vector{String}
end

# ------------------------------------------------------------------
# Estimation Configuration
# ------------------------------------------------------------------

"""
Configuration for parameter estimation.

Fields:
- method: Estimation method (FOCEIMethod, SAEMMethod, or LaplacianMethod)
- theta_init: Initial values for fixed effects (population parameters)
- theta_lower: Lower bounds for theta
- theta_upper: Upper bounds for theta
- theta_names: Names of theta parameters
- omega_init: Initial omega matrix (variance of random effects)
- omega_structure: Structure of omega matrix
- omega_names: Names of eta parameters
- sigma_init: Initial residual error specification
- max_iter: Maximum outer iterations
- tol: Convergence tolerance
- compute_se: Whether to compute standard errors
- compute_ci: Whether to compute confidence intervals
- ci_level: Confidence level (default: 0.95)
- verbose: Print progress during estimation
- seed: Random seed for reproducibility
- covariate_on_iiv: Covariate effects on IIV (optional)
- iov_specs: Inter-occasion variability specifications (optional)
- variance_log_lower: Lower bound for log-scale variance parameters (default: -10.0)
- variance_log_upper: Upper bound for log-scale variance parameters (default: 5.0)
- blq_config: BLQ/censoring handling configuration (optional)
"""
struct EstimationConfig{M<:EstimationMethod}
    method::M
    theta_init::Vector{Float64}
    theta_lower::Vector{Float64}
    theta_upper::Vector{Float64}
    theta_names::Vector{Symbol}
    omega_init::Matrix{Float64}
    omega_structure::OmegaStructure
    omega_names::Vector{Symbol}
    sigma_init::ResidualErrorSpec
    max_iter::Int
    tol::Float64
    compute_se::Bool
    compute_ci::Bool
    ci_level::Float64
    verbose::Bool
    seed::UInt64
    covariate_on_iiv::Vector{CovariateOnIIV}
    iov_specs::Vector{EstimationIOVSpec}
    variance_log_lower::Float64
    variance_log_upper::Float64
    blq_config::Union{Nothing,BLQConfig}

    function EstimationConfig(
        method::M;
        theta_init::Vector{Float64},
        theta_lower::Vector{Float64}=fill(-Inf, length(theta_init)),
        theta_upper::Vector{Float64}=fill(Inf, length(theta_init)),
        theta_names::Vector{Symbol}=Symbol[],
        omega_init::Matrix{Float64},
        omega_structure::OmegaStructure=DiagonalOmega(),
        omega_names::Vector{Symbol}=Symbol[],
        sigma_init::ResidualErrorSpec,
        max_iter::Int=500,
        tol::Float64=1e-4,
        compute_se::Bool=true,
        compute_ci::Bool=false,
        ci_level::Float64=0.95,
        verbose::Bool=false,
        seed::UInt64=UInt64(12345),
        covariate_on_iiv::Vector{CovariateOnIIV}=CovariateOnIIV[],
        iov_specs::Vector{EstimationIOVSpec}=EstimationIOVSpec[],
        variance_log_lower::Float64=-10.0,
        variance_log_upper::Float64=5.0,
        blq_config::Union{Nothing,BLQConfig}=nothing
    ) where {M<:EstimationMethod}
        n_theta = length(theta_init)
        n_eta = size(omega_init, 1)

        # Validate dimensions
        @assert length(theta_lower) == n_theta "theta_lower length mismatch"
        @assert length(theta_upper) == n_theta "theta_upper length mismatch"
        @assert size(omega_init, 1) == size(omega_init, 2) "omega_init must be square"
        @assert issymmetric(omega_init) "omega_init must be symmetric"
        @assert all(theta_lower .<= theta_init .<= theta_upper) "theta_init out of bounds"
        @assert 0.0 < ci_level < 1.0 "ci_level must be in (0, 1)"

        # Default names if not provided
        if isempty(theta_names)
            theta_names = [Symbol("theta_$i") for i in 1:n_theta]
        end
        if isempty(omega_names)
            omega_names = [Symbol("eta_$i") for i in 1:n_eta]
        end

        new{M}(
            method, theta_init, theta_lower, theta_upper, theta_names,
            omega_init, omega_structure, omega_names, sigma_init,
            max_iter, tol, compute_se, compute_ci, ci_level, verbose, seed,
            covariate_on_iiv, iov_specs, variance_log_lower, variance_log_upper,
            blq_config
        )
    end
end

# ------------------------------------------------------------------
# Individual Estimates
# ------------------------------------------------------------------

"""
Estimated quantities for a single individual.

Fields:
- subject_id: Subject identifier
- eta: Estimated random effects (empirical Bayes estimates)
- eta_se: Standard errors of eta (if computed)
- ipred: Individual predictions at observation times
- pred: Population predictions at observation times
- cwres: Conditional weighted residuals
- iwres: Individual weighted residuals
- wres: Weighted residuals (population)
- ofv_contribution: Subject's contribution to objective function
"""
struct IndividualEstimate
    subject_id::String
    eta::Vector{Float64}
    eta_se::Union{Nothing,Vector{Float64}}
    ipred::Vector{Float64}
    pred::Vector{Float64}
    cwres::Vector{Float64}
    iwres::Vector{Float64}
    wres::Vector{Float64}
    ofv_contribution::Float64
end

# ------------------------------------------------------------------
# Estimation Result
# ------------------------------------------------------------------

"""
Result of parameter estimation.

Fields:
- config: EstimationConfig used
- theta: Estimated fixed effects
- theta_se: Standard errors of theta (from Hessian)
- theta_se_robust: Robust standard errors (sandwich estimator)
- theta_rse: Relative standard errors (CV%)
- theta_ci_lower: Lower confidence bounds
- theta_ci_upper: Upper confidence bounds
- omega: Estimated omega matrix
- omega_se: Standard errors of omega elements
- omega_corr: Correlation matrix derived from omega
- sigma: Estimated residual error specification
- sigma_se: Standard errors of sigma parameters
- ofv: Objective function value (-2LL)
- aic: Akaike Information Criterion
- bic: Bayesian Information Criterion
- individuals: Per-subject estimates
- convergence: Did optimization converge?
- n_iterations: Number of iterations used
- gradient_norm: Final gradient norm
- condition_number: Condition number of Hessian
- eigenvalue_ratio: Ratio of smallest to largest eigenvalue
- covariance_step_successful: Did SE computation succeed?
- messages: Any warnings or messages
- runtime_seconds: Total runtime
- blq_summary: Summary of BLQ observations (if BLQ handling enabled)
"""
struct EstimationResult{M<:EstimationMethod}
    config::EstimationConfig{M}

    # Fixed effects
    theta::Vector{Float64}
    theta_se::Union{Nothing,Vector{Float64}}
    theta_se_robust::Union{Nothing,Vector{Float64}}  # Sandwich estimator SEs
    theta_rse::Union{Nothing,Vector{Float64}}
    theta_ci_lower::Union{Nothing,Vector{Float64}}
    theta_ci_upper::Union{Nothing,Vector{Float64}}

    # Random effects variance
    omega::Matrix{Float64}
    omega_se::Union{Nothing,Matrix{Float64}}
    omega_corr::Matrix{Float64}

    # Residual error
    sigma::ResidualErrorSpec
    sigma_se::Union{Nothing,Vector{Float64}}

    # Model fit statistics
    ofv::Float64
    aic::Float64
    bic::Float64

    # Individual estimates
    individuals::Vector{IndividualEstimate}

    # Convergence information
    convergence::Bool
    n_iterations::Int
    gradient_norm::Float64
    condition_number::Float64
    eigenvalue_ratio::Float64
    covariance_step_successful::Bool

    # Messages and runtime
    messages::Vector{String}
    runtime_seconds::Float64

    # BLQ summary
    blq_summary::Union{Nothing,BLQSummary}
end

# ------------------------------------------------------------------
# Utility Functions for Results
# ------------------------------------------------------------------

"""
Get number of estimated parameters (for AIC/BIC calculation).
"""
function n_parameters(result::EstimationResult)::Int
    n_theta = length(result.theta)
    n_omega = _count_omega_params(result.config.omega_structure, size(result.omega, 1))
    n_sigma = _count_sigma_params(result.sigma)
    return n_theta + n_omega + n_sigma
end

function _count_omega_params(::DiagonalOmega, n::Int)::Int
    return n  # Just the diagonal elements
end

function _count_omega_params(structure::BlockOmega, n::Int)::Int
    total = 0
    for block_size in structure.block_sizes
        total += block_size * (block_size + 1) ÷ 2
    end
    return total
end

function _count_omega_params(::FullOmega, n::Int)::Int
    return n * (n + 1) ÷ 2  # Lower triangle including diagonal
end

function _count_sigma_params(spec::ResidualErrorSpec{AdditiveError})::Int
    return 1  # sigma
end

function _count_sigma_params(spec::ResidualErrorSpec{ProportionalError})::Int
    return 1  # sigma
end

function _count_sigma_params(spec::ResidualErrorSpec{CombinedError})::Int
    return 2  # sigma_add, sigma_prop
end

function _count_sigma_params(spec::ResidualErrorSpec{ExponentialError})::Int
    return 1  # sigma
end

"""
Get number of observations across all subjects.
"""
function n_observations(result::EstimationResult)::Int
    return sum(length(ind.ipred) for ind in result.individuals)
end

"""
Get number of subjects.
"""
function n_subjects(result::EstimationResult)::Int
    return length(result.individuals)
end

export n_parameters, n_observations, n_subjects

# ------------------------------------------------------------------
# Pretty Printing
# ------------------------------------------------------------------

function Base.show(io::IO, result::EstimationResult)
    println(io, "EstimationResult")
    println(io, "  Method: ", typeof(result.config.method))
    println(io, "  Convergence: ", result.convergence ? "Yes" : "No")
    println(io, "  OFV: ", round(result.ofv, digits=3))
    println(io, "  AIC: ", round(result.aic, digits=3))
    println(io, "  BIC: ", round(result.bic, digits=3))
    println(io, "  Subjects: ", n_subjects(result))
    println(io, "  Observations: ", n_observations(result))
    println(io, "  Iterations: ", result.n_iterations)
    if result.theta_se !== nothing
        println(io, "\n  Fixed Effects (theta):")
        for (i, (name, val, se)) in enumerate(zip(
            result.config.theta_names, result.theta, result.theta_se
        ))
            rse = result.theta_rse !== nothing ? result.theta_rse[i] : NaN
            println(io, "    $name: $(round(val, digits=4)) (SE=$(round(se, digits=4)), RSE=$(round(rse, digits=1))%)")
        end
    end
end
