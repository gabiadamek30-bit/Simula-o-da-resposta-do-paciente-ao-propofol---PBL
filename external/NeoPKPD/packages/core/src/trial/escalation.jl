# Dose Escalation Simulation
# Implements 3+3, mTPI, and CRM dose escalation algorithms

using Distributions
using StableRNGs

export simulate_escalation, EscalationResult, CohortResult
export simulate_3plus3, simulate_mtpi, simulate_crm, simulate_boin
export boin_boundaries, boin_isotonic_regression

# ============================================================================
# Result Types
# ============================================================================

"""
Result for a single cohort in dose escalation.
"""
struct CohortResult
    cohort_number::Int
    dose_level::Int
    dose_amount::Float64
    n_subjects::Int
    n_dlt::Int
    subject_ids::Vector{Int}
    dlt_flags::Vector{Bool}
    pk_exposures::Vector{Float64}  # Cmax or AUC for each subject
    decision::Symbol  # :escalate, :expand, :deescalate, :stop, :mtd
end

"""
Result of dose escalation simulation.
"""
struct EscalationResult
    design::DoseEscalationDesign
    cohorts::Vector{CohortResult}
    mtd_level::Union{Nothing, Int}
    mtd_dose::Union{Nothing, Float64}
    total_subjects::Int
    total_dlt::Int
    dlt_rate_by_dose::Dict{Int, Float64}
    n_subjects_by_dose::Dict{Int, Int}
    completed::Bool
    termination_reason::String
    seed::UInt64
end

# ============================================================================
# Main Simulation Function
# ============================================================================

"""
    simulate_escalation(design, pk_model_spec, grid, solver; dlt_model, seed, verbose)

Simulate a dose escalation trial.

# Arguments
- `design::DoseEscalationDesign`: Dose escalation design specification
- `pk_model_spec::ModelSpec`: PK model for simulating exposures
- `grid::SimGrid`: Simulation grid
- `solver::SolverSpec`: ODE solver specification
- `dlt_model`: Model for DLT probability (function of dose or exposure)
- `seed::UInt64`: Random seed for reproducibility
- `verbose::Bool`: Print progress

# DLT Models
The dlt_model can be:
- A function `f(dose, exposure) -> probability`
- A tuple `(:threshold, threshold_value)` for exposure-based DLT
- A tuple `(:logistic, α, β)` for logistic dose-response

# Returns
- `EscalationResult` with cohort-by-cohort results and MTD determination
"""
function simulate_escalation(
    design::DoseEscalationDesign,
    pk_model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec;
    dlt_model = (:threshold, 100.0),  # Default: Cmax > 100 causes DLT
    seed::UInt64 = UInt64(12345),
    verbose::Bool = false
)::EscalationResult
    # Dispatch based on escalation rule type
    if design.escalation_rule isa ThreePlusThree
        return simulate_3plus3(design, pk_model_spec, grid, solver;
            dlt_model=dlt_model, seed=seed, verbose=verbose)
    elseif design.escalation_rule isa mTPI
        return simulate_mtpi(design, pk_model_spec, grid, solver;
            dlt_model=dlt_model, seed=seed, verbose=verbose)
    elseif design.escalation_rule isa CRM
        return simulate_crm(design, pk_model_spec, grid, solver;
            dlt_model=dlt_model, seed=seed, verbose=verbose)
    elseif design.escalation_rule isa BOIN
        return simulate_boin(design, pk_model_spec, grid, solver;
            dlt_model=dlt_model, seed=seed, verbose=verbose)
    else
        error("Unsupported escalation rule type: $(typeof(design.escalation_rule))")
    end
end

# ============================================================================
# 3+3 Design
# ============================================================================

"""
    simulate_3plus3(design, pk_model_spec, grid, solver; dlt_model, seed, verbose)

Simulate a 3+3 dose escalation trial.

Traditional 3+3 rules:
- Start with 3 subjects at current dose
- 0/3 DLT: Escalate to next dose
- 1/3 DLT: Expand cohort to 6
- ≥2/3 DLT: De-escalate or declare MTD
- 1/6 DLT: Escalate to next dose (MTD is current dose)
- ≥2/6 DLT: De-escalate (MTD is previous dose)
"""
function simulate_3plus3(
    design::DoseEscalationDesign,
    pk_model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec;
    dlt_model = (:threshold, 100.0),
    seed::UInt64 = UInt64(12345),
    verbose::Bool = false
)::EscalationResult
    rng = StableRNG(seed)

    n_doses = length(design.dose_levels)
    cohort_size = design.cohort_size  # Typically 3

    # Find starting dose level index
    current_level = findfirst(d -> d == design.starting_dose, design.dose_levels)

    # Track results
    cohorts = CohortResult[]
    n_per_dose = Dict{Int, Int}()
    n_dlt_per_dose = Dict{Int, Int}()
    for i in 1:n_doses
        n_per_dose[i] = 0
        n_dlt_per_dose[i] = 0
    end

    subject_counter = 0
    cohort_counter = 0
    mtd_level = nothing
    mtd_dose = nothing
    completed = false
    termination_reason = ""

    while subject_counter < design.max_subjects
        cohort_counter += 1
        current_dose = design.dose_levels[current_level]

        if verbose
            println("Cohort $cohort_counter: Dose level $current_level ($(current_dose) mg)")
        end

        # Determine cohort size for this cohort
        n_at_current = n_per_dose[current_level]
        if n_at_current == 0
            # New dose level, enroll 3
            n_enroll = cohort_size
        elseif n_at_current == cohort_size && n_dlt_per_dose[current_level] == 1
            # Expansion cohort, enroll 3 more
            n_enroll = cohort_size
        else
            # Something went wrong or we shouldn't be here
            n_enroll = cohort_size
        end

        # Cap at max subjects
        if subject_counter + n_enroll > design.max_subjects
            n_enroll = design.max_subjects - subject_counter
        end

        if n_enroll <= 0
            termination_reason = "Maximum subjects reached"
            break
        end

        # Simulate cohort
        subject_ids = Int[]
        dlt_flags = Bool[]
        pk_exposures = Float64[]

        for i in 1:n_enroll
            subject_counter += 1
            push!(subject_ids, subject_counter)

            # Simulate PK for this subject
            exposure = _simulate_subject_exposure(
                pk_model_spec, current_dose, grid, solver, rng
            )
            push!(pk_exposures, exposure)

            # Determine DLT based on model
            dlt = _evaluate_dlt(current_dose, exposure, dlt_model, rng)
            push!(dlt_flags, dlt)
        end

        n_dlt = sum(dlt_flags)
        n_per_dose[current_level] += n_enroll
        n_dlt_per_dose[current_level] += n_dlt

        # Apply 3+3 decision rules
        total_at_dose = n_per_dose[current_level]
        total_dlt_at_dose = n_dlt_per_dose[current_level]

        decision = _3plus3_decision(total_at_dose, total_dlt_at_dose, current_level, n_doses)

        if verbose
            println("  Enrolled: $n_enroll, DLTs: $n_dlt (total: $total_dlt_at_dose/$total_at_dose)")
            println("  Decision: $decision")
        end

        push!(cohorts, CohortResult(
            cohort_counter, current_level, current_dose,
            n_enroll, n_dlt, subject_ids, dlt_flags, pk_exposures, decision
        ))

        # Execute decision
        if decision == :escalate
            if current_level < n_doses
                current_level += 1
            else
                # At highest dose with acceptable toxicity
                mtd_level = current_level
                mtd_dose = design.dose_levels[current_level]
                completed = true
                termination_reason = "MTD at highest dose level"
                break
            end
        elseif decision == :expand
            # Will expand in next iteration
        elseif decision == :deescalate
            if current_level > 1
                # MTD is one level below current
                mtd_level = current_level - 1
                mtd_dose = design.dose_levels[mtd_level]
                completed = true
                termination_reason = "MTD determined (dose level $(mtd_level))"
            else
                # At lowest dose with excessive toxicity
                mtd_level = nothing
                mtd_dose = nothing
                completed = true
                termination_reason = "Excessive toxicity at lowest dose"
            end
            break
        elseif decision == :mtd
            mtd_level = current_level
            mtd_dose = design.dose_levels[current_level]
            completed = true
            termination_reason = "MTD determined (dose level $(mtd_level))"
            break
        elseif decision == :stop
            completed = true
            termination_reason = "Trial stopped"
            break
        end
    end

    if !completed
        termination_reason = "Maximum subjects reached"
        # Determine MTD from available data
        mtd_level, mtd_dose = _determine_mtd_from_data(design, n_per_dose, n_dlt_per_dose)
    end

    # Compute DLT rates
    dlt_rates = Dict{Int, Float64}()
    for i in 1:n_doses
        if n_per_dose[i] > 0
            dlt_rates[i] = n_dlt_per_dose[i] / n_per_dose[i]
        end
    end

    return EscalationResult(
        design,
        cohorts,
        mtd_level,
        mtd_dose,
        subject_counter,
        sum(values(n_dlt_per_dose)),
        dlt_rates,
        n_per_dose,
        completed,
        termination_reason,
        seed
    )
end

"""
Apply 3+3 decision rules.
"""
function _3plus3_decision(n_treated::Int, n_dlt::Int, current_level::Int, n_doses::Int)::Symbol
    if n_treated == 3
        if n_dlt == 0
            return :escalate
        elseif n_dlt == 1
            return :expand
        else  # n_dlt >= 2
            return :deescalate
        end
    elseif n_treated == 6
        if n_dlt <= 1
            return :mtd  # MTD is current dose
        else  # n_dlt >= 2
            return :deescalate
        end
    else
        # Shouldn't happen in standard 3+3
        return :escalate
    end
end

# ============================================================================
# mTPI Design
# ============================================================================

"""
    simulate_mtpi(design, pk_model_spec, grid, solver; dlt_model, seed, verbose)

Simulate a modified Toxicity Probability Interval (mTPI) dose escalation trial.

mTPI uses Bayesian posterior probability to make decisions based on
whether the true DLT rate falls within an equivalence interval around the target.
"""
function simulate_mtpi(
    design::DoseEscalationDesign,
    pk_model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec;
    dlt_model = (:threshold, 100.0),
    seed::UInt64 = UInt64(12345),
    verbose::Bool = false
)::EscalationResult
    rng = StableRNG(seed)

    rule = design.escalation_rule::mTPI
    n_doses = length(design.dose_levels)
    cohort_size = design.cohort_size

    current_level = findfirst(d -> d == design.starting_dose, design.dose_levels)

    cohorts = CohortResult[]
    n_per_dose = Dict{Int, Int}()
    n_dlt_per_dose = Dict{Int, Int}()
    for i in 1:n_doses
        n_per_dose[i] = 0
        n_dlt_per_dose[i] = 0
    end

    subject_counter = 0
    cohort_counter = 0
    mtd_level = nothing
    mtd_dose = nothing
    completed = false
    termination_reason = ""

    # mTPI parameters
    target = rule.target_dlt_rate
    eps1, eps2 = rule.equivalence_interval
    max_cohorts = design.max_subjects ÷ cohort_size

    while cohort_counter < max_cohorts
        cohort_counter += 1
        current_dose = design.dose_levels[current_level]

        if verbose
            println("Cohort $cohort_counter: Dose level $current_level ($(current_dose) mg)")
        end

        n_enroll = min(cohort_size, design.max_subjects - subject_counter)
        if n_enroll <= 0
            break
        end

        # Simulate cohort
        subject_ids = Int[]
        dlt_flags = Bool[]
        pk_exposures = Float64[]

        for i in 1:n_enroll
            subject_counter += 1
            push!(subject_ids, subject_counter)

            exposure = _simulate_subject_exposure(pk_model_spec, current_dose, grid, solver, rng)
            push!(pk_exposures, exposure)

            dlt = _evaluate_dlt(current_dose, exposure, dlt_model, rng)
            push!(dlt_flags, dlt)
        end

        n_dlt = sum(dlt_flags)
        n_per_dose[current_level] += n_enroll
        n_dlt_per_dose[current_level] += n_dlt

        # Compute mTPI decision
        n_total = n_per_dose[current_level]
        n_dlt_total = n_dlt_per_dose[current_level]

        decision = _mtpi_decision(n_total, n_dlt_total, target, eps1, eps2,
                                   current_level, n_doses, rng)

        if verbose
            println("  Enrolled: $n_enroll, DLTs: $n_dlt (total: $n_dlt_total/$n_total)")
            println("  Decision: $decision")
        end

        push!(cohorts, CohortResult(
            cohort_counter, current_level, current_dose,
            n_enroll, n_dlt, subject_ids, dlt_flags, pk_exposures, decision
        ))

        # Execute decision
        if decision == :escalate
            if current_level < n_doses
                current_level += 1
            else
                mtd_level = current_level
                mtd_dose = design.dose_levels[mtd_level]
                completed = true
                termination_reason = "MTD at highest dose"
                break
            end
        elseif decision == :stay
            # Continue at same dose
        elseif decision == :deescalate
            if current_level > 1
                current_level -= 1
            else
                # Can't go lower
                completed = true
                termination_reason = "Excessive toxicity at lowest dose"
                break
            end
        elseif decision == :stop
            # Determine MTD from data
            mtd_level, mtd_dose = _determine_mtd_from_data(design, n_per_dose, n_dlt_per_dose)
            completed = true
            termination_reason = "mTPI stopping rule"
            break
        end
    end

    if !completed
        mtd_level, mtd_dose = _determine_mtd_from_data(design, n_per_dose, n_dlt_per_dose)
        termination_reason = "Maximum cohorts reached"
        completed = true
    end

    dlt_rates = Dict{Int, Float64}()
    for i in 1:n_doses
        if n_per_dose[i] > 0
            dlt_rates[i] = n_dlt_per_dose[i] / n_per_dose[i]
        end
    end

    return EscalationResult(
        design, cohorts, mtd_level, mtd_dose,
        subject_counter, sum(values(n_dlt_per_dose)),
        dlt_rates, n_per_dose, completed, termination_reason, seed
    )
end

"""
mTPI decision based on posterior probability.
"""
function _mtpi_decision(n::Int, y::Int, target::Float64, eps1::Float64, eps2::Float64,
                        current_level::Int, n_doses::Int, rng)::Symbol
    # Use Beta-Binomial conjugate prior (Beta(1,1) = Uniform)
    α_post = 1 + y
    β_post = 1 + n - y

    # Compute posterior probabilities for each interval
    # Under-dosing: [0, eps1)
    # Target: [eps1, eps2]
    # Over-dosing: (eps2, 1]
    p_under = cdf(Beta(α_post, β_post), eps1)
    p_target = cdf(Beta(α_post, β_post), eps2) - p_under
    p_over = 1 - cdf(Beta(α_post, β_post), eps2)

    # mTPI decision rule
    if p_under > p_target && p_under > p_over
        # Most likely under-dosing, escalate
        return :escalate
    elseif p_over > p_target && p_over > p_under
        # Most likely over-dosing, de-escalate
        return :deescalate
    else
        # Most likely at target, stay
        return :stay
    end
end

# ============================================================================
# CRM Design
# ============================================================================

"""
    simulate_crm(design, pk_model_spec, grid, solver; dlt_model, seed, verbose)

Simulate a Continual Reassessment Method (CRM) dose escalation trial.

CRM uses Bayesian updating of a dose-toxicity model to recommend
the dose closest to the target DLT rate.
"""
function simulate_crm(
    design::DoseEscalationDesign,
    pk_model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec;
    dlt_model = (:threshold, 100.0),
    seed::UInt64 = UInt64(12345),
    verbose::Bool = false
)::EscalationResult
    rng = StableRNG(seed)

    rule = design.escalation_rule::CRM
    n_doses = length(design.dose_levels)
    cohort_size = design.cohort_size

    current_level = findfirst(d -> d == design.starting_dose, design.dose_levels)

    # CRM skeleton (prior DLT probabilities)
    skeleton = rule.skeleton
    if length(skeleton) != n_doses
        # Interpolate or use default
        skeleton = _generate_skeleton(n_doses, rule.target_dlt_rate)
    end

    cohorts = CohortResult[]
    n_per_dose = Dict{Int, Int}()
    n_dlt_per_dose = Dict{Int, Int}()
    for i in 1:n_doses
        n_per_dose[i] = 0
        n_dlt_per_dose[i] = 0
    end

    # CRM model parameter (for power model: p_d = skeleton[d]^exp(β))
    β_prior_mean = 0.0
    β_prior_var = 1.0  # Prior variance

    subject_counter = 0
    cohort_counter = 0
    mtd_level = nothing
    mtd_dose = nothing
    completed = false
    termination_reason = ""

    max_cohorts = design.max_subjects ÷ cohort_size

    while cohort_counter < max_cohorts
        cohort_counter += 1
        current_dose = design.dose_levels[current_level]

        if verbose
            println("Cohort $cohort_counter: Dose level $current_level ($(current_dose) mg)")
        end

        n_enroll = min(cohort_size, design.max_subjects - subject_counter)
        if n_enroll <= 0
            break
        end

        # Simulate cohort
        subject_ids = Int[]
        dlt_flags = Bool[]
        pk_exposures = Float64[]

        for i in 1:n_enroll
            subject_counter += 1
            push!(subject_ids, subject_counter)

            exposure = _simulate_subject_exposure(pk_model_spec, current_dose, grid, solver, rng)
            push!(pk_exposures, exposure)

            dlt = _evaluate_dlt(current_dose, exposure, dlt_model, rng)
            push!(dlt_flags, dlt)
        end

        n_dlt = sum(dlt_flags)
        n_per_dose[current_level] += n_enroll
        n_dlt_per_dose[current_level] += n_dlt

        # Update CRM model and get recommendation
        β_post_mean, β_post_var = _update_crm_model(
            n_per_dose, n_dlt_per_dose, skeleton, β_prior_mean, β_prior_var
        )

        recommended_level = _crm_recommendation(
            skeleton, β_post_mean, rule.target_dlt_rate, current_level, n_doses
        )

        # Determine decision
        if recommended_level > current_level
            decision = :escalate
        elseif recommended_level < current_level
            decision = :deescalate
        else
            decision = :stay
        end

        if verbose
            estimated_probs = [skeleton[d]^exp(β_post_mean) for d in 1:n_doses]
            println("  Enrolled: $n_enroll, DLTs: $n_dlt")
            println("  β estimate: $(round(β_post_mean, digits=3))")
            println("  Estimated DLT probs: $(round.(estimated_probs, digits=3))")
            println("  Recommended level: $recommended_level, Decision: $decision")
        end

        push!(cohorts, CohortResult(
            cohort_counter, current_level, current_dose,
            n_enroll, n_dlt, subject_ids, dlt_flags, pk_exposures, decision
        ))

        # Execute decision (with safety constraint: can only escalate one level)
        if decision == :escalate
            if current_level < n_doses
                current_level = min(recommended_level, current_level + 1)
            end
        elseif decision == :deescalate
            current_level = max(recommended_level, 1)
        end
    end

    # Final MTD determination
    if subject_counter > 0
        β_final, _ = _update_crm_model(n_per_dose, n_dlt_per_dose, skeleton, β_prior_mean, β_prior_var)
        mtd_level = _crm_recommendation(skeleton, β_final, rule.target_dlt_rate, 1, n_doses)
        mtd_dose = design.dose_levels[mtd_level]
    end

    completed = true
    termination_reason = "CRM completed"

    dlt_rates = Dict{Int, Float64}()
    for i in 1:n_doses
        if n_per_dose[i] > 0
            dlt_rates[i] = n_dlt_per_dose[i] / n_per_dose[i]
        end
    end

    return EscalationResult(
        design, cohorts, mtd_level, mtd_dose,
        subject_counter, sum(values(n_dlt_per_dose)),
        dlt_rates, n_per_dose, completed, termination_reason, seed
    )
end

"""
Update CRM model parameter using Bayesian updating.
"""
function _update_crm_model(
    n_per_dose::Dict{Int, Int},
    n_dlt_per_dose::Dict{Int, Int},
    skeleton::Vector{Float64},
    prior_mean::Float64,
    prior_var::Float64
)
    # Simplified conjugate update for illustration
    # In practice, use numerical integration or MCMC

    # Collect all data
    doses_treated = [d for d in keys(n_per_dose) if n_per_dose[d] > 0]

    if isempty(doses_treated)
        return prior_mean, prior_var
    end

    # Grid-based approximation for posterior
    β_grid = range(-3, 3, length=101)
    log_posterior = zeros(length(β_grid))

    for (i, β) in enumerate(β_grid)
        # Prior
        log_prior = -0.5 * (β - prior_mean)^2 / prior_var

        # Likelihood
        log_lik = 0.0
        for d in doses_treated
            p_d = skeleton[d]^exp(β)
            n = n_per_dose[d]
            y = n_dlt_per_dose[d]
            # Binomial log-likelihood
            log_lik += y * log(max(p_d, 1e-10)) + (n - y) * log(max(1 - p_d, 1e-10))
        end

        log_posterior[i] = log_prior + log_lik
    end

    # Normalize
    log_posterior .-= maximum(log_posterior)
    posterior = exp.(log_posterior)
    posterior ./= sum(posterior)

    # Compute posterior mean and variance
    post_mean = sum(β_grid .* posterior)
    post_var = sum((β_grid .- post_mean).^2 .* posterior)

    return post_mean, post_var
end

"""
Get CRM dose recommendation.
"""
function _crm_recommendation(
    skeleton::Vector{Float64},
    β::Float64,
    target::Float64,
    current_level::Int,
    n_doses::Int
)::Int
    # Compute estimated DLT probabilities
    estimated_probs = [skeleton[d]^exp(β) for d in 1:n_doses]

    # Find dose with probability closest to target
    distances = abs.(estimated_probs .- target)
    recommended = argmin(distances)

    return recommended
end

"""
Generate default skeleton for CRM.
"""
function _generate_skeleton(n_doses::Int, target::Float64)::Vector{Float64}
    # Generate skeleton such that probabilities increase from near 0 to near target*2
    return [target * 2 * i / n_doses for i in 1:n_doses]
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
Simulate a single subject's PK exposure at a given dose.
"""
function _simulate_subject_exposure(
    pk_model_spec::ModelSpec,
    dose::Float64,
    grid::SimGrid,
    solver::SolverSpec,
    rng
)::Float64
    # Create a dose event at time 0
    doses = [DoseEvent(0.0, dose)]

    # Create a new ModelSpec with the updated dose
    new_spec = ModelSpec(
        pk_model_spec.kind,
        pk_model_spec.name,
        pk_model_spec.params,
        doses
    )

    # Simulate
    result = simulate(new_spec, grid, solver)

    # Return Cmax as exposure metric
    conc = result.observations[:conc]
    return maximum(conc)
end

"""
Evaluate DLT based on dose/exposure model.

dlt_model can be:
- A function taking (exposure) and returning probability
- A function taking (dose, exposure) and returning probability
- A tuple (:threshold, value) for threshold-based DLT
- A tuple (:logistic, α, β) for logistic dose-response
- A tuple (:probability, Dict/function) for dose-specific probabilities
"""
function _evaluate_dlt(dose::Float64, exposure::Float64, dlt_model, rng)::Bool
    if dlt_model isa Tuple && dlt_model[1] == :threshold
        threshold = dlt_model[2]
        return exposure > threshold
    elseif dlt_model isa Tuple && dlt_model[1] == :logistic
        α, β = dlt_model[2], dlt_model[3]
        prob = 1 / (1 + exp(-(α + β * dose)))
        return rand(rng) < prob
    elseif dlt_model isa Tuple && dlt_model[1] == :probability
        prob_per_dose = dlt_model[2]  # Dict or function
        if prob_per_dose isa Dict
            prob = get(prob_per_dose, dose, 0.0)
        else
            prob = prob_per_dose(dose)
        end
        return rand(rng) < prob
    elseif dlt_model isa Function
        # Try to call with (exposure) first, then (dose, exposure)
        prob = try
            dlt_model(exposure)  # Single argument: exposure-based
        catch
            dlt_model(dose, exposure)  # Two arguments: dose and exposure
        end
        return rand(rng) < prob
    else
        # Default: no DLT
        return false
    end
end

"""
Determine MTD from accumulated data.
"""
function _determine_mtd_from_data(
    design::DoseEscalationDesign,
    n_per_dose::Dict{Int, Int},
    n_dlt_per_dose::Dict{Int, Int}
)::Tuple{Union{Nothing, Int}, Union{Nothing, Float64}}
    n_doses = length(design.dose_levels)

    # Find highest dose with DLT rate ≤ 33%
    mtd_level = nothing
    for d in n_doses:-1:1
        if n_per_dose[d] >= 3  # Need at least 3 subjects
            rate = n_dlt_per_dose[d] / n_per_dose[d]
            if rate <= 1/3
                mtd_level = d
                break
            end
        end
    end

    if mtd_level !== nothing
        return mtd_level, design.dose_levels[mtd_level]
    else
        return nothing, nothing
    end
end

# ============================================================================
# BOIN Design (Bayesian Optimal Interval)
# ============================================================================

"""
    boin_boundaries(target::Float64; p1_ratio::Float64 = 0.6, p2_ratio::Float64 = 1.4)

Compute BOIN escalation and de-escalation boundaries.

# Arguments
- `target`: Target DLT rate
- `p1_ratio`: Ratio of p1 to target (default 0.6)
- `p2_ratio`: Ratio of p2 to target (default 1.4)

# Returns
- `(λe, λd)`: Escalation and de-escalation boundaries
"""
function boin_boundaries(target::Float64; p1_ratio::Float64 = 0.6, p2_ratio::Float64 = 1.4)
    p1 = p1_ratio * target
    p2 = p2_ratio * target

    # Optimal boundaries from Bayesian decision theory
    λe = log((1 - p1) / (1 - target)) / log(target * (1 - p1) / (p1 * (1 - target)))
    λd = log((1 - target) / (1 - p2)) / log(p2 * (1 - target) / (target * (1 - p2)))

    return λe, λd
end

"""
    boin_decision(n_dlt::Int, n::Int, boin::BOIN)

Make BOIN dose escalation decision.

# Arguments
- `n_dlt`: Number of DLTs at current dose
- `n`: Number of subjects at current dose
- `boin`: BOIN design parameters

# Returns
- `:escalate`, `:stay`, `:deescalate`, or `:eliminate`
"""
function boin_decision(n_dlt::Int, n::Int, boin::BOIN)
    p_hat = n_dlt / n

    if p_hat <= boin.cutoff_e
        return :escalate
    elseif p_hat >= boin.cutoff_d
        # Check for elimination
        if n >= boin.n_eliminate
            # Use posterior probability for elimination
            # P(p > target | data) using Beta-Binomial
            posterior = Beta(1 + n_dlt, 1 + n - n_dlt)
            prob_toxic = 1 - cdf(posterior, boin.target_dlt_rate)
            if prob_toxic >= boin.eliminate_threshold
                return :eliminate
            end
        end
        return :deescalate
    else
        return :stay
    end
end

"""
    boin_isotonic_regression(n_per_dose::Dict{Int, Int}, n_dlt_per_dose::Dict{Int, Int})

Apply isotonic regression to DLT rates for MTD selection.

BOIN uses isotonic regression to ensure DLT rates are monotonically increasing
with dose when selecting the MTD at the end of the trial.

# Returns
- Vector of isotonic DLT rates
"""
function boin_isotonic_regression(n_per_dose::Dict{Int, Int}, n_dlt_per_dose::Dict{Int, Int})
    n_doses = maximum(keys(n_per_dose))

    # Extract raw rates for doses with data
    rates = Float64[]
    weights = Float64[]

    for d in 1:n_doses
        if n_per_dose[d] > 0
            push!(rates, n_dlt_per_dose[d] / n_per_dose[d])
            push!(weights, Float64(n_per_dose[d]))
        else
            push!(rates, NaN)
            push!(weights, 0.0)
        end
    end

    # Pool-adjacent-violators algorithm (PAVA) for isotonic regression
    isotonic_rates = copy(rates)
    n = length(isotonic_rates)

    # Forward pass - ensure non-decreasing
    for i in 2:n
        if isnan(isotonic_rates[i])
            continue
        end
        for j in (i-1):-1:1
            if isnan(isotonic_rates[j])
                continue
            end
            if isotonic_rates[j] > isotonic_rates[i]
                # Pool j and i
                total_dlt = 0.0
                total_n = 0.0
                for k in j:i
                    if !isnan(rates[k]) && weights[k] > 0
                        total_dlt += rates[k] * weights[k]
                        total_n += weights[k]
                    end
                end
                pooled_rate = total_n > 0 ? total_dlt / total_n : 0.0
                for k in j:i
                    if !isnan(isotonic_rates[k])
                        isotonic_rates[k] = pooled_rate
                    end
                end
            end
            break
        end
    end

    return isotonic_rates
end

"""
    simulate_boin(design, pk_model_spec, grid, solver; dlt_model, seed, verbose)

Simulate a BOIN dose escalation trial.

BOIN (Bayesian Optimal Interval) design uses:
1. Pre-calculated optimal boundaries for dose decisions
2. Posterior probability for dose elimination
3. Isotonic regression for MTD selection

# Decision Rules
- If observed DLT rate ≤ λe: Escalate
- If observed DLT rate ≥ λd: De-escalate (or eliminate if P(p > target) ≥ threshold)
- Otherwise: Stay at current dose

# Reference
Liu S, Yuan Y (2015). Bayesian Optimal Interval Designs for Phase I Clinical Trials.
"""
function simulate_boin(
    design::DoseEscalationDesign,
    pk_model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec;
    dlt_model = (:threshold, 100.0),
    seed::UInt64 = UInt64(12345),
    verbose::Bool = false
)::EscalationResult
    rng = StableRNG(seed)
    boin = design.escalation_rule::BOIN

    n_doses = length(design.dose_levels)
    cohort_size = design.cohort_size

    # Find starting dose level index
    current_level = findfirst(d -> d == design.starting_dose, design.dose_levels)

    # Track results
    cohorts = CohortResult[]
    n_per_dose = Dict{Int, Int}()
    n_dlt_per_dose = Dict{Int, Int}()
    eliminated = Set{Int}()

    for i in 1:n_doses
        n_per_dose[i] = 0
        n_dlt_per_dose[i] = 0
    end

    subject_counter = 0
    cohort_counter = 0
    mtd_level = nothing
    mtd_dose = nothing
    completed = false
    termination_reason = ""

    if verbose
        println("BOIN Dose Escalation Simulation")
        println("Target DLT rate: $(boin.target_dlt_rate)")
        println("Escalation cutoff (λe): $(round(boin.cutoff_e, digits=3))")
        println("De-escalation cutoff (λd): $(round(boin.cutoff_d, digits=3))")
        println()
    end

    while subject_counter < design.max_subjects
        # Check if current dose is eliminated
        while current_level in eliminated && current_level > 1
            current_level -= 1
        end

        if current_level in eliminated
            completed = true
            termination_reason = "All doses eliminated"
            break
        end

        cohort_counter += 1
        current_dose = design.dose_levels[current_level]

        if verbose
            println("Cohort $cohort_counter: Dose level $current_level ($(current_dose) mg)")
        end

        # Enroll cohort
        n_enroll = min(cohort_size, design.max_subjects - subject_counter)

        # Simulate subjects
        subject_ids = Int[]
        dlt_flags = Bool[]
        pk_exposures = Float64[]

        for _ in 1:n_enroll
            subject_counter += 1
            push!(subject_ids, subject_counter)

            # Simulate PK
            exposure = _simulate_subject_exposure(pk_model_spec, current_dose, grid, solver, rng)
            push!(pk_exposures, exposure)

            # Determine DLT
            has_dlt = _evaluate_dlt(current_dose, exposure, dlt_model, rng)
            push!(dlt_flags, has_dlt)
        end

        n_dlt = sum(dlt_flags)
        n_per_dose[current_level] += n_enroll
        n_dlt_per_dose[current_level] += n_dlt

        if verbose
            println("  Enrolled: $n_enroll, DLTs: $n_dlt")
            println("  Cumulative at dose: $(n_per_dose[current_level]) subjects, $(n_dlt_per_dose[current_level]) DLTs")
        end

        # Make BOIN decision based on cumulative data at current dose
        decision = boin_decision(n_dlt_per_dose[current_level], n_per_dose[current_level], boin)

        if verbose
            println("  Decision: $decision")
        end

        # Record cohort result
        push!(cohorts, CohortResult(
            cohort_counter,
            current_level,
            current_dose,
            n_enroll,
            n_dlt,
            subject_ids,
            dlt_flags,
            pk_exposures,
            decision
        ))

        # Apply decision
        if decision == :escalate
            if current_level < n_doses && !(current_level + 1 in eliminated)
                current_level += 1
            end
        elseif decision == :deescalate
            if current_level > 1
                current_level -= 1
            else
                # Can't go lower, stay at dose 1
            end
        elseif decision == :eliminate
            push!(eliminated, current_level)
            # Also eliminate all higher doses
            for d in (current_level + 1):n_doses
                push!(eliminated, d)
            end

            if current_level > 1
                current_level -= 1
            else
                completed = true
                termination_reason = "Lowest dose eliminated"
                break
            end
        end
        # :stay - keep current_level unchanged
    end

    if !completed
        termination_reason = "Maximum subjects reached"
        completed = true
    end

    # Determine MTD using isotonic regression
    isotonic_rates = boin_isotonic_regression(n_per_dose, n_dlt_per_dose)

    # Find MTD: highest dose where isotonic rate ≤ target and not eliminated
    for d in n_doses:-1:1
        if n_per_dose[d] > 0 && !(d in eliminated)
            if !isnan(isotonic_rates[d]) && isotonic_rates[d] <= boin.target_dlt_rate
                mtd_level = d
                mtd_dose = design.dose_levels[d]
                break
            end
        end
    end

    if verbose
        println()
        println("Trial completed: $termination_reason")
        if mtd_level !== nothing
            println("MTD: Level $mtd_level ($(mtd_dose) mg)")
        else
            println("MTD: Not determined")
        end
    end

    # Compute DLT rates
    dlt_rates = Dict{Int, Float64}()
    for d in 1:n_doses
        if n_per_dose[d] > 0
            dlt_rates[d] = n_dlt_per_dose[d] / n_per_dose[d]
        else
            dlt_rates[d] = NaN
        end
    end

    return EscalationResult(
        design,
        cohorts,
        mtd_level,
        mtd_dose,
        subject_counter,
        sum(values(n_dlt_per_dose)),
        dlt_rates,
        n_per_dose,
        completed,
        termination_reason,
        seed
    )
end
