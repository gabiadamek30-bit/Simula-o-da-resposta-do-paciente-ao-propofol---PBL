# Stepwise Covariate Modeling (SCM)
# Automatic covariate selection using forward selection and backward elimination

using Distributions
using LinearAlgebra

export SCMResult, SCMSpec, CovariateRelationship
export run_scm, run_scm_native, forward_selection, backward_elimination
# Note: likelihood_ratio_test is defined in diagnostics.jl

# ============================================================================
# Types and Specifications
# ============================================================================

"""
    CovariateRelationship

Defines a potential covariate-parameter relationship to test.

Fields:
- param: Parameter affected by covariate (e.g., :CL)
- covariate: Covariate name (e.g., :WT, :AGE)
- relationship: Type of relationship (:linear, :power, :exponential)
- reference: Reference value for centering (e.g., median weight)
"""
struct CovariateRelationship
    param::Symbol
    covariate::Symbol
    relationship::Symbol  # :linear, :power, :exponential
    reference::Float64

    function CovariateRelationship(
        param::Symbol,
        covariate::Symbol;
        relationship::Symbol=:power,
        reference::Float64=1.0
    )
        relationship in (:linear, :power, :exponential) ||
            error("relationship must be :linear, :power, or :exponential")
        new(param, covariate, relationship, reference)
    end
end

"""
    SCMSpec

Specification for stepwise covariate modeling.

Fields:
- relationships: Vector of covariate relationships to test
- forward_p_value: P-value threshold for forward inclusion (default: 0.05)
- backward_p_value: P-value threshold for backward elimination (default: 0.01)
- max_iterations: Maximum iterations for forward/backward steps
- parallel: Run covariate tests in parallel
"""
struct SCMSpec
    relationships::Vector{CovariateRelationship}
    forward_p_value::Float64
    backward_p_value::Float64
    max_iterations::Int
    parallel::Bool

    function SCMSpec(
        relationships::Vector{CovariateRelationship};
        forward_p_value::Float64=0.05,
        backward_p_value::Float64=0.01,
        max_iterations::Int=100,
        parallel::Bool=false
    )
        0 < forward_p_value < 1 || error("forward_p_value must be in (0, 1)")
        0 < backward_p_value < 1 || error("backward_p_value must be in (0, 1)")
        forward_p_value >= backward_p_value ||
            @warn "backward_p_value should typically be <= forward_p_value"

        new(relationships, forward_p_value, backward_p_value, max_iterations, parallel)
    end
end

"""
    SCMStepResult

Result of a single SCM step (adding or removing one covariate).
"""
struct SCMStepResult
    relationship::CovariateRelationship
    ofv_change::Float64
    p_value::Float64
    coefficient::Float64
    coefficient_se::Float64
    included::Bool
    step_type::Symbol  # :forward or :backward
end

"""
    SCMResult

Complete result of stepwise covariate modeling.

Fields:
- final_model: Vector of included covariate relationships
- forward_steps: Results from forward selection
- backward_steps: Results from backward elimination
- base_ofv: OFV of base model (no covariates)
- final_ofv: OFV of final model
- all_tested: All relationships that were tested
- convergence: Whether SCM converged
"""
struct SCMResult
    final_model::Vector{CovariateRelationship}
    forward_steps::Vector{SCMStepResult}
    backward_steps::Vector{SCMStepResult}
    base_ofv::Float64
    final_ofv::Float64
    all_tested::Vector{CovariateRelationship}
    convergence::Bool
    n_forward_iterations::Int
    n_backward_iterations::Int
end

export SCMStepResult

# ============================================================================
# Covariate Model Building
# Note: likelihood_ratio_test is defined in diagnostics.jl
# ============================================================================

"""
Build covariate effect function based on relationship type.

Returns a function: (theta_pop, cov_value, coefficient) -> theta_individual
"""
function build_covariate_effect(rel::CovariateRelationship)::Function
    if rel.relationship == :linear
        # θ_i = θ_pop * (1 + β * (COV - REF))
        return (theta_pop, cov, beta) -> theta_pop * (1 + beta * (cov - rel.reference))
    elseif rel.relationship == :power
        # θ_i = θ_pop * (COV / REF)^β
        return (theta_pop, cov, beta) -> theta_pop * (cov / rel.reference)^beta
    elseif rel.relationship == :exponential
        # θ_i = θ_pop * exp(β * (COV - REF))
        return (theta_pop, cov, beta) -> theta_pop * exp(beta * (cov - rel.reference))
    else
        error("Unknown relationship type: $(rel.relationship)")
    end
end

"""
Apply covariate effect to a population parameter value.

Arguments:
- theta_pop: Population parameter value
- cov_value: Individual covariate value
- beta: Covariate coefficient
- rel: Covariate relationship specification

Returns:
- Adjusted individual parameter value
"""
function apply_covariate_effect(
    theta_pop::Float64,
    cov_value::Float64,
    beta::Float64,
    rel::CovariateRelationship
)::Float64
    if rel.relationship == :linear
        # θ_i = θ_pop * (1 + β * (COV - REF))
        return theta_pop * (1.0 + beta * (cov_value - rel.reference))
    elseif rel.relationship == :power
        # θ_i = θ_pop * (COV / REF)^β
        # Guard against negative/zero covariate values
        ratio = max(cov_value, 1e-10) / max(rel.reference, 1e-10)
        return theta_pop * ratio^beta
    elseif rel.relationship == :exponential
        # θ_i = θ_pop * exp(β * (COV - REF))
        return theta_pop * exp(beta * (cov_value - rel.reference))
    else
        return theta_pop  # No effect for unknown relationship
    end
end

"""
Get covariate values for a specific covariate from all subjects.

Returns a vector of covariate values, one per subject.
Missing values are replaced with the reference value.
"""
function get_covariate_values(
    observed_data::ObservedData,
    covariate_name::Symbol,
    reference::Float64
)::Vector{Float64}
    values = Float64[]
    for subj in observed_data.subjects
        if haskey(subj.covariates, covariate_name)
            val = subj.covariates[covariate_name]
            push!(values, Float64(val))
        else
            # Use reference value for missing covariates
            push!(values, reference)
        end
    end
    return values
end

"""
Map parameter name to theta index based on model kind.

Returns the 1-based index into theta for the parameter, or nothing if not found.
"""
function get_param_index(model_kind, param::Symbol)::Union{Int, Nothing}
    # Standard PK parameter mappings
    param_maps = Dict(
        :OneCompIVBolus => Dict(:CL => 1, :V => 2),
        :OneCompOralFirstOrder => Dict(:CL => 1, :V => 2, :Ka => 3),
        :TwoCompIVBolus => Dict(:CL => 1, :V1 => 2, :Q => 3, :V2 => 4),
        :TwoCompOralFirstOrder => Dict(:CL => 1, :V1 => 2, :Q => 3, :V2 => 4, :Ka => 5),
        :ThreeCompIVBolus => Dict(:CL => 1, :V1 => 2, :Q2 => 3, :V2 => 4, :Q3 => 5, :V3 => 6),
    )

    model_name = string(typeof(model_kind).name.name)
    if haskey(param_maps, Symbol(model_name))
        param_map = param_maps[Symbol(model_name)]
        if haskey(param_map, param)
            return param_map[param]
        end
        # Also check common aliases
        if param == :V && haskey(param_map, :V1)
            return param_map[:V1]
        end
    end
    return nothing
end

# ============================================================================
# Forward Selection
# ============================================================================

"""
    forward_selection(estimate_fn, base_theta, candidates, spec) -> (selected, steps)

Perform forward covariate selection.

Arguments:
- estimate_fn: Function (relationships) -> (ofv, theta, se) for fitting model
- base_theta: Initial parameter estimates
- candidates: Vector of candidate covariate relationships
- spec: SCM specification

Returns:
- selected: Vector of selected relationships
- steps: Vector of SCMStepResult for each step
"""
function forward_selection(
    estimate_fn::Function,
    base_theta::Vector{Float64},
    candidates::Vector{CovariateRelationship},
    spec::SCMSpec
)::Tuple{Vector{CovariateRelationship}, Vector{SCMStepResult}}
    selected = CovariateRelationship[]
    remaining = copy(candidates)
    steps = SCMStepResult[]

    # Get base model OFV
    base_ofv, _, _ = estimate_fn(selected)

    iteration = 0
    while !isempty(remaining) && iteration < spec.max_iterations
        iteration += 1

        # Test each remaining covariate
        best_candidate = nothing
        best_p_value = 1.0
        best_ofv = base_ofv
        best_coef = 0.0
        best_se = 0.0

        for rel in remaining
            # Fit model with this covariate added
            test_model = vcat(selected, [rel])
            try
                new_ofv, theta, se = estimate_fn(test_model)

                # LRT with 1 df (one additional parameter)
                _, p_value = likelihood_ratio_test(new_ofv, base_ofv, 1)

                if p_value < best_p_value
                    best_candidate = rel
                    best_p_value = p_value
                    best_ofv = new_ofv
                    best_coef = length(theta) > 0 ? theta[end] : 0.0
                    best_se = length(se) > 0 ? se[end] : 0.0
                end
            catch e
                @warn "Failed to fit model with $(rel.param)-$(rel.covariate): $e"
                continue
            end
        end

        # Check if best candidate meets threshold
        if best_candidate !== nothing && best_p_value < spec.forward_p_value
            push!(selected, best_candidate)
            filter!(r -> r != best_candidate, remaining)
            base_ofv = best_ofv

            push!(steps, SCMStepResult(
                best_candidate,
                base_ofv - best_ofv,
                best_p_value,
                best_coef,
                best_se,
                true,
                :forward
            ))
        else
            # No more significant covariates
            break
        end
    end

    return selected, steps
end

# ============================================================================
# Backward Elimination
# ============================================================================

"""
    backward_elimination(estimate_fn, included, spec) -> (final, steps)

Perform backward covariate elimination.

Arguments:
- estimate_fn: Function (relationships) -> (ofv, theta, se) for fitting model
- included: Vector of currently included relationships
- spec: SCM specification

Returns:
- final: Vector of relationships remaining after elimination
- steps: Vector of SCMStepResult for each elimination step
"""
function backward_elimination(
    estimate_fn::Function,
    included::Vector{CovariateRelationship},
    spec::SCMSpec
)::Tuple{Vector{CovariateRelationship}, Vector{SCMStepResult}}
    current = copy(included)
    steps = SCMStepResult[]

    if isempty(current)
        return current, steps
    end

    # Get full model OFV
    full_ofv, _, _ = estimate_fn(current)

    iteration = 0
    while length(current) > 0 && iteration < spec.max_iterations
        iteration += 1

        # Test removing each covariate
        worst_candidate = nothing
        worst_p_value = 0.0
        worst_ofv = full_ofv
        worst_coef = 0.0
        worst_se = 0.0

        for (i, rel) in enumerate(current)
            # Fit model without this covariate
            test_model = vcat(current[1:i-1], current[i+1:end])

            try
                reduced_ofv, theta, se = estimate_fn(test_model)

                # LRT: is the reduced model significantly worse?
                _, p_value = likelihood_ratio_test(full_ofv, reduced_ofv, 1)

                # We want to remove covariates with HIGH p-values (not significant)
                if p_value > worst_p_value
                    worst_candidate = rel
                    worst_p_value = p_value
                    worst_ofv = reduced_ofv
                    # Get coefficient from full model (before removal)
                    worst_coef = i <= length(theta) ? theta[i] : 0.0
                    worst_se = i <= length(se) ? se[i] : 0.0
                end
            catch e
                @warn "Failed to fit reduced model without $(rel.param)-$(rel.covariate): $e"
                continue
            end
        end

        # Check if worst candidate should be removed
        if worst_candidate !== nothing && worst_p_value > spec.backward_p_value
            filter!(r -> r != worst_candidate, current)
            full_ofv = worst_ofv

            push!(steps, SCMStepResult(
                worst_candidate,
                worst_ofv - full_ofv,
                worst_p_value,
                worst_coef,
                worst_se,
                false,
                :backward
            ))
        else
            # All remaining covariates are significant
            break
        end
    end

    return current, steps
end

# ============================================================================
# Main SCM Function
# ============================================================================

"""
    run_scm(estimate_fn, base_theta, spec) -> SCMResult

Run full stepwise covariate modeling.

This performs:
1. Forward selection to add significant covariates
2. Backward elimination to remove non-significant ones

Arguments:
- estimate_fn: Function (relationships) -> (ofv, theta, se)
  Takes a vector of covariate relationships and returns:
  - ofv: Objective function value
  - theta: Parameter estimates
  - se: Standard errors
- base_theta: Initial parameter estimates for base model
- spec: SCM specification

Returns:
- SCMResult with final model and all tested relationships

Example:
```julia
# Define candidate relationships
candidates = [
    CovariateRelationship(:CL, :WT, relationship=:power, reference=70.0),
    CovariateRelationship(:CL, :AGE, relationship=:linear, reference=40.0),
    CovariateRelationship(:V, :WT, relationship=:power, reference=70.0),
    CovariateRelationship(:CL, :CRCL, relationship=:power, reference=100.0)
]

spec = SCMSpec(candidates, forward_p_value=0.05, backward_p_value=0.01)

# estimate_fn should fit the model and return (ofv, theta, se)
result = run_scm(estimate_fn, base_theta, spec)
```
"""
function run_scm(
    estimate_fn::Function,
    base_theta::Vector{Float64},
    spec::SCMSpec
)::SCMResult
    # Get base model OFV
    base_ofv, _, _ = estimate_fn(CovariateRelationship[])

    # Forward selection
    selected, forward_steps = forward_selection(
        estimate_fn, base_theta, spec.relationships, spec
    )

    # Backward elimination
    final_model, backward_steps = backward_elimination(
        estimate_fn, selected, spec
    )

    # Get final OFV
    final_ofv = if isempty(final_model)
        base_ofv
    else
        ofv, _, _ = estimate_fn(final_model)
        ofv
    end

    return SCMResult(
        final_model,
        forward_steps,
        backward_steps,
        base_ofv,
        final_ofv,
        spec.relationships,
        true,  # convergence
        length(forward_steps),
        length(backward_steps)
    )
end

export run_scm, forward_selection, backward_elimination

# ============================================================================
# Utility Functions
# ============================================================================

"""
Print SCM result summary.
"""
function summarize_scm(result::SCMResult)::String
    lines = String[]

    push!(lines, "=== Stepwise Covariate Modeling Results ===")
    push!(lines, "")
    push!(lines, "Base OFV: $(round(result.base_ofv, digits=3))")
    push!(lines, "Final OFV: $(round(result.final_ofv, digits=3))")
    push!(lines, "ΔOFV: $(round(result.base_ofv - result.final_ofv, digits=3))")
    push!(lines, "")

    push!(lines, "Forward Selection Steps: $(result.n_forward_iterations)")
    for step in result.forward_steps
        rel = step.relationship
        push!(lines, "  + $(rel.param) ~ $(rel.covariate) ($(rel.relationship)): p=$(round(step.p_value, digits=4)), β=$(round(step.coefficient, digits=4))")
    end

    push!(lines, "")
    push!(lines, "Backward Elimination Steps: $(result.n_backward_iterations)")
    for step in result.backward_steps
        rel = step.relationship
        push!(lines, "  - $(rel.param) ~ $(rel.covariate) ($(rel.relationship)): p=$(round(step.p_value, digits=4))")
    end

    push!(lines, "")
    push!(lines, "Final Model Covariates:")
    if isempty(result.final_model)
        push!(lines, "  (none)")
    else
        for rel in result.final_model
            push!(lines, "  $(rel.param) ~ $(rel.covariate) ($(rel.relationship), ref=$(rel.reference))")
        end
    end

    return join(lines, "\n")
end

export summarize_scm

"""
Create standard covariate relationships for common PK parameters.

Returns relationships for:
- Weight effects on CL, V (allometric scaling)
- Age effects on CL
- Creatinine clearance effects on CL
- Sex effects on CL, V
"""
function standard_pk_covariates(;
    weight_ref::Float64=70.0,
    age_ref::Float64=40.0,
    crcl_ref::Float64=100.0
)::Vector{CovariateRelationship}
    return [
        # Weight effects (allometric)
        CovariateRelationship(:CL, :WT, relationship=:power, reference=weight_ref),
        CovariateRelationship(:V, :WT, relationship=:power, reference=weight_ref),
        CovariateRelationship(:V1, :WT, relationship=:power, reference=weight_ref),
        CovariateRelationship(:V2, :WT, relationship=:power, reference=weight_ref),
        CovariateRelationship(:Q, :WT, relationship=:power, reference=weight_ref),

        # Age effects
        CovariateRelationship(:CL, :AGE, relationship=:linear, reference=age_ref),

        # Renal function
        CovariateRelationship(:CL, :CRCL, relationship=:power, reference=crcl_ref),
        CovariateRelationship(:CL, :eGFR, relationship=:power, reference=crcl_ref),

        # Sex effects (typically as proportional)
        CovariateRelationship(:CL, :SEX, relationship=:linear, reference=0.0),
        CovariateRelationship(:V, :SEX, relationship=:linear, reference=0.0)
    ]
end

export standard_pk_covariates


# ============================================================================
# Native SCM Function (For Python Integration)
# ============================================================================

"""
    run_scm_native(observed_data, model_spec, est_config, scm_spec; grid, solver, verbose)

Run stepwise covariate modeling without requiring a callback function.

This is the preferred method for calling from Python via JuliaCall, as it
creates the estimation function internally and handles all covariate effect
computations.

# Arguments
- `observed_data`: ObservedData containing subjects with covariates
- `model_spec`: ModelSpec defining the base PK/PD model
- `est_config`: EstimationConfig with estimation settings
- `scm_spec`: SCMSpec with covariate relationships to test
- `grid`: SimGrid for simulation time points
- `solver`: SolverSpec for ODE solver settings
- `verbose`: Print progress information

# Returns
- SCMResult with selected covariates, steps, and OFV changes

# Example
```julia
# Define covariate relationships to test
relationships = [
    CovariateRelationship(:CL, :WT, relationship=:power, reference=70.0),
    CovariateRelationship(:V, :WT, relationship=:power, reference=70.0),
    CovariateRelationship(:CL, :AGE, relationship=:linear, reference=40.0),
]

scm_spec = SCMSpec(relationships, forward_p_value=0.05, backward_p_value=0.01)

result = run_scm_native(
    observed_data,
    model_spec,
    est_config,
    scm_spec,
    grid=sim_grid,
    solver=solver_spec
)
```
"""
function run_scm_native(
    observed_data::ObservedData,
    model_spec::ModelSpec,
    est_config::EstimationConfig,
    scm_spec::SCMSpec;
    grid::SimGrid,
    solver::SolverSpec,
    verbose::Bool=false
)::SCMResult
    # Extract base theta from estimation config
    base_theta = Vector{Float64}(est_config.theta_init)
    n_base_params = length(base_theta)

    # Pre-extract covariate values for each relationship
    covariate_cache = Dict{Symbol, Vector{Float64}}()
    for rel in scm_spec.relationships
        if !haskey(covariate_cache, rel.covariate)
            covariate_cache[rel.covariate] = get_covariate_values(
                observed_data, rel.covariate, rel.reference
            )
        end
    end

    # Create the estimation function that handles covariate effects
    function estimate_with_covariates(relationships::Vector{CovariateRelationship})
        n_covariates = length(relationships)

        if n_covariates == 0
            # Base model: just run standard estimation
            try
                result = estimate(observed_data, model_spec, est_config; grid=grid, solver=solver)
                theta_se = result.theta_se !== nothing ? result.theta_se : zeros(n_base_params)
                return (result.ofv, Vector{Float64}(result.theta), Vector{Float64}(theta_se))
            catch e
                if verbose
                    println("Base estimation failed: $e")
                end
                return (Inf, base_theta, zeros(n_base_params))
            end
        end

        # Extended model with covariate effects
        # theta_extended = [base_params..., beta_1, beta_2, ..., beta_n]
        n_extended = n_base_params + n_covariates
        theta_init_extended = vcat(base_theta, zeros(n_covariates))

        # Build mapping from relationship to parameter index and covariate values
        rel_info = []
        for (i, rel) in enumerate(relationships)
            param_idx = get_param_index(model_spec.kind, rel.param)
            if param_idx === nothing
                if verbose
                    println("Warning: Parameter $(rel.param) not found in model, skipping")
                end
                continue
            end
            cov_values = covariate_cache[rel.covariate]
            push!(rel_info, (rel=rel, param_idx=param_idx, cov_values=cov_values, beta_idx=n_base_params + i))
        end

        if isempty(rel_info)
            # No valid relationships, return base model results
            try
                result = estimate(observed_data, model_spec, est_config; grid=grid, solver=solver)
                theta_se = result.theta_se !== nothing ? result.theta_se : zeros(n_base_params)
                return (result.ofv, Vector{Float64}(result.theta), Vector{Float64}(theta_se))
            catch e
                return (Inf, base_theta, zeros(n_base_params))
            end
        end

        # Compute OFV for extended model by modifying individual predictions
        # For each subject, apply covariate effects to the relevant parameters
        function compute_extended_ofv(theta_ext::Vector{Float64})::Float64
            theta_base = theta_ext[1:n_base_params]
            betas = theta_ext[n_base_params+1:end]

            total_ll = 0.0
            n_obs = 0

            for (subj_idx, subj) in enumerate(observed_data.subjects)
                # Start with population parameters
                theta_individual = copy(theta_base)

                # Apply covariate effects for this subject
                for (j, info) in enumerate(rel_info)
                    beta = betas[j]
                    cov_val = info.cov_values[subj_idx]
                    theta_pop = theta_individual[info.param_idx]
                    theta_individual[info.param_idx] = apply_covariate_effect(
                        theta_pop, cov_val, beta, info.rel
                    )
                end

                # Ensure all parameters are positive
                if any(theta_individual .<= 0)
                    return Inf
                end

                # Create individual model spec with adjusted parameters
                try
                    ind_params = theta_to_params(theta_individual, model_spec)
                    ind_model = ModelSpec(model_spec.kind, model_spec.name, ind_params, subj.doses)

                    # Simulate
                    sim_result = simulate(ind_model, grid, solver)

                    # Compute residuals
                    for (i, t) in enumerate(subj.times)
                        idx = searchsortedfirst(sim_result.t, t)
                        if idx > length(sim_result.t)
                            idx = length(sim_result.t)
                        elseif idx > 1 && abs(sim_result.t[idx] - t) > abs(sim_result.t[idx-1] - t)
                            idx = idx - 1
                        end

                        pred = sim_result.observations[:conc][idx]
                        obs = subj.observations[i]

                        # Proportional error model: SD = sigma * pred
                        sigma = est_config.sigma_init.params.sigma
                        sd = sigma * max(pred, 1e-10)

                        # Log-likelihood contribution
                        resid = (obs - pred) / sd
                        total_ll -= 0.5 * (resid^2 + log(2π) + 2*log(sd))
                        n_obs += 1
                    end
                catch e
                    return Inf
                end
            end

            # Return -2LL (OFV)
            return -2.0 * total_ll
        end

        # Optimize the extended model using simple gradient descent
        # (A full implementation would use the same optimizer as the main estimation)
        theta_opt = copy(theta_init_extended)
        ofv_opt = compute_extended_ofv(theta_opt)

        # Simple optimization: coordinate descent with line search
        step_sizes = vcat(fill(0.1, n_base_params), fill(0.01, n_covariates))
        max_opt_iter = 100
        tol = 1e-4

        for iter in 1:max_opt_iter
            improved = false
            for p in 1:n_extended
                # Try positive and negative steps
                for direction in [1.0, -1.0]
                    theta_test = copy(theta_opt)
                    theta_test[p] += direction * step_sizes[p]

                    # Ensure positivity for base params
                    if p <= n_base_params && theta_test[p] <= 0
                        continue
                    end

                    ofv_test = compute_extended_ofv(theta_test)
                    if ofv_test < ofv_opt - tol
                        theta_opt = theta_test
                        ofv_opt = ofv_test
                        improved = true
                        step_sizes[p] *= 1.2  # Increase step size on success
                    else
                        step_sizes[p] *= 0.5  # Decrease on failure
                    end
                end
            end

            if !improved
                break
            end
        end

        # Compute approximate SEs using finite differences of Hessian diagonal
        theta_se = zeros(n_extended)
        h = 1e-5
        ofv_center = compute_extended_ofv(theta_opt)

        for p in 1:n_extended
            theta_plus = copy(theta_opt)
            theta_minus = copy(theta_opt)
            theta_plus[p] += h * max(abs(theta_opt[p]), 1.0)
            theta_minus[p] -= h * max(abs(theta_opt[p]), 1.0)

            ofv_plus = compute_extended_ofv(theta_plus)
            ofv_minus = compute_extended_ofv(theta_minus)

            # Second derivative approximation
            d2 = (ofv_plus - 2*ofv_center + ofv_minus) / (h * max(abs(theta_opt[p]), 1.0))^2

            if d2 > 0
                theta_se[p] = sqrt(2.0 / d2)  # SE from curvature of -2LL
            else
                theta_se[p] = 0.0
            end
        end

        return (ofv_opt, theta_opt, theta_se)
    end

    # Run SCM with the covariate-aware estimation function
    if verbose
        println("Starting SCM with $(length(scm_spec.relationships)) candidate relationships")
    end

    result = run_scm(estimate_with_covariates, base_theta, scm_spec)

    if verbose
        println("\nSCM completed:")
        println("  Forward steps: $(result.n_forward_iterations)")
        println("  Backward steps: $(result.n_backward_iterations)")
        println("  Final model: $(length(result.final_model)) covariates")
        println("  ΔOFV: $(round(result.base_ofv - result.final_ofv, digits=2))")
    end

    return result
end

export run_scm_native
