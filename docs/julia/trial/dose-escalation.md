# Dose Escalation Designs

Comprehensive guide for Phase I dose escalation trial simulation with multiple escalation rules.

---

## Overview

Dose escalation designs are used in Phase I trials to identify the maximum tolerated dose (MTD) while minimizing patient exposure to toxic doses.

```julia
using NeoPKPD

# 3+3 design
design = ThreePlusThree(
    doses = [10.0, 25.0, 50.0, 100.0, 200.0],
    target_dlt_rate = 0.33,
    starting_dose_index = 1
)
```

---

## Escalation Algorithms

### 3+3 Design

The classic rule-based design for oncology Phase I trials.

```julia
struct ThreePlusThree <: DoseEscalationDesign
    doses::Vector{Float64}           # Available dose levels
    target_dlt_rate::Float64         # Target DLT rate (typically 0.33)
    starting_dose_index::Int         # Index of starting dose
    max_subjects_per_dose::Int       # Maximum per cohort (typically 6)
end
```

#### Creating 3+3 Design

```julia
# Standard 3+3
design = ThreePlusThree(
    doses = [10.0, 25.0, 50.0, 100.0, 200.0, 400.0],
    target_dlt_rate = 0.33,
    starting_dose_index = 1
)

# With custom parameters
design = ThreePlusThree(
    doses = [1.0, 3.0, 10.0, 30.0, 100.0],  # Log-spaced
    target_dlt_rate = 0.33,
    starting_dose_index = 1,
    max_subjects_per_dose = 6
)
```

#### Decision Rules

```julia
# Automatic decision based on DLT count
function three_plus_three_decision(n_dlt::Int, n_patients::Int)
    if n_patients == 3
        if n_dlt == 0
            return :escalate       # 0/3 → escalate
        elseif n_dlt == 1
            return :expand         # 1/3 → expand to 6
        else
            return :deescalate     # 2-3/3 → de-escalate
        end
    elseif n_patients == 6
        if n_dlt <= 1
            return :escalate       # 0-1/6 → escalate
        else
            return :deescalate     # 2+/6 → de-escalate, MTD = previous
        end
    end
end
```

---

### Modified Toxicity Probability Interval (mTPI)

Model-assisted design that uses probability intervals.

```julia
struct MTPI <: DoseEscalationDesign
    doses::Vector{Float64}
    target_toxicity::Float64         # Target DLT probability (e.g., 0.25)
    epsilon1::Float64                # Lower equivalence margin
    epsilon2::Float64                # Upper equivalence margin
    prior_alpha::Float64             # Beta prior alpha
    prior_beta::Float64              # Beta prior beta
    starting_dose_index::Int
    max_sample_size::Int
end
```

#### Creating mTPI Design

```julia
# Standard mTPI
design = MTPI(
    doses = [10.0, 25.0, 50.0, 100.0, 200.0],
    target_toxicity = 0.25,
    epsilon1 = 0.05,
    epsilon2 = 0.05,
    prior_alpha = 1.0,
    prior_beta = 1.0,
    starting_dose_index = 1,
    max_sample_size = 36
)

# mTPI-2 (improved version)
design = MTPI2(
    doses = [10.0, 25.0, 50.0, 100.0, 200.0],
    target_toxicity = 0.25,
    epsilon1 = 0.05,
    epsilon2 = 0.05,
    starting_dose_index = 1,
    max_sample_size = 36
)
```

#### mTPI Decision Making

```julia
# Decision based on posterior probability intervals
function mtpi_decision(n_dlt::Int, n_patients::Int, design::MTPI)
    # Posterior Beta distribution
    alpha_post = design.prior_alpha + n_dlt
    beta_post = design.prior_beta + n_patients - n_dlt

    # Calculate unit probability mass (UPM)
    p_target = design.target_toxicity
    e1, e2 = design.epsilon1, design.epsilon2

    # Underdosing interval: (0, p_target - e1)
    upm_under = cdf(Beta(alpha_post, beta_post), p_target - e1) / (p_target - e1)

    # Target interval: (p_target - e1, p_target + e2)
    prob_target = cdf(Beta(alpha_post, beta_post), p_target + e2) -
                  cdf(Beta(alpha_post, beta_post), p_target - e1)
    upm_target = prob_target / (e1 + e2)

    # Overdosing interval: (p_target + e2, 1)
    upm_over = (1 - cdf(Beta(alpha_post, beta_post), p_target + e2)) / (1 - p_target - e2)

    # Decision based on maximum UPM
    if upm_under > upm_target && upm_under > upm_over
        return :escalate
    elseif upm_over > upm_target && upm_over > upm_under
        return :deescalate
    else
        return :stay
    end
end
```

---

### Continual Reassessment Method (CRM)

Model-based design using a dose-toxicity model.

```julia
struct CRM <: DoseEscalationDesign
    doses::Vector{Float64}
    target_toxicity::Float64         # Target DLT rate
    skeleton::Vector{Float64}        # Prior toxicity probabilities
    model::Symbol                    # :power, :logistic, :empiric
    prior_mean::Float64              # Prior for model parameter
    prior_sd::Float64                # Prior SD
    starting_dose_index::Int
    max_sample_size::Int
    cohort_size::Int
end
```

#### Creating CRM Design

```julia
# Power model CRM
design = CRM(
    doses = [10.0, 25.0, 50.0, 100.0, 200.0],
    target_toxicity = 0.25,
    skeleton = [0.05, 0.10, 0.25, 0.40, 0.55],  # Prior p(DLT)
    model = :power,
    prior_mean = 0.0,
    prior_sd = 1.34,
    starting_dose_index = 1,
    max_sample_size = 30,
    cohort_size = 1
)

# Logistic model CRM
design = CRM(
    doses = [10.0, 25.0, 50.0, 100.0, 200.0],
    target_toxicity = 0.25,
    skeleton = [0.05, 0.10, 0.25, 0.40, 0.55],
    model = :logistic,
    prior_mean = 0.0,
    prior_sd = 2.0,
    starting_dose_index = 2,  # Start at dose 2 for safety
    max_sample_size = 30,
    cohort_size = 3
)
```

#### CRM Model Specification

```julia
# Power model: p(DLT|dose) = skeleton[i]^exp(a)
function power_model(skeleton::Vector{Float64}, a::Float64)
    return skeleton .^ exp(a)
end

# Logistic model: logit(p) = a + b * log(dose)
function logistic_model(doses::Vector{Float64}, a::Float64, b::Float64)
    logit_p = a .+ b .* log.(doses)
    return 1.0 ./ (1.0 .+ exp.(-logit_p))
end
```

#### CRM Decision Algorithm

```julia
function crm_next_dose(
    design::CRM,
    dose_history::Vector{Int},
    dlt_history::Vector{Bool}
)
    # Fit model using Bayesian update
    posterior = fit_crm_model(
        design.skeleton,
        dose_history,
        dlt_history,
        design.prior_mean,
        design.prior_sd
    )

    # Estimate toxicity at each dose
    estimated_toxicity = compute_toxicity_estimates(
        design.skeleton,
        posterior
    )

    # Select dose closest to target
    distances = abs.(estimated_toxicity .- design.target_toxicity)
    next_dose_index = argmin(distances)

    # Safety rule: don't skip more than 1 dose
    current_dose = isempty(dose_history) ? design.starting_dose_index : dose_history[end]
    if next_dose_index > current_dose + 1
        next_dose_index = current_dose + 1
    end

    return next_dose_index
end
```

---

### Bayesian Optimal Interval (BOIN)

Model-assisted design with optimal boundaries.

```julia
struct BOIN <: DoseEscalationDesign
    doses::Vector{Float64}
    target_toxicity::Float64
    p1::Float64                      # Highest acceptable DLT rate
    p2::Float64                      # Lowest unacceptable DLT rate
    starting_dose_index::Int
    max_sample_size::Int
    cohort_size::Int
end
```

#### Creating BOIN Design

```julia
# Standard BOIN
design = BOIN(
    doses = [10.0, 25.0, 50.0, 100.0, 200.0],
    target_toxicity = 0.30,
    p1 = 0.25,                       # Acceptable toxicity
    p2 = 0.35,                       # Unacceptable toxicity
    starting_dose_index = 1,
    max_sample_size = 36,
    cohort_size = 3
)

# Calculate optimal boundaries
lambda_e = log((1 - design.p1) / (1 - design.target_toxicity)) /
           log(design.target_toxicity * (1 - design.p1) / (design.p1 * (1 - design.target_toxicity)))

lambda_d = log((1 - design.target_toxicity) / (1 - design.p2)) /
           log(design.p2 * (1 - design.target_toxicity) / (design.target_toxicity * (1 - design.p2)))

println("Escalation boundary (λe): $(round(lambda_e, digits=3))")
println("De-escalation boundary (λd): $(round(lambda_d, digits=3))")
```

#### BOIN Decision Rules

```julia
function boin_decision(
    n_dlt::Int,
    n_patients::Int,
    lambda_e::Float64,
    lambda_d::Float64
)
    observed_rate = n_dlt / n_patients

    if observed_rate <= lambda_e
        return :escalate
    elseif observed_rate >= lambda_d
        return :deescalate
    else
        return :stay
    end
end
```

---

## Dose-Toxicity Modeling

### Exposure-Toxicity Relationship

```julia
# Link PK exposure to DLT probability
struct ExposureToxicityModel
    pk_model::PKModel
    pk_params::PKParams
    exposure_metric::Symbol          # :cmax, :auc, :cmin
    toxicity_model::Symbol           # :logistic, :probit
    toxicity_params::Vector{Float64} # [intercept, slope]
end

# Create exposure-toxicity model
exp_tox = ExposureToxicityModel(
    pk_model = TwoCompOral(),
    pk_params = TwoCompOralParams(Ka=1.5, CL=10.0, V1=50.0, Q=5.0, V2=100.0),
    exposure_metric = :cmax,
    toxicity_model = :logistic,
    toxicity_params = [-3.0, 0.05]   # P(DLT) = logistic(-3 + 0.05*Cmax)
)
```

### Simulating DLT Events

```julia
# Simulate DLT based on exposure
function simulate_dlt(
    dose::Float64,
    exp_tox::ExposureToxicityModel;
    seed::Int = nothing
)
    if seed !== nothing
        Random.seed!(seed)
    end

    # Simulate PK
    pk_result = simulate_pk(
        exp_tox.pk_model,
        exp_tox.pk_params,
        dose,
        times = 0.0:0.1:72.0
    )

    # Calculate exposure metric
    exposure = if exp_tox.exposure_metric == :cmax
        maximum(pk_result.concentrations)
    elseif exp_tox.exposure_metric == :auc
        compute_auc(pk_result)
    else
        minimum(pk_result.concentrations[pk_result.times .> 0])
    end

    # Calculate DLT probability
    a, b = exp_tox.toxicity_params
    logit_p = a + b * exposure
    p_dlt = 1.0 / (1.0 + exp(-logit_p))

    # Simulate DLT
    return rand() < p_dlt
end
```

---

## Trial Simulation

### Running Dose Escalation Trial

```julia
# Simulate complete Phase I trial
function simulate_dose_escalation(
    design::DoseEscalationDesign,
    exp_tox::ExposureToxicityModel;
    n_simulations::Int = 1000,
    seed::Int = 42
)
    Random.seed!(seed)

    mtd_selections = zeros(Int, length(design.doses) + 1)  # +1 for "none"
    n_patients_per_sim = Int[]
    n_dlts_per_sim = Int[]

    for sim in 1:n_simulations
        result = run_single_escalation(design, exp_tox)
        push!(n_patients_per_sim, result.n_patients)
        push!(n_dlts_per_sim, result.n_dlts)

        if result.mtd_index === nothing
            mtd_selections[end] += 1
        else
            mtd_selections[result.mtd_index] += 1
        end
    end

    return DoseEscalationSimResult(
        design = design,
        n_simulations = n_simulations,
        mtd_selection_prob = mtd_selections ./ n_simulations,
        mean_patients = mean(n_patients_per_sim),
        mean_dlts = mean(n_dlts_per_sim)
    )
end
```

### Single Trial Run

```julia
function run_single_escalation(
    design::ThreePlusThree,
    exp_tox::ExposureToxicityModel
)
    current_dose_idx = design.starting_dose_index
    dose_history = Int[]
    dlt_history = Bool[]
    dose_data = Dict{Int, Tuple{Int, Int}}()  # dose_idx => (n_patients, n_dlt)

    while true
        # Treat cohort
        n_treat = 3
        n_dlt = sum(simulate_dlt(design.doses[current_dose_idx], exp_tox) for _ in 1:n_treat)

        # Update history
        append!(dose_history, fill(current_dose_idx, n_treat))
        append!(dlt_history, [n_dlt > i for i in 0:n_treat-1][1:n_treat])

        # Update dose data
        if haskey(dose_data, current_dose_idx)
            prev = dose_data[current_dose_idx]
            dose_data[current_dose_idx] = (prev[1] + n_treat, prev[2] + n_dlt)
        else
            dose_data[current_dose_idx] = (n_treat, n_dlt)
        end

        total_n, total_dlt = dose_data[current_dose_idx]

        # Make decision
        decision = three_plus_three_decision(total_dlt, total_n)

        if decision == :escalate
            if current_dose_idx < length(design.doses)
                current_dose_idx += 1
            else
                # At max dose, declare MTD
                return (mtd_index = current_dose_idx,
                        n_patients = length(dose_history),
                        n_dlts = sum(dlt_history))
            end
        elseif decision == :expand
            continue  # Treat 3 more at same dose
        else  # :deescalate
            mtd_idx = current_dose_idx > 1 ? current_dose_idx - 1 : nothing
            return (mtd_index = mtd_idx,
                    n_patients = length(dose_history),
                    n_dlts = sum(dlt_history))
        end

        # Safety stop
        if length(dose_history) >= 50
            return (mtd_index = nothing,
                    n_patients = length(dose_history),
                    n_dlts = sum(dlt_history))
        end
    end
end
```

---

## Operating Characteristics

### Evaluating Design Performance

```julia
# Calculate operating characteristics
function operating_characteristics(
    design::DoseEscalationDesign,
    true_toxicity::Vector{Float64};
    n_simulations::Int = 1000,
    seed::Int = 42
)
    # Find true MTD
    true_mtd = findfirst(x -> x >= design.target_toxicity, true_toxicity)

    results = simulate_dose_escalation_with_true_tox(
        design, true_toxicity, n_simulations, seed
    )

    # Calculate metrics
    correct_selection = results.mtd_selection_prob[true_mtd]
    overdose_rate = sum(results.mtd_selection_prob[true_mtd+1:end])

    # Average DLT rate experienced
    avg_dlt_rate = results.mean_dlts / results.mean_patients

    println("=== Operating Characteristics ===")
    println("True MTD: Dose $(true_mtd) ($(design.doses[true_mtd]))")
    println("Correct MTD selection: $(round(correct_selection * 100, digits=1))%")
    println("Overdose selection rate: $(round(overdose_rate * 100, digits=1))%")
    println("Average patients per trial: $(round(results.mean_patients, digits=1))")
    println("Average DLT rate: $(round(avg_dlt_rate * 100, digits=1))%")

    return results
end

# Example: Compare designs
true_toxicity = [0.05, 0.10, 0.25, 0.40, 0.55]

# 3+3 design
design_3p3 = ThreePlusThree(
    doses = [10.0, 25.0, 50.0, 100.0, 200.0],
    target_dlt_rate = 0.33
)
oc_3p3 = operating_characteristics(design_3p3, true_toxicity)

# CRM design
design_crm = CRM(
    doses = [10.0, 25.0, 50.0, 100.0, 200.0],
    target_toxicity = 0.25,
    skeleton = [0.05, 0.10, 0.25, 0.40, 0.55]
)
oc_crm = operating_characteristics(design_crm, true_toxicity)
```

### Design Comparison

```julia
# Compare multiple designs
function compare_escalation_designs(
    designs::Vector{<:DoseEscalationDesign},
    true_toxicity::Vector{Float64};
    n_simulations::Int = 1000
)
    results = []

    for design in designs
        oc = operating_characteristics(design, true_toxicity, n_simulations=n_simulations)
        push!(results, (design=design, oc=oc))
    end

    println("\n=== Design Comparison ===")
    println("Design          Correct%  Overdose%  Avg N  Avg DLT%")
    println("-" ^ 55)

    for r in results
        name = string(typeof(r.design).name.name)[1:12]
        correct = r.oc.correct_selection * 100
        overdose = r.oc.overdose_rate * 100
        avg_n = r.oc.mean_patients
        avg_dlt = r.oc.mean_dlts / r.oc.mean_patients * 100

        @printf("%-15s %6.1f    %6.1f     %5.1f  %6.1f\n",
                name, correct, overdose, avg_n, avg_dlt)
    end

    return results
end
```

---

## Safety Constraints

### Stopping Rules

```julia
# Define safety stopping rules
struct SafetyStoppingRules
    max_total_dlt::Int               # Stop if total DLTs exceed
    max_dlt_rate::Float64            # Stop if DLT rate at dose 1 exceeds
    min_patients_dose1::Int          # Min patients at dose 1 before stopping
    posterior_prob_threshold::Float64 # P(toxicity > target) threshold
end

rules = SafetyStoppingRules(
    max_total_dlt = 10,
    max_dlt_rate = 0.50,
    min_patients_dose1 = 6,
    posterior_prob_threshold = 0.95
)

# Check if should stop for safety
function check_safety_stop(
    dose_data::Dict{Int, Tuple{Int, Int}},
    rules::SafetyStoppingRules,
    design::DoseEscalationDesign
)
    # Total DLT check
    total_dlt = sum(d[2] for d in values(dose_data))
    if total_dlt >= rules.max_total_dlt
        return true, "Exceeded maximum total DLTs"
    end

    # Dose 1 toxicity check
    if haskey(dose_data, 1)
        n, dlt = dose_data[1]
        if n >= rules.min_patients_dose1 && dlt / n >= rules.max_dlt_rate
            return true, "Dose 1 toxicity too high"
        end
    end

    return false, ""
end
```

### Dose Limiting Rules

```julia
# Prevent unsafe escalation
struct EscalationConstraints
    skip_dose_allowed::Bool          # Can skip intermediate doses?
    max_skip::Int                    # Maximum doses to skip
    min_patients_before_escalate::Int # Min at current dose before escalating
    require_dlt_free::Bool           # Must have 0 DLTs to escalate?
end

constraints = EscalationConstraints(
    skip_dose_allowed = false,
    max_skip = 0,
    min_patients_before_escalate = 3,
    require_dlt_free = false
)
```

---

## Complete Example

```julia
using NeoPKPD

# =======================================
# Phase I Dose Escalation Trial Simulation
# =======================================

# 1. Define dose levels
doses = [10.0, 25.0, 50.0, 100.0, 200.0, 400.0]  # mg
true_toxicity = [0.02, 0.05, 0.15, 0.30, 0.50, 0.70]

println("=== Phase I Trial Simulation ===")
println("\nDose Levels and True Toxicity:")
for (i, (d, t)) in enumerate(zip(doses, true_toxicity))
    println("  Dose $i: $d mg → P(DLT) = $(t * 100)%")
end

# 2. Create designs
design_3p3 = ThreePlusThree(
    doses = doses,
    target_dlt_rate = 0.33,
    starting_dose_index = 1
)

design_crm = CRM(
    doses = doses,
    target_toxicity = 0.25,
    skeleton = [0.05, 0.10, 0.20, 0.35, 0.50, 0.65],
    model = :power,
    starting_dose_index = 1,
    max_sample_size = 30,
    cohort_size = 1
)

design_boin = BOIN(
    doses = doses,
    target_toxicity = 0.25,
    p1 = 0.20,
    p2 = 0.30,
    starting_dose_index = 1,
    max_sample_size = 36,
    cohort_size = 3
)

# 3. Run simulations
println("\n--- Running 1000 Simulations ---")

results_3p3 = simulate_with_true_toxicity(design_3p3, true_toxicity, n_sim=1000)
results_crm = simulate_with_true_toxicity(design_crm, true_toxicity, n_sim=1000)
results_boin = simulate_with_true_toxicity(design_boin, true_toxicity, n_sim=1000)

# 4. Display results
println("\n--- MTD Selection Probabilities ---")
println("Dose    True P(DLT)   3+3      CRM      BOIN")
println("-" ^ 55)

for i in 1:length(doses)
    @printf("%3d mg    %5.1f%%    %5.1f%%   %5.1f%%   %5.1f%%\n",
            doses[i], true_toxicity[i]*100,
            results_3p3.mtd_prob[i]*100,
            results_crm.mtd_prob[i]*100,
            results_boin.mtd_prob[i]*100)
end

# 5. Summary statistics
println("\n--- Summary Statistics ---")
println("                     3+3      CRM      BOIN")
println("-" ^ 45)

# True MTD is dose 4 (30% DLT rate, closest to 25% target)
true_mtd_idx = 4
@printf("Correct MTD (%%):    %5.1f    %5.1f    %5.1f\n",
        results_3p3.mtd_prob[true_mtd_idx]*100,
        results_crm.mtd_prob[true_mtd_idx]*100,
        results_boin.mtd_prob[true_mtd_idx]*100)

overdose_3p3 = sum(results_3p3.mtd_prob[5:end]) * 100
overdose_crm = sum(results_crm.mtd_prob[5:end]) * 100
overdose_boin = sum(results_boin.mtd_prob[5:end]) * 100
@printf("Overdose (%%):       %5.1f    %5.1f    %5.1f\n",
        overdose_3p3, overdose_crm, overdose_boin)

@printf("Avg patients:       %5.1f    %5.1f    %5.1f\n",
        results_3p3.mean_n, results_crm.mean_n, results_boin.mean_n)

@printf("Avg DLTs:           %5.1f    %5.1f    %5.1f\n",
        results_3p3.mean_dlt, results_crm.mean_dlt, results_boin.mean_dlt)

# 6. Toxicity at recommended MTD
println("\n--- Safety Analysis ---")
for (name, result) in [("3+3", results_3p3), ("CRM", results_crm), ("BOIN", results_boin)]
    avg_tox = sum(result.mtd_prob[i] * true_toxicity[i] for i in 1:length(doses))
    println("$name: Average toxicity at selected MTD = $(round(avg_tox*100, digits=1))%")
end
```

---

## See Also

- [Parallel Design](parallel.md) - Phase II parallel designs
- [Crossover Design](crossover.md) - Crossover study designs
- [Power Analysis](power.md) - Sample size calculation
- [Population Simulation](../population/index.md) - IIV in PK
