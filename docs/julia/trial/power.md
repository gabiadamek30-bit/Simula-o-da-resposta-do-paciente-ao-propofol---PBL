# Power Analysis

Comprehensive guide for sample size calculation and power analysis in clinical trial simulation.

---

## Overview

Power analysis enables determination of the sample size required to detect a clinically meaningful effect with adequate statistical power.

```julia
using NeoPKPD

# Calculate sample size for parallel design
result = power_analysis(
    design = :parallel,
    effect_size = 0.5,
    alpha = 0.05,
    power = 0.80,
    cv = 0.30
)

println("Required N per arm: $(result.n_per_arm)")
```

---

## Key Concepts

### Statistical Power

Power is the probability of detecting an effect when one truly exists:

- **Power = 1 - β** (β = Type II error rate)
- Typical target: 80% or 90%
- Higher power requires larger sample size

### Effect Size

The magnitude of the difference to detect:

```julia
# Cohen's d for continuous outcomes
cohens_d = (mean1 - mean2) / pooled_sd

# Geometric Mean Ratio (GMR) for PK
gmr = exp(mean_log_test - mean_log_reference)

# Relative difference
relative_diff = (mean_test - mean_reference) / mean_reference
```

---

## Parallel Design Power

### Two-Arm Comparison

```julia
struct ParallelPowerSpec
    effect_size::Float64        # Expected difference / SD
    alpha::Float64              # Type I error rate
    power::Float64              # Target power
    cv::Float64                 # Coefficient of variation
    allocation_ratio::Float64   # N_treatment / N_control
    test::Symbol                # :two_sided, :one_sided
end

# Create specification
spec = ParallelPowerSpec(
    effect_size = 0.5,          # Medium effect
    alpha = 0.05,
    power = 0.80,
    cv = 0.30,
    allocation_ratio = 1.0,
    test = :two_sided
)

# Calculate sample size
result = calculate_sample_size(spec)
println("N per arm: $(result.n_per_arm)")
println("Total N: $(result.n_total)")
```

### Sample Size Formula

```julia
function sample_size_parallel(
    effect_size::Float64,
    alpha::Float64,
    power::Float64;
    test::Symbol = :two_sided
)
    # Z values
    z_alpha = test == :two_sided ? quantile(Normal(), 1 - alpha/2) : quantile(Normal(), 1 - alpha)
    z_beta = quantile(Normal(), power)

    # Sample size per group
    n = 2 * ((z_alpha + z_beta) / effect_size)^2

    return ceil(Int, n)
end

# Example
n = sample_size_parallel(0.5, 0.05, 0.80)
println("Required n per arm: $n")  # 64 per arm
```

### Multi-Arm Designs

```julia
# Sample size for multiple treatment arms
function sample_size_multiarm(
    n_arms::Int,
    effect_size::Float64,
    alpha::Float64,
    power::Float64;
    comparison::Symbol = :dunnett  # :dunnett, :bonferroni, :none
)
    # Adjust alpha for multiple comparisons
    adjusted_alpha = if comparison == :bonferroni
        alpha / (n_arms - 1)
    elseif comparison == :dunnett
        dunnett_alpha(alpha, n_arms - 1)
    else
        alpha
    end

    # Per-arm sample size
    z_alpha = quantile(Normal(), 1 - adjusted_alpha/2)
    z_beta = quantile(Normal(), power)

    n_per_arm = ceil(Int, 2 * ((z_alpha + z_beta) / effect_size)^2)

    return (
        n_per_arm = n_per_arm,
        n_total = n_per_arm * n_arms,
        adjusted_alpha = adjusted_alpha
    )
end
```

---

## Crossover Design Power

### 2×2 Crossover

```julia
struct CrossoverPowerSpec
    gmr::Float64                # Expected geometric mean ratio
    cv_within::Float64          # Within-subject CV
    theta1::Float64             # Lower BE limit (e.g., 0.80)
    theta2::Float64             # Upper BE limit (e.g., 1.25)
    alpha::Float64              # Type I error rate
    power::Float64              # Target power
end

# BE study sample size
spec = CrossoverPowerSpec(
    gmr = 0.95,                 # Expected GMR of 95%
    cv_within = 0.25,           # 25% within-subject CV
    theta1 = 0.80,
    theta2 = 1.25,
    alpha = 0.05,
    power = 0.80
)

result = calculate_be_sample_size(spec)
println("Required N: $(result.n_total)")
```

### Sample Size Calculation

```julia
function sample_size_crossover_be(
    gmr::Float64,
    cv_within::Float64,
    theta1::Float64,
    theta2::Float64,
    alpha::Float64,
    power::Float64
)
    # Convert CV to log-scale SD
    sigma_w = sqrt(log(1 + cv_within^2))

    # Z values
    z_alpha = quantile(Normal(), 1 - alpha)
    z_beta = quantile(Normal(), power)

    # Distance from GMR to nearest limit
    delta = min(log(gmr) - log(theta1), log(theta2) - log(gmr))

    # Sample size (2×2 crossover)
    n = 2 * ((z_alpha + z_beta) * sigma_w / delta)^2

    return ceil(Int, n)
end

# Example: 25% CV, GMR = 0.95
n = sample_size_crossover_be(0.95, 0.25, 0.80, 1.25, 0.05, 0.80)
println("Required N: $n")  # ~24 subjects
```

### Replicate Crossover Designs

```julia
# Sample size for 4-period replicate design
function sample_size_replicate_be(
    gmr::Float64,
    cv_within::Float64,
    theta1::Float64,
    theta2::Float64,
    alpha::Float64,
    power::Float64;
    design::Symbol = :trtr_rtrt  # :trtr_rtrt, :trt_rtr, etc.
)
    # Design efficiency factors
    efficiency = Dict(
        :trtr_rtrt => 2.0,      # Full replicate
        :trt_rtr => 1.5,        # Partial replicate (TRR, RTR, RRT)
        :tr_rt => 1.0           # Standard 2×2
    )

    # Adjusted sample size
    n_2x2 = sample_size_crossover_be(gmr, cv_within, theta1, theta2, alpha, power)
    n_rep = ceil(Int, n_2x2 / efficiency[design])

    return n_rep
end
```

---

## Highly Variable Drugs

### RSABE Sample Size (FDA)

```julia
function sample_size_rsabe(
    gmr::Float64,
    cv_reference::Float64;
    theta_s::Float64 = 0.8928,  # Scaling factor
    sigma_w0::Float64 = 0.25,   # Regulatory cutoff
    alpha::Float64 = 0.05,
    power::Float64 = 0.80
)
    sigma_wr = sqrt(log(1 + cv_reference^2))

    # Use scaled limits if CV > 30%
    if sigma_wr > sigma_w0
        # Scaled limits
        scaled_limit = exp(theta_s * sigma_wr)
        n = sample_size_scaled_be(gmr, cv_reference, scaled_limit, alpha, power)
    else
        # Standard limits
        n = sample_size_crossover_be(gmr, cv_reference, 0.80, 1.25, alpha, power)
    end

    return n
end

# HVD with 45% CV
n = sample_size_rsabe(0.95, 0.45)
println("RSABE sample size: $n")
```

### ABEL Sample Size (EMA)

```julia
function sample_size_abel(
    gmr::Float64,
    cv_reference::Float64;
    cv_cutoff::Float64 = 0.30,
    max_widening::Float64 = 0.50,  # Up to 69.84%-143.19%
    alpha::Float64 = 0.05,
    power::Float64 = 0.80
)
    if cv_reference > cv_cutoff
        # Calculate widened limits
        k = 0.760  # Regulatory constant
        sigma_wr = sqrt(log(1 + cv_reference^2))
        widening = min(k * sigma_wr, max_widening)

        lower = exp(log(0.80) - widening)
        upper = exp(log(1.25) + widening)

        n = sample_size_crossover_be(gmr, cv_reference, lower, upper, alpha, power)
    else
        n = sample_size_crossover_be(gmr, cv_reference, 0.80, 1.25, alpha, power)
    end

    return n
end
```

---

## Power Curves

### Generating Power Curves

```julia
function power_curve(
    design_spec::PowerSpec,
    n_range::UnitRange{Int}
)
    powers = Float64[]

    for n in n_range
        pow = calculate_power(design_spec, n)
        push!(powers, pow)
    end

    return (n = collect(n_range), power = powers)
end

# Example
spec = ParallelPowerSpec(effect_size=0.5, alpha=0.05, cv=0.30)
curve = power_curve(spec, 10:100)

# Find minimum N for 80% power
min_n = findfirst(p -> p >= 0.80, curve.power)
println("Minimum N for 80% power: $(curve.n[min_n])")
```

### Effect Size Sensitivity

```julia
function effect_sensitivity(
    n::Int,
    alpha::Float64,
    effect_sizes::Vector{Float64}
)
    powers = Float64[]

    for es in effect_sizes
        spec = ParallelPowerSpec(effect_size=es, alpha=alpha)
        pow = calculate_power(spec, n)
        push!(powers, pow)
    end

    return (effect_size = effect_sizes, power = powers)
end

# What effect sizes can we detect with N=50 per arm?
sens = effect_sensitivity(50, 0.05, 0.1:0.1:1.0)
```

---

## Dropout Adjustment

### Adjusting for Anticipated Dropout

```julia
function adjust_for_dropout(
    n_calculated::Int,
    dropout_rate::Float64
)
    n_adjusted = ceil(Int, n_calculated / (1 - dropout_rate))
    return n_adjusted
end

# Example: 20% anticipated dropout
n_base = 64
n_adjusted = adjust_for_dropout(n_base, 0.20)
println("Enroll: $n_adjusted to retain $n_base")  # 80 subjects
```

### Per-Period Dropout (Crossover)

```julia
function adjust_crossover_dropout(
    n_calculated::Int,
    dropout_per_period::Float64,
    n_periods::Int
)
    # Expected completers
    retention = (1 - dropout_per_period)^n_periods

    n_adjusted = ceil(Int, n_calculated / retention)
    return n_adjusted
end

# 2×2 crossover with 5% dropout per period
n = adjust_crossover_dropout(24, 0.05, 2)  # ~27 subjects
```

---

## Simulation-Based Power

### Monte Carlo Power Estimation

```julia
function simulate_power(
    trial_spec::TrialSpec,
    n_subjects::Int,
    n_simulations::Int = 1000;
    alpha::Float64 = 0.05,
    seed::Int = 42
)
    Random.seed!(seed)
    significant_count = 0

    for _ in 1:n_simulations
        # Simulate trial
        result = simulate_trial(trial_spec, n=n_subjects)

        # Perform test
        pvalue = perform_test(result)

        if pvalue < alpha
            significant_count += 1
        end
    end

    power = significant_count / n_simulations
    se = sqrt(power * (1 - power) / n_simulations)

    return (power = power, se = se, ci_lower = power - 1.96*se, ci_upper = power + 1.96*se)
end
```

### Power for BE Studies

```julia
function simulate_be_power(
    pk_model::PKModel,
    pk_params::PKParams,
    omega::Matrix{Float64},
    n_subjects::Int;
    formulation_effect::Dict = Dict(:Ka => 1.05),
    n_simulations::Int = 1000,
    seed::Int = 42
)
    Random.seed!(seed)
    be_pass = 0

    for _ in 1:n_simulations
        # Simulate crossover trial
        result = simulate_crossover_be(
            pk_model, pk_params, omega,
            n = n_subjects,
            formulation_effect = formulation_effect
        )

        # Assess BE
        be = assess_bioequivalence(result.test_auc, result.ref_auc)

        if be.is_bioequivalent
            be_pass += 1
        end
    end

    return be_pass / n_simulations
end

# Example: Power analysis for BE study
power = simulate_be_power(
    TwoCompOral(),
    TwoCompOralParams(Ka=1.5, CL=10.0, V1=50.0, Q=5.0, V2=100.0),
    [0.09 0.0; 0.0 0.04],
    24
)
println("BE power with N=24: $(round(power*100, digits=1))%")
```

---

## Special Populations

### Pediatric Trial Power

```julia
function pediatric_power_adjustment(
    adult_n::Int,
    age_groups::Vector{Tuple{Int, Int}},  # (min_age, max_age)
    group_weights::Vector{Float64};
    variability_inflation::Float64 = 1.2
)
    # Adjust for higher variability in pediatrics
    adjusted_n = ceil(Int, adult_n * variability_inflation)

    # Distribute across age groups
    group_n = [ceil(Int, adjusted_n * w) for w in group_weights]

    return (
        total_n = sum(group_n),
        group_allocation = Dict(zip(age_groups, group_n))
    )
end

# Example: Pediatric study
result = pediatric_power_adjustment(
    64,  # Adult study N
    [(2, 6), (6, 12), (12, 18)],  # Age groups
    [0.3, 0.4, 0.3]  # Allocation weights
)
```

### Renal/Hepatic Impairment

```julia
function impairment_study_size(
    base_effect::Float64,
    cv::Float64;
    impairment_groups::Vector{Symbol} = [:mild, :moderate, :severe],
    alpha::Float64 = 0.05,
    power::Float64 = 0.80
)
    # Typically 6-8 per group for PK studies
    # Calculate based on detecting specified fold-change

    n_per_group = sample_size_parallel(base_effect, alpha, power) ÷ 2

    # Ensure minimum of 6 per group (FDA guidance)
    n_per_group = max(n_per_group, 6)

    return (
        n_per_group = n_per_group,
        n_total = n_per_group * (length(impairment_groups) + 1),  # +1 for normal
        groups = [:normal; impairment_groups]
    )
end
```

---

## Adaptive Sample Size

### Sample Size Re-estimation

```julia
struct AdaptiveSampleSizeSpec
    initial_n::Int              # Initial sample size
    interim_fraction::Float64   # Fraction for interim (e.g., 0.5)
    conditional_power_target::Float64  # Target conditional power
    max_n::Int                  # Maximum total N
    min_increase::Int           # Minimum increase if needed
end

function sample_size_reestimation(
    spec::AdaptiveSampleSizeSpec,
    interim_result::InterimResult
)
    # Calculate observed effect size
    observed_es = interim_result.effect_estimate / interim_result.se

    # Calculate conditional power at current N
    n_remaining = spec.initial_n - interim_result.n
    cp_current = conditional_power(observed_es, interim_result.n, n_remaining)

    if cp_current >= spec.conditional_power_target
        return spec.initial_n  # No increase needed
    end

    # Calculate new sample size for target conditional power
    new_n = find_n_for_conditional_power(
        observed_es,
        interim_result.n,
        spec.conditional_power_target
    )

    new_n = min(new_n, spec.max_n)
    new_n = max(new_n, spec.initial_n + spec.min_increase)

    return new_n
end
```

---

## Complete Example

```julia
using NeoPKPD

# ============================================
# Power Analysis for Phase III Parallel Study
# ============================================

println("=== Power Analysis ===\n")

# 1. Study Parameters
println("--- Study Parameters ---")
effect_size = 0.40            # Expected treatment effect (Cohen's d)
cv = 0.35                     # Between-subject CV
alpha = 0.05                  # Type I error
target_power = 0.80           # Target power
dropout_rate = 0.15           # Expected dropout

println("Effect size (Cohen's d): $effect_size")
println("Between-subject CV: $(cv * 100)%")
println("Alpha: $alpha")
println("Target power: $(target_power * 100)%")
println("Expected dropout: $(dropout_rate * 100)%")

# 2. Calculate base sample size
println("\n--- Sample Size Calculation ---")
n_per_arm = sample_size_parallel(effect_size, alpha, target_power)
println("Base sample size: $n_per_arm per arm")

# 3. Adjust for dropout
n_adjusted = adjust_for_dropout(n_per_arm, dropout_rate)
println("Adjusted for dropout: $n_adjusted per arm")
println("Total enrollment: $(2 * n_adjusted)")

# 4. Power curve
println("\n--- Power by Sample Size ---")
println("N/arm    Power")
println("-" ^ 20)

for n in [30, 40, 50, 60, 70, 80, 90, 100]
    pow = calculate_power_parallel(effect_size, alpha, n)
    marker = pow >= 0.80 ? "*" : ""
    @printf("%4d    %5.1f%%%s\n", n, pow * 100, marker)
end

# 5. Effect size sensitivity
println("\n--- Detectable Effect Sizes ---")
println("With N = $n_per_arm per arm:")

for power_target in [0.70, 0.80, 0.90]
    min_es = minimum_detectable_effect(n_per_arm, alpha, power_target)
    @printf("  For %d%% power: Cohen's d = %.3f\n", Int(power_target*100), min_es)
end

# 6. Simulation-based verification
println("\n--- Simulation Verification ---")
println("Running 1000 simulations...")

# Define trial
design = ParallelDesign(
    n_arms = 2,
    arm_names = ["Placebo", "Treatment"],
    randomization_ratio = [1, 1]
)

pk_model = OneCompOral()
pk_params = OneCompOralParams(Ka=1.5, CL=10.0, V=60.0)
omega = [0.09 0.0; 0.0 0.04]

# Simulate power
simulated_power = simulate_trial_power(
    design = design,
    pk_model = pk_model,
    pk_params = pk_params,
    omega = omega,
    n_per_arm = n_per_arm,
    treatment_effect = effect_size,
    n_simulations = 1000
)

println("Simulated power: $(round(simulated_power.power * 100, digits=1))%")
println("95% CI: [$(round(simulated_power.ci_lower * 100, digits=1))%, " *
        "$(round(simulated_power.ci_upper * 100, digits=1))%]")

# 7. Summary
println("\n" * "=" ^ 50)
println("RECOMMENDATION")
println("=" ^ 50)
println("\nEnroll $(2 * n_adjusted) subjects ($(n_adjusted) per arm)")
println("Expected completers: ~$(2 * n_per_arm)")
println("Achieves $(Int(target_power * 100))% power to detect effect size = $effect_size")
```

---

## Power Tables

### Parallel Design Reference

| Effect Size | CV | N per Arm (80% power) | N per Arm (90% power) |
|-------------|-----|----------------------|----------------------|
| 0.3 | 30% | 176 | 235 |
| 0.4 | 30% | 99 | 132 |
| 0.5 | 30% | 64 | 85 |
| 0.6 | 30% | 44 | 59 |
| 0.3 | 40% | 176 | 235 |
| 0.4 | 40% | 99 | 132 |
| 0.5 | 40% | 64 | 85 |

### Crossover BE Reference

| CV Within | GMR | N (80% power) | N (90% power) |
|-----------|-----|---------------|---------------|
| 15% | 0.95 | 10 | 14 |
| 20% | 0.95 | 16 | 22 |
| 25% | 0.95 | 24 | 32 |
| 30% | 0.95 | 36 | 48 |
| 25% | 1.00 | 18 | 24 |
| 25% | 0.90 | 40 | 54 |

---

## See Also

- [Parallel Design](parallel.md) - Parallel trial simulation
- [Crossover Design](crossover.md) - Crossover trial simulation
- [Dose Escalation](dose-escalation.md) - Phase I designs
- [Bioequivalence](../nca/bioequivalence.md) - BE analysis methods
