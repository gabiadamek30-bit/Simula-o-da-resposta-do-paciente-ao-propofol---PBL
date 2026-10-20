# Parallel Design

Comprehensive guide for parallel group clinical trial simulation with multiple treatment arms.

---

## Overview

Parallel designs assign subjects to independent treatment groups, enabling comparison of interventions without within-subject correlation.

```julia
using NeoPKPD

design = ParallelDesign(
    n_arms = 3,
    arm_names = ["Placebo", "Low Dose", "High Dose"],
    randomization_ratio = [1, 1, 1]
)
```

---

## ParallelDesign Structure

```julia
struct ParallelDesign <: TrialDesign
    n_arms::Int                           # Number of treatment arms
    arm_names::Vector{String}             # Names for each arm
    randomization_ratio::Vector{Int}      # Allocation ratio
    stratification_factors::Vector{Symbol} # Stratification variables
    block_size::Int                       # Block randomization size
end
```

### Creating Designs

```julia
# Two-arm parallel (1:1)
design = ParallelDesign(
    n_arms = 2,
    arm_names = ["Placebo", "Treatment"],
    randomization_ratio = [1, 1]
)

# Three-arm with 2:1:1 randomization
design = ParallelDesign(
    n_arms = 3,
    arm_names = ["Placebo", "Low", "High"],
    randomization_ratio = [2, 1, 1]
)

# Stratified by sex and age group
design = ParallelDesign(
    n_arms = 2,
    arm_names = ["Control", "Active"],
    randomization_ratio = [1, 1],
    stratification_factors = [:sex, :age_group],
    block_size = 4
)
```

---

## Dosing Regimens

### Standard Frequencies

```julia
# Once daily (QD)
qd_regimen = dosing_qd(
    dose = 100.0,
    days = 28,
    loading_dose = nothing
)

# Twice daily (BID)
bid_regimen = dosing_bid(
    dose = 50.0,
    days = 28,
    loading_dose = 100.0  # Double first dose
)

# Three times daily (TID)
tid_regimen = dosing_tid(
    dose = 25.0,
    days = 14
)

# Four times daily (QID)
qid_regimen = dosing_qid(
    dose = 20.0,
    days = 7
)
```

### Custom Regimens

```julia
# Custom dose times (hours after midnight)
custom_regimen = dosing_custom(
    dose = 100.0,
    days = 28,
    dose_times = [8.0, 20.0]  # 8 AM and 8 PM
)

# Variable doses
variable_regimen = DosingRegimen(
    doses = [100.0, 50.0, 50.0],
    times = [0.0, 12.0, 24.0],
    repeat_days = 2
)
```

### Titration Regimens

```julia
# Gradual dose escalation
titration = TitrationRegimen(
    start_dose = 25.0,
    target_dose = 100.0,
    n_steps = 4,
    days_per_step = 7,
    frequency = :qd,
    loading_dose = nothing,
    maintenance_days = 14
)

# Example schedule:
# Week 1: 25 mg QD
# Week 2: 50 mg QD
# Week 3: 75 mg QD
# Week 4+: 100 mg QD (maintenance)
```

---

## Virtual Population

### Demographic Specification

```julia
# Healthy volunteer population
spec = DemographicSpec(
    age_mean = 35.0,
    age_sd = 10.0,
    age_min = 18.0,
    age_max = 55.0,
    weight_mean = 75.0,
    weight_sd = 12.0,
    weight_min = 50.0,
    weight_max = 100.0,
    female_fraction = 0.5,
    race_distribution = Dict(
        "Caucasian" => 0.70,
        "Asian" => 0.15,
        "Black" => 0.10,
        "Other" => 0.05
    )
)

population = generate_virtual_population(spec, 150)
```

### Pre-Built Specifications

```julia
# Healthy volunteers (Phase I typical)
spec = healthy_volunteer_spec()
# Age: 30 ± 8 years (18-45)
# Weight: 72 ± 10 kg (55-90)
# 50% female

# Patient population
spec = patient_population_spec(:diabetes)
# Age: 55 ± 12 years
# Weight: 90 ± 18 kg
# Includes HbA1c baseline

# Renal impairment
spec = patient_population_spec(:renal_impairment)
# Includes eGFR distribution

# Hepatic impairment
spec = patient_population_spec(:hepatic_impairment)
# Includes Child-Pugh score
```

### Population Summary

```julia
population = generate_virtual_population(spec, 100)
summary = summarize_population(population)

println("Age: $(summary.age_mean) ± $(summary.age_sd) years")
println("Weight: $(summary.weight_mean) ± $(summary.weight_sd) kg")
println("Female: $(summary.female_fraction * 100)%")
println("Race distribution: $(summary.race_distribution)")
```

---

## Trial Specification

### Complete Trial Setup

```julia
using NeoPKPD

# Design
design = ParallelDesign(
    n_arms = 3,
    arm_names = ["Placebo", "50 mg", "100 mg"],
    randomization_ratio = [1, 1, 1]
)

# Dosing regimens per arm
regimens = [
    dosing_qd(dose=0.0, days=28),     # Placebo
    dosing_qd(dose=50.0, days=28),    # Low dose
    dosing_qd(dose=100.0, days=28)    # High dose
]

# Population
spec = healthy_volunteer_spec()
population = generate_virtual_population(spec, 150)

# PK model
model = TwoCompOral()
params = TwoCompOralParams(
    Ka = 1.5,
    CL = 10.0,
    V1 = 50.0,
    Q = 5.0,
    V2 = 100.0
)

# Inter-individual variability
omega = OmegaMatrix([
    0.09 0.03 0.0;    # IIV on CL (30% CV)
    0.03 0.04 0.0;    # IIV on V1 (20% CV), correlated with CL
    0.0  0.0  0.16    # IIV on Ka (40% CV)
])

# Residual error (proportional)
sigma = 0.1  # 10% CV

# Create trial specification
trial = TrialSpec(
    name = "Phase 2 Dose Finding Study",
    design = design,
    regimens = regimens,
    population = population,
    pk_model = model,
    pk_params = params,
    omega = omega,
    sigma = sigma,
    observation_times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
    endpoints = [:cmax, :auc_0_24, :auc_0_inf]
)
```

---

## Trial Simulation

### Running Simulation

```julia
# Simulate trial
result = simulate_trial(trial, seed=12345)

# Access results
println("Study: $(result.name)")
println("Enrolled: $(result.n_enrolled)")
println("Completed: $(result.n_completed)")
println("Dropout rate: $(result.dropout_rate * 100)%")
```

### Arm Results

```julia
for (arm_name, arm) in result.arms
    println("\n$(arm_name):")
    println("  N completed: $(arm.n_completed)")
    println("  AUC: $(arm.mean_auc) ± $(arm.sd_auc)")
    println("  Cmax: $(arm.mean_cmax) ± $(arm.sd_cmax)")
    println("  CV%: $(arm.cv_auc * 100)%")
end
```

### Individual Subject Data

```julia
# Access individual-level data
for subject in result.subjects
    println("Subject $(subject.id):")
    println("  Arm: $(subject.arm)")
    println("  AUC: $(subject.auc)")
    println("  Cmax: $(subject.cmax)")
    println("  Completed: $(subject.completed)")
end
```

---

## Compliance Modeling

### Compliance Patterns

```julia
# Random compliance (probabilistic missed doses)
compliance = ComplianceSpec(
    pattern = :random,
    rate = 0.85  # 85% compliance
)

# Weekend-miss pattern
compliance = ComplianceSpec(
    pattern = :weekend_miss,
    weekday_rate = 0.95,
    weekend_rate = 0.70
)

# Decay pattern (decreasing over time)
compliance = ComplianceSpec(
    pattern = :decay,
    initial_rate = 0.95,
    final_rate = 0.75,
    half_time_days = 14
)

# Early-good pattern
compliance = ComplianceSpec(
    pattern = :early_good,
    initial_rate = 0.95,
    later_rate = 0.80,
    transition_day = 14
)
```

### Adding Compliance to Trial

```julia
trial = TrialSpec(
    name = "Phase 2 with Compliance",
    design = design,
    regimens = regimens,
    population = population,
    pk_model = model,
    pk_params = params,
    omega = omega,
    sigma = sigma,
    compliance = ComplianceSpec(pattern=:random, rate=0.85)
)
```

---

## Dropout Modeling

### Dropout Specification

```julia
# Constant dropout rate
dropout = DropoutSpec(
    rate = 0.15,  # 15% dropout
    pattern = :constant
)

# Time-dependent dropout
dropout = DropoutSpec(
    pattern = :increasing,
    initial_rate = 0.05,
    final_rate = 0.20,
    shape = :linear
)

# Exposure-dependent dropout (adverse events)
dropout = DropoutSpec(
    pattern = :exposure_dependent,
    threshold = 100.0,  # If Cmax > threshold
    rate_below = 0.10,
    rate_above = 0.30
)
```

### Adding Dropout to Trial

```julia
trial = TrialSpec(
    ...,
    dropout = DropoutSpec(rate=0.15, pattern=:constant)
)
```

---

## Statistical Analysis

### Arm Comparisons

```julia
# Compare active to placebo
comparison = compare_arms(
    result.arms["100 mg"].auc_values,
    result.arms["Placebo"].auc_values,
    arm1_name = "100 mg",
    arm2_name = "Placebo",
    alpha = 0.05
)

println("Mean difference: $(comparison.difference)")
println("95% CI: [$(comparison.ci_lower), $(comparison.ci_upper)]")
println("p-value: $(comparison.pvalue)")
println("Significant: $(comparison.significant)")
```

### Dose-Response Analysis

```julia
# Analyze dose-response
doses = [0.0, 50.0, 100.0]
means = [
    result.arms["Placebo"].mean_auc,
    result.arms["50 mg"].mean_auc,
    result.arms["100 mg"].mean_auc
]

# Linear trend test
trend = dose_response_trend(doses, means)
println("Slope: $(trend.slope)")
println("p-value (trend): $(trend.pvalue)")
```

### Responder Analysis

```julia
# Count responders (AUC > threshold)
threshold = 50.0  # Target exposure

for (arm_name, arm) in result.arms
    resp = responder_analysis(
        arm.auc_values,
        threshold = threshold,
        direction = :above,
        confidence = 0.95
    )

    println("$(arm_name):")
    println("  Responders: $(resp.n_responders)/$(resp.n_total)")
    println("  Rate: $(resp.rate * 100)% [$(resp.ci_lower*100), $(resp.ci_upper*100)]")
end
```

---

## Complete Example

```julia
using NeoPKPD

# ======================
# Phase 2 Dose Finding Study
# ======================

# 1. Define design
design = ParallelDesign(
    n_arms = 4,
    arm_names = ["Placebo", "25 mg", "50 mg", "100 mg"],
    randomization_ratio = [1, 1, 1, 1],
    stratification_factors = [:sex],
    block_size = 4
)

# 2. Define dosing regimens
regimens = [
    dosing_qd(dose=0.0, days=28),
    dosing_qd(dose=25.0, days=28),
    dosing_qd(dose=50.0, days=28),
    dosing_qd(dose=100.0, days=28)
]

# 3. Generate population
spec = healthy_volunteer_spec()
population = generate_virtual_population(spec, 200)

# 4. Define PK model
model = OneCompOralFirstOrder()
params = OneCompOralFirstOrderParams(
    Ka = 1.2,
    CL = 8.0,
    V = 60.0
)

omega = OmegaMatrix([
    0.09 0.02;    # IIV CL
    0.02 0.04     # IIV V
])

# 5. Create trial specification
trial = TrialSpec(
    name = "Phase 2 Dose Finding",
    design = design,
    regimens = regimens,
    population = population,
    pk_model = model,
    pk_params = params,
    omega = omega,
    sigma = 0.1,
    observation_times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
    endpoints = [:cmax, :auc_0_24],
    compliance = ComplianceSpec(pattern=:random, rate=0.90),
    dropout = DropoutSpec(rate=0.10, pattern=:constant)
)

# 6. Run simulation
println("Running trial simulation...")
result = simulate_trial(trial, seed=42)

# 7. Report results
println("\n" * "=" ^ 60)
println("PHASE 2 DOSE FINDING STUDY RESULTS")
println("=" ^ 60)

println("\n--- Enrollment Summary ---")
println("Total enrolled: $(result.n_enrolled)")
println("Total completed: $(result.n_completed)")
println("Dropout rate: $(round(result.dropout_rate * 100, digits=1))%")

println("\n--- PK Results by Arm ---")
println("Arm          N    AUC Mean±SD        Cmax Mean±SD")
println("-" ^ 60)

for arm_name in ["Placebo", "25 mg", "50 mg", "100 mg"]
    arm = result.arms[arm_name]
    auc_str = "$(round(arm.mean_auc, digits=1)) ± $(round(arm.sd_auc, digits=1))"
    cmax_str = "$(round(arm.mean_cmax, digits=2)) ± $(round(arm.sd_cmax, digits=2))"
    println("$(rpad(arm_name, 12)) $(arm.n_completed)    $(rpad(auc_str, 18)) $(cmax_str)")
end

# 8. Statistical comparisons
println("\n--- Statistical Comparisons vs Placebo ---")
placebo_auc = result.arms["Placebo"].auc_values

for dose in ["25 mg", "50 mg", "100 mg"]
    comp = compare_arms(
        result.arms[dose].auc_values,
        placebo_auc,
        arm1_name = dose,
        arm2_name = "Placebo"
    )

    sig = comp.significant ? "*" : ""
    println("$(dose) vs Placebo: Δ=$(round(comp.difference, digits=1)), " *
            "95% CI [$(round(comp.ci_lower, digits=1)), $(round(comp.ci_upper, digits=1))], " *
            "p=$(round(comp.pvalue, digits=4))$sig")
end

# 9. Dose-response
println("\n--- Dose-Response Analysis ---")
doses = [0.0, 25.0, 50.0, 100.0]
mean_aucs = [result.arms[arm].mean_auc for arm in ["Placebo", "25 mg", "50 mg", "100 mg"]]
trend = dose_response_trend(doses, mean_aucs)
println("Linear trend slope: $(round(trend.slope, digits=3)) per mg")
println("Trend p-value: $(round(trend.pvalue, digits=4))")

# 10. Responder analysis
println("\n--- Responder Analysis (AUC > 30) ---")
for arm_name in ["Placebo", "25 mg", "50 mg", "100 mg"]
    arm = result.arms[arm_name]
    resp = responder_analysis(arm.auc_values, threshold=30.0, direction=:above)
    println("$(arm_name): $(resp.n_responders)/$(resp.n_total) " *
            "($(round(resp.rate*100, digits=1))%)")
end
```

---

## See Also

- [Crossover Design](crossover.md) - Within-subject designs
- [Dose Escalation](dose-escalation.md) - Phase I designs
- [Power Analysis](power.md) - Sample size calculation
- [Population Generation](../population/index.md) - IIV modeling

