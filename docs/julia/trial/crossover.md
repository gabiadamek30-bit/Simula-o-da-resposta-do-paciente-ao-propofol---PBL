# Crossover Design

Comprehensive guide for crossover clinical trial simulation with within-subject comparisons.

---

## Overview

Crossover designs allow each subject to receive multiple treatments in different periods, enabling within-subject comparisons with reduced variability.

```julia
using NeoPKPD

design = CrossoverDesign(
    n_periods = 2,
    sequences = [["Test", "Reference"], ["Reference", "Test"]],
    washout_days = 14
)
```

---

## CrossoverDesign Structure

```julia
struct CrossoverDesign <: TrialDesign
    n_periods::Int                    # Number of treatment periods
    n_sequences::Int                  # Number of sequences
    sequences::Vector{Vector{String}} # Treatment sequences
    washout_days::Int                 # Washout period between treatments
    period_duration_days::Int         # Duration of each period
end
```

---

## Standard Designs

### 2×2 Crossover

```julia
# Classic AB/BA design
design = crossover_2x2(
    treatments = ["Test", "Reference"],
    washout_days = 14
)

# Equivalent to:
design = CrossoverDesign(
    n_periods = 2,
    sequences = [
        ["Test", "Reference"],    # Sequence 1: TR
        ["Reference", "Test"]     # Sequence 2: RT
    ],
    washout_days = 14
)
```

### 3×3 Latin Square

```julia
# Three treatments, three periods
design = crossover_3x3(
    treatments = ["A", "B", "C"],
    washout_days = 7
)

# Sequences: ABC, BCA, CAB (Latin square)
```

### Williams Design

```julia
# Balanced for first-order carryover
design = williams_design(
    n_treatments = 2,
    washout_days = 14
)
# Returns 2×2 design

design = williams_design(
    n_treatments = 3,
    washout_days = 7
)
# Returns 6 sequences (3×6 design)

design = williams_design(
    n_treatments = 4,
    washout_days = 7
)
# Returns 4 sequences (4×4 design)
```

### Replicate Designs

```julia
# 2×4 replicate (TRTR/RTRT)
design = replicate_crossover_2x4(
    treatments = ["Test", "Reference"],
    washout_days = 7
)

# Sequences:
# TRTR, RTRT

# 3×3 partial replicate
design = partial_replicate_3x3(
    treatments = ["Test", "Reference"],
    washout_days = 7
)

# Sequences:
# TRR, RTR, RRT
```

---

## Trial Specification

### Complete Crossover Setup

```julia
using NeoPKPD

# Design
design = crossover_2x2(
    treatments = ["Test", "Reference"],
    washout_days = 14
)

# Dosing regimens for each treatment
regimens = Dict(
    "Test" => dosing_single(dose=100.0),
    "Reference" => dosing_single(dose=100.0)
)

# Population (equal allocation to sequences)
spec = healthy_volunteer_spec()
population = generate_virtual_population(spec, 24)  # 12 per sequence

# PK model
model = TwoCompOral()
params = TwoCompOralParams(
    Ka = 1.5,
    CL = 10.0,
    V1 = 50.0,
    Q = 5.0,
    V2 = 100.0
)

# IIV
omega = OmegaMatrix([
    0.09 0.0;
    0.0  0.04
])

# Formulation effect (Test may have different Ka)
formulation_effects = Dict(
    "Test" => Dict(:Ka => 1.0),      # No change
    "Reference" => Dict(:Ka => 1.0)  # Baseline
)

# Trial spec
trial = CrossoverTrialSpec(
    name = "BE Crossover Study",
    design = design,
    regimens = regimens,
    population = population,
    pk_model = model,
    pk_params = params,
    omega = omega,
    sigma = 0.1,
    formulation_effects = formulation_effects,
    observation_times = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0]
)
```

---

## Simulation

### Running Crossover Trial

```julia
result = simulate_crossover_trial(trial, seed=42)

# Access by period
for period in 1:result.n_periods
    println("Period $period:")
    for treatment in result.treatments
        period_data = result.period_data[period, treatment]
        println("  $treatment: mean AUC = $(period_data.mean_auc)")
    end
end

# Access by subject
for subject in result.subjects
    println("Subject $(subject.id):")
    for (period, treatment) in enumerate(subject.sequence)
        println("  Period $period ($treatment): AUC = $(subject.auc[period])")
    end
end
```

---

## Washout Period

### Configuring Washout

```julia
# Standard washout (5 half-lives recommended)
design = crossover_2x2(
    treatments = ["Test", "Reference"],
    washout_days = 14  # Adequate for drug with t½ ~ 3 days
)

# Long-acting drug
design = crossover_2x2(
    treatments = ["Test", "Reference"],
    washout_days = 28  # For drug with t½ ~ 5-6 days
)
```

### Washout Verification

```julia
# Check washout adequacy
t_half = 3.0  # days
recommended_washout = 5 * t_half
actual_washout = design.washout_days

if actual_washout >= recommended_washout
    println("✓ Washout adequate ($(actual_washout) ≥ $(recommended_washout) days)")
else
    println("⚠ Washout may be insufficient")
end
```

---

## Period and Sequence Effects

### Testing Period Effect

```julia
result = simulate_crossover_trial(trial, seed=42)

period_effect = test_period_effect(result)

println("Period effect p-value: $(period_effect.pvalue)")
if period_effect.pvalue < 0.05
    println("⚠ Significant period effect detected")
    println("  Period 1 mean: $(period_effect.period1_mean)")
    println("  Period 2 mean: $(period_effect.period2_mean)")
else
    println("✓ No significant period effect")
end
```

### Testing Sequence Effect

```julia
sequence_effect = test_sequence_effect(result)

println("Sequence effect p-value: $(sequence_effect.pvalue)")
if sequence_effect.pvalue < 0.05
    println("⚠ Significant sequence (carryover) effect detected")
else
    println("✓ No significant sequence effect")
end
```

---

## Within-Subject Variability

### Calculating Intra-Subject CV

```julia
result = simulate_crossover_trial(trial, seed=42)

# For replicate designs
cv_within = compute_within_subject_cv(result, :auc)

println("Within-subject CV (AUC): $(cv_within * 100)%")

# Classification
if cv_within < 0.30
    println("Standard variability")
elseif cv_within < 0.40
    println("Moderate variability")
else
    println("Highly variable drug (HVD)")
end
```

### Subject-by-Formulation Interaction

```julia
# For replicate designs
sbf = test_subject_by_formulation(result)

println("Subject × Formulation p-value: $(sbf.pvalue)")
if sbf.pvalue < 0.05
    println("⚠ Significant subject-by-formulation interaction")
end
```

---

## Bioequivalence Analysis

### Standard BE Assessment

```julia
result = simulate_crossover_trial(trial, seed=42)

# Extract paired data
test_auc = result.treatment_data["Test"].auc_values
ref_auc = result.treatment_data["Reference"].auc_values

# Assess BE
be_auc = assess_bioequivalence(
    test = test_auc,
    reference = ref_auc,
    theta1 = 0.80,
    theta2 = 1.25,
    alpha = 0.05
)

println("=== AUC Bioequivalence ===")
println("GMR: $(round(be_auc.gmr * 100, digits=2))%")
println("90% CI: [$(round(be_auc.ci_lower * 100, digits=2))%, $(round(be_auc.ci_upper * 100, digits=2))%]")
println("Within-subject CV: $(round(be_auc.cv_within * 100, digits=1))%")
println("BE demonstrated: $(be_auc.is_bioequivalent)")
```

### Multiple Endpoints

```julia
# Assess both AUC and Cmax
endpoints = [:auc, :cmax]
be_results = Dict()

for endpoint in endpoints
    test_vals = getfield(result.treatment_data["Test"], endpoint * :_values)
    ref_vals = getfield(result.treatment_data["Reference"], endpoint * :_values)

    be_results[endpoint] = assess_bioequivalence(
        test = test_vals,
        reference = ref_vals
    )
end

# Overall BE requires both to pass
overall_be = all(r.is_bioequivalent for r in values(be_results))
println("Overall BE: $(overall_be ? "PASS" : "FAIL")")
```

---

## Highly Variable Drugs

### Reference-Scaled BE (FDA RSABE)

```julia
# For HVD with CV > 30%
rsabe_result = rsabe_analysis(
    result,
    regulatory = :FDA
)

println("Reference CV: $(rsabe_result.cv_reference * 100)%")
println("Scaling applied: $(rsabe_result.scaling_applied)")
println("RSABE criterion met: $(rsabe_result.criterion_met)")
println("Point estimate constraint: $(rsabe_result.point_estimate_ok)")
println("RSABE conclusion: $(rsabe_result.is_bioequivalent)")
```

### ABEL (EMA)

```julia
# Average Bioequivalence with Expanding Limits
abel_result = abel_analysis(
    result,
    regulatory = :EMA
)

println("Reference CV: $(abel_result.cv_reference * 100)%")
println("Widened limits: [$(abel_result.lower_limit*100)%, $(abel_result.upper_limit*100)%]")
println("ABEL conclusion: $(abel_result.is_bioequivalent)")
```

---

## Complete Example

```julia
using NeoPKPD

# ===================================
# Bioequivalence Crossover Study
# ===================================

# 1. Design: 2×2 crossover
design = crossover_2x2(
    treatments = ["Test", "Reference"],
    washout_days = 14
)

# 2. Dosing: Single dose 500 mg
regimens = Dict(
    "Test" => dosing_single(dose=500.0),
    "Reference" => dosing_single(dose=500.0)
)

# 3. Population: 24 healthy volunteers
spec = healthy_volunteer_spec()
population = generate_virtual_population(spec, 24)

# 4. PK Model
model = TwoCompOral()
params = TwoCompOralParams(
    Ka = 1.2,
    CL = 8.0,
    V1 = 50.0,
    Q = 3.0,
    V2 = 80.0
)

# 5. Variability
omega = OmegaMatrix([
    0.09 0.02 0.0;    # CL
    0.02 0.04 0.0;    # V1
    0.0  0.0  0.16    # Ka
])

# Formulation effect: Test has 5% higher Ka (faster absorption)
formulation_effects = Dict(
    "Test" => Dict(:Ka => 1.05),
    "Reference" => Dict(:Ka => 1.0)
)

# 6. Trial spec
trial = CrossoverTrialSpec(
    name = "BE Study - Drug X",
    design = design,
    regimens = regimens,
    population = population,
    pk_model = model,
    pk_params = params,
    omega = omega,
    sigma = 0.1,
    formulation_effects = formulation_effects,
    observation_times = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0],
    endpoints = [:cmax, :auc_0_t, :auc_0_inf, :tmax]
)

# 7. Simulate
println("Simulating crossover trial...")
result = simulate_crossover_trial(trial, seed=12345)

# 8. Report
println("\n" * "=" ^ 60)
println("BIOEQUIVALENCE STUDY RESULTS")
println("=" ^ 60)

println("\n--- Study Summary ---")
println("Subjects enrolled: $(result.n_enrolled)")
println("Subjects completed: $(result.n_completed)")
println("Design: $(result.design.n_periods)×$(result.design.n_sequences) crossover")
println("Washout: $(result.design.washout_days) days")

println("\n--- PK Results by Treatment ---")
for treatment in ["Test", "Reference"]
    data = result.treatment_data[treatment]
    println("\n$treatment:")
    println("  Cmax:     $(round(data.mean_cmax, digits=2)) ± $(round(data.sd_cmax, digits=2))")
    println("  Tmax:     $(round(data.mean_tmax, digits=2)) h")
    println("  AUC0-t:   $(round(data.mean_auc_0_t, digits=1)) ± $(round(data.sd_auc_0_t, digits=1))")
    println("  AUC0-inf: $(round(data.mean_auc_0_inf, digits=1)) ± $(round(data.sd_auc_0_inf, digits=1))")
end

# 9. Period/Sequence effects
println("\n--- Effect Tests ---")
period = test_period_effect(result)
sequence = test_sequence_effect(result)
println("Period effect p-value: $(round(period.pvalue, digits=4))")
println("Sequence effect p-value: $(round(sequence.pvalue, digits=4))")

# 10. Within-subject CV
cv_auc = compute_within_subject_cv(result, :auc_0_inf)
cv_cmax = compute_within_subject_cv(result, :cmax)
println("\n--- Within-Subject Variability ---")
println("CV (AUC): $(round(cv_auc * 100, digits=1))%")
println("CV (Cmax): $(round(cv_cmax * 100, digits=1))%")

# 11. Bioequivalence assessment
println("\n--- Bioequivalence Assessment ---")

be_auc = assess_bioequivalence(
    result.treatment_data["Test"].auc_0_inf_values,
    result.treatment_data["Reference"].auc_0_inf_values
)

be_cmax = assess_bioequivalence(
    result.treatment_data["Test"].cmax_values,
    result.treatment_data["Reference"].cmax_values
)

println("\nAUC0-inf:")
println("  GMR: $(round(be_auc.gmr * 100, digits=2))%")
println("  90% CI: [$(round(be_auc.ci_lower * 100, digits=2))%, $(round(be_auc.ci_upper * 100, digits=2))%]")
println("  BE: $(be_auc.is_bioequivalent ? "PASS" : "FAIL")")

println("\nCmax:")
println("  GMR: $(round(be_cmax.gmr * 100, digits=2))%")
println("  90% CI: [$(round(be_cmax.ci_lower * 100, digits=2))%, $(round(be_cmax.ci_upper * 100, digits=2))%]")
println("  BE: $(be_cmax.is_bioequivalent ? "PASS" : "FAIL")")

# 12. Overall conclusion
overall = be_auc.is_bioequivalent && be_cmax.is_bioequivalent
println("\n" * "=" ^ 60)
println("OVERALL CONCLUSION: Bioequivalence $(overall ? "DEMONSTRATED" : "NOT DEMONSTRATED")")
println("=" ^ 60)
```

---

## See Also

- [Parallel Design](parallel.md) - Independent group designs
- [Power Analysis](power.md) - Sample size for crossover
- [NCA](../nca/index.md) - Non-compartmental analysis

