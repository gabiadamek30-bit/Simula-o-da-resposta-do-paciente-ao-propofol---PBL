# Bioequivalence Analysis

Comprehensive documentation for bioequivalence (BE) assessment following FDA, EMA, and Health Canada guidance.

---

## Overview

Bioequivalence analysis determines whether two drug formulations have comparable rate and extent of absorption within predefined limits (typically 80-125%).

---

## Quick Start

```julia
using NeoPKPD

# Test vs Reference NCA results
test_auc = [45.2, 52.1, 38.7, 61.3, 49.8]
ref_auc = [48.1, 55.3, 41.2, 58.9, 52.4]
test_cmax = [8.5, 9.2, 7.1, 10.8, 8.9]
ref_cmax = [9.1, 9.8, 7.8, 10.2, 9.5]

# Calculate 90% CI for AUC
ci_auc = bioequivalence_90ci(test_auc, ref_auc)
println("AUC GMR: $(ci_auc.gmr)")
println("90% CI: [$(ci_auc.lower), $(ci_auc.upper)]")
println("BE met: $(ci_auc.be_met)")
```

---

## Study Designs

### 2x2 Crossover (Standard)

```julia
design = Crossover2x2()

# Periods: 1, 2
# Sequences: TR, RT
# Each subject receives both test and reference
```

### 2x4 Replicate Crossover

```julia
design = Crossover2x4()

# Periods: 1, 2, 3, 4
# Sequences: TRTR, RTRT
# For highly variable drugs
```

### 3-Way Crossover

```julia
design = Crossover3Way()

# For three treatments (e.g., test, reference, fasted reference)
```

### Parallel Group

```julia
design = ParallelGroup()

# Used when crossover not feasible (long half-life drugs)
```

---

## 90% Confidence Intervals

### Basic 90% CI Calculation

```julia
# Log-transformed data analysis
ci = bioequivalence_90ci(test_values, ref_values)

println("Geometric Mean Ratio: $(round(ci.gmr * 100, digits=2))%")
println("90% CI: [$(round(ci.lower * 100, digits=2))%, $(round(ci.upper * 100, digits=2))%]")
println("BE criteria (80-125%): $(ci.be_met ? "MET" : "NOT MET")")
```

### CI Result Structure

```julia
struct BioequivalenceResult
    gmr::Float64           # Geometric mean ratio (Test/Reference)
    lower::Float64         # 90% CI lower bound
    upper::Float64         # 90% CI upper bound
    be_met::Bool           # Within 80-125%
    cv_within::Float64     # Within-subject CV%
    n_subjects::Int        # Number of subjects
    power::Float64         # Statistical power achieved
end
```

### Crossover Analysis with Period Effects

```julia
# Full crossover analysis with period and sequence effects
result = analyze_crossover_be(
    data,
    treatment_col = :formulation,
    subject_col = :subject_id,
    period_col = :period,
    sequence_col = :sequence,
    response_col = :auc
)

println("Period effect p-value: $(result.period_pvalue)")
println("Sequence effect p-value: $(result.sequence_pvalue)")
println("Carryover effect p-value: $(result.carryover_pvalue)")
```

---

## TOST Analysis

Two One-Sided Tests for equivalence:

```julia
# TOST with custom bounds
tost = tost_analysis(
    test_values,
    ref_values,
    lower_bound = 0.80,
    upper_bound = 1.25,
    alpha = 0.05
)

println("TOST p-value (lower): $(tost.p_lower)")
println("TOST p-value (upper): $(tost.p_upper)")
println("Overall BE conclusion: $(tost.be_concluded)")
```

### Regulatory Acceptance Limits

| Parameter | FDA | EMA | Health Canada |
|-----------|-----|-----|---------------|
| AUC | 80-125% | 80-125% | 80-125% |
| Cmax | 80-125% | 80-125% | 80-125% |
| AUC (HVD) | 80-125% or scaled | Widened | 80-125% |
| Cmax (HVD) | 80-125% | Widened | 80-125% |

---

## Geometric Mean Ratio

```julia
# Calculate GMR
gmr = geometric_mean_ratio(test_values, ref_values)
println("GMR: $(round(gmr * 100, digits=2))%")

# Point estimate
point_estimate = exp(mean(log.(test_values)) - mean(log.(ref_values)))
```

---

## Within-Subject Variability

### CV% Calculation

```julia
# From crossover data
cv_within = within_subject_cv(data, :auc, :subject_id)
println("Within-subject CV: $(round(cv_within * 100, digits=1))%")

# Classification
if cv_within < 0.30
    println("Standard variability drug")
elseif cv_within < 0.40
    println("Moderately variable drug")
else
    println("Highly variable drug (HVD)")
end
```

### Intra-Subject Variability from Replicate Design

```julia
# From 2x4 replicate crossover
result = replicate_be_analysis(data, design=Crossover2x4())

println("Reference CV: $(result.cv_reference)%")
println("Test CV: $(result.cv_test)%")
println("Subject-by-formulation interaction: $(result.sbf_interaction)")
```

---

## Reference-Scaled Average Bioequivalence

### RSABE (FDA Approach)

For highly variable drugs (CV > 30%):

```julia
# RSABE analysis
rsabe = rsabe_analysis(
    test_values,
    ref_values,
    design = FullReplicate2x4(),
    regulatory = FDAGuidance()
)

println("Within-subject CV: $(rsabe.cv_within)%")
println("Scaling applied: $(rsabe.scaling_applied)")
println("Scaled criterion: $(rsabe.scaled_criterion)")
println("Upper bound: $(rsabe.upper_bound)")
println("RSABE conclusion: $(rsabe.be_met)")
```

#### FDA RSABE Criterion

For CV > 30%:
- Scaled upper bound: $\sqrt{(\ln GMR)^2 + \theta \cdot s_{WR}^2} \leq \theta \cdot \sigma_{W0}$
- Where $\sigma_{W0} = 0.25$ (regulatory constant)
- $\theta = (\ln 1.25)^2 / \sigma_{W0}^2$

```julia
# FDA scaling parameters
const FDA_SIGMA_W0 = 0.25
const FDA_THETA = (log(1.25))^2 / FDA_SIGMA_W0^2
```

### ABEL (EMA Approach)

Average Bioequivalence with Expanding Limits:

```julia
# ABEL analysis
abel = abel_analysis(
    test_values,
    ref_values,
    design = FullReplicate2x4(),
    regulatory = EMAGuidance()
)

println("Reference CV: $(abel.cv_reference)%")
println("Widened limits: [$(abel.lower_limit)%, $(abel.upper_limit)%]")
println("GMR constraint (80-125%): $(abel.gmr_constraint_met)")
println("ABEL conclusion: $(abel.be_met)")
```

#### EMA ABEL Widening

For CV > 30%:
- Lower limit: $\exp(-k \cdot s_{WR})$
- Upper limit: $\exp(k \cdot s_{WR})$
- Maximum widening: 69.84% - 143.19%
- GMR must remain within 80-125%

```julia
# EMA widening parameters
const EMA_K = log(1.25) / 0.25  # Regulatory constant
const EMA_MAX_LOWER = 0.6984    # Maximum widened lower
const EMA_MAX_UPPER = 1.4319    # Maximum widened upper
```

---

## Replicate Study Designs

### Partial Replicate 3x3

```julia
design = PartialReplicate3x3()

# Sequences: TRR, RTR, RRT
# Reference replicated, Test single
# Used when Test formulation is limited

result = replicate_be_analysis(data, design=design)
```

### Full Replicate 2x4

```julia
design = FullReplicate2x4()

# Sequences: TRTR, RTRT
# Both formulations replicated
# Gold standard for HVD

result = replicate_be_analysis(data, design=design)
println("Subject-by-formulation variance: $(result.var_sbf)")
```

### Full Replicate 2x3

```julia
design = FullReplicate2x3()

# Sequences: TRT, RTR
# Shorter than 2x4 but still allows scaling
```

---

## Sample Size Calculation

### Standard BE Study

```julia
# Calculate required sample size
n = be_sample_size(
    cv = 0.25,          # Expected CV (25%)
    gmr = 0.95,         # Expected GMR
    power = 0.80,       # Target power
    alpha = 0.05,       # Significance level
    design = Crossover2x2()
)

println("Required subjects: $n per sequence")
println("Total subjects: $(2 * n)")
```

### Power Calculation

```julia
# Calculate power for given sample size
power = be_power(
    n = 24,
    cv = 0.25,
    gmr = 0.95,
    design = Crossover2x2()
)

println("Expected power: $(round(power * 100, digits=1))%")
```

### Sample Size Table

| CV | GMR=0.95 | GMR=0.90 | GMR=1.00 |
|----|----------|----------|----------|
| 15% | 10 | 14 | 8 |
| 20% | 16 | 24 | 12 |
| 25% | 24 | 36 | 18 |
| 30% | 34 | 50 | 26 |
| 35% | 46 | 68 | 36 |

---

## Regulatory Guidance

### FDA Guidance

```julia
config = BEConfig(
    regulatory = FDAGuidance(),
    acceptance_lower = 0.80,
    acceptance_upper = 1.25,
    alpha = 0.05
)

# FDA requires:
# - Fasted and fed studies (where applicable)
# - AUC0-t, AUC0-inf, Cmax
# - Log transformation
```

### EMA Guidance

```julia
config = BEConfig(
    regulatory = EMAGuidance(),
    acceptance_lower = 0.80,
    acceptance_upper = 1.25,
    alpha = 0.05,
    tmax_analysis = true  # EMA includes Tmax
)

# EMA requires:
# - Usually fasted only
# - AUC0-t, Cmax (and Tmax as supportive)
# - Widened limits for HVD Cmax
```

### Health Canada

```julia
config = BEConfig(
    regulatory = HealthCanadaGuidance(),
    acceptance_lower = 0.80,
    acceptance_upper = 1.25
)
```

---

## Example: Complete BE Analysis

```julia
using NeoPKPD, DataFrames

# Crossover study data
data = DataFrame(
    subject = repeat(1:24, inner=2),
    period = repeat([1, 2], 24),
    sequence = repeat(["TR", "RT"], inner=24),
    formulation = vcat(
        repeat(["T", "R"], 12),  # TR sequence
        repeat(["R", "T"], 12)   # RT sequence
    ),
    auc = [45.2, 48.1, 52.1, 55.3, ...],  # AUC values
    cmax = [8.5, 9.1, 9.2, 9.8, ...]       # Cmax values
)

# Configure analysis
config = BEConfig(
    regulatory = FDAGuidance(),
    design = Crossover2x2(),
    alpha = 0.05
)

# Analyze AUC
auc_result = analyze_be(data, :auc, config)
println("=== AUC Bioequivalence ===")
println("GMR: $(round(auc_result.gmr * 100, digits=2))%")
println("90% CI: [$(round(auc_result.lower * 100, digits=2))%, $(round(auc_result.upper * 100, digits=2))%]")
println("Within-subject CV: $(round(auc_result.cv_within * 100, digits=1))%")
println("BE conclusion: $(auc_result.be_met ? "PASS" : "FAIL")")

# Analyze Cmax
cmax_result = analyze_be(data, :cmax, config)
println("\n=== Cmax Bioequivalence ===")
println("GMR: $(round(cmax_result.gmr * 100, digits=2))%")
println("90% CI: [$(round(cmax_result.lower * 100, digits=2))%, $(round(cmax_result.upper * 100, digits=2))%]")
println("Within-subject CV: $(round(cmax_result.cv_within * 100, digits=1))%")
println("BE conclusion: $(cmax_result.be_met ? "PASS" : "FAIL")")

# Overall conclusion
overall_be = auc_result.be_met && cmax_result.be_met
println("\n=== Overall Conclusion ===")
println("Bioequivalence: $(overall_be ? "ESTABLISHED" : "NOT ESTABLISHED")")

# Generate regulatory report
report = be_regulatory_report(
    auc_result,
    cmax_result,
    regulatory = FDAGuidance()
)
println(report)
```

---

## Highly Variable Drug Example

```julia
# HVD with CV > 30%
data_hvd = load_hvd_study_data()

# Check variability
cv = within_subject_cv(data_hvd, :auc, :subject_id)
println("Within-subject CV: $(round(cv * 100, digits=1))%")

if cv > 0.30
    println("HVD criteria met - applying reference scaling")

    # Use replicate design analysis
    result = rsabe_analysis(
        data_hvd,
        design = FullReplicate2x4(),
        regulatory = FDAGuidance()
    )

    println("\n=== RSABE Analysis ===")
    println("GMR: $(round(result.gmr * 100, digits=2))%")
    println("Reference CV: $(round(result.cv_reference * 100, digits=1))%")
    println("Scaling applied: $(result.scaling_applied)")
    println("Scaled criterion met: $(result.scaled_criterion_met)")
    println("Point estimate constraint: $(result.point_estimate_met)")
    println("RSABE conclusion: $(result.be_met ? "PASS" : "FAIL")")
else
    # Standard ABE
    result = analyze_be(data_hvd, :auc, BEConfig())
end
```

---

## Formulas Summary

| Parameter | Formula |
|-----------|---------|
| GMR | $\exp(\bar{X}_T - \bar{X}_R)$ where X = ln(value) |
| 90% CI | $GMR \cdot \exp(\pm t_{0.95,df} \cdot SE)$ |
| Within-subject CV | $\sqrt{\exp(MSE) - 1}$ |
| RSABE criterion | $\sqrt{(\ln GMR)^2 + \theta \cdot s_{WR}^2}$ |
| ABEL limits | $\exp(\pm k \cdot s_{WR})$ |

---

## See Also

- [Exposure Metrics](exposure-metrics.md) - AUC and Cmax calculation
- [Population NCA](population-nca.md) - Multi-subject analysis
- [Terminal Phase](terminal-phase.md) - Lambda_z and half-life

