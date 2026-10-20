# Population NCA

Comprehensive documentation for multi-subject non-compartmental analysis including summary statistics and stratification.

---

## Overview

Population NCA applies non-compartmental analysis to multiple subjects simultaneously, generating individual results and summary statistics across the population.

---

## Quick Start

```julia
using NeoPKPD, DataFrames

# Multi-subject PK data
data = DataFrame(
    subject_id = repeat(1:12, inner=8),
    time = repeat([0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0], 12),
    conc = [...],  # Concentration values
    dose = repeat([100.0], 96)
)

# Run population NCA
pop_result = run_population_nca(data; dose_col=:dose)

# View summary
summary = summarize_population_nca(pop_result)
println(summary)
```

---

## Data Format

### Required Columns

```julia
data = DataFrame(
    subject_id = [1, 1, 1, 2, 2, 2, ...],  # Subject identifier
    time = [0.0, 1.0, 4.0, 0.0, 1.0, 4.0, ...],  # Time points
    conc = [0.0, 5.2, 2.1, 0.0, 4.8, 1.9, ...],  # Concentrations
    dose = [100.0, 100.0, 100.0, 150.0, ...]      # Dose per subject
)
```

### Optional Columns

```julia
# With covariates and grouping
data = DataFrame(
    subject_id = [...],
    time = [...],
    conc = [...],
    dose = [...],
    weight = [70.0, 85.0, ...],        # Body weight
    sex = ["M", "F", ...],              # Sex
    formulation = ["Test", "Ref", ...], # Treatment group
    period = [1, 1, 2, 2, ...]          # Study period
)
```

---

## Running Population NCA

### Basic Usage

```julia
pop_result = run_population_nca(
    data;
    subject_col = :subject_id,
    time_col = :time,
    conc_col = :conc,
    dose_col = :dose
)

# Access individual results
for (subj, result) in pop_result.individual_results
    println("Subject $subj: AUC=$(result.auc_0_inf), Cmax=$(result.cmax)")
end
```

### With Configuration

```julia
config = NCAConfig(
    method = LinLogMixedMethod(),
    lambda_z_min_points = 3,
    lambda_z_r2_threshold = 0.9,
    lloq = 0.05,
    blq_handling = BLQMissing()
)

pop_result = run_population_nca(
    data;
    config = config,
    route = :extravascular
)
```

### Multiple Dose Analysis

```julia
pop_result = run_population_nca(
    data;
    dosing_type = :multiple,
    tau = 24.0,  # Dosing interval
    dose_col = :dose
)
```

---

## Summary Statistics

### Generate Summary

```julia
summary = summarize_population_nca(pop_result)

# Access statistics
println("=== AUC0-inf Summary ===")
println("N: $(summary.auc_0_inf.n)")
println("Mean: $(summary.auc_0_inf.mean)")
println("SD: $(summary.auc_0_inf.sd)")
println("CV%: $(summary.auc_0_inf.cv_pct)")
println("Median: $(summary.auc_0_inf.median)")
println("Min: $(summary.auc_0_inf.min)")
println("Max: $(summary.auc_0_inf.max)")
println("Geometric Mean: $(summary.auc_0_inf.geomean)")
println("Geometric CV%: $(summary.auc_0_inf.geocv_pct)")
```

### Summary Table Format

```julia
# Get formatted summary table
table = population_nca_summary_table(pop_result)

# Output format:
# Parameter     N    Mean    SD     CV%    Median   Min    Max    GeoMean  GeoCV%
# ---------------------------------------------------------------------------
# Cmax          12   5.23    1.12   21.4   5.15     3.42   7.81   5.12     22.1
# Tmax          12   1.25    0.35   28.0   1.00     0.50   2.00   -        -
# AUC0-t        12   45.2    8.7    19.3   44.1     31.2   62.5   44.5     19.8
# AUC0-inf      12   48.1    9.2    19.1   47.2     33.5   68.2   47.4     19.5
# t1/2          12   4.52    0.85   18.8   4.35     3.12   6.21   4.45     19.2
# CL/F          12   2.12    0.42   19.8   2.08     1.47   3.01   2.08     20.1
```

### Custom Statistics

```julia
# Select specific parameters
summary = summarize_population_nca(
    pop_result;
    parameters = [:cmax, :auc_0_inf, :t_half, :cl_f]
)

# Include percentiles
summary = summarize_population_nca(
    pop_result;
    percentiles = [5, 25, 50, 75, 95]
)
```

---

## Stratified Analysis

### By Single Variable

```julia
# Stratify by formulation
stratified = stratified_population_nca(
    data,
    strat_col = :formulation
)

for (stratum, result) in stratified
    println("=== $stratum ===")
    summary = summarize_population_nca(result)
    println("AUC0-inf mean: $(summary.auc_0_inf.mean)")
end
```

### By Multiple Variables

```julia
# Stratify by formulation and fed/fasted
stratified = stratified_population_nca(
    data,
    strat_cols = [:formulation, :fed_state]
)

# Results grouped by combination
# e.g., ("Test", "Fed"), ("Test", "Fasted"), ("Ref", "Fed"), ("Ref", "Fasted")
```

### Comparison Between Strata

```julia
# Compare Test vs Reference
test_summary = summarize_population_nca(stratified["Test"])
ref_summary = summarize_population_nca(stratified["Ref"])

# Geometric mean ratio
gmr_auc = test_summary.auc_0_inf.geomean / ref_summary.auc_0_inf.geomean
gmr_cmax = test_summary.cmax.geomean / ref_summary.cmax.geomean

println("AUC GMR: $(round(gmr_auc * 100, digits=2))%")
println("Cmax GMR: $(round(gmr_cmax * 100, digits=2))%")
```

---

## Individual Results Access

### Accessing All Results

```julia
pop_result = run_population_nca(data)

# Iterate through individuals
for (subject_id, result) in pop_result.individual_results
    println("Subject $subject_id:")
    println("  Cmax: $(result.cmax) at Tmax: $(result.tmax)")
    println("  AUC0-t: $(result.auc_0_t)")
    println("  AUC0-inf: $(result.auc_0_inf)")
    println("  t1/2: $(result.t_half)")
    println("  Lambda_z R²: $(result.lambda_z_result.r_squared)")
    println("")
end
```

### Export to DataFrame

```julia
# Convert results to DataFrame
results_df = population_nca_to_dataframe(pop_result)

# Columns: subject_id, cmax, tmax, auc_0_t, auc_0_inf, t_half, cl_f, vz_f, ...
```

### Filter by Quality

```julia
# Only include subjects with good lambda_z fit
quality_results = filter(pop_result.individual_results) do (_, result)
    !isnothing(result.lambda_z_result) &&
    result.lambda_z_result.r_squared >= 0.9 &&
    result.auc_extra_pct <= 20.0
end

println("Subjects meeting quality criteria: $(length(quality_results))/$(length(pop_result.individual_results))")
```

---

## Dose Normalization

### Normalize by Dose

```julia
# Dose-normalized population NCA
pop_result = run_population_nca(
    data;
    dose_normalize = true,
    dose_col = :dose
)

# Access normalized values
for (subj, result) in pop_result.individual_results
    println("Subject $subj: AUC/D = $(result.auc_dn), Cmax/D = $(result.cmax_dn)")
end
```

### Weight-Normalized Dose

```julia
# Dose per kg body weight
data.dose_per_kg = data.dose ./ data.weight

pop_result = run_population_nca(
    data;
    dose_col = :dose_per_kg,
    dose_normalize = true
)
```

---

## Handling Missing Data

### Subjects with Insufficient Data

```julia
pop_result = run_population_nca(data; config=config)

# Check for failed subjects
if !isempty(pop_result.failed_subjects)
    println("Subjects with NCA failures:")
    for (subj, reason) in pop_result.failed_subjects
        println("  Subject $subj: $reason")
    end
end
```

### Partial Results

```julia
# Some parameters may be missing for individual subjects
for (subj, result) in pop_result.individual_results
    if isnan(result.auc_0_inf)
        println("Subject $subj: AUC0-inf not calculable (lambda_z failed)")
    end
end
```

---

## Covariate Analysis

### Summary by Covariate

```julia
# Summarize by weight category
data.wt_cat = ifelse.(data.weight .< 70, "Low", "High")

stratified = stratified_population_nca(data, strat_col=:wt_cat)

for (cat, result) in stratified
    summary = summarize_population_nca(result)
    println("$cat weight: CL/F = $(summary.cl_f.mean) ± $(summary.cl_f.sd)")
end
```

### Correlation Analysis

```julia
# Extract PK parameters and covariates
results_df = population_nca_to_dataframe(pop_result)

# Merge with covariate data
merged = leftjoin(results_df, unique(data[:, [:subject_id, :weight, :age]]), on=:subject_id)

# Calculate correlations
using Statistics
r_cl_wt = cor(merged.cl_f, merged.weight)
r_vz_wt = cor(merged.vz_f, merged.weight)

println("CL/F vs Weight correlation: $(round(r_cl_wt, digits=3))")
println("Vz/F vs Weight correlation: $(round(r_vz_wt, digits=3))")
```

---

## Example: Complete Population NCA

```julia
using NeoPKPD, DataFrames

# Load multi-subject PK data
data = DataFrame(
    subject_id = repeat(1:24, inner=10),
    time = repeat([0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0], 24),
    conc = [...],  # 240 concentration values
    dose = repeat([500.0], 240),
    weight = repeat([65.0, 72.0, 80.0, ...], inner=10),  # 24 weights
    formulation = repeat(["Test", "Reference"], inner=120)
)

# Configure NCA
config = NCAConfig(
    method = LinLogMixedMethod(),
    lambda_z_min_points = 3,
    lambda_z_r2_threshold = 0.9,
    extrapolation_max_pct = 20.0,
    lloq = 0.1
)

# Run population NCA
pop_result = run_population_nca(
    data;
    config = config,
    route = :extravascular,
    dose_col = :dose
)

# Overall summary
println("=== Population NCA Summary ===")
println("Total subjects: $(length(pop_result.individual_results))")
println("Failed subjects: $(length(pop_result.failed_subjects))")

summary = summarize_population_nca(pop_result)
println("\nParameter    N    Mean±SD         CV%    GeoMean  GeoCV%")
println("-" ^ 60)
for param in [:cmax, :tmax, :auc_0_t, :auc_0_inf, :t_half, :cl_f, :vz_f]
    s = getfield(summary, param)
    if param == :tmax
        @printf("%-10s  %d   %.2f±%.2f    %.1f%%   -        -\n",
            string(param), s.n, s.mean, s.sd, s.cv_pct)
    else
        @printf("%-10s  %d   %.2f±%.2f    %.1f%%   %.2f     %.1f%%\n",
            string(param), s.n, s.mean, s.sd, s.cv_pct, s.geomean, s.geocv_pct)
    end
end

# Stratified by formulation
println("\n=== By Formulation ===")
stratified = stratified_population_nca(data, strat_col=:formulation)

for form in ["Test", "Reference"]
    s = summarize_population_nca(stratified[form])
    println("\n$form:")
    println("  AUC0-inf: $(round(s.auc_0_inf.geomean, digits=2)) (GeoCV $(round(s.auc_0_inf.geocv_pct, digits=1))%)")
    println("  Cmax:     $(round(s.cmax.geomean, digits=2)) (GeoCV $(round(s.cmax.geocv_pct, digits=1))%)")
end

# GMR calculation
test_s = summarize_population_nca(stratified["Test"])
ref_s = summarize_population_nca(stratified["Reference"])
gmr_auc = test_s.auc_0_inf.geomean / ref_s.auc_0_inf.geomean
gmr_cmax = test_s.cmax.geomean / ref_s.cmax.geomean

println("\n=== Geometric Mean Ratios (Test/Reference) ===")
println("AUC GMR: $(round(gmr_auc * 100, digits=2))%")
println("Cmax GMR: $(round(gmr_cmax * 100, digits=2))%")

# Quality assessment
println("\n=== Quality Assessment ===")
n_good_lz = count(r -> r.lambda_z_result.r_squared >= 0.9, values(pop_result.individual_results))
n_low_extrap = count(r -> r.auc_extra_pct <= 20.0, values(pop_result.individual_results))
total = length(pop_result.individual_results)

println("Lambda_z R² ≥ 0.9: $n_good_lz/$total ($(round(100*n_good_lz/total, digits=1))%)")
println("AUC extrap ≤ 20%: $n_low_extrap/$total ($(round(100*n_low_extrap/total, digits=1))%)")

# Export for further analysis
results_df = population_nca_to_dataframe(pop_result)
println("\nResults exported to DataFrame with $(nrow(results_df)) subjects")
```

---

## Output Formats

### Summary Statistics Structure

```julia
struct NCAParameterSummary
    n::Int              # Number of subjects
    mean::Float64       # Arithmetic mean
    sd::Float64         # Standard deviation
    cv_pct::Float64     # CV% (100 * SD/mean)
    median::Float64     # Median
    min::Float64        # Minimum
    max::Float64        # Maximum
    geomean::Float64    # Geometric mean
    geocv_pct::Float64  # Geometric CV%
    percentiles::Dict{Int,Float64}  # Requested percentiles
end
```

### Export Formats

```julia
# CSV export
population_nca_to_csv(pop_result, "nca_results.csv")

# JSON export
population_nca_to_json(pop_result, "nca_results.json")

# DataFrame (in-memory)
df = population_nca_to_dataframe(pop_result)
```

---

## See Also

- [Exposure Metrics](exposure-metrics.md) - Individual NCA metrics
- [Terminal Phase](terminal-phase.md) - Lambda_z estimation
- [Bioequivalence](bioequivalence.md) - BE analysis from population NCA
- [Multiple Dose](multiple-dose.md) - Steady-state population analysis

