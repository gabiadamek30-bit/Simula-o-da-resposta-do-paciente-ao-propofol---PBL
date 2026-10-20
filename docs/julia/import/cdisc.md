# CDISC Data Import

Comprehensive guide for importing CDISC/SDTM and ADaM formatted data into NeoPKPD.

---

## Overview

NeoPKPD supports importing data from CDISC (Clinical Data Interchange Standards Consortium) formats, including SDTM (Study Data Tabulation Model) and ADaM (Analysis Data Model) datasets.

```julia
using NeoPKPD

# Import CDISC data
data = import_cdisc(
    pc_path = "pc.csv",
    ex_path = "ex.csv",
    dm_path = "dm.csv"
)
```

---

## Supported Domains

### SDTM Domains

| Domain | Full Name | Description | Support |
|--------|-----------|-------------|---------|
| PC | Pharmacokinetic Concentrations | Drug concentration measurements | ✅ Full |
| PP | Pharmacokinetic Parameters | Derived PK parameters | ✅ Full |
| EX | Exposure | Dosing information | ✅ Full |
| DM | Demographics | Subject characteristics | ✅ Full |
| VS | Vital Signs | Weight, height, etc. | ⚠️ Partial |
| LB | Laboratory | Lab test results | ⚠️ Partial |

### ADaM Datasets

| Dataset | Description | Support |
|---------|-------------|---------|
| ADPC | Analysis PK Concentrations | ✅ Full |
| ADPP | Analysis PK Parameters | ✅ Full |
| ADSL | Subject-Level Analysis | ✅ Full |
| ADEX | Analysis Exposure | ⚠️ Partial |

---

## PC Domain (Pharmacokinetic Concentrations)

### Required Variables

| Variable | Label | Type | Description |
|----------|-------|------|-------------|
| STUDYID | Study Identifier | Char | Unique study ID |
| USUBJID | Unique Subject ID | Char | Unique subject identifier |
| PCTESTCD | Test Code | Char | Short test name (e.g., "DRUG") |
| PCTEST | Test Name | Char | Full test name |
| PCORRES | Original Result | Char | Result as collected |
| PCSTRESN | Numeric Result | Num | Standardized numeric result |
| PCSTRESU | Unit | Char | Units (ng/mL, µg/L, etc.) |
| PCELTM | Elapsed Time | Char | Time from reference (ISO 8601) |

### Optional Variables

| Variable | Label | Description |
|----------|-------|-------------|
| PCLLOQ | Lower LOQ | Lower limit of quantification |
| PCSTAT | Status | Completion status (null or "NOT DONE") |
| PCREASND | Reason Not Done | Reason if PCSTAT="NOT DONE" |
| PCBLFL | Baseline Flag | "Y" for baseline records |
| PCSPEC | Specimen | Specimen type (PLASMA, SERUM, etc.) |
| PCMETHOD | Method | Analytical method |
| VISITNUM | Visit Number | Planned visit number |
| VISIT | Visit Name | Visit description |
| PCDTC | Collection Datetime | ISO 8601 datetime |
| PCDY | Study Day | Study day of collection |

### Example PC Data

```csv
STUDYID,USUBJID,PCTESTCD,PCTEST,PCSTRESN,PCSTRESU,PCELTM,PCDTC,PCSPEC,PCLLOQ
STUDY01,SUBJ001,DRUGA,Drug A Concentration,0.0,ng/mL,PT0H,2024-01-15T08:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A Concentration,125.3,ng/mL,PT1H,2024-01-15T09:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A Concentration,89.7,ng/mL,PT2H,2024-01-15T10:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A Concentration,52.1,ng/mL,PT4H,2024-01-15T12:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A Concentration,28.4,ng/mL,PT8H,2024-01-15T16:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A Concentration,12.1,ng/mL,PT12H,2024-01-15T20:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A Concentration,3.2,ng/mL,PT24H,2024-01-16T08:00,PLASMA,0.5
```

---

## EX Domain (Exposure)

### Required Variables

| Variable | Label | Type | Description |
|----------|-------|------|-------------|
| STUDYID | Study Identifier | Char | Unique study ID |
| USUBJID | Unique Subject ID | Char | Unique subject identifier |
| EXDOSE | Dose | Num | Administered dose amount |
| EXDOSU | Dose Units | Char | Units (mg, µg, etc.) |
| EXROUTE | Route | Char | Route of administration |
| EXSTDTC | Start Datetime | Char | Dose start (ISO 8601) |

### Optional Variables

| Variable | Label | Description |
|----------|-------|-------------|
| EXENDTC | End Datetime | Dose end (for infusions) |
| EXDOSFRM | Dose Form | TABLET, CAPSULE, SOLUTION, etc. |
| EXDOSFRQ | Dosing Frequency | QD, BID, Q12H, etc. |
| EXDUR | Duration | Infusion duration (ISO 8601) |
| EXTRT | Treatment | Treatment name |
| EXSEQ | Sequence | Sequence number |
| VISITNUM | Visit Number | Planned visit number |
| EXDY | Study Day | Study day of dose |

### Example EX Data

```csv
STUDYID,USUBJID,EXTRT,EXDOSE,EXDOSU,EXROUTE,EXDOSFRM,EXSTDTC,EXDY
STUDY01,SUBJ001,Drug A,100,mg,ORAL,TABLET,2024-01-15T08:00,1
STUDY01,SUBJ001,Drug A,100,mg,ORAL,TABLET,2024-01-16T08:00,2
STUDY01,SUBJ002,Drug A,200,mg,ORAL,TABLET,2024-01-15T08:00,1
STUDY01,SUBJ002,Drug A,200,mg,ORAL,TABLET,2024-01-16T08:00,2
```

---

## DM Domain (Demographics)

### Required Variables

| Variable | Label | Type | Description |
|----------|-------|------|-------------|
| STUDYID | Study Identifier | Char | Unique study ID |
| USUBJID | Unique Subject ID | Char | Unique subject identifier |
| SUBJID | Subject ID | Char | Subject identifier within study |
| RFSTDTC | Reference Start | Char | Reference start date (ISO 8601) |
| RFENDTC | Reference End | Char | Reference end date |

### Optional Variables

| Variable | Label | Description |
|----------|-------|-------------|
| AGE | Age | Age at screening |
| AGEU | Age Units | YEARS, MONTHS, etc. |
| SEX | Sex | M, F, U |
| RACE | Race | Racial designation |
| ETHNIC | Ethnicity | Ethnic group |
| ARMCD | Arm Code | Treatment arm code |
| ARM | Arm Description | Treatment arm description |
| COUNTRY | Country | Country of participation |
| SITEID | Site Identifier | Study site ID |
| DMDTC | Collection Date | Demographics collection date |

### Example DM Data

```csv
STUDYID,USUBJID,SUBJID,AGE,AGEU,SEX,RACE,ETHNIC,ARM,ARMCD,RFSTDTC
STUDY01,SUBJ001,001,45,YEARS,M,WHITE,NOT HISPANIC OR LATINO,Treatment A,TRT_A,2024-01-15
STUDY01,SUBJ002,002,38,YEARS,F,ASIAN,NOT HISPANIC OR LATINO,Treatment A,TRT_A,2024-01-15
STUDY01,SUBJ003,003,52,YEARS,M,BLACK OR AFRICAN AMERICAN,NOT HISPANIC OR LATINO,Treatment B,TRT_B,2024-01-15
```

---

## PP Domain (PK Parameters)

### Required Variables

| Variable | Label | Type | Description |
|----------|-------|------|-------------|
| STUDYID | Study Identifier | Char | Unique study ID |
| USUBJID | Unique Subject ID | Char | Unique subject identifier |
| PPTESTCD | Parameter Code | Char | Short parameter name |
| PPTEST | Parameter Name | Char | Full parameter name |
| PPSTRESN | Numeric Result | Num | Standardized result |
| PPSTRESU | Units | Char | Result units |

### Common PK Parameter Codes

| PPTESTCD | Description | Units |
|----------|-------------|-------|
| CMAX | Maximum Concentration | ng/mL |
| TMAX | Time to Maximum | h |
| AUCLST | AUC to Last | ng·h/mL |
| AUCIFO | AUC Extrapolated | ng·h/mL |
| LAMZ | Lambda z | 1/h |
| LAMZHL | Terminal Half-life | h |
| CLFO | Clearance | L/h |
| VZFO | Volume of Distribution | L |

### Example PP Data

```csv
STUDYID,USUBJID,PPTESTCD,PPTEST,PPSTRESN,PPSTRESU,PPCAT
STUDY01,SUBJ001,CMAX,Maximum Concentration,125.3,ng/mL,PK PARAMETERS
STUDY01,SUBJ001,TMAX,Time to Maximum,1.0,h,PK PARAMETERS
STUDY01,SUBJ001,AUCLST,AUC to Last Observation,845.2,ng.h/mL,PK PARAMETERS
STUDY01,SUBJ001,AUCIFO,AUC Extrapolated to Infinity,892.4,ng.h/mL,PK PARAMETERS
STUDY01,SUBJ001,LAMZHL,Terminal Half-life,8.5,h,PK PARAMETERS
```

---

## Import Functions

### Basic Import

```julia
using NeoPKPD

# Import from CSV files
data = import_cdisc(
    pc_path = "pc.csv",
    ex_path = "ex.csv",
    dm_path = "dm.csv"
)

# Access imported data
println("Subjects: $(length(data.subjects))")
println("Total observations: $(n_observations(data))")
```

### Individual Domain Import

```julia
# Import PC domain only
pc_data = read_cdisc_pc("pc.csv")

# Import EX domain only
ex_data = read_cdisc_ex("ex.csv")

# Import DM domain only
dm_data = read_cdisc_dm("dm.csv")

# Combine manually
data = cdisc_to_observed_data(pc_data, ex_data, dm_data)
```

### Import from SAS Transport (XPT)

```julia
# Import XPT files
data = import_cdisc(
    pc_path = "pc.xpt",
    ex_path = "ex.xpt",
    dm_path = "dm.xpt",
    format = :xpt
)
```

---

## Data Structures

### ObservedData

```julia
struct ObservedData
    subjects::Vector{SubjectData}  # All subjects
    study_id::String               # Study identifier
    analyte::String                # Drug/analyte name
    units::String                  # Concentration units
    time_units::String             # Time units
end
```

### SubjectData

```julia
struct SubjectData
    subject_id::String             # Unique subject ID
    times::Vector{Float64}         # Time points (h from first dose)
    observations::Vector{Float64}  # Concentrations
    doses::Vector{DoseEvent}       # Dose events
    covariates::Dict{Symbol,Any}   # Age, weight, sex, etc.
    lloq::Float64                  # Lower limit of quantification
    blq_flags::Vector{Bool}        # Below LOQ indicators
end
```

### Accessing Subject Data

```julia
data = import_cdisc(pc_path="pc.csv", ex_path="ex.csv", dm_path="dm.csv")

# Iterate through subjects
for subj in data.subjects
    println("Subject: $(subj.subject_id)")
    println("  Time points: $(length(subj.times))")
    println("  Doses: $(length(subj.doses))")
    println("  Age: $(subj.covariates[:age])")
    println("  Weight: $(subj.covariates[:weight])")
end

# Get specific subject
subj = get_subject(data, "SUBJ001")
```

---

## Time Handling

### Elapsed Time Parsing

```julia
# PCELTM format (ISO 8601 duration)
# PT0H = 0 hours (pre-dose)
# PT1H = 1 hour post-dose
# PT2H30M = 2.5 hours post-dose
# P1DT2H = 1 day + 2 hours = 26 hours

# NeoPKPD converts to numeric hours from first dose
```

### Reference Time

```julia
# Time is calculated relative to:
# 1. First dose in EX domain (EXSTDTC)
# 2. Or reference start date in DM (RFSTDTC)

# Custom reference time
data = import_cdisc(
    pc_path = "pc.csv",
    ex_path = "ex.csv",
    dm_path = "dm.csv",
    reference_time = :first_dose  # or :rfstdtc
)
```

---

## BLQ Handling

### Detection

```julia
# BLQ detected from:
# 1. PCSTRESN == 0 or missing
# 2. PCSTAT == "NOT DONE"
# 3. PCORRES contains "<", "BLQ", "BLOQ"
# 4. PCSTRESN < PCLLOQ

data = import_cdisc(pc_path="pc.csv", ex_path="ex.csv", dm_path="dm.csv")

# Access BLQ flags
for subj in data.subjects
    n_blq = sum(subj.blq_flags)
    println("Subject $(subj.subject_id): $(n_blq) BLQ samples")
end
```

### BLQ Handling Options

```julia
# Default: BLQ values set to 0
data = import_cdisc(pc_path="pc.csv", ex_path="ex.csv", dm_path="dm.csv",
    blq_handling = :zero
)

# Set BLQ to LLOQ/2
data = import_cdisc(pc_path="pc.csv", ex_path="ex.csv", dm_path="dm.csv",
    blq_handling = :lloq_half
)

# Keep as missing (NaN)
data = import_cdisc(pc_path="pc.csv", ex_path="ex.csv", dm_path="dm.csv",
    blq_handling = :missing
)
```

---

## Validation

### Automatic Validation

```julia
# Validation checks performed automatically:
# 1. Required columns present
# 2. Data types correct
# 3. Subject IDs consistent across domains
# 4. Dose information available for PK subjects
# 5. No negative concentrations
# 6. Time sequence valid

data = import_cdisc(pc_path="pc.csv", ex_path="ex.csv", dm_path="dm.csv")

# Check validation results
if !isempty(data.validation_warnings)
    println("Validation warnings:")
    for w in data.validation_warnings
        println("  ⚠️ $w")
    end
end
```

### Manual Validation

```julia
# Validate before import
validation = validate_cdisc_dataset(
    pc_path = "pc.csv",
    ex_path = "ex.csv",
    dm_path = "dm.csv"
)

if validation.valid
    println("✓ Dataset passes validation")
    data = import_cdisc(pc_path="pc.csv", ex_path="ex.csv", dm_path="dm.csv")
else
    println("❌ Validation failed:")
    for error in validation.errors
        println("  - $error")
    end
end
```

---

## Dose Reconstruction

### From EX Domain

```julia
# Doses automatically extracted from EX domain
# Route mapping:
# ORAL → :oral
# INTRAVENOUS → :iv_bolus
# INTRAVENOUS INFUSION → :iv_infusion
# SUBCUTANEOUS → :subcutaneous
# INTRAMUSCULAR → :intramuscular

data = import_cdisc(pc_path="pc.csv", ex_path="ex.csv", dm_path="dm.csv")

for subj in data.subjects
    for dose in subj.doses
        println("Dose: $(dose.amount) $(dose.units) at $(dose.time)h via $(dose.route)")
    end
end
```

### Multiple Dose Handling

```julia
# Multiple doses per subject are supported
# Dose times relative to first dose

# For steady-state data:
data = import_cdisc(
    pc_path = "pc_ss.csv",
    ex_path = "ex_ss.csv",
    dm_path = "dm.csv",
    dosing_type = :steady_state,
    tau = 24.0  # Dosing interval
)
```

---

## Covariate Extraction

### From DM Domain

```julia
# Standard covariates extracted:
# - age (AGE)
# - sex (SEX → :M, :F)
# - weight (if in DM or VS)
# - race (RACE)
# - ethnicity (ETHNIC)

data = import_cdisc(pc_path="pc.csv", ex_path="ex.csv", dm_path="dm.csv")

# Access covariates
for subj in data.subjects
    println("Subject $(subj.subject_id):")
    println("  Age: $(subj.covariates[:age])")
    println("  Sex: $(subj.covariates[:sex])")
    println("  Weight: $(get(subj.covariates, :weight, "N/A"))")
end
```

### Additional Covariates from VS/LB

```julia
# Include vital signs
data = import_cdisc(
    pc_path = "pc.csv",
    ex_path = "ex.csv",
    dm_path = "dm.csv",
    vs_path = "vs.csv",  # Include vital signs
    covariates = [:weight, :height, :bmi]
)
```

---

## Complete Example

### Sample Dataset Files

**pc.csv**:
```csv
STUDYID,USUBJID,PCTESTCD,PCTEST,PCSTRESN,PCSTRESU,PCELTM,PCDTC,PCSPEC,PCLLOQ
STUDY01,SUBJ001,DRUGA,Drug A,0.0,ng/mL,PT0H,2024-01-15T08:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A,125.3,ng/mL,PT1H,2024-01-15T09:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A,89.7,ng/mL,PT2H,2024-01-15T10:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A,52.1,ng/mL,PT4H,2024-01-15T12:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A,28.4,ng/mL,PT8H,2024-01-15T16:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A,12.1,ng/mL,PT12H,2024-01-15T20:00,PLASMA,0.5
STUDY01,SUBJ001,DRUGA,Drug A,3.2,ng/mL,PT24H,2024-01-16T08:00,PLASMA,0.5
STUDY01,SUBJ002,DRUGA,Drug A,0.0,ng/mL,PT0H,2024-01-15T08:00,PLASMA,0.5
STUDY01,SUBJ002,DRUGA,Drug A,198.5,ng/mL,PT1H,2024-01-15T09:00,PLASMA,0.5
STUDY01,SUBJ002,DRUGA,Drug A,142.3,ng/mL,PT2H,2024-01-15T10:00,PLASMA,0.5
```

**ex.csv**:
```csv
STUDYID,USUBJID,EXTRT,EXDOSE,EXDOSU,EXROUTE,EXDOSFRM,EXSTDTC,EXDY
STUDY01,SUBJ001,Drug A,100,mg,ORAL,TABLET,2024-01-15T08:00,1
STUDY01,SUBJ002,Drug A,200,mg,ORAL,TABLET,2024-01-15T08:00,1
```

**dm.csv**:
```csv
STUDYID,USUBJID,SUBJID,AGE,AGEU,SEX,RACE,ETHNIC,ARM,ARMCD,RFSTDTC
STUDY01,SUBJ001,001,45,YEARS,M,WHITE,NOT HISPANIC OR LATINO,100mg,TRT100,2024-01-15
STUDY01,SUBJ002,002,38,YEARS,F,ASIAN,NOT HISPANIC OR LATINO,200mg,TRT200,2024-01-15
```

### Julia Import Code

```julia
using NeoPKPD

# Import CDISC data
data = import_cdisc(
    pc_path = "pc.csv",
    ex_path = "ex.csv",
    dm_path = "dm.csv"
)

# Display import results
println("=" ^ 50)
println("CDISC Data Import Results")
println("=" ^ 50)

println("\n--- Study Information ---")
println("Study ID: $(data.study_id)")
println("Analyte: $(data.analyte)")
println("Units: $(data.units)")

println("\n--- Dataset Summary ---")
println("Total subjects: $(length(data.subjects))")
println("Total observations: $(n_observations(data))")

println("\n--- Subject Details ---")
for subj in data.subjects
    println("\nSubject $(subj.subject_id):")
    println("  Observations: $(length(subj.times))")
    println("  Time range: $(minimum(subj.times)) to $(maximum(subj.times)) h")
    println("  Doses: $(length(subj.doses))")

    if !isempty(subj.doses)
        dose = subj.doses[1]
        println("  First dose: $(dose.amount) mg at $(dose.time)h")
    end

    println("  Covariates:")
    println("    Age: $(subj.covariates[:age]) years")
    println("    Sex: $(subj.covariates[:sex])")

    n_blq = sum(subj.blq_flags)
    if n_blq > 0
        println("  BLQ samples: $n_blq")
    end
end

# Convert to DataFrame for analysis
println("\n--- Export to DataFrame ---")
df = to_dataframe(data)
println("DataFrame created with $(nrow(df)) rows")

# Use with NCA
println("\n--- NCA Analysis ---")
for subj in data.subjects
    result = run_nca(subj.times, subj.observations, subj.doses[1].amount)
    println("Subject $(subj.subject_id):")
    println("  Cmax: $(round(result.cmax, digits=2)) $(data.units)")
    println("  Tmax: $(result.tmax) h")
    println("  AUC0-inf: $(round(result.auc_0_inf, digits=2)) $(data.units)·h")
end
```

### Expected Output

```
==================================================
CDISC Data Import Results
==================================================

--- Study Information ---
Study ID: STUDY01
Analyte: Drug A
Units: ng/mL

--- Dataset Summary ---
Total subjects: 2
Total observations: 10

--- Subject Details ---

Subject SUBJ001:
  Observations: 7
  Time range: 0.0 to 24.0 h
  Doses: 1
  First dose: 100.0 mg at 0.0h
  Covariates:
    Age: 45 years
    Sex: M

Subject SUBJ002:
  Observations: 3
  Time range: 0.0 to 2.0 h
  Doses: 1
  First dose: 200.0 mg at 0.0h
  Covariates:
    Age: 38 years
    Sex: F

--- Export to DataFrame ---
DataFrame created with 10 rows

--- NCA Analysis ---
Subject SUBJ001:
  Cmax: 125.3 ng/mL
  Tmax: 1.0 h
  AUC0-inf: 523.45 ng/mL·h
Subject SUBJ002:
  Cmax: 198.5 ng/mL
  Tmax: 1.0 h
  AUC0-inf: N/A (insufficient data)
```

---

## ADaM Import

### ADPC Dataset

```julia
# Import ADaM PK concentrations
adpc = import_adam_adpc("adpc.csv")

# ADaM-specific variables:
# AVAL - Analysis value
# ATPT - Analysis timepoint
# ATPTN - Analysis timepoint (numeric)
# TRTP - Planned treatment
# TRTA - Actual treatment
```

### ADSL Dataset

```julia
# Import subject-level dataset
adsl = import_adam_adsl("adsl.csv")

# Subject-level variables merged automatically
# with concentration data
```

### Combined ADaM Import

```julia
data = import_adam(
    adpc_path = "adpc.csv",
    adsl_path = "adsl.csv"
)
```

---

## Helper Functions

### Data Access

```julia
# Get all subject IDs
ids = subject_ids(data)

# Get number of subjects
n = n_subjects(data)

# Get total observations
n_obs = n_observations(data)

# Get pooled times and observations
all_t = all_times(data)
all_c = all_observations(data)
```

### Data Filtering

```julia
# Filter by treatment arm
arm_a = filter_by_arm(data, "TRT_A")

# Filter by sex
males = filter_by_covariate(data, :sex, :M)

# Filter by dose
high_dose = filter_by_dose(data, dose -> dose >= 200)
```

### Export

```julia
# Export to DataFrame
df = to_dataframe(data)

# Export to NONMEM format
export_nonmem_data(data, "nm_data.csv")

# Export to Monolix format
export_monolix_data(data, "mlx_data.csv")
```

---

## See Also

- [NONMEM Import](nonmem.md) - Import NONMEM control files
- [Monolix Import](monolix.md) - Import Monolix projects
- [Population NCA](../nca/population-nca.md) - Population NCA analysis

