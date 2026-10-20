# CDISC Data Import

Comprehensive guide for importing CDISC/SDTM and ADaM formatted data in Python.

---

## Overview

The `neopkpd.data` module provides utilities for importing data from CDISC (Clinical Data Interchange Standards Consortium) formats.

```python
from neopkpd.data import import_cdisc

data = import_cdisc(
    pc_path="pc.csv",
    ex_path="ex.csv",
    dm_path="dm.csv"
)
```

---

## Quick Start

### Basic Import

```python
from neopkpd.data import import_cdisc

# Import from CSV files
data = import_cdisc(
    pc_path="pc.csv",      # Pharmacokinetic Concentrations
    ex_path="ex.csv",      # Exposure (Dosing)
    dm_path="dm.csv"       # Demographics
)

# Access data
print(f"Subjects: {data.n_subjects}")
print(f"Observations: {data.n_observations}")
```

### Import from SAS Transport (XPT)

```python
data = import_cdisc(
    pc_path="pc.xpt",
    ex_path="ex.xpt",
    dm_path="dm.xpt",
    format="xpt"
)
```

---

## Supported Domains

### PC Domain (Pharmacokinetic Concentrations)

#### Required Variables

| Variable | Label | Type | Description |
|----------|-------|------|-------------|
| STUDYID | Study ID | str | Study identifier |
| USUBJID | Subject ID | str | Unique subject identifier |
| PCTESTCD | Test Code | str | Short test name |
| PCSTRESN | Numeric Result | float | Standardized numeric result |
| PCSTRESU | Unit | str | Units (ng/mL, µg/L, etc.) |
| PCELTM | Elapsed Time | str | Time from reference (ISO 8601) |

#### Optional Variables

| Variable | Label | Description |
|----------|-------|-------------|
| PCLLOQ | Lower LOQ | Lower limit of quantification |
| PCSTAT | Status | Completion status |
| PCBLFL | Baseline Flag | "Y" for baseline |
| PCSPEC | Specimen | PLASMA, SERUM, etc. |
| PCDTC | Datetime | ISO 8601 datetime |
| PCDY | Study Day | Study day |

#### Example Data

```python
import pandas as pd

pc_data = pd.DataFrame({
    'STUDYID': ['STUDY01'] * 7,
    'USUBJID': ['SUBJ001'] * 7,
    'PCTESTCD': ['DRUGA'] * 7,
    'PCSTRESN': [0.0, 125.3, 89.7, 52.1, 28.4, 12.1, 3.2],
    'PCSTRESU': ['ng/mL'] * 7,
    'PCELTM': ['PT0H', 'PT1H', 'PT2H', 'PT4H', 'PT8H', 'PT12H', 'PT24H'],
    'PCLLOQ': [0.5] * 7,
    'PCSPEC': ['PLASMA'] * 7
})
```

---

### EX Domain (Exposure)

#### Required Variables

| Variable | Label | Type | Description |
|----------|-------|------|-------------|
| STUDYID | Study ID | str | Study identifier |
| USUBJID | Subject ID | str | Unique subject identifier |
| EXDOSE | Dose | float | Administered dose |
| EXDOSU | Units | str | Dose units (mg, µg) |
| EXROUTE | Route | str | Route of administration |
| EXSTDTC | Start Datetime | str | Dose start (ISO 8601) |

#### Optional Variables

| Variable | Label | Description |
|----------|-------|-------------|
| EXENDTC | End Datetime | Dose end (infusions) |
| EXDOSFRM | Dose Form | TABLET, CAPSULE, etc. |
| EXDOSFRQ | Frequency | QD, BID, Q12H |
| EXDUR | Duration | ISO 8601 duration |
| EXTRT | Treatment | Treatment name |
| EXDY | Study Day | Study day |

#### Route Mapping

| EX Route | NeoPKPD Route |
|----------|----------------|
| ORAL | `"oral"` |
| INTRAVENOUS | `"iv_bolus"` |
| INTRAVENOUS INFUSION | `"iv_infusion"` |
| SUBCUTANEOUS | `"subcutaneous"` |
| INTRAMUSCULAR | `"intramuscular"` |

---

### DM Domain (Demographics)

#### Required Variables

| Variable | Label | Type | Description |
|----------|-------|------|-------------|
| STUDYID | Study ID | str | Study identifier |
| USUBJID | Subject ID | str | Unique subject identifier |
| RFSTDTC | Reference Start | str | Reference start date |

#### Optional Variables

| Variable | Label | Description |
|----------|-------|-------------|
| AGE | Age | Age at screening |
| AGEU | Age Units | YEARS, MONTHS |
| SEX | Sex | M, F, U |
| RACE | Race | Racial designation |
| ETHNIC | Ethnicity | Ethnic group |
| ARMCD | Arm Code | Treatment arm code |
| ARM | Arm | Arm description |
| COUNTRY | Country | Country |
| SITEID | Site ID | Study site |

---

## Import Functions

### import_cdisc()

```python
from neopkpd.data import import_cdisc

def import_cdisc(
    pc_path: str,
    ex_path: str,
    dm_path: str,
    *,
    format: str = "csv",
    blq_handling: str = "zero",
    lloq: float | None = None,
    time_reference: str = "first_dose"
) -> CDISCData:
    """
    Import CDISC data from PC, EX, and DM domains.

    Parameters
    ----------
    pc_path : str
        Path to PC domain file
    ex_path : str
        Path to EX domain file
    dm_path : str
        Path to DM domain file
    format : str
        "csv" or "xpt"
    blq_handling : str
        "zero", "lloq_half", or "missing"
    lloq : float, optional
        Lower limit of quantification
    time_reference : str
        "first_dose" or "rfstdtc"

    Returns
    -------
    CDISCData
        Imported data object
    """
```

### Individual Domain Functions

```python
from neopkpd.data import read_cdisc_pc, read_cdisc_ex, read_cdisc_dm

# Read individual domains
pc = read_cdisc_pc("pc.csv")
ex = read_cdisc_ex("ex.csv")
dm = read_cdisc_dm("dm.csv")

# Combine manually
from neopkpd.data import cdisc_to_observed_data
data = cdisc_to_observed_data(pc, ex, dm)
```

---

## CDISCData Class

### Structure

```python
@dataclass
class CDISCData:
    """Container for imported CDISC data."""

    # Study information
    study_id: str
    analyte: str
    units: str
    time_units: str

    # Subject data
    subjects: list[SubjectData]

    # Validation
    validation_warnings: list[str]

    # Properties
    @property
    def n_subjects(self) -> int: ...

    @property
    def n_observations(self) -> int: ...

    @property
    def subject_ids(self) -> list[str]: ...

    # Methods
    def get_subject(self, subject_id: str) -> SubjectData: ...
    def to_dataframe(self) -> pd.DataFrame: ...
    def to_estimation_data(self) -> EstimationData: ...
```

### SubjectData Class

```python
@dataclass
class SubjectData:
    """Individual subject data."""

    subject_id: str
    times: list[float]           # Hours from first dose
    observations: list[float]    # Concentrations
    doses: list[DoseEvent]       # Dose events
    covariates: dict[str, Any]   # Demographics
    lloq: float                  # LLOQ value
    blq_flags: list[bool]        # BLQ indicators
```

### Accessing Data

```python
from neopkpd.data import import_cdisc

data = import_cdisc("pc.csv", "ex.csv", "dm.csv")

# Iterate through subjects
for subj in data.subjects:
    print(f"Subject: {subj.subject_id}")
    print(f"  Time points: {len(subj.times)}")
    print(f"  Max conc: {max(subj.observations):.2f}")
    print(f"  Age: {subj.covariates.get('age', 'N/A')}")
    print(f"  Doses: {len(subj.doses)}")

# Get specific subject
subj = data.get_subject("SUBJ001")

# Convert to DataFrame
df = data.to_dataframe()
print(df.head())
```

---

## Time Handling

### Elapsed Time Parsing

```python
# PCELTM format (ISO 8601 duration)
# PT0H     → 0.0 hours
# PT1H     → 1.0 hours
# PT2H30M  → 2.5 hours
# P1DT2H   → 26.0 hours

# Parsed automatically during import
data = import_cdisc("pc.csv", "ex.csv", "dm.csv")

# Times are relative to first dose (in hours)
for subj in data.subjects:
    print(f"Time range: {min(subj.times)} to {max(subj.times)} h")
```

### Reference Time Options

```python
# Times relative to first dose (default)
data = import_cdisc("pc.csv", "ex.csv", "dm.csv",
    time_reference="first_dose"
)

# Times relative to reference start date (RFSTDTC)
data = import_cdisc("pc.csv", "ex.csv", "dm.csv",
    time_reference="rfstdtc"
)
```

---

## BLQ Handling

### Detection

BLQ values are detected from:
- `PCSTRESN == 0` or missing
- `PCSTAT == "NOT DONE"`
- `PCORRES` contains "<", "BLQ", "BLOQ"
- `PCSTRESN < PCLLOQ`

```python
# Check BLQ flags
for subj in data.subjects:
    n_blq = sum(subj.blq_flags)
    if n_blq > 0:
        print(f"Subject {subj.subject_id}: {n_blq} BLQ samples")
```

### Handling Options

```python
# Set BLQ to zero (default)
data = import_cdisc("pc.csv", "ex.csv", "dm.csv",
    blq_handling="zero"
)

# Set BLQ to LLOQ/2
data = import_cdisc("pc.csv", "ex.csv", "dm.csv",
    blq_handling="lloq_half",
    lloq=0.5
)

# Keep as NaN (missing)
data = import_cdisc("pc.csv", "ex.csv", "dm.csv",
    blq_handling="missing"
)
```

---

## Validation

### Automatic Validation

```python
data = import_cdisc("pc.csv", "ex.csv", "dm.csv")

# Check for warnings
if data.validation_warnings:
    print("Validation warnings:")
    for warning in data.validation_warnings:
        print(f"  ⚠️ {warning}")
```

### Manual Validation

```python
from neopkpd.data import validate_cdisc_dataset

validation = validate_cdisc_dataset(
    pc_path="pc.csv",
    ex_path="ex.csv",
    dm_path="dm.csv"
)

if validation.valid:
    print("✓ Dataset passes validation")
else:
    print("❌ Validation failed:")
    for error in validation.errors:
        print(f"  - {error}")
```

### Validation Checks

| Check | Description |
|-------|-------------|
| Required columns | All required variables present |
| Data types | Correct types for each column |
| Subject consistency | IDs match across domains |
| Dose availability | All PK subjects have dosing data |
| Concentration values | No negative concentrations |
| Time sequence | Valid time progression |

---

## Covariate Extraction

### Standard Covariates

```python
# Covariates extracted from DM domain:
# - age (AGE)
# - sex (SEX → "M", "F")
# - race (RACE)
# - ethnicity (ETHNIC)
# - treatment arm (ARM, ARMCD)

data = import_cdisc("pc.csv", "ex.csv", "dm.csv")

for subj in data.subjects:
    print(f"Subject {subj.subject_id}:")
    print(f"  Age: {subj.covariates.get('age', 'N/A')}")
    print(f"  Sex: {subj.covariates.get('sex', 'N/A')}")
    print(f"  Weight: {subj.covariates.get('weight', 'N/A')}")
```

### Additional Covariates from VS

```python
# Include vital signs for weight, height
data = import_cdisc(
    pc_path="pc.csv",
    ex_path="ex.csv",
    dm_path="dm.csv",
    vs_path="vs.csv",  # Optional vital signs
    covariates=["weight", "height", "bmi"]
)
```

---

## Export Options

### To DataFrame

```python
df = data.to_dataframe()

# DataFrame columns:
# SUBJID, TIME, DV, AMT, EVID, MDV, CMT, AGE, SEX, WT, ...
print(df.head())
```

### To NONMEM Format

```python
from neopkpd.data import export_nonmem_data

export_nonmem_data(data, "nonmem_data.csv")

# Creates NONMEM-compatible dataset with:
# ID, TIME, DV, AMT, EVID, MDV, CMT, ...
```

### To NCA Input

```python
# Direct use with NCA
from neopkpd.nca import run_population_nca

pop_result = run_population_nca(
    data.to_dataframe(),
    subject_col='SUBJID',
    time_col='TIME',
    conc_col='DV',
    dose_col='AMT'
)
```

---

## Complete Example

```python
import pandas as pd
from neopkpd.data import import_cdisc, validate_cdisc_dataset

# Step 1: Validate data
print("Validating CDISC data...")
validation = validate_cdisc_dataset(
    pc_path="pc.csv",
    ex_path="ex.csv",
    dm_path="dm.csv"
)

if not validation.valid:
    print("Validation failed:")
    for error in validation.errors:
        print(f"  - {error}")
    exit(1)

print("✓ Validation passed")

# Step 2: Import data
print("\nImporting CDISC data...")
data = import_cdisc(
    pc_path="pc.csv",
    ex_path="ex.csv",
    dm_path="dm.csv",
    blq_handling="lloq_half",
    lloq=0.5
)

# Step 3: Display summary
print("\n" + "=" * 50)
print("CDISC Data Import Summary")
print("=" * 50)

print(f"\nStudy: {data.study_id}")
print(f"Analyte: {data.analyte}")
print(f"Units: {data.units}")

print(f"\nSubjects: {data.n_subjects}")
print(f"Total observations: {data.n_observations}")

# Step 4: Subject details
print("\n--- Subject Details ---")
for subj in data.subjects[:5]:  # First 5 subjects
    n_obs = len(subj.times)
    n_doses = len(subj.doses)
    age = subj.covariates.get('age', 'N/A')
    sex = subj.covariates.get('sex', 'N/A')
    n_blq = sum(subj.blq_flags)

    print(f"\n{subj.subject_id}:")
    print(f"  Observations: {n_obs}")
    print(f"  Time range: {min(subj.times):.1f} - {max(subj.times):.1f} h")
    print(f"  Doses: {n_doses}")
    print(f"  Age: {age}, Sex: {sex}")
    if n_blq > 0:
        print(f"  BLQ samples: {n_blq}")

# Step 5: Check warnings
if data.validation_warnings:
    print("\n--- Warnings ---")
    for w in data.validation_warnings:
        print(f"  ⚠️ {w}")

# Step 6: Use with NCA
print("\n--- Running Population NCA ---")
from neopkpd.nca import run_population_nca, summarize_population_nca

df = data.to_dataframe()
pop_result = run_population_nca(
    df,
    subject_col='SUBJID',
    time_col='TIME',
    conc_col='DV',
    dose_col='AMT'
)

summary = summarize_population_nca(pop_result)
print(f"\nCmax: {summary['cmax']['mean']:.2f} ± {summary['cmax']['sd']:.2f} {data.units}")
print(f"AUC: {summary['auc_0_inf']['mean']:.2f} ± {summary['auc_0_inf']['sd']:.2f} {data.units}·h")

# Step 7: Export
print("\n--- Exporting Data ---")
df.to_csv("analysis_dataset.csv", index=False)
print("Saved to analysis_dataset.csv")
```

---

## ADaM Import

### ADPC Dataset

```python
from neopkpd.data import import_adam_adpc

adpc = import_adam_adpc("adpc.csv")

# ADaM-specific variables:
# AVAL - Analysis value
# ATPT - Analysis timepoint
# ATPTN - Numeric timepoint
# TRTP - Planned treatment
# TRTA - Actual treatment
```

### ADSL Dataset

```python
from neopkpd.data import import_adam_adsl

adsl = import_adam_adsl("adsl.csv")
```

### Combined Import

```python
from neopkpd.data import import_adam

data = import_adam(
    adpc_path="adpc.csv",
    adsl_path="adsl.csv"
)
```

---

## PP Domain (PK Parameters)

Import pre-calculated PK parameters:

```python
from neopkpd.data import read_cdisc_pp

pp = read_cdisc_pp("pp.csv")

# Common parameter codes:
# CMAX - Maximum concentration
# TMAX - Time to max
# AUCLST - AUC to last
# AUCIFO - AUC extrapolated
# LAMZHL - Half-life
```

---

## See Also

- [Data Import Overview](index.md) - All import options
- [Population NCA](../nca/population.md) - NCA analysis
- [Estimation](../estimation/index.md) - Population modeling

