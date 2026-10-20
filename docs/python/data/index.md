# Data Import

The `neopkpd.data` module provides utilities for importing and preparing PK/PD data from various formats.

---

## Supported Formats

<div class="grid cards" markdown>

-   :material-file-table:{ .lg .middle } **CDISC Formats**

    ---

    SDTM PC, EX, DM domains

    [:octicons-arrow-right-24: CDISC Import](cdisc.md)

-   :material-file-code:{ .lg .middle } **NONMEM Import**

    ---

    Parse NONMEM control files

    [:octicons-arrow-right-24: NONMEM](nonmem.md)

-   :material-file-document:{ .lg .middle } **Monolix Import**

    ---

    Parse Monolix project files

    [:octicons-arrow-right-24: Monolix](monolix.md)

-   :material-file-delimited:{ .lg .middle } **CSV Files**

    ---

    Standard CSV data import

    [:octicons-arrow-right-24: CSV Import](#csv-import)

</div>

---

## Quick Start

### Import CDISC Data

```python
from neopkpd.data import import_cdisc

# Import from CSV files
data = import_cdisc(
    pc_path="pc.csv",      # Pharmacokinetic Concentrations
    ex_path="ex.csv",      # Exposure (Dosing)
    dm_path="dm.csv",      # Demographics
    format="csv"
)

# Access standardized data
print(f"Subjects: {len(set(data.usubjid))}")
print(f"Observations: {len(data.dv)}")
```

### Import from XPT

```python
from neopkpd.data import import_cdisc

data = import_cdisc(
    pc_path="pc.xpt",
    ex_path="ex.xpt",
    dm_path="dm.xpt",
    format="xpt"
)
```

---

## CDISC Domains

### PC (Pharmacokinetic Concentrations)

Required variables:

| Variable | Description | Required |
|----------|-------------|----------|
| USUBJID | Unique subject identifier | Yes |
| PCDTC | Collection datetime | Yes |
| PCSTRESN | Numeric result | Yes |
| PCTEST | Analyte name | Yes |
| PCSPEC | Specimen type | No |
| PCBLFL | Baseline flag | No |

### EX (Exposure)

Required variables:

| Variable | Description | Required |
|----------|-------------|----------|
| USUBJID | Unique subject identifier | Yes |
| EXSTDTC | Start datetime | Yes |
| EXDOSE | Dose amount | Yes |
| EXDOSU | Dose units | Yes |
| EXROUTE | Route of administration | No |
| EXDUR | Duration (for infusions) | No |

### DM (Demographics)

Common variables:

| Variable | Description | Required |
|----------|-------------|----------|
| USUBJID | Unique subject identifier | Yes |
| AGE | Age | No |
| SEX | Sex | No |
| RACE | Race | No |
| WEIGHT | Weight (from VS or custom) | No |

---

## CSV Import

### Standard Format

```python
from neopkpd.data import read_pk_data

# Expected columns: ID, TIME, DV, AMT, EVID, CMT
data = read_pk_data("pk_data.csv")

# Custom column mapping
data = read_pk_data(
    "pk_data.csv",
    column_map={
        "subject_id": "ID",
        "time_hr": "TIME",
        "concentration": "DV",
        "dose_mg": "AMT"
    }
)
```

### NONMEM Format

```python
from neopkpd.data import read_nonmem_data

# Standard NONMEM dataset
data = read_nonmem_data("data.csv")

# Access as DataFrame or structured data
print(data.to_dataframe())
```

---

## XPT Import

```python
from neopkpd.data import read_xpt

# Single file
pc_data = read_xpt("pc.xpt")

# Multiple files
data = {
    "pc": read_xpt("pc.xpt"),
    "ex": read_xpt("ex.xpt"),
    "dm": read_xpt("dm.xpt")
}
```

---

## Data Preparation

### Time Calculation

```python
from neopkpd.data import calculate_time_after_dose

data = calculate_time_after_dose(
    observations=obs_df,
    dosing=dose_df,
    time_column="PCDTC",
    dose_time_column="EXSTDTC"
)
```

### BLQ Handling

```python
from neopkpd.data import handle_blq

# Replace BLQ with LLOQ/2
data = handle_blq(
    data,
    lloq=0.1,
    method="lloq_half"
)

# Or set to zero
data = handle_blq(data, lloq=0.1, method="zero")

# Or drop
data = handle_blq(data, lloq=0.1, method="drop")
```

### Merge Demographics

```python
from neopkpd.data import merge_demographics

merged = merge_demographics(
    pk_data=pc_data,
    demographics=dm_data,
    on="USUBJID"
)
```

---

## Data Classes

### PKData

```python
class PKData:
    ids: list[str]           # Subject IDs
    times: list[float]       # Time points
    dv: list[float]          # Observations
    amt: list[float]         # Doses
    evid: list[int]          # Event IDs
    cmt: list[int]           # Compartments
    covariates: dict         # Covariate data

    def to_dataframe(self) -> pd.DataFrame: ...
    def to_dict(self) -> dict: ...
    def subset(self, ids: list) -> PKData: ...
```

### CDISCData

```python
class CDISCData:
    usubjid: list[str]
    subjid: list[int]
    times: list[float]
    dv: list[float]
    doses: list[DoseEvent]
    demographics: dict

    def to_estimation_data(self) -> EstimationData: ...
```

---

## Examples

### Complete Workflow

```python
from neopkpd.data import import_cdisc, handle_blq, calculate_time_after_dose
import neopkpd

# Import CDISC data
data = import_cdisc(
    pc_path="pc.csv",
    ex_path="ex.csv",
    dm_path="dm.csv"
)

# Handle BLQ
data = handle_blq(data, lloq=0.1, method="lloq_half")

# Calculate time after first dose
data = calculate_time_after_dose(data)

# Convert to estimation format
est_data = data.to_estimation_data()

# Run NCA
from neopkpd.nca import run_population_nca
nca_results = run_population_nca(est_data)
```

---

## Next Steps

- [CDISC Import Details](cdisc.md) - SDTM and ADaM data
- [NONMEM Import](nonmem.md) - Import control stream files
- [Monolix Import](monolix.md) - Import mlxtran projects
- [NCA Module](../nca/index.md) - Non-compartmental analysis
- [Estimation](../estimation/index.md) - Population modeling
