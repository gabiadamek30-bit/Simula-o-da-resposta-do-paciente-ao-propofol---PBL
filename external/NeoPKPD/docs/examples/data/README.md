# Data Import Examples

Examples demonstrating import of clinical data in CDISC and custom formats.

## CDISC Domains

| Domain | Description | Use |
|--------|-------------|-----|
| PC | Pharmacokinetic Concentrations | Observed PK data |
| EX | Exposure | Dosing records |
| DM | Demographics | Subject characteristics |
| VS | Vital Signs | Covariates (weight, etc.) |
| LB | Lab Results | Biomarkers, renal function |

## Examples

| Example | Description | Directory |
|---------|-------------|-----------|
| Basic PC/EX | Standard PK data import | [01_basic_pc_ex](cdisc/01_basic_pc_ex/README.md) |
| With Infusion | Infusion data from EXDUR | [02_with_infusion](cdisc/02_with_infusion/README.md) |
| XPT Format | SAS transport files | [03_xpt_format](cdisc/03_xpt_format/README.md) |

## Basic Usage

```python
from neopkpd.data import load_cdisc

# Load from CSV files
data = load_cdisc(
    pc="pc.csv",
    ex="ex.csv",
    dm="dm.csv"
)

# Access as NeoPKPD format
print(f"Subjects: {len(data.subjects)}")
print(f"Observations: {len(data.observations)}")
print(f"Doses: {len(data.doses)}")
```

## CDISC to NeoPKPD Mapping

### PC Domain

| CDISC Variable | NeoPKPD |
|----------------|----------|
| USUBJID | subject_id |
| PCTPTNUM / PCELTM | time |
| PCSTRESN | dv |
| PCSPEC | compartment |
| PCLLOQ | lloq |

### EX Domain

| CDISC Variable | NeoPKPD |
|----------------|----------|
| USUBJID | subject_id |
| EXSTDTC + EXSTTM | dose_time |
| EXDOSE | amount |
| EXDUR | duration (for infusion) |
| EXROUTE | route |

### DM Domain

| CDISC Variable | NeoPKPD |
|----------------|----------|
| USUBJID | subject_id |
| AGE | covariates.age |
| SEX | covariates.sex |
| RACE | covariates.race |

## See Also

- [Model Import](../import/README.md) - NONMEM/Monolix import
- [Population Examples](../population/README.md) - Using covariates
