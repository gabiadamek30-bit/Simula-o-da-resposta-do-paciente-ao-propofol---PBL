# Basic CDISC PC/EX Import

Import pharmacokinetic data from standard CDISC PC (concentrations) and EX (exposure) domains.

## Data Structure

### PC Domain (Pharmacokinetic Concentrations)
- USUBJID: Subject identifier
- PCTPTNUM: Nominal time
- PCSTRESN: Numeric result
- PCSTAT: Status (if BLQ)
- PCLLOQ: Lower limit of quantification

### EX Domain (Exposure)
- USUBJID: Subject identifier
- EXSTDTC: Start datetime
- EXDOSE: Dose amount
- EXROUTE: Route of administration

## Files

| File | Description |
|------|-------------|
| [pc.csv](pc.csv) | PC domain data |
| [ex.csv](ex.csv) | EX domain data |
| [dm.csv](dm.csv) | DM domain data |
| [load.py](load.py) | Python loader |
| [load.jl](load.jl) | Julia loader |

## Usage

```python
from neopkpd.data import load_cdisc

data = load_cdisc(pc="pc.csv", ex="ex.csv", dm="dm.csv")
```
