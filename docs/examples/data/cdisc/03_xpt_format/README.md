# CDISC XPT Format Import

Import pharmacokinetic data from SAS transport (.xpt) files.

## XPT Format

XPT (SAS Transport) is the standard format for CDISC submission packages:
- Binary format specified by FDA
- Contains datasets in SAS transport format
- Common in regulatory submissions

## Files

| File | Description |
|------|-------------|
| [load.py](load.py) | Python XPT loader |

## Usage

```python
from neopkpd.data import load_cdisc_xpt

# Load from XPT files
data = load_cdisc_xpt(
    pc="pc.xpt",
    ex="ex.xpt",
    dm="dm.xpt"
)
```

## Dependencies

- `pyreadstat` or `pandas` with SAS support
- XPT files must be valid SAS transport format
