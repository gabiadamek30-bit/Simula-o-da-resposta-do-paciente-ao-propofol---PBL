# NONMEM ADVAN4 Import

Import a two-compartment IV bolus model from NONMEM.

## NONMEM Model

- **ADVAN**: 4 (Two-compartment)
- **TRANS**: 4 (CL, V1, Q, V2 parameterization)
- **Parameters**: CL=5.0 L/h, V1=10.0 L, Q=2.0 L/h, V2=20.0 L
- **IIV**: CL (30% CV), V1 (25% CV)
- **Error**: Proportional (15% CV)

## Files

| File | Description |
|------|-------------|
| [run003.ctl](run003.ctl) | NONMEM control file |
| [convert.py](convert.py) | Python conversion |
| [expected.json](expected.json) | Expected NeoPKPD spec |

## Parameter Mapping

| NONMEM | NeoPKPD |
|--------|----------|
| THETA(1) = CL | params.CL |
| THETA(2) = V1 | params.V1 |
| THETA(3) = Q | params.Q |
| THETA(4) = V2 | params.V2 |
