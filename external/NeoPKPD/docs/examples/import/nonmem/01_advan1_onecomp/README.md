# NONMEM ADVAN1 Import

Import a one-compartment IV bolus model from NONMEM control file.

## NONMEM Model

- **ADVAN**: 1 (One-compartment)
- **TRANS**: 2 (CL, V parameterization)
- **Parameters**: CL=5.0 L/h, V=50.0 L
- **IIV**: CL (30% CV), V (25% CV)
- **Error**: Proportional (10% CV)

## Files

| File | Description |
|------|-------------|
| [run001.ctl](run001.ctl) | NONMEM control file |
| [convert.py](convert.py) | Python conversion |
| [convert.jl](convert.jl) | Julia conversion |
| [expected.json](expected.json) | Expected NeoPKPD spec |

## Conversion

```python
from neopkpd import import_nonmem

model = import_nonmem("run001.ctl")
print(model.spec)  # NeoPKPD ModelSpec
```

## Parameter Mapping

| NONMEM | NeoPKPD |
|--------|----------|
| THETA(1) = CL | params.CL |
| THETA(2) = V | params.V |
| OMEGA(1,1) | iiv.omegas.CL |
| OMEGA(2,2) | iiv.omegas.V |
| SIGMA(1,1) | error.proportional |
