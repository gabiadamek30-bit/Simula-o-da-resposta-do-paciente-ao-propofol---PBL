# NONMEM ADVAN2 Import

Import a one-compartment oral model with first-order absorption from NONMEM.

## NONMEM Model

- **ADVAN**: 2 (One-compartment with first-order absorption)
- **TRANS**: 2 (CL, V, Ka parameterization)
- **Parameters**: Ka=1.5 1/h, CL=5.0 L/h, V=50.0 L
- **IIV**: Ka (40% CV), CL (30% CV), V (25% CV)
- **Error**: Combined (additive + proportional)

## Files

| File | Description |
|------|-------------|
| [run002.ctl](run002.ctl) | NONMEM control file |
| [convert.py](convert.py) | Python conversion |
| [expected.json](expected.json) | Expected NeoPKPD spec |

## Parameter Mapping

| NONMEM | NeoPKPD |
|--------|----------|
| THETA(1) = Ka | params.Ka |
| THETA(2) = CL | params.CL |
| THETA(3) = V | params.V |
| OMEGA(1,1) | iiv.omegas.Ka |
| OMEGA(2,2) | iiv.omegas.CL |
| OMEGA(3,3) | iiv.omegas.V |
