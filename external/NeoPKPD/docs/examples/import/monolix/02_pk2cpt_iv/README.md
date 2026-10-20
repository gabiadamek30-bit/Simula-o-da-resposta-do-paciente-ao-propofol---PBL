# Monolix pk2cpt Import

Import a two-compartment IV model from Monolix project.

## Monolix Model

- **Structural**: lib:bolus2_2cpt_ClV1QV2
- **Parameters**: Cl=5.0 L/h, V1=10.0 L, Q=2.0 L/h, V2=20.0 L
- **IIV**: Cl, V1 (log-normal)
- **Error**: Proportional

## Files

| File | Description |
|------|-------------|
| [project.mlxtran](project.mlxtran) | Monolix project file |
| [convert.py](convert.py) | Python conversion |
| [expected.json](expected.json) | Expected NeoPKPD spec |

## Parameter Mapping

| Monolix | NeoPKPD |
|---------|----------|
| Cl_pop | params.CL |
| V1_pop | params.V1 |
| Q_pop | params.Q |
| V2_pop | params.V2 |
