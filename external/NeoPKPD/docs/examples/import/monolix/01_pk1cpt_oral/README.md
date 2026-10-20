# Monolix pk1cpt Import

Import a one-compartment oral model from Monolix project.

## Monolix Model

- **Structural**: lib:oral1_1cpt_kaVCl
- **Parameters**: ka=1.5 1/h, V=50.0 L, Cl=5.0 L/h
- **IIV**: ka, V, Cl (log-normal)
- **Error**: Combined

## Files

| File | Description |
|------|-------------|
| [project.mlxtran](project.mlxtran) | Monolix project file |
| [convert.py](convert.py) | Python conversion |
| [expected.json](expected.json) | Expected NeoPKPD spec |

## Conversion

```python
from neopkpd import import_monolix

model = import_monolix("project.mlxtran")
print(model.spec)
```

## Parameter Mapping

| Monolix | NeoPKPD |
|---------|----------|
| ka_pop | params.Ka |
| V_pop | params.V |
| Cl_pop | params.CL |
| omega_ka | iiv.omegas.Ka |
| omega_V | iiv.omegas.V |
| omega_Cl | iiv.omegas.CL |
