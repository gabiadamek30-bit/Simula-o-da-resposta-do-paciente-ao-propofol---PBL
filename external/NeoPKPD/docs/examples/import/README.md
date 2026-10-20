# Model Import Examples

Examples demonstrating import of models from NONMEM and Monolix.

## Supported Formats

| Source | Format | Import Function |
|--------|--------|-----------------|
| NONMEM | .ctl, .mod | `import_nonmem()` |
| Monolix | .mlxtran | `import_monolix()` |

## NONMEM Examples

| Example | ADVAN | Description | Directory |
|---------|-------|-------------|-----------|
| ADVAN1 | 1-comp IV | One-compartment IV bolus | [01_advan1_onecomp](nonmem/01_advan1_onecomp/README.md) |
| ADVAN2 | 1-comp Oral | One-compartment first-order absorption | [02_advan2_oral](nonmem/02_advan2_oral/README.md) |
| ADVAN4 | 2-comp IV | Two-compartment IV bolus | [03_advan4_twocomp](nonmem/03_advan4_twocomp/README.md) |

## Monolix Examples

| Example | Model | Description | Directory |
|---------|-------|-------------|-----------|
| pk1cpt | 1-comp Oral | One-compartment with Ka | [01_pk1cpt_oral](monolix/01_pk1cpt_oral/README.md) |
| pk2cpt | 2-comp IV | Two-compartment IV | [02_pk2cpt_iv](monolix/02_pk2cpt_iv/README.md) |

## Import Workflow

```python
from neopkpd import import_nonmem, import_monolix

# Import NONMEM model
model = import_nonmem("run001.ctl")

# Import Monolix project
model = import_monolix("project.mlxtran")

# Use imported model
result = simulate(model.spec, t_end=24.0)
```

## What Gets Imported

| Component | NONMEM | Monolix |
|-----------|--------|---------|
| Structural model | ADVAN/TRANS | Structural model |
| Fixed effects (THETA) | Yes | Yes |
| Random effects (OMEGA) | Yes | Yes |
| Residual error (SIGMA) | Yes | Yes |
| Covariate model | $COV block | Covariate section |
| Initial estimates | Yes | Yes |

## Validation

Each import example includes:
1. Original source file
2. Import script
3. Expected NeoPKPD specification
4. Validation by comparing simulations

## See Also

- [Data Import](../data/README.md) - CDISC data import
- [Models](../models/README.md) - NeoPKPD model specifications
