# NONMEM Migration Use Case

Complete workflow for migrating NONMEM models to NeoPKPD, demonstrating
how to convert an existing NONMEM population PK model and validate that
simulation results match.

## Source Model

This example migrates a one-compartment oral PK model from NONMEM ADVAN2:

```
$PROBLEM THEOPHYLLINE PK - ORAL ONE COMPARTMENT
$DATA theo_sd.csv IGNORE=@
$INPUT ID TIME DV AMT EVID CMT WT
$SUBROUTINE ADVAN2 TRANS2

$PK
TVCL = THETA(1) * (WT/70)**0.75
TVV  = THETA(2) * (WT/70)
TVKA = THETA(3)
CL = TVCL * EXP(ETA(1))
V  = TVV  * EXP(ETA(2))
KA = TVKA * EXP(ETA(3))
S2 = V

$ERROR
IPRED = F
Y = F * (1 + ERR(1))

$THETA
(0, 2.8)    ; CL
(0, 35)     ; V
(0, 1.5)    ; KA

$OMEGA
0.09        ; IIV CL
0.04        ; IIV V
0.16        ; IIV KA

$SIGMA
0.04        ; Proportional error
```

## Workflow Steps

| Step | Description |
|------|-------------|
| 1 | Parse NONMEM control file structure |
| 2 | Map ADVAN2 to OneCompOralFirstOrder |
| 3 | Extract THETA, OMEGA values |
| 4 | Convert covariate model |
| 5 | Simulate population and compare |
| 6 | Validate numerical equivalence |

## Quick Run

```bash
# Run migration workflow
julia --project=packages/core docs/examples/use_cases/nonmem_migration/run.jl

# Validate outputs
julia --project=packages/core docs/examples/use_cases/nonmem_migration/validate.jl
```

## Migration Mapping

### Model Type
| NONMEM | NeoPKPD |
|--------|----------|
| ADVAN1 | OneCompIVBolus |
| ADVAN2 | OneCompOralFirstOrder |
| ADVAN3 | TwoCompIVBolus |
| ADVAN4 | TwoCompOral |

### Parameters
| NONMEM | NeoPKPD |
|--------|----------|
| THETA(1) | typical_cl |
| THETA(2) | typical_v |
| THETA(3) | typical_ka |
| OMEGA(1,1) | omega_cl |
| OMEGA(2,2) | omega_v |
| OMEGA(3,3) | omega_ka |

### Covariates
| NONMEM | NeoPKPD |
|--------|----------|
| (WT/70)**0.75 | PowerCovariate(exponent=0.75, reference=70.0) |
| (WT/70) | PowerCovariate(exponent=1.0, reference=70.0) |

## Validation Criteria

The migration is considered successful if:
- Simulated concentrations match within 0.1% relative error
- Population summary statistics (mean, percentiles) align
- Covariate effects produce identical parameter adjustments

## Files

| File | Description |
|------|-------------|
| `run001.ctl` | Original NONMEM control file |
| `run.jl` | Migration and validation script |
| `output/` | Comparison results |

## See Also

- [Model Import Examples](../../import/README.md)
- [Estimation Examples](../../estimation/README.md)
- [Population Examples](../../population/README.md)
