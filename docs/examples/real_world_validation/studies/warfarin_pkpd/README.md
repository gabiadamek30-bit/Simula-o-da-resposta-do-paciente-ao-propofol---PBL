# Warfarin PKPD Real-World Validation

## Dataset

- **Source**: `warfarin` from nlmixr2data (GPL >= 3)
- **Drug**: Warfarin (oral anticoagulant)
- **Subjects**: 32 subjects
- **Endpoints**:
  - PK: Plasma concentration (`cp`)
  - PD: Prothrombin complex activity (`pca`)

## Columns

| Column | Description |
|--------|-------------|
| id | Subject identifier |
| time | Time (hours) |
| amt | Dose amount (mg) |
| dv | Dependent variable |
| dvid | Endpoint type (cp=concentration, pca=PD response) |
| evid | Event ID (1=dose, 0=observation) |
| wt | Body weight (kg) |
| age | Age (years) |
| sex | Sex (male/female) |

## Model

### PK Component
- **Model**: One-compartment oral (ADVAN2)
- **Parameters**:
  | Parameter | Value | Unit |
  |-----------|-------|------|
  | Ka | 1.2 | 1/h |
  | CL | 3.0 | L/h |
  | V | 35.0 | L |

### PD Component
- **Model**: Indirect response turnover (inhibition of production)
- **Equation**: dR/dt = Kin * (1 - Emax * C / (EC50 + C)) - Kout * R
- **Parameters**:
  | Parameter | Value | Unit |
  |-----------|-------|------|
  | Kin | 10.0 | 1/h |
  | Kout | 0.5 | 1/h |
  | R0 | 20.0 | % activity |
  | Emax | 0.8 | - |
  | EC50 | 2.0 | mg/L |

## Validation Approach

This is **simulation-only validation** using fixed literature-derived parameters:
1. Parse warfarin dataset with both PK and PD observations
2. Simulate coupled PKPD using NeoPKPD
3. Compute RMSE separately for PK (concentration) and PD (response) endpoints
4. Compare against expected output (deterministic replay)

## Outputs

Per subject (32 subjects total):
- `subj_<id>.json` - Coupled PKPD simulation artifact

Summary:
- `metrics.json` - Per-subject RMSE for PK and PD endpoints

## Running

```bash
# Run simulation (from repository root)
julia --project=packages/core \
  docs/examples/real_world_validation/datasets/warfarin_nlmixr2data/run.jl

# Copy to expected (first time or when regenerating)
rm -f docs/examples/real_world_validation/studies/warfarin_pkpd/expected/*.json
cp docs/examples/real_world_validation/studies/warfarin_pkpd/output/*.json \
   docs/examples/real_world_validation/studies/warfarin_pkpd/expected/

# Validate
julia --project=packages/core \
  docs/examples/real_world_validation/datasets/warfarin_nlmixr2data/validate.jl
```

## Validation Criteria

- Deterministic replay of all artifacts
- Metrics matched with strict tolerance (1e-12)
- Both PK and PD predictions validated

## Key Features Demonstrated

1. **Coupled PKPD Simulation**: PK drives PD through concentration-effect relationship
2. **Indirect Response Model**: Turnover dynamics with inhibition of production
3. **Multi-endpoint Handling**: Separate PK (cp) and PD (pca) observations
4. **Real-World Data**: Actual clinical trial data structure with covariates

## See Also

- [Theophylline SD Study](../theophylline_theo_sd/README.md)
- [Theophylline MD Study](../theophylline_theo_md/README.md)
