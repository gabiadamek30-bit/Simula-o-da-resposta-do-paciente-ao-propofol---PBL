# Theophylline theo_sd Real-World Validation

## Dataset

- Source: `theo_sd` from nlmixr2data (GPL >= 3)
- Columns: ID, TIME, DV, AMT, EVID, CMT, WT

## Design Realism Addressed

1. **Missing DV rows**: Ignored in RMSE calculation (kept for audit, excluded from metric).

2. **Duplicate observation times**: Supported via unique grid with mapping back to all observations.

3. **AMT ambiguity**: Handled explicitly with deterministic rule:
   - If AMT < 50 (typical mg/kg scale): `dose_mg = AMT * WT`
   - Otherwise: `dose_mg = AMT`
   - Rule is recorded in artifact metadata as `dose_unit_rule`

4. **WT covariate scaling scenario**: Two simulation modes:
   - **Fixed**: CL = 2.75 L/hr, V = 31.8 L (population typical)
   - **WT-scaled**:
     - `CL_i = CL_ref * (WT/70)^0.75`
     - `V_i = V_ref * (WT/70)^1.0`

   This is simulation-only validation (no fitting).

## Model Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| Ka | 1.59 | 1/hr |
| CL_ref | 2.75 | L/hr |
| V_ref | 31.8 | L |
| WT_ref | 70.0 | kg |

## Outputs

Per subject (12 subjects total):
- `subj_<id>_fixed.json` - Fixed parameter artifact
- `subj_<id>_wt_scaled.json` - WT-scaled parameter artifact

Metrics JSON with per-subject:
- `id`, `wt`, `dose_mg`
- `dose_unit_rule` - How AMT was interpreted
- `rmse_fixed` - RMSE for fixed params
- `rmse_wt_scaled` - RMSE for WT-scaled params

## Running

```bash
# From repository root
julia docs/examples/real_world_validation/datasets/theophylline_theo_sd/run.jl

# Copy to expected (first time or when regenerating)
rm -f docs/examples/real_world_validation/studies/theophylline_theo_sd/expected/*.json
cp docs/examples/real_world_validation/studies/theophylline_theo_sd/output/*.json \
   docs/examples/real_world_validation/studies/theophylline_theo_sd/expected/

# Validate
julia docs/examples/real_world_validation/datasets/theophylline_theo_sd/validate.jl
```

## Validation

- Metrics compared with strict tolerance (1e-12)
- All artifacts replayed and compared
- Both fixed and WT-scaled scenarios validated
