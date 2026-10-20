# Theophylline theo_sd real-world validation

Dataset:
- `theo_sd` from nlmixr2data (GPL >= 3)
- Columns: ID, TIME, DV, AMT, EVID, CMT, WT

Design realism addressed:
1. Missing DV rows are ignored (kept for audit, excluded from RMSE).
2. Duplicate observation times are supported.
3. AMT ambiguity is handled explicitly:
   - If AMT is small (typical mg/kg scale), dose_mg = AMT * WT
   - Else dose_mg = AMT
   This rule is recorded in artifact metadata.
4. Optional WT covariate scaling scenario:
   - CL_i = CL_ref * (WT/70)^0.75
   - V_i = V_ref * (WT/70)^1.0
   This is simulation-only, no fitting.

Outputs:
- Per subject artifacts for:
  - fixed params (no WT scaling)
  - WT scaled params
- Metrics JSON capturing RMSE for each scenario.
