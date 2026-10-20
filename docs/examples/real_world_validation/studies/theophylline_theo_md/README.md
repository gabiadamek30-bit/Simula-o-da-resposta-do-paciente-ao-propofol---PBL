# Theophylline theo_md Multiple-Dose Validation

## Dataset

- Source: `theo_md` from nlmixr2data (GPL >= 3)
- Columns: ID, TIME, DV, AMT, EVID, CMT, WT

## Purpose

Stress test event handling and replay under multiple-dose designs:
- Multiple dosing records per subject
- Mixed dose and observation records
- Long time horizons
- Duplicated times and missing DV tolerance

## Model

One-compartment oral PK with first-order absorption.

Fixed parameters (no fitting):
| Parameter | Value | Unit |
|-----------|-------|------|
| Ka | 1.59 | 1/hr |
| CL | 2.75 | L/hr |
| V | 31.8 | L |

## Dose Handling

- All records with AMT > 0 are interpreted as dosing events
- Dose unit rule is explicit and recorded in metadata:
  - If AMT < 50: treat as mg/kg and multiply by WT
  - Else: treat as mg

## Outputs

- One execution artifact per subject: `subj_<id>.json`
- Metrics JSON with per-subject: id, wt, dose_unit_rule, dose_events, rmse

## Running

```bash
# From repository root
julia docs/examples/real_world_validation/datasets/theophylline_theo_md/run.jl
julia docs/examples/real_world_validation/datasets/theophylline_theo_md/validate.jl
```
