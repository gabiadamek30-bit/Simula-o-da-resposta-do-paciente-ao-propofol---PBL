# Warfarin PK/PD real-world validation

Dataset:
- nlmixr2data::warfarin (GPL >= 3)
- dv holds either PK concentration or PD endpoint value
- dvid indicates endpoint identity

Purpose:
- Stress mixed PK/PD observations, dosing records, and deterministic replay.
- Validate that coupled PKPD simulation can be compared against interleaved real observations.

Model:
- PK: 1-comp oral first-order
- PD: indirect response turnover driven by concentration

Approach:
- No fitting.
- Fixed parameters chosen solely to generate stable, reproducible predictions.
- Compare predictions to observed rows by endpoint (dvid).
- Report descriptive RMSE per subject and endpoint.

Outputs:
- One coupled artifact per subject
- Metrics JSON per subject: rmse_pk, rmse_pd
