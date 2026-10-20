# PKPD Biomarker Turnover

This use case models exposure-driven biomarker suppression using a coupled PKPD indirect response turnover model.

## Purpose

- Validate coupled PKPD, regimen comparison, and IOV behavior under realistic workflows
- Provide a reproducible, replayable contract for response metrics used in decision making

## Models

**PK Model:** One-compartment IV bolus (CL=10 L/h, V=50 L typical)

**PD Model:** Indirect response turnover
- Baseline R0 = Kin / Kout = 100
- Drug effect: inhibitory Emax (Imax=0.9, IC50=5)

## Scenarios

| Regimen | Doses |
|---------|-------|
| QD | 100 mg at 0, 24, 48 h |
| BID | 50 mg at 0, 12, 24, 36, 48, 60 h |

Each scenario runs with and without IOV on CL across occasions.

## Outputs

- Population artifacts with mean and quantile response summaries
- Decision metrics: Emin, time below 80% baseline, suppression AUC

## Running

```bash
julia docs/examples/use_cases/pkpd_biomarker_turnover/run.jl
julia docs/examples/use_cases/pkpd_biomarker_turnover/validate.jl
```
