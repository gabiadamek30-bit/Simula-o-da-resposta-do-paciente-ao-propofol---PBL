# theo_md provenance

This dataset is `theo_md` distributed by the R package `nlmixr2data`.

## Description

- Starts with day 1 theophylline concentrations from the nlme/NONMEM teaching dataset
- Later observations are simulated under a once-daily regimen as described by nlmixr2data

## Columns

- ID: Subject identifier
- TIME: Time since first dose (hours)
- DV: Dependent variable (concentration, mg/L)
- AMT: Dose amount
- EVID: Event ID (0=observation, 1=dose)
- CMT: Compartment
- WT: Body weight (kg)
