# Glossary

Common pharmacometrics and NeoPKPD terminology.

---

## A

**Absorption Rate Constant (Ka)**
: First-order rate constant describing the rate of drug absorption from the gut to systemic circulation. Units: 1/h.

**AUC (Area Under the Curve)**
: Total drug exposure over time, calculated as the integral of concentration over time. Units: mg·h/L.

**AUC0-t**
: Area under the concentration-time curve from time zero to the last measurable concentration.

**AUC0-inf**
: Area under the concentration-time curve extrapolated to infinity.

---

## B

**Bioavailability (F)**
: Fraction of administered dose that reaches systemic circulation unchanged. Range: 0-1.

**Bioequivalence (BE)**
: Demonstration that two formulations have similar bioavailability within defined limits (typically 80-125% for the 90% CI of the geometric mean ratio).

**BLQ (Below Limit of Quantification)**
: Concentration measurement below the lower limit of quantification of the analytical assay.

---

## C

**Clearance (CL)**
: Volume of plasma from which drug is completely removed per unit time. Units: L/h.

**Cmax**
: Maximum observed concentration. Units: mg/L.

**Compartment**
: Theoretical space in which drug is assumed to distribute uniformly and instantaneously.

**CWRES (Conditional Weighted Residuals)**
: Residuals weighted by the conditional variance, used for model diagnostics.

---

## D

**Distribution**
: Process by which drug reversibly transfers between blood and tissues.

**Dose Event**
: A single drug administration with time, amount, and optional duration (for infusions).

---

## E

**EBE (Empirical Bayes Estimate)**
: Individual parameter estimates derived using Bayesian estimation with population parameters as priors.

**EC50**
: Concentration producing 50% of maximum effect. Units: mg/L.

**Effect Compartment**
: Hypothetical compartment representing the site of drug action, used to model PK-PD delays.

**Emax**
: Maximum achievable drug effect.

**Eta (η)**
: Random effect representing individual deviation from the typical population parameter value.

---

## F

**FOCE (First-Order Conditional Estimation)**
: NLME estimation method that linearizes the model around individual parameter estimates.

**FOCE-I (FOCE with Interaction)**
: FOCE method that accounts for the interaction between random effects and residual error.

---

## G

**Gamma (γ)**
: Hill coefficient in the sigmoid Emax model, controlling the steepness of the concentration-response curve.

---

## H

**Half-life (t½)**
: Time required for drug concentration to decrease by 50%. Calculated as ln(2)/λz.

**Hill Coefficient**
: See Gamma.

---

## I

**IIV (Inter-Individual Variability)**
: Random variability in pharmacokinetic or pharmacodynamic parameters between subjects in a population.

**Indirect Response Model**
: PD model where drug affects the rate of production or elimination of a response variable rather than the response directly.

**IOV (Inter-Occasion Variability)**
: Random variability in parameters within the same subject across different occasions.

---

## K

**Km (Michaelis Constant)**
: Concentration at which the elimination rate is half of Vmax in saturable kinetics. Units: mg/L.

---

## L

**Lambda_z (λz)**
: Terminal elimination rate constant, estimated from the terminal log-linear phase of the concentration-time curve. Units: 1/h.

**LLOQ (Lower Limit of Quantification)**
: Lowest concentration that can be reliably measured with acceptable precision and accuracy.

---

## M

**MRT (Mean Residence Time)**
: Average time a drug molecule spends in the body. Calculated as AUMC/AUC.

**Michaelis-Menten**
: Nonlinear elimination kinetics characterized by Vmax and Km, where elimination rate saturates at high concentrations.

---

## N

**NCA (Non-Compartmental Analysis)**
: Model-independent method for calculating PK parameters directly from concentration-time data.

**NLME (Nonlinear Mixed Effects)**
: Statistical framework for analyzing population data with both fixed and random effects.

**NPDE (Normalized Prediction Distribution Errors)**
: Diagnostic metric based on the cumulative distribution of observations relative to simulations.

---

## O

**OFV (Objective Function Value)**
: Negative twice the log-likelihood; minimized during parameter estimation. Lower is better for nested models.

**Omega (Ω)**
: Variance-covariance matrix of inter-individual random effects.

---

## P

**pcVPC (Prediction-Corrected VPC)**
: VPC variant that corrects predictions for differences in independent variables across bins.

**Population PK**
: Approach to characterize typical PK parameters and their variability in a population.

---

## Q

**Q (Inter-compartmental Clearance)**
: Rate of drug transfer between compartments. Units: L/h.

---

## R

**Residual Error**
: Unexplained variability between model predictions and observations, including measurement error and model misspecification.

**RSE (Relative Standard Error)**
: Standard error expressed as a percentage of the parameter estimate.

---

## S

**SAEM (Stochastic Approximation EM)**
: NLME estimation algorithm using stochastic approximation for robust parameter estimation.

**Shrinkage**
: Phenomenon where individual parameter estimates are pulled toward population values due to limited individual data.

**Sigma (σ)**
: Variance of residual error.

---

## T

**Theta (θ)**
: Fixed-effect (typical population) parameter values.

**Tmax**
: Time at which Cmax occurs. Units: h.

**TMDD (Target-Mediated Drug Disposition)**
: PK behavior where drug binding to target significantly affects its disposition.

**Transit Compartment**
: Series of compartments used to model delayed or complex absorption processes.

---

## V

**Vmax**
: Maximum rate of saturable elimination. Units: mg/h.

**Volume of Distribution (V)**
: Apparent volume into which drug distributes at steady state. Units: L.

**VPC (Visual Predictive Check)**
: Graphical method for assessing model adequacy by comparing observed data percentiles to simulated prediction intervals.

**Vss (Volume of Distribution at Steady State)**
: Volume relating amount of drug in body to plasma concentration at steady state.

---

## W

**Washout Period**
: Time between treatment periods in crossover studies to allow drug elimination before the next period.

---

## Abbreviations

| Abbreviation | Full Term |
|--------------|-----------|
| AUC | Area Under the Curve |
| BE | Bioequivalence |
| BID | Twice Daily |
| BLQ | Below Limit of Quantification |
| CDISC | Clinical Data Interchange Standards Consortium |
| CI | Confidence Interval |
| CL | Clearance |
| CV | Coefficient of Variation |
| EBE | Empirical Bayes Estimate |
| FDA | Food and Drug Administration |
| FOCE | First-Order Conditional Estimation |
| GOF | Goodness of Fit |
| IIV | Inter-Individual Variability |
| IOV | Inter-Occasion Variability |
| IV | Intravenous |
| NCA | Non-Compartmental Analysis |
| NLME | Nonlinear Mixed Effects |
| NPDE | Normalized Prediction Distribution Error |
| ODE | Ordinary Differential Equation |
| OFV | Objective Function Value |
| PD | Pharmacodynamics |
| PK | Pharmacokinetics |
| QD | Once Daily |
| RSE | Relative Standard Error |
| SAEM | Stochastic Approximation EM |
| SDTM | Study Data Tabulation Model |
| SE | Standard Error |
| TID | Three Times Daily |
| TMDD | Target-Mediated Drug Disposition |
| TOST | Two One-Sided Tests |
| VPC | Visual Predictive Check |
