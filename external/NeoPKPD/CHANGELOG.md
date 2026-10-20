# Changelog

All notable changes to NeoPKPD will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2025-01-11

### Fixed
- CI segmentation fault in Python tests by adding proper juliacall signal handling
- Added `conftest.py` with session-scoped Julia initialization fixture
- Set `PYTHON_JULIACALL_HANDLE_SIGNALS=yes` environment variable in CI

### Added
- SoftwareX publication figures and diagrams
- TikZ architecture diagrams for publication
- Mermaid architecture diagrams

### Changed
- Updated all version references to 0.1.1

## [0.1.0] - 2025-01-10

### PK Models
- One-compartment IV bolus (`OneCompIVBolus`)
- One-compartment oral first-order absorption (`OneCompOralFirstOrder`)
- Two-compartment IV bolus (`TwoCompIVBolus`)
- Two-compartment oral (`TwoCompOral`)
- Three-compartment IV bolus (`ThreeCompIVBolus`)
- Transit absorption compartment model (`TransitAbsorption`)
- Michaelis-Menten elimination (`MichaelisMentenElimination`)
- IV infusion support with duration and overlapping infusions

### PD Models
- Direct Emax (`DirectEmax`)
- Sigmoid Emax (`SigmoidEmax`)
- Biophase equilibration (`BiophaseEquilibration`)
- Indirect response turnover (`IndirectResponseTurnover`)
- Indirect response IRM1-4 models
- Transit compartment PD (`TransitCompartmentPD`)
- Disease progression models
- Tolerance and counter-regulation

### TMDD Models (Target-Mediated Drug Disposition)
- Full TMDD model for biologics
- Quasi-steady-state (QSS) approximation
- Michaelis-Menten approximation
- TMDD population simulation support
- TMDD-specific analysis functions

### Drug Interaction Models
- Bliss independence
- Competitive inhibition
- Receptor regulation

### Population Simulation
- Inter-individual variability (IIV) with log-normal distribution
- Inter-occasion variability (IOV)
- Static covariates (linear, power, exponential effects)
- Time-varying covariates (step and linear interpolation)
- Below limit of quantification (BLQ) handling

### Parameter Estimation (NLME)
- First-Order Conditional Estimation with Interaction (FOCE-I)
- Stochastic Approximation Expectation Maximization (SAEM)
- Laplacian estimation method
- Bayesian estimation with AdvancedHMC
- Standard error computation (sandwich estimator, bootstrap)
- Covariance matrix estimation
- Objective function value (OFV) computation
- Eta shrinkage calculation

### Advanced Estimation Features
- Bootstrap analysis with confidence intervals
- Mixture models for subpopulation identification
- Model averaging (AIC, BIC, stacked weighting with NNLS)
- Stepwise covariate modeling (SCM) - forward/backward selection
- Parallel MCMC chains
- Optimizer fallback strategies

### Non-Compartmental Analysis (NCA)
- FDA/EMA-compliant NCA metrics
- AUC (linear, log-linear, linear-up/log-down trapezoidal)
- Cmax, Tmax, half-life, clearance, volume of distribution
- Bioavailability and accumulation ratios
- Lambda-z estimation with R-squared criteria
- Partial AUC calculations

### Visual Predictive Check (VPC)
- Standard VPC with prediction intervals
- Prediction-corrected VPC (pcVPC)
- Stratification by covariates
- Multiple binning strategies (quantile, equal width, K-means, Jenks)
- Bootstrap confidence intervals
- LLOQ handling in VPC

### Diagnostics
- Conditional weighted residuals (CWRES)
- Individual weighted residuals (IWRES)
- Normalized prediction distribution errors (NPDE)
- Goodness-of-fit plots
- Estimation diagnostics

### Sensitivity Analysis
- Local sensitivity analysis (single parameter perturbation)
- Population-level sensitivity
- Sobol' indices (first-order, second-order, total-order)
- Morris screening method (elementary effects)
- Sensitivity metrics (AUC of delta, max absolute delta, RMSE)

### Clinical Trial Simulation
- Parallel group designs
- Crossover designs (2x2, 3x3, Williams)
- Dose-escalation designs (3+3, CRM, mTPI, BOIN)
- Bioequivalence studies (ABE, RSABE, NTID)
- Adaptive trial designs
- Virtual population generation
- Power analysis
- Treatment arm specification
- Dosing regimen builders (QD, BID, TID, loading dose)

### Adaptive Trial Features
- Response-adaptive randomization (Thall-Wathen, DBCD)
- Sample size re-estimation
- Treatment selection/arm dropping
- Biomarker enrichment designs
- Interim analysis with alpha spending (O'Brien-Fleming, Pocock)
- Futility boundaries

### Model Import
- NONMEM control stream (.ctl) parsing
- Monolix model (.mlxtran) parsing
- Automatic model structure detection
- Parameter extraction

### Data Import
- CDISC/SDTM format support
- PC (Pharmacokinetic Concentrations) domain
- EX (Exposure) domain
- DM (Demographics) domain
- XPT (SAS Transport) file reader

### Residual Error Models
- Additive error
- Proportional error
- Combined (additive + proportional) error
- Exponential error

### Compliance (FDA 21 CFR Part 11)
- Digital signatures for artifacts
- Audit trail generation
- Environment capture (Julia version, package versions, OS)
- Validation reports
- Schema validation with integrity checks
- Content hashing (SHA-256)

### Serialization & Reproducibility
- JSON artifact format with semantic versioning
- Deterministic replay system
- Golden artifact validation
- Schema migration support (1.0.0 → 1.1.0)
- Compliance metadata in artifacts

### Interfaces
- Julia API (NeoPKPD)
- Python bindings (neopkpd) with JuliaCall
- Command-line interface (CLI)

### CLI Commands
- `simulate` - Run PK/PD simulations
- `nca` - Non-compartmental analysis
- `trial` - Clinical trial simulation
- `vpc` - Visual Predictive Checks
- `import` - Import NONMEM/Monolix models
- `replay` - Replay artifacts deterministically
- `validate-golden` - Validate golden artifacts

### Python Package Features
- Full simulation wrappers for all PK/PD models
- NCA calculations (pure Python + Julia integration)
- Trial simulation interface
- VPC computation and visualization
- Estimation interface
- Visualization module (matplotlib/plotly backends)
- CDISC data utilities

### Visualization (Python)
- Concentration-time plots
- VPC plots with prediction intervals
- Goodness-of-fit plots (DV vs PRED, DV vs IPRED)
- Residual plots (CWRES, IWRES, NPDE)
- Parameter distribution plots
- Sensitivity tornado plots
- Bootstrap forest plots
- Estimation diagnostic dashboards

### Documentation
- MkDocs with Material theme
- Architecture overview
- Sensitivity analysis guide
- Reproducibility guide
- Numerical semantics documentation
- CLI reference

### Testing
- 5,427 Julia unit tests
- 315 Python tests
- Golden artifact validation
- Real-world data validation (Theophylline, Warfarin)
- Cross-platform CI (Ubuntu, macOS, Windows)

### Semantics Versions
- Event Semantics: 1.0.0
- Solver Semantics: 1.0.0
- Artifact Schema: 1.0.0

[0.1.1]: https://github.com/shramish2057/neopkpd/releases/tag/v0.1.1
[0.1.0]: https://github.com/shramish2057/neopkpd/releases/tag/v0.1.0
