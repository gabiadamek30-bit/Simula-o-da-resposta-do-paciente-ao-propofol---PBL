# Changelog

All notable changes to the neopkpd Python package are documented here.

## [0.1.1] - 2025-01-11

### Fixed
- CI segmentation fault by adding proper juliacall signal handling in `conftest.py`
- Session-scoped Julia initialization fixture for stable test execution

### Added
- `conftest.py` with `PYTHON_JULIACALL_HANDLE_SIGNALS=yes` configuration
- `julia_initialized` and `init` fixtures for backward compatibility

## [0.1.0] - 2025-01-10

### Simulation
- Full PK model simulation wrappers (1/2/3-compartment, transit, Michaelis-Menten)
- PD model simulation (Emax, sigmoid Emax, indirect response)
- Coupled PK-PD simulation
- IV bolus and infusion support
- Population simulation with IIV/IOV

### Non-Compartmental Analysis (NCA)
- Pure Python NCA implementation
- Julia-integrated NCA for validated results
- AUC, Cmax, Tmax, half-life, clearance calculations
- Multiple AUC methods (linear, log-linear, linear-up/log-down)
- Lambda-z estimation

### Visual Predictive Check (VPC)
- VPC computation via Julia backend
- Standard and prediction-corrected VPC
- Stratification support
- Multiple binning strategies

### Clinical Trial Simulation
- Trial design specification
- Parallel and crossover designs
- Dose-escalation trials
- Power analysis interface

### Parameter Estimation
- FOCE-I, SAEM, Laplacian interface
- Bootstrap analysis
- Estimation diagnostics

### Model Import
- NONMEM control stream parsing
- Monolix model parsing

### Data Utilities
- CDISC/SDTM data handling
- DataFrame integration

### Visualization
- Matplotlib backend (default)
- Plotly backend (interactive)
- Concentration-time plots
- VPC visualization with confidence bands
- Goodness-of-fit plots (DV vs PRED, DV vs IPRED)
- Residual diagnostic plots
- Sensitivity tornado plots
- Bootstrap forest plots
- Estimation summary dashboards

### Infrastructure
- JuliaCall integration for Julia backend
- Automatic Julia initialization
- Bridge module for Julia interop
- Type conversion utilities

### Testing
- 315 unit tests
- Integration tests with Julia backend
- Visualization tests

[0.1.1]: https://github.com/shramish2057/neopkpd/releases/tag/v0.1.1
[0.1.0]: https://github.com/shramish2057/neopkpd/releases/tag/v0.1.0
