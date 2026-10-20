# TMDD Module
#
# Target-Mediated Drug Disposition models for biologics and monoclonal antibodies.
#
# This module provides industry-standard TMDD model implementations:
# - Full TMDD (Mager-Jusko)
# - QSS (Quasi-Steady-State) approximation
# - QE (Quasi-Equilibrium) approximation
# - Michaelis-Menten approximation
# - Rapid binding (Wagner) approximation
# - Irreversible binding
# - Two-compartment variants
# - Soluble target models
# - Receptor internalization/recycling
#
# References:
# - Mager DE, Jusko WJ. J Pharmacokinet Pharmacodyn. 2001;28(6):507-532
# - Gibiansky L, et al. J Pharmacokinet Pharmacodyn. 2008;35(5):573-591
# - Peletier LA, Gabrielsson J. J Theor Biol. 2012;293:39-58

# Include model implementations
include("tmdd_models.jl")
include("tmdd_solve.jl")
include("tmdd_population.jl")
include("tmdd_analysis.jl")

# Re-export all public symbols
# Type definitions are exported from tmdd_specs.jl (included in specs)
# Model functions are exported from tmdd_models.jl
# Solve functions are exported from tmdd_solve.jl
