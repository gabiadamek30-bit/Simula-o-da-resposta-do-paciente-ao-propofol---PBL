module NeoPKPD

# ------------------------------------------------------------------
# External dependencies
# ------------------------------------------------------------------
using SciMLBase
using DifferentialEquations

const NEOPKPD_VERSION = "0.1.1"
export NEOPKPD_VERSION

# ------------------------------------------------------------------
# Core specs and shared types
# ------------------------------------------------------------------
include("specs/time_covariates.jl")
include("specs/specs.jl")
include("specs/error_models.jl")
include("specs/sensitivity.jl")
include("specs/tmdd_specs.jl")

# ------------------------------------------------------------------
# PK model definitions
# ------------------------------------------------------------------
include("pk/onecomp_iv_bolus.jl")
include("pk/onecomp_oral_first_order.jl")
include("pk/twocomp_iv_bolus.jl")
include("pk/twocomp_oral.jl")
include("pk/threecomp_iv_bolus.jl")
include("pk/transit_absorption.jl")
include("pk/michaelis_menten.jl")
include("pk/pk_interface.jl")
include("pk/custom.jl")

# ------------------------------------------------------------------
# PD model definitions
# ------------------------------------------------------------------
include("pd/pd_helpers.jl")
include("pd/direct_emax.jl")
include("pd/sigmoid_emax.jl")
include("pd/indirect_response_turnover.jl")
include("pd/indirect_response_irm1.jl")
include("pd/indirect_response_irm2.jl")
include("pd/indirect_response_irm4.jl")
include("pd/transit_compartment_pd.jl")
include("pd/disease_progression.jl")
include("pd/bliss_independence.jl")
include("pd/competitive_inhibition.jl")
include("pd/drug_interaction.jl")
include("pd/tolerance_counter_regulation.jl")
include("pd/receptor_regulation.jl")
include("pd/biophase_equilibration.jl")

# ------------------------------------------------------------------
# TMDD (Target-Mediated Drug Disposition) Models
# Industry-standard models for biologics and monoclonal antibodies
# ------------------------------------------------------------------
include("tmdd/tmdd.jl")

# ------------------------------------------------------------------
# Numerical semantics 
# These define versioned scientific meaning
# ------------------------------------------------------------------
include("engine/semantics.jl")
include("engine/solver_semantics.jl")
include("engine/semantics_fingerprint.jl")

# ------------------------------------------------------------------
# Perturbation + sensitivity core (local)
# ------------------------------------------------------------------
include("engine/perturb.jl")
include("engine/sensitivity_metrics.jl")
include("engine/sensitivity.jl")
include("engine/sensitivity_population.jl")

# ------------------------------------------------------------------
# Core simulation engine
# ------------------------------------------------------------------
include("engine/events.jl")
include("engine/callbacks.jl")
include("engine/infusion.jl")
include("engine/dose_modifiers.jl")  # ALAG and bioavailability support
include("engine/residual_error.jl")
include("engine/solve.jl")

# ------------------------------------------------------------------
# PK–PD execution layers
# ------------------------------------------------------------------
include("engine/pkpd.jl")
include("engine/pkpd_coupled.jl")

# ------------------------------------------------------------------
# Population engine
# Defines PopulationSpec, PopulationResult, IIV, etc.
# ------------------------------------------------------------------
include("engine/iov.jl")
include("engine/segment_sim.jl")
include("engine/segment_sim_pkpd.jl")
include("engine/time_covariates.jl")
include("engine/covariates.jl")
include("engine/population.jl")

# ------------------------------------------------------------------
# Serialization
# ------------------------------------------------------------------
include("serialization/schema.jl")
include("serialization/serialize.jl")
include("serialization/deserialize.jl")
include("serialization/serialize_population.jl")
include("serialization/deserialize_population.jl")
include("serialization/serialize_sensitivity.jl")
include("serialization/deserialize_sensitivity.jl")
include("serialization/serialize_error.jl")
include("serialization/deserialize_error.jl")
include("serialization/serialize_tmdd.jl")
include("serialization/deserialize_tmdd.jl")

# ------------------------------------------------------------------
# FDA 21 CFR Part 11 Compliance Module
# Audit trails, data integrity, environment capture, validation reports
# ------------------------------------------------------------------
include("compliance/compliance.jl")

# Compliance-aware serialization (extends base serialization with compliance_metadata)
include("serialization/serialize_compliance.jl")
include("serialization/deserialize_compliance.jl")

include("analysis/exposure.jl")
include("analysis/response_metrics.jl")

# ------------------------------------------------------------------
# Model Import (NONMEM, Monolix)
# ------------------------------------------------------------------
include("import/import.jl")

# ------------------------------------------------------------------
# Data Module (CDISC/SDTM support)
# ------------------------------------------------------------------
include("data/data.jl")

# ------------------------------------------------------------------
# Visual Predictive Check (VPC)
# Requires data module for ObservedData
# ------------------------------------------------------------------
include("analysis/vpc_types.jl")
include("analysis/binning.jl")
include("analysis/vpc.jl")

# ------------------------------------------------------------------
# Parameter Estimation (NLME)
# FOCE-I, SAEM, Laplacian methods
# ------------------------------------------------------------------
include("estimation/estimation.jl")

# ------------------------------------------------------------------
# Global Sensitivity Analysis (Sobol', Morris)
# Requires estimation/parallel.jl for ParallelConfig
# ------------------------------------------------------------------
include("engine/sensitivity_sampling.jl")
include("engine/sensitivity_sobol.jl")
include("engine/sensitivity_morris.jl")

# GSA Serialization (after GSA types are defined)
include("serialization/serialize_gsa.jl")
include("serialization/deserialize_gsa.jl")

# ------------------------------------------------------------------
# Residuals (CWRES, IWRES, NPDE)
# Requires estimation module for EstimationResult
# ------------------------------------------------------------------
include("analysis/residuals.jl")
include("analysis/npde.jl")

# ------------------------------------------------------------------
# NCA (Non-Compartmental Analysis)
# FDA/EMA compliant NCA metrics
# ------------------------------------------------------------------
include("nca/nca.jl")

# ------------------------------------------------------------------
# Clinical Trial Simulation
# Parallel, crossover, dose-escalation, adaptive designs
# Virtual population, power analysis, bioequivalence
# ------------------------------------------------------------------
include("trial/trial.jl")

end
