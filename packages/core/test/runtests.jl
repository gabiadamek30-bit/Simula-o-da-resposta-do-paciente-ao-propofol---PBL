# NeoPKPD Test Suite
#
# This file serves as the entry point for all tests. Individual test files
# are organized by functionality in separate files.

using Test
using NeoPKPD
using StableRNGs

# Include test helper functions
include("test_helpers.jl")

@testset "NeoPKPD" begin
    # One-compartment PK models
    @testset "One-Compartment PK" begin
        include("test_onecomp.jl")
    end

    # Two-compartment PK models
    @testset "Two-Compartment PK" begin
        include("test_twocomp.jl")
    end

    # Three-compartment PK models
    @testset "Three-Compartment PK" begin
        include("test_threecomp.jl")
    end

    # Transit absorption models
    @testset "Transit Absorption" begin
        include("test_transit.jl")
    end

    # Michaelis-Menten elimination
    @testset "Michaelis-Menten Elimination" begin
        include("test_michaelis_menten.jl")
    end

    # IV Infusion administration
    @testset "IV Infusion" begin
        include("test_infusion.jl")
    end

    # Residual Error Models
    @testset "Residual Error Models" begin
        include("test_residual_error.jl")
    end

    # PKPD coupling tests
    @testset "PKPD Coupling" begin
        include("test_pkpd.jl")
    end

    # Sigmoid Emax PD model
    @testset "Sigmoid Emax PD" begin
        include("test_sigmoid_emax.jl")
    end

    # Biophase equilibration PD model
    @testset "Biophase Equilibration PD" begin
        include("test_biophase.jl")
    end

    # Event semantics
    @testset "Event Semantics" begin
        include("test_event_semantics.jl")
    end

    # Serialization and replay
    @testset "Serialization" begin
        include("test_serialization.jl")
    end

    # Population simulations
    @testset "Population Simulations" begin
        include("test_population.jl")
    end

    # Sensitivity analysis (Local)
    @testset "Sensitivity Analysis" begin
        include("test_sensitivity.jl")
    end

    # Global Sensitivity Analysis (Sobol', Morris)
    @testset "Global Sensitivity Analysis" begin
        include("test_gsa.jl")
    end

    # IOV and covariates
    @testset "IOV and Covariates" begin
        include("test_iov_covariates.jl")
    end

    # PK/PD metrics
    @testset "Metrics" begin
        include("test_metrics.jl")
    end

    # Integration tests (PK + PD combinations)
    @testset "Integration" begin
        include("test_integration.jl")
    end

    # NCA (Non-Compartmental Analysis)
    @testset "NCA" begin
        include("test_nca.jl")
    end

    # Clinical Trial Simulation
    @testset "Trial" begin
        include("test_trial.jl")
    end

    # Model Import (NONMEM, Monolix)
    @testset "Model Import" begin
        include("test_import.jl")
    end

    # CDISC/SDTM Data Support
    @testset "CDISC Data" begin
        include("test_cdisc.jl")
    end

    # XPT (SAS Transport) Format Support
    @testset "XPT Reader" begin
        include("test_xpt.jl")
    end

    # Visual Predictive Check (VPC)
    @testset "VPC" begin
        include("test_vpc.jl")
    end

    # NPDE (Normalized Prediction Distribution Errors)
    @testset "NPDE" begin
        include("test_npde.jl")
    end

    # Parameter Estimation (NLME)
    @testset "Estimation" begin
        include("test_estimation.jl")
    end

    # Parallel MCMC Chains for SAEM
    @testset "Parallel MCMC" begin
        include("test_parallel_mcmc.jl")
    end

    # Bootstrap for Parameter Uncertainty
    @testset "Bootstrap" begin
        include("test_bootstrap.jl")
    end

    # Dose Escalation Simulation
    @testset "Dose Escalation" begin
        include("test_escalation.jl")
    end

    # TMDD (Target-Mediated Drug Disposition) Models
    @testset "TMDD Models" begin
        include("test_tmdd.jl")
    end

    # MCP-Mod (Phase II Dose-Response)
    @testset "MCP-Mod" begin
        include("test_mcpmod.jl")
    end

    # Exposure-Response Analysis
    @testset "Exposure-Response" begin
        include("test_exposure_response.jl")
    end

    # Probability of Success (Go/No-Go Decisions)
    @testset "Probability of Success" begin
        include("test_pos.jl")
    end

    # Oncology Case Study Integration
    @testset "Oncology Case Study" begin
        include("test_oncology_case_study.jl")
    end

    # Combination Synergy Study (2D Dose-Finding)
    @testset "Combination Synergy" begin
        include("test_combination_synergy.jl")
    end
end
