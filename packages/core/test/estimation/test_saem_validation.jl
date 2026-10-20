# SAEM Validation Tests
# Validates SAEM estimation on simulated datasets with known parameters

using Test
using NeoPKPD
using LinearAlgebra
using StableRNGs
using Statistics

@testset "SAEM Validation" begin

    @testset "SAEM Basic Estimation" begin
        # Create test data with known parameters
        doses = [DoseEvent(0.0, 100.0)]

        subj1 = SubjectData("S1", [0.5, 1.0, 2.0, 4.0, 8.0], [1.8, 1.6, 1.3, 0.9, 0.4], doses)
        subj2 = SubjectData("S2", [0.5, 1.0, 2.0, 4.0, 8.0], [2.2, 1.9, 1.5, 1.0, 0.5], doses)
        subj3 = SubjectData("S3", [0.5, 1.0, 2.0, 4.0, 8.0], [1.9, 1.7, 1.4, 0.95, 0.45], doses)

        observed = ObservedData([subj1, subj2, subj3])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk", model_params, doses)
        grid = SimGrid(0.0, 12.0, collect(0.0:0.5:12.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        # Use reduced iterations for faster testing
        config = EstimationConfig(
            SAEMMethod(
                n_burn=50,
                n_iter=30,
                n_chains=2,
                n_mcmc_steps=20,
                adapt_proposal=true,
                track_diagnostics=true,
                use_all_chains=true
            );
            theta_init=[10.0, 50.0],
            theta_lower=[1.0, 10.0],
            theta_upper=[100.0, 200.0],
            theta_names=[:CL, :V],
            omega_init=diagm([0.09, 0.04]),
            omega_names=[:eta_CL, :eta_V],
            sigma_init=sigma_spec,
            compute_se=false,
            verbose=false,
            seed=UInt64(12345)
        )

        result = NeoPKPD.estimate(observed, model_spec, config; grid=grid, solver=solver)

        @test result isa EstimationResult
        @test length(result.theta) == 2
        @test result.theta[1] > 0  # CL > 0
        @test result.theta[2] > 0  # V > 0
        @test isfinite(result.ofv)
        @test length(result.individuals) == 3
        @test result.n_iterations > 0

        # Check messages contain SAEM info
        @test any(contains.(result.messages, "SAEM"))
        @test any(contains.(result.messages, "acceptance"))
    end

    @testset "SAEM with Diagonal Omega" begin
        doses = [DoseEvent(0.0, 100.0)]

        subj1 = SubjectData("S1", [0.5, 1.0, 2.0, 4.0], [1.8, 1.6, 1.3, 0.9], doses)
        subj2 = SubjectData("S2", [0.5, 1.0, 2.0, 4.0], [2.2, 1.9, 1.5, 1.0], doses)

        observed = ObservedData([subj1, subj2])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk", model_params, doses)
        grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        config = EstimationConfig(
            SAEMMethod(n_burn=30, n_iter=20, n_mcmc_steps=15);
            theta_init=[10.0, 50.0],
            omega_init=diagm([0.09, 0.04]),
            sigma_init=sigma_spec,
            omega_structure=DiagonalOmega(),
            compute_se=false,
            verbose=false,
            seed=UInt64(42)
        )

        result = NeoPKPD.estimate(observed, model_spec, config; grid=grid, solver=solver)

        @test result isa EstimationResult
        @test isfinite(result.ofv)

        # Check omega remains diagonal
        omega_est = result.omega
        for i in 1:size(omega_est, 1)
            for j in 1:size(omega_est, 2)
                if i != j
                    @test omega_est[i, j] ≈ 0.0 atol=1e-10
                end
            end
        end
    end

    @testset "SAEM Reproducibility with Seed" begin
        doses = [DoseEvent(0.0, 100.0)]
        subj1 = SubjectData("S1", [0.5, 1.0, 2.0, 4.0], [1.8, 1.6, 1.3, 0.9], doses)
        subj2 = SubjectData("S2", [0.5, 1.0, 2.0, 4.0], [2.2, 1.9, 1.5, 1.0], doses)

        observed = ObservedData([subj1, subj2])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk", model_params, doses)
        grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        # Run twice with same seed
        config = EstimationConfig(
            SAEMMethod(n_burn=20, n_iter=15, n_mcmc_steps=10);
            theta_init=[10.0, 50.0],
            omega_init=diagm([0.09, 0.04]),
            sigma_init=sigma_spec,
            compute_se=false,
            verbose=false,
            seed=UInt64(99999)
        )

        result1 = NeoPKPD.estimate(observed, model_spec, config; grid=grid, solver=solver)
        result2 = NeoPKPD.estimate(observed, model_spec, config; grid=grid, solver=solver)

        # Same seed should give same results
        @test result1.theta ≈ result2.theta atol=1e-10
        @test result1.ofv ≈ result2.ofv atol=1e-10
        @test result1.omega ≈ result2.omega atol=1e-10
    end

    @testset "SAEM All Chains vs Single Chain" begin
        doses = [DoseEvent(0.0, 100.0)]
        subj1 = SubjectData("S1", [0.5, 1.0, 2.0, 4.0], [1.8, 1.6, 1.3, 0.9], doses)
        subj2 = SubjectData("S2", [0.5, 1.0, 2.0, 4.0], [2.2, 1.9, 1.5, 1.0], doses)

        observed = ObservedData([subj1, subj2])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk", model_params, doses)
        grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        # All chains (default)
        config_all = EstimationConfig(
            SAEMMethod(n_burn=20, n_iter=15, n_mcmc_steps=10, use_all_chains=true);
            theta_init=[10.0, 50.0],
            omega_init=diagm([0.09, 0.04]),
            sigma_init=sigma_spec,
            compute_se=false,
            verbose=false,
            seed=UInt64(111)
        )

        # Single chain
        config_single = EstimationConfig(
            SAEMMethod(n_burn=20, n_iter=15, n_mcmc_steps=10, use_all_chains=false);
            theta_init=[10.0, 50.0],
            omega_init=diagm([0.09, 0.04]),
            sigma_init=sigma_spec,
            compute_se=false,
            verbose=false,
            seed=UInt64(111)
        )

        result_all = NeoPKPD.estimate(observed, model_spec, config_all; grid=grid, solver=solver)
        result_single = NeoPKPD.estimate(observed, model_spec, config_single; grid=grid, solver=solver)

        # Both should produce valid results
        @test isfinite(result_all.ofv)
        @test isfinite(result_single.ofv)
        @test result_all.theta[1] > 0
        @test result_single.theta[1] > 0
    end

    @testset "SAEM Adaptive Proposal" begin
        doses = [DoseEvent(0.0, 100.0)]
        subj1 = SubjectData("S1", [0.5, 1.0, 2.0, 4.0], [1.8, 1.6, 1.3, 0.9], doses)
        subj2 = SubjectData("S2", [0.5, 1.0, 2.0, 4.0], [2.2, 1.9, 1.5, 1.0], doses)

        observed = ObservedData([subj1, subj2])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk", model_params, doses)
        grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        # With adaptation
        config_adapt = EstimationConfig(
            SAEMMethod(
                n_burn=30,
                n_iter=20,
                n_mcmc_steps=15,
                adapt_proposal=true,
                target_acceptance=0.3,
                adaptation_interval=10
            );
            theta_init=[10.0, 50.0],
            omega_init=diagm([0.09, 0.04]),
            sigma_init=sigma_spec,
            compute_se=false,
            verbose=false,
            seed=UInt64(222)
        )

        result = NeoPKPD.estimate(observed, model_spec, config_adapt; grid=grid, solver=solver)

        @test result isa EstimationResult
        @test isfinite(result.ofv)
        @test any(contains.(result.messages, "Adaptive Metropolis-Hastings"))
    end

    @testset "SAEM Step Size Schedules" begin
        doses = [DoseEvent(0.0, 100.0)]
        subj1 = SubjectData("S1", [0.5, 1.0, 2.0], [1.8, 1.6, 1.3], doses)
        subj2 = SubjectData("S2", [0.5, 1.0, 2.0], [2.2, 1.9, 1.5], doses)

        observed = ObservedData([subj1, subj2])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk", model_params, doses)
        grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        # Harmonic schedule
        config_harmonic = EstimationConfig(
            SAEMMethod(n_burn=15, n_iter=10, n_mcmc_steps=10, step_size_schedule=:harmonic);
            theta_init=[10.0, 50.0],
            omega_init=diagm([0.09, 0.04]),
            sigma_init=sigma_spec,
            compute_se=false,
            verbose=false
        )

        # Constant schedule
        config_constant = EstimationConfig(
            SAEMMethod(n_burn=15, n_iter=10, n_mcmc_steps=10, step_size_schedule=:constant);
            theta_init=[10.0, 50.0],
            omega_init=diagm([0.09, 0.04]),
            sigma_init=sigma_spec,
            compute_se=false,
            verbose=false
        )

        result_harmonic = NeoPKPD.estimate(observed, model_spec, config_harmonic; grid=grid, solver=solver)
        result_constant = NeoPKPD.estimate(observed, model_spec, config_constant; grid=grid, solver=solver)

        @test isfinite(result_harmonic.ofv)
        @test isfinite(result_constant.ofv)
    end

end
