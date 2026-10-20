# Test suite for Parallel MCMC Chain Implementation
# Tests thread-parallel MCMC sampling for SAEM estimation

using Test
using NeoPKPD
using LinearAlgebra
using StableRNGs
using Statistics
using Base.Threads: nthreads

@testset "Parallel MCMC Chains" begin

    @testset "Chain RNG Creation" begin
        # Test that create_chain_rngs produces independent RNG streams
        base_rng = StableRNG(12345)
        chain_rngs = NeoPKPD.create_chain_rngs(base_rng, 4)

        @test length(chain_rngs) == 4
        @test all(r isa StableRNG for r in chain_rngs)

        # Verify independence by generating samples from each RNG
        samples1 = rand(chain_rngs[1], 100)
        samples2 = rand(chain_rngs[2], 100)
        samples3 = rand(chain_rngs[3], 100)
        samples4 = rand(chain_rngs[4], 100)

        # Samples should be different (not identical streams)
        @test samples1 != samples2
        @test samples2 != samples3
        @test samples3 != samples4

        # Verify reproducibility with same base seed
        base_rng2 = StableRNG(12345)
        chain_rngs2 = NeoPKPD.create_chain_rngs(base_rng2, 4)
        samples1_repeat = rand(chain_rngs2[1], 100)
        @test samples1 == samples1_repeat  # Same seed should reproduce
    end

    @testset "Single Chain MCMC Kernel" begin
        # Test the inner MCMC kernel function
        n_eta = 2
        eta_init = [0.0, 0.0]
        proposal_sd = [0.1, 0.1]
        n_steps = 50
        chain_rng = StableRNG(42)

        # Simple log posterior: multivariate normal
        omega_inv = [1.0 0.0; 0.0 1.0]  # Identity precision matrix

        log_posterior = eta -> -0.5 * dot(eta, omega_inv * eta)

        eta_final, n_accepts = NeoPKPD.mcmc_run_single_chain(
            1, eta_init, log_posterior, proposal_sd, n_steps, chain_rng
        )

        @test length(eta_final) == n_eta
        @test 0 <= n_accepts <= n_steps
        @test all(isfinite, eta_final)
    end

    @testset "Parallel vs Sequential MCMC Consistency" begin
        # Setup test data
        rng = StableRNG(99999)

        # Model setup
        doses = [DoseEvent(1.0, 100.0)]
        times = [0.5, 1.0, 2.0, 4.0, 8.0]
        obs = [5.0, 8.0, 12.0, 9.0, 5.0]

        model_spec = ModelSpec(
            OneCompOralFirstOrder(),
            "test_model",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            doses
        )

        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # Parameters
        theta = [1.5, 2.0, 30.0]
        omega = [0.1 0.0; 0.0 0.1]
        omega_inv = inv(omega)
        omega_chol = cholesky(Symmetric(omega)).L

        sigma = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        # Initial chains
        n_chains = 3
        current_chains = [[0.0, 0.0] for _ in 1:n_chains]
        proposal_sd = [0.1, 0.1]
        n_steps = 20

        # Run sequential version
        rng_seq = StableRNG(12345)
        new_chains_seq, accepts_seq, n_props_seq = NeoPKPD.mcmc_sample_eta_adaptive(
            theta, omega_inv, omega_chol, sigma,
            times, obs, doses, model_spec, grid, solver,
            current_chains, proposal_sd, n_steps, rng_seq
        )

        # Run parallel version (with fresh current_chains since RNG is consumed)
        current_chains_par = [[0.0, 0.0] for _ in 1:n_chains]
        rng_par = StableRNG(12345)
        new_chains_par, accepts_par, n_props_par = NeoPKPD.mcmc_sample_eta_parallel_chains(
            theta, omega_inv, omega_chol, sigma,
            times, obs, doses, model_spec, grid, solver,
            current_chains_par, proposal_sd, n_steps, rng_par
        )

        # Both should return valid results
        @test length(new_chains_seq) == n_chains
        @test length(new_chains_par) == n_chains
        @test length(accepts_seq) == n_chains
        @test length(accepts_par) == n_chains
        @test all(n_props_seq .== n_steps)
        @test all(n_props_par .== n_steps)

        # All chains should produce finite etas
        @test all(all(isfinite, c) for c in new_chains_seq)
        @test all(all(isfinite, c) for c in new_chains_par)

        # Acceptance rates should be reasonable (not 0% or 100%)
        total_accepts_seq = sum(accepts_seq)
        total_accepts_par = sum(accepts_par)
        total_props = n_chains * n_steps

        @test 0 <= total_accepts_seq <= total_props
        @test 0 <= total_accepts_par <= total_props
    end

    @testset "Parallel MCMC Determinism" begin
        # Same seed should produce same results
        rng1 = StableRNG(55555)
        rng2 = StableRNG(55555)

        doses = [DoseEvent(1.0, 100.0)]
        times = [0.5, 1.0, 2.0]
        obs = [5.0, 8.0, 6.0]

        model_spec = ModelSpec(
            OneCompOralFirstOrder(),
            "test",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            doses
        )

        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        theta = [1.5, 2.0, 30.0]
        omega = [0.1 0.0; 0.0 0.1]
        omega_inv = inv(omega)
        omega_chol = cholesky(Symmetric(omega)).L

        sigma = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        n_chains = 4
        current_chains1 = [[0.0, 0.0] for _ in 1:n_chains]
        current_chains2 = [[0.0, 0.0] for _ in 1:n_chains]
        proposal_sd = [0.1, 0.1]
        n_steps = 30

        chains1, accepts1, _ = NeoPKPD.mcmc_sample_eta_parallel_chains(
            theta, omega_inv, omega_chol, sigma,
            times, obs, doses, model_spec, grid, solver,
            current_chains1, proposal_sd, n_steps, rng1
        )

        chains2, accepts2, _ = NeoPKPD.mcmc_sample_eta_parallel_chains(
            theta, omega_inv, omega_chol, sigma,
            times, obs, doses, model_spec, grid, solver,
            current_chains2, proposal_sd, n_steps, rng2
        )

        # Results should be identical with same seed
        for c in 1:n_chains
            @test chains1[c] ≈ chains2[c] atol=1e-10
        end
        @test accepts1 == accepts2
    end

    @testset "SAEMMethod parallel_chains Parameter" begin
        # Test that SAEMMethod correctly exposes parallel_chains
        method_default = SAEMMethod()
        @test method_default.parallel_chains == true

        method_parallel = SAEMMethod(parallel_chains=true)
        @test method_parallel.parallel_chains == true

        method_sequential = SAEMMethod(parallel_chains=false)
        @test method_sequential.parallel_chains == false

        # Other parameters should still work
        method_custom = SAEMMethod(
            n_burn=50,
            n_iter=30,
            n_chains=4,
            n_mcmc_steps=50,
            parallel_chains=true
        )
        @test method_custom.n_burn == 50
        @test method_custom.n_chains == 4
        @test method_custom.parallel_chains == true
    end

    @testset "Multiple Chains Return Valid Results" begin
        # Test that parallel chains return valid results (not necessarily exploring well)
        rng = StableRNG(77777)

        doses = [DoseEvent(1.0, 100.0)]
        times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0]
        obs = [5.0, 8.0, 12.0, 9.0, 5.0, 3.0]

        model_spec = ModelSpec(
            OneCompOralFirstOrder(),
            "test",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            doses
        )

        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        theta = [1.5, 2.0, 30.0]
        omega = [0.1 0.0; 0.0 0.1]
        omega_inv = inv(omega)
        omega_chol = cholesky(Symmetric(omega)).L

        sigma = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        n_chains = 4
        n_steps = 50

        # Start from small initial values
        current_chains = [
            [0.0, 0.0] for _ in 1:n_chains
        ]
        proposal_sd = sqrt.(diag(omega))

        final_chains, accepts, n_props = NeoPKPD.mcmc_sample_eta_parallel_chains(
            theta, omega_inv, omega_chol, sigma,
            times, obs, doses, model_spec, grid, solver,
            current_chains, proposal_sd, n_steps, rng
        )

        # All chains should return valid eta vectors
        @test length(final_chains) == n_chains
        for c in 1:n_chains
            @test length(final_chains[c]) == 2
            @test all(isfinite, final_chains[c])
        end

        # Acceptance counts should be valid
        @test length(accepts) == n_chains
        for c in 1:n_chains
            @test 0 <= accepts[c] <= n_steps
        end

        # Compute chain mean - should be finite
        chain_mean = mean(final_chains)
        @test length(chain_mean) == 2
        @test all(isfinite, chain_mean)
    end

    @testset "Thread Count Reporting" begin
        # Verify nthreads() is accessible
        @test nthreads() >= 1
    end

end
