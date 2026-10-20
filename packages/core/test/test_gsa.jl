using Test
using NeoPKPD
using Statistics
using StableRNGs

@testset "Global Sensitivity Analysis" begin

    @testset "Sampling Functions" begin
        @testset "Sobol sequence generation" begin
            rng = StableRNG(12345)
            samples = generate_sobol_sequence(100, 3, rng)

            @test size(samples) == (100, 3)
            @test all(0 .<= samples .<= 1)
        end

        @testset "Saltelli sampling" begin
            bounds = ParameterBounds([:CL, :V], [0.5, 5.0], [2.0, 20.0])
            rng = StableRNG(12345)

            samples = generate_saltelli_samples(bounds, 64, rng)

            @test size(samples.A) == (64, 2)
            @test size(samples.B) == (64, 2)
            @test length(samples.AB) == 2
            @test length(samples.BA) == 2

            # Check bounds are respected
            @test all(samples.A[:, 1] .>= 0.5)
            @test all(samples.A[:, 1] .<= 2.0)
            @test all(samples.A[:, 2] .>= 5.0)
            @test all(samples.A[:, 2] .<= 20.0)

            # Check AB matrix construction
            @test samples.AB[1][:, 2] == samples.A[:, 2]  # Column 2 unchanged
            @test samples.AB[1][:, 1] == samples.B[:, 1]  # Column 1 from B
        end

        @testset "Morris trajectories" begin
            bounds = ParameterBounds([:CL, :V], [0.5, 5.0], [2.0, 20.0])
            rng = StableRNG(12345)

            trajectories = generate_morris_trajectories(bounds, 10, 4, 0.667, rng)

            @test length(trajectories) == 10
            @test all(size(t) == (3, 2) for t in trajectories)  # d+1 points, d dims

            # Check bounds are respected
            for traj in trajectories
                @test all(traj[:, 1] .>= 0.5)
                @test all(traj[:, 1] .<= 2.0)
                @test all(traj[:, 2] .>= 5.0)
                @test all(traj[:, 2] .<= 20.0)
            end
        end
    end

    @testset "Sobol' Analysis" begin
        @testset "Type construction" begin
            method = SobolMethod(base_sample_size=128)
            @test method.base_sample_size == 128
            @test method.compute_second_order == false
            @test method.bootstrap_samples == 1000
            @test method.bootstrap_ci_level == 0.95

            method2 = SobolMethod(base_sample_size=256, bootstrap_samples=500)
            @test method2.base_sample_size == 256
            @test method2.bootstrap_samples == 500

            # Validation
            @test_throws Exception SobolMethod(base_sample_size=32)  # Too small
            @test_throws Exception SobolMethod(bootstrap_ci_level=1.5)  # Invalid CI
        end

        @testset "Parameter bounds" begin
            bounds = ParameterBounds([:CL, :V, :ka], [0.5, 5.0, 0.1], [2.0, 20.0, 2.0])
            @test length(bounds) == 3
            @test bounds.params == [:CL, :V, :ka]

            # Dict constructor
            bounds2 = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
            @test length(bounds2) == 2
        end

        @testset "Basic Sobol' analysis" begin
            # Simple one-compartment model
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
            gsa_spec = GlobalSensitivitySpec(
                SobolMethod(base_sample_size=64, bootstrap_samples=0),
                bounds
            )

            result = run_sobol_sensitivity(spec, grid, solver, gsa_spec)

            @test result isa SobolResult
            @test length(result.indices) == 2
            @test :CL in keys(result.indices)
            @test :V in keys(result.indices)

            # Indices should be in valid range
            for (param, idx) in result.indices
                @test 0.0 <= idx.Si <= 1.0
                @test 0.0 <= idx.STi <= 1.0
                # Total should be >= first-order
                @test idx.STi >= idx.Si - 0.1  # Allow small numerical tolerance
            end

            # Convergence metric should be reasonable
            @test 0.0 <= result.convergence_metric <= 2.0

            # Check metadata
            @test result.metadata["method"] == "Sobol"
            @test result.metadata["base_sample_size"] == 64
        end

        @testset "Sobol' with bootstrap CI" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:2.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
            gsa_spec = GlobalSensitivitySpec(
                SobolMethod(base_sample_size=64, bootstrap_samples=100, bootstrap_ci_level=0.90),
                bounds
            )

            result = run_sobol_sensitivity(spec, grid, solver, gsa_spec)

            # CI bounds should be ordered correctly
            for (param, idx) in result.indices
                @test idx.Si_ci_lower <= idx.Si <= idx.Si_ci_upper
                @test idx.STi_ci_lower <= idx.STi <= idx.STi_ci_upper
            end
        end
    end

    @testset "Morris Analysis" begin
        @testset "Type construction" begin
            method = MorrisMethod(n_trajectories=20)
            @test method.n_trajectories == 20
            @test method.n_levels == 4
            @test method.delta > 0

            method2 = MorrisMethod(n_trajectories=30, n_levels=6)
            @test method2.n_trajectories == 30
            @test method2.n_levels == 6

            # Validation
            @test_throws Exception MorrisMethod(n_trajectories=2)  # Too few
            @test_throws Exception MorrisMethod(n_levels=1)  # Too few
        end

        @testset "Basic Morris analysis" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
            gsa_spec = GlobalSensitivitySpec(
                MorrisMethod(n_trajectories=10),
                bounds
            )

            result = run_morris_sensitivity(spec, grid, solver, gsa_spec)

            @test result isa MorrisResult
            @test length(result.indices) == 2
            @test :CL in keys(result.indices)
            @test :V in keys(result.indices)

            # Morris indices should be computed
            for (param, idx) in result.indices
                @test idx.mu_star >= 0  # Absolute mean is non-negative
                @test idx.sigma >= 0    # Std dev is non-negative
            end

            # Elementary effects should be stored
            @test length(result.elementary_effects) == 2
            for (param, ees) in result.elementary_effects
                @test length(ees) >= 1
            end

            # Check metadata
            @test result.metadata["method"] == "Morris"
            @test result.metadata["n_trajectories"] == 10
        end

        @testset "Morris ranking" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
            gsa_spec = GlobalSensitivitySpec(
                MorrisMethod(n_trajectories=15),
                bounds
            )

            result = run_morris_sensitivity(spec, grid, solver, gsa_spec)

            # Test ranking function
            rankings = rank_parameters(result; by=:mu_star)
            @test length(rankings) == 2
            @test rankings[1][2] >= rankings[2][2]  # Descending order

            # Test importance identification
            important = identify_important_parameters(result; threshold=0.1)
            @test length(important) >= 1
        end
    end

    @testset "Serialization" begin
        @testset "Sobol' serialization" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:2.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
            gsa_spec = GlobalSensitivitySpec(
                SobolMethod(base_sample_size=64, bootstrap_samples=0),
                bounds
            )

            result = run_sobol_sensitivity(spec, grid, solver, gsa_spec)

            # Serialize
            artifact = serialize_sobol_result(result)

            @test artifact["artifact_type"] == "sobol_sensitivity"
            @test haskey(artifact, "indices")
            @test haskey(artifact, "n_evaluations")
            @test haskey(artifact, "convergence_metric")

            # Deserialize
            result2 = deserialize_sobol_result(artifact)

            @test result2 isa SobolResult
            @test result2.n_evaluations == result.n_evaluations
            @test length(result2.indices) == length(result.indices)
        end

        @testset "Morris serialization" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:2.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
            gsa_spec = GlobalSensitivitySpec(
                MorrisMethod(n_trajectories=10),
                bounds
            )

            result = run_morris_sensitivity(spec, grid, solver, gsa_spec)

            # Serialize
            artifact = serialize_morris_result(result)

            @test artifact["artifact_type"] == "morris_sensitivity"
            @test haskey(artifact, "indices")
            @test haskey(artifact, "elementary_effects")

            # Deserialize
            result2 = deserialize_morris_result(artifact)

            @test result2 isa MorrisResult
            @test result2.n_evaluations == result.n_evaluations
            @test length(result2.indices) == length(result.indices)
        end
    end

    @testset "Index Computation" begin
        @testset "Sobol' point estimates" begin
            # Test with known values
            rng = StableRNG(12345)
            N = 100

            # Create synthetic outputs
            Y_A = randn(rng, N)
            Y_B = randn(rng, N)
            Y_AB = [randn(rng, N) for _ in 1:2]

            indices = compute_sobol_indices(Y_A, Y_B, Y_AB; bootstrap_samples=0)

            @test length(indices) == 2
            for (j, idx) in indices
                @test 0.0 <= idx.Si <= 1.0 || idx.Si < 0.1  # Allow small negative due to estimation
                @test 0.0 <= idx.STi <= 1.0 || idx.STi < 0.1
            end
        end

        @testset "Morris indices from EEs" begin
            # Create synthetic elementary effects
            ees = Dict(
                :CL => [1.0, 2.0, 1.5, 0.5, 2.5],
                :V => [-0.5, 0.5, -0.3, 0.2, -0.4]
            )

            indices = compute_morris_indices(ees)

            @test indices[:CL].mu ≈ mean(ees[:CL])
            @test indices[:CL].mu_star ≈ mean(abs.(ees[:CL]))
            @test indices[:CL].sigma ≈ std(ees[:CL])

            @test indices[:V].mu ≈ mean(ees[:V])
            @test indices[:V].mu_star ≈ mean(abs.(ees[:V]))
        end
    end

    @testset "Reproducibility" begin
        @testset "Sobol' reproducibility with seed" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:2.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
            gsa_spec1 = GlobalSensitivitySpec(
                SobolMethod(base_sample_size=64, bootstrap_samples=0),
                bounds;
                seed=UInt64(42)
            )
            gsa_spec2 = GlobalSensitivitySpec(
                SobolMethod(base_sample_size=64, bootstrap_samples=0),
                bounds;
                seed=UInt64(42)
            )

            result1 = run_sobol_sensitivity(spec, grid, solver, gsa_spec1)
            result2 = run_sobol_sensitivity(spec, grid, solver, gsa_spec2)

            @test result1.indices[:CL].Si ≈ result2.indices[:CL].Si
            @test result1.indices[:V].Si ≈ result2.indices[:V].Si
        end

        @testset "Morris reproducibility with seed" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:2.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
            gsa_spec1 = GlobalSensitivitySpec(
                MorrisMethod(n_trajectories=10),
                bounds;
                seed=UInt64(42)
            )
            gsa_spec2 = GlobalSensitivitySpec(
                MorrisMethod(n_trajectories=10),
                bounds;
                seed=UInt64(42)
            )

            result1 = run_morris_sensitivity(spec, grid, solver, gsa_spec1)
            result2 = run_morris_sensitivity(spec, grid, solver, gsa_spec2)

            @test result1.indices[:CL].mu_star ≈ result2.indices[:CL].mu_star
            @test result1.indices[:V].mu_star ≈ result2.indices[:V].mu_star
        end
    end

end
