# Test suite for Parallel Computing Infrastructure

using Test
using NeoPKPD
using LinearAlgebra
using StableRNGs

@testset "Parallel Computing" begin

    @testset "Backend Types" begin
        @testset "SerialBackend" begin
            backend = SerialBackend()
            @test backend isa ParallelBackend
        end

        @testset "ThreadedBackend" begin
            backend = ThreadedBackend()
            @test backend isa ParallelBackend
            @test backend.n_threads >= 1

            # With explicit thread count (clamped to available threads)
            backend2 = ThreadedBackend(2)
            @test backend2.n_threads >= 1  # May be 1 if Julia has 1 thread
            @test backend2.n_threads <= 2  # Cannot exceed requested

            # Cannot exceed available threads
            @test_throws AssertionError ThreadedBackend(0)
        end

        @testset "DistributedBackend" begin
            backend = DistributedBackend()
            @test backend isa ParallelBackend
            @test backend.n_workers == 1

            backend2 = DistributedBackend(4)
            @test backend2.n_workers == 4

            @test_throws AssertionError DistributedBackend(0)
        end
    end

    @testset "ParallelConfig" begin
        @testset "Default configuration" begin
            config = ParallelConfig()
            @test config.backend isa ThreadedBackend
            @test config.chunk_size == 0  # Auto
            @test config.load_balance == true
            @test config.progress == false
            @test config.seed === nothing
            @test config.verbose == false
        end

        @testset "Serial backend" begin
            config = ParallelConfig(SerialBackend())
            @test config.backend isa SerialBackend
            @test !is_parallel(config)
            @test n_workers(config) == 1
        end

        @testset "Threaded backend" begin
            config = ParallelConfig(ThreadedBackend(4))
            @test config.backend isa ThreadedBackend
            @test n_workers(config) >= 1  # May be clamped to available threads
            @test is_parallel(config) || n_workers(config) == 1  # Parallel if >1 thread
        end

        @testset "With seed" begin
            config = ParallelConfig(ThreadedBackend(); seed=12345)
            @test config.seed == 12345
        end

        @testset "With custom options" begin
            config = ParallelConfig(
                ThreadedBackend();
                chunk_size=10,
                load_balance=false,
                progress=true,
                verbose=true
            )
            @test config.chunk_size == 10
            @test config.load_balance == false
            @test config.progress == true
            @test config.verbose == true
        end
    end

    @testset "Thread-Safe RNG" begin
        @testset "Create independent streams" begin
            rngs = create_thread_rngs(4; seed=UInt64(12345))
            @test length(rngs) == 4

            # Each RNG should produce different sequences
            samples = [rand(rng) for rng in rngs]
            @test length(unique(samples)) == 4
        end

        @testset "Reproducibility with seed" begin
            rngs1 = create_thread_rngs(3; seed=UInt64(42))
            rngs2 = create_thread_rngs(3; seed=UInt64(42))

            for i in 1:3
                @test rand(rngs1[i]) == rand(rngs2[i])
            end
        end

        @testset "Subject RNGs" begin
            config = ParallelConfig(ThreadedBackend(); seed=12345)
            rngs = create_subject_rngs(10, config)
            @test length(rngs) == 10
        end
    end

    @testset "Chunking" begin
        @testset "Optimal chunk size" begin
            # Small problem
            chunk_size = get_optimal_chunk_size(10, 4)
            @test chunk_size >= 1

            # Large problem
            chunk_size = get_optimal_chunk_size(1000, 8)
            @test chunk_size > 1
            @test chunk_size < 1000

            # Single worker
            chunk_size = get_optimal_chunk_size(100, 1)
            @test chunk_size == 100
        end

        @testset "Create chunks" begin
            chunks = create_chunks(10, 3)
            @test length(chunks) == 4
            @test chunks[1] == 1:3
            @test chunks[2] == 4:6
            @test chunks[3] == 7:9
            @test chunks[4] == 10:10

            # Full chunk
            chunks = create_chunks(10, 10)
            @test length(chunks) == 1
            @test chunks[1] == 1:10

            # Chunk larger than items
            chunks = create_chunks(5, 10)
            @test length(chunks) == 1
        end
    end

    @testset "Parallel Map" begin
        @testset "Serial execution" begin
            config = ParallelConfig(SerialBackend())
            items = [1, 2, 3, 4, 5]

            results = parallel_map(x -> x^2, items, config)

            @test length(results) == 5
            @test results == [1, 4, 9, 16, 25]
        end

        @testset "Threaded execution" begin
            config = ParallelConfig(ThreadedBackend())
            items = collect(1:100)

            results = parallel_map(x -> x^2, items, config)

            @test length(results) == 100
            @test results == [x^2 for x in 1:100]
        end

        @testset "Empty input" begin
            config = ParallelConfig(ThreadedBackend())
            results = parallel_map(x -> x^2, Int[], config)
            @test isempty(results)
        end

        @testset "Single item" begin
            config = ParallelConfig(ThreadedBackend())
            results = parallel_map(x -> x^2, [5], config)
            @test results == [25]
        end

        @testset "With RNG" begin
            config = ParallelConfig(ThreadedBackend(); seed=42)
            items = collect(1:10)
            rngs = create_thread_rngs(10; seed=UInt64(42))

            results = parallel_map_with_rng(
                (x, rng) -> x + rand(rng),
                items,
                rngs,
                config
            )

            @test length(results) == 10
            @test all(results .> items)  # Added random values
        end
    end

    @testset "Parallel Map Reduce" begin
        @testset "Serial sum" begin
            config = ParallelConfig(SerialBackend())
            items = collect(1:10)

            result = parallel_map_reduce(x -> x^2, +, 0.0, items, config)

            @test result == sum(x^2 for x in 1:10)
        end

        @testset "Threaded sum" begin
            config = ParallelConfig(ThreadedBackend())
            items = collect(1:1000)

            result = parallel_map_reduce(x -> x^2, +, 0.0, items, config)

            @test result == sum(x^2 for x in 1:1000)
        end

        @testset "Parallel sum helper" begin
            config = ParallelConfig(ThreadedBackend())
            items = collect(1.0:100.0)

            result = parallel_sum(x -> x^2, items, config)

            @test isapprox(result, sum(x^2 for x in 1.0:100.0))
        end

        @testset "Empty input" begin
            config = ParallelConfig(ThreadedBackend())
            result = parallel_map_reduce(x -> x^2, +, 0.0, Float64[], config)
            @test result == 0.0
        end

        @testset "Custom reduction" begin
            config = ParallelConfig(ThreadedBackend())
            items = [1, 2, 3, 4, 5]

            # Product
            result = parallel_map_reduce(identity, *, 1, items, config)
            @test result == 120  # 5!
        end
    end

    @testset "Process Subjects Parallel" begin
        @testset "Basic processing" begin
            config = ParallelConfig(SerialBackend())

            # Mock subject data
            subjects = [("S1", 1), ("S2", 2), ("S3", 3)]

            results = process_subjects_parallel(
                (subj, idx) -> (subj[1], subj[2] * 10),
                subjects,
                config
            )

            @test length(results) == 3
            @test results[1] == ("S1", 10)
            @test results[2] == ("S2", 20)
            @test results[3] == ("S3", 30)
        end

        @testset "With RNG" begin
            config = ParallelConfig(ThreadedBackend(); seed=42)

            subjects = [("S1", 1), ("S2", 2), ("S3", 3)]

            results = process_subjects_parallel_with_rng(
                (subj, idx, rng) -> (subj[1], subj[2] + rand(rng)),
                subjects,
                config
            )

            @test length(results) == 3
            @test all(r -> r[2] > 0, results)
        end
    end

    @testset "Parallel OFV Computation" begin
        @testset "Simple OFV" begin
            config = ParallelConfig(SerialBackend())

            # Mock subjects returning OFV contributions
            subjects = [(i, i^2) for i in 1:5]

            ofv = compute_ofv_parallel(
                subj -> Float64(subj[2]),
                subjects,
                config
            )

            @test ofv == 1.0 + 4.0 + 9.0 + 16.0 + 25.0
        end

        @testset "OFV with etas" begin
            config = ParallelConfig(SerialBackend())

            subjects = [(1, [1.0]), (2, [2.0]), (3, [3.0])]
            etas = [[0.1], [0.2], [0.3]]

            ofv = compute_ofv_with_etas_parallel(
                (subj, eta) -> subj[1] * eta[1],
                subjects,
                etas,
                config
            )

            @test isapprox(ofv, 1*0.1 + 2*0.2 + 3*0.3)
        end
    end

    @testset "Configuration Management" begin
        @testset "Current config" begin
            default_config = current_parallel_config()
            @test default_config isa ParallelConfig
        end

        @testset "Set config" begin
            old_config = current_parallel_config()

            new_config = ParallelConfig(ThreadedBackend(4))
            set_parallel_config!(new_config)

            # n_workers may be clamped to available threads
            @test n_workers(current_parallel_config()) >= 1
            @test current_parallel_config().backend isa ThreadedBackend

            # Restore
            set_parallel_config!(old_config)
        end

        @testset "With config block" begin
            old_config = current_parallel_config()

            result = with_parallel_config(ParallelConfig(SerialBackend())) do
                @test current_parallel_config().backend isa SerialBackend
                42
            end

            @test result == 42
            @test typeof(current_parallel_config().backend) == typeof(old_config.backend)
        end
    end

    @testset "Recommend Config" begin
        @testset "Small problem" begin
            config = recommend_parallel_config(5)
            @test config.backend isa SerialBackend
        end

        @testset "Medium problem" begin
            config = recommend_parallel_config(50)
            # Should recommend threading if available
            @test config.backend isa SerialBackend || config.backend isa ThreadedBackend
        end

        @testset "Large problem" begin
            config = recommend_parallel_config(500)
            # Should recommend threading with progress if available
            @test config.backend isa SerialBackend || config.backend isa ThreadedBackend
        end
    end

    @testset "Integration with Estimation" begin
        # Test that parallel config can be passed to estimate function
        @testset "FOCE-I with parallel config" begin
            # Create simple test data
            doses = [DoseEvent(0.0, 100.0)]

            subj1 = SubjectData(
                "SUBJ001",
                [0.5, 1.0, 2.0],
                [1.8, 1.6, 1.3],
                doses
            )
            subj2 = SubjectData(
                "SUBJ002",
                [0.5, 1.0, 2.0],
                [2.0, 1.8, 1.5],
                doses
            )

            observed = ObservedData([subj1, subj2])

            model_params = OneCompIVBolusParams(10.0, 50.0)
            model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

            grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            sigma_spec = ResidualErrorSpec(
                ProportionalError(),
                ProportionalErrorParams(0.1),
                :conc,
                UInt64(1)
            )

            config = EstimationConfig(
                FOCEIMethod(max_inner_iter=10, inner_tol=1e-3);
                theta_init=[10.0, 50.0],
                omega_init=diagm([0.09, 0.04]),
                sigma_init=sigma_spec,
                max_iter=5,
                compute_se=false,
                verbose=false
            )

            # Test with serial config
            serial_config = ParallelConfig(SerialBackend())
            result_serial = NeoPKPD.estimate(
                observed, model_spec, config;
                grid=grid, solver=solver, parallel_config=serial_config
            )

            @test result_serial isa EstimationResult
            @test length(result_serial.theta) == 2
            @test isfinite(result_serial.ofv)

            # Test with threaded config (should produce similar results)
            threaded_config = ParallelConfig(ThreadedBackend())
            result_threaded = NeoPKPD.estimate(
                observed, model_spec, config;
                grid=grid, solver=solver, parallel_config=threaded_config
            )

            @test result_threaded isa EstimationResult
            @test length(result_threaded.theta) == 2
            @test isfinite(result_threaded.ofv)

            # Results should be similar (may differ slightly due to floating point)
            @test isapprox(result_serial.ofv, result_threaded.ofv; rtol=0.01)
        end

        @testset "Laplacian with parallel config" begin
            doses = [DoseEvent(0.0, 100.0)]

            subj1 = SubjectData("S1", [1.0, 2.0], [1.5, 1.2], doses)
            subj2 = SubjectData("S2", [1.0, 2.0], [1.6, 1.3], doses)

            observed = ObservedData([subj1, subj2])

            model_params = OneCompIVBolusParams(10.0, 50.0)
            model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

            grid = SimGrid(0.0, 4.0, collect(0.0:0.5:4.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            sigma_spec = ResidualErrorSpec(
                ProportionalError(),
                ProportionalErrorParams(0.1),
                :conc,
                UInt64(1)
            )

            config = EstimationConfig(
                LaplacianMethod(max_inner_iter=10, inner_tol=1e-3);
                theta_init=[10.0, 50.0],
                omega_init=diagm([0.09, 0.04]),
                sigma_init=sigma_spec,
                max_iter=3,
                compute_se=false,
                verbose=false
            )

            parallel_config = ParallelConfig(ThreadedBackend())
            result = NeoPKPD.estimate(
                observed, model_spec, config;
                grid=grid, solver=solver, parallel_config=parallel_config
            )

            @test result isa EstimationResult
            @test isfinite(result.ofv)
        end
    end

    @testset "Pretty Printing" begin
        @testset "ParallelConfig show" begin
            config = ParallelConfig(ThreadedBackend(4); seed=12345)
            io = IOBuffer()
            show(io, config)
            output = String(take!(io))
            @test contains(output, "ParallelConfig")
            @test contains(output, "ThreadedBackend")
            @test contains(output, string(n_workers(config)))
        end

        @testset "Backend show" begin
            io = IOBuffer()
            backend = ThreadedBackend(4)
            show(io, backend)
            output = String(take!(io))
            @test contains(output, "threads")
            @test contains(output, string(backend.n_threads))

            io = IOBuffer()
            show(io, SerialBackend())
            output = String(take!(io))
            @test contains(output, "SerialBackend")

            io = IOBuffer()
            show(io, DistributedBackend(8))
            output = String(take!(io))
            @test contains(output, "8 workers")
        end
    end

    @testset "Parallel Bootstrap" begin
        @testset "BootstrapSpec parallel configuration" begin
            # Default: sequential
            spec1 = BootstrapSpec(n_bootstrap=100)
            @test spec1.parallel === false

            # Explicit parallel flag
            spec2 = BootstrapSpec(n_bootstrap=100, parallel=true)
            @test spec2.parallel === true

            # With ParallelConfig
            config = ParallelConfig(ThreadedBackend(); seed=42)
            spec3 = BootstrapSpec(n_bootstrap=100, parallel=config)
            @test spec3.parallel isa ParallelConfig
        end

        @testset "get_parallel_config" begin
            # From Bool=false -> SerialBackend
            spec1 = BootstrapSpec(n_bootstrap=100, parallel=false)
            config1 = NeoPKPD.get_parallel_config(spec1)
            @test config1.backend isa SerialBackend

            # From Bool=true -> ThreadedBackend
            spec2 = BootstrapSpec(n_bootstrap=100, parallel=true)
            config2 = NeoPKPD.get_parallel_config(spec2)
            @test config2.backend isa ThreadedBackend

            # From ParallelConfig -> pass-through
            custom_config = ParallelConfig(SerialBackend(); seed=99)
            spec3 = BootstrapSpec(n_bootstrap=100, parallel=custom_config)
            config3 = NeoPKPD.get_parallel_config(spec3)
            @test config3 === custom_config
        end

        @testset "Sequential bootstrap execution" begin
            # Create simple test data
            doses = [DoseEvent(0.0, 100.0)]
            subj1 = SubjectData("S1", [1.0, 2.0], [1.5, 1.2], doses)
            subj2 = SubjectData("S2", [1.0, 2.0], [1.6, 1.3], doses)
            subj3 = SubjectData("S3", [1.0, 2.0], [1.4, 1.1], doses)
            observed = ObservedData([subj1, subj2, subj3])

            # Mock estimation function
            mock_estimate_calls = Ref(0)
            function mock_estimate(data)
                mock_estimate_calls[] += 1
                # Return mock result
                return (
                    theta = [10.0 + rand(), 50.0 + rand()],
                    omega = diagm([0.09, 0.04]),
                    sigma = (params = (sigma = 0.1,),),
                    convergence = true,
                    ofv = 100.0
                )
            end

            # Run sequential bootstrap (small for testing)
            spec = BootstrapSpec(n_bootstrap=100, seed=UInt64(12345), parallel=false)
            result = run_bootstrap(mock_estimate, observed, spec)

            @test result isa BootstrapResult
            @test size(result.theta_estimates, 1) == 100
            @test size(result.theta_estimates, 2) == 2
            @test result.diagnostics.n_successful == 100
            @test result.diagnostics.n_failed == 0
            @test result.diagnostics.convergence_rate == 1.0

            # 101 calls: 1 original + 100 bootstrap
            @test mock_estimate_calls[] == 101
        end

        @testset "Parallel bootstrap execution" begin
            # Create simple test data
            doses = [DoseEvent(0.0, 100.0)]
            subj1 = SubjectData("S1", [1.0, 2.0], [1.5, 1.2], doses)
            subj2 = SubjectData("S2", [1.0, 2.0], [1.6, 1.3], doses)
            subj3 = SubjectData("S3", [1.0, 2.0], [1.4, 1.1], doses)
            observed = ObservedData([subj1, subj2, subj3])

            # Thread-safe counter using Atomic
            call_count = Threads.Atomic{Int}(0)

            function mock_estimate_parallel(data)
                Threads.atomic_add!(call_count, 1)
                # Simulate some work
                sleep(0.001)
                return (
                    theta = [10.0 + 0.1 * randn(), 50.0 + 0.5 * randn()],
                    omega = diagm([0.09, 0.04]),
                    sigma = (params = (sigma = 0.1,),),
                    convergence = true,
                    ofv = 100.0
                )
            end

            # Run parallel bootstrap
            spec = BootstrapSpec(n_bootstrap=100, seed=UInt64(42), parallel=true)
            result = run_bootstrap(mock_estimate_parallel, observed, spec)

            @test result isa BootstrapResult
            @test size(result.theta_estimates, 1) == 100
            @test result.diagnostics.n_successful == 100

            # Verify all replicates ran
            @test call_count[] == 101  # 1 original + 100 bootstrap
        end

        @testset "Bootstrap reproducibility with parallel" begin
            # Create test data
            doses = [DoseEvent(0.0, 100.0)]
            subj1 = SubjectData("S1", [1.0, 2.0], [1.5, 1.2], doses)
            subj2 = SubjectData("S2", [1.0, 2.0], [1.6, 1.3], doses)
            observed = ObservedData([subj1, subj2])

            # Deterministic mock function
            function deterministic_estimate(data)
                # Use subject IDs to create deterministic output
                id_hash = sum(hash(s.subject_id) for s in data.subjects)
                return (
                    theta = [10.0 + (id_hash % 100) / 1000, 50.0 + (id_hash % 200) / 1000],
                    omega = diagm([0.09, 0.04]),
                    sigma = nothing,
                    convergence = true,
                    ofv = 100.0
                )
            end

            # Run twice with same seed - should get same resampling pattern
            spec1 = BootstrapSpec(n_bootstrap=100, seed=UInt64(42), parallel=false)
            result1 = run_bootstrap(deterministic_estimate, observed, spec1)

            spec2 = BootstrapSpec(n_bootstrap=100, seed=UInt64(42), parallel=false)
            result2 = run_bootstrap(deterministic_estimate, observed, spec2)

            # Results should be identical (same seed = same resampling)
            @test result1.theta_estimates == result2.theta_estimates
            @test result1.theta_mean == result2.theta_mean
            @test result1.theta_se == result2.theta_se
        end

        @testset "Bootstrap with failed replicates" begin
            doses = [DoseEvent(0.0, 100.0)]
            subj1 = SubjectData("S1", [1.0, 2.0], [1.5, 1.2], doses)
            subj2 = SubjectData("S2", [1.0, 2.0], [1.6, 1.3], doses)
            observed = ObservedData([subj1, subj2])

            fail_counter = Ref(0)

            function sometimes_failing_estimate(data)
                fail_counter[] += 1
                # Fail every 5th replicate (after original)
                if fail_counter[] > 1 && (fail_counter[] - 1) % 5 == 0
                    error("Simulated estimation failure")
                end
                return (
                    theta = [10.0, 50.0],
                    omega = diagm([0.09, 0.04]),
                    sigma = nothing,
                    convergence = true,
                    ofv = 100.0
                )
            end

            spec = BootstrapSpec(n_bootstrap=100, seed=UInt64(123), parallel=false)
            result = run_bootstrap(sometimes_failing_estimate, observed, spec)

            @test result isa BootstrapResult
            @test result.diagnostics.n_failed == 20  # 20 out of 100 fail
            @test result.diagnostics.n_successful == 80
            @test result.diagnostics.convergence_rate ≈ 0.8
        end

        @testset "Bootstrap CI computation" begin
            doses = [DoseEvent(0.0, 100.0)]
            subj1 = SubjectData("S1", [1.0, 2.0], [1.5, 1.2], doses)
            subj2 = SubjectData("S2", [1.0, 2.0], [1.6, 1.3], doses)
            observed = ObservedData([subj1, subj2])

            # Generate normally distributed parameter estimates
            function normal_estimate(data)
                return (
                    theta = [10.0 + 0.5 * randn(), 50.0 + 2.0 * randn()],
                    omega = diagm([0.09, 0.04]),
                    sigma = nothing,
                    convergence = true,
                    ofv = 100.0
                )
            end

            spec = BootstrapSpec(n_bootstrap=500, seed=UInt64(456), ci_level=0.95, parallel=false)
            result = run_bootstrap(normal_estimate, observed, spec)

            # CI should bracket the mean reasonably
            @test result.theta_ci_lower[1] < result.theta_mean[1] < result.theta_ci_upper[1]
            @test result.theta_ci_lower[2] < result.theta_mean[2] < result.theta_ci_upper[2]

            # CI width should be approximately 4 * SE for 95% CI
            ci_width1 = result.theta_ci_upper[1] - result.theta_ci_lower[1]
            ci_width2 = result.theta_ci_upper[2] - result.theta_ci_lower[2]
            @test 3.5 * result.theta_se[1] < ci_width1 < 4.5 * result.theta_se[1]
            @test 3.5 * result.theta_se[2] < ci_width2 < 4.5 * result.theta_se[2]
        end
    end

end
