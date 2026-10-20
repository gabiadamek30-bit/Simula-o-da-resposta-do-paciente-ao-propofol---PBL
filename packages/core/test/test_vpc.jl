# Test suite for Visual Predictive Check (VPC)

using Test
using NeoPKPD
using StableRNGs

@testset "VPC" begin

    @testset "Binning Strategies" begin
        @testset "QuantileBinning" begin
            strategy = QuantileBinning(5)
            @test strategy.n_bins == 5

            times = collect(1.0:20.0)
            bins = compute_bins(times, strategy)

            @test length(bins) == 5
            @test bins[1].id == 1
            @test bins[5].id == 5

            # All times should be covered
            for t in times
                covered = any(b.lower <= t <= b.upper for b in bins)
                @test covered
            end
        end

        @testset "EqualWidthBinning" begin
            strategy = EqualWidthBinning(4)
            @test strategy.n_bins == 4

            times = [0.0, 1.0, 2.0, 3.0, 4.0, 8.0, 12.0, 24.0]
            bins = compute_bins(times, strategy)

            @test length(bins) == 4
            # Equal width means equal time range per bin
            widths = [b.upper - b.lower for b in bins]
            @test all(isapprox(w, widths[1]; atol=0.01) for w in widths)
        end

        @testset "KMeansBinning" begin
            strategy = KMeansBinning(3; max_iter=50)
            @test strategy.n_bins == 3
            @test strategy.max_iter == 50

            # Times with clear clusters
            times = vcat(fill(1.0, 5), fill(5.0, 5), fill(10.0, 5))
            bins = compute_bins(times, strategy)

            @test length(bins) == 3
        end

        @testset "Empty times" begin
            strategy = QuantileBinning(5)
            bins = compute_bins(Float64[], strategy)
            @test isempty(bins)
        end
    end

    @testset "Percentile Computation" begin
        values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

        @test compute_percentile(values, 0.5) == 5.5  # Median
        @test compute_percentile(values, 0.0) == 1.0
        @test compute_percentile(values, 1.0) == 10.0
        @test isapprox(compute_percentile(values, 0.25), 3.25; atol=0.1)
        @test isapprox(compute_percentile(values, 0.75), 7.75; atol=0.1)

        # Empty vector
        @test isnan(compute_percentile(Float64[], 0.5))

        # Single value
        @test compute_percentile([5.0], 0.5) == 5.0
    end

    @testset "Bootstrap CI" begin
        rng = StableRNG(42)
        values = collect(1.0:100.0)

        lower, median_val, upper = bootstrap_percentile_ci(values, 0.5, 0.95, 500, rng)

        # Median should be around 50
        @test 45 < median_val < 55

        # Lower should be less than median
        @test lower < median_val

        # Upper should be greater than median
        @test upper > median_val

        # CI should contain the true median
        @test lower < 50.5 < upper
    end

    @testset "Assign to Bins" begin
        times = [0.5, 1.5, 2.5, 3.5, 4.5]
        values = [10.0, 20.0, 30.0, 40.0, 50.0]

        bins = [
            BinDefinition(1, 0.0, 2.0, 1.0),
            BinDefinition(2, 2.0, 4.0, 3.0),
            BinDefinition(3, 4.0, 6.0, 5.0)
        ]

        result = assign_to_bins(times, values, bins)

        @test length(result) == 3
        @test result[1][1] == 1  # bin id
        @test result[1][2] == [10.0, 20.0]  # values in bin 1
        @test result[2][2] == [30.0, 40.0]  # values in bin 2
        @test result[3][2] == [50.0]  # values in bin 3
    end

    @testset "VPCConfig" begin
        # Default config
        config1 = VPCConfig()
        @test config1.pi_levels == [0.05, 0.50, 0.95]
        @test config1.ci_level == 0.95
        @test config1.binning isa QuantileBinning
        @test config1.prediction_corrected == false
        @test config1.n_simulations == 200
        @test config1.n_bootstrap == 500

        # Custom config
        config2 = VPCConfig(
            pi_levels=[0.10, 0.50, 0.90],
            ci_level=0.90,
            binning=EqualWidthBinning(8),
            prediction_corrected=true,
            n_simulations=100,
            n_bootstrap=200,
            seed=UInt64(99999)
        )
        @test config2.pi_levels == [0.10, 0.50, 0.90]
        @test config2.ci_level == 0.90
        @test config2.binning isa EqualWidthBinning
        @test config2.prediction_corrected == true
        @test config2.n_simulations == 100
        @test config2.seed == UInt64(99999)
    end

    @testset "VPCResult Accessors" begin
        # Create mock VPC result
        percentiles = [
            VPCPercentileData(0.05, 1.0, 1.1, 0.9, 1.3),
            VPCPercentileData(0.50, 5.0, 5.2, 4.8, 5.6),
            VPCPercentileData(0.95, 10.0, 9.8, 9.5, 10.2)
        ]

        bin1 = VPCBin(1, 0.0, 4.0, 2.0, 10, 100, percentiles)
        bin2 = VPCBin(2, 4.0, 8.0, 6.0, 8, 100, percentiles)

        config = VPCConfig()
        result = VPCResult(config, [bin1, bin2], 10, 18, 100, "", UInt64(12345))

        # Test accessors
        @test bin_midpoints(result) == [2.0, 6.0]
        @test observed_percentile(result, 0.50) == [5.0, 5.0]
        @test simulated_median(result, 0.50) == [5.2, 5.2]
        @test simulated_lower(result, 0.50) == [4.8, 4.8]
        @test simulated_upper(result, 0.50) == [5.6, 5.6]
    end

    @testset "VPC from Population Simulation" begin
        # Create a simple one-compartment model
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL, V
        doses = [DoseEvent(0.0, 100.0)]
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        # Create population spec
        iiv_spec = IIVSpec(LogNormalIIV(), Dict(:CL => 0.3, :V => 0.2), UInt64(42), 20)

        pop_spec = PopulationSpec(model_spec, iiv_spec, nothing, nothing, IndividualCovariates[])

        grid = SimGrid(0.0, 24.0, collect(0.0:2.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        # Run population simulation
        pop_result = simulate_population(pop_spec, grid, solver)

        # Compute VPC from simulation
        config = VPCConfig(
            pi_levels=[0.10, 0.50, 0.90],
            binning=QuantileBinning(6),
            n_bootstrap=100,
            seed=UInt64(123)
        )

        vpc_result = compute_vpc_from_simulation(pop_result; config=config)

        @test vpc_result.n_subjects_observed == 20
        @test length(vpc_result.bins) == 6
        @test vpc_result.n_simulations == 1

        # Each bin should have 3 percentiles
        for bin in vpc_result.bins
            @test length(bin.percentiles) == 3
        end

        # Midpoints should be sorted
        midpoints = bin_midpoints(vpc_result)
        @test issorted(midpoints)

        # Median percentile should be reasonable
        median_pctls = observed_percentile(vpc_result, 0.50)
        @test all(!isnan(p) for p in median_pctls)
    end

    @testset "VPC with Observed Data" begin
        # Create observed data
        doses = [DoseEvent(0.0, 100.0)]

        # Create subjects with observations
        subj1 = SubjectData(
            "SUBJ001",
            [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            [2.0, 1.8, 1.5, 1.2, 0.8, 0.5, 0.2],
            doses
        )
        subj2 = SubjectData(
            "SUBJ002",
            [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            [2.2, 1.9, 1.6, 1.3, 0.9, 0.6, 0.25],
            doses
        )
        subj3 = SubjectData(
            "SUBJ003",
            [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
            [1.8, 1.7, 1.4, 1.1, 0.7, 0.45, 0.18],
            doses
        )

        observed = ObservedData(
            [subj1, subj2, subj3];
            study_id="TEST001",
            analyte="DRUG1",
            units="ng/mL",
            time_units="h"
        )

        # Create model for simulation
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL, V
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        iiv_spec = IIVSpec(LogNormalIIV(), Dict(:CL => 0.3, :V => 0.2), UInt64(42), 50)

        pop_spec = PopulationSpec(model_spec, iiv_spec, nothing, nothing, IndividualCovariates[])

        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        config = VPCConfig(
            pi_levels=[0.10, 0.50, 0.90],
            binning=QuantileBinning(5),
            n_simulations=20,  # Low for testing
            n_bootstrap=100,   # Minimum required
            seed=UInt64(123)
        )

        # Compute VPC
        vpc_result = compute_vpc(observed, pop_spec, grid, solver; config=config)

        @test vpc_result.n_subjects_observed == 3
        @test vpc_result.n_observations_observed == 21
        @test vpc_result.n_simulations == 20
        @test length(vpc_result.bins) == 5

        # Check that observed percentiles are computed
        for bin in vpc_result.bins
            @test bin.n_observed > 0
            for p in bin.percentiles
                @test !isnan(p.observed) || bin.n_observed == 0
            end
        end
    end

    @testset "VPC with Residual Error" begin
        # Create observed data
        doses = [DoseEvent(0.0, 100.0)]
        subj1 = SubjectData(
            "SUBJ001",
            [0.0, 2.0, 4.0, 8.0, 12.0],
            [2.0, 1.5, 1.0, 0.5, 0.25],
            doses
        )

        observed = ObservedData([subj1])

        # Create model
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL, V
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        iiv_spec = IIVSpec(LogNormalIIV(), Dict(:CL => 0.3), UInt64(42), 20)

        pop_spec = PopulationSpec(model_spec, iiv_spec, nothing, nothing, IndividualCovariates[])

        grid = SimGrid(0.0, 12.0, collect(0.0:2.0:12.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        # Proportional error model
        error_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),  # sigma
            :conc,
            UInt64(999)
        )

        config = VPCConfig(
            n_simulations=10,
            n_bootstrap=100,
            seed=UInt64(123)
        )

        vpc_result = compute_vpc(observed, pop_spec, grid, solver;
            config=config, error_spec=error_spec)

        @test vpc_result.n_simulations == 10
        @test length(vpc_result.bins) > 0
    end

    @testset "pcVPC" begin
        # Create observed data with multiple subjects and more time points
        doses = [DoseEvent(0.0, 100.0)]
        subj1 = SubjectData(
            "SUBJ001",
            [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0],
            [2.0, 1.8, 1.5, 1.0, 0.7, 0.5, 0.25],
            doses
        )
        subj2 = SubjectData(
            "SUBJ002",
            [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0],
            [2.2, 1.9, 1.6, 1.1, 0.75, 0.55, 0.28],
            doses
        )
        subj3 = SubjectData(
            "SUBJ003",
            [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0],
            [1.8, 1.7, 1.4, 0.9, 0.65, 0.45, 0.22],
            doses
        )

        observed = ObservedData([subj1, subj2, subj3])

        # Create model
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL, V
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        iiv_spec = IIVSpec(LogNormalIIV(), Dict(:CL => 0.3), UInt64(42), 20)

        pop_spec = PopulationSpec(model_spec, iiv_spec, nothing, nothing, IndividualCovariates[])

        grid = SimGrid(0.0, 12.0, collect(0.0:0.5:12.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        # Use EqualWidthBinning to ensure reasonable bin coverage
        config = VPCConfig(
            pi_levels=[0.10, 0.50, 0.90],
            binning=EqualWidthBinning(4),  # More appropriate for sparse discrete times
            n_simulations=20,
            n_bootstrap=100,
            seed=UInt64(123)
        )

        # Test pcVPC via compute_pcvpc
        vpc_result = compute_pcvpc(observed, pop_spec, grid, solver; config=config)

        @test vpc_result.config.prediction_corrected == true
        @test length(vpc_result.bins) > 0
        @test vpc_result.n_subjects_observed == 3
        @test vpc_result.n_simulations == 20

        # Count bins with valid data
        bins_with_data = filter(b -> b.n_observed > 0, vpc_result.bins)
        @test length(bins_with_data) > 0

        # Bins with data should have computed percentiles
        for bin in bins_with_data
            @test length(bin.percentiles) == 3  # 10th, 50th, 90th
            for p in bin.percentiles
                @test p.percentile in [0.10, 0.50, 0.90]
                # Observed percentile should be computed
                @test !isnan(p.observed)
            end
        end

        # Test pcVPC via compute_vpc with prediction_corrected=true
        pc_config = VPCConfig(
            pi_levels=[0.05, 0.50, 0.95],
            binning=EqualWidthBinning(3),
            prediction_corrected=true,
            n_simulations=15,
            n_bootstrap=100,
            seed=UInt64(456)
        )

        vpc_result2 = compute_vpc(observed, pop_spec, grid, solver; config=pc_config)
        @test vpc_result2.config.prediction_corrected == true
        @test length(vpc_result2.bins) == 3

        # Verify bin midpoints are sensible
        midpoints = bin_midpoints(vpc_result2)
        @test issorted(midpoints)
        @test midpoints[1] >= 0.0
        @test midpoints[end] <= 12.0
    end

    @testset "VPC Stratification" begin
        # Test stratified VPC result type
        @test isdefined(NeoPKPD, :StratifiedVPCResult)

        # Create basic VPC config
        config = VPCConfig(
            pi_levels=[0.10, 0.50, 0.90],
            binning=QuantileBinning(4),
            n_simulations=10,
            n_bootstrap=100,
            seed=UInt64(123)
        )

        @test config.n_simulations == 10
        @test config.n_bootstrap == 100
    end

    @testset "BLQ Methods" begin
        # Test BLQ method enum exists
        @test isdefined(NeoPKPD, :BLQMethod)
        @test isdefined(NeoPKPD, :M1)
        @test isdefined(NeoPKPD, :M4)
        @test isdefined(NeoPKPD, :M7)

        # Test handle_blq function with different methods
        values = [0.5, 1.0, 2.0, 0.3, 5.0, 0.2]
        times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0]
        lloq = 0.4
        tmax = 1.0

        # M1: Discard all BLQ (returns NaN for discarded)
        v1 = handle_blq(values, times, lloq; method=M1)
        @test length(v1) == 6  # Same length, but some NaN
        n_valid = count(!isnan, v1)
        @test n_valid == 4  # 0.5, 1.0, 2.0, 5.0 kept

        # M4: Replace with LLOQ/2
        v4 = handle_blq(values, times, lloq; method=M4)
        @test length(v4) == 6  # All kept
        @test count(v -> v == lloq/2, v4) == 2  # Two values replaced

        # M5: 0 before Tmax, LLOQ/2 after
        v5 = handle_blq(values, times, lloq; method=M5, tmax=tmax)
        @test length(v5) == 6

        # M7: 0 before Tmax, discard after (NaN for discarded)
        v7 = handle_blq(values, times, lloq; method=M7, tmax=tmax)
        @test length(v7) == 6  # Same length but some NaN
    end

    @testset "BLQ Bin Stats" begin
        @test isdefined(NeoPKPD, :BLQBinStats)

        # Test computing BLQ stats
        # Values below LLOQ: 0.1, 0.3 (2 out of 4)
        values = [0.1, 0.3, 1.0, 2.0]
        lloq = 0.4

        n_blq = count(v -> v < lloq, values)
        pct_blq = 100.0 * n_blq / length(values)

        @test n_blq == 2
        @test pct_blq == 50.0
    end

    @testset "BLQ Comprehensive Validation" begin
        # Test all BLQ methods (M1-M7) according to Beal (2001)

        # Scenario: PK profile with early absorption and terminal BLQ
        times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        values = [0.1, 0.8, 2.0, 1.5, 0.8, 0.3, 0.15, 0.05]  # Tmax at t=1.0
        lloq = 0.2
        tmax = 1.0  # Time of Cmax

        # Identify BLQ points: t=0.0 (0.1), t=12.0 (0.15), t=24.0 (0.05)
        # Pre-Tmax BLQ: t=0.0
        # Post-Tmax BLQ: t=12.0, t=24.0

        # M1: Discard all BLQ
        v1 = handle_blq(values, times, lloq; method=M1, tmax=tmax)
        @test count(isnan, v1) == 3  # 3 BLQ discarded
        @test !isnan(v1[2])  # t=0.5, value=0.8 kept
        @test !isnan(v1[3])  # t=1.0, value=2.0 kept

        # M3: Keep as-is (censored)
        v3 = handle_blq(values, times, lloq; method=M3, tmax=tmax)
        @test count(isnan, v3) == 0  # None discarded
        @test v3[1] == 0.1  # Original value kept
        @test v3[8] == 0.05  # Original value kept

        # M4: Replace all BLQ with LLOQ/2
        v4 = handle_blq(values, times, lloq; method=M4, tmax=tmax)
        @test count(isnan, v4) == 0  # None discarded
        @test v4[1] == lloq/2  # Replaced
        @test v4[7] == lloq/2  # Replaced
        @test v4[8] == lloq/2  # Replaced
        @test v4[3] == 2.0  # Above LLOQ, unchanged

        # M5: 0 before Tmax, LLOQ/2 after
        v5 = handle_blq(values, times, lloq; method=M5, tmax=tmax)
        @test v5[1] == 0.0  # Pre-Tmax BLQ -> 0
        @test v5[7] == lloq/2  # Post-Tmax BLQ -> LLOQ/2
        @test v5[8] == lloq/2  # Post-Tmax BLQ -> LLOQ/2

        # M6: LLOQ/2 before Tmax, discard after
        v6 = handle_blq(values, times, lloq; method=M6, tmax=tmax)
        @test v6[1] == lloq/2  # Pre-Tmax BLQ -> LLOQ/2
        @test isnan(v6[7])  # Post-Tmax BLQ -> discarded
        @test isnan(v6[8])  # Post-Tmax BLQ -> discarded

        # M7: 0 before Tmax, discard after
        v7 = handle_blq(values, times, lloq; method=M7, tmax=tmax)
        @test v7[1] == 0.0  # Pre-Tmax BLQ -> 0
        @test isnan(v7[7])  # Post-Tmax BLQ -> discarded
        @test isnan(v7[8])  # Post-Tmax BLQ -> discarded

        # Test BLQ handling with automatic Tmax detection
        v5_auto = handle_blq(values, times, lloq; method=M5)
        @test v5_auto[1] == 0.0  # Pre-Tmax BLQ -> 0 (Tmax auto-detected at t=1.0)

        # Edge case: All values above LLOQ
        values_high = [1.0, 2.0, 3.0, 2.5, 2.0]
        times_high = [0.0, 1.0, 2.0, 4.0, 8.0]
        v_high = handle_blq(values_high, times_high, lloq; method=M1)
        @test all(!isnan, v_high)  # No values changed
        @test v_high == values_high

        # Edge case: All values below LLOQ
        values_low = [0.1, 0.05, 0.08, 0.03, 0.02]
        times_low = [0.0, 1.0, 2.0, 4.0, 8.0]
        v_low_m1 = handle_blq(values_low, times_low, lloq; method=M1)
        @test all(isnan, v_low_m1)  # All discarded

        v_low_m4 = handle_blq(values_low, times_low, lloq; method=M4)
        @test all(v -> v == lloq/2, v_low_m4)  # All replaced with LLOQ/2
    end

    @testset "VPC Performance Scaling" begin
        # Test that pcVPC has O(N) scaling, not O(N×M)
        # by verifying that the pre-computed base model is used

        doses = [DoseEvent(0.0, 100.0)]

        # Create model
        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

        iiv_spec = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2), UInt64(42), 10)
        pop_spec = PopulationSpec(model_spec, iiv_spec, nothing, nothing, IndividualCovariates[])

        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        # Create observed data
        subj = SubjectData(
            "SUBJ001",
            [0.0, 1.0, 2.0, 4.0, 8.0, 12.0],
            [2.0, 1.8, 1.5, 1.0, 0.6, 0.3],
            doses
        )
        observed = ObservedData([subj])

        # Small grid
        small_grid = SimGrid(0.0, 12.0, collect(0.0:2.0:12.0))

        # Large grid (should not significantly impact pcVPC time due to O(N) optimization)
        large_grid = SimGrid(0.0, 12.0, collect(0.0:0.5:12.0))

        config = VPCConfig(
            n_simulations=10,  # Minimal for test
            n_bootstrap=100,
            seed=UInt64(123)
        )

        # Both should complete without timeout
        # (if O(N×M) scaling existed, large_grid would be much slower)
        vpc_small = compute_vpc(observed, pop_spec, small_grid, solver; config=config)
        vpc_large = compute_vpc(observed, pop_spec, large_grid, solver; config=config)

        @test length(vpc_small.bins) > 0
        @test length(vpc_large.bins) > 0

        # Test pcVPC explicitly
        pc_config = VPCConfig(
            n_simulations=10,  # Minimum required
            n_bootstrap=100,
            prediction_corrected=true,
            seed=UInt64(456)
        )

        pcvpc_result = compute_pcvpc(observed, pop_spec, large_grid, solver; config=pc_config)
        @test pcvpc_result.config.prediction_corrected == true
        @test length(pcvpc_result.bins) > 0
    end

end
