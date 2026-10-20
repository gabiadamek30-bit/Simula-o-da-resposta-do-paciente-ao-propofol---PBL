# Test suite for Bootstrap Parameter Uncertainty Estimation
# Tests regulatory-compliant bootstrap implementation

using Test
using NeoPKPD
using LinearAlgebra
using StableRNGs
using Statistics

@testset "Bootstrap" begin

    @testset "Bootstrap Type Definitions" begin
        @testset "CaseBootstrap" begin
            bt = CaseBootstrap()
            @test bt isa BootstrapType
        end

        @testset "ParametricBootstrap" begin
            bt = ParametricBootstrap()
            @test bt isa BootstrapType
            @test bt.n_simulations_per_subject == 1

            bt2 = ParametricBootstrap(5)
            @test bt2.n_simulations_per_subject == 5
        end

        @testset "ResidualBootstrap" begin
            bt = ResidualBootstrap()
            @test bt isa BootstrapType
            @test bt.standardize == true

            bt2 = ResidualBootstrap(false)
            @test bt2.standardize == false
        end
    end

    @testset "Bootstrap CI Methods" begin
        @testset "PercentileCI" begin
            ci = PercentileCI()
            @test ci isa BootstrapCIMethod
            @test NeoPKPD.ci_method_name(ci) == "Percentile"
        end

        @testset "BCCI" begin
            ci = BCCI()
            @test ci isa BootstrapCIMethod
            @test ci.acceleration == 0.0
            @test NeoPKPD.ci_method_name(ci) == "BCa"

            ci2 = BCCI(0.05)
            @test ci2.acceleration == 0.05
        end

        @testset "BasicCI" begin
            ci = BasicCI()
            @test ci isa BootstrapCIMethod
            @test NeoPKPD.ci_method_name(ci) == "Basic"
        end

        @testset "StudentizedCI" begin
            ci = StudentizedCI()
            @test ci isa BootstrapCIMethod
            @test NeoPKPD.ci_method_name(ci) == "Studentized"
        end
    end

    @testset "BootstrapSpec" begin
        @testset "Default values" begin
            spec = BootstrapSpec()
            @test spec.n_bootstrap == 1000
            @test spec.bootstrap_type isa CaseBootstrap
            @test isempty(spec.stratify_by)
            @test spec.ci_level == 0.95
            @test spec.ci_method isa PercentileCI
            @test spec.compute_omega_ci == true
            @test spec.compute_sigma_ci == true
            @test spec.min_success_rate == 0.8
        end

        @testset "Custom values" begin
            spec = BootstrapSpec(
                n_bootstrap=500,
                bootstrap_type=ParametricBootstrap(),
                stratify_by=[:study, :dose],
                ci_level=0.90,
                ci_method=BCCI(),
                compute_omega_ci=false,
                min_success_rate=0.9
            )
            @test spec.n_bootstrap == 500
            @test spec.bootstrap_type isa ParametricBootstrap
            @test :study in spec.stratify_by
            @test spec.ci_level == 0.90
            @test spec.ci_method isa BCCI
            @test spec.compute_omega_ci == false
            @test spec.min_success_rate == 0.9
        end

        @testset "Validation" begin
            @test_throws AssertionError BootstrapSpec(n_bootstrap=50)  # Must be >= 100
            @test_throws AssertionError BootstrapSpec(ci_level=0.0)   # Must be in (0, 1)
            @test_throws AssertionError BootstrapSpec(ci_level=1.0)   # Must be in (0, 1)
            @test_throws AssertionError BootstrapSpec(min_success_rate=0.0)  # Must be in (0, 1]
        end
    end

    @testset "Stratified Resampling" begin
        @testset "Simple resampling" begin
            rng = StableRNG(12345)
            indices = stratified_resample(10, rng)
            @test length(indices) == 10
            @test all(1 .<= indices .<= 10)
        end

        @testset "Stratified resampling" begin
            rng = StableRNG(12345)
            subject_ids = ["S001", "S002", "S003", "S004", "S005", "S006"]
            strata = [:A, :A, :B, :B, :B, :A]

            indices = stratified_resample(subject_ids, strata, rng)
            @test length(indices) == 6

            # Check that resampling preserves stratum counts
            # Original: 3 from stratum A (indices 1, 2, 6), 3 from stratum B (indices 3, 4, 5)
            stratum_a_indices = [1, 2, 6]
            stratum_b_indices = [3, 4, 5]

            n_from_a = sum(i in stratum_a_indices for i in indices)
            n_from_b = sum(i in stratum_b_indices for i in indices)

            @test n_from_a == 3  # Should have 3 from stratum A
            @test n_from_b == 3  # Should have 3 from stratum B
        end

        @testset "Reproducibility" begin
            rng1 = StableRNG(99999)
            rng2 = StableRNG(99999)

            indices1 = stratified_resample(20, rng1)
            indices2 = stratified_resample(20, rng2)

            @test indices1 == indices2
        end
    end

    @testset "Bootstrap CI Computation" begin
        # Create synthetic bootstrap samples
        rng = StableRNG(42)
        n_boot = 1000
        n_params = 3

        # Generate samples from known distribution
        true_mean = [10.0, 50.0, 1.5]
        true_sd = [1.0, 5.0, 0.2]
        estimates = zeros(n_boot, n_params)

        for j in 1:n_params
            estimates[:, j] = true_mean[j] .+ true_sd[j] .* randn(rng, n_boot)
        end

        original = true_mean

        @testset "Percentile CI" begin
            lower, upper = compute_bootstrap_ci(estimates, original, 0.95, PercentileCI())
            @test length(lower) == n_params
            @test length(upper) == n_params
            @test all(lower .< upper)

            # 95% CI should contain the mean for well-behaved distributions
            for j in 1:n_params
                @test lower[j] < true_mean[j] < upper[j]
            end
        end

        @testset "BCa CI" begin
            lower, upper = compute_bootstrap_ci(estimates, original, 0.95, BCCI())
            @test length(lower) == n_params
            @test length(upper) == n_params
            @test all(lower .< upper)
        end

        @testset "Basic CI" begin
            lower, upper = compute_bootstrap_ci(estimates, original, 0.95, BasicCI())
            @test length(lower) == n_params
            @test length(upper) == n_params
            @test all(lower .< upper)
        end

        @testset "Different CI levels" begin
            lower_90, upper_90 = compute_bootstrap_ci(estimates, original, 0.90, PercentileCI())
            lower_95, upper_95 = compute_bootstrap_ci(estimates, original, 0.95, PercentileCI())
            lower_99, upper_99 = compute_bootstrap_ci(estimates, original, 0.99, PercentileCI())

            # Wider CI should have lower lower bound and higher upper bound
            for j in 1:n_params
                @test lower_99[j] <= lower_95[j] <= lower_90[j]
                @test upper_90[j] <= upper_95[j] <= upper_99[j]
            end
        end
    end

    @testset "Omega Summary Computation" begin
        # Create synthetic omega estimates
        n_boot = 100
        n_eta = 2

        omega_estimates = [
            [0.1 + 0.01*randn() 0.0; 0.0 0.05 + 0.005*randn()]
            for _ in 1:n_boot
        ]

        summary = NeoPKPD._compute_omega_summary(omega_estimates, 0.95, PercentileCI())

        @test summary isa OmegaBootstrapSummary
        @test size(summary.mean) == (n_eta, n_eta)
        @test size(summary.se) == (n_eta, n_eta)
        @test size(summary.ci_lower) == (n_eta, n_eta)
        @test size(summary.ci_upper) == (n_eta, n_eta)

        # Mean should be close to 0.1 and 0.05 for diagonal
        @test abs(summary.mean[1, 1] - 0.1) < 0.02
        @test abs(summary.mean[2, 2] - 0.05) < 0.02

        # CI should contain the mean
        @test summary.ci_lower[1, 1] < summary.mean[1, 1] < summary.ci_upper[1, 1]
        @test summary.ci_lower[2, 2] < summary.mean[2, 2] < summary.ci_upper[2, 2]
    end

    @testset "Sigma Summary Computation" begin
        # Create synthetic sigma estimates
        rng = StableRNG(123)
        sigma_estimates = 0.1 .+ 0.01 .* randn(rng, 100)

        summary = NeoPKPD._compute_sigma_summary(sigma_estimates, 0.95)

        @test summary isa SigmaBootstrapSummary
        @test abs(summary.mean - 0.1) < 0.02
        @test summary.se > 0
        @test summary.ci_lower < summary.mean < summary.ci_upper
    end

    @testset "Bootstrap Diagnostics" begin
        diagnostics = BootstrapDiagnostics(
            950,                                # n_successful
            50,                                 # n_failed
            0.95,                               # convergence_rate
            25.0,                               # median_iterations
            [10, 25, 47],                       # outlier_indices
            Dict(:minimization => 30, :boundary => 20),  # failure_reasons
            [5.0, 4.5, 6.2]                     # rse_stability
        )

        @test diagnostics.n_successful == 950
        @test diagnostics.n_failed == 50
        @test diagnostics.convergence_rate == 0.95
        @test diagnostics.median_iterations == 25.0
        @test length(diagnostics.outlier_indices) == 3
        @test diagnostics.failure_reasons[:minimization] == 30
        @test length(diagnostics.rse_stability) == 3
    end

    @testset "Parallel Config" begin
        @testset "Bool parallel" begin
            spec_serial = BootstrapSpec(parallel=false)
            config_serial = NeoPKPD.get_parallel_config(spec_serial)
            @test config_serial.backend isa SerialBackend

            spec_parallel = BootstrapSpec(parallel=true)
            config_parallel = NeoPKPD.get_parallel_config(spec_parallel)
            @test config_parallel.backend isa ThreadedBackend
        end
    end

    @testset "Bootstrap Stability Analysis" begin
        # Create mock bootstrap result
        rng = StableRNG(555)
        n_boot = 500
        n_params = 3

        theta_estimates = zeros(n_boot, n_params)
        theta_estimates[:, 1] = 10.0 .+ randn(rng, n_boot)         # Low CV
        theta_estimates[:, 2] = 50.0 .+ 5.0 .* randn(rng, n_boot)  # Moderate CV
        theta_estimates[:, 3] = 1.5 .+ 0.1 .* randn(rng, n_boot)   # Low CV

        theta_mean = vec(mean(theta_estimates, dims=1))
        theta_se = vec(std(theta_estimates, dims=1))
        theta_rse = theta_se ./ abs.(theta_mean) .* 100

        result = BootstrapResult(
            theta_estimates,
            theta_mean,
            theta_se,
            theta_rse,
            theta_mean .- 1.96 .* theta_se,  # ci_lower
            theta_mean .+ 1.96 .* theta_se,  # ci_upper
            [10.0, 50.0, 1.5],                # original_estimate
            theta_mean .- [10.0, 50.0, 1.5],  # bias
            [10.0, 50.0, 1.5],                # bias_corrected
            nothing,                           # omega_summary
            nothing,                           # sigma_summary
            nothing,                           # eta_shrinkage
            BootstrapDiagnostics(n_boot, 0, 1.0, 20.0, Int[], Dict{Symbol,Int}(), Float64[]),
            0.95,
            "Percentile",
            nothing,
            nothing
        )

        stability = compute_bootstrap_stability(result)

        @test haskey(stability, :cv_percent)
        @test haskey(stability, :skewness)
        @test haskey(stability, :excess_kurtosis)
        @test haskey(stability, :stability_grade)

        @test length(stability[:cv_percent]) == n_params
        @test all(stability[:cv_percent] .> 0)
        @test stability[:stability_grade] in ["Excellent", "Good", "Acceptable", "Poor"]
    end

    @testset "Influential Subject Analysis" begin
        # Create mock bootstrap result with some outliers
        rng = StableRNG(777)
        n_boot = 200
        n_params = 2

        theta_estimates = zeros(n_boot, n_params)
        theta_estimates[:, 1] = 10.0 .+ randn(rng, n_boot)
        theta_estimates[:, 2] = 50.0 .+ 5.0 .* randn(rng, n_boot)

        # Add a few outliers
        theta_estimates[5, 1] = 20.0   # Outlier
        theta_estimates[10, 2] = 100.0 # Outlier

        theta_mean = vec(mean(theta_estimates, dims=1))
        theta_se = vec(std(theta_estimates, dims=1))
        theta_rse = theta_se ./ abs.(theta_mean) .* 100

        result = BootstrapResult(
            theta_estimates,
            theta_mean,
            theta_se,
            theta_rse,
            theta_mean .- 1.96 .* theta_se,
            theta_mean .+ 1.96 .* theta_se,
            [10.0, 50.0],
            theta_mean .- [10.0, 50.0],
            [10.0, 50.0],
            nothing,
            nothing,
            nothing,
            BootstrapDiagnostics(n_boot, 0, 1.0, 20.0, Int[], Dict{Symbol,Int}(), Float64[]),
            0.95,
            "Percentile",
            nothing,
            nothing
        )

        analysis = influential_subject_analysis(result)

        @test haskey(analysis, :high_influence_replicates)
        @test haskey(analysis, :stability_percent)
        @test haskey(analysis, :assessment)

        @test analysis[:assessment] in ["Stable", "Acceptable", "Unstable"]
    end

    @testset "Bootstrap Coverage" begin
        # Create mock result
        result = BootstrapResult(
            zeros(100, 2),         # theta_estimates
            [10.0, 50.0],          # theta_mean
            [1.0, 5.0],            # theta_se
            [10.0, 10.0],          # theta_rse
            [8.0, 40.0],           # theta_ci_lower
            [12.0, 60.0],          # theta_ci_upper
            [10.0, 50.0],          # original_estimate
            [0.0, 0.0],            # bias
            [10.0, 50.0],          # bias_corrected
            nothing,
            nothing,
            nothing,
            BootstrapDiagnostics(100, 0, 1.0, 20.0, Int[], Dict{Symbol,Int}(), Float64[]),
            0.95,
            "Percentile",
            nothing,
            nothing
        )

        # True values within CI
        coverage = compute_bootstrap_coverage(result, [10.5, 45.0])
        @test coverage[:coverage_rate] == 1.0  # Both parameters covered

        # One true value outside CI
        coverage2 = compute_bootstrap_coverage(result, [10.5, 65.0])
        @test coverage2[:coverage_rate] == 0.5  # Only first parameter covered
    end

    @testset "SE Comparison" begin
        # Create mock result
        result = BootstrapResult(
            zeros(100, 3),
            [10.0, 50.0, 1.5],
            [1.0, 5.0, 0.15],      # bootstrap SE
            [10.0, 10.0, 10.0],
            [8.0, 40.0, 1.2],
            [12.0, 60.0, 1.8],
            [10.0, 50.0, 1.5],
            [0.0, 0.0, 0.0],
            [10.0, 50.0, 1.5],
            nothing,
            nothing,
            nothing,
            BootstrapDiagnostics(100, 0, 1.0, 20.0, Int[], Dict{Symbol,Int}(), Float64[]),
            0.95,
            "Percentile",
            nothing,
            nothing
        )

        analytical_se = [0.95, 5.5, 0.14]

        comparison = compare_se_methods(result, analytical_se)

        @test haskey(comparison, :se_ratio)
        @test haskey(comparison, :relative_diff_percent)
        @test haskey(comparison, :discrepancy_flags)
        @test haskey(comparison, :mean_ratio)

        @test length(comparison[:se_ratio]) == 3
    end

end
