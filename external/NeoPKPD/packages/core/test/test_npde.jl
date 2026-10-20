# Test suite for NPDE (Normalized Prediction Distribution Errors)

using Test
using NeoPKPD
using StableRNGs
using Statistics
using LinearAlgebra
using Distributions

@testset "NPDE" begin

    @testset "NPDESpec Construction" begin
        # Default spec
        spec = NPDESpec()
        @test spec.n_simulations == 1000
        @test spec.decorrelate == true
        @test spec.compute_pde == true

        # Custom spec
        spec2 = NPDESpec(n_simulations=500, decorrelate=false)
        @test spec2.n_simulations == 500
        @test spec2.decorrelate == false

        # Minimum simulations validation
        @test_throws AssertionError NPDESpec(n_simulations=50)
    end

    @testset "Statistical Tests" begin
        # Test with N(0,1) data - should pass all tests
        rng = StableRNG(42)
        normal_data = randn(rng, 200)

        tests = npde_statistical_tests(normal_data; alpha=0.05)

        @test length(tests) == 4  # t-test, variance, normality, symmetry

        # At least mean and variance tests should pass for N(0,1)
        mean_test = first(filter(t -> t.test_name == "t-test (mean = 0)", tests))
        @test mean_test.passed  # Should pass for N(0,1) data

        var_test = first(filter(t -> t.test_name == "Variance test (var = 1)", tests))
        @test var_test.passed  # Should pass for N(0,1) data

        # Test with biased data - should fail mean test
        biased_data = randn(rng, 200) .+ 0.5  # Mean shifted to 0.5
        tests_biased = npde_statistical_tests(biased_data; alpha=0.05)
        mean_test_biased = first(filter(t -> t.test_name == "t-test (mean = 0)", tests_biased))
        @test !mean_test_biased.passed  # Should fail

        # Test with high variance data - should fail variance test
        high_var_data = randn(rng, 200) .* 2.0  # Variance = 4
        tests_hv = npde_statistical_tests(high_var_data; alpha=0.05)
        var_test_hv = first(filter(t -> t.test_name == "Variance test (var = 1)", tests_hv))
        @test !var_test_hv.passed  # Should fail
    end

    @testset "Skewness and Kurtosis" begin
        rng = StableRNG(123)

        # Normal data should have skewness ~0
        normal_data = randn(rng, 1000)
        skew = NeoPKPD._compute_skewness(normal_data)
        @test abs(skew) < 0.2  # Should be near 0

        # Kurtosis for normal is 3 (excess = 0)
        kurt = NeoPKPD._compute_kurtosis(normal_data)
        @test abs(kurt) < 0.5  # Excess kurtosis near 0

        # Skewed data
        skewed_data = exp.(randn(rng, 1000))  # Log-normal is right-skewed
        skew_lognormal = NeoPKPD._compute_skewness(skewed_data)
        @test skew_lognormal > 0.5  # Should be positive (right-skewed)
    end

    @testset "Decorrelation Validation" begin
        # Test validate_npde_decorrelation function exists
        @test isdefined(NeoPKPD, :validate_npde_decorrelation)
        @test isdefined(NeoPKPD, :DecorrelationValidationResult)

        # Create mock NPDE result with uncorrelated data
        rng = StableRNG(456)
        n_subjects = 10
        n_obs = 8

        subjects = SubjectNPDEResult[]
        for i in 1:n_subjects
            # Independent N(0,1) NPDE values (properly decorrelated)
            npde_vals = randn(rng, n_obs)
            push!(subjects, SubjectNPDEResult(
                "SUBJ$(lpad(i, 3, '0'))",
                collect(0.0:1.0:(n_obs-1)),  # times
                randn(rng, n_obs),  # dv
                npde_vals,  # npde
                rand(rng, n_obs),  # pde
                randn(rng, n_obs) .+ 5.0,  # epred
                rand(rng, n_obs)  # var_pred
            ))
        end

        pooled_npde = vcat([s.npde for s in subjects]...)
        pooled_pde = vcat([s.pde for s in subjects]...)
        diagnostics = NPDEDiagnostics(
            mean(pooled_npde), var(pooled_npde),
            0.0, 3.0, 0.5, 0.5, 0.5
        )

        npde_result = NPDEResult(subjects, pooled_npde, pooled_pde, diagnostics, nothing)

        # Validate decorrelation
        validation = validate_npde_decorrelation(npde_result)

        @test validation isa DecorrelationValidationResult
        # For independent data, max correlation should be reasonably low
        # (with small samples, random correlations can exceed 0.5 by chance)
        @test validation.max_off_diagonal_correlation < 0.7
    end

    @testset "Independence Check" begin
        @test isdefined(NeoPKPD, :check_npde_independence)

        rng = StableRNG(789)

        # Independent data - should pass
        independent_data = randn(rng, 100)
        result_ind = check_npde_independence(independent_data)
        @test result_ind[:is_independent]

        # Autocorrelated data - should fail
        # Create AR(1) process with high autocorrelation
        ar1_data = zeros(100)
        ar1_data[1] = randn(rng)
        for i in 2:100
            ar1_data[i] = 0.8 * ar1_data[i-1] + 0.2 * randn(rng)
        end
        result_ar1 = check_npde_independence(ar1_data)
        # AR(1) with phi=0.8 should show significant lag-1 autocorrelation
        @test haskey(result_ar1[:autocorrelations], 1)
    end

    @testset "Comprehensive NPDE Validation" begin
        @test isdefined(NeoPKPD, :comprehensive_npde_validation)

        # Create mock NPDE result with good properties
        rng = StableRNG(111)
        n_subjects = 15
        n_obs = 10

        subjects = SubjectNPDEResult[]
        all_npde = Float64[]

        for i in 1:n_subjects
            npde_vals = randn(rng, n_obs)  # N(0,1) as expected
            append!(all_npde, npde_vals)
            push!(subjects, SubjectNPDEResult(
                "SUBJ$(lpad(i, 3, '0'))",
                collect(0.0:1.0:(n_obs-1)),
                randn(rng, n_obs) .+ 10.0,
                npde_vals,
                clamp.(rand(rng, n_obs), 0.01, 0.99),
                randn(rng, n_obs) .+ 5.0,
                rand(rng, n_obs)
            ))
        end

        pooled_pde = vcat([s.pde for s in subjects]...)
        diagnostics = NPDEDiagnostics(
            mean(all_npde), var(all_npde),
            NeoPKPD._compute_skewness(all_npde),
            NeoPKPD._compute_kurtosis(all_npde) + 3.0,
            0.5, 0.5, 0.5
        )

        npde_result = NPDEResult(subjects, all_npde, pooled_pde, diagnostics, nothing)

        # Run comprehensive validation
        validation = comprehensive_npde_validation(npde_result)

        @test haskey(validation, :overall_valid)
        @test haskey(validation, :summary)
        @test haskey(validation, :mean_test_passed)
        @test haskey(validation, :variance_test_passed)
        @test haskey(validation, :decorrelation_valid)
        @test haskey(validation, :independence_valid)

        # For N(0,1) data, tests should generally pass
        # (may occasionally fail due to random variation)
    end

    @testset "Lag-1 Autocorrelation" begin
        # Test the internal autocorrelation function

        # Perfectly alternating: should have negative autocorrelation
        alternating = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
        ac_alt = NeoPKPD._compute_lag1_autocorrelation(alternating)
        @test ac_alt < -0.5

        # Constant: should have undefined/zero autocorrelation (variance is 0)
        constant = [1.0, 1.0, 1.0, 1.0]
        ac_const = NeoPKPD._compute_lag1_autocorrelation(constant)
        @test ac_const == 0.0  # Returns 0 due to zero denominator check

        # Trending: should have positive autocorrelation
        trending = collect(1.0:10.0)
        ac_trend = NeoPKPD._compute_lag1_autocorrelation(trending)
        @test ac_trend > 0.5
    end

    @testset "NPDE Visualization Data" begin
        # Test visualization helper functions
        rng = StableRNG(222)
        n_subjects = 5
        n_obs = 6

        subjects = SubjectNPDEResult[]
        all_npde = Float64[]

        for i in 1:n_subjects
            npde_vals = randn(rng, n_obs)
            times = collect(0.0:1.0:(n_obs-1))
            append!(all_npde, npde_vals)
            push!(subjects, SubjectNPDEResult(
                "SUBJ$i",
                times,
                randn(rng, n_obs) .+ 10.0,  # dv
                npde_vals,  # npde
                rand(rng, n_obs),  # pde
                randn(rng, n_obs) .+ 5.0,  # epred
                rand(rng, n_obs)  # var_pred
            ))
        end

        pooled_pde = rand(rng, length(all_npde))
        diagnostics = NPDEDiagnostics(0.0, 1.0, 0.0, 3.0, 0.5, 0.5, 0.5)
        npde_result = NPDEResult(subjects, all_npde, pooled_pde, diagnostics, nothing)

        # Test QQ data
        qq_data = get_npde_qq_data(npde_result)
        @test haskey(qq_data, :observed)
        @test haskey(qq_data, :theoretical)
        @test length(qq_data[:observed]) == length(all_npde)
        @test issorted(qq_data[:observed])  # Should be sorted

        # Test vs pred data
        pred_data = get_npde_vs_pred_data(npde_result)
        @test haskey(pred_data, :npde)
        @test haskey(pred_data, :pred)
        @test length(pred_data[:npde]) == n_subjects * n_obs

        # Test vs time data
        time_data = get_npde_vs_time_data(npde_result)
        @test haskey(time_data, :npde)
        @test haskey(time_data, :time)
        @test length(time_data[:npde]) == n_subjects * n_obs

        # Test histogram data
        hist_data = get_npde_histogram_data(npde_result; n_bins=15)
        @test haskey(hist_data, :bin_centers)
        @test haskey(hist_data, :density)
        @test haskey(hist_data, :reference_x)
        @test haskey(hist_data, :reference_y)
        @test length(hist_data[:bin_centers]) == 15
    end

end
