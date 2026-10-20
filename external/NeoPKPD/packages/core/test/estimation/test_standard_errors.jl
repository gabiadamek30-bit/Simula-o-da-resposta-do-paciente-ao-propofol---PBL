# Test suite for Standard Error Computation
# Tests marginal likelihood-based SEs with proper integration over random effects

using Test
using NeoPKPD
using LinearAlgebra
using StableRNGs

@testset "Standard Errors" begin

    @testset "SE Method Types" begin
        # Test LaplaceSE
        @test LaplaceSE() isa SEMethod

        # Test MonteCarloSE
        mc = MonteCarloSE()
        @test mc isa SEMethod
        @test mc.n_samples == 500  # Default

        mc2 = MonteCarloSE(n_samples=1000)
        @test mc2.n_samples == 1000

        # Test LouisSE
        louis = LouisSE()
        @test louis isa SEMethod
        @test louis.n_importance_samples == 200  # Default

        # Test SandwichSE
        @test SandwichSE() isa SEMethod
    end

    @testset "compute_se_from_hessian" begin
        # Test with positive definite Hessian
        hessian = [4.0 0.0; 0.0 4.0]
        se, cov, success = compute_se_from_hessian(hessian)

        @test success == true
        @test se !== nothing
        @test cov !== nothing
        @test isapprox(se[1], 0.5, atol=1e-10)
        @test isapprox(se[2], 0.5, atol=1e-10)
        @test isapprox(cov[1,1], 0.25, atol=1e-10)
        @test isapprox(cov[2,2], 0.25, atol=1e-10)

        # Test with non-positive definite Hessian (standard method)
        bad_hessian = [1.0 2.0; 2.0 1.0]
        se_bad, cov_bad, success_bad = compute_se_from_hessian(bad_hessian)
        @test success_bad == false
        @test se_bad === nothing

        # Test with non-positive definite Hessian (robust method)
        se_robust, cov_robust, success_robust = compute_se_from_hessian(bad_hessian; method=:robust)
        @test success_robust == true
        @test se_robust !== nothing

        # Test with correlated Hessian
        corr_hessian = [4.0 1.0; 1.0 4.0]
        se_corr, cov_corr, success_corr = compute_se_from_hessian(corr_hessian)
        @test success_corr == true
        @test se_corr !== nothing
        # Off-diagonal covariance should be non-zero
        @test cov_corr[1,2] != 0.0
    end

    @testset "compute_se_sandwich" begin
        # Simple test with identity-like Hessian and gradients
        hessian = [1.0 0.0; 0.0 1.0]
        gradients = [0.1 0.2; -0.1 -0.2; 0.05 0.1]

        se, cov, success = compute_se_sandwich(hessian, gradients)

        @test success == true
        @test se !== nothing
        @test length(se) == 2

        # Sandwich SEs should capture variance from gradients
        @test se[1] > 0
        @test se[2] > 0
    end

    @testset "compute_ci" begin
        estimates = [10.0, 50.0, 0.5]
        se = [1.0, 5.0, 0.1]

        # Standard Wald CI
        lower, upper = compute_ci(estimates, se)
        @test all(lower .< estimates)
        @test all(upper .> estimates)
        # Symmetric around estimate
        @test isapprox(estimates .- lower, upper .- estimates, atol=1e-10)

        # Check 95% CI width (approximately 1.96 * 2 * SE)
        widths = upper .- lower
        @test isapprox(widths[1], 2 * 1.96 * 1.0, atol=0.1)

        # Log-transformed CI (for positive parameters)
        lower_log, upper_log = compute_ci(estimates[1:2], se[1:2]; transform=:log)
        @test all(lower_log .> 0)  # Must be positive
        @test all(lower_log .< estimates[1:2])
        @test all(upper_log .> estimates[1:2])

        # Logit-transformed CI (for 0-1 parameters)
        lower_logit, upper_logit = compute_ci([0.5], [0.1]; transform=:logit)
        @test lower_logit[1] > 0
        @test upper_logit[1] < 1
    end

    @testset "compute_rse" begin
        estimates = [10.0, 50.0, 100.0]
        se = [1.0, 5.0, 20.0]

        rse = compute_rse(estimates, se)

        @test isapprox(rse[1], 10.0, atol=1e-10)  # 1/10 * 100 = 10%
        @test isapprox(rse[2], 10.0, atol=1e-10)  # 5/50 * 100 = 10%
        @test isapprox(rse[3], 20.0, atol=1e-10)  # 20/100 * 100 = 20%
    end

    @testset "cov_to_corr_matrix" begin
        # Covariance matrix
        cov = [0.09 0.03; 0.03 0.04]
        corr = cov_to_corr_matrix(cov)

        # Diagonal should be 1
        @test isapprox(corr[1,1], 1.0, atol=1e-10)
        @test isapprox(corr[2,2], 1.0, atol=1e-10)

        # Off-diagonal should be correlation
        expected_corr = 0.03 / (sqrt(0.09) * sqrt(0.04))
        @test isapprox(corr[1,2], expected_corr, atol=1e-10)
        @test isapprox(corr[2,1], expected_corr, atol=1e-10)

        # Correlation should be between -1 and 1
        @test abs(corr[1,2]) <= 1.0
    end

    @testset "bootstrap_se" begin
        rng = StableRNG(12345)
        data = randn(rng, 100) .+ 10.0  # Normal data with mean 10

        # Bootstrap SE of the mean
        mean_fn(x) = sum(x) / length(x)
        se_bootstrap = bootstrap_se(data, mean_fn, 1000, StableRNG(42))

        # Theoretical SE of mean is σ/√n ≈ 1/10 = 0.1
        @test isapprox(se_bootstrap, 0.1, atol=0.05)
    end

    @testset "StandardErrorResult Structure" begin
        # Test that all fields can be accessed
        result = StandardErrorResult(
            [1.0, 2.0],  # theta_se
            [0.01 0.0; 0.0 0.02],  # omega_se
            [0.05],  # sigma_se
            [1.0 0.1 0.0; 0.1 2.0 0.0; 0.0 0.0 0.5],  # covariance
            [1.0 0.5 0.0; 0.5 1.0 0.0; 0.0 0.0 1.0],  # correlation
            [8.0, 46.0],  # theta_ci_lower
            [12.0, 54.0],  # theta_ci_upper
            [0.005 0.0; 0.0 0.01],  # omega_ci_lower
            [0.02 0.0; 0.0 0.04],  # omega_ci_upper
            100.0,  # condition_number
            0.01,  # min_eigenvalue
            true,  # hessian_pd
            :LaplaceSE,  # method
            true,  # successful
            "Success"  # message
        )

        @test result.theta_se == [1.0, 2.0]
        @test result.successful == true
        @test result.method == :LaplaceSE
        @test result.hessian_pd == true
    end

    @testset "Hessian Computation" begin
        # Test central difference Hessian computation
        f(x) = x[1]^2 + x[2]^2 + x[1]*x[2]
        x = [1.0, 2.0]

        # Use the internal function
        hessian, success = NeoPKPD.compute_hessian_central(f, x)

        @test success == true
        # Expected Hessian for f(x) = x1² + x2² + x1*x2
        # ∂²f/∂x1² = 2, ∂²f/∂x2² = 2, ∂²f/∂x1∂x2 = 1
        @test isapprox(hessian[1,1], 2.0, atol=1e-4)
        @test isapprox(hessian[2,2], 2.0, atol=1e-4)
        @test isapprox(hessian[1,2], 1.0, atol=1e-4)
        @test isapprox(hessian[2,1], 1.0, atol=1e-4)

        # Test with non-finite function
        g(x) = x[1] < 0 ? Inf : x[1]^2
        y = [-1.0, 0.0]
        hessian_bad, success_bad = NeoPKPD.compute_hessian_central(g, y)
        @test success_bad == false
    end

    @testset "CI Width vs RSE Consistency" begin
        # For 95% CI, width ≈ 2 * 1.96 * SE
        # RSE = 100 * SE / estimate
        # So CI width / estimate ≈ 2 * 1.96 * RSE / 100

        estimates = [10.0, 100.0]
        se = [1.0, 10.0]  # Both have 10% RSE

        lower, upper = compute_ci(estimates, se; level=0.95)
        rse = compute_rse(estimates, se)

        # Check RSE
        @test all(isapprox.(rse, 10.0, atol=1e-10))

        # Check CI width is consistent with RSE
        ci_width = upper .- lower
        expected_width = 2 * 1.96 .* se
        @test isapprox(ci_width, expected_width, rtol=1e-4)  # Use relative tolerance
    end

    @testset "Delta Method for Log Transform" begin
        # For ω with SE(ω), the SE of log(ω) is SE(ω)/ω
        # Conversely, for SE on log scale, SE(ω) = ω * SE(log(ω))

        omega = 0.09  # Variance
        se_log_omega = 0.2  # SE on log scale (20% CV on log scale)

        # Delta method: SE(ω) = ω × SE(log(ω))
        se_omega = omega * se_log_omega
        @test isapprox(se_omega, 0.018, atol=1e-10)

        # Verify: RSE of omega
        rse_omega = 100 * se_omega / omega
        @test isapprox(rse_omega, 20.0, atol=1e-10)  # Should equal CV on log scale * 100
    end

    @testset "Covariance Matrix Properties" begin
        # A valid covariance matrix should be symmetric and positive semi-definite

        # Create a valid covariance matrix
        L = [1.0 0.0 0.0; 0.3 0.95 0.0; 0.1 0.2 0.97]
        cov = L * L'

        # Check symmetry
        @test issymmetric(cov) || isapprox(cov, cov', atol=1e-10)

        # Check positive definiteness
        eigenvalues = eigvals(Symmetric(cov))
        @test all(eigenvalues .> 0)

        # Convert to correlation
        corr = cov_to_corr_matrix(cov)

        # All correlations should be between -1 and 1
        for i in 1:3
            for j in 1:3
                @test abs(corr[i,j]) <= 1.0 + 1e-10
            end
        end

        # Diagonal of correlation should be 1
        for i in 1:3
            @test isapprox(corr[i,i], 1.0, atol=1e-10)
        end
    end

end
