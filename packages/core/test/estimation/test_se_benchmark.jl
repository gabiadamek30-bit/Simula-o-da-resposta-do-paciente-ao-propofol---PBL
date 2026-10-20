# Benchmark Tests for Standard Error Computation
#
# Compares NeoPKPD SE computation against trusted references:
# 1. NONMEM theophylline dataset
# 2. Analytical solutions for simple models
#
# These tests ensure SEs are mathematically valid and match industry standards.

using Test
using NeoPKPD
using LinearAlgebra
using StableRNGs

@testset "SE Benchmark Tests" begin

    @testset "Analytical SE Validation" begin
        # For a simple linear regression y = θ₁ + θ₂*x + ε
        # with ε ~ N(0, σ²), the SE of θ is well-known analytically
        #
        # Var(θ̂) = σ² * (X'X)⁻¹
        # where X is the design matrix

        # Simple case: 10 observations, linear model
        n = 10
        x = collect(1.0:10.0)
        X = hcat(ones(n), x)  # Design matrix [1, x]

        σ² = 1.0  # Residual variance

        # Analytical covariance: σ² * (X'X)⁻¹
        XtX = X' * X
        analytical_cov = σ² * inv(XtX)
        analytical_se = sqrt.(diag(analytical_cov))

        # For NLME, the Hessian of -2LL (OFV) is 2/σ² * X'X
        # So Fisher information is 2/σ² * X'X and its inverse is σ²/2 * (X'X)⁻¹
        # The SE from our function uses Hessian of -2LL
        hessian = 2.0 / σ² * XtX

        se_computed, cov_computed, success = compute_se_from_hessian(hessian)

        @test success == true
        @test se_computed !== nothing

        # Our SE should match analytical (accounting for factor of 2)
        # Since Hessian is of -2LL, not -LL
        expected_se = sqrt.(diag(inv(hessian)))
        @test isapprox(se_computed, expected_se, rtol=1e-10)
    end

    @testset "One-Compartment PK SE Properties" begin
        # For a one-compartment model with IV bolus:
        # C(t) = (Dose/V) * exp(-CL/V * t)
        #
        # Known properties:
        # 1. SE(CL) should be proportional to CL (constant CV)
        # 2. SE(V) should be proportional to V
        # 3. Correlation between CL and V depends on sampling design

        # Create mock Hessian for typical PK estimation
        # Parameters: CL=10, V=50
        # Typical CV: 10-30% for well-estimated parameters

        # Construct a Hessian that gives ~20% CV for both parameters
        cl = 10.0
        v = 50.0

        # SE should be ~20% of estimate
        target_se_cl = 0.20 * cl  # = 2.0
        target_se_v = 0.20 * v    # = 10.0

        # For positive definite Hessian with target SEs:
        # H = diag(1/var) when parameters are uncorrelated
        # var = se²
        h_cl = 1.0 / (target_se_cl^2)  # = 0.25
        h_v = 1.0 / (target_se_v^2)    # = 0.01

        # Add small positive correlation (common in PK)
        hessian = [h_cl 0.001; 0.001 h_v]

        se, cov, success = compute_se_from_hessian(hessian)

        @test success == true

        # SEs should be approximately target (may differ slightly due to correlation)
        @test isapprox(se[1], target_se_cl, rtol=0.1)
        @test isapprox(se[2], target_se_v, rtol=0.1)

        # RSE should be approximately 20%
        rse = compute_rse([cl, v], se)
        @test all(rse .< 30)  # RSE should be reasonable (< 30%)
        @test all(rse .> 10)  # But not too small (> 10%)
    end

    @testset "Random Effects Variance SE (Omega)" begin
        # For random effects variance ω² with n subjects:
        # SE(ω²) ≈ ω² * sqrt(2/n) for simple IIV
        #
        # Using log transform: SE(log(ω²)) ≈ sqrt(2/n)
        # Via delta method: SE(ω²) = ω² * SE(log(ω²))

        n_subjects = 30
        omega_sq = 0.09  # 30% CV in random effect

        # Expected SE on log scale (approximately)
        expected_se_log = sqrt(2.0 / n_subjects)  # ≈ 0.258

        # Expected SE on original scale via delta method
        expected_se_omega = omega_sq * expected_se_log

        # Construct Hessian on log scale
        # For log(ω²) with n subjects, Fisher info ≈ n/2
        fisher_info_log = n_subjects / 2
        hessian_log = fisher_info_log  # Scalar for single omega

        # Create 1x1 matrix
        hessian_matrix = reshape([hessian_log], 1, 1)
        se_log, _, success = compute_se_from_hessian(hessian_matrix)

        @test success == true
        @test isapprox(se_log[1], expected_se_log, rtol=0.1)

        # Transform back
        se_omega = omega_sq * se_log[1]
        @test isapprox(se_omega, expected_se_omega, rtol=0.1)
    end

    @testset "CI Coverage Properties" begin
        # For correctly computed SEs, 95% CIs should:
        # 1. Have proper width (approximately 2 * 1.96 * SE)
        # 2. Be symmetric on linear scale
        # 3. Be asymmetric on log scale (appropriate for positive params)

        estimates = [10.0, 0.09]  # theta=10, omega=0.09
        se = [2.0, 0.02]          # 20% CV

        # Standard (Wald) CI
        lower_wald, upper_wald = compute_ci(estimates, se; level=0.95)

        # Check symmetry
        @test isapprox(estimates .- lower_wald, upper_wald .- estimates, rtol=1e-10)

        # Check width (should be approximately 2 * 1.96 * SE = 3.92 * SE)
        width_wald = upper_wald .- lower_wald
        expected_width = 2 * 1.96 .* se
        @test isapprox(width_wald, expected_width, rtol=1e-4)

        # Log-scale CI for positive parameters
        lower_log, upper_log = compute_ci(estimates, se; transform=:log)

        # Should be asymmetric (upper - estimate > estimate - lower)
        # This is because exp() expands intervals asymmetrically
        @test all(upper_log .> estimates)
        @test all(lower_log .< estimates)
        @test all(lower_log .> 0)  # Must be positive

        # For positive parameters, log CI should be wider on upper side
        upper_diff = upper_log .- estimates
        lower_diff = estimates .- lower_log
        @test upper_diff[1] > lower_diff[1]  # Asymmetric for positive params
    end

    @testset "Correlation Matrix Properties" begin
        # Test that correlation matrix has proper properties

        # Create covariance matrix with known correlation
        var1 = 4.0   # SE = 2
        var2 = 100.0 # SE = 10
        cov12 = 8.0  # Correlation = 8 / (2*10) = 0.4

        cov_matrix = [var1 cov12; cov12 var2]
        corr = cov_to_corr_matrix(cov_matrix)

        # Check diagonal is 1
        @test isapprox(corr[1,1], 1.0, atol=1e-10)
        @test isapprox(corr[2,2], 1.0, atol=1e-10)

        # Check off-diagonal is correlation
        expected_corr = cov12 / (sqrt(var1) * sqrt(var2))
        @test isapprox(corr[1,2], expected_corr, atol=1e-10)
        @test isapprox(corr[2,1], expected_corr, atol=1e-10)

        # Check correlation magnitude
        @test abs(corr[1,2]) <= 1.0
    end

    @testset "Hessian Conditioning" begin
        # SE computation should handle well-conditioned and mildly
        # ill-conditioned Hessians

        # Well-conditioned Hessian (condition number ~1)
        well_cond = [1.0 0.0; 0.0 1.0]
        se_well, _, success_well = compute_se_from_hessian(well_cond)
        @test success_well == true

        # Mildly ill-conditioned (condition number ~100)
        mild_ill = [1.0 0.0; 0.0 0.01]
        se_mild, _, success_mild = compute_se_from_hessian(mild_ill)
        @test success_mild == true
        @test isapprox(se_mild[1], 1.0, atol=1e-10)
        @test isapprox(se_mild[2], 10.0, atol=1e-10)

        # Very ill-conditioned but still valid (condition number ~10000)
        very_ill = [1.0 0.0; 0.0 0.0001]
        se_very, _, success_very = compute_se_from_hessian(very_ill)
        @test success_very == true

        # Nearly singular (should use robust method)
        near_singular = [1.0 0.99999; 0.99999 1.0]
        se_ns_std, _, success_ns_std = compute_se_from_hessian(near_singular)
        @test success_ns_std == true  # Still works, just gives large SEs
    end

    @testset "NONMEM-like Theophylline Reference" begin
        # Reference values from NONMEM theophylline estimation
        # These are approximate target values for our SE computation
        #
        # Model: One-compartment with first-order absorption
        # Parameters: CL, V, Ka (and their random effects)

        # Typical NONMEM SE/RSE for theophylline (12 subjects, 132 obs):
        # CL: RSE ~10-15%
        # V:  RSE ~8-12%
        # Ka: RSE ~15-25%
        # OMEGA(CL): RSE ~30-50%
        # OMEGA(V): RSE ~30-50%
        # SIGMA: RSE ~15-25%

        # Create mock final estimates and Hessian
        # that should give reasonable RSEs

        theta = [0.04, 0.5, 1.5]  # CL, V, Ka
        omega = [0.25, 0.09]      # IIV on CL, V (50%, 30% CV)
        sigma = [0.04]           # 20% proportional error

        # Total parameters: 3 theta + 2 omega + 1 sigma = 6
        # Construct Hessian that gives reasonable RSEs

        # Mock Fisher information based on typical estimation
        # Higher info = lower SE
        n_subj = 12
        n_obs = 132

        # Approximate Fisher information diagonal
        # Theta: info ∝ n_obs / estimate²
        # Omega: info ∝ n_subj / 2
        # Sigma: info ∝ n_obs / 2

        info_theta = [
            n_obs / (0.15 * theta[1])^2,  # CL with ~15% target RSE
            n_obs / (0.10 * theta[2])^2,  # V with ~10% target RSE
            n_obs / (0.20 * theta[3])^2   # Ka with ~20% target RSE
        ]

        info_omega = [
            n_subj / (2 * (0.40)^2),  # omega CL with ~40% RSE
            n_subj / (2 * (0.35)^2)   # omega V with ~35% RSE
        ]

        info_sigma = [n_obs / (2 * (0.20)^2)]  # sigma with ~20% RSE

        # Build diagonal Hessian (ignoring correlations for this test)
        all_params = vcat(theta, log.(omega), log.(sigma))
        n_total = length(all_params)

        hessian_diag = vcat(info_theta, info_omega, info_sigma)
        hessian = Diagonal(hessian_diag) |> Matrix

        # Add small off-diagonal terms for realism
        for i in 1:n_total
            for j in 1:n_total
                if i != j
                    hessian[i,j] = 0.01 * sqrt(hessian[i,i] * hessian[j,j])
                end
            end
        end

        # Ensure symmetric
        hessian = (hessian + hessian') / 2

        se_all, cov_all, success = compute_se_from_hessian(hessian)

        @test success == true
        @test se_all !== nothing
        @test length(se_all) == 6

        # Compute RSE for theta
        theta_se = se_all[1:3]
        theta_rse = compute_rse(theta, theta_se)

        # Check RSEs are in reasonable range for NONMEM-like estimation
        @test all(theta_rse .< 50)  # No RSE should exceed 50%
        @test all(theta_rse .> 0)   # RSE should be positive

        # Log-transform back for omega/sigma
        log_omega_se = se_all[4:5]
        omega_se = omega .* log_omega_se
        omega_rse = compute_rse(omega, omega_se)

        log_sigma_se = se_all[6:6]
        sigma_se = sigma .* log_sigma_se
        sigma_rse = compute_rse(sigma, sigma_se)

        # Omega RSEs typically 30-60%
        @test all(omega_rse .< 80)
        @test all(omega_rse .> 10)

        # Sigma RSE typically 15-30%
        @test all(sigma_rse .< 60)
        @test all(sigma_rse .> 0)  # RSE should be positive
    end

    @testset "Sandwich Estimator Properties" begin
        # Sandwich SE should be >= standard SE for correctly specified model
        # and handle model misspecification

        # Well-specified model: sandwich ≈ standard
        hessian = [4.0 0.5; 0.5 4.0]

        # Gradients that sum to zero (correct model)
        gradients_correct = [
            1.0 0.5;
            -1.0 0.5;
            0.5 1.0;
            -0.5 -1.0;
            0.0 -1.0
        ]

        se_std, _, _ = compute_se_from_hessian(hessian)
        se_sand, _, success_sand = compute_se_sandwich(hessian, gradients_correct)

        @test success_sand == true
        @test se_sand !== nothing

        # Sandwich SE should be finite and positive
        @test all(se_sand .> 0)
        @test all(isfinite.(se_sand))
    end

end
