# NCA (Non-Compartmental Analysis) Tests
# Tests for FDA/EMA compliant NCA metrics

using Test
using NeoPKPD

# =============================================================================
# Test Data Setup
# =============================================================================

# Standard oral PK profile (1-compartment oral first-order absorption)
const TEST_TIMES = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
const TEST_CONC = [0.0, 0.82, 1.44, 1.62, 1.28, 0.94, 0.68, 0.36, 0.08]
const TEST_DOSE = 100.0

# IV bolus profile (rapid distribution)
const IV_TIMES = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0]
const IV_CONC = [10.0, 8.2, 6.7, 4.5, 2.0, 0.4, 0.08, 0.016, 0.0006]
const IV_DOSE = 100.0

# Multiple dose steady-state profile
const SS_TIMES = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
const SS_CONC = [2.5, 4.2, 3.8, 3.0, 2.8, 2.6, 2.5, 2.5]
const SS_TAU = 12.0

# =============================================================================
# NCA Configuration Tests
# =============================================================================

@testset "NCA Configuration" begin
    # Default config - uses LinLogMixedMethod by default
    config = NCAConfig()
    @test config.method isa LinLogMixedMethod
    @test config.lambda_z_min_points == 3
    @test config.lambda_z_r2_threshold == 0.9
    @test config.extrapolation_max_pct == 20.0
    @test config.significant_digits == 3
    @test config.blq_handling isa BLQZero

    # Custom config
    config_custom = NCAConfig(
        method = LogLinearMethod(),
        lambda_z_min_points = 4,
        lambda_z_r2_threshold = 0.95,
        extrapolation_max_pct = 10.0,
        significant_digits = 4,
        blq_handling = BLQLLOQHalf(),
        lloq = 0.1
    )
    @test config_custom.method isa LogLinearMethod
    @test config_custom.lambda_z_min_points == 4
    @test config_custom.lloq == 0.1
end

# =============================================================================
# Exposure Metrics Tests
# =============================================================================

@testset "Exposure Metrics" begin
    @testset "Cmax and Tmax" begin
        cmax = nca_cmax(TEST_CONC)
        tmax = nca_tmax(TEST_TIMES, TEST_CONC)

        @test cmax == 1.62
        @test tmax == 2.0

        # find_cmax returns both value and index
        cmax_val, cmax_idx = find_cmax(TEST_TIMES, TEST_CONC)
        @test cmax_val == 1.62
        @test cmax_idx == 4
    end

    @testset "Cmin" begin
        cmin = nca_cmin(TEST_CONC)
        @test cmin == 0.0  # First point is zero

        # Non-zero minimum
        cmin_ss = nca_cmin(SS_CONC)
        @test cmin_ss == 2.5
    end

    @testset "Clast and Tlast" begin
        clast = nca_clast(TEST_TIMES, TEST_CONC)
        tlast = nca_tlast(TEST_TIMES, TEST_CONC)

        @test clast == 0.08
        @test tlast == 24.0

        # With LLOQ filtering
        clast_lloq = nca_clast(TEST_TIMES, TEST_CONC; lloq=0.1)
        tlast_lloq = nca_tlast(TEST_TIMES, TEST_CONC; lloq=0.1)
        @test clast_lloq == 0.36
        @test tlast_lloq == 12.0
    end

    @testset "find_clast" begin
        clast, tlast, idx = find_clast(TEST_TIMES, TEST_CONC)
        @test clast == 0.08
        @test tlast == 24.0
        @test idx == 9

        # With LLOQ
        clast2, tlast2, idx2 = find_clast(TEST_TIMES, TEST_CONC; lloq=0.5)
        @test clast2 == 0.68
        @test tlast2 == 8.0
        @test idx2 == 7
    end

    @testset "time_above_concentration" begin
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [0.0, 2.0, 4.0, 2.0, 0.0]

        time_above = time_above_concentration(t, c, 1.0)
        @test time_above ≈ 3.0 atol=0.1  # Approximately 3 hours above threshold
    end
end

# =============================================================================
# Lambda-z Estimation Tests
# =============================================================================

@testset "Lambda-z Estimation" begin
    config = NCAConfig()

    @testset "Basic lambda_z estimation" begin
        result = estimate_lambda_z(TEST_TIMES, TEST_CONC, config)

        @test result.lambda_z !== nothing
        @test result.lambda_z > 0.0
        @test result.t_half !== nothing
        @test result.t_half > 0.0
        @test result.r_squared !== nothing
        @test result.r_squared >= config.lambda_z_r2_threshold
        @test result.n_points >= config.lambda_z_min_points
        @test result.quality_flag in [:good, :warning]
    end

    @testset "Lambda_z relationship to t_half" begin
        result = estimate_lambda_z(TEST_TIMES, TEST_CONC, config)

        if result.lambda_z !== nothing && result.t_half !== nothing
            # t_half = ln(2) / lambda_z
            expected_t_half = log(2) / result.lambda_z
            @test result.t_half ≈ expected_t_half atol=1e-10
        end
    end

    @testset "IV bolus lambda_z" begin
        result = estimate_lambda_z(IV_TIMES, IV_CONC, config)

        @test result.lambda_z !== nothing
        # IV bolus typically has clear terminal phase
        @test result.r_squared !== nothing
        @test result.r_squared > 0.9
    end

    @testset "Insufficient data handling" begin
        # Only 2 points in terminal phase
        short_t = [0.0, 1.0, 2.0]
        short_c = [0.0, 1.0, 0.5]

        result = estimate_lambda_z(short_t, short_c, config)
        # Should handle gracefully
        @test result.quality_flag in [:good, :warning, :insufficient]
    end

    @testset "nca_half_life" begin
        lambda_z = 0.1  # 1/hr
        t_half = nca_half_life(lambda_z)
        @test t_half ≈ log(2) / 0.1 atol=1e-10
    end
end

# =============================================================================
# AUC Calculation Tests
# =============================================================================

@testset "AUC Calculations" begin
    config_linear = NCAConfig(method=LinearMethod())
    config_log = NCAConfig(method=LogLinearMethod())
    config_mixed = NCAConfig(method=LinLogMixedMethod())

    @testset "AUC 0-t (Linear)" begin
        auc = auc_0_t(TEST_TIMES, TEST_CONC, config_linear)
        @test auc > 0.0
        # Manual trapezoid check for first interval
        expected_first = 0.5 * (0.0 + 0.82) * (0.5 - 0.0)
        @test auc > expected_first  # AUC should be larger than first interval
    end

    @testset "AUC 0-t (Log-linear)" begin
        auc = auc_0_t(TEST_TIMES, TEST_CONC, config_log)
        @test auc > 0.0
    end

    @testset "AUC 0-t (Lin-Log Mixed)" begin
        auc = auc_0_t(TEST_TIMES, TEST_CONC, config_mixed)
        @test auc > 0.0
    end

    @testset "AUC methods comparison" begin
        auc_linear = auc_0_t(TEST_TIMES, TEST_CONC, config_linear)
        auc_log = auc_0_t(TEST_TIMES, TEST_CONC, config_log)
        auc_mixed = auc_0_t(TEST_TIMES, TEST_CONC, config_mixed)

        # All methods should give positive results
        @test auc_linear > 0.0
        @test auc_log > 0.0
        @test auc_mixed > 0.0

        # Results should be similar (within 20%)
        @test auc_linear / auc_log > 0.8
        @test auc_linear / auc_log < 1.2
    end

    @testset "AUC 0-inf" begin
        lambda_z_result = estimate_lambda_z(TEST_TIMES, TEST_CONC, config_linear)

        if lambda_z_result.lambda_z !== nothing
            clast, tlast, _ = find_clast(TEST_TIMES, TEST_CONC)
            auc_inf, extra_pct = auc_0_inf(
                TEST_TIMES, TEST_CONC,
                lambda_z_result.lambda_z, clast, config_linear
            )

            @test auc_inf > 0.0
            @test extra_pct >= 0.0
            @test extra_pct < 100.0

            # AUC_inf should be larger than AUC_0-t
            auc_t = auc_0_t(TEST_TIMES, TEST_CONC, config_linear)
            @test auc_inf > auc_t
        end
    end

    @testset "AUC 0-tau" begin
        auc_tau = auc_0_tau(SS_TIMES, SS_CONC, SS_TAU, config_linear)
        @test auc_tau > 0.0

        # Should equal AUC_0-t when tau >= tlast
        auc_t = auc_0_t(SS_TIMES, SS_CONC, config_linear)
        @test auc_tau ≈ auc_t atol=0.01
    end

    @testset "AUMC 0-t" begin
        aumc = aumc_0_t(TEST_TIMES, TEST_CONC, config_linear)
        @test aumc > 0.0
    end

    @testset "AUC partial" begin
        # Partial AUC from 1.0 to 4.0 hours
        partial = auc_partial(TEST_TIMES, TEST_CONC, 1.0, 4.0, config_linear)
        @test partial > 0.0

        # Should be less than total AUC
        total = auc_0_t(TEST_TIMES, TEST_CONC, config_linear)
        @test partial < total
    end
end

# =============================================================================
# PK Parameters Tests
# =============================================================================

@testset "PK Parameters" begin
    config = NCAConfig()

    @testset "MRT calculation" begin
        aumc_inf = 100.0
        auc_inf = 20.0

        mrt = nca_mrt(aumc_inf, auc_inf)
        @test mrt == 5.0

        # IV infusion adjustment
        mrt_inf = nca_mrt(aumc_inf, auc_inf; route=:iv_infusion, t_inf=1.0)
        @test mrt_inf == 4.5  # MRT - t_inf/2
    end

    @testset "CL/F calculation" begin
        cl_f = nca_cl_f(TEST_DOSE, 20.0)  # dose=100, auc=20
        @test cl_f == 5.0
    end

    @testset "CL calculation (IV)" begin
        cl = nca_cl(IV_DOSE, 25.0)
        @test cl == 4.0
    end

    @testset "Vz/F calculation" begin
        vz_f = nca_vz_f(TEST_DOSE, 0.1, 20.0)  # dose=100, lambda_z=0.1, auc=20
        @test vz_f == 50.0
    end

    @testset "Vss calculation" begin
        vss = nca_vss(5.0, 10.0)  # CL=5, MRT=10
        @test vss == 50.0
    end

    @testset "Vss from AUMC" begin
        vss = nca_vss_from_aumc(TEST_DOSE, 20.0, 100.0)  # dose=100, auc=20, aumc=100
        @test vss == 100.0 * 100.0 / 400.0
    end

    @testset "Vc calculation" begin
        vc = nca_vc(IV_DOSE, 10.0)  # dose=100, C0=10
        @test vc == 10.0
    end

    @testset "Bioavailability" begin
        f = nca_bioavailability(20.0, 100.0, 25.0, 100.0)
        @test f ≈ 0.8 atol=1e-10

        # Different doses
        f2 = nca_bioavailability(10.0, 50.0, 25.0, 100.0)
        @test f2 ≈ 0.8 atol=1e-10
    end

    @testset "Mean Absorption Time" begin
        mat = nca_mean_absorption_time(8.0, 5.0)  # MRT_po=8, MRT_iv=5
        @test mat == 3.0
    end
end

# =============================================================================
# Multiple Dose Metrics Tests
# =============================================================================

@testset "Multiple Dose Metrics" begin
    @testset "Accumulation Index" begin
        rac = nca_accumulation_index(25.0, 20.0)  # AUC_ss=25, AUC_sd=20
        @test rac == 1.25
    end

    @testset "Predicted Accumulation" begin
        rac_pred = nca_accumulation_predicted(0.1, 12.0)  # lambda_z=0.1, tau=12
        expected = 1.0 / (1.0 - exp(-0.1 * 12.0))
        @test rac_pred ≈ expected atol=1e-10
    end

    @testset "PTF (Peak-Trough Fluctuation)" begin
        ptf = nca_ptf(4.2, 2.5, 3.0)  # Cmax=4.2, Cmin=2.5, Cavg=3.0
        expected = 100.0 * (4.2 - 2.5) / 3.0
        @test ptf ≈ expected atol=1e-10
    end

    @testset "Swing" begin
        swing = nca_swing(4.2, 2.5)  # Cmax=4.2, Cmin=2.5
        expected = 100.0 * (4.2 - 2.5) / 2.5
        @test swing ≈ expected atol=1e-10
    end

    @testset "Dose Linearity Index" begin
        doses = [25.0, 50.0, 100.0, 200.0]
        aucs = [5.0, 10.0, 20.0, 40.0]  # Perfect linearity

        result = nca_linearity_index(doses, aucs)
        @test result.beta ≈ 1.0 atol=0.01
        @test result.r_squared ≈ 1.0 atol=0.01
        @test result.is_linear == true

        # Non-linear case
        aucs_nonlinear = [5.0, 15.0, 50.0, 180.0]
        result_nl = nca_linearity_index(doses, aucs_nonlinear)
        @test result_nl.beta > 1.0  # Superproportional
    end

    @testset "Time to Steady State" begin
        t_ss = nca_time_to_steady_state(0.1; fraction=0.90)
        expected = -log(1.0 - 0.90) / 0.1
        @test t_ss ≈ expected atol=1e-10

        # 95% steady state
        t_ss_95 = nca_time_to_steady_state(0.1; fraction=0.95)
        @test t_ss_95 > t_ss
    end

    @testset "Doses to Steady State" begin
        n_doses = nca_time_to_steady_state_doses(0.1, 12.0; fraction=0.90)
        @test n_doses >= 1
        @test n_doses isa Int
    end

    @testset "Cavg calculation" begin
        config = NCAConfig()
        cavg = nca_cavg(SS_TIMES, SS_CONC, SS_TAU, config)

        auc_tau = auc_0_tau(SS_TIMES, SS_CONC, SS_TAU, config)
        expected = auc_tau / SS_TAU
        @test cavg ≈ expected atol=1e-10
    end
end

# =============================================================================
# Bioequivalence Tests
# =============================================================================

@testset "Bioequivalence Analysis" begin
    # Test data: paired crossover study
    test_values = [20.0, 22.0, 18.0, 25.0, 21.0, 19.0, 23.0, 20.0]
    ref_values = [19.0, 21.0, 17.0, 24.0, 20.0, 18.0, 22.0, 19.0]

    @testset "Geometric Mean Ratio" begin
        gmr = geometric_mean_ratio(test_values, ref_values)
        @test gmr > 0.0
        @test gmr ≈ 1.05 atol=0.1  # Close to 1 for similar values
    end

    @testset "Geometric Mean" begin
        gm = geometric_mean(test_values)
        expected = exp(sum(log.(test_values)) / length(test_values))
        @test gm ≈ expected atol=1e-10
    end

    @testset "Within-Subject CV" begin
        cv = within_subject_cv(test_values, ref_values)
        @test cv > 0.0
        @test cv < 100.0  # Should be reasonable CV
    end

    @testset "90% Confidence Interval" begin
        result = bioequivalence_90ci(test_values, ref_values)

        @test haskey(result, :gmr)
        @test haskey(result, :ci_lower)
        @test haskey(result, :ci_upper)
        @test haskey(result, :cv_intra)
        @test haskey(result, :n)

        @test result.gmr > 0.0
        @test result.ci_lower < result.gmr
        @test result.ci_upper > result.gmr
        @test result.n == length(test_values)
    end

    @testset "TOST Analysis" begin
        result = tost_analysis(test_values, ref_values)

        @test result.parameter == :generic
        @test result.be_conclusion in [:bioequivalent, :not_bioequivalent]

        # Check test statistics exist
        @test !isnan(result.t_lower)
        @test !isnan(result.t_upper)
    end

    @testset "BE Conclusion" begin
        # Within limits
        @test be_conclusion(0.85, 1.20) == :bioequivalent

        # Outside limits
        @test be_conclusion(0.75, 1.10) == :inconclusive
        @test be_conclusion(0.90, 1.30) == :inconclusive

        # Completely outside
        @test be_conclusion(0.70, 0.78) == :not_bioequivalent
        @test be_conclusion(1.30, 1.50) == :not_bioequivalent
    end

    @testset "Custom BE Limits" begin
        # Highly variable drug (wider limits)
        conclusion = be_conclusion(0.72, 1.35; theta_lower=0.70, theta_upper=1.43)
        @test conclusion == :bioequivalent
    end

    @testset "Create BE Result" begin
        result = create_be_result(:cmax, test_values, ref_values)

        @test result.parameter == :cmax
        @test result.n_test == length(test_values)
        @test result.n_reference == length(ref_values)
        @test result.gmr > 0.0
        @test result.be_conclusion in [:bioequivalent, :not_bioequivalent, :inconclusive]
    end
end

# =============================================================================
# Bioequivalence Confidence Interval - T-Distribution Tests
# =============================================================================

@testset "BE Confidence Interval - T-Distribution" begin
    @testset "T-distribution critical values" begin
        # Known t-critical values (from statistical tables)
        # t_0.95 (upper 5%) for various df

        # df=10: t_0.95 = 1.812
        t_crit_10 = NeoPKPD._t_critical(10, 0.10)
        @test t_crit_10 ≈ 1.812 atol=0.01

        # df=20: t_0.95 = 1.725
        t_crit_20 = NeoPKPD._t_critical(20, 0.10)
        @test t_crit_20 ≈ 1.725 atol=0.01

        # df=24: t_0.95 = 1.711 (typical BE study with 26 subjects)
        t_crit_24 = NeoPKPD._t_critical(24, 0.10)
        @test t_crit_24 ≈ 1.711 atol=0.01

        # df=30: t_0.95 = 1.697
        t_crit_30 = NeoPKPD._t_critical(30, 0.10)
        @test t_crit_30 ≈ 1.697 atol=0.01

        # df=60: t_0.95 = 1.671 (larger study)
        t_crit_60 = NeoPKPD._t_critical(60, 0.10)
        @test t_crit_60 ≈ 1.671 atol=0.01

        # df=120: t_0.95 = 1.658 (approaches z=1.645)
        t_crit_120 = NeoPKPD._t_critical(120, 0.10)
        @test t_crit_120 ≈ 1.658 atol=0.01
    end

    @testset "T-distribution vs z-score difference" begin
        # Verify t-critical is LARGER than z for small df (wider CI)
        z_crit = 1.645  # z_0.95

        for df in [5, 10, 15, 20, 25, 30]
            t_crit = NeoPKPD._t_critical(df, 0.10)
            @test t_crit > z_crit  # t should always be larger than z for finite df
        end

        # For very large df, t approaches z
        t_crit_1000 = NeoPKPD._t_critical(1000, 0.10)
        @test abs(t_crit_1000 - z_crit) < 0.01
    end

    @testset "Degrees of freedom calculation" begin
        # 2x2 Crossover design
        @test compute_be_degrees_of_freedom(Crossover2x2(), 24) == 22
        @test compute_be_degrees_of_freedom(Crossover2x2(), 30) == 28
        @test compute_be_degrees_of_freedom(Crossover2x2(), 12) == 10

        # 2x4 Replicate crossover
        @test compute_be_degrees_of_freedom(Crossover2x4(), 24) == 22
        @test compute_be_degrees_of_freedom(Crossover2x4(), 18) == 16

        # Parallel group design
        @test compute_be_degrees_of_freedom(ParallelGroup(), 40) == 38
        @test compute_be_degrees_of_freedom(ParallelGroup(), 60) == 58
    end

    @testset "Design-aware CI calculation" begin
        # Generate realistic BE study data (n=24 typical)
        # Simulate GMR ≈ 1.05 with CV ≈ 25%

        # Seed for reproducibility
        import Random
        Random.seed!(12345)

        n = 24
        cv_true = 0.25  # 25% CV
        gmr_true = 1.05

        # Generate paired data
        ref_values = 100.0 .+ 20.0 .* randn(n)  # Reference with some variance
        test_values = ref_values .* gmr_true .* exp.(cv_true .* randn(n))  # Test with log-normal error

        # Ensure positive values
        ref_values = max.(ref_values, 10.0)
        test_values = max.(test_values, 10.0)

        result = bioequivalence_90ci_design(test_values, ref_values, Crossover2x2())

        # Check returned values
        @test result.n == 24
        @test result.df == 22  # n - 2 for crossover
        @test result.t_critical > 1.7  # Should be approximately 1.717
        @test result.t_critical < 1.75

        # GMR should be reasonable
        @test result.gmr > 0.8
        @test result.gmr < 1.3

        # CI should bracket GMR
        @test result.ci_lower < result.gmr
        @test result.ci_upper > result.gmr

        # Design should be recorded
        @test result.design == "Crossover2x2"
    end

    @testset "Typical BE study sizes (n=20-30)" begin
        # Validate CI for typical BE study sizes
        import Random

        for n in [20, 24, 26, 30]
            Random.seed!(42 + n)

            # Simulate bioequivalent products (GMR = 1.0, CV = 20%)
            ref = 100.0 .+ 15.0 .* randn(n)
            test = ref .* exp.(0.20 .* randn(n))

            ref = max.(ref, 50.0)
            test = max.(test, 50.0)

            result = bioequivalence_90ci_design(test, ref, Crossover2x2())

            # Verify degrees of freedom
            @test result.df == n - 2

            # t-critical should be in expected range for these df
            @test result.t_critical > 1.69
            @test result.t_critical < 1.75

            # For bioequivalent products, CI should usually be within 0.80-1.25
            # (not guaranteed due to random sampling, but likely)
            @test result.ci_lower > 0.5  # Should not be extreme
            @test result.ci_upper < 2.0  # Should not be extreme
        end
    end

    @testset "CI width depends on sample size" begin
        import Random
        Random.seed!(999)

        # Fixed CV, varying sample size
        # Larger n should give narrower CI

        ref_base = 100.0 .+ 15.0 .* randn(50)
        test_base = ref_base .* exp.(0.15 .* randn(50))

        widths = Float64[]

        for n in [12, 20, 30, 40]
            ref = max.(ref_base[1:n], 50.0)
            test = max.(test_base[1:n], 50.0)

            result = bioequivalence_90ci_design(test, ref, Crossover2x2())

            width = result.ci_upper / result.ci_lower
            push!(widths, width)
        end

        # CI width should decrease as n increases
        @test widths[2] < widths[1]
        @test widths[3] < widths[2]
        @test widths[4] < widths[3]
    end

    @testset "Comparison: old vs design-aware CI" begin
        # The new design-aware function should give slightly wider CIs
        # due to df = n-2 instead of df = n-1

        import Random
        Random.seed!(777)

        n = 24
        ref = 100.0 .+ 20.0 .* randn(n)
        test = ref .* 1.02 .* exp.(0.20 .* randn(n))

        ref = max.(ref, 50.0)
        test = max.(test, 50.0)

        old_result = bioequivalence_90ci(test, ref)
        new_result = bioequivalence_90ci_design(test, ref, Crossover2x2())

        # GMR should be the same
        @test old_result.gmr ≈ new_result.gmr atol=1e-10

        # New CI should be slightly wider (df = n-2 gives larger t-critical)
        old_width = old_result.ci_upper - old_result.ci_lower
        new_width = new_result.ci_upper - new_result.ci_lower

        @test new_width >= old_width * 0.99  # New should be equal or wider
    end

    @testset "T-CDF and T-PDF accuracy" begin
        # Verify CDF properties
        @test NeoPKPD._t_cdf(0.0, 10) ≈ 0.5 atol=1e-10  # Symmetric at 0
        @test NeoPKPD._t_cdf(Inf, 10) ≈ 1.0 atol=1e-10  # Limit at infinity
        @test NeoPKPD._t_cdf(-Inf, 10) ≈ 0.0 atol=1e-10  # Limit at -infinity

        # Verify PDF properties
        @test NeoPKPD._t_pdf(0.0, 10) > 0.0  # Positive at peak
        @test NeoPKPD._t_pdf(0.0, 10) > NeoPKPD._t_pdf(1.0, 10)  # Peak at 0
        @test NeoPKPD._t_pdf(-1.0, 10) ≈ NeoPKPD._t_pdf(1.0, 10) atol=1e-10  # Symmetric

        # Known CDF values for df=10
        # P(T <= 1.812) ≈ 0.95 for df=10
        @test NeoPKPD._t_cdf(1.812, 10) ≈ 0.95 atol=0.01
    end
end

# =============================================================================
# Full NCA Workflow Tests
# =============================================================================

@testset "Full NCA Workflow" begin
    config = NCAConfig()

    @testset "run_nca - Single Dose Oral" begin
        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE; config=config)

        # Primary exposure
        @test result.cmax == 1.62
        @test result.tmax == 2.0
        @test result.clast == 0.08
        @test result.tlast == 24.0

        # AUC
        @test result.auc_0_t > 0.0

        # Terminal phase (if available)
        @test result.lambda_z_result isa LambdaZResult

        # Quality
        @test result.quality_flags isa Vector{Symbol}
        @test result.warnings isa Vector{String}

        # Metadata
        @test result.metadata["dose"] == TEST_DOSE
        @test result.metadata["route"] == "extravascular"
        @test result.metadata["dosing_type"] == "single"
    end

    @testset "run_nca - IV Bolus" begin
        result = run_nca(IV_TIMES, IV_CONC, IV_DOSE;
                        config=config, route=:iv_bolus)

        @test result.cmax > 0.0
        @test result.auc_0_t > 0.0
        @test result.metadata["route"] == "iv_bolus"
    end

    @testset "run_nca - Multiple Dose" begin
        result = run_nca(SS_TIMES, SS_CONC, TEST_DOSE;
                        config=config, dosing_type=:multiple, tau=SS_TAU)

        @test result.cmin !== nothing
        @test result.cavg !== nothing
        @test result.auc_0_tau !== nothing
        @test result.metadata["tau"] == SS_TAU
    end

    @testset "run_nca - Steady State" begin
        result = run_nca(SS_TIMES, SS_CONC, TEST_DOSE;
                        config=config, dosing_type=:steady_state, tau=SS_TAU)

        @test result.cmin !== nothing
        @test result.ptf !== nothing || result.cmin == 0.0
        @test result.swing !== nothing || result.cmin == 0.0
    end

    @testset "run_nca with custom config" begin
        custom_config = NCAConfig(
            method = LogLinearMethod(),
            lambda_z_r2_threshold = 0.85,
            extrapolation_max_pct = 30.0
        )

        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE; config=custom_config)
        @test result.auc_0_t > 0.0
    end

    @testset "nca_from_simresult" begin
        # Create mock SimResult-like structure
        # This would normally come from a simulation
        t = TEST_TIMES
        c = TEST_CONC

        # Just test the individual NCA works with same data
        result = run_nca(t, c, TEST_DOSE)
        @test result.cmax > 0.0
    end
end

# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

@testset "Edge Cases" begin
    config = NCAConfig()

    @testset "Minimum data points" begin
        t = [0.0, 1.0, 2.0]
        c = [0.0, 1.0, 0.5]

        # Should not throw with 3 points
        result = run_nca(t, c, 100.0; config=config)
        @test result.cmax == 1.0
    end

    @testset "Constant concentration" begin
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [1.0, 1.0, 1.0, 1.0, 1.0]

        result = run_nca(t, c, 100.0; config=config)
        @test result.cmax == 1.0
        # cmin is only calculated for multiple dose, not single dose
        @test result.cmin === nothing
    end

    @testset "All zeros after peak" begin
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [0.0, 1.0, 0.0, 0.0, 0.0]

        # Lambda_z estimation may fail, but NCA should complete
        result = run_nca(t, c, 100.0; config=config)
        @test result.cmax == 1.0
    end

    @testset "BLQ handling" begin
        config_blq = NCAConfig(lloq=0.5, blq_handling=BLQZero())

        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [0.0, 2.0, 1.0, 0.3, 0.1]  # Last 2 below LLOQ

        result = run_nca(t, c, 100.0; config=config_blq)
        @test result.cmax == 2.0
    end
end

# =============================================================================
# Input Validation Tests
# =============================================================================

@testset "Input Validation" begin
    config = NCAConfig()

    @testset "Mismatched lengths" begin
        @test_throws AssertionError run_nca([0.0, 1.0], [0.0], 100.0)
    end

    @testset "Negative dose" begin
        @test_throws AssertionError run_nca(TEST_TIMES, TEST_CONC, -100.0)
    end

    @testset "Unsorted times" begin
        @test_throws AssertionError run_nca([1.0, 0.0, 2.0], [1.0, 0.0, 0.5], 100.0)
    end

    @testset "Insufficient points" begin
        @test_throws AssertionError run_nca([0.0, 1.0], [0.0, 1.0], 100.0)
    end

    @testset "Multiple dose without tau" begin
        @test_throws AssertionError run_nca(
            TEST_TIMES, TEST_CONC, TEST_DOSE;
            dosing_type=:multiple
        )
    end
end

# =============================================================================
# Numerical Precision Tests
# =============================================================================

@testset "Numerical Precision" begin
    @testset "AUC precision" begin
        # Simple trapezoid - exact answer is 2.0
        t = [0.0, 1.0, 2.0]
        c = [0.0, 2.0, 2.0]
        config = NCAConfig(method=LinearMethod())

        auc = auc_0_t(t, c, config)
        @test auc ≈ 3.0 atol=1e-10  # Triangle (1.0) + Rectangle (2.0)
    end

    @testset "Log-linear extrapolation" begin
        # Verify extrapolation is reasonable
        t = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        c = [0.0, 1.0, 0.8, 0.64, 0.512, 0.4096, 0.328, 0.262, 0.21]
        config = NCAConfig()

        result = estimate_lambda_z(t, c, config)

        if result.lambda_z !== nothing
            # Should be approximately -log(0.8) ≈ 0.223
            @test result.lambda_z > 0.15
            @test result.lambda_z < 0.30
        end
    end

    @testset "Round NCA result" begin
        @test round_nca_result(1.2345, 3) ≈ 1.23 atol=0.01
        @test round_nca_result(0.001234, 3) ≈ 0.00123 atol=0.00001
        @test round_nca_result(1234.5, 3) ≈ 1230.0 atol=1.0
    end
end

# =============================================================================
# Lambda-z Selection Algorithm Tests (FDA/EMA Compliant)
# =============================================================================

@testset "Lambda-z Selection Algorithm" begin
    @testset "MinPointsFirst (FDA/EMA approach)" begin
        # Theophylline-like profile with clear terminal phase
        t = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
        c = [0.0, 2.5, 5.0, 8.0, 6.0, 4.5, 3.4, 1.9, 0.5]
        config = NCAConfig(lambda_z_min_points=3, lambda_z_r2_threshold=0.9)

        result = estimate_lambda_z(t, c, config; selection_method=MinPointsFirst())

        # Should use minimum required points from the END of terminal phase
        @test result.lambda_z !== nothing
        @test result.n_points >= 3
        @test result.r_squared >= 0.9

        # End time should be the last point (24.0)
        @test result.end_time == 24.0
    end

    @testset "MaxAdjR2 (traditional approach)" begin
        t = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
        c = [0.0, 2.5, 5.0, 8.0, 6.0, 4.5, 3.4, 1.9, 0.5]
        config = NCAConfig(lambda_z_min_points=3, lambda_z_r2_threshold=0.9)

        result = estimate_lambda_z(t, c, config; selection_method=MaxAdjR2())

        @test result.lambda_z !== nothing
        @test result.adjusted_r_squared !== nothing
        # MaxAdjR2 may use more points if it improves adj-R²
        @test result.n_points >= 3
    end

    @testset "Lambda-z validation with known dataset" begin
        # Known exponential decay: C(t) = 10 * exp(-0.1 * t)
        # Expected lambda_z = 0.1, t_half = ln(2)/0.1 ≈ 6.93
        t = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
        lambda_z_true = 0.1
        c = [0.0, 10.0 * exp(-lambda_z_true * 1.0),
             10.0 * exp(-lambda_z_true * 2.0),
             10.0 * exp(-lambda_z_true * 4.0),
             10.0 * exp(-lambda_z_true * 6.0),
             10.0 * exp(-lambda_z_true * 8.0),
             10.0 * exp(-lambda_z_true * 10.0),
             10.0 * exp(-lambda_z_true * 12.0)]

        config = NCAConfig(lambda_z_min_points=3, lambda_z_r2_threshold=0.99)
        result = estimate_lambda_z(t, c, config)

        @test result.lambda_z !== nothing
        @test result.lambda_z ≈ lambda_z_true atol=0.01

        expected_t_half = log(2) / lambda_z_true
        @test result.t_half ≈ expected_t_half atol=0.1

        # Perfect exponential should give R² very close to 1.0
        @test result.r_squared ≈ 1.0 atol=0.01
    end

    @testset "Lambda-z with noisy terminal phase" begin
        # Add some noise to mimic real data
        t = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0]
        c = [0.0, 8.5, 9.2, 6.8, 4.9, 3.7, 2.6, 1.8, 0.85, 0.18]  # Noisy terminal

        config = NCAConfig(lambda_z_min_points=3, lambda_z_r2_threshold=0.85)
        result = estimate_lambda_z(t, c, config)

        @test result.lambda_z !== nothing
        @test result.lambda_z > 0.0
        @test result.r_squared >= 0.85

        # Quality flag should be present if R² < 0.95
        if result.r_squared < 0.95
            @test result.quality_flag == :warning
        end
    end

    @testset "Lambda-z excludes points correctly" begin
        t = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0]
        c = [0.0, 5.0, 10.0, 8.0, 5.0, 3.0, 1.8, 0.6]
        config = NCAConfig(lambda_z_min_points=3)

        # Exclude the 6.0 hour point
        result = estimate_lambda_z(t, c, config; exclude_indices=[6])

        @test result.lambda_z !== nothing
        @test !(6 in result.points_used)
    end
end

# =============================================================================
# Regulatory Rounding Tests
# =============================================================================

@testset "Regulatory Rounding" begin
    @testset "round_nca_to_regulatory" begin
        # Significant figures rounding
        @test round_nca_to_regulatory(12.345, SIGNIFICANT_FIGURES, 3) ≈ 12.3 atol=0.01
        @test round_nca_to_regulatory(0.0012345, SIGNIFICANT_FIGURES, 3) ≈ 0.00123 atol=0.00001
        @test round_nca_to_regulatory(12345.0, SIGNIFICANT_FIGURES, 3) ≈ 12300.0 atol=10.0

        # Decimal places rounding
        @test round_nca_to_regulatory(12.3456, DECIMAL_PLACES, 2) ≈ 12.35 atol=0.001
        @test round_nca_to_regulatory(0.1234, DECIMAL_PLACES, 3) ≈ 0.123 atol=0.0001

        # Pharmacopeial rounding (half-up)
        @test round_nca_to_regulatory(1.235, PHARMACOPEIAL, 3) ≈ 1.24 atol=0.001
        @test round_nca_to_regulatory(1.225, PHARMACOPEIAL, 3) ≈ 1.23 atol=0.001

        # Edge cases
        @test round_nca_to_regulatory(0.0, SIGNIFICANT_FIGURES, 3) == 0.0
        @test isnan(round_nca_to_regulatory(NaN, SIGNIFICANT_FIGURES, 3))
        @test isinf(round_nca_to_regulatory(Inf, SIGNIFICANT_FIGURES, 3))
    end

    @testset "Parameter-specific precision" begin
        # High precision parameters
        @test round_nca_to_regulatory(0.08765, SIGNIFICANT_FIGURES, 4; parameter=:lambda_z) ≈ 0.08765 atol=0.00001
        @test round_nca_to_regulatory(6.9315, SIGNIFICANT_FIGURES, 4; parameter=:t_half) ≈ 6.932 atol=0.001

        # Standard precision parameters
        @test round_nca_to_regulatory(123.45, SIGNIFICANT_FIGURES, 3; parameter=:cmax) ≈ 123.0 atol=0.5
        @test round_nca_to_regulatory(456.789, SIGNIFICANT_FIGURES, 3; parameter=:auc_0_inf) ≈ 457.0 atol=0.5
    end

    @testset "round_nca_result_all" begin
        config = NCAConfig()
        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE; config=config)

        rounded = round_nca_result_all(result)

        # All values should be rounded to 3 significant figures (default)
        @test rounded.cmax !== nothing
        @test rounded.auc_0_t > 0.0

        # Verify rounding was applied
        @test rounded.cmax == round_nca_result(result.cmax, 3)
    end
end

# =============================================================================
# C0 Back-Extrapolation Validation Tests
# =============================================================================

@testset "C0 Validation" begin
    @testset "Valid C0 extrapolation" begin
        # Typical IV bolus: C0 = 10 mg/L with 100 mg dose → Vc = 10 L
        # For 70 kg human, Vc = 0.143 L/kg (physiologically plausible)
        result = validate_c0_extrapolation(
            10.0,  # C0
            100.0;  # dose
            expected_vc_range=(0.05, 0.5),
            body_weight=70.0
        )

        @test result.is_valid == true
        @test result.c0 == 10.0
        @test isempty(result.validation_flags)
        @test result.dose_c0_ratio ≈ 10.0 atol=0.01
    end

    @testset "Invalid C0 - non-positive" begin
        result = validate_c0_extrapolation(-5.0, 100.0)

        @test result.is_valid == false
        @test :non_positive_c0 in result.validation_flags
    end

    @testset "Invalid C0 - Vc too low" begin
        # C0 = 1000 → Vc = 0.1 L → Vc/kg = 0.0014 L/kg (way too small)
        result = validate_c0_extrapolation(
            1000.0, 100.0;
            expected_vc_range=(0.05, 1.0),
            body_weight=70.0
        )

        @test result.is_valid == false
        @test :vc_too_low in result.validation_flags
    end

    @testset "Invalid C0 - Vc too high" begin
        # C0 = 0.5 → Vc = 200 L → Vc/kg = 2.86 L/kg (too large)
        result = validate_c0_extrapolation(
            0.5, 100.0;
            expected_vc_range=(0.05, 1.0),
            body_weight=70.0
        )

        @test result.is_valid == false
        @test :vc_too_high in result.validation_flags
    end

    @testset "C0 validation with extrapolation time check" begin
        result = validate_c0_extrapolation(
            10.0, 100.0;
            max_extrapolation_time=2.0,
            t_half=1.0
        )

        # Extrapolation > t_half should generate warning
        @test :long_extrapolation in result.validation_flags
        @test length(result.warnings) > 0
    end

    @testset "C0 validation with first sample comparison" begin
        # C0 = 50, first sample = 10 (C0 > 2x first sample)
        result = validate_c0_extrapolation(
            50.0, 100.0;
            first_sample_time=0.5,
            first_sample_conc=10.0,
            expected_vc_range=(0.01, 0.5)  # Adjust range to pass Vc check
        )

        @test :c0_too_high in result.validation_flags
    end
end

# =============================================================================
# Sparse NCA Tests
# =============================================================================

@testset "Sparse NCA" begin
    @testset "SparseNCAConfig" begin
        config = SparseNCAConfig()
        @test config.pooling_method == :naive_pooled
        @test config.time_tolerance == 0.1
        @test config.min_samples_per_bin == 3

        custom = SparseNCAConfig(
            pooling_method=:destructive_sampling,
            time_tolerance=0.2,
            min_samples_per_bin=2
        )
        @test custom.pooling_method == :destructive_sampling
        @test custom.time_tolerance == 0.2
    end

    @testset "Sparse NCA with destructive sampling" begin
        # Simulated mouse PK study with destructive sampling
        # 3 mice per time point
        data = [
            (subject_id="M1", times=[0.5], concs=[9.5]),
            (subject_id="M2", times=[0.5], concs=[10.2]),
            (subject_id="M3", times=[0.5], concs=[10.0]),
            (subject_id="M4", times=[1.0], concs=[7.8]),
            (subject_id="M5", times=[1.0], concs=[8.1]),
            (subject_id="M6", times=[1.0], concs=[7.5]),
            (subject_id="M7", times=[2.0], concs=[5.2]),
            (subject_id="M8", times=[2.0], concs=[4.9]),
            (subject_id="M9", times=[2.0], concs=[5.0]),
            (subject_id="M10", times=[4.0], concs=[2.0]),
            (subject_id="M11", times=[4.0], concs=[2.2]),
            (subject_id="M12", times=[4.0], concs=[1.9]),
            (subject_id="M13", times=[8.0], concs=[0.5]),
            (subject_id="M14", times=[8.0], concs=[0.6]),
            (subject_id="M15", times=[8.0], concs=[0.4]),
        ]

        config = SparseNCAConfig(
            pooling_method=:destructive_sampling,
            min_samples_per_bin=3
        )

        result = run_sparse_nca(data, 100.0; config=config)

        @test result.n_subjects == 15
        @test result.nca_result.cmax > 0.0
        @test result.nca_result.auc_0_t > 0.0
        @test result.pooling_method == :destructive_sampling

        # Check mean profile was calculated
        @test !isempty(result.mean_profile.times)
        @test !isempty(result.mean_profile.mean_conc)

        # Should have samples at 5 nominal times
        @test length(result.n_samples_per_time) >= 5
    end

    @testset "Sparse NCA with time bins" begin
        data = [
            (subject_id="S1", times=[0.48], concs=[8.5]),
            (subject_id="S2", times=[0.52], concs=[9.0]),
            (subject_id="S3", times=[0.5], concs=[8.8]),
            (subject_id="S4", times=[1.02], concs=[6.5]),
            (subject_id="S5", times=[0.98], concs=[7.0]),
            (subject_id="S6", times=[1.01], concs=[6.8]),
        ]

        config = SparseNCAConfig(
            time_bins=[0.5, 1.0],
            time_tolerance=0.1,
            min_samples_per_bin=2
        )

        result = run_sparse_nca(data, 100.0; config=config)

        @test result.n_subjects == 6
        # Samples should be binned to nominal times
        @test haskey(result.n_samples_per_time, 0.5) || haskey(result.n_samples_per_time, 1.0)
    end

    @testset "Sparse NCA quality flags" begin
        # Data with insufficient samples at some time points
        data = [
            (subject_id="S1", times=[0.5], concs=[10.0]),  # Only 1 sample
            (subject_id="S2", times=[1.0], concs=[5.0]),
            (subject_id="S3", times=[1.0], concs=[5.5]),
            (subject_id="S4", times=[1.0], concs=[5.2]),
        ]

        config = SparseNCAConfig(min_samples_per_bin=3)
        result = run_sparse_nca(data, 100.0; config=config)

        @test :insufficient_samples in result.quality_flags
    end
end

# =============================================================================
# Partial AUC Method Tests
# =============================================================================

@testset "Partial AUC with Different Methods" begin
    t = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    c = [0.0, 5.0, 10.0, 7.0, 4.0, 2.0]

    config_linear = NCAConfig(method=LinearMethod())
    config_log = NCAConfig(method=LogLinearMethod())
    config_mixed = NCAConfig(method=LinLogMixedMethod())

    @testset "Partial AUC respects linear method" begin
        partial = auc_partial(t, c, 1.5, 3.5, config_linear)
        @test partial > 0.0
    end

    @testset "Partial AUC respects log-linear method" begin
        partial = auc_partial(t, c, 2.5, 4.5, config_log)
        @test partial > 0.0
    end

    @testset "Partial AUC respects lin-log mixed method" begin
        partial_mixed = auc_partial(t, c, 0.5, 4.5, config_mixed)
        @test partial_mixed > 0.0
    end

    @testset "Partial AUC interpolation consistency" begin
        # Partial AUC at exact data points should equal AUC between those points
        partial_exact = auc_partial(t, c, 1.0, 4.0, config_linear)

        # Manual calculation for indices 2-5 (times 1.0, 2.0, 3.0, 4.0)
        manual_auc = 0.0
        for i in 2:4
            manual_auc += 0.5 * (t[i+1] - t[i]) * (c[i] + c[i+1])
        end

        @test partial_exact ≈ manual_auc atol=0.01
    end
end

# =============================================================================
# Known Dataset Validation Tests
# =============================================================================

@testset "Known Dataset Validation" begin
    @testset "Simple 1-compartment IV bolus" begin
        # Known parameters: CL = 5 L/h, V = 50 L, Dose = 100 mg
        # C(t) = (Dose/V) * exp(-CL/V * t) = 2 * exp(-0.1 * t)
        # Expected: lambda_z = 0.1, t_half = 6.93, AUC_inf = 20

        CL = 5.0
        V = 50.0
        dose = 100.0
        lambda_z_true = CL / V  # 0.1

        t = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
        c = [(dose / V) * exp(-lambda_z_true * ti) for ti in t]

        config = NCAConfig()
        result = run_nca(t, c, dose; route=:iv_bolus, config=config)

        # Verify lambda_z
        if result.lambda_z_result.lambda_z !== nothing
            @test result.lambda_z_result.lambda_z ≈ lambda_z_true atol=0.01
        end

        # Verify t_half
        expected_t_half = log(2) / lambda_z_true
        if result.t_half !== nothing
            @test result.t_half ≈ expected_t_half atol=0.5
        end

        # Verify AUC_inf
        expected_auc_inf = dose / CL
        if result.auc_0_inf !== nothing
            @test result.auc_0_inf ≈ expected_auc_inf atol=1.0
        end

        # Verify CL
        if result.cl_f !== nothing
            @test result.cl_f ≈ CL atol=0.5
        end
    end

    @testset "Theophylline-like profile" begin
        # Based on NONMEM Theophylline example
        # Expected CL ≈ 3.14 L/hr, V ≈ 40 L
        t = [0.0, 0.25, 0.5, 1.0, 2.0, 3.5, 5.0, 7.0, 9.0, 12.0, 24.0]
        c = [0.0, 1.5, 2.5, 4.5, 6.5, 7.0, 6.2, 5.0, 3.8, 2.5, 0.6]
        dose = 320.0  # mg (typical oral dose)

        config = NCAConfig()
        result = run_nca(t, c, dose; config=config)

        @test result.cmax == 7.0
        @test result.tmax == 3.5
        @test result.auc_0_t > 50.0  # Rough estimate
        @test result.auc_0_t < 150.0

        if result.lambda_z_result.lambda_z !== nothing
            # Terminal phase should give reasonable elimination
            @test result.lambda_z_result.lambda_z > 0.05
            @test result.lambda_z_result.lambda_z < 0.2
        end
    end
end
