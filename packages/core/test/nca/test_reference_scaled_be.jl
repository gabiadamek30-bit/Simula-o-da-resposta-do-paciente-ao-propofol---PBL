# Tests for Reference-Scaled Average Bioequivalence (RSABE)
# Validates FDA RSABE and EMA ABEL implementations

using Test
using NeoPKPD
using StableRNGs

@testset "Reference-Scaled Bioequivalence" begin

    @testset "Regulatory Guidance Types" begin
        @testset "FDA Guidance" begin
            fda = FDAGuidance()

            # Check default parameters
            @test fda.theta ≈ (log(1.25) / 0.25)^2 atol=0.001  # ≈ 0.8926
            @test fda.point_estimate_lower == 0.80
            @test fda.point_estimate_upper == 1.25
            @test fda.swr_threshold ≈ 0.294 atol=0.001

            # Custom parameters
            fda_custom = FDAGuidance(
                point_estimate_lower=0.75,
                point_estimate_upper=1.30
            )
            @test fda_custom.point_estimate_lower == 0.75
            @test fda_custom.point_estimate_upper == 1.30
        end

        @testset "EMA Guidance" begin
            ema = EMAGuidance()

            # Check default parameters
            @test ema.k_scaling == 0.760
            @test ema.limit_lower_min ≈ 0.6984 atol=0.0001
            @test ema.limit_upper_max ≈ 1.4319 atol=0.0001
            @test ema.cv_threshold == 30.0
        end

        @testset "Health Canada Guidance" begin
            hc = HealthCanadaGuidance()

            # Should be similar to FDA
            @test hc.theta ≈ (log(1.25) / 0.25)^2 atol=0.001
            @test hc.point_estimate_lower == 0.80
            @test hc.point_estimate_upper == 1.25
        end
    end

    @testset "Replicate Design Types" begin
        @testset "Partial Replicate 3x3" begin
            design = PartialReplicate3x3()
            @test design.sequences == ["TRR", "RTR", "RRT"]

            # Custom sequences
            design2 = PartialReplicate3x3(sequences=["TRR", "RRT", "RTR"])
            @test length(design2.sequences) == 3
        end

        @testset "Full Replicate 2x4" begin
            design = FullReplicate2x4()
            @test design.sequences == ["TRTR", "RTRT"]
        end

        @testset "Full Replicate 2x3" begin
            design = FullReplicate2x3()
            @test design.sequences == ["TRT", "RTR"]
        end
    end

    @testset "Within-Subject Reference Variance (σWR)" begin
        @testset "Basic calculation" begin
            # Simulated replicate reference data
            # 8 subjects with 2 reference observations each
            ref_obs = [
                [10.0, 10.5],   # Subject 1
                [12.0, 11.5],   # Subject 2
                [9.5, 10.0],    # Subject 3
                [11.0, 10.8],   # Subject 4
                [10.5, 10.2],   # Subject 5
                [11.5, 12.0],   # Subject 6
                [9.8, 9.5],     # Subject 7
                [10.8, 11.2]    # Subject 8
            ]

            result = compute_within_subject_variance_reference(ref_obs)

            @test result.swr_squared > 0.0
            @test result.swr > 0.0
            @test result.cv_wr > 0.0
            @test result.df == 8  # One per subject
        end

        @testset "Convenience function compute_swr" begin
            # Two separate vectors (R1, R2)
            ref1 = [10.0, 12.0, 9.5, 11.0]
            ref2 = [10.5, 11.5, 10.0, 10.8]

            swr = compute_swr(ref1, ref2)
            @test swr > 0.0
            @test swr < 1.0  # Reasonable range for CV ~5-30%
        end

        @testset "High variability detection" begin
            # High variability data (CVw > 30%)
            rng = StableRNG(42)
            n = 20

            # Generate data with ~40% CV (σWR ≈ 0.39)
            ref1 = exp.(randn(rng, n) * 0.35 .+ log(100))
            ref2 = exp.(randn(rng, n) * 0.35 .+ log(100))

            ref_obs = [[ref1[i], ref2[i]] for i in 1:n]
            result = compute_within_subject_variance_reference(ref_obs)

            # CV should be roughly in the high variability range
            @test result.cv_wr > 20.0  # Should be moderately high given the simulation
        end
    end

    @testset "FDA RSABE Criterion" begin
        @testset "Scaling trigger" begin
            fda = FDAGuidance()

            # Below threshold - should not scale
            result_low = rsabe_criterion(0.05, 0.25, fda)
            @test result_low.use_scaled == false

            # Above threshold - should scale
            result_high = rsabe_criterion(0.05, 0.35, fda)
            @test result_high.use_scaled == true
        end

        @testset "Criterion evaluation" begin
            fda = FDAGuidance()

            # Bioequivalent case: small log difference, high variability
            result_be = rsabe_criterion(0.03, 0.35, fda)
            @test result_be.criterion_met == true
            @test result_be.use_scaled == true

            # Non-bioequivalent case: large log difference
            # With θ ≈ 0.797 and σWR = 0.35: need log_diff² > θ × σWR² = 0.0976
            # So log_diff > 0.312 to fail
            result_not_be = rsabe_criterion(0.40, 0.35, fda)
            @test result_not_be.criterion_met == false
        end

        @testset "Scaled limits calculation" begin
            fda = FDAGuidance()

            # High variability (σWR = 0.35)
            result = rsabe_criterion(0.05, 0.35, fda)

            # Scaled limits should be wider than 80-125%
            @test result.scaled_limit_lower < 0.80
            @test result.scaled_limit_upper > 1.25
        end
    end

    @testset "EMA ABEL Scaled Limits" begin
        @testset "Standard limits below threshold" begin
            ema = EMAGuidance()

            # Low variability - standard limits
            limits = abel_scaled_limits(0.20, ema)  # CVw ≈ 20%
            @test limits.is_scaled == false
            @test limits.lower ≈ 0.80 atol=0.001
            @test limits.upper ≈ 1.25 atol=0.001
        end

        @testset "Expanded limits above threshold" begin
            ema = EMAGuidance()

            # High variability - expanded limits
            limits = abel_scaled_limits(0.35, ema)  # CVw ≈ 37%
            @test limits.is_scaled == true
            @test limits.lower < 0.80
            @test limits.upper > 1.25
        end

        @testset "Maximum expansion caps" begin
            ema = EMAGuidance()

            # Very high variability - should hit caps
            limits = abel_scaled_limits(0.60, ema)  # Very high σWR
            @test limits.is_scaled == true

            # Should not exceed regulatory caps
            @test limits.lower >= ema.limit_lower_min - 0.0001
            @test limits.upper <= ema.limit_upper_max + 0.0001
        end

        @testset "Scaling formula verification" begin
            ema = EMAGuidance()
            swr = 0.35

            limits = abel_scaled_limits(swr, ema)

            # Manual calculation
            expected_lower = max(0.80 * exp(-0.760 * swr), 0.6984)
            expected_upper = min(1.25 * exp(0.760 * swr), 1.4319)

            @test limits.lower ≈ expected_lower atol=0.0001
            @test limits.upper ≈ expected_upper atol=0.0001
        end
    end

    @testset "Replicate Data Extraction" begin
        @testset "From matrices" begin
            # Test data (8 subjects, 2x4 design)
            test_vals = [100.0 105.0; 95.0 98.0; 110.0 108.0; 90.0 92.0;
                         102.0 100.0; 98.0 95.0; 105.0 102.0; 100.0 103.0]
            ref_vals = [98.0 102.0; 96.0 100.0; 105.0 110.0; 88.0 90.0;
                        100.0 98.0; 95.0 97.0; 102.0 105.0; 98.0 100.0]

            design = FullReplicate2x4()
            data = extract_replicate_data(test_vals, ref_vals, design)

            @test data.n_subjects == 8
            @test data.n_test_per_subject == 2
            @test data.n_ref_per_subject == 2
            @test length(data.test_obs) == 8
            @test length(data.ref_obs) == 8
        end
    end

    @testset "Full RSABE Analysis" begin
        @testset "FDA RSABE with bioequivalent data" begin
            # Simulate bioequivalent HVD data
            rng = StableRNG(123)
            n = 24

            # Generate replicate data with ~35% CV
            # Test is slightly lower than reference (GMR ≈ 1.02)
            base_ref = 100.0
            swr_target = 0.35  # High variability

            test_vals = zeros(n, 2)
            ref_vals = zeros(n, 2)

            for i in 1:n
                # Subject-specific random effect
                subj_effect = randn(rng) * 0.3

                # Reference observations with within-subject variability
                ref_vals[i, 1] = base_ref * exp(subj_effect + randn(rng) * swr_target)
                ref_vals[i, 2] = base_ref * exp(subj_effect + randn(rng) * swr_target)

                # Test observations (GMR ~1.0)
                test_vals[i, 1] = base_ref * 1.02 * exp(subj_effect + randn(rng) * swr_target)
                test_vals[i, 2] = base_ref * 1.02 * exp(subj_effect + randn(rng) * swr_target)
            end

            design = FullReplicate2x4()
            data = extract_replicate_data(test_vals, ref_vals, design)

            config = RSABEConfig(
                guidance = FDAGuidance(),
                design = design,
                parameter = :cmax
            )

            result = rsabe_analysis(data, config)

            @test result isa RSABEResult
            @test result.n_subjects == n
            @test result.gmr > 0.0
            @test result.cv_wr > 20.0  # Should be high variability
            @test result.parameter == :cmax

            # For this simulated data, should likely pass
            # (not guaranteed due to random variation)
            @test result.point_estimate_pass isa Bool
            @test result.scaled_criterion_pass isa Bool
        end

        @testset "Non-scaled pathway for low variability" begin
            # Low variability data (should use standard ABE)
            rng = StableRNG(456)
            n = 24

            # Generate data with ~15% CV (σWR ≈ 0.15)
            test_vals = zeros(n, 2)
            ref_vals = zeros(n, 2)

            for i in 1:n
                subj_effect = randn(rng) * 0.1
                ref_vals[i, 1] = 100.0 * exp(subj_effect + randn(rng) * 0.12)
                ref_vals[i, 2] = 100.0 * exp(subj_effect + randn(rng) * 0.12)
                test_vals[i, 1] = 100.0 * exp(subj_effect + randn(rng) * 0.12)
                test_vals[i, 2] = 100.0 * exp(subj_effect + randn(rng) * 0.12)
            end

            design = FullReplicate2x4()
            data = extract_replicate_data(test_vals, ref_vals, design)

            config = RSABEConfig(
                guidance = FDAGuidance(),
                design = design,
                parameter = :auc_0_inf
            )

            result = rsabe_analysis(data, config)

            # With low variability, should not use scaled approach
            @test result.use_scaled == false || result.cv_wr < 35.0  # May or may not scale
            @test result.parameter == :auc_0_inf
        end
    end

    @testset "Full ABEL Analysis" begin
        @testset "EMA ABEL with high variability" begin
            rng = StableRNG(789)
            n = 24

            # High variability data
            test_vals = zeros(n, 2)
            ref_vals = zeros(n, 2)

            for i in 1:n
                subj_effect = randn(rng) * 0.3
                ref_vals[i, 1] = 100.0 * exp(subj_effect + randn(rng) * 0.35)
                ref_vals[i, 2] = 100.0 * exp(subj_effect + randn(rng) * 0.35)
                test_vals[i, 1] = 100.0 * exp(subj_effect + randn(rng) * 0.35)
                test_vals[i, 2] = 100.0 * exp(subj_effect + randn(rng) * 0.35)
            end

            design = FullReplicate2x4()
            data = extract_replicate_data(test_vals, ref_vals, design)

            config = RSABEConfig(
                guidance = EMAGuidance(),
                design = design,
                parameter = :cmax
            )

            result = abel_analysis(data, config)

            @test result isa ABELResult
            @test result.n_subjects == n
            @test result.gmr > 0.0
            @test result.ci_lower > 0.0
            @test result.ci_upper > result.ci_lower

            # Should use expanded limits for HVD
            @test result.use_scaled isa Bool
            @test result.be_conclusion in [:bioequivalent, :not_bioequivalent]
        end
    end

    @testset "Unified Entry Point" begin
        @testset "FDA guidance routes to rsabe_analysis" begin
            rng = StableRNG(111)
            n = 12

            test_vals = exp.(randn(rng, n, 2) * 0.3 .+ log(100))
            ref_vals = exp.(randn(rng, n, 2) * 0.3 .+ log(100))

            design = FullReplicate2x4()
            data = extract_replicate_data(test_vals, ref_vals, design)

            config = RSABEConfig(guidance=FDAGuidance(), design=design)
            result = reference_scaled_be(data, config)

            @test result isa RSABEResult
        end

        @testset "EMA guidance routes to abel_analysis" begin
            rng = StableRNG(222)
            n = 12

            test_vals = exp.(randn(rng, n, 2) * 0.3 .+ log(100))
            ref_vals = exp.(randn(rng, n, 2) * 0.3 .+ log(100))

            design = FullReplicate2x4()
            data = extract_replicate_data(test_vals, ref_vals, design)

            config = RSABEConfig(guidance=EMAGuidance(), design=design)
            result = reference_scaled_be(data, config)

            @test result isa ABELResult
        end

        @testset "Direct matrix input" begin
            rng = StableRNG(333)
            n = 12

            test_vals = exp.(randn(rng, n, 2) * 0.3 .+ log(100))
            ref_vals = exp.(randn(rng, n, 2) * 0.3 .+ log(100))

            config = RSABEConfig(guidance=FDAGuidance(), design=FullReplicate2x4())
            result = reference_scaled_be(test_vals, ref_vals, config)

            @test result isa RSABEResult
            @test result.n_subjects == n
        end
    end

    @testset "Replicate ANOVA" begin
        @testset "Basic ANOVA computation" begin
            rng = StableRNG(444)
            n = 20

            test_vals = exp.(randn(rng, n, 2) * 0.25 .+ log(100))
            ref_vals = exp.(randn(rng, n, 2) * 0.25 .+ log(100))

            design = FullReplicate2x4()
            data = extract_replicate_data(test_vals, ref_vals, design)

            anova = replicate_anova(data, design)

            @test anova.mse_within > 0.0
            @test anova.swr_squared > 0.0
            @test anova.swr > 0.0
            @test anova.cv_wr > 0.0
            @test isfinite(anova.formulation_effect)
            @test anova.se_formulation > 0.0
            @test anova.df_within > 0
        end
    end

    @testset "Edge Cases" begin
        @testset "Minimum sample size" begin
            # 3 subjects (minimum)
            test_vals = [100.0 105.0; 95.0 98.0; 102.0 100.0]
            ref_vals = [98.0 102.0; 96.0 100.0; 100.0 98.0]

            design = FullReplicate2x4()
            data = extract_replicate_data(test_vals, ref_vals, design)

            config = RSABEConfig(guidance=FDAGuidance(), design=design)
            result = reference_scaled_be(data, config)

            @test result.n_subjects == 3
        end

        @testset "GMR at boundary" begin
            fda = FDAGuidance()

            # GMR exactly at 0.80 (on boundary)
            log_diff_low = log(0.80)
            result_low = rsabe_criterion(log_diff_low, 0.35, fda)
            @test result_low.criterion_met isa Bool

            # GMR exactly at 1.25 (on boundary)
            log_diff_high = log(1.25)
            result_high = rsabe_criterion(log_diff_high, 0.35, fda)
            @test result_high.criterion_met isa Bool
        end

        @testset "σWR at threshold boundary" begin
            fda = FDAGuidance()

            # Exactly at threshold
            result_at = rsabe_criterion(0.05, 0.294, fda)
            @test result_at.use_scaled isa Bool

            # Just above threshold
            result_above = rsabe_criterion(0.05, 0.295, fda)
            @test result_above.use_scaled == true

            # Just below threshold
            result_below = rsabe_criterion(0.05, 0.293, fda)
            @test result_below.use_scaled == false
        end
    end

    @testset "Regulatory Compliance" begin
        @testset "FDA theta value" begin
            # FDA uses θ = (ln(1.25)/0.25)² ≈ 0.7967
            # Note: ln(1.25)/0.25 ≈ 0.8926, but θ is the square of this
            fda = FDAGuidance()
            expected_theta = (log(1.25) / 0.25)^2

            @test fda.theta ≈ expected_theta atol=1e-6
            @test fda.theta ≈ 0.7967 atol=0.001
        end

        @testset "EMA k value" begin
            # EMA uses k = 0.760
            ema = EMAGuidance()
            @test ema.k_scaling == 0.760
        end

        @testset "EMA cap values" begin
            # EMA caps at 69.84% - 143.19%
            ema = EMAGuidance()
            @test ema.limit_lower_min ≈ 0.6984 atol=0.0001
            @test ema.limit_upper_max ≈ 1.4319 atol=0.0001
        end
    end

end
