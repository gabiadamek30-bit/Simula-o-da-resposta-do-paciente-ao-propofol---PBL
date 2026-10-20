# Probability of Success Test Suite
# Tests for PoS/Assurance calculations

using Test
using NeoPKPD
using Statistics
using Random

@testset "Prior Types" begin
    @testset "NormalPoSPrior" begin
        prior = NormalPoSPrior(0.3, 0.04)
        @test prior.mean == 0.3
        @test prior.variance == 0.04
        @test prior.source == :expert

        prior2 = NormalPoSPrior(0.5, 0.1; source = :meta_analysis)
        @test prior2.source == :meta_analysis

        @test_throws AssertionError NormalPoSPrior(0.3, -0.1)  # Negative variance
    end

    @testset "BetaPoSPrior" begin
        prior = BetaPoSPrior(10.0, 5.0)
        @test prior.alpha == 10.0
        @test prior.beta == 5.0

        @test_throws AssertionError BetaPoSPrior(0.0, 5.0)  # Zero alpha
    end

    @testset "HistoricalDataPrior" begin
        effects = [0.25, 0.30, 0.28]
        ses = [0.05, 0.06, 0.04]

        prior = HistoricalDataPrior(effects, ses; discount_factor = 0.5)
        @test length(prior.effect_estimates) == 3
        @test prior.discount_factor == 0.5
        @test prior.source == :historical

        # Convert to normal
        normal_prior = to_normal_prior(prior)
        @test normal_prior isa NormalPoSPrior
        @test 0.25 < normal_prior.mean < 0.32  # Should be pooled estimate
    end
end

@testset "Success Criteria" begin
    @testset "SuperioritySuccess" begin
        crit = SuperioritySuccess()
        @test crit.margin == 0.0
        @test crit.alpha == 0.025
        @test crit.one_sided == true

        crit2 = SuperioritySuccess(margin = 0.1, alpha = 0.05)
        @test crit2.margin == 0.1
    end

    @testset "NonInferioritySuccess" begin
        crit = NonInferioritySuccess(margin = 0.15)
        @test crit.margin == 0.15

        @test_throws AssertionError NonInferioritySuccess(margin = -0.1)
    end

    @testset "EquivalenceSuccess" begin
        crit = EquivalenceSuccess(lower_margin = -0.2, upper_margin = 0.2)
        @test crit.lower_margin == -0.2
        @test crit.upper_margin == 0.2
    end
end

@testset "PoSConfig" begin
    prior = NormalPoSPrior(0.3, 0.04)
    criterion = SuperioritySuccess()

    config = PoSConfig(
        prior = prior,
        criterion = criterion,
        n_subjects_planned = 200,
        n_samples = 1000
    )

    @test config.n_subjects_planned == 200
    @test config.n_samples == 1000
    @test config.include_uncertainty == true
end

@testset "Power Calculations" begin
    @testset "Superiority power" begin
        crit = SuperioritySuccess(margin = 0.0, alpha = 0.025)

        # High effect size should have high power
        power_high = compute_power(0.5, 0.1, crit)
        @test power_high > 0.95

        # Low effect size should have low power
        power_low = compute_power(0.1, 0.1, crit)
        @test power_low < 0.5

        # Zero effect should have ~alpha
        power_zero = compute_power(0.0, 0.1, crit)
        @test 0.02 < power_zero < 0.05
    end

    @testset "Non-inferiority power" begin
        crit = NonInferioritySuccess(margin = 0.1, alpha = 0.025)

        # Effect at -margin should have ~alpha power (null boundary)
        power_at_null = compute_power(-0.1, 0.1, crit)
        @test 0.02 < power_at_null < 0.05

        # Effect at 0 has limited power with high SE
        # (effect+margin)/se = 0.1/0.1 = 1, need ~2 to have high power
        power_zero = compute_power(0.0, 0.1, crit)
        @test power_zero > 0.1  # Modest power

        # Large positive effect with smaller SE should have high power
        power_positive = compute_power(0.3, 0.05, crit)  # Larger effect, smaller SE
        @test power_positive > 0.9
    end
end

@testset "Assurance Calculation" begin
    @testset "Basic assurance" begin
        prior = NormalPoSPrior(0.3, 0.04)
        config = PoSConfig(
            prior = prior,
            criterion = SuperioritySuccess(),
            n_subjects_planned = 200,
            n_samples = 5000,
            seed = 12345
        )

        assurance = compute_assurance(config)

        @test 0 <= assurance.assurance <= 1
        @test assurance.assurance_ci[1] < assurance.assurance_ci[2]
        @test assurance.n_samples == 5000
        @test isfinite(assurance.predictive_mean)
        @test isfinite(assurance.predictive_sd)
        @test haskey(assurance.predictive_quantiles, 0.5)
    end

    @testset "Strong prior should give high assurance" begin
        # Strong effect, low variance
        prior = NormalPoSPrior(0.5, 0.01)  # Very confident about 0.5 effect
        config = PoSConfig(
            prior = prior,
            criterion = SuperioritySuccess(),
            n_subjects_planned = 400,
            n_samples = 2000,
            seed = 12345
        )

        assurance = compute_assurance(config)
        @test assurance.assurance > 0.8
    end

    @testset "Weak prior should give low assurance" begin
        # Effect near zero, high variance
        prior = NormalPoSPrior(0.0, 0.25)
        config = PoSConfig(
            prior = prior,
            criterion = SuperioritySuccess(),
            n_subjects_planned = 100,
            n_samples = 2000,
            seed = 12345
        )

        assurance = compute_assurance(config)
        @test assurance.assurance < 0.5
    end
end

@testset "Complete PoS Analysis" begin
    prior = NormalPoSPrior(0.25, 0.04; source = :previous_trial)
    criterion = SuperioritySuccess(margin = 0.0, alpha = 0.025)

    config = PoSConfig(
        prior = prior,
        criterion = criterion,
        n_subjects_planned = 300,
        n_samples = 2000,
        seed = 42
    )

    result = compute_pos(config; run_sensitivity = false)

    @testset "Result structure" begin
        @test result.config === config
        @test result.assurance isa AssuranceResult
        @test length(result.posterior_samples) == 2000
        @test !isempty(result.power_curve)
        @test result.go_nogo_recommendation in [:go, :nogo, :uncertain]
        @test !isempty(result.decision_rationale)
    end

    @testset "Power curve" begin
        # Power should increase with effect size
        effects = sort(collect(keys(result.power_curve)))
        powers = [result.power_curve[e] for e in effects]

        # Check overall increasing trend
        @test powers[end] > powers[1]
    end
end

@testset "Sensitivity Analysis" begin
    prior = NormalPoSPrior(0.3, 0.04)
    config = PoSConfig(
        prior = prior,
        criterion = SuperioritySuccess(),
        n_subjects_planned = 200,
        n_samples = 1000,
        seed = 12345
    )

    sensitivity = sensitivity_analysis_pos(config)

    @test !isempty(sensitivity)
    @test haskey(sensitivity, "Prior mean 80%")
    @test haskey(sensitivity, "Prior mean 120%")
    @test haskey(sensitivity, "Sample size 75%")
    @test haskey(sensitivity, "Sample size 125%")

    # All results should be valid AssuranceResult
    for (label, result) in sensitivity
        @test result isa AssuranceResult
        @test 0 <= result.assurance <= 1
    end
end

@testset "Prior from Historical Trials" begin
    # Three historical trials
    effects = [0.22, 0.28, 0.25]
    standard_errors = [0.06, 0.05, 0.07]

    prior = create_pos_prior_from_trials(effects, standard_errors; discount = 0.5)

    @test prior isa NormalPoSPrior
    @test prior.source == :historical

    # Mean should be close to weighted average
    @test 0.22 < prior.mean < 0.28

    # Variance should be inflated due to discount
    @test prior.variance > 0
end

@testset "Recommendation Logic" begin
    criterion = SuperioritySuccess()

    @testset "High assurance -> Go" begin
        # Very confident prior
        prior = NormalPoSPrior(0.6, 0.02)
        config = PoSConfig(
            prior = prior,
            criterion = criterion,
            n_subjects_planned = 500,
            n_samples = 2000
        )

        result = compute_pos(config)
        @test result.go_nogo_recommendation == :go
    end

    @testset "Low assurance -> NoGo" begin
        # Poor prior
        prior = NormalPoSPrior(-0.1, 0.1)
        config = PoSConfig(
            prior = prior,
            criterion = criterion,
            n_subjects_planned = 100,
            n_samples = 2000
        )

        result = compute_pos(config)
        @test result.go_nogo_recommendation == :nogo
    end
end

@testset "Decision Boundary" begin
    prior = NormalPoSPrior(0.3, 0.04)
    config = PoSConfig(
        prior = prior,
        criterion = SuperioritySuccess(),
        n_subjects_planned = 200,
        n_samples = 500,  # Fewer for speed
        seed = 12345
    )

    min_n = pos_decision_boundary(config; target_assurance = 0.80)

    @test min_n > 0
    @test min_n < 2000  # Should find a solution
end

@testset "Predictive Probability at Interim" begin
    prior = NormalPoSPrior(0.3, 0.04)
    config = PoSConfig(
        prior = prior,
        criterion = SuperioritySuccess(),
        n_subjects_planned = 200,
        n_samples = 2000,
        seed = 12345
    )

    # Simulate interim data
    current_data = Dict{Symbol, Float64}(
        :n_current => 100.0,
        :effect_current => 0.35,
        :se_current => 0.08
    )

    pp = compute_predictive_probability(current_data, 200, config)

    @test 0 <= pp <= 1

    # Good interim result should give reasonable PP
    @test pp > 0.3
end

@testset "Summary Report" begin
    prior = NormalPoSPrior(0.3, 0.04)
    config = PoSConfig(
        prior = prior,
        criterion = SuperioritySuccess(),
        n_subjects_planned = 200,
        n_samples = 1000
    )

    result = compute_pos(config)
    summary = summarize_pos_result(result)

    @test summary isa String
    @test occursin("Probability of Success", summary)
    @test occursin("Assurance", summary)
    @test occursin("Recommendation", summary)
end
