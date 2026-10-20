# Exposure-Response Test Suite
# Tests for E-R analysis module

using Test
using NeoPKPD
using Statistics
using Random

@testset "Exposure Metric Types" begin
    @testset "CmaxMetric" begin
        metric = CmaxMetric()
        @test metric.log_transform == false

        metric_log = CmaxMetric(log_transform = true)
        @test metric_log.log_transform == true
    end

    @testset "AUCMetric" begin
        metric = AUCMetric()
        @test metric.type == :auc_0_inf

        metric_tau = AUCMetric(type = :auc_0_tau)
        @test metric_tau.type == :auc_0_tau

        @test_throws AssertionError AUCMetric(type = :invalid)
    end

    @testset "CtroughMetric" begin
        metric = CtroughMetric(log_transform = true)
        @test metric.log_transform == true
    end

    @testset "CavgMetric" begin
        metric = CavgMetric()
        @test metric.log_transform == false
    end
end

@testset "Response Type Types" begin
    @testset "BinaryResponse" begin
        resp = BinaryResponse()
        @test resp.success_label == :responder

        resp2 = BinaryResponse(success_label = :success)
        @test resp2.success_label == :success
    end

    @testset "ContinuousResponse" begin
        resp = ContinuousResponse()
        @test resp.baseline_corrected == true
    end

    @testset "CountResponse" begin
        resp = CountResponse()
        @test resp.offset_variable === nothing
    end
end

@testset "E-R Models" begin
    @testset "LogisticERModel initialization" begin
        model = LogisticERModel()
        @test model.intercept === nothing
        @test model.slope === nothing
        @test isempty(model.covariates)
    end

    @testset "LinearERModel initialization" begin
        model = LinearERModel()
        @test model.intercept === nothing
        @test model.residual_variance === nothing
    end

    @testset "EmaxERModel initialization" begin
        model = EmaxERModel()
        @test model.e0 === nothing
        @test model.emax === nothing
        @test model.ec50 === nothing
    end
end

@testset "ExposureResponseConfig" begin
    @testset "Default configuration" begin
        config = ExposureResponseConfig()
        @test config.exposure_metric isa AUCMetric
        @test config.response_type isa BinaryResponse
        @test config.model isa LogisticERModel
        @test config.confidence_level == 0.95
    end

    @testset "Custom configuration" begin
        config = ExposureResponseConfig(
            exposure_metric = CmaxMetric(log_transform = true),
            response_type = ContinuousResponse(),
            model = LinearERModel(),
            confidence_level = 0.90,
            n_bootstrap = 500
        )

        @test config.exposure_metric isa CmaxMetric
        @test config.exposure_metric.log_transform == true
        @test config.response_type isa ContinuousResponse
        @test config.confidence_level == 0.90
    end

    @testset "Invalid confidence level" begin
        @test_throws AssertionError ExposureResponseConfig(confidence_level = 1.5)
        @test_throws AssertionError ExposureResponseConfig(confidence_level = 0.0)
    end
end

@testset "Logistic E-R Model Fitting" begin
    Random.seed!(42)

    # Generate data with known relationship
    n = 200
    exposure = rand(n) * 100  # Exposure 0-100

    # True logistic: P(Y=1) = logistic(-2 + 0.04*exposure)
    true_intercept = -2.0
    true_slope = 0.04
    prob = 1.0 ./ (1.0 .+ exp.(-(true_intercept .+ true_slope .* exposure)))
    response = Float64.(rand(n) .< prob)

    config = ExposureResponseConfig(
        model = LogisticERModel(),
        n_bootstrap = 100,
        seed = 12345
    )

    result = run_exposure_response(exposure, response, config)

    @testset "Basic results" begin
        @test result.n_subjects == n
        @test length(result.exposure_values) == n
        @test length(result.predicted_responses) == n
        @test all(0 .<= result.predicted_responses .<= 1)
    end

    @testset "Parameter estimates" begin
        @test haskey(result.parameter_estimates, :intercept)
        @test haskey(result.parameter_estimates, :slope)

        # Should be close to true values
        @test isapprox(result.parameter_estimates[:intercept], true_intercept, atol = 0.5)
        @test isapprox(result.parameter_estimates[:slope], true_slope, atol = 0.02)
    end

    @testset "Standard errors" begin
        @test haskey(result.parameter_se, :intercept)
        @test haskey(result.parameter_se, :slope)
        @test result.parameter_se[:slope] > 0
    end

    @testset "Odds ratios" begin
        @test haskey(result.odds_ratios, :per_unit)
        @test result.odds_ratios[:per_unit] > 1  # Positive slope
    end

    @testset "Model fit statistics" begin
        @test haskey(result.model_fit, :log_likelihood)
        @test haskey(result.model_fit, :mcfadden_r2)
        @test haskey(result.model_fit, :aic)
        @test haskey(result.model_fit, :brier_score)

        @test result.model_fit[:mcfadden_r2] >= 0
        @test result.model_fit[:brier_score] <= 0.25  # Should be better than random
    end

    @testset "Quantile predictions" begin
        @test haskey(result.exposure_quantiles, 0.5)
        @test haskey(result.response_at_quantiles, 0.5)
    end
end

@testset "Linear E-R Model Fitting" begin
    Random.seed!(42)

    # Generate data with known relationship
    n = 100
    exposure = rand(n) * 50  # Exposure 0-50

    # True linear: Y = 10 + 0.5*exposure + error
    true_intercept = 10.0
    true_slope = 0.5
    response = true_intercept .+ true_slope .* exposure .+ randn(n) * 2.0

    config = ExposureResponseConfig(
        response_type = ContinuousResponse(),
        model = LinearERModel(),
        n_bootstrap = 100
    )

    result = run_exposure_response(exposure, response, config)

    @testset "Parameter estimates" begin
        @test isapprox(result.parameter_estimates[:intercept], true_intercept, atol = 2.0)
        @test isapprox(result.parameter_estimates[:slope], true_slope, atol = 0.2)
        @test result.parameter_estimates[:residual_variance] > 0
    end

    @testset "Model fit" begin
        @test haskey(result.model_fit, :r_squared)
        @test haskey(result.model_fit, :rmse)
        @test result.model_fit[:r_squared] > 0.5  # Should explain significant variance
    end
end

@testset "Emax E-R Model Fitting" begin
    Random.seed!(42)

    # Generate Emax data
    n = 100
    exposure = rand(n) * 100

    # True Emax: Y = 5 + 50 * x / (20 + x) + error
    true_e0 = 5.0
    true_emax = 50.0
    true_ec50 = 20.0
    response = true_e0 .+ true_emax .* exposure ./ (true_ec50 .+ exposure) .+ randn(n) * 3.0

    config = ExposureResponseConfig(
        response_type = ContinuousResponse(),
        model = EmaxERModel(),
        n_bootstrap = 50
    )

    result = run_exposure_response(exposure, response, config)

    @testset "Parameter estimates" begin
        @test haskey(result.parameter_estimates, :e0)
        @test haskey(result.parameter_estimates, :emax)
        @test haskey(result.parameter_estimates, :ec50)

        @test result.parameter_estimates[:ec50] > 0
    end
end

@testset "Prediction Functions" begin
    Random.seed!(42)

    # Generate cleaner data with clear positive relationship
    n = 50
    exposure = collect(range(10.0, 100.0, length=n))
    true_intercept = -2.0
    true_slope = 0.04
    prob = 1.0 ./ (1.0 .+ exp.(-(true_intercept .+ true_slope .* exposure)))
    response = Float64.(rand(n) .< prob)

    config = ExposureResponseConfig(model = LogisticERModel(), n_bootstrap = 50)
    result = run_exposure_response(exposure, response, config)

    @testset "predict_response_at_dose" begin
        pred50 = predict_response_at_dose(result, 50.0)
        @test 0 <= pred50 <= 1
    end

    @testset "find_target_exposure" begin
        # Find exposure for 50% response probability
        target_exp = find_target_exposure(result, 0.5)

        # With a positive slope model, target should be finite
        @test isfinite(target_exp)

        # Verify prediction at found exposure
        pred = predict_response_at_dose(result, target_exp)
        @test isapprox(pred, 0.5, atol = 0.02)
    end
end

@testset "Exposure Quantile Analysis" begin
    Random.seed!(42)

    n = 100
    exposure = rand(n) * 100
    prob = 0.2 .+ 0.6 .* (exposure ./ 100)
    response = Float64.(rand(n) .< prob)

    result = compute_exposure_quantiles(exposure, response; n_groups = 4)

    @testset "Group structure" begin
        @test result[:n_groups] == 4
        @test length(result[:boundaries]) == 5
        @test length(result[:group_stats]) == 4
    end

    @testset "Response trend" begin
        # Higher exposure should have higher response rate
        rates = [s.response_rate for s in result[:group_stats]]
        # Not strictly monotonic due to randomness, but overall trend should be positive
        @test rates[4] > rates[1]  # Highest quartile > lowest quartile
    end
end

@testset "Bootstrap CI" begin
    Random.seed!(42)

    n = 100
    exposure = rand(n) * 50
    response = Float64.(exposure .> 25)  # Perfect separation at 25

    config = ExposureResponseConfig(
        model = LogisticERModel(),
        n_bootstrap = 100,
        confidence_level = 0.90
    )

    result = run_exposure_response(exposure, response, config)

    @testset "CI structure" begin
        if !isempty(result.parameter_ci)
            for (param, ci) in result.parameter_ci
                @test ci[1] < ci[2]  # Lower < upper
            end
        end
    end
end

@testset "Missing Data Handling" begin
    exposure = [10.0, 20.0, NaN, 40.0, 50.0]
    response = [0.0, 1.0, 1.0, NaN, 1.0]

    config = ExposureResponseConfig(n_bootstrap = 10)
    result = run_exposure_response(exposure, response, config)

    # Should handle missing values gracefully
    @test result.n_subjects == 3  # Only 3 valid pairs
    @test length(result.exposure_values) == 3
end
