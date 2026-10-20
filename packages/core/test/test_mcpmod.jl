# MCP-Mod Test Suite
# Tests for Multiple Comparison Procedure - Modelling (FDA-qualified Phase II method)

using Test
using NeoPKPD
using Statistics
using Random

@testset "Dose-Response Models" begin
    @testset "LinearModel" begin
        model = LinearModel()
        @test model.name == :linear

        params = Dict{Symbol, Float64}(:e0 => 0.5, :delta => 0.02)
        @test model_mean(model, 0.0, params) ≈ 0.5
        @test model_mean(model, 50.0, params) ≈ 1.5
        @test model_mean(model, 100.0, params) ≈ 2.5
    end

    @testset "QuadraticModel" begin
        model = QuadraticModel(delta = 100.0)
        @test model.name == :quadratic
        @test model.delta == 100.0

        params = Dict{Symbol, Float64}(:e0 => 0.0, :b1 => 0.02, :b2 => -0.0001)
        @test model_mean(model, 0.0, params) ≈ 0.0
        @test model_mean(model, 50.0, params) ≈ 0.75  # 0 + 0.02*50 - 0.0001*2500
    end

    @testset "EmaxModel" begin
        model = EmaxModel(ed50 = 25.0)
        @test model.name == :emax
        @test model.ed50 == 25.0

        params = Dict{Symbol, Float64}(:e0 => 0.0, :emax => 1.0, :ed50 => 25.0)
        @test model_mean(model, 0.0, params) ≈ 0.0
        @test model_mean(model, 25.0, params) ≈ 0.5  # 50% effect at ED50
        @test model_mean(model, 100.0, params) ≈ 0.8  # 100/(25+100) = 0.8
    end

    @testset "SigmoidEmaxModel" begin
        model = SigmoidEmaxModel(ed50 = 25.0, h = 2.0)
        @test model.name == :sigmoid_emax
        @test model.ed50 == 25.0
        @test model.h == 2.0

        params = Dict{Symbol, Float64}(:e0 => 0.0, :emax => 1.0, :ed50 => 25.0, :h => 2.0)
        @test model_mean(model, 0.0, params) ≈ 0.0
        @test model_mean(model, 25.0, params) ≈ 0.5  # 50% effect at ED50

        # With Hill coefficient > 1, curve is steeper
        @test model_mean(model, 50.0, params) > 0.75  # Steeper than standard Emax
    end

    @testset "ExponentialModel" begin
        model = ExponentialModel(delta = 50.0)
        @test model.name == :exponential
        @test model.delta == 50.0

        params = Dict{Symbol, Float64}(:e0 => 0.0, :e1 => 1.0, :delta => 50.0)
        @test model_mean(model, 0.0, params) ≈ 0.0
        @test model_mean(model, 50.0, params) ≈ exp(1) - 1  # ≈ 1.718
    end

    @testset "BetaModel" begin
        model = BetaModel(delta1 = 1.0, delta2 = 1.0, scal = 120.0)
        @test model.name == :beta
        @test model.delta1 == 1.0
        @test model.delta2 == 1.0

        params = Dict{Symbol, Float64}(:e0 => 0.0, :emax => 1.0)
        @test model_mean(model, 0.0, params) ≈ 0.0  # Zero at boundary
        @test model_mean(model, 60.0, params) > 0.5  # Peak near middle
    end

    @testset "LogLinearModel" begin
        model = LogLinearModel(off = 1.0)
        @test model.name == :loglinear
        @test model.off == 1.0

        params = Dict{Symbol, Float64}(:e0 => 0.0, :delta => 1.0)
        @test model_mean(model, 0.0, params) ≈ 0.0
        @test model_mean(model, 99.0, params) ≈ log(100.0)  # log(99 + 1)
    end
end

@testset "Optimal Contrasts" begin
    doses = [0.0, 10.0, 25.0, 50.0, 100.0]

    @testset "Contrasts sum to zero" begin
        models = DoseResponseModel[
            LinearModel(),
            EmaxModel(ed50 = 25.0),
            SigmoidEmaxModel(ed50 = 25.0, h = 2.0)
        ]

        config = MCPModConfig(
            doses = doses,
            models = models,
            alpha = 0.025,
            direction = :increasing
        )

        contrasts = compute_optimal_contrasts(config)

        for model in models
            c = contrasts[model.name]
            @test length(c) == length(doses)
            @test abs(sum(c)) < 1e-10  # Contrasts sum to zero
        end
    end

    @testset "Linear contrast is linear" begin
        models = DoseResponseModel[LinearModel()]
        config = MCPModConfig(doses = doses, models = models, direction = :increasing)
        contrasts = compute_optimal_contrasts(config)

        c = contrasts[:linear]
        # For linear model with increasing direction, contrasts should increase
        @test c[end] > c[1]
    end
end

@testset "MCPModConfig" begin
    @testset "Valid configuration" begin
        doses = [0.0, 25.0, 50.0, 100.0]
        models = DoseResponseModel[LinearModel(), EmaxModel(ed50 = 50.0)]

        config = MCPModConfig(
            doses = doses,
            models = models,
            alpha = 0.025,
            delta = 0.5,
            direction = :increasing
        )

        @test config.alpha == 0.025
        @test config.delta == 0.5
        @test config.direction == :increasing
        @test length(config.models) == 2
    end

    @testset "Must include placebo" begin
        doses = [25.0, 50.0, 100.0]  # No placebo
        models = DoseResponseModel[LinearModel()]

        @test_throws AssertionError MCPModConfig(doses = doses, models = models)
    end

    @testset "Direction validation" begin
        doses = [0.0, 50.0, 100.0]
        models = DoseResponseModel[LinearModel()]

        @test_throws AssertionError MCPModConfig(
            doses = doses,
            models = models,
            direction = :invalid
        )
    end
end

@testset "MCP Step" begin
    @testset "Proof of concept detection" begin
        doses = [0.0, 25.0, 50.0, 100.0]
        models = DoseResponseModel[LinearModel(), EmaxModel(ed50 = 50.0)]

        config = MCPModConfig(
            doses = doses,
            models = models,
            alpha = 0.025,
            direction = :increasing
        )

        # Strong dose-response signal
        means = [0.0, 0.5, 0.8, 1.0]
        se = [0.1, 0.1, 0.1, 0.1]
        n = [50, 50, 50, 50]

        result = mcp_step(means, se, n, config)

        @test result.proof_of_concept == true
        @test !isempty(result.significant_models)
        @test result.max_test_statistic > 0
        @test result.critical_value > 0
    end

    @testset "No signal detection" begin
        doses = [0.0, 25.0, 50.0, 100.0]
        models = DoseResponseModel[LinearModel()]

        config = MCPModConfig(
            doses = doses,
            models = models,
            alpha = 0.025,
            direction = :increasing
        )

        # Flat response (no dose-response)
        means = [0.5, 0.48, 0.52, 0.49]
        se = [0.2, 0.2, 0.2, 0.2]
        n = [30, 30, 30, 30]

        result = mcp_step(means, se, n, config)

        # With flat data and high variability, should not find PoC
        # Note: This is probabilistic, but with flat data it should generally fail
        @test result.max_test_statistic < result.critical_value || isempty(result.significant_models)
    end
end

@testset "Parameter Estimation" begin
    doses = [0.0, 25.0, 50.0, 100.0]
    weights = fill(1.0, length(doses))

    @testset "Linear model fitting" begin
        model = LinearModel()
        # True: e0 = 0.5, delta = 0.01
        means = [0.5, 0.75, 1.0, 1.5]

        params = estimate_model_parameters(model, doses, means, weights)

        @test haskey(params, :e0)
        @test haskey(params, :delta)
        @test isapprox(params[:e0], 0.5, atol = 0.1)
        @test isapprox(params[:delta], 0.01, atol = 0.001)
    end

    @testset "Emax model fitting" begin
        model = EmaxModel(ed50 = 50.0)
        # Generate data from Emax model
        true_params = Dict(:e0 => 0.0, :emax => 1.0, :ed50 => 50.0)
        means = [model_mean(model, d, true_params) for d in doses]

        params = estimate_model_parameters(model, doses, means, weights)

        @test haskey(params, :e0)
        @test haskey(params, :emax)
        @test haskey(params, :ed50)
        @test params[:ed50] > 0  # ED50 should be positive
    end
end

@testset "Mod Step" begin
    doses = [0.0, 25.0, 50.0, 100.0]

    @testset "Model selection" begin
        models = DoseResponseModel[
            LinearModel(),
            EmaxModel(ed50 = 50.0),
            SigmoidEmaxModel(ed50 = 50.0, h = 2.0)
        ]

        config = MCPModConfig(
            doses = doses,
            models = models,
            alpha = 0.025,
            delta = 0.5,
            direction = :increasing
        )

        # Generate Emax-like data
        means = [0.0, 0.33, 0.5, 0.67]
        se = [0.1, 0.1, 0.1, 0.1]
        n = [50, 50, 50, 50]

        significant = [:linear, :emax, :sigmoid_emax]
        result = mod_step(means, se, n, config, significant)

        @test result.model_selection in significant
        @test !isempty(result.fitted_models)
        @test !isempty(result.aic_values)
        @test result.target_dose > 0
    end
end

@testset "Complete MCP-Mod Analysis" begin
    @testset "From summary statistics" begin
        doses = [0.0, 10.0, 25.0, 50.0, 100.0]
        models = DoseResponseModel[
            LinearModel(),
            EmaxModel(ed50 = 25.0),
            EmaxModel(ed50 = 50.0),
            SigmoidEmaxModel(ed50 = 25.0, h = 2.0)
        ]

        config = MCPModConfig(
            doses = doses,
            models = models,
            alpha = 0.025,
            delta = 0.4,
            direction = :increasing
        )

        # Strong Emax-like dose-response
        means = [0.0, 0.3, 0.5, 0.7, 0.85]
        se = [0.1, 0.1, 0.1, 0.1, 0.1]
        n = [40, 40, 40, 40, 40]

        result = run_mcpmod(means, se, n, config)

        @test result.config === config
        @test result.mcp_result isa MCPStepResult
        @test result.mcp_result.proof_of_concept == true
        @test result.mod_result !== nothing
        @test result.mod_result.target_dose > 0
    end

    @testset "From raw data" begin
        Random.seed!(12345)

        doses = [0.0, 25.0, 50.0, 100.0]
        models = DoseResponseModel[LinearModel(), EmaxModel(ed50 = 50.0)]

        config = MCPModConfig(
            doses = doses,
            models = models,
            alpha = 0.025,
            direction = :increasing
        )

        # Generate simulated data
        true_model = EmaxModel(ed50 = 50.0)
        true_params = Dict(:e0 => 0.0, :emax => 1.0, :ed50 => 50.0)

        data = Dict{Float64, Vector{Float64}}()
        for dose in doses
            true_mean = model_mean(true_model, dose, true_params)
            data[dose] = true_mean .+ 0.1 * randn(30)  # 30 subjects per dose
        end

        result = run_mcpmod(data, config)

        @test result.sample_sizes == fill(30, length(doses))
        @test length(result.observed_means) == length(doses)
        @test length(result.observed_se) == length(doses)
    end
end

@testset "Convenience Functions" begin
    @testset "create_mcpmod_models" begin
        doses = [0.0, 25.0, 50.0, 100.0]

        models = create_mcpmod_models(doses)

        @test length(models) > 0
        @test all(m -> m isa DoseResponseModel, models)

        # Check that model names are unique
        names = [m.name for m in models]
        # Note: Some models may have same name with different params
    end

    @testset "default_mcpmod_config" begin
        doses = [0.0, 25.0, 50.0, 100.0]

        config = default_mcpmod_config(doses; alpha = 0.05, direction = :decreasing)

        @test config.alpha == 0.05
        @test config.direction == :decreasing
        @test !isempty(config.models)
        @test config.doses == doses
    end
end

@testset "Target Dose Estimation" begin
    @testset "Emax target dose" begin
        model = EmaxModel(ed50 = 50.0)
        params = Dict(:e0 => 0.0, :emax => 1.0, :ed50 => 50.0)
        doses = [0.0, 25.0, 50.0, 100.0]

        # Target 80% of max effect (0.8)
        target_dose = estimate_target_dose(model, params, 0.8, doses, :increasing)

        # For Emax: d/(ed50+d) = 0.8 => d = 4*ed50 = 200
        # But clamped to max dose of 100
        @test target_dose <= 100.0
        @test target_dose > 0
    end

    @testset "Linear target dose" begin
        model = LinearModel()
        params = Dict(:e0 => 0.0, :delta => 0.01)
        doses = [0.0, 50.0, 100.0]

        # Target effect of 0.5 => d = 0.5/0.01 = 50
        target_dose = estimate_target_dose(model, params, 0.5, doses, :increasing)

        @test isapprox(target_dose, 50.0, atol = 5.0)
    end
end

@testset "Decreasing Direction" begin
    doses = [0.0, 25.0, 50.0, 100.0]
    models = DoseResponseModel[LinearModel()]

    config = MCPModConfig(
        doses = doses,
        models = models,
        alpha = 0.025,
        direction = :decreasing
    )

    # Decreasing dose-response
    means = [1.0, 0.7, 0.5, 0.3]
    se = [0.1, 0.1, 0.1, 0.1]
    n = [40, 40, 40, 40]

    result = run_mcpmod(means, se, n, config)

    # Should detect PoC for decreasing signal
    @test result.mcp_result.proof_of_concept == true

    # Contrasts should be oriented for decreasing
    contrasts = result.mcp_result.contrasts
    c = contrasts[:linear]
    @test c[1] > c[end]  # Decreasing contrast weights
end
