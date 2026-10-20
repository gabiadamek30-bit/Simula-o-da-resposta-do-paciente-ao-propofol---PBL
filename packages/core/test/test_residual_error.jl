# Test suite for Residual Error Models

using Test
using NeoPKPD

@testset "Residual Error Models" begin
    # Test predictions
    F = [1.0, 2.0, 5.0, 10.0, 20.0]

    @testset "AdditiveError" begin
        spec = ResidualErrorSpec(AdditiveError(), AdditiveErrorParams(0.5), :conc, 12345)

        # Test SD is constant
        sds = residual_sd.(F, Ref(spec))
        @test all(sds .== 0.5)

        # Test variance
        vars = residual_variance.(F, Ref(spec))
        @test all(vars .== 0.25)

        # Test apply_residual_error returns correct length
        Y = apply_residual_error(F, spec)
        @test length(Y) == length(F)

        # Test reproducibility with same seed
        spec2 = ResidualErrorSpec(AdditiveError(), AdditiveErrorParams(0.5), :conc, 12345)
        Y2 = apply_residual_error(F, spec2)
        @test Y ≈ Y2

        # Test log-likelihood is finite
        ll = residual_log_likelihood(Y, F, spec)
        @test isfinite(ll)

        # Test IWRES
        iwres = compute_iwres(Y, F, spec)
        @test length(iwres) == length(F)
    end

    @testset "ProportionalError" begin
        spec = ResidualErrorSpec(ProportionalError(), ProportionalErrorParams(0.1), :conc, 12345)

        # Test SD is proportional to F
        sds = residual_sd.(F, Ref(spec))
        @test sds ≈ 0.1 .* F

        # Test variance
        vars = residual_variance.(F, Ref(spec))
        @test vars ≈ (0.1 .* F).^2

        # Test apply_residual_error
        Y = apply_residual_error(F, spec)
        @test length(Y) == length(F)

        # Test log-likelihood
        ll = residual_log_likelihood(Y, F, spec)
        @test isfinite(ll)

        # Test IWRES
        iwres = compute_iwres(Y, F, spec)
        @test length(iwres) == length(F)
    end

    @testset "CombinedError" begin
        spec = ResidualErrorSpec(CombinedError(), CombinedErrorParams(0.1, 0.1), :conc, 12345)

        # Test SD is combined
        sds = residual_sd.(F, Ref(spec))
        expected_sds = sqrt.(0.1^2 .+ (0.1 .* F).^2)
        @test sds ≈ expected_sds

        # Test apply_residual_error
        Y = apply_residual_error(F, spec)
        @test length(Y) == length(F)

        # Test log-likelihood
        ll = residual_log_likelihood(Y, F, spec)
        @test isfinite(ll)

        # Test IWRES
        iwres = compute_iwres(Y, F, spec)
        @test length(iwres) == length(F)
    end

    @testset "ExponentialError" begin
        spec = ResidualErrorSpec(ExponentialError(), ExponentialErrorParams(0.2), :conc, 12345)

        # Test apply_residual_error returns positive values
        Y = apply_residual_error(F, spec)
        @test length(Y) == length(F)
        @test all(Y .> 0)

        # Test log-likelihood
        ll = residual_log_likelihood(Y, F, spec)
        @test isfinite(ll)

        # Test IWRES
        iwres = compute_iwres(Y, F, spec)
        @test length(iwres) == length(F)
    end

    @testset "Parameter Validation" begin
        # Additive: sigma must be positive
        @test_throws Exception AdditiveErrorParams(0.0)
        @test_throws Exception AdditiveErrorParams(-0.1)

        # Proportional: sigma must be positive
        @test_throws Exception ProportionalErrorParams(0.0)
        @test_throws Exception ProportionalErrorParams(-0.1)

        # Combined: at least one must be positive
        @test_throws Exception CombinedErrorParams(0.0, 0.0)
        @test_throws Exception CombinedErrorParams(-0.1, 0.1)

        # This should work (one is positive)
        @test_nowarn CombinedErrorParams(0.1, 0.0)
        @test_nowarn CombinedErrorParams(0.0, 0.1)

        # Exponential: sigma must be positive
        @test_throws Exception ExponentialErrorParams(0.0)
        @test_throws Exception ExponentialErrorParams(-0.1)
    end

    @testset "Serialization Round-Trip" begin
        # Test additive
        add_spec = ResidualErrorSpec(AdditiveError(), AdditiveErrorParams(0.5), :conc, 12345)
        add_dict = serialize_error_spec(add_spec)
        add_spec2 = deserialize_error_spec(add_dict)
        @test add_spec2.params.sigma == add_spec.params.sigma
        @test add_spec2.observation == add_spec.observation
        @test add_spec2.seed == add_spec.seed

        # Test proportional
        prop_spec = ResidualErrorSpec(ProportionalError(), ProportionalErrorParams(0.1), :conc, 54321)
        prop_dict = serialize_error_spec(prop_spec)
        prop_spec2 = deserialize_error_spec(prop_dict)
        @test prop_spec2.params.sigma == prop_spec.params.sigma

        # Test combined
        comb_spec = ResidualErrorSpec(CombinedError(), CombinedErrorParams(0.1, 0.2), :effect, 99999)
        comb_dict = serialize_error_spec(comb_spec)
        comb_spec2 = deserialize_error_spec(comb_dict)
        @test comb_spec2.params.sigma_add == comb_spec.params.sigma_add
        @test comb_spec2.params.sigma_prop == comb_spec.params.sigma_prop

        # Test exponential
        exp_spec = ResidualErrorSpec(ExponentialError(), ExponentialErrorParams(0.3), :conc, 11111)
        exp_dict = serialize_error_spec(exp_spec)
        exp_spec2 = deserialize_error_spec(exp_dict)
        @test exp_spec2.params.sigma == exp_spec.params.sigma
    end

    @testset "Integration with PK Simulation" begin
        # Run a simulation and apply error
        spec = ModelSpec(
            OneCompIVBolus(),
            "test",
            OneCompIVBolusParams(10.0, 50.0),
            [DoseEvent(0.0, 100.0)]
        )
        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        result = simulate(spec, grid, solver)
        F = result.observations[:conc]

        # Apply combined error
        error_spec = ResidualErrorSpec(CombinedError(), CombinedErrorParams(0.01, 0.1), :conc, 42)
        Y = apply_residual_error(F, error_spec)

        @test length(Y) == length(F)

        # Compute log-likelihood
        ll = residual_log_likelihood(Y, F, error_spec)
        @test isfinite(ll)

        # Compute IWRES
        iwres = compute_iwres(Y, F, error_spec)
        @test length(iwres) == length(F)
    end
end
