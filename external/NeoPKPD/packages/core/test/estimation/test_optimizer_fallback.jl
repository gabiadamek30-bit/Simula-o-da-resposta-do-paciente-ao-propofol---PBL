using Test
using NeoPKPD
using Optim
using LinearAlgebra

@testset "Optimizer Fallback" begin

    @testset "OptimizerConfig" begin
        @testset "Default configuration" begin
            config = OptimizerConfig()
            @test config.primary == BFGS_OPTIMIZER
            @test LBFGS_OPTIMIZER in config.fallback_chain
            @test NELDER_MEAD_OPTIMIZER in config.fallback_chain
            @test config.max_attempts_per_optimizer >= 1
            @test config.scale_initial_step == 1.0
        end

        @testset "Custom configuration" begin
            config = OptimizerConfig(
                primary=LBFGS_OPTIMIZER,
                fallback_chain=[BFGS_OPTIMIZER],
                max_attempts_per_optimizer=3,
                verbose=true
            )
            @test config.primary == LBFGS_OPTIMIZER
            @test config.fallback_chain == [BFGS_OPTIMIZER]
            @test config.max_attempts_per_optimizer == 3
            @test config.verbose == true
        end

        @testset "No fallback configuration" begin
            config = OptimizerConfig(
                fallback_chain=NeoPKPD.OptimizerType[]
            )
            @test isempty(config.fallback_chain)
        end
    end

    @testset "Unbounded Optimization" begin
        @testset "Simple quadratic (all optimizers should work)" begin
            # f(x) = (x[1]-2)^2 + (x[2]-3)^2
            f(x) = (x[1] - 2.0)^2 + (x[2] - 3.0)^2
            x0 = [0.0, 0.0]

            result = optimize_with_fallback(f, x0)

            @test result.converged
            @test isapprox(result.minimizer[1], 2.0, atol=1e-4)
            @test isapprox(result.minimizer[2], 3.0, atol=1e-4)
            @test result.minimum < 1e-6
            @test result.optimizer_used == BFGS_OPTIMIZER  # Primary should succeed
            @test result.fallback_count == 0
        end

        @testset "Rosenbrock function (harder problem)" begin
            # Rosenbrock: f(x,y) = (a-x)^2 + b(y-x^2)^2 with a=1, b=100
            f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
            x0 = [-1.0, 1.0]

            config = OptimizerConfig()
            options = Optim.Options(iterations=1000, g_tol=1e-8)

            result = optimize_with_fallback(f, x0, config; options=options)

            @test isapprox(result.minimizer[1], 1.0, atol=1e-3)
            @test isapprox(result.minimizer[2], 1.0, atol=1e-3)
            @test result.minimum < 1e-4
        end

        @testset "Different optimizers" begin
            f(x) = sum(x.^2)
            x0 = [5.0, 5.0]

            # Test BFGS
            config_bfgs = OptimizerConfig(primary=BFGS_OPTIMIZER, fallback_chain=OptimizerType[])
            result_bfgs = optimize_with_fallback(f, x0, config_bfgs)
            @test result_bfgs.converged
            @test result_bfgs.minimum < 1e-6

            # Test L-BFGS
            config_lbfgs = OptimizerConfig(primary=LBFGS_OPTIMIZER, fallback_chain=OptimizerType[])
            result_lbfgs = optimize_with_fallback(f, x0, config_lbfgs)
            @test result_lbfgs.converged
            @test result_lbfgs.minimum < 1e-6

            # Test Nelder-Mead
            config_nm = OptimizerConfig(primary=NELDER_MEAD_OPTIMIZER, fallback_chain=OptimizerType[])
            result_nm = optimize_with_fallback(f, x0, config_nm)
            @test result_nm.converged
            @test result_nm.minimum < 1e-6
        end
    end

    @testset "Bounded Optimization" begin
        @testset "Simple bounded quadratic" begin
            f(x) = (x[1] - 5.0)^2 + (x[2] - 5.0)^2
            lower = [0.0, 0.0]
            upper = [3.0, 3.0]  # Optimum at (5,5) but constrained to (3,3)
            x0 = [1.0, 1.0]

            result = optimize_bounded_with_fallback(f, lower, upper, x0)

            @test result.converged
            # Optimum should be at the upper bounds
            @test isapprox(result.minimizer[1], 3.0, atol=1e-4)
            @test isapprox(result.minimizer[2], 3.0, atol=1e-4)
        end

        @testset "Initial point clamping" begin
            f(x) = sum(x.^2)
            lower = [1.0, 1.0]
            upper = [10.0, 10.0]
            x0 = [0.0, 0.0]  # Outside bounds!

            result = optimize_bounded_with_fallback(f, lower, upper, x0)

            # Should still work (clamped to bounds)
            @test all(result.minimizer .>= lower)
            @test all(result.minimizer .<= upper)
        end

        @testset "Interior optimum" begin
            f(x) = (x[1] - 2.0)^2 + (x[2] - 2.0)^2
            lower = [0.0, 0.0]
            upper = [10.0, 10.0]
            x0 = [5.0, 5.0]

            result = optimize_bounded_with_fallback(f, lower, upper, x0)

            @test result.converged
            @test isapprox(result.minimizer[1], 2.0, atol=1e-4)
            @test isapprox(result.minimizer[2], 2.0, atol=1e-4)
        end
    end

    @testset "Fallback Behavior" begin
        @testset "Fallback triggers on bad objective" begin
            # Create an objective that fails for BFGS but works for Nelder-Mead
            # (e.g., non-differentiable at optimum)
            function tricky_objective(x)
                val = sum(abs.(x .- 1.0))  # Non-differentiable at x=1
                return val
            end

            x0 = [0.0, 0.0]

            # With full fallback chain
            config = OptimizerConfig(verbose=false)
            result = optimize_with_fallback(tricky_objective, x0, config)

            # Should still find a reasonable solution
            @test result.minimum < 0.1
        end

        @testset "Multiple attempts" begin
            f(x) = sum(x.^2)
            x0 = [10.0, 10.0]

            config = OptimizerConfig(max_attempts_per_optimizer=3)
            result = optimize_with_fallback(f, x0, config)

            @test result.converged
            @test result.minimum < 1e-6
        end
    end

    @testset "Simple Wrappers" begin
        @testset "robust_optimize" begin
            f(x) = sum(x.^2)
            x0 = [5.0, 5.0]

            result = robust_optimize(f, x0; max_iter=100, g_tol=1e-6)

            @test result.converged
            @test result.minimum < 1e-6
        end

        @testset "robust_optimize_bounded" begin
            f(x) = sum(x.^2)
            lower = [1.0, 1.0]
            upper = [10.0, 10.0]
            x0 = [5.0, 5.0]

            result = robust_optimize_bounded(f, lower, upper, x0)

            # Minimum at bounds since optimum is at 0
            @test isapprox(result.minimizer[1], 1.0, atol=1e-4)
            @test isapprox(result.minimizer[2], 1.0, atol=1e-4)
        end
    end

    @testset "OptimizationResult" begin
        @testset "Result structure" begin
            f(x) = sum(x.^2)
            x0 = [5.0, 5.0]

            result = optimize_with_fallback(f, x0)

            @test length(result.minimizer) == 2
            @test result.minimum isa Float64
            @test result.converged isa Bool
            @test result.iterations isa Int
            @test result.optimizer_used isa OptimizerType
            @test result.fallback_count isa Int
            @test result.messages isa Vector{String}
        end
    end

    @testset "Edge Cases" begin
        @testset "Single dimension" begin
            f(x) = (x[1] - 3.0)^2
            x0 = [0.0]

            result = optimize_with_fallback(f, x0)

            @test result.converged
            @test isapprox(result.minimizer[1], 3.0, atol=1e-4)
        end

        @testset "High dimension" begin
            n = 20
            f(x) = sum((x .- collect(1.0:n)).^2)
            x0 = zeros(n)

            result = optimize_with_fallback(f, x0)

            @test result.converged
            for i in 1:n
                @test isapprox(result.minimizer[i], Float64(i), atol=1e-3)
            end
        end

        @testset "Flat objective" begin
            f(x) = 0.0
            x0 = [1.0, 2.0]

            result = optimize_with_fallback(f, x0)

            @test result.minimum == 0.0
            @test isfinite(result.minimizer[1])
            @test isfinite(result.minimizer[2])
        end
    end

end
