using Test
using NeoPKPD
using LinearAlgebra
using Random

@testset "Mixture Models" begin

    @testset "MixtureComponent" begin
        @testset "Default construction" begin
            comp = MixtureComponent()
            @test comp.name == "Component"
            @test comp.theta === nothing
            @test comp.omega === nothing
            @test isempty(comp.theta_indices)
            @test isempty(comp.omega_indices)
        end

        @testset "Named component with theta" begin
            comp = MixtureComponent(
                name = "Fast Metabolizers",
                theta = [5.0, 50.0],
                theta_indices = [1]
            )
            @test comp.name == "Fast Metabolizers"
            @test comp.theta == [5.0, 50.0]
            @test comp.theta_indices == [1]
        end

        @testset "Component with omega" begin
            omega = diagm([0.09, 0.04])
            comp = MixtureComponent(
                name = "High Variability",
                omega = omega,
                omega_indices = [1, 2]
            )
            @test comp.omega == omega
            @test comp.omega_indices == [1, 2]
        end
    end

    @testset "MixtureSpec" begin
        @testset "Default two-component mixture" begin
            spec = MixtureSpec(n_components = 2)
            @test spec.n_components == 2
            @test length(spec.components) == 2
            @test length(spec.mixing_probabilities) == 2
            @test isapprox(sum(spec.mixing_probabilities), 1.0)
            @test spec.mixing_probabilities == [0.5, 0.5]
        end

        @testset "Three-component mixture" begin
            spec = MixtureSpec(
                n_components = 3,
                mixing_probabilities = [0.5, 0.3, 0.2]
            )
            @test spec.n_components == 3
            @test spec.mixing_probabilities == [0.5, 0.3, 0.2]
        end

        @testset "Custom components" begin
            spec = MixtureSpec(
                n_components = 2,
                components = [
                    MixtureComponent(name="Slow", theta=[2.0, 50.0]),
                    MixtureComponent(name="Fast", theta=[8.0, 50.0])
                ],
                mixing_probabilities = [0.7, 0.3]
            )
            @test spec.components[1].name == "Slow"
            @test spec.components[2].name == "Fast"
            @test spec.components[1].theta == [2.0, 50.0]
            @test spec.components[2].theta == [8.0, 50.0]
        end

        @testset "Parameterization types" begin
            spec_theta = MixtureSpec(n_components=2, parameterization=ThetaOnlyMixture())
            @test spec_theta.parameterization isa ThetaOnlyMixture

            spec_omega = MixtureSpec(n_components=2, parameterization=OmegaOnlyMixture())
            @test spec_omega.parameterization isa OmegaOnlyMixture

            spec_full = MixtureSpec(n_components=2, parameterization=FullMixture())
            @test spec_full.parameterization isa FullMixture
        end

        @testset "Probability bounds" begin
            spec = MixtureSpec(
                n_components = 2,
                probability_bounds = (0.05, 0.95)
            )
            @test spec.probability_bounds == (0.05, 0.95)
        end

        @testset "Invalid specifications" begin
            # Need at least 2 components
            @test_throws AssertionError MixtureSpec(n_components = 1)

            # Probabilities must sum to 1
            @test_throws AssertionError MixtureSpec(
                n_components = 2,
                mixing_probabilities = [0.6, 0.6]
            )

            # Probabilities must be in [0, 1]
            @test_throws AssertionError MixtureSpec(
                n_components = 2,
                mixing_probabilities = [1.5, -0.5]
            )
        end
    end

    @testset "MixtureEstimationMethod" begin
        @testset "MixtureEM defaults" begin
            em = MixtureEM()
            @test em.max_iter == 100
            @test em.tol == 1e-4
            @test em.n_init == 5
            @test em.verbose == false
            @test em.inner_method isa FOCEIMethod
        end

        @testset "MixtureEM custom" begin
            em = MixtureEM(
                max_iter = 200,
                tol = 1e-6,
                n_init = 10,
                verbose = true
            )
            @test em.max_iter == 200
            @test em.tol == 1e-6
            @test em.n_init == 10
            @test em.verbose == true
        end

        @testset "MixtureSAEM defaults" begin
            saem = MixtureSAEM()
            @test saem.n_burn == 200
            @test saem.n_iter == 300
            @test saem.n_mcmc_steps == 50
        end
    end

    @testset "MixtureConfig" begin
        @testset "Basic configuration" begin
            spec = MixtureSpec(n_components = 2)
            sigma = ResidualErrorSpec(
                ProportionalError(),
                ProportionalErrorParams(0.1),
                :conc,
                UInt64(1)
            )

            config = MixtureConfig(
                spec,
                MixtureEM();
                base_theta = [3.0, 50.0],
                base_omega = diagm([0.09, 0.04]),
                sigma = sigma
            )

            @test config.mixture_spec === spec
            @test config.method isa MixtureEM
            @test config.base_theta == [3.0, 50.0]
            @test config.base_omega == diagm([0.09, 0.04])
            @test config.compute_posteriors == true
            @test config.classification_threshold == 0.5
        end

        @testset "With bounds" begin
            spec = MixtureSpec(n_components = 2)
            sigma = ResidualErrorSpec(
                AdditiveError(),
                AdditiveErrorParams(1.0),
                :conc,
                UInt64(1)
            )

            config = MixtureConfig(
                spec,
                MixtureEM();
                base_theta = [3.0, 50.0],
                base_omega = diagm([0.09, 0.04]),
                sigma = sigma,
                theta_lower = [0.1, 10.0],
                theta_upper = [10.0, 100.0]
            )

            @test config.theta_lower == [0.1, 10.0]
            @test config.theta_upper == [10.0, 100.0]
        end
    end

    @testset "Component Log-Likelihood" begin
        sigma = ResidualErrorSpec(
            AdditiveError(),
            AdditiveErrorParams(1.0),
            :conc,
            UInt64(1)
        )

        @testset "Perfect prediction" begin
            y = [1.0, 2.0, 3.0]
            f = [1.0, 2.0, 3.0]  # Perfect match
            ll = NeoPKPD.component_log_likelihood(y, f, sigma)
            # With sigma=1.0, -2LL = n*log(2π) for perfect fit
            expected_ll = -0.5 * 3 * log(2π)
            @test isapprox(ll, expected_ll, atol=1e-6)
        end

        @testset "With residuals" begin
            y = [1.0, 2.0, 3.0]
            f = [1.5, 2.5, 3.5]  # 0.5 residual each
            ll = NeoPKPD.component_log_likelihood(y, f, sigma)
            # Higher residuals -> lower (more negative) LL
            @test ll < -0.5 * 3 * log(2π)
        end
    end

    @testset "Posterior Probabilities" begin
        @testset "Two-component equal priors" begin
            # Observations clearly belonging to component 1
            observations = [[1.0, 1.1, 0.9]]

            # Component 1 predicts ~1.0, Component 2 predicts ~5.0
            predictions_by_component = [
                [[1.0, 1.0, 1.0]],  # Component 1 - good fit
                [[5.0, 5.0, 5.0]]   # Component 2 - poor fit
            ]

            sigma = ResidualErrorSpec(
                AdditiveError(),
                AdditiveErrorParams(0.5),
                :conc,
                UInt64(1)
            )

            mixing_probs = [0.5, 0.5]

            posteriors = compute_posterior_probabilities(
                observations,
                predictions_by_component,
                [sigma],
                mixing_probs
            )

            @test size(posteriors) == (1, 2)
            @test posteriors[1, 1] > 0.99  # Strong preference for component 1
            @test posteriors[1, 2] < 0.01
            @test isapprox(sum(posteriors[1, :]), 1.0)
        end

        @testset "Unequal priors" begin
            observations = [[3.0, 3.0, 3.0]]

            predictions_by_component = [
                [[3.0, 3.0, 3.0]],  # Component 1 - perfect fit
                [[3.0, 3.0, 3.0]]   # Component 2 - perfect fit (same)
            ]

            sigma = ResidualErrorSpec(
                AdditiveError(),
                AdditiveErrorParams(0.5),
                :conc,
                UInt64(1)
            )

            mixing_probs = [0.8, 0.2]

            posteriors = compute_posterior_probabilities(
                observations,
                predictions_by_component,
                [sigma],
                mixing_probs
            )

            # With equal likelihood, posteriors should match priors
            @test isapprox(posteriors[1, 1], 0.8, atol=1e-6)
            @test isapprox(posteriors[1, 2], 0.2, atol=1e-6)
        end

        @testset "Multiple subjects" begin
            observations = [
                [1.0, 1.0],  # Subject 1 - fits component 1
                [5.0, 5.0]   # Subject 2 - fits component 2
            ]

            predictions_by_component = [
                [[1.0, 1.0], [1.0, 1.0]],  # Component 1 predictions
                [[5.0, 5.0], [5.0, 5.0]]   # Component 2 predictions
            ]

            sigma = ResidualErrorSpec(
                AdditiveError(),
                AdditiveErrorParams(0.5),
                :conc,
                UInt64(1)
            )

            mixing_probs = [0.5, 0.5]

            posteriors = compute_posterior_probabilities(
                observations,
                predictions_by_component,
                [sigma],
                mixing_probs
            )

            @test size(posteriors) == (2, 2)
            @test posteriors[1, 1] > 0.99  # Subject 1 -> Component 1
            @test posteriors[2, 2] > 0.99  # Subject 2 -> Component 2
        end
    end

    @testset "Subject Classification" begin
        @testset "Clear classification" begin
            posteriors = [
                0.95 0.05;
                0.10 0.90;
                0.80 0.20
            ]

            assignments = classify_subjects(posteriors)

            @test assignments == [1, 2, 1]
        end

        @testset "Three components" begin
            posteriors = [
                0.1 0.8 0.1;
                0.7 0.2 0.1;
                0.1 0.1 0.8
            ]

            assignments = classify_subjects(posteriors)

            @test assignments == [2, 1, 3]
        end
    end

    @testset "Classification Entropy" begin
        @testset "Perfect separation" begin
            posteriors = [
                1.0 0.0;
                0.0 1.0
            ]
            entropy = NeoPKPD.classification_entropy(posteriors)
            @test isapprox(entropy, 0.0, atol=1e-10)
        end

        @testset "Maximum uncertainty" begin
            posteriors = [
                0.5 0.5;
                0.5 0.5
            ]
            entropy = NeoPKPD.classification_entropy(posteriors)
            # Maximum entropy for 2 classes is log(2)
            @test isapprox(entropy, log(2), atol=1e-6)
        end

        @testset "Partial separation" begin
            posteriors = [
                0.9 0.1;
                0.1 0.9
            ]
            entropy = NeoPKPD.classification_entropy(posteriors)
            @test entropy > 0.0
            @test entropy < log(2)
        end
    end

    @testset "EM M-Step Mixing Probabilities" begin
        @testset "Equal posteriors" begin
            posteriors = [
                0.5 0.5;
                0.5 0.5;
                0.5 0.5;
                0.5 0.5
            ]
            bounds = (0.01, 0.99)

            new_probs = NeoPKPD.em_m_step_mixing_probs(posteriors, bounds)

            @test isapprox(new_probs[1], 0.5, atol=1e-6)
            @test isapprox(new_probs[2], 0.5, atol=1e-6)
            @test isapprox(sum(new_probs), 1.0)
        end

        @testset "Skewed posteriors" begin
            posteriors = [
                0.9 0.1;
                0.8 0.2;
                0.9 0.1;
                0.8 0.2
            ]
            bounds = (0.01, 0.99)

            new_probs = NeoPKPD.em_m_step_mixing_probs(posteriors, bounds)

            # Average of column 1: (0.9+0.8+0.9+0.8)/4 = 0.85
            @test isapprox(new_probs[1], 0.85, atol=1e-6)
            @test isapprox(new_probs[2], 0.15, atol=1e-6)
        end

        @testset "Bounds enforcement" begin
            # All subjects assigned to component 1
            posteriors = [
                1.0 0.0;
                1.0 0.0
            ]
            bounds = (0.05, 0.95)

            new_probs = NeoPKPD.em_m_step_mixing_probs(posteriors, bounds)

            # Should be clamped and renormalized
            @test new_probs[1] <= bounds[2]
            @test new_probs[2] >= bounds[1]
            @test isapprox(sum(new_probs), 1.0)
        end
    end

    @testset "Convenience Functions" begin
        @testset "create_metabolizer_mixture" begin
            spec = create_metabolizer_mixture(2.0, 8.0, 50.0; slow_fraction=0.65)

            @test spec.n_components == 2
            @test spec.components[1].name == "Slow Metabolizers"
            @test spec.components[2].name == "Fast Metabolizers"
            @test spec.components[1].theta == [2.0, 50.0]
            @test spec.components[2].theta == [8.0, 50.0]
            @test spec.mixing_probabilities == [0.65, 0.35]
            @test spec.parameterization isa ThetaOnlyMixture
        end

        @testset "create_responder_mixture" begin
            spec = create_responder_mixture(responder_fraction=0.7)

            @test spec.n_components == 2
            @test spec.components[1].name == "Responders"
            @test spec.components[2].name == "Non-Responders"
            @test isapprox(spec.mixing_probabilities[1], 0.7, atol=1e-10)
            @test isapprox(spec.mixing_probabilities[2], 0.3, atol=1e-10)
            @test spec.parameterization isa FullMixture
        end
    end

    @testset "MixtureResult Structure" begin
        # Create a mock result to test structure
        spec = MixtureSpec(n_components = 2)
        sigma = ResidualErrorSpec(
            AdditiveError(),
            AdditiveErrorParams(1.0),
            :conc,
            UInt64(1)
        )

        config = MixtureConfig(
            spec,
            MixtureEM();
            base_theta = [3.0, 50.0],
            base_omega = diagm([0.09, 0.04]),
            sigma = sigma
        )

        subject_result = MixtureSubjectResult(
            "S001",
            [0.8, 0.2],
            1,
            0.8,
            [-10.0, -20.0],
            [0.0, 0.0],
            [1.0, 2.0, 3.0]
        )

        result = MixtureResult(
            config,
            true,  # converged
            50,    # n_iterations
            -100.0, # log_likelihood
            210.0,  # aic
            220.0,  # bic
            [0.7, 0.3],  # mixing_probabilities
            [[3.0, 50.0], [5.0, 50.0]],  # component_theta
            [diagm([0.09, 0.04]), diagm([0.09, 0.04])],  # component_omega
            sigma,
            [subject_result],
            Dict(1 => 1, 2 => 0),  # classification_summary
            0.5,  # entropy
            String[],  # messages
            1.5  # runtime_seconds
        )

        @test result.converged == true
        @test result.n_iterations == 50
        @test result.log_likelihood == -100.0
        @test result.aic == 210.0
        @test result.bic == 220.0
        @test result.mixing_probabilities == [0.7, 0.3]
        @test length(result.component_theta) == 2
        @test result.classification_summary[1] == 1
        @test result.entropy == 0.5
    end

    @testset "MixtureSubjectResult Structure" begin
        result = MixtureSubjectResult(
            "SUBJ001",
            [0.9, 0.1],
            1,
            0.9,
            [-15.0, -30.0],
            [0.1, -0.05],
            [5.0, 4.0, 3.0, 2.0]
        )

        @test result.subject_id == "SUBJ001"
        @test result.posterior_probabilities == [0.9, 0.1]
        @test result.most_likely_component == 1
        @test result.classification_confidence == 0.9
        @test result.component_likelihoods == [-15.0, -30.0]
        @test result.eta == [0.1, -0.05]
        @test result.ipred == [5.0, 4.0, 3.0, 2.0]
    end

    @testset "Mixture Log-Likelihood" begin
        @testset "Single component dominance" begin
            observations = [[1.0, 1.0, 1.0]]

            predictions_by_component = [
                [[1.0, 1.0, 1.0]],  # Perfect fit
                [[10.0, 10.0, 10.0]]  # Poor fit
            ]

            sigma = ResidualErrorSpec(
                AdditiveError(),
                AdditiveErrorParams(1.0),
                :conc,
                UInt64(1)
            )

            mixing_probs = [0.5, 0.5]

            ll = NeoPKPD.mixture_log_likelihood(
                observations,
                predictions_by_component,
                [sigma],
                mixing_probs
            )

            # Should be close to log-likelihood of component 1 (good fit)
            @test isfinite(ll)
            @test ll < 0  # Log-likelihood is negative
        end

        @testset "Multiple subjects" begin
            observations = [
                [1.0, 1.0],
                [5.0, 5.0]
            ]

            predictions_by_component = [
                [[1.0, 1.0], [1.0, 1.0]],
                [[5.0, 5.0], [5.0, 5.0]]
            ]

            sigma = ResidualErrorSpec(
                AdditiveError(),
                AdditiveErrorParams(1.0),
                :conc,
                UInt64(1)
            )

            mixing_probs = [0.5, 0.5]

            ll = NeoPKPD.mixture_log_likelihood(
                observations,
                predictions_by_component,
                [sigma],
                mixing_probs
            )

            @test isfinite(ll)
            @test ll < 0
        end
    end

    @testset "Pretty Printing" begin
        @testset "MixtureSpec show" begin
            spec = create_metabolizer_mixture(2.0, 8.0, 50.0)
            io = IOBuffer()
            show(io, spec)
            output = String(take!(io))
            @test contains(output, "MixtureSpec")
            @test contains(output, "Components: 2")
            @test contains(output, "Slow Metabolizers")
            @test contains(output, "Fast Metabolizers")
        end
    end

end
