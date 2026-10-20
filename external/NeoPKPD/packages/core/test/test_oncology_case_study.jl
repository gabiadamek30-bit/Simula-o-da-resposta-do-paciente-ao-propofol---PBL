# Oncology Case Study Integration Tests
# Validates that all trial simulation modules work together correctly

using Test
using NeoPKPD
using Statistics
using Random

@testset "Oncology Case Study Integration" begin

    # =========================================================================
    # PHASE I: Dose Escalation Integration
    # =========================================================================
    @testset "Phase I Dose Escalation" begin
        # Define common parameters
        dose_levels = [7.0, 21.0, 70.0, 140.0, 280.0, 420.0]
        starting_dose = 7.0

        # PK model for exposure simulation
        pk_params = TwoCompIVBolusParams(0.2, 3.0, 0.5, 2.0)  # CL, V1, Q, V2
        pk_model = ModelSpec(
            TwoCompIVBolus(),
            "test_mab",
            pk_params,
            [DoseEvent(0.0, starting_dose)]
        )

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

        # DLT model with known properties
        dlt_model = (:logistic, -3.0, 0.01)

        @testset "3+3 Design" begin
            design = DoseEscalationDesign(
                dose_levels;
                starting_dose = starting_dose,
                max_subjects = 24,
                cohort_size = 3,
                escalation_rule = ThreePlusThree()
            )

            result = simulate_escalation(
                design, pk_model, grid, solver;
                dlt_model = dlt_model,
                seed = UInt64(12345)
            )

            @test result.total_subjects > 0
            @test result.total_subjects <= 24
            @test length(result.cohorts) > 0
            @test result.completed
        end

        @testset "BOIN Design" begin
            boin_rule = BOIN(target_dlt_rate = 0.25)

            @test boin_rule.target_dlt_rate == 0.25
            @test 0 < boin_rule.cutoff_e < boin_rule.cutoff_d

            design = DoseEscalationDesign(
                dose_levels;
                starting_dose = starting_dose,
                max_subjects = 24,
                cohort_size = 3,
                escalation_rule = boin_rule
            )

            result = simulate_escalation(
                design, pk_model, grid, solver;
                dlt_model = dlt_model,
                seed = UInt64(12346)
            )

            @test result.total_subjects > 0
            @test result.completed
        end

        @testset "CRM Design" begin
            skeleton = [0.05, 0.10, 0.20, 0.30, 0.45, 0.60]
            crm_rule = CRM(
                target_dlt_rate = 0.25,
                skeleton = skeleton
            )

            design = DoseEscalationDesign(
                dose_levels;
                starting_dose = starting_dose,
                max_subjects = 24,
                cohort_size = 3,
                escalation_rule = crm_rule
            )

            result = simulate_escalation(
                design, pk_model, grid, solver;
                dlt_model = dlt_model,
                seed = UInt64(12347)
            )

            @test result.total_subjects > 0
            @test result.completed
        end
    end

    # =========================================================================
    # PHASE II: MCP-Mod Integration
    # =========================================================================
    @testset "Phase II MCP-Mod" begin
        # Simulated dose-response data
        doses = [0.0, 25.0, 50.0, 100.0, 150.0]
        n_per_dose = 30

        # Simulate data with true Emax response
        Random.seed!(42)
        true_e0 = -5.0
        true_emax = -30.0
        true_ed50 = 40.0

        observed_means = Float64[]
        observed_sds = Float64[]
        group_n = Int[]

        for dose in doses
            true_mean = true_e0 + true_emax * dose / (true_ed50 + dose)
            responses = true_mean .+ randn(n_per_dose) * 12.0

            push!(observed_means, mean(responses))
            push!(observed_sds, std(responses) / sqrt(n_per_dose))  # SE, not SD
            push!(group_n, n_per_dose)
        end

        @testset "MCP-Mod Configuration" begin
            models = DoseResponseModel[
                LinearModel(),
                EmaxModel(ed50 = 30.0),
                EmaxModel(ed50 = 60.0)
            ]

            config = MCPModConfig(
                doses = doses,
                models = models,
                alpha = 0.025,
                direction = :decreasing  # More negative is better
            )

            @test config.doses == doses
            @test length(config.models) == 3
        end

        @testset "MCP-Mod Analysis" begin
            models = DoseResponseModel[
                LinearModel(),
                EmaxModel(ed50 = 40.0),
                EmaxModel(ed50 = 80.0)
            ]

            config = MCPModConfig(
                doses = doses,
                models = models,
                alpha = 0.05,
                direction = :decreasing
            )

            result = run_mcpmod(observed_means, observed_sds, group_n, config)

            @test result isa MCPModResult
            @test result.mcp_result isa MCPStepResult
            @test result.mcp_result.proof_of_concept isa Bool
            @test !isempty(result.mcp_result.p_values)
        end
    end

    # =========================================================================
    # EXPOSURE-RESPONSE Integration
    # =========================================================================
    @testset "Exposure-Response Analysis" begin
        Random.seed!(123)

        # Simulate E-R data
        n = 150
        exposures = rand(n) .* 1000  # AUC: 0-1000
        true_e0 = -5.0
        true_emax = -35.0
        true_ec50 = 400.0

        true_responses = true_e0 .+ true_emax .* exposures ./ (true_ec50 .+ exposures)
        responses = true_responses .+ randn(n) .* 10.0

        @testset "E-R Configuration" begin
            config = ExposureResponseConfig(
                exposure_metric = AUCMetric(),
                response_type = ContinuousResponse(),
                model = EmaxERModel()
            )

            @test config.exposure_metric isa AUCMetric
            @test config.model isa EmaxERModel
        end

        @testset "E-R Model Fitting" begin
            config = ExposureResponseConfig(
                exposure_metric = AUCMetric(),
                response_type = ContinuousResponse(),
                model = EmaxERModel(),
                n_bootstrap = 50
            )

            result = run_exposure_response(exposures, responses, config)

            @test haskey(result.parameter_estimates, :e0)
            @test haskey(result.parameter_estimates, :emax)
            @test haskey(result.parameter_estimates, :ec50)

            # Check parameter estimates are reasonable
            @test -20.0 < result.parameter_estimates[:e0] < 10.0
            @test result.parameter_estimates[:emax] < 0  # Negative effect
            @test result.parameter_estimates[:ec50] > 0
        end
    end

    # =========================================================================
    # PROBABILITY OF SUCCESS Integration
    # =========================================================================
    @testset "Probability of Success" begin
        @testset "Prior Specification" begin
            prior = NormalPoSPrior(-20.0, 25.0; source = :previous_trial)

            @test prior.mean == -20.0
            @test prior.variance == 25.0
            @test prior.source == :previous_trial
        end

        @testset "Assurance Calculation" begin
            prior = NormalPoSPrior(-15.0, 16.0)

            config = PoSConfig(
                prior = prior,
                criterion = SuperioritySuccess(margin = 0.0, alpha = 0.025),
                n_subjects_planned = 400,
                n_samples = 2000,
                seed = 54321
            )

            assurance = compute_assurance(config)

            @test 0 <= assurance.assurance <= 1
            # CI bounds might be equal if assurance is very high or low
            @test assurance.assurance_ci[1] <= assurance.assurance_ci[2]
            @test isfinite(assurance.predictive_mean)
        end

        @testset "Full PoS Analysis" begin
            prior = NormalPoSPrior(-20.0, 20.0; source = :meta_analysis)

            config = PoSConfig(
                prior = prior,
                criterion = SuperioritySuccess(),
                n_subjects_planned = 300,
                n_samples = 1000,
                seed = 67890
            )

            result = compute_pos(config; run_sensitivity = true)

            @test result.assurance isa AssuranceResult
            @test result.go_nogo_recommendation in [:go, :nogo, :uncertain]
            @test !isempty(result.decision_rationale)
            @test !isempty(result.sensitivity_results)
        end

        @testset "Decision Boundary" begin
            prior = NormalPoSPrior(-25.0, 16.0)

            config = PoSConfig(
                prior = prior,
                criterion = SuperioritySuccess(),
                n_subjects_planned = 200,
                n_samples = 500,
                seed = 11111
            )

            min_n = pos_decision_boundary(config; target_assurance = 0.75)

            @test min_n > 0
            @test min_n <= 2000  # Could be at max if no solution found
        end
    end

    # =========================================================================
    # FULL PIPELINE Integration
    # =========================================================================
    @testset "Full Pipeline Workflow" begin
        # This tests the complete flow from Phase I to Phase III decision

        # Step 1: Phase I determines RP2D
        pipeline_dose_levels = [10.0, 30.0, 100.0, 200.0, 400.0]
        pk_params = OneCompIVBolusParams(5.0, 50.0)  # CL, V
        pk_model = ModelSpec(
            OneCompIVBolus(),
            "pipeline_test",
            pk_params,
            [DoseEvent(0.0, 10.0)]
        )

        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

        design = DoseEscalationDesign(
            pipeline_dose_levels;
            starting_dose = 10.0,
            max_subjects = 30,
            cohort_size = 3,
            escalation_rule = ThreePlusThree()
        )

        phase1_result = simulate_escalation(
            design, pk_model, grid, solver;
            dlt_model = (:logistic, -2.5, 0.01),
            seed = UInt64(99999)
        )

        @test phase1_result.completed

        # Step 2: Phase II with MCP-Mod
        doses = [0.0, 30.0, 100.0, 200.0]
        Random.seed!(88888)

        means = [-5.0, -12.0, -20.0, -25.0]
        ses = [1.8, 2.0, 2.2, 2.0]  # Standard errors
        ns = [40, 40, 40, 40]

        models = DoseResponseModel[
            LinearModel(),
            EmaxModel(ed50 = 50.0)
        ]

        mcpmod_config = MCPModConfig(
            doses = doses,
            models = models,
            alpha = 0.05,
            direction = :decreasing
        )

        phase2_result = run_mcpmod(means, ses, ns, mcpmod_config)

        @test phase2_result isa MCPModResult
        @test phase2_result.mcp_result.proof_of_concept isa Bool

        # Step 3: PoS for Phase III
        prior = NormalPoSPrior(-18.0, 16.0; source = :previous_trial)

        pos_config = PoSConfig(
            prior = prior,
            criterion = SuperioritySuccess(margin = 0.0, alpha = 0.025),
            n_subjects_planned = 400,
            n_samples = 1000,
            seed = 77777
        )

        pos_result = compute_pos(pos_config)

        @test pos_result.go_nogo_recommendation in [:go, :nogo, :uncertain]

        # Verify pipeline completed successfully
        @test phase1_result.completed && pos_result.go_nogo_recommendation !== nothing
    end
end
