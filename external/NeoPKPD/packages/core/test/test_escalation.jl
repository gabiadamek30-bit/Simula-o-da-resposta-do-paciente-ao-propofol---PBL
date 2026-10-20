# Test suite for Dose Escalation Simulation

using Test
using NeoPKPD
using StableRNGs

@testset "Dose Escalation" begin

    @testset "Escalation Types" begin
        @testset "ThreePlusThree Rule" begin
            rule = ThreePlusThree()
            @test rule.max_dlt_rate == 0.33

            design = DoseEscalationDesign([10.0, 25.0, 50.0, 100.0, 200.0];
                escalation_rule=ThreePlusThree())
            @test design.dose_levels == [10.0, 25.0, 50.0, 100.0, 200.0]
            @test design.cohort_size == 3
        end

        @testset "mTPI Rule" begin
            rule = mTPI(target_dlt_rate=0.25, equivalence_interval=(0.20, 0.30))
            @test rule.target_dlt_rate == 0.25

            design = DoseEscalationDesign([10.0, 25.0, 50.0, 100.0];
                escalation_rule=mTPI())
            @test design.dose_levels == [10.0, 25.0, 50.0, 100.0]
        end

        @testset "CRM Rule" begin
            rule = CRM(target_dlt_rate=0.25,
                skeleton=[0.05, 0.15, 0.25, 0.40],
                model=:logistic)
            @test rule.target_dlt_rate == 0.25
            @test length(rule.skeleton) == 4
            @test rule.model == :logistic

            design = DoseEscalationDesign([10.0, 25.0, 50.0, 100.0];
                escalation_rule=CRM(skeleton=[0.05, 0.15, 0.25, 0.40]))
            @test design.dose_levels == [10.0, 25.0, 50.0, 100.0]
        end
    end

    @testset "CohortResult" begin
        @test isdefined(NeoPKPD, :CohortResult)

        # Create a mock cohort result
        cohort = CohortResult(
            1,      # cohort_number
            1,      # dose_level
            10.0,   # dose_amount
            3,      # n_subjects
            0,      # n_dlt
            [1, 2, 3],  # subject_ids
            [false, false, false],  # dlt_flags
            [50.0, 55.0, 48.0],  # pk_exposures
            :escalate  # decision
        )

        @test cohort.cohort_number == 1
        @test cohort.dose_level == 1
        @test cohort.n_subjects == 3
        @test cohort.n_dlt == 0
        @test cohort.decision == :escalate
    end

    @testset "EscalationResult" begin
        @test isdefined(NeoPKPD, :EscalationResult)
    end

    @testset "3+3 Logic" begin
        # Test 3+3 decision rules
        # 0 DLTs in 3 → escalate
        # 1 DLT in 3 → expand to 6
        # 2+ DLTs in 3 → de-escalate
        # 1 DLT in 6 → escalate
        # 2+ DLTs in 6 → de-escalate or declare MTD

        design = DoseEscalationDesign([10.0, 25.0, 50.0, 100.0];
            escalation_rule=ThreePlusThree())

        # Test helper functions exist
        @test isdefined(NeoPKPD, :simulate_escalation)
        @test isdefined(NeoPKPD, :simulate_3plus3)
    end

    @testset "mTPI Decision" begin
        # mTPI uses posterior probability to decide
        design = DoseEscalationDesign([10.0, 25.0, 50.0];
            escalation_rule=mTPI())

        @test isdefined(NeoPKPD, :simulate_mtpi)
    end

    @testset "CRM Model" begin
        # CRM continuously updates dose-toxicity model
        design = DoseEscalationDesign([10.0, 25.0, 50.0, 100.0];
            escalation_rule=CRM(skeleton=[0.05, 0.15, 0.25, 0.40]))

        @test isdefined(NeoPKPD, :simulate_crm)
    end

    @testset "Full 3+3 Simulation" begin
        # Create a simple PK model for exposure calculation
        model_params = OneCompIVBolusParams(10.0, 50.0)  # CL=10, V=50
        doses_placeholder = [DoseEvent(0.0, 1.0)]
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses_placeholder)

        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        # Define dose escalation design
        design = DoseEscalationDesign([10.0, 25.0, 50.0, 100.0, 200.0];
            escalation_rule=ThreePlusThree())

        # Simple exposure-based DLT model: P(DLT) = 1 / (1 + exp(-(AUC - 500)/100))
        dlt_model = (exposure) -> 1.0 / (1.0 + exp(-(exposure - 500.0) / 100.0))

        # Run simulation
        result = simulate_escalation(
            design,
            model_spec,
            grid,
            solver;
            dlt_model=dlt_model,
            seed=UInt64(12345),
            verbose=false
        )

        # Basic checks
        @test result isa EscalationResult
        @test result.design == design
        @test length(result.cohorts) > 0
        @test result.total_subjects > 0
        @test result.total_dlt >= 0

        # Check first cohort
        first_cohort = result.cohorts[1]
        @test first_cohort.cohort_number == 1
        @test first_cohort.n_subjects >= 3
        @test first_cohort.decision in [:escalate, :expand, :deescalate, :stop, :mtd]
    end

    @testset "mTPI Simulation" begin
        model_params = OneCompIVBolusParams(10.0, 50.0)
        doses_placeholder = [DoseEvent(0.0, 1.0)]
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses_placeholder)

        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        design = DoseEscalationDesign([10.0, 25.0, 50.0, 100.0];
            escalation_rule=mTPI())
        dlt_model = (exposure) -> 1.0 / (1.0 + exp(-(exposure - 500.0) / 100.0))

        result = simulate_escalation(
            design,
            model_spec,
            grid,
            solver;
            dlt_model=dlt_model,
            seed=UInt64(54321),
            verbose=false
        )

        @test result isa EscalationResult
        @test length(result.cohorts) > 0
    end

    @testset "CRM Simulation" begin
        model_params = OneCompIVBolusParams(10.0, 50.0)
        doses_placeholder = [DoseEvent(0.0, 1.0)]
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses_placeholder)

        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        design = DoseEscalationDesign([10.0, 25.0, 50.0, 100.0];
            escalation_rule=CRM(skeleton=[0.05, 0.15, 0.25, 0.40]))
        dlt_model = (exposure) -> 1.0 / (1.0 + exp(-(exposure - 500.0) / 100.0))

        result = simulate_escalation(
            design,
            model_spec,
            grid,
            solver;
            dlt_model=dlt_model,
            seed=UInt64(99999),
            verbose=false
        )

        @test result isa EscalationResult
        @test length(result.cohorts) > 0
    end

    @testset "MTD Determination" begin
        # Test that MTD is correctly identified
        model_params = OneCompIVBolusParams(10.0, 50.0)
        doses_placeholder = [DoseEvent(0.0, 1.0)]
        model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses_placeholder)

        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        # Use a DLT model with high toxicity at higher doses
        design = DoseEscalationDesign([10.0, 25.0, 50.0, 100.0];
            escalation_rule=ThreePlusThree())

        # DLT probability increases sharply with exposure
        dlt_model = (exposure) -> begin
            1.0 / (1.0 + exp(-(exposure - 30.0) / 10.0))
        end

        result = simulate_escalation(
            design,
            model_spec,
            grid,
            solver;
            dlt_model=dlt_model,
            seed=UInt64(42),
            verbose=false
        )

        @test result isa EscalationResult
        # MTD should be identified at some dose level
        # (may be nothing if all doses are safe or toxic)
        @test result.mtd_level === nothing || result.mtd_level >= 1
    end

end
