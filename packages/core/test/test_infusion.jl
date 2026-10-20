# Test suite for IV Infusion functionality
# Tests zero-order infusion (constant rate) administration

using Test
using NeoPKPD

@testset "IV Infusion Tests" begin
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

    @testset "DoseEvent with duration" begin
        # Bolus dose (default)
        d_bolus = DoseEvent(0.0, 100.0)
        @test d_bolus.duration == 0.0
        @test is_bolus(d_bolus)
        @test !is_infusion(d_bolus)
        @test infusion_rate(d_bolus) == Inf

        # Infusion dose
        d_infusion = DoseEvent(0.0, 100.0, 1.0)
        @test d_infusion.duration == 1.0
        @test !is_bolus(d_infusion)
        @test is_infusion(d_infusion)
        @test infusion_rate(d_infusion) ≈ 100.0

        # Dose end time
        @test dose_end_time(d_bolus) == 0.0
        @test dose_end_time(d_infusion) == 1.0
    end

    @testset "OneCompIVBolus with infusion" begin
        params = OneCompIVBolusParams(10.0, 50.0)  # CL=10, V=50

        # Bolus simulation
        spec_bolus = ModelSpec(OneCompIVBolus(), "bolus", params, [DoseEvent(0.0, 100.0)])
        result_bolus = simulate(spec_bolus, grid, solver)

        # Infusion simulation (1 hour)
        spec_infusion = ModelSpec(OneCompIVBolus(), "infusion", params, [DoseEvent(0.0, 100.0, 1.0)])
        result_infusion = simulate(spec_infusion, grid, solver)

        # Infusion should have lower Cmax than bolus (same dose)
        @test maximum(result_infusion.observations[:conc]) < maximum(result_bolus.observations[:conc])

        # At t=0, infusion should start at 0 concentration
        @test result_infusion.observations[:conc][1] ≈ 0.0 atol=1e-10

        # Bolus Cmax should be at t=0
        @test result_bolus.observations[:conc][1] ≈ maximum(result_bolus.observations[:conc])
    end

    @testset "TwoCompIVBolus with infusion" begin
        params = TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0)  # CL, V1, Q, V2

        # 2-hour infusion
        spec = ModelSpec(TwoCompIVBolus(), "2comp_infusion", params, [DoseEvent(0.0, 500.0, 2.0)])
        result = simulate(spec, grid, solver)

        @test length(result.observations[:conc]) == length(grid.saveat)
        @test result.observations[:conc][1] ≈ 0.0 atol=1e-10
        @test maximum(result.observations[:conc]) > 0
    end

    @testset "ThreeCompIVBolus with infusion" begin
        params = ThreeCompIVBolusParams(15.0, 30.0, 8.0, 60.0, 3.0, 120.0)

        # 30-minute infusion
        spec = ModelSpec(ThreeCompIVBolus(), "3comp_infusion", params, [DoseEvent(0.0, 200.0, 0.5)])
        result = simulate(spec, grid, solver)

        @test length(result.observations[:conc]) == length(grid.saveat)
        @test maximum(result.observations[:conc]) > 0
    end

    @testset "MichaelisMentenElimination with infusion" begin
        params = MichaelisMentenEliminationParams(100.0, 5.0, 50.0)  # Vmax, Km, V

        # 1-hour infusion
        spec = ModelSpec(MichaelisMentenElimination(), "mm_infusion", params, [DoseEvent(0.0, 200.0, 1.0)])
        result = simulate(spec, grid, solver)

        @test length(result.observations[:conc]) == length(grid.saveat)
        @test result.observations[:conc][1] ≈ 0.0 atol=1e-10
    end

    @testset "InfusionSchedule" begin
        doses = [
            DoseEvent(0.0, 100.0),           # bolus at t=0
            DoseEvent(1.0, 200.0, 2.0),      # infusion from t=1 to t=3
            DoseEvent(6.0, 150.0, 1.0),      # infusion from t=6 to t=7
        ]

        schedule = build_infusion_schedule(doses, 0.0, 24.0)

        @test schedule.initial_amount == 100.0
        @test length(schedule.start_times) == 2
        @test length(schedule.bolus_times) == 0  # t=0 bolus goes to initial_amount

        # Check infusion rates at different times
        @test compute_infusion_rate_at_time(schedule, 0.5) == 0.0
        @test compute_infusion_rate_at_time(schedule, 1.5) ≈ 100.0  # 200/2
        @test compute_infusion_rate_at_time(schedule, 6.5) ≈ 150.0  # 150/1
        @test compute_infusion_rate_at_time(schedule, 10.0) == 0.0
    end

    @testset "Serialization with duration" begin
        params = OneCompIVBolusParams(10.0, 50.0)
        spec = ModelSpec(OneCompIVBolus(), "test", params, [DoseEvent(0.0, 100.0, 2.0)])
        result = simulate(spec, grid, solver)

        # Serialize
        artifact = serialize_execution(
            model_spec=spec,
            grid=grid,
            solver=solver,
            result=result
        )

        # Check duration is in serialized output
        @test artifact["model_spec"]["doses"][1]["duration"] == 2.0

        # Deserialize and replay
        parsed = deserialize_execution(artifact)
        @test parsed.model_spec.doses[1].duration == 2.0

        replayed = simulate(parsed.model_spec, parsed.grid, parsed.solver)
        @test maximum(abs.(result.observations[:conc] .- replayed.observations[:conc])) < 1e-10
    end

    @testset "Backward compatibility (no duration in JSON)" begin
        # Simulate old artifact format without duration field
        artifact = Dict(
            "artifact_schema_version" => "1.0.0",
            "model_spec" => Dict(
                "kind" => "NeoPKPD.OneCompIVBolus",
                "name" => "old_format",
                "params" => Dict("CL" => 10.0, "V" => 50.0),
                "doses" => [Dict("time" => 0.0, "amount" => 100.0)]  # No duration field
            ),
            "grid" => Dict("t0" => 0.0, "t1" => 24.0, "saveat" => collect(0.0:1.0:24.0)),
            "solver" => Dict("alg" => "Tsit5", "reltol" => 1e-6, "abstol" => 1e-8, "maxiters" => 10000)
        )

        parsed = deserialize_execution(artifact)
        @test parsed.model_spec.doses[1].duration == 0.0  # Should default to bolus
        @test is_bolus(parsed.model_spec.doses[1])
    end
end
