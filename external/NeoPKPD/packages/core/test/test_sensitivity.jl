# Sensitivity analysis tests

@testset "Single-run sensitivity: increasing CL reduces concentration" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "sens_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    out = run_sensitivity(spec, grid, solver; plan=plan, observation=:conc)

    # Concentration should be lower at positive times after increasing CL
    # We avoid t=0 because it reflects initial condition and bolus application.
    @test out.pert_metric_series[2] < out.base_metric_series[2]
    @test out.metrics.max_abs_delta > 0.0
end

@testset "Population sensitivity: increasing CL reduces mean concentration" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "sens_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(9999), 20)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    out = run_population_sensitivity(pop, grid, solver; plan=plan, observation=:conc)

    @test out.pert_summary_mean[2] < out.base_summary_mean[2]
    @test out.metrics_mean.max_abs_delta > 0.0
end

@testset "Sensitivity artifact replay matches stored series (single)" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "sens_art_single",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    r1 = run_sensitivity(spec, grid, solver; plan=plan, observation=:conc)

    art = serialize_sensitivity_execution(
        model_spec=spec, grid=grid, solver=solver, result=r1
    )

    r2 = replay_sensitivity_execution(art)

    for i in eachindex(r1.base_metric_series)
        @test isapprox(
            r2.base_metric_series[i], r1.base_metric_series[i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            r2.pert_metric_series[i], r1.pert_metric_series[i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "Sensitivity artifact replay matches stored series (population)" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "sens_art_pop",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(9999), 20)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    r1 = run_population_sensitivity(
        pop, grid, solver; plan=plan, observation=:conc, probs=[0.05, 0.95]
    )

    art = serialize_population_sensitivity_execution(
        population_spec=pop, grid=grid, solver=solver, result=r1
    )

    r2 = replay_population_sensitivity_execution(art)

    for i in eachindex(r1.base_summary_mean)
        @test isapprox(
            r2.base_summary_mean[i], r1.base_summary_mean[i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            r2.pert_summary_mean[i], r1.pert_summary_mean[i]; rtol=1e-12, atol=1e-12
        )
    end
end
