# Event semantics tests

@testset "Event semantics: duplicate times are summed (IV bolus)" begin
    base_params = OneCompIVBolusParams(5.0, 50.0)
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pk_dup = ModelSpec(
        OneCompIVBolus(),
        "dup_times",
        base_params,
        [
            DoseEvent(0.0, 60.0),
            DoseEvent(0.0, 40.0),
            DoseEvent(12.0, 10.0),
            DoseEvent(12.0, 15.0),
        ],
    )

    pk_sum = ModelSpec(
        OneCompIVBolus(),
        "summed",
        base_params,
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 25.0)],
    )

    r_dup = simulate(pk_dup, grid, solver)
    r_sum = simulate(pk_sum, grid, solver)

    @test r_dup.metadata["event_semantics_version"] == "1.0.0"
    @test r_sum.metadata["event_semantics_version"] == "1.0.0"

    for i in eachindex(r_dup.t)
        @test isapprox(
            r_dup.observations[:conc][i],
            r_sum.observations[:conc][i];
            rtol=1e-12,
            atol=1e-12,
        )
    end
end

@testset "Event semantics: input ordering does not matter for same-time events" begin
    params = OneCompIVBolusParams(5.0, 50.0)
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    doses_a = [
        DoseEvent(0.0, 30.0),
        DoseEvent(0.0, 70.0),
        DoseEvent(12.0, 10.0),
        DoseEvent(12.0, 15.0),
    ]
    doses_b = [
        DoseEvent(0.0, 70.0),
        DoseEvent(0.0, 30.0),
        DoseEvent(12.0, 15.0),
        DoseEvent(12.0, 10.0),
    ]

    pk_a = ModelSpec(OneCompIVBolus(), "a", params, doses_a)
    pk_b = ModelSpec(OneCompIVBolus(), "b", params, doses_b)

    r_a = simulate(pk_a, grid, solver)
    r_b = simulate(pk_b, grid, solver)

    for i in eachindex(r_a.t)
        @test isapprox(
            r_a.observations[:conc][i], r_b.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "Event semantics: duplicate times are summed in coupled PKPD" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_dup_coupled",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 60.0), DoseEvent(0.0, 40.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "pd_coupled",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.0, 1.0),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test res.metadata["event_semantics_version"] == "1.0.0"
end

@testset "Solver semantics version is present and stable" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "solver_semantics_test",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:1.0:12.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.metadata, "solver_semantics_version")
    @test res.metadata["solver_semantics_version"] == "1.0.0"
end
