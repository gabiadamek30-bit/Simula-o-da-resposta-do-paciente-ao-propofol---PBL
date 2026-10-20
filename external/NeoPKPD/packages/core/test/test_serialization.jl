# Serialization and artifact replay tests

@testset "Execution artifact serialization is complete and self-consistent" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "serialize_test",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:1.0:12.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    res = simulate(pk, grid, solver)

    artifact = serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)

    @test artifact["artifact_schema_version"] == ARTIFACT_SCHEMA_VERSION
    @test haskey(artifact, "model_spec")
    @test haskey(artifact, "grid")
    @test haskey(artifact, "solver")
    @test haskey(artifact, "result")

    @test artifact["result"]["metadata"]["event_semantics_version"] == "1.0.0"
    @test artifact["result"]["metadata"]["solver_semantics_version"] == "1.0.0"

    @test artifact["model_spec"]["name"] == "serialize_test"
    @test artifact["solver"]["alg"] == "Tsit5"
end

@testset "Artifact replay: PK only matches original outputs" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "replay_pk",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 60.0), DoseEvent(0.0, 40.0), DoseEvent(12.0, 25.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res1 = simulate(pk, grid, solver)

    artifact = serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res1)

    res2 = replay_execution(artifact)

    @test res2.t == res1.t
    for i in eachindex(res1.t)
        @test isapprox(
            res2.observations[:conc][i], res1.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "Artifact replay: PK then PD (DirectEmax) matches original outputs" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "replay_pkpd_direct",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "emax", DirectEmaxParams(10.0, 40.0, 0.8), :conc, :effect)

    res1 = simulate_pkpd(pk, pd, grid, solver)

    artifact = serialize_execution(
        model_spec=pk, grid=grid, solver=solver, result=res1, pd_spec=pd
    )

    res2 = replay_execution(artifact)

    @test res2.t == res1.t
    for i in eachindex(res1.t)
        @test isapprox(
            res2.observations[:conc][i], res1.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            res2.observations[:effect][i],
            res1.observations[:effect][i];
            rtol=1e-12,
            atol=1e-12,
        )
    end
end

@testset "Artifact replay: coupled PKPD (IndirectResponseTurnover) matches original outputs" begin
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "replay_coupled_oral",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    res1 = simulate_pkpd_coupled(pk, pd, grid, solver)

    artifact = serialize_execution(
        model_spec=pk, grid=grid, solver=solver, result=res1, pd_spec=pd
    )

    res2 = replay_execution(artifact)

    @test res2.t == res1.t
    for i in eachindex(res1.t)
        @test isapprox(
            res2.observations[:conc][i], res1.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            res2.observations[:response][i],
            res1.observations[:response][i];
            rtol=1e-12,
            atol=1e-12,
        )
    end
end

@testset "Semantics fingerprint is complete" begin
    fp = semantics_fingerprint()

    @test haskey(fp, "event_semantics_version")
    @test haskey(fp, "solver_semantics_version")
    @test haskey(fp, "artifact_schema_version")
end
