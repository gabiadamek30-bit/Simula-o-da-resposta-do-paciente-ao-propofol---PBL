# Population simulation tests

@testset "Population simulation is deterministic with StableRNG seed" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "pop_pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(12345), 20)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    r1 = simulate_population(pop, grid, solver)
    r2 = simulate_population(pop, grid, solver)

    @test r1.metadata["n"] == 20
    @test r2.metadata["n"] == 20

    # Realized individual parameters must match exactly run-to-run
    for i in 1:20
        @test r1.params[i][:CL] == r2.params[i][:CL]
        @test r1.params[i][:V] == r2.params[i][:V]
    end

    # Concentration series must match exactly to tight tolerance
    for i in 1:20
        c1 = r1.individuals[i].observations[:conc]
        c2 = r2.individuals[i].observations[:conc]
        @test length(c1) == length(c2)
        for j in eachindex(c1)
            @test isapprox(c1[j], c2[j]; rtol=1e-12, atol=1e-12)
        end
    end
end

@testset "Population simulation without IIV produces one individual matching base simulation" begin
    base = ModelSpec(
        OneCompOralFirstOrder(),
        "pop_pk_oral_base",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pop = PopulationSpec(base, nothing, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    base_res = simulate(base, grid, solver)
    pop_res = simulate_population(pop, grid, solver)

    @test length(pop_res.individuals) == 1

    c_base = base_res.observations[:conc]
    c_pop = pop_res.individuals[1].observations[:conc]

    for j in eachindex(c_base)
        @test isapprox(c_pop[j], c_base[j]; rtol=1e-12, atol=1e-12)
    end
end

@testset "Population artifact serialization and replay matches original outputs" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "pop_artifact_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(424242), 5)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    r1 = simulate_population(pop, grid, solver)

    artifact = serialize_population_execution(
        population_spec=pop, grid=grid, solver=solver, result=r1
    )

    r2 = replay_population_execution(artifact)

    @test r2.metadata["n"] == r1.metadata["n"]
    @test length(r2.individuals) == length(r1.individuals)

    for i in eachindex(r1.individuals)
        c1 = r1.individuals[i].observations[:conc]
        c2 = r2.individuals[i].observations[:conc]
        @test length(c1) == length(c2)
        for j in eachindex(c1)
            @test isapprox(c2[j], c1[j]; rtol=1e-12, atol=1e-12)
        end
    end

    @test haskey(r1.summaries, :conc)
    @test haskey(r2.summaries, :conc)

    s1 = r1.summaries[:conc]
    s2 = r2.summaries[:conc]

    @test s1.probs == s2.probs

    for i in eachindex(s1.mean)
        @test isapprox(s2.mean[i], s1.mean[i]; rtol=1e-12, atol=1e-12)
        @test isapprox(s2.median[i], s1.median[i]; rtol=1e-12, atol=1e-12)
    end

    for p in s1.probs
        q1 = s1.quantiles[p]
        q2 = s2.quantiles[p]
        for i in eachindex(q1)
            @test isapprox(q2[i], q1[i]; rtol=1e-12, atol=1e-12)
        end
    end
end
