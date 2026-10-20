# Integration tests: PK models with PD models

@testset "TwoCompIVBolus with DirectEmax PKPD" begin
    pk = ModelSpec(
        TwoCompIVBolus(),
        "2c_pkpd",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "emax", DirectEmaxParams(0.0, 100.0, 2.0), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)

    # Effect should track concentration
    for (i, t) in enumerate(res.t)
        C = res.observations[:conc][i]
        e_ref = direct_emax(C, 0.0, 100.0, 2.0)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10)
    end
end

@testset "TransitAbsorption with SigmoidEmax PKPD" begin
    pk = ModelSpec(
        TransitAbsorption(),
        "transit_pkpd",
        TransitAbsorptionParams(5, 2.0, 1.0, 10.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(SigmoidEmax(), "semax", SigmoidEmaxParams(0.0, 100.0, 2.0, 2.0), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :effect)

    # Effect should follow sigmoid Emax of concentration
    for (i, t) in enumerate(res.t)
        C = res.observations[:conc][i]
        e_ref = sigmoid_emax_ref(C, 0.0, 100.0, 2.0, 2.0)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10)
    end
end

@testset "MichaelisMentenElimination with DirectEmax PKPD" begin
    pk = ModelSpec(
        MichaelisMentenElimination(),
        "mm_pkpd",
        MichaelisMentenEliminationParams(100.0, 5.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "emax", DirectEmaxParams(10.0, 50.0, 3.0), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :effect)

    for (i, t) in enumerate(res.t)
        C = res.observations[:conc][i]
        e_ref = direct_emax(C, 10.0, 50.0, 3.0)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10)
    end
end

@testset "ThreeCompIVBolus with BiophaseEquilibration PKPD" begin
    pk = ModelSpec(
        ThreeCompIVBolus(),
        "3c_biophase",
        ThreeCompIVBolusParams(10.0, 50.0, 10.0, 80.0, 2.0, 200.0),
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(BiophaseEquilibration(), "biophase", BiophaseEquilibrationParams(0.5, 0.0, 100.0, 5.0), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)
end
