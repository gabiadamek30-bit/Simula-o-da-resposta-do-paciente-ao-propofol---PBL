# Biophase equilibration PD model tests

@testset "BiophaseEquilibration basic evaluation (quasi-steady state)" begin
    pd = PDSpec(
        BiophaseEquilibration(),
        "biophase",
        BiophaseEquilibrationParams(0.5, 10.0, 40.0, 0.8),  # ke0=0.5, E0=10, Emax=40, EC50=0.8
        :conc,
        :effect,
    )

    concentrations = [0.0, 0.5, 1.0, 2.0, 5.0]
    effects = evaluate(pd, concentrations)

    # Quasi-steady state: Ce = Cp, so effect is direct Emax
    for (i, C) in enumerate(concentrations)
        e_ref = direct_emax(C, 10.0, 40.0, 0.8)
        @test isapprox(effects[i], e_ref; rtol=1e-10)
    end
end

@testset "BiophaseEquilibration PKPD coupling" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_biophase",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        BiophaseEquilibration(),
        "biophase_test",
        BiophaseEquilibrationParams(0.5, 10.0, 40.0, 0.8),
        :conc,
        :effect,
    )

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)
    @test res.metadata["pd_model"] == "BiophaseEquilibration"
end
