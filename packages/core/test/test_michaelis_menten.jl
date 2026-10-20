# Michaelis-Menten elimination model tests

@testset "MichaelisMentenElimination basic simulation" begin
    spec = ModelSpec(
        MichaelisMentenElimination(),
        "mm_elim",
        MichaelisMentenEliminationParams(100.0, 5.0, 50.0),  # Vmax=100, Km=5, V=50
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_central)
    @test haskey(res.observations, :conc)

    # Initial concentration should be Dose/V
    @test isapprox(res.observations[:conc][1], 500.0 / 50.0; rtol=0.01)

    # Concentration should decline
    @test res.observations[:conc][end] < res.observations[:conc][1]

    # Metadata should be correct
    @test res.metadata["model"] == "MichaelisMentenElimination"
end

@testset "MichaelisMentenElimination shows nonlinear kinetics" begin
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    # Low dose - approximately linear
    spec_low = ModelSpec(
        MichaelisMentenElimination(),
        "mm_low",
        MichaelisMentenEliminationParams(100.0, 10.0, 50.0),  # Km=10
        [DoseEvent(0.0, 50.0)],  # C0 = 1, << Km
    )

    # High dose - saturated elimination
    spec_high = ModelSpec(
        MichaelisMentenElimination(),
        "mm_high",
        MichaelisMentenEliminationParams(100.0, 10.0, 50.0),
        [DoseEvent(0.0, 5000.0)],  # C0 = 100, >> Km
    )

    res_low = simulate(spec_low, grid, solver)
    res_high = simulate(spec_high, grid, solver)

    # AUC should increase disproportionately with dose
    auc_low = auc_trapezoid(res_low.t, res_low.observations[:conc])
    auc_high = auc_trapezoid(res_high.t, res_high.observations[:conc])

    # AUC ratio should be > dose ratio (100x) due to saturable kinetics
    @test auc_high / auc_low > 100.0
end
