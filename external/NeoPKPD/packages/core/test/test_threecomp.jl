# Three-compartment PK model tests

@testset "ThreeCompIVBolus basic simulation" begin
    spec = ModelSpec(
        ThreeCompIVBolus(),
        "3c_iv_bolus",
        ThreeCompIVBolusParams(10.0, 50.0, 10.0, 80.0, 2.0, 200.0),  # CL=10, V1=50, Q2=10, V2=80, Q3=2, V3=200
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_central)
    @test haskey(res.states, :A_periph1)
    @test haskey(res.states, :A_periph2)
    @test haskey(res.observations, :conc)

    # Initial concentration should be Dose/V1
    @test isapprox(res.observations[:conc][1], 1000.0 / 50.0; rtol=0.01)

    # Metadata should be correct
    @test res.metadata["model"] == "ThreeCompIVBolus"
end

@testset "ThreeCompIVBolus shows tri-exponential decline" begin
    spec = ModelSpec(
        ThreeCompIVBolus(),
        "3c_triexp",
        ThreeCompIVBolusParams(5.0, 50.0, 20.0, 100.0, 1.0, 500.0),
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 200.0, collect(0.0:1.0:200.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    conc = res.observations[:conc]

    # Concentration should be monotonically decreasing (or nearly so)
    for i in 2:length(conc)
        @test conc[i] <= conc[i-1] + 1e-10
    end
end
