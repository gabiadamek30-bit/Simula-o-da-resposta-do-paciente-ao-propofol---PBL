# Two-compartment PK model tests

@testset "TwoCompIVBolus basic simulation" begin
    spec = ModelSpec(
        TwoCompIVBolus(),
        "2c_iv_bolus",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),  # CL=10, V1=50, Q=5, V2=100
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_central)
    @test haskey(res.states, :A_peripheral)
    @test haskey(res.observations, :conc)
    @test length(res.observations[:conc]) == length(res.t)

    # Initial concentration should be Dose/V1
    @test isapprox(res.observations[:conc][1], 500.0 / 50.0; rtol=0.01)

    # Concentration should decline over time
    @test res.observations[:conc][end] < res.observations[:conc][1]

    # Metadata should be correct
    @test res.metadata["model"] == "TwoCompIVBolus"
end

@testset "TwoCompIVBolus mass distribution" begin
    spec = ModelSpec(
        TwoCompIVBolus(),
        "2c_iv_mass",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),  # Normal elimination
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 100.0, collect(0.0:1.0:100.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    # Total mass should decline over time due to elimination
    initial_total = res.states[:A_central][1] + res.states[:A_peripheral][1]
    final_total = res.states[:A_central][end] + res.states[:A_peripheral][end]
    @test final_total < initial_total

    # Mass should be distributed to peripheral at equilibrium
    # Initially all drug in central, later some moves to peripheral
    @test res.states[:A_peripheral][1] ≈ 0.0 atol=1e-10
    @test res.states[:A_peripheral][10] > 0.0  # After some time, peripheral has drug
end

@testset "TwoCompOral basic simulation" begin
    spec = ModelSpec(
        TwoCompOral(),
        "2c_oral",
        TwoCompOralParams(1.0, 10.0, 50.0, 5.0, 100.0),  # Ka=1, CL=10, V1=50, Q=5, V2=100
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_gut)
    @test haskey(res.states, :A_central)
    @test haskey(res.states, :A_peripheral)
    @test haskey(res.observations, :conc)

    # Initial concentration should be zero (drug in gut)
    @test isapprox(res.observations[:conc][1], 0.0; atol=1e-10)

    # Concentration should rise then fall
    conc = res.observations[:conc]
    max_idx = argmax(conc)
    @test max_idx > 1  # Peak occurs after t=0

    # Metadata should be correct
    @test res.metadata["model"] == "TwoCompOral"
end
