# Transit absorption model tests

@testset "TransitAbsorption basic simulation" begin
    spec = ModelSpec(
        TransitAbsorption(),
        "transit_abs",
        TransitAbsorptionParams(5, 2.0, 1.0, 10.0, 50.0),  # N=5, Ktr=2, Ka=1, CL=10, V=50
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_central)
    @test haskey(res.states, :Transit_1)
    @test haskey(res.states, :Transit_5)
    @test haskey(res.observations, :conc)

    # Initial concentration should be zero
    @test isapprox(res.observations[:conc][1], 0.0; atol=1e-10)

    # Delayed peak due to transit chain
    conc = res.observations[:conc]
    max_idx = argmax(conc)
    tmax = res.t[max_idx]
    @test tmax > 0.5  # Peak should be delayed

    # Metadata should be correct
    @test res.metadata["model"] == "TransitAbsorption"
    @test res.metadata["N_transit"] == 5
end

@testset "TransitAbsorption: more compartments delay peak" begin
    grid = SimGrid(0.0, 24.0, collect(0.0:0.1:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    spec_n2 = ModelSpec(
        TransitAbsorption(),
        "transit_n2",
        TransitAbsorptionParams(2, 2.0, 1.0, 10.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    spec_n8 = ModelSpec(
        TransitAbsorption(),
        "transit_n8",
        TransitAbsorptionParams(8, 2.0, 1.0, 10.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    res_n2 = simulate(spec_n2, grid, solver)
    res_n8 = simulate(spec_n8, grid, solver)

    tmax_n2 = res_n2.t[argmax(res_n2.observations[:conc])]
    tmax_n8 = res_n8.t[argmax(res_n8.observations[:conc])]

    # More transit compartments should delay the peak
    @test tmax_n8 > tmax_n2
end
